#!/usr/bin/env Rscript

# Initialize and validate the normalized DuckVEP annotation-model relation pack.
#
# The SQL file beside this script is the physical DuckDB/Parquet contract.
# Validation always emits one deterministic JSON document.  Exit status is 0
# for the selected profile passing, 1 for an invalid pack, and 2 for invocation
# or tool failure.

CONTRACT <- "duckvep_model"
SCHEMA_NAME <- "duckvep_model"
FORMAT_VERSION <- 2L
VALIDATOR_VERSION <- 1L
COMPATIBILITY_TARGET <- "ensembl-vep/116.0"
COMPATIBILITY_COMMIT <- "57ea5c52340acc1f156267f810ad162e26597082"
ENSEMBL_CORE_COMMIT <- "c0cf13daa961d80584bad797b2eb0ff3a7500ef3"
ENSEMBL_CORE_DDL_SHA256 <- "ec7bb8dd2fcd6a7012bfd67aff8065b912915d541450fa4f3f05334a92944e8c"
PROFILES <- c("structural", "model_candidate", "conformance_gate")

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    return(normalizePath("scripts/duckvep_model.R", mustWork = FALSE))
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
}

SCHEMA_PATH <- file.path(dirname(script_path()), "duckvep_model_schema.sql")

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("required R package is not installed: %s", package), call. = FALSE)
  }
}

for (package in c("DBI", "duckdb", "jsonlite")) {
  require_package(package)
}

sql_ident <- function(value) {
  paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
}

relation_sql <- function(name) {
  paste0(sql_ident(SCHEMA_NAME), ".", sql_ident(name))
}

sql_literal <- function(value) {
  paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

read_schema <- function() {
  if (!file.exists(SCHEMA_PATH)) {
    stop(sprintf("schema SQL does not exist: %s", SCHEMA_PATH), call. = FALSE)
  }
  paste(readLines(SCHEMA_PATH, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

json_scalar <- function(value) {
  if (inherits(value, "integer64")) {
    return(as.character(value))
  }
  if (inherits(value, "POSIXt")) {
    if (is.na(value)) return(NA_character_)
    return(format(value, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  }
  if (is.raw(value)) {
    return(paste(format(value), collapse = ""))
  }
  if (is.list(value) && length(value) == 1L) {
    return(json_scalar(value[[1L]]))
  }
  if (!length(value)) return(NA_character_)
  value[[1L]]
}

rows_as_list <- function(data) {
  if (!nrow(data)) return(list())
  lapply(seq_len(nrow(data)), function(row_index) {
    row <- lapply(data, function(column) json_scalar(column[row_index]))
    names(row) <- names(data)
    row
  })
}

emit_json <- function(payload, pretty = FALSE) {
  text <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    digits = NA,
    pretty = pretty
  )
  cat(text, "\n", sep = "")
}

tool_error <- function(code, message) {
  condition <- structure(
    list(message = message, call = NULL, code = code),
    class = c("duckvep_model_tool_error", "error", "condition")
  )
  stop(condition)
}

connect_duckdb <- function(dbdir = ":memory:", read_only = FALSE) {
  DBI::dbConnect(duckdb::duckdb(), dbdir = dbdir, read_only = read_only)
}

disconnect_duckdb <- function(connection) {
  try(DBI::dbDisconnect(connection, shutdown = TRUE), silent = TRUE)
}

expected_schema <- function() {
  connection <- connect_duckdb()
  on.exit(disconnect_duckdb(connection), add = TRUE)
  DBI::dbExecute(connection, read_schema())

  tables <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT table_name FROM information_schema.tables",
      "WHERE table_schema = ? AND table_type = 'BASE TABLE'",
      "ORDER BY table_name"
    ),
    params = list(SCHEMA_NAME)
  )$table_name

  columns <- setNames(vector("list", length(tables)), tables)
  for (table in tables) {
    columns[[table]] <- DBI::dbGetQuery(
      connection,
      paste(
        "SELECT column_name, data_type, is_nullable = 'YES' AS nullable",
        "FROM information_schema.columns",
        "WHERE table_schema = ? AND table_name = ?",
        "ORDER BY ordinal_position"
      ),
      params = list(SCHEMA_NAME, table)
    )
  }

  constraints <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT table_name, constraint_column_names",
      "FROM duckdb_constraints()",
      "WHERE schema_name = ? AND constraint_type = 'PRIMARY KEY'",
      "ORDER BY table_name"
    ),
    params = list(SCHEMA_NAME)
  )
  keys <- setNames(vector("list", nrow(constraints)), constraints$table_name)
  if (nrow(constraints)) {
    for (index in seq_len(nrow(constraints))) {
      keys[[constraints$table_name[[index]]]] <-
        as.character(constraints$constraint_column_names[[index]])
    }
  }
  list(columns = columns, keys = keys)
}

new_validator <- function(connection, expected, profile, max_samples, pack_kind) {
  validator <- new.env(parent = emptyenv())
  validator$connection <- connection
  validator$expected <- expected
  validator$profile <- profile
  validator$max_samples <- max_samples
  validator$pack_kind <- pack_kind
  validator$issues <- list()
  validator$ready <- setNames(rep(FALSE, length(expected$columns)), names(expected$columns))
  validator$row_counts <- setNames(as.list(rep(NA_real_, length(expected$columns))),
                                   names(expected$columns))
  validator
}

add_issue <- function(validator, code, relation, message, count = 1, sample = list()) {
  validator$issues[[length(validator$issues) + 1L]] <- list(
    code = code,
    count = as.numeric(count),
    message = message,
    relation = relation,
    sample = sample
  )
}

check_query <- function(validator, code, relation, message, query, order_by) {
  count_query <- paste0("SELECT count(*) AS n FROM (", query, ") AS failures")
  count <- DBI::dbGetQuery(validator$connection, count_query)$n[[1L]]
  if (as.numeric(count) == 0) return(invisible(NULL))

  sample <- list()
  if (validator$max_samples > 0L) {
    sample_query <- paste0(
      "SELECT * FROM (", query, ") AS failures ORDER BY ", order_by,
      " LIMIT ", validator$max_samples
    )
    sample <- rows_as_list(DBI::dbGetQuery(validator$connection, sample_query))
  }
  add_issue(validator, code, relation, message, count, sample)
  invisible(NULL)
}

relations_ready <- function(validator, relations) {
  all(vapply(relations, function(name) isTRUE(validator$ready[[name]]), logical(1)))
}

inspect_shape <- function(validator) {
  actual_relations <- DBI::dbGetQuery(
    validator$connection,
    paste(
      "SELECT table_name,table_type FROM information_schema.tables",
      "WHERE table_schema = ? ORDER BY table_name"
    ),
    params = list(SCHEMA_NAME)
  )
  actual <- actual_relations$table_name

  unexpected <- setdiff(actual, names(validator$expected$columns))
  for (table in unexpected) {
    add_issue(
      validator, "DVMA1000", table,
      sprintf("unexpected relation %s.%s is not in format version %d",
              SCHEMA_NAME, table, FORMAT_VERSION)
    )
  }

  for (table in names(validator$expected$columns)) {
    expected_columns <- validator$expected$columns[[table]]
    if (!table %in% actual) {
      add_issue(
        validator, "DVMA1001", table,
        sprintf("required relation %s.%s is missing", SCHEMA_NAME, table)
      )
      next
    }

    actual_columns <- DBI::dbGetQuery(
      validator$connection,
      paste(
        "SELECT column_name, data_type FROM information_schema.columns",
        "WHERE table_schema = ? AND table_name = ? ORDER BY ordinal_position"
      ),
      params = list(SCHEMA_NAME, table)
    )
    good <- TRUE
    expected_table_type <- if (validator$pack_kind == "duckdb_database")
      "BASE TABLE" else "VIEW"
    table_type <- actual_relations$table_type[match(table, actual_relations$table_name)]
    if (!identical(table_type, expected_table_type)) {
      good <- FALSE
      add_issue(
        validator, "DVMA1006", table,
        sprintf(
          "%s pack relation has table_type %s; expected %s",
          validator$pack_kind, table_type, expected_table_type
        )
      )
    }
    for (index in seq_len(nrow(expected_columns))) {
      column <- expected_columns$column_name[[index]]
      match_index <- match(column, actual_columns$column_name)
      if (is.na(match_index)) {
        good <- FALSE
        add_issue(validator, "DVMA1002", table,
                  sprintf("required column %s is missing", column))
      } else if (actual_columns$data_type[[match_index]] !=
                 expected_columns$data_type[[index]]) {
        good <- FALSE
        add_issue(
          validator, "DVMA1003", table,
          sprintf("column %s has type %s; expected %s", column,
                  actual_columns$data_type[[match_index]],
                  expected_columns$data_type[[index]])
        )
      }
    }
    extra <- setdiff(actual_columns$column_name, expected_columns$column_name)
    for (column in extra) {
      good <- FALSE
      add_issue(validator, "DVMA1004", table,
                sprintf("unexpected column %s is not in format version %d",
                        column, FORMAT_VERSION))
    }
    if (!identical(actual_columns$column_name, expected_columns$column_name)) {
      good <- FALSE
      add_issue(
        validator, "DVMA1005", table,
        "column order differs from the versioned schema",
        sample = list(list(
          actual = as.list(actual_columns$column_name),
          expected = as.list(expected_columns$column_name)
        ))
      )
    }
    validator$ready[[table]] <- good
    validator$row_counts[[table]] <- as.numeric(DBI::dbGetQuery(
      validator$connection,
      sprintf("SELECT count(*) AS n FROM %s", relation_sql(table))
    )$n[[1L]])
  }
}

check_generic_constraints <- function(validator) {
  for (table in names(validator$expected$columns)) {
    if (!isTRUE(validator$ready[[table]])) next
    columns <- validator$expected$columns[[table]]
    keys <- validator$expected$keys[[table]]
    if (is.null(keys) || !length(keys)) {
      add_issue(validator, "DVMA1100", table,
                "versioned relation has no declared primary identity")
      next
    }
    key_sql <- paste(vapply(keys, sql_ident, character(1)), collapse = ", ")
    required <- columns$column_name[!columns$nullable]
    if (length(required)) {
      predicate <- paste(
        paste0(vapply(required, sql_ident, character(1)), " IS NULL"),
        collapse = " OR "
      )
      check_query(
        validator, "DVMA1101", table, "NOT NULL contract violated",
        sprintf("SELECT %s FROM %s WHERE %s", key_sql, relation_sql(table), predicate),
        key_sql
      )
    }
    check_query(
      validator, "DVMA1102", table, "primary identity is not unique",
      sprintf(
        "SELECT %s, count(*) AS row_count FROM %s GROUP BY %s HAVING count(*) > 1",
        key_sql, relation_sql(table), key_sql
      ),
      key_sql
    )
  }
}

model_relations <- function(expected) {
  names(expected$columns)[vapply(expected$columns, function(columns) {
    "model_id" %in% columns$column_name
  }, logical(1))]
}

check_manifest_sources_build <- function(validator) {
  strict <- validator$profile != "structural"
  if (validator$ready[["model_manifest"]]) {
    comparator <- if (strict) "<> 1" else "> 1"
    check_query(
      validator, "DVMA2001", "model_manifest",
      if (strict) "a model candidate must contain exactly one manifest row" else
        "a structural pack may contain at most one manifest row",
      sprintf(
        "SELECT count(*) AS row_count FROM %s HAVING count(*) %s",
        relation_sql("model_manifest"), comparator
      ),
      "row_count"
    )
    check_query(
      validator, "DVMA2002", "model_manifest",
      "manifest is not pinned to the annotation-model format and VEP target",
      sprintf(
        paste(
          "SELECT model_id, format_version, compatibility_target,",
          "compatibility_commit, canonicalization_version",
          "FROM %s WHERE format_version <> %d",
          "OR compatibility_target <> %s OR compatibility_commit <> %s"
        ),
        relation_sql("model_manifest"), FORMAT_VERSION,
        sql_literal(COMPATIBILITY_TARGET), sql_literal(COMPATIBILITY_COMMIT)
      ),
      "model_id"
    )
    check_query(
      validator, "DVMA2003", "model_manifest",
      "manifest identifiers and provenance strings must be non-empty canonical values",
      sprintf(
        paste(
          "SELECT model_id FROM %s WHERE",
          "NOT regexp_full_match(model_id, '[0-9a-f]{64}')",
          "OR length(trim(species)) = 0 OR taxonomy_id = 0",
          "OR length(trim(assembly_name)) = 0",
          "OR length(trim(assembly_accession)) = 0",
          "OR length(trim(selection_policy_version)) = 0",
          "OR length(trim(canonicalization_version)) = 0"
        ),
        relation_sql("model_manifest")
      ),
      "model_id"
    )
  }

  if (relations_ready(validator, c("model_manifest", "model_source"))) {
    check_query(
      validator, "DVMA2010", "model_source",
      "source namespace has invalid provenance, kind, or priority",
      sprintf(
        paste(
          "SELECT model_id, source_namespace FROM %s WHERE",
          "length(trim(source_namespace)) = 0 OR length(trim(provider)) = 0",
          "OR source_kind NOT IN ('core','refseq','merged','overlay','crosswalk')",
          "OR length(trim(selection_role)) = 0 OR length(trim(annotation_set)) = 0",
          "OR length(trim(annotation_source_release)) = 0 OR source_priority < 0",
          "OR (schema_sha256 IS NOT NULL AND",
          "NOT regexp_full_match(schema_sha256, '[0-9a-f]{64}'))",
          "OR ((schema_repository IS NULL) <> (schema_commit IS NULL))",
          "OR ((schema_repository IS NULL) <> (schema_sha256 IS NULL))",
          "OR (source_kind = 'core' AND (source_database IS NULL",
          "OR schema_repository IS NULL))"
        ),
        relation_sql("model_source")
      ),
      "model_id, source_namespace"
    )
    check_query(
      validator, "DVMA2011", "model_source",
      "source priorities must be unique within one model",
      sprintf(
        paste(
          "SELECT model_id, source_priority, count(*) AS source_count FROM %s",
          "GROUP BY model_id, source_priority HAVING count(*) > 1"
        ),
        relation_sql("model_source")
      ),
      "model_id, source_priority"
    )
    check_query(
      validator, "DVMA2012", "model_source",
      "Ensembl core source is not pinned to the VEP-116 API/schema receipt",
      sprintf(
        paste(
          "SELECT model_id,source_namespace,schema_commit,schema_sha256 FROM %s",
          "WHERE provider='Ensembl' AND source_kind='core'",
          "AND (schema_commit<>%s OR schema_sha256<>%s)"
        ),
        relation_sql("model_source"), sql_literal(ENSEMBL_CORE_COMMIT),
        sql_literal(ENSEMBL_CORE_DDL_SHA256)
      ),
      "model_id, source_namespace"
    )
  }

  if (validator$ready[["model_build"]]) {
    comparator <- if (strict) "<> 1" else "> 1"
    check_query(
      validator, "DVMA2020", "model_build",
      if (strict) "a model candidate must select exactly one build receipt" else
        "a structural pack may contain at most one build receipt",
      sprintf(
        "SELECT count(*) AS row_count FROM %s HAVING count(*) %s",
        relation_sql("model_build"), comparator
      ),
      "row_count"
    )
    check_query(
      validator, "DVMA2021", "model_build",
      "build receipt hashes or repositories are malformed",
      sprintf(
        paste(
          "SELECT model_id, model_build_id FROM %s WHERE",
          "NOT regexp_full_match(model_build_id, '[0-9a-f]{64}')",
          "OR NOT regexp_full_match(importer_sha256, '[0-9a-f]{64}')",
          "OR NOT regexp_full_match(exporter_sha256, '[0-9a-f]{64}')",
          "OR NOT regexp_full_match(invocation_sha256, '[0-9a-f]{64}')",
          "OR NOT regexp_full_match(environment_sha256, '[0-9a-f]{64}')",
          "OR NOT regexp_full_match(importer_commit, '[0-9a-f]{40,64}')",
          "OR NOT regexp_full_match(exporter_commit, '[0-9a-f]{40,64}')",
          "OR length(trim(importer_repository)) = 0",
          "OR length(trim(exporter_repository)) = 0"
        ),
        relation_sql("model_build")
      ),
      "model_id, model_build_id"
    )
  }

  if (relations_ready(validator, c("model_manifest", "model_build"))) {
    check_query(
      validator, "DVMA2022", "model_build",
      "build receipt does not belong to the sole manifest model",
      sprintf(
        paste(
          "SELECT b.model_id, b.model_build_id FROM %s b",
          "LEFT JOIN %s m USING (model_id) WHERE m.model_id IS NULL"
        ),
        relation_sql("model_build"), relation_sql("model_manifest")
      ),
      "model_id, model_build_id"
    )
  }

  if (relations_ready(validator, c("model_source", "model_selection_audit"))) {
    check_query(
      validator, "DVMA2030", "model_selection_audit",
      "selection audit row has no source or has invalid categorical data",
      sprintf(
        paste(
          "SELECT a.model_id, a.source_namespace, a.object_kind, a.decision,",
          "a.reason_code FROM %s a LEFT JOIN %s s",
          "USING (model_id, source_namespace) WHERE s.source_namespace IS NULL",
          "OR length(trim(a.object_kind)) = 0",
          "OR a.decision NOT IN ('included','excluded','quarantined')",
          "OR length(trim(a.reason_code)) = 0"
        ),
        relation_sql("model_selection_audit"), relation_sql("model_source")
      ),
      "model_id, source_namespace, object_kind, decision, reason_code"
    )
    check_query(
      validator, "DVMA2033", "model_selection_audit",
      "one selection reason cannot carry multiple mutually exclusive decisions",
      sprintf(
        paste(
          "SELECT model_id,source_namespace,object_kind,reason_code,",
          "count(DISTINCT decision) AS decision_count FROM %s",
          "GROUP BY model_id,source_namespace,object_kind,reason_code",
          "HAVING count(DISTINCT decision)>1"
        ),
        relation_sql("model_selection_audit")
      ),
      "model_id, source_namespace, object_kind, reason_code"
    )
    if (validator$profile != "structural") {
      check_query(
        validator, "DVMA2031", "model_selection_audit",
        "every declared source namespace requires explicit selection-audit coverage",
        sprintf(
          paste(
            "SELECT s.model_id,s.source_namespace FROM %s s LEFT JOIN %s a",
            "USING (model_id,source_namespace) GROUP BY s.model_id,s.source_namespace",
            "HAVING count(a.object_kind)=0"
          ),
          relation_sql("model_source"), relation_sql("model_selection_audit")
        ),
        "model_id, source_namespace"
      )
    }
  }

  selection_relations <- list(
    seq_region = c("model_seq_region", "source_namespace"),
    gene = c("model_gene", "source_namespace"),
    transcript = c("model_transcript", "source_namespace"),
    exon = c("model_exon", "source_namespace"),
    transcript_exon = c("model_transcript_exon", "source_namespace"),
    translation = c("model_translation", "source_namespace"),
    attribute_type = c("model_attribute_type", "source_namespace"),
    seq_region_attribute = c("model_seq_region_attribute", "source_namespace"),
    gene_attribute = c("model_gene_attribute", "source_namespace"),
    transcript_attribute = c("model_transcript_attribute", "source_namespace"),
    translation_attribute = c("model_translation_attribute", "source_namespace"),
    external_db = c("model_external_db", "source_namespace"),
    xref = c("model_xref", "source_namespace"),
    object_xref = c("model_object_xref", "object_source_namespace"),
    seq_region_synonym = c("model_seq_region_synonym", "source_namespace")
  )
  selection_requirements <- c(
    "model_source", "model_selection_audit",
    vapply(selection_relations, `[[`, character(1), 1L)
  )
  if (relations_ready(validator, selection_requirements)) {
    actual_queries <- vapply(names(selection_relations), function(kind) {
      spec <- selection_relations[[kind]]
      sprintf(
        paste(
          "SELECT model_id,%s AS source_namespace,%s AS object_kind,",
          "count(*) AS actual_count FROM %s GROUP BY model_id,%s"
        ),
        sql_ident(spec[[2L]]), sql_literal(kind), relation_sql(spec[[1L]]),
        sql_ident(spec[[2L]])
      )
    }, character(1))
    kinds <- paste(
      sprintf("(%s)", vapply(names(selection_relations), sql_literal, character(1))),
      collapse = ","
    )
    check_query(
      validator, "DVMA2032", "model_selection_audit",
      "included selection-audit row counts disagree with retained source relations",
      sprintf(
        paste(
          "WITH kinds(object_kind) AS (VALUES %s), actual AS (%s),",
          "audited AS (SELECT model_id,source_namespace,object_kind,sum(row_count) AS audited_count",
          "FROM %s WHERE decision='included'",
          "GROUP BY model_id,source_namespace,object_kind),",
          "expected AS (SELECT s.model_id,s.source_namespace,k.object_kind",
          "FROM %s s CROSS JOIN kinds k)",
          "SELECT e.model_id,e.source_namespace,e.object_kind,",
          "coalesce(a.actual_count,0) AS actual_count,",
          "coalesce(u.audited_count,0) AS audited_count FROM expected e",
          "LEFT JOIN actual a USING (model_id,source_namespace,object_kind)",
          "LEFT JOIN audited u USING (model_id,source_namespace,object_kind)",
          "WHERE coalesce(a.actual_count,0)<>coalesce(u.audited_count,0)"
        ),
        kinds, paste(actual_queries, collapse = " UNION ALL "),
        relation_sql("model_selection_audit"), relation_sql("model_source")
      ),
      "model_id, source_namespace, object_kind"
    )
  }
}

check_artifacts_and_pack_references <- function(validator) {
  if (validator$ready[["model_artifact"]]) {
    check_query(
      validator, "DVMA2040", "model_artifact",
      "artifact receipt has an invalid role, locator, hash, or count",
      sprintf(
        paste(
          "SELECT model_id, model_build_id, role_class, role, logical_name",
          "FROM %s WHERE role_class NOT IN",
          "('source_input','software_input','generated_output','audit_output','oracle_output')",
          "OR length(trim(role)) = 0 OR length(trim(logical_name)) = 0",
          "OR length(trim(locator)) = 0",
          "OR NOT regexp_full_match(sha256, '[0-9a-f]{64}')",
          "OR byte_count<0",
          "OR (role_class='generated_output' AND row_count IS NULL)"
        ),
        relation_sql("model_artifact")
      ),
      "model_id, model_build_id, role_class, role, logical_name"
    )
    check_query(
      validator, "DVMA2045", "model_artifact",
      "core_ddl artifact does not carry the exact pinned Ensembl-116 DDL digest",
      sprintf(
        "SELECT model_id,model_build_id,logical_name,sha256 FROM %s WHERE role='core_ddl' AND sha256<>%s",
        relation_sql("model_artifact"), sql_literal(ENSEMBL_CORE_DDL_SHA256)
      ),
      "model_id, model_build_id, logical_name"
    )
  }

  if (relations_ready(validator, c("model_artifact", "model_build"))) {
    check_query(
      validator, "DVMA2041", "model_artifact",
      "artifact does not reference the selected model/build receipt",
      sprintf(
        paste(
          "SELECT a.model_id, a.model_build_id, a.role_class, a.role, a.logical_name",
          "FROM %s a LEFT JOIN %s b USING (model_id, model_build_id)",
          "WHERE b.model_build_id IS NULL"
        ),
        relation_sql("model_artifact"), relation_sql("model_build")
      ),
      "model_id, model_build_id, role_class, role, logical_name"
    )
    check_query(
      validator, "DVMA2047", "model_artifact",
      "importer/exporter software artifact digest disagrees with model_build",
      sprintf(
        paste(
          "SELECT a.model_id,a.model_build_id,a.role,a.logical_name,a.sha256",
          "FROM %s a JOIN %s b USING (model_id,model_build_id)",
          "WHERE (a.role='importer' AND a.sha256<>b.importer_sha256)",
          "OR (a.role='sequence_state_exporter' AND a.sha256<>b.exporter_sha256)"
        ),
        relation_sql("model_artifact"), relation_sql("model_build")
      ),
      "model_id, model_build_id, role, logical_name"
    )
  }

  if (relations_ready(validator, c("model_artifact", "model_source"))) {
    check_query(
      validator, "DVMA2042", "model_artifact",
      "artifact source_namespace does not reference a declared source",
      sprintf(
        paste(
          "SELECT a.model_id, a.model_build_id, a.role, a.logical_name,",
          "a.source_namespace FROM %s a LEFT JOIN %s s",
          "ON s.model_id=a.model_id AND s.source_namespace=a.source_namespace",
          "WHERE a.source_namespace IS NOT NULL AND s.source_namespace IS NULL"
        ),
        relation_sql("model_artifact"), relation_sql("model_source")
      ),
      "model_id, model_build_id, role, logical_name"
    )
    check_query(
      validator, "DVMA2046", "model_artifact",
      "every declared source namespace requires a forward required source-input receipt",
      sprintf(
        paste(
          "SELECT s.model_id,s.source_namespace,count(a.logical_name) AS receipt_count",
          "FROM %s s LEFT JOIN %s a ON a.model_id=s.model_id",
          "AND a.source_namespace=s.source_namespace",
          "AND a.role_class='source_input' AND a.required",
          "GROUP BY s.model_id,s.source_namespace HAVING count(a.logical_name)=0"
        ),
        relation_sql("model_source"), relation_sql("model_artifact")
      ),
      "model_id, source_namespace"
    )
    check_query(
      validator, "DVMA2048", "model_artifact",
      "core_ddl artifact digest disagrees with its model_source schema receipt",
      sprintf(
        paste(
          "SELECT s.model_id,s.source_namespace,s.schema_sha256,a.sha256 AS artifact_sha256",
          "FROM %s s LEFT JOIN %s a ON a.model_id=s.model_id",
          "AND a.source_namespace=s.source_namespace AND a.role_class='source_input'",
          "AND a.role='core_ddl' AND a.required",
          "WHERE s.source_kind='core' AND",
          "(a.logical_name IS NULL OR a.sha256<>s.schema_sha256)"
        ),
        relation_sql("model_source"), relation_sql("model_artifact")
      ),
      "model_id, source_namespace"
    )
  }

  if (relations_ready(validator, c(
      "model_selection_audit", "model_artifact", "model_build"
  ))) {
    check_query(
      validator, "DVMA2044", "model_selection_audit",
      "named row-level rejection ledger is not a checksummed model artifact",
      sprintf(
        paste(
          "SELECT a.model_id,a.source_namespace,a.object_kind,a.decision,a.reason_code,",
          "a.rejected_rows_artifact FROM %s a JOIN %s b ON b.model_id=a.model_id",
          "LEFT JOIN %s r ON r.model_id=a.model_id AND r.model_build_id=b.model_build_id",
          "AND r.logical_name=a.rejected_rows_artifact AND r.role_class='audit_output'",
          "AND r.source_namespace=a.source_namespace AND r.required",
          "WHERE a.rejected_rows_artifact IS NOT NULL AND r.logical_name IS NULL"
        ),
        relation_sql("model_selection_audit"), relation_sql("model_build"),
        relation_sql("model_artifact")
      ),
      "model_id, source_namespace, object_kind, decision, reason_code"
    )
  }

  if (validator$profile != "structural" &&
      relations_ready(validator, c("model_manifest", "model_artifact"))) {
    required <- data.frame(
      role_class = c(
        "source_input", "source_input", "source_input",
        "software_input", "software_input", "software_input"
      ),
      role = c(
        "core_ddl", "reference_fasta", "reference_fasta_index",
        "importer", "sequence_state_exporter", "vep_api_modules"
      ),
      stringsAsFactors = FALSE
    )
    values <- paste(
      sprintf("(%s,%s)", vapply(required$role_class, sql_literal, character(1)),
              vapply(required$role, sql_literal, character(1))),
      collapse = ","
    )
    check_query(
      validator, "DVMA2043", "model_artifact",
      "required annotation-model receipt role is absent",
      sprintf(
        paste(
          "WITH required(role_class, role) AS (VALUES %s),",
          "models AS (SELECT model_id FROM %s)",
          "SELECT m.model_id, r.role_class, r.role FROM models m CROSS JOIN required r",
          "LEFT JOIN %s a ON a.model_id=m.model_id",
          "AND a.role_class=r.role_class AND a.role=r.role AND a.required",
          "WHERE a.model_id IS NULL"
        ),
        values, relation_sql("model_manifest"), relation_sql("model_artifact")
      ),
      "model_id, role_class, role"
    )
  }

  if (relations_ready(validator, c("model_manifest"))) {
    for (table in setdiff(model_relations(validator$expected), "model_manifest")) {
      if (!isTRUE(validator$ready[[table]])) next
      check_query(
        validator, "DVMA2050", table,
        "model_id does not reference the sole pack manifest",
        sprintf(
          paste(
            "SELECT DISTINCT r.model_id FROM %s r LEFT JOIN %s m",
            "USING (model_id) WHERE m.model_id IS NULL"
          ),
          relation_sql(table), relation_sql("model_manifest")
        ),
        "model_id"
      )
    }
  }

  source_columns <- list(
    model_seq_region = "source_namespace",
    model_gene = "source_namespace",
    model_transcript = "source_namespace",
    model_exon = "source_namespace",
    model_transcript_exon = "source_namespace",
    model_translation = "source_namespace",
    model_attribute_type = "source_namespace",
    model_seq_region_attribute = "source_namespace",
    model_gene_attribute = "source_namespace",
    model_transcript_attribute = "source_namespace",
    model_translation_attribute = "source_namespace",
    model_transcript_edit = "source_namespace",
    model_translation_edit = "source_namespace",
    model_mature_mirna = "source_namespace",
    model_external_db = "source_namespace",
    model_xref = "source_namespace",
    model_seq_region_synonym = "source_namespace",
    model_sequence_state = "source_namespace"
  )
  if (validator$ready[["model_source"]]) {
    for (table in names(source_columns)) {
      if (!validator$ready[[table]]) next
      column <- source_columns[[table]]
      check_query(
        validator, "DVMA2051", table,
        "source namespace does not reference model_source",
        sprintf(
          paste(
            "SELECT DISTINCT r.model_id, r.%s AS source_namespace FROM %s r",
            "LEFT JOIN %s s ON s.model_id=r.model_id",
            "AND s.source_namespace=r.%s WHERE s.source_namespace IS NULL"
          ),
          sql_ident(column), relation_sql(table), relation_sql("model_source"),
          sql_ident(column)
        ),
        "model_id, source_namespace"
      )
    }
  }
}

check_nonempty_candidate <- function(validator) {
  if (validator$profile == "structural") return(invisible(NULL))
  required <- c(
    "model_source", "model_selection_audit", "model_build", "model_artifact",
    "model_seq_region", "model_gene", "model_transcript", "model_exon",
    "model_transcript_exon", "model_translation", "model_attribute_type",
    "model_transcript_attribute", "model_external_db", "model_xref",
    "model_object_xref", "model_xref_identity", "model_sequence_blob",
    "model_sequence_state"
  )
  for (table in required) {
    if (isTRUE(validator$ready[[table]]) && validator$row_counts[[table]] == 0) {
      add_issue(
        validator, "DVMA2060", table,
        "model candidate relation must be nonempty for consequence conformance"
      )
    }
  }
}

stable_version_predicate <- function(alias) {
  paste0(
    "(", alias, ".stable_id IS NULL AND ", alias,
    ".stable_id_version IS NOT NULL) OR (", alias,
    ".stable_id IS NOT NULL AND length(trim(", alias,
    ".stable_id)) = 0) OR (", alias,
    ".stable_id_version IS NOT NULL AND ", alias, ".stable_id_version < 1)"
  )
}

check_geometry <- function(validator) {
  if (validator$ready[["model_seq_region"]]) {
    check_query(
      validator, "DVMA3001", "model_seq_region",
      "sequence-region identity, geometry, or reference receipt is invalid",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, seq_region_key FROM %s WHERE",
          "seq_region_key = 0 OR length(trim(source_internal_id)) = 0",
          "OR length(trim(name)) = 0 OR length <= 0",
          "OR length(trim(coord_system_name)) = 0",
          "OR length(trim(assembly_accession)) = 0",
          "OR length(trim(reference_accession)) = 0",
          "OR length(trim(reference_version)) = 0",
          "OR NOT regexp_full_match(sequence_sha256, '[0-9a-f]{64}')"
        ),
        relation_sql("model_seq_region")
      ),
      "model_id, source_namespace, seq_region_key"
    )
    if (validator$ready[["model_manifest"]]) {
      check_query(
        validator, "DVMA3002", "model_seq_region",
        "sequence region assembly accession disagrees with the model manifest",
        sprintf(
          paste(
            "SELECT r.model_id,r.source_namespace,r.seq_region_key,r.assembly_accession",
            "FROM %s r JOIN %s m USING (model_id)",
            "WHERE r.assembly_accession<>m.assembly_accession"
          ),
          relation_sql("model_seq_region"), relation_sql("model_manifest")
        ),
        "model_id, source_namespace, seq_region_key"
      )
    }
  }

  if (relations_ready(validator, c("model_gene", "model_seq_region"))) {
    check_query(
      validator, "DVMA3010", "model_gene",
      "gene identity, stable-ID metadata, or geometry is invalid",
      sprintf(
        paste(
          "SELECT g.model_id, g.source_namespace, g.gene_key FROM %s g",
          "LEFT JOIN %s r USING (model_id, source_namespace, seq_region_key)",
          "WHERE g.gene_key = 0 OR length(trim(g.source_internal_id)) = 0",
          "OR %s OR g.start1 < 1 OR g.end1 < g.start1",
          "OR g.strand NOT IN (-1,1) OR length(trim(g.source)) = 0",
          "OR length(trim(g.biotype)) = 0 OR r.seq_region_key IS NULL",
          "OR g.end1 > r.length"
        ),
        relation_sql("model_gene"), relation_sql("model_seq_region"),
        stable_version_predicate("g")
      ),
      "model_id, source_namespace, gene_key"
    )
  }

  if (relations_ready(validator, c("model_transcript", "model_seq_region", "model_gene"))) {
    check_query(
      validator, "DVMA3020", "model_transcript",
      "transcript identity, owner, stable-ID metadata, or geometry is invalid",
      sprintf(
        paste(
          "SELECT t.model_id, t.source_namespace, t.transcript_key FROM %s t",
          "LEFT JOIN %s r USING (model_id, source_namespace, seq_region_key)",
          "LEFT JOIN %s g ON g.model_id=t.model_id",
          "AND g.source_namespace=t.source_namespace AND g.gene_key=t.gene_key",
          "WHERE t.transcript_key = 0 OR length(trim(t.source_internal_id)) = 0",
          "OR %s OR t.start1 < 1 OR t.end1 < t.start1",
          "OR t.strand NOT IN (-1,1) OR length(trim(t.source)) = 0",
          "OR length(trim(t.biotype)) = 0 OR r.seq_region_key IS NULL",
          "OR t.end1 > r.length OR (t.gene_key IS NOT NULL AND g.gene_key IS NULL)",
          "OR (g.gene_key IS NOT NULL AND (g.seq_region_key <> t.seq_region_key",
          "OR g.strand <> t.strand OR t.start1 < g.start1 OR t.end1 > g.end1))"
        ),
        relation_sql("model_transcript"), relation_sql("model_seq_region"),
        relation_sql("model_gene"), stable_version_predicate("t")
      ),
      "model_id, source_namespace, transcript_key"
    )
    if (validator$ready[["model_selection_audit"]]) {
      check_query(
        validator, "DVMA3021", "model_transcript",
        "a gene-orphan transcript requires an explicit source-policy audit count",
        sprintf(
          paste(
            "WITH orphans AS (SELECT model_id,source_namespace,count(*) AS n FROM %s",
            "WHERE gene_key IS NULL GROUP BY model_id,source_namespace),",
            "admitted AS (SELECT model_id,source_namespace,sum(row_count) AS n FROM %s",
            "WHERE object_kind='transcript' AND decision='included'",
            "AND reason_code='orphan_transcript_admitted'",
            "GROUP BY model_id,source_namespace)",
            "SELECT o.model_id,o.source_namespace,o.n AS orphan_count,",
            "coalesce(a.n,0) AS audited_count FROM orphans o LEFT JOIN admitted a",
            "USING (model_id,source_namespace) WHERE coalesce(a.n,0)<>o.n"
          ),
          relation_sql("model_transcript"), relation_sql("model_selection_audit")
        ),
        "model_id, source_namespace"
      )
    }
    check_query(
      validator, "DVMA3022", "model_transcript",
      "gene-orphan transcript cannot be marked as its gene's canonical transcript",
      sprintf(
        "SELECT model_id,source_namespace,transcript_key FROM %s WHERE gene_key IS NULL AND is_canonical",
        relation_sql("model_transcript")
      ),
      "model_id, source_namespace, transcript_key"
    )
  }

  if (relations_ready(validator, c("model_exon", "model_seq_region"))) {
    check_query(
      validator, "DVMA3030", "model_exon",
      "exon identity, stable-ID metadata, phase, or geometry is invalid",
      sprintf(
        paste(
          "SELECT e.model_id, e.source_namespace, e.exon_key FROM %s e",
          "LEFT JOIN %s r USING (model_id, source_namespace, seq_region_key)",
          "WHERE e.exon_key = 0 OR length(trim(e.source_internal_id)) = 0",
          "OR %s OR e.start1 < 1 OR e.end1 < e.start1",
          "OR e.strand NOT IN (-1,1) OR e.phase NOT IN (-1,0,1,2)",
          "OR e.end_phase NOT IN (-1,0,1,2) OR r.seq_region_key IS NULL",
          "OR e.end1 > r.length"
        ),
        relation_sql("model_exon"), relation_sql("model_seq_region"),
        stable_version_predicate("e")
      ),
      "model_id, source_namespace, exon_key"
    )
  }
}

check_membership_and_translation <- function(validator) {
  if (relations_ready(validator, c(
      "model_transcript_exon", "model_transcript", "model_exon"
  ))) {
    check_query(
      validator, "DVMA3100", "model_transcript_exon",
      "transcript/exon membership has a missing owner or inconsistent geometry",
      sprintf(
        paste(
          "SELECT m.model_id, m.source_namespace, m.transcript_key, m.exon_key,",
          "m.exon_rank FROM %s m",
          "LEFT JOIN %s t USING (model_id, source_namespace, transcript_key)",
          "LEFT JOIN %s e USING (model_id, source_namespace, exon_key)",
          "WHERE t.transcript_key IS NULL OR e.exon_key IS NULL",
          "OR m.exon_rank < 1 OR m.raw_cdna_start1 < 1",
          "OR m.raw_cdna_end1 < m.raw_cdna_start1",
          "OR m.raw_cdna_end1-m.raw_cdna_start1 <> e.end1-e.start1",
          "OR e.seq_region_key <> t.seq_region_key OR e.strand <> t.strand",
          "OR e.start1 < t.start1 OR e.end1 > t.end1"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_transcript"),
        relation_sql("model_exon")
      ),
      "model_id, source_namespace, transcript_key, exon_rank, exon_key"
    )
    check_query(
      validator, "DVMA3101", "model_transcript_exon",
      "exon ranks must be contiguous from one and one membership must occupy each rank",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, transcript_key, min(exon_rank) AS min_rank,",
          "max(exon_rank) AS max_rank, count(*) AS membership_count,",
          "count(DISTINCT exon_rank) AS rank_count FROM %s",
          "GROUP BY model_id, source_namespace, transcript_key",
          "HAVING min(exon_rank) <> 1 OR max(exon_rank) <> count(*)",
          "OR count(DISTINCT exon_rank) <> count(*)"
        ),
        relation_sql("model_transcript_exon")
      ),
      "model_id, source_namespace, transcript_key"
    )
    check_query(
      validator, "DVMA3102", "model_transcript_exon",
      "repeated exon membership must be quarantined before model-candidate validation",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, transcript_key, exon_key,",
          "count(*) AS membership_count FROM %s",
          "GROUP BY model_id, source_namespace, transcript_key, exon_key",
          "HAVING count(*) > 1"
        ),
        relation_sql("model_transcript_exon")
      ),
      "model_id, source_namespace, transcript_key, exon_key"
    )
    check_query(
      validator, "DVMA3103", "model_transcript_exon",
      "exons must be non-overlapping in transcript order with contiguous raw cDNA coordinates",
      sprintf(
        paste(
          "WITH ordered AS (SELECT m.*, t.strand, e.start1, e.end1,",
          "lag(e.start1) OVER w AS previous_start1,",
          "lag(e.end1) OVER w AS previous_end1,",
          "lag(m.raw_cdna_end1) OVER w AS previous_cdna_end1",
          "FROM %s m JOIN %s t USING (model_id, source_namespace, transcript_key)",
          "JOIN %s e USING (model_id, source_namespace, exon_key)",
          "WINDOW w AS (PARTITION BY m.model_id, m.source_namespace, m.transcript_key",
          "ORDER BY m.exon_rank))",
          "SELECT model_id, source_namespace, transcript_key, exon_rank FROM ordered",
          "WHERE (exon_rank=1 AND raw_cdna_start1<>1)",
          "OR (exon_rank>1 AND raw_cdna_start1-1<>previous_cdna_end1)",
          "OR (exon_rank>1 AND strand=1 AND start1<=previous_end1)",
          "OR (exon_rank>1 AND strand=-1 AND end1>=previous_start1)"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_transcript"),
        relation_sql("model_exon")
      ),
      "model_id, source_namespace, transcript_key, exon_rank"
    )
    check_query(
      validator, "DVMA3104", "model_transcript",
      "every retained transcript must have at least one exon membership",
      sprintf(
        paste(
          "SELECT t.model_id, t.source_namespace, t.transcript_key FROM %s t",
          "LEFT JOIN %s m USING (model_id, source_namespace, transcript_key)",
          "WHERE m.transcript_key IS NULL"
        ),
        relation_sql("model_transcript"), relation_sql("model_transcript_exon")
      ),
      "model_id, source_namespace, transcript_key"
    )
    check_query(
      validator, "DVMA3105", "model_transcript",
      "transcript source bounds must equal the genomic envelope of retained exon memberships",
      sprintf(
        paste(
          "WITH envelope AS (SELECT m.model_id,m.source_namespace,m.transcript_key,",
          "min(e.start1) AS exon_start1,max(e.end1) AS exon_end1 FROM %s m JOIN %s e",
          "USING (model_id,source_namespace,exon_key)",
          "GROUP BY m.model_id,m.source_namespace,m.transcript_key)",
          "SELECT t.model_id,t.source_namespace,t.transcript_key,t.start1,t.end1,",
          "e.exon_start1,e.exon_end1 FROM %s t JOIN envelope e",
          "USING (model_id,source_namespace,transcript_key)",
          "WHERE t.start1<>e.exon_start1 OR t.end1<>e.exon_end1"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_exon"),
        relation_sql("model_transcript")
      ),
      "model_id, source_namespace, transcript_key"
    )
  }

  if (relations_ready(validator, c(
      "model_translation", "model_transcript", "model_transcript_exon", "model_exon"
  ))) {
    check_query(
      validator, "DVMA3200", "model_translation",
      "translation identity, stable-ID metadata, coding bounds, phase, or codon table is invalid",
      sprintf(
        paste(
          "SELECT x.model_id, x.source_namespace, x.translation_key FROM %s x",
          "LEFT JOIN %s t USING (model_id, source_namespace, transcript_key)",
          "WHERE x.translation_key=0 OR length(trim(x.source_internal_id))=0",
          "OR %s OR t.transcript_key IS NULL",
          "OR x.start_offset1<1 OR x.end_offset1<1",
          "OR x.raw_cdna_coding_start1<1",
          "OR x.raw_cdna_coding_end1<x.raw_cdna_coding_start1",
          "OR x.edited_cdna_coding_start1<1",
          "OR x.edited_cdna_coding_end1<x.edited_cdna_coding_start1",
          "OR x.start_phase_padding NOT BETWEEN 0 AND 2",
          "OR x.codon_table NOT IN (1,2)"
        ),
        relation_sql("model_translation"), relation_sql("model_transcript"),
        stable_version_predicate("x")
      ),
      "model_id, source_namespace, translation_key"
    )
    check_query(
      validator, "DVMA3201", "model_translation",
      "translation boundary exon or exon-relative offset is invalid",
      sprintf(
        paste(
          "SELECT x.model_id, x.source_namespace, x.translation_key FROM %s x",
          "LEFT JOIN %s sm ON sm.model_id=x.model_id",
          "AND sm.source_namespace=x.source_namespace",
          "AND sm.transcript_key=x.transcript_key AND sm.exon_key=x.start_exon_key",
          "LEFT JOIN %s se ON se.model_id=x.model_id",
          "AND se.source_namespace=x.source_namespace AND se.exon_key=x.start_exon_key",
          "LEFT JOIN %s em ON em.model_id=x.model_id",
          "AND em.source_namespace=x.source_namespace",
          "AND em.transcript_key=x.transcript_key AND em.exon_key=x.end_exon_key",
          "LEFT JOIN %s ee ON ee.model_id=x.model_id",
          "AND ee.source_namespace=x.source_namespace AND ee.exon_key=x.end_exon_key",
          "WHERE sm.exon_key IS NULL OR em.exon_key IS NULL",
          "OR x.start_offset1>se.end1-se.start1+1",
          "OR x.end_offset1>ee.end1-ee.start1+1 OR sm.exon_rank>em.exon_rank",
          "OR x.raw_cdna_coding_start1-sm.raw_cdna_start1<>x.start_offset1-1",
          "OR x.raw_cdna_coding_end1-em.raw_cdna_start1<>x.end_offset1-1"
        ),
        relation_sql("model_translation"), relation_sql("model_transcript_exon"),
        relation_sql("model_exon"), relation_sql("model_transcript_exon"),
        relation_sql("model_exon")
      ),
      "model_id, source_namespace, translation_key"
    )
    check_query(
      validator, "DVMA3202", "model_translation",
      "raw cDNA coding bounds exceed the transcript's raw spliced length",
      sprintf(
        paste(
          "WITH lengths AS (SELECT model_id, source_namespace, transcript_key,",
          "max(raw_cdna_end1) AS raw_length FROM %s",
          "GROUP BY model_id, source_namespace, transcript_key)",
          "SELECT x.model_id, x.source_namespace, x.translation_key FROM %s x",
          "JOIN lengths l USING (model_id, source_namespace, transcript_key)",
          "WHERE x.raw_cdna_coding_end1>l.raw_length"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_translation")
      ),
      "model_id, source_namespace, translation_key"
    )
  }

  if (relations_ready(validator, c("model_gene", "model_transcript"))) {
    check_query(
      validator, "DVMA3210", "model_gene",
      "canonical transcript pointer and transcript is_canonical flag disagree",
      sprintf(
        paste(
          "WITH canonical AS (SELECT model_id, source_namespace, gene_key,",
          "min(transcript_key) AS transcript_key, count(*) AS n FROM %s",
          "WHERE is_canonical GROUP BY model_id, source_namespace, gene_key)",
          "SELECT g.model_id, g.source_namespace, g.gene_key,",
          "g.canonical_transcript_key, c.transcript_key AS flagged_key, c.n",
          "FROM %s g LEFT JOIN canonical c USING (model_id, source_namespace, gene_key)",
          "WHERE (g.canonical_transcript_key IS NULL) <> (c.transcript_key IS NULL)",
          "OR g.canonical_transcript_key IS DISTINCT FROM c.transcript_key",
          "OR coalesce(c.n,0)>1"
        ),
        relation_sql("model_transcript"), relation_sql("model_gene")
      ),
      "model_id, source_namespace, gene_key"
    )
  }

  if (relations_ready(validator, c("model_transcript", "model_translation"))) {
    check_query(
      validator, "DVMA3211", "model_transcript",
      "canonical translation pointer and translation is_canonical flag disagree",
      sprintf(
        paste(
          "WITH canonical AS (SELECT model_id, source_namespace, transcript_key,",
          "min(translation_key) AS translation_key, count(*) AS n FROM %s",
          "WHERE is_canonical GROUP BY model_id, source_namespace, transcript_key)",
          "SELECT t.model_id, t.source_namespace, t.transcript_key,",
          "t.canonical_translation_key, c.translation_key AS flagged_key, c.n",
          "FROM %s t LEFT JOIN canonical c",
          "USING (model_id, source_namespace, transcript_key)",
          "WHERE (t.canonical_translation_key IS NULL) <> (c.translation_key IS NULL)",
          "OR t.canonical_translation_key IS DISTINCT FROM c.translation_key",
          "OR coalesce(c.n,0)>1"
        ),
        relation_sql("model_translation"), relation_sql("model_transcript")
      ),
      "model_id, source_namespace, transcript_key"
    )
  }
}

ATTRIBUTE_RELATIONS <- list(
  model_seq_region_attribute = list(owner_table = "model_seq_region", owner_key = "seq_region_key"),
  model_gene_attribute = list(owner_table = "model_gene", owner_key = "gene_key"),
  model_transcript_attribute = list(owner_table = "model_transcript", owner_key = "transcript_key"),
  model_translation_attribute = list(owner_table = "model_translation", owner_key = "translation_key")
)

check_attributes <- function(validator) {
  if (validator$ready[["model_attribute_type"]]) {
    check_query(
      validator, "DVMA3300", "model_attribute_type",
      "attribute-type source identity or named definition is invalid",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, attrib_type_key FROM %s WHERE",
          "attrib_type_key=0 OR length(trim(source_internal_id))=0",
          "OR length(trim(code))=0 OR length(trim(name))=0"
        ),
        relation_sql("model_attribute_type")
      ),
      "model_id, source_namespace, attrib_type_key"
    )
  }

  for (table in names(ATTRIBUTE_RELATIONS)) {
    spec <- ATTRIBUTE_RELATIONS[[table]]
    requirements <- c(table, spec$owner_table, "model_attribute_type")
    if (!relations_ready(validator, requirements)) next
    owner_key <- sql_ident(spec$owner_key)
    check_query(
      validator, "DVMA3310", table,
      "raw attribute has an invalid owner/type reference or source locator",
      sprintf(
        paste(
          "SELECT a.model_id, a.source_namespace, a.attribute_key,",
          "a.%s AS owner_key FROM %s a",
          "LEFT JOIN %s o ON o.model_id=a.model_id",
          "AND o.source_namespace=a.source_namespace AND o.%s=a.%s",
          "LEFT JOIN %s t ON t.model_id=a.model_id",
          "AND t.source_namespace=a.source_namespace",
          "AND t.attrib_type_key=a.attrib_type_key",
          "WHERE a.attribute_key=0 OR a.duplicate_ordinal<1",
          "OR length(trim(a.source_row_locator))=0",
          "OR o.%s IS NULL OR t.attrib_type_key IS NULL"
        ),
        owner_key, relation_sql(table), relation_sql(spec$owner_table),
        owner_key, owner_key, relation_sql("model_attribute_type"), owner_key
      ),
      "model_id, source_namespace, attribute_key"
    )
    check_query(
      validator, "DVMA3311", table,
      "duplicate ordinals must be contiguous from one for each complete raw tuple",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, %s AS owner_key, attrib_type_key, value,",
          "min(duplicate_ordinal) AS min_ordinal,",
          "max(duplicate_ordinal) AS max_ordinal, count(*) AS row_count",
          "FROM %s GROUP BY model_id, source_namespace, %s, attrib_type_key, value",
          "HAVING min(duplicate_ordinal)<>1 OR max(duplicate_ordinal)<>count(*)",
          "OR count(DISTINCT duplicate_ordinal)<>count(*)"
        ),
        owner_key, relation_sql(table), owner_key
      ),
      "model_id, source_namespace, owner_key, attrib_type_key, value"
    )
    check_query(
      validator, "DVMA3312", table,
      "pack-local attribute_key is not unique within its source namespace",
      sprintf(
        paste(
          "SELECT model_id,source_namespace,attribute_key,count(*) AS row_count FROM %s",
          "GROUP BY model_id,source_namespace,attribute_key HAVING count(*)>1"
        ),
        relation_sql(table)
      ),
      "model_id, source_namespace, attribute_key"
    )
  }

  if (relations_ready(validator, c(
      "model_transcript_attribute", "model_attribute_type", "model_selection_audit",
      "model_source"
  ))) {
    check_query(
      validator, "DVMA3320", "model_transcript_attribute",
      "translation-boundary override is malformed or conflicting for one transcript",
      sprintf(
        paste(
          "SELECT a.model_id,a.source_namespace,a.transcript_key,t.code,",
          "count(DISTINCT a.value) AS distinct_values,",
          "min(try_cast(a.value AS BIGINT)) AS parsed_value FROM %s a JOIN %s t",
          "USING (model_id,source_namespace,attrib_type_key)",
          "WHERE t.code IN ('_transl_start','_transl_end')",
          "GROUP BY a.model_id,a.source_namespace,a.transcript_key,t.code",
          "HAVING count(DISTINCT a.value)<>1 OR min(try_cast(a.value AS BIGINT))<1",
          "OR min(try_cast(a.value AS BIGINT)) IS NULL"
        ),
        relation_sql("model_transcript_attribute"), relation_sql("model_attribute_type")
      ),
      "model_id, source_namespace, transcript_key, code"
    )
    if (validator$profile != "structural") {
      check_query(
        validator, "DVMA3321", "model_selection_audit",
        "model must preserve translation-boundary overrides or carry a zero-row absence proof",
        sprintf(
          paste(
            "WITH override_rows AS (SELECT a.model_id,a.source_namespace,count(*) AS n",
            "FROM %s a JOIN %s t",
            "USING (model_id,source_namespace,attrib_type_key)",
            "WHERE t.code IN ('_transl_start','_transl_end')",
            "GROUP BY a.model_id,a.source_namespace),",
            "proof AS (SELECT model_id,source_namespace,count(*) AS n FROM %s",
            "WHERE object_kind='transcript_attribute'",
            "AND reason_code='translation_boundary_overrides_absent' AND row_count=0",
            "GROUP BY model_id,source_namespace)",
            "SELECT s.model_id,s.source_namespace,coalesce(o.n,0) AS override_rows,",
            "coalesce(p.n,0) AS proof_rows FROM %s s LEFT JOIN override_rows o",
            "USING (model_id,source_namespace) LEFT JOIN proof p",
            "USING (model_id,source_namespace)",
            "WHERE coalesce(o.n,0)=0 AND coalesce(p.n,0)=0"
          ),
          relation_sql("model_transcript_attribute"), relation_sql("model_attribute_type"),
          relation_sql("model_selection_audit"), relation_sql("model_source")
        ),
        "model_id, source_namespace"
      )
    }
  }

  if (relations_ready(validator, c(
      "model_transcript_attribute", "model_attribute_type", "model_translation"
  ))) {
    check_query(
      validator, "DVMA3322", "model_translation",
      "materialized edited cDNA coding bounds disagree with transcript override attributes",
      sprintf(
        paste(
          "WITH overrides AS (SELECT a.model_id,a.source_namespace,a.transcript_key,",
          "max(try_cast(a.value AS BIGINT)) FILTER (WHERE t.code='_transl_start') AS start1,",
          "max(try_cast(a.value AS BIGINT)) FILTER (WHERE t.code='_transl_end') AS end1",
          "FROM %s a JOIN %s t USING (model_id,source_namespace,attrib_type_key)",
          "WHERE t.code IN ('_transl_start','_transl_end')",
          "GROUP BY a.model_id,a.source_namespace,a.transcript_key)",
          "SELECT x.model_id,x.source_namespace,x.translation_key,o.start1,o.end1,",
          "x.edited_cdna_coding_start1,x.edited_cdna_coding_end1 FROM %s x JOIN overrides o",
          "USING (model_id,source_namespace,transcript_key)",
          "WHERE (o.start1 IS NOT NULL AND o.start1<>x.edited_cdna_coding_start1)",
          "OR (o.end1 IS NOT NULL AND o.end1<>x.edited_cdna_coding_end1)"
        ),
        relation_sql("model_transcript_attribute"), relation_sql("model_attribute_type"),
        relation_sql("model_translation")
      ),
      "model_id, source_namespace, translation_key"
    )
  }

  if (relations_ready(validator, c(
      "model_transcript_attribute", "model_attribute_type",
      "model_translation", "model_transcript_edit"
  ))) {
    check_query(
      validator, "DVMA3323", "model_translation",
      "edited coding bounds disagree with the raw-to-edited nucleotide-edit coordinate map",
      sprintf(
        paste(
          "WITH overrides AS (SELECT a.model_id,a.source_namespace,a.transcript_key,",
          "max(try_cast(a.value AS BIGINT)) FILTER (WHERE t.code='_transl_start') AS start1,",
          "max(try_cast(a.value AS BIGINT)) FILTER (WHERE t.code='_transl_end') AS end1",
          "FROM %s a JOIN %s t USING (model_id,source_namespace,attrib_type_key)",
          "WHERE t.code IN ('_transl_start','_transl_end')",
          "GROUP BY a.model_id,a.source_namespace,a.transcript_key),",
          "mapped AS (SELECT x.model_id,x.source_namespace,x.translation_key,x.transcript_key,",
          "x.raw_cdna_coding_start1,x.raw_cdna_coding_end1,",
          "x.edited_cdna_coding_start1,x.edited_cdna_coding_end1,",
          "o.start1 AS override_start1,o.end1 AS override_end1,",
          "coalesce(sum(CASE WHEN e.status='applied' AND",
          "((e.start1-1=e.end1 AND e.start1<=x.raw_cdna_coding_start1) OR",
          "(e.start1<=e.end1 AND e.end1<x.raw_cdna_coding_start1))",
          "THEN octet_length(e.alt_seq)-octet_length(e.preapply_ref_seq) ELSE 0 END),0)",
          "AS start_delta,",
          "coalesce(sum(CASE WHEN e.status='applied' AND",
          "((e.start1-1=e.end1 AND e.start1<=x.raw_cdna_coding_end1) OR",
          "(e.start1<=e.end1 AND e.end1<x.raw_cdna_coding_end1))",
          "THEN octet_length(e.alt_seq)-octet_length(e.preapply_ref_seq) ELSE 0 END),0)",
          "AS end_delta,",
          "coalesce(bool_or(e.status='applied' AND e.start1<=e.end1",
          "AND octet_length(e.alt_seq)<>octet_length(e.preapply_ref_seq)",
          "AND e.start1<=x.raw_cdna_coding_start1",
          "AND e.end1>=x.raw_cdna_coding_start1),false) AS start_overlap,",
          "coalesce(bool_or(e.status='applied' AND e.start1<=e.end1",
          "AND octet_length(e.alt_seq)<>octet_length(e.preapply_ref_seq)",
          "AND e.start1<=x.raw_cdna_coding_end1",
          "AND e.end1>=x.raw_cdna_coding_end1),false) AS end_overlap",
          "FROM %s x LEFT JOIN %s e ON e.model_id=x.model_id",
          "AND e.source_namespace=x.source_namespace AND e.transcript_key=x.transcript_key",
          "LEFT JOIN overrides o ON o.model_id=x.model_id",
          "AND o.source_namespace=x.source_namespace AND o.transcript_key=x.transcript_key",
          "GROUP BY x.model_id,x.source_namespace,x.translation_key,x.transcript_key,",
          "x.raw_cdna_coding_start1,x.raw_cdna_coding_end1,",
          "x.edited_cdna_coding_start1,x.edited_cdna_coding_end1,o.start1,o.end1)",
          "SELECT model_id,source_namespace,translation_key,",
          "raw_cdna_coding_start1,raw_cdna_coding_end1,",
          "edited_cdna_coding_start1,edited_cdna_coding_end1,",
          "override_start1,override_end1,start_delta,end_delta,start_overlap,end_overlap",
          "FROM mapped WHERE (start_overlap AND override_start1 IS NULL)",
          "OR (end_overlap AND override_end1 IS NULL)",
          "OR edited_cdna_coding_start1<>",
          "coalesce(override_start1,raw_cdna_coding_start1+start_delta)",
          "OR edited_cdna_coding_end1<>",
          "coalesce(override_end1,raw_cdna_coding_end1+end_delta)"
        ),
        relation_sql("model_transcript_attribute"), relation_sql("model_attribute_type"),
        relation_sql("model_translation"), relation_sql("model_transcript_edit")
      ),
      "model_id, source_namespace, translation_key"
    )
  }

  if (relations_ready(validator, c(
      "model_seq_region_attribute", "model_attribute_type"
  ))) {
    check_query(
      validator, "DVMA3330", "model_seq_region_attribute",
      "codon_table source attribute is malformed or conflicting",
      sprintf(
        paste(
          "SELECT a.model_id,a.source_namespace,a.seq_region_key,",
          "count(DISTINCT a.value) AS distinct_values,",
          "min(try_cast(a.value AS INTEGER)) AS codon_table FROM %s a JOIN %s t",
          "USING (model_id,source_namespace,attrib_type_key) WHERE t.code='codon_table'",
          "GROUP BY a.model_id,a.source_namespace,a.seq_region_key",
          "HAVING count(DISTINCT a.value)<>1",
          "OR min(try_cast(a.value AS INTEGER)) NOT IN (1,2)",
          "OR min(try_cast(a.value AS INTEGER)) IS NULL"
        ),
        relation_sql("model_seq_region_attribute"), relation_sql("model_attribute_type")
      ),
      "model_id, source_namespace, seq_region_key"
    )
  }

  if (relations_ready(validator, c(
      "model_seq_region_attribute", "model_attribute_type",
      "model_translation", "model_transcript"
  ))) {
    check_query(
      validator, "DVMA3331", "model_translation",
      "translation codon table disagrees with seq_region_attrib or its explicit table-1 default",
      sprintf(
        paste(
          "WITH codon AS (SELECT a.model_id,a.source_namespace,a.seq_region_key,",
          "max(try_cast(a.value AS INTEGER)) AS codon_table FROM %s a JOIN %s ty",
          "USING (model_id,source_namespace,attrib_type_key) WHERE ty.code='codon_table'",
          "GROUP BY a.model_id,a.source_namespace,a.seq_region_key)",
          "SELECT x.model_id,x.source_namespace,x.translation_key,x.codon_table,",
          "coalesce(c.codon_table,1) AS expected_codon_table FROM %s x JOIN %s t",
          "USING (model_id,source_namespace,transcript_key) LEFT JOIN codon c",
          "ON c.model_id=t.model_id AND c.source_namespace=t.source_namespace",
          "AND c.seq_region_key=t.seq_region_key",
          "WHERE x.codon_table<>coalesce(c.codon_table,1)"
        ),
        relation_sql("model_seq_region_attribute"), relation_sql("model_attribute_type"),
        relation_sql("model_translation"), relation_sql("model_transcript")
      ),
      "model_id, source_namespace, translation_key"
    )
  }
}

check_edit_relation <- function(validator, table, owner_table, owner_key,
                                attribute_table, expected_basis, allowed_applied_codes) {
  requirements <- c(table, owner_table, attribute_table, "model_attribute_type")
  if (!relations_ready(validator, requirements)) return(invisible(NULL))
  owner_key_sql <- sql_ident(owner_key)
  code_values <- paste(vapply(allowed_applied_codes, sql_literal, character(1)), collapse = ",")

  check_query(
    validator, "DVMA3400", table,
    "typed edit does not resolve to its owner, raw attribute, and pinned attribute code",
    sprintf(
      paste(
        "SELECT e.model_id, e.source_namespace, e.edit_key, e.%s AS owner_key",
        "FROM %s e LEFT JOIN %s o ON o.model_id=e.model_id",
        "AND o.source_namespace=e.source_namespace AND o.%s=e.%s",
        "LEFT JOIN %s a ON a.model_id=e.model_id",
        "AND a.source_namespace=e.source_namespace AND a.attribute_key=e.attribute_key",
        "LEFT JOIN %s t ON t.model_id=a.model_id",
        "AND t.source_namespace=a.source_namespace AND t.attrib_type_key=a.attrib_type_key",
        "WHERE e.edit_key=0 OR o.%s IS NULL OR a.attribute_key IS NULL",
        "OR a.%s IS DISTINCT FROM e.%s OR t.code IS DISTINCT FROM e.code"
      ),
      owner_key_sql, relation_sql(table), relation_sql(owner_table),
      owner_key_sql, owner_key_sql, relation_sql(attribute_table),
      relation_sql("model_attribute_type"), owner_key_sql, owner_key_sql, owner_key_sql
    ),
    "model_id, source_namespace, edit_key"
  )
  check_query(
    validator, "DVMA3401", table,
    "typed edit basis, status, applied-code allowlist, or ordinal coupling is invalid",
    sprintf(
      paste(
        "SELECT model_id, source_namespace, edit_key, status, apply_ordinal, code",
        "FROM %s WHERE coordinate_basis<>%s",
        "OR status NOT IN ('applied','source_ignored',",
        "'quarantined_ambiguous_order','unsupported_code')",
        "OR ((status='applied') <> (apply_ordinal IS NOT NULL))",
        "OR (apply_ordinal IS NOT NULL AND apply_ordinal<1)",
        "OR (status='applied' AND code NOT IN (%s))",
        "OR (status='unsupported_code' AND code IN (%s))"
      ),
      relation_sql(table), sql_literal(expected_basis), code_values, code_values
    ),
    "model_id, source_namespace, edit_key"
  )
  check_query(
    validator, "DVMA3402", table,
    "typed edit coordinates or original-basis reference snapshot are invalid",
    sprintf(
      paste(
        "SELECT model_id, source_namespace, edit_key, start1, end1,",
        "octet_length(basis_ref_seq) AS basis_ref_length",
        "FROM %s WHERE start1<1 OR end1<0 OR start1-1>end1",
        "OR octet_length(basis_ref_seq)<>greatest(end1-start1+1,0)"
      ),
      relation_sql(table)
    ),
    "model_id, source_namespace, edit_key"
  )
  check_query(
    validator, "DVMA3403", table,
    "applied edit ordinals must be contiguous and follow descending start coordinates",
    sprintf(
      paste(
        "WITH ordered AS (SELECT *,",
        "lag(start1) OVER (PARTITION BY model_id, source_namespace, %s",
        "ORDER BY apply_ordinal) AS previous_start,",
        "min(apply_ordinal) OVER (PARTITION BY model_id, source_namespace, %s) AS min_ordinal,",
        "max(apply_ordinal) OVER (PARTITION BY model_id, source_namespace, %s) AS max_ordinal,",
        "count(*) OVER (PARTITION BY model_id, source_namespace, %s) AS edit_count",
        "FROM %s WHERE status='applied')",
        "SELECT model_id, source_namespace, edit_key, apply_ordinal, start1",
        "FROM ordered WHERE min_ordinal<>1 OR max_ordinal<>edit_count",
        "OR start1>previous_start"
      ),
      owner_key_sql, owner_key_sql, owner_key_sql, owner_key_sql, relation_sql(table)
    ),
    "model_id, source_namespace, owner_key, apply_ordinal"
  )
}

raw_identical <- function(left, right) {
  identical(as.raw(left), as.raw(right))
}

blob_value <- function(value) {
  if (is.raw(value)) return(value)
  if (is.list(value) && length(value) == 1L && is.raw(value[[1L]])) return(value[[1L]])
  if (inherits(value, "blob") && length(value) == 1L) return(value[[1L]])
  as.raw(value)
}

raw_slice <- function(value, start1, width) {
  if (width <= 0L || start1 > length(value)) return(raw())
  last <- min(length(value), start1 + width - 1L)
  value[seq.int(start1, last)]
}

raw_replace <- function(value, start1, width, replacement) {
  before <- if (start1 <= 1L) raw() else value[seq_len(min(length(value), start1 - 1L))]
  after_start <- start1 + width
  after <- if (after_start > length(value)) raw() else value[seq.int(after_start, length(value))]
  c(before, replacement, after)
}

check_edit_replay_one <- function(validator, table, owner_key, basis_role, result_role) {
  if (!relations_ready(validator, c(table, "model_sequence_state", "model_sequence_blob"))) {
    return(invisible(NULL))
  }
  owner_kind <- if (owner_key == "transcript_key") "transcript" else "translation"
  query <- sprintf(
    paste(
      "SELECT e.model_id, e.source_namespace, e.%s AS owner_key, e.edit_key,",
      "e.start1, e.end1, e.apply_ordinal, e.basis_ref_seq,",
      "e.preapply_ref_seq, e.alt_seq,",
      "bb.sequence_bytes AS basis_bytes, rb.sequence_bytes AS result_bytes",
      "FROM %s e",
      "LEFT JOIN %s bs ON bs.model_id=e.model_id",
      "AND bs.source_namespace=e.source_namespace AND bs.owner_kind=%s",
      "AND bs.owner_key=e.%s AND bs.role=%s AND bs.status='present'",
      "LEFT JOIN %s bb ON bb.sequence_sha256=bs.sequence_sha256",
      "AND bb.alphabet=bs.alphabet",
      "LEFT JOIN %s rs ON rs.model_id=e.model_id",
      "AND rs.source_namespace=e.source_namespace AND rs.owner_kind=%s",
      "AND rs.owner_key=e.%s AND rs.role=%s AND rs.status='present'",
      "LEFT JOIN %s rb ON rb.sequence_sha256=rs.sequence_sha256",
      "AND rb.alphabet=rs.alphabet",
      "WHERE e.status='applied' ORDER BY e.model_id, e.source_namespace,",
      "e.%s, e.apply_ordinal, e.edit_key"
    ),
    sql_ident(owner_key), relation_sql(table), relation_sql("model_sequence_state"),
    sql_literal(owner_kind), sql_ident(owner_key), sql_literal(basis_role),
    relation_sql("model_sequence_blob"), relation_sql("model_sequence_state"),
    sql_literal(owner_kind), sql_ident(owner_key), sql_literal(result_role),
    relation_sql("model_sequence_blob"), sql_ident(owner_key)
  )
  rows <- DBI::dbGetQuery(validator$connection, query)
  if (!nrow(rows)) return(invisible(NULL))

  group_id <- paste(rows$model_id, rows$source_namespace, as.character(rows$owner_key), sep = "\r")
  failures <- list()
  for (indices in split(seq_len(nrow(rows)), group_id)) {
    first <- indices[[1L]]
    if (is.na(rows$basis_bytes[first]) || is.na(rows$result_bytes[first])) next
    basis <- blob_value(rows$basis_bytes[first])
    current <- basis
    replay_valid <- TRUE
    for (index in indices) {
      start1 <- as.numeric(rows$start1[[index]])
      end1 <- as.numeric(rows$end1[[index]])
      width <- max(end1 - start1 + 1, 0)
      representable <- all(is.finite(c(start1, end1, width))) &&
        start1 == floor(start1) && end1 == floor(end1) && width == floor(width) &&
        start1 >= 1 && end1 >= 0 && start1 - 1 <= end1 &&
        start1 <= length(basis) + 1
      if (!representable || start1 > length(current) + 1) {
        failures[[length(failures) + 1L]] <- list(
          model_id = rows$model_id[[index]], source_namespace = rows$source_namespace[[index]],
          owner_key = as.character(rows$owner_key[[index]]),
          edit_key = as.character(rows$edit_key[[index]]),
          failure = "replay_coordinate_unrepresentable"
        )
        replay_valid <- FALSE
        break
      }
      expected_basis <- raw_slice(basis, start1, width)
      actual_preapply <- raw_slice(current, start1, width)
      if (!raw_identical(blob_value(rows$basis_ref_seq[index]), expected_basis)) {
        failures[[length(failures) + 1L]] <- list(
          model_id = rows$model_id[[index]], source_namespace = rows$source_namespace[[index]],
          owner_key = as.character(rows$owner_key[[index]]),
          edit_key = as.character(rows$edit_key[[index]]), failure = "basis_ref_mismatch"
        )
      }
      if (!raw_identical(blob_value(rows$preapply_ref_seq[index]), actual_preapply)) {
        failures[[length(failures) + 1L]] <- list(
          model_id = rows$model_id[[index]], source_namespace = rows$source_namespace[[index]],
          owner_key = as.character(rows$owner_key[[index]]),
          edit_key = as.character(rows$edit_key[[index]]), failure = "preapply_ref_mismatch"
        )
      }
      current <- raw_replace(current, start1, width, blob_value(rows$alt_seq[index]))
    }
    if (replay_valid && !raw_identical(current, blob_value(rows$result_bytes[first]))) {
      failures[[length(failures) + 1L]] <- list(
        model_id = rows$model_id[[first]], source_namespace = rows$source_namespace[[first]],
        owner_key = as.character(rows$owner_key[[first]]), edit_key = NA_character_,
        failure = "result_state_mismatch"
      )
    }
  }
  if (length(failures)) {
    ord <- order(
      vapply(failures, `[[`, character(1), "model_id"),
      vapply(failures, `[[`, character(1), "source_namespace"),
      vapply(failures, `[[`, character(1), "owner_key"),
      vapply(failures, `[[`, character(1), "failure")
    )
    failures <- failures[ord]
    add_issue(
      validator, "DVMA3404", table,
      "applied edit replay disagrees with reference snapshots or final sequence state",
      length(failures), head(failures, validator$max_samples)
    )
  }
}

check_edits_and_annotations <- function(validator) {
  check_edit_relation(
    validator, "model_transcript_edit", "model_transcript", "transcript_key",
    "model_transcript_attribute", "raw_spliced", c("_rna_edit")
  )
  check_edit_relation(
    validator, "model_translation_edit", "model_translation", "translation_key",
    "model_translation_attribute", "peptide_unedited",
    c("initial_met", "_selenocysteine", "amino_acid_sub", "_stop_codon_rt")
  )
  check_edit_replay_one(
    validator, "model_transcript_edit", "transcript_key",
    "raw_spliced", "edited_spliced"
  )
  check_edit_replay_one(
    validator, "model_translation_edit", "translation_key",
    "peptide_unedited", "peptide_final"
  )

  if (relations_ready(validator, c(
      "model_mature_mirna", "model_transcript", "model_transcript_attribute",
      "model_attribute_type", "model_transcript_exon"
  ))) {
    check_query(
      validator, "DVMA3410", "model_mature_mirna",
      "mature-miRNA span has invalid ownership, basis, coordinates, or source attribute",
      sprintf(
        paste(
          "WITH lengths AS (SELECT model_id, source_namespace, transcript_key,",
          "max(raw_cdna_end1) AS raw_length FROM %s",
          "GROUP BY model_id, source_namespace, transcript_key)",
          "SELECT m.model_id, m.source_namespace, m.feature_key FROM %s m",
          "LEFT JOIN %s t USING (model_id, source_namespace, transcript_key)",
          "LEFT JOIN %s a ON a.model_id=m.model_id",
          "AND a.source_namespace=m.source_namespace",
          "AND a.attribute_key=m.source_attribute_key",
          "LEFT JOIN %s ty ON ty.model_id=a.model_id",
          "AND ty.source_namespace=a.source_namespace AND ty.attrib_type_key=a.attrib_type_key",
          "LEFT JOIN lengths l USING (model_id, source_namespace, transcript_key)",
          "WHERE m.feature_key=0 OR t.transcript_key IS NULL OR a.attribute_key IS NULL",
          "OR a.transcript_key<>m.transcript_key OR ty.code<>'miRNA'",
          "OR m.coordinate_basis<>'raw_spliced' OR m.cdna_start1<1",
          "OR m.cdna_end1<m.cdna_start1 OR m.cdna_end1>l.raw_length"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_mature_mirna"),
        relation_sql("model_transcript"), relation_sql("model_transcript_attribute"),
        relation_sql("model_attribute_type")
      ),
      "model_id, source_namespace, feature_key"
    )
  }
}

check_xrefs <- function(validator) {
  if (validator$ready[["model_external_db"]]) {
    check_query(
      validator, "DVMA3500", "model_external_db",
      "external-database source identity or metadata is invalid",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, external_db_key FROM %s WHERE",
          "external_db_key=0 OR length(trim(source_internal_id))=0",
          "OR length(trim(db_name))=0 OR priority<0",
          "OR status NOT IN ('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO')",
          "OR type NOT IN",
          "('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT',",
          "'PRIMARY_DB_SYNONYM','ENSEMBL')"
        ),
        relation_sql("model_external_db")
      ),
      "model_id, source_namespace, external_db_key"
    )
  }

  if (relations_ready(validator, c("model_xref", "model_external_db"))) {
    check_query(
      validator, "DVMA3510", "model_xref",
      "xref source identity, accession, or external-database reference is invalid",
      sprintf(
        paste(
          "SELECT x.model_id, x.source_namespace, x.xref_key FROM %s x",
          "LEFT JOIN %s d USING (model_id, source_namespace, external_db_key)",
          "WHERE x.xref_key=0 OR length(trim(x.source_internal_id))=0",
          "OR length(trim(x.accession))=0 OR d.external_db_key IS NULL",
          "OR (x.accession_version IS NOT NULL AND",
          "length(trim(x.accession_version))=0)",
          "OR x.info_type NOT IN",
          "('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH',",
          "'INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM')"
        ),
        relation_sql("model_xref"), relation_sql("model_external_db")
      ),
      "model_id, source_namespace, xref_key"
    )
  }

  if (relations_ready(validator, c(
      "model_object_xref", "model_source", "model_xref",
      "model_gene", "model_transcript", "model_translation"
  ))) {
    check_query(
      validator, "DVMA3520", "model_object_xref",
      "object-xref source namespace, xref, or typed object owner is invalid",
      sprintf(
        paste(
          "SELECT o.model_id, o.object_source_namespace, o.xref_source_namespace,",
          "o.object_xref_key FROM %s o",
          "LEFT JOIN %s os ON os.model_id=o.model_id",
          "AND os.source_namespace=o.object_source_namespace",
          "LEFT JOIN %s xs ON xs.model_id=o.model_id",
          "AND xs.source_namespace=o.xref_source_namespace",
          "LEFT JOIN %s x ON x.model_id=o.model_id",
          "AND x.source_namespace=o.xref_source_namespace AND x.xref_key=o.xref_key",
          "LEFT JOIN %s g ON o.object_kind='gene' AND g.model_id=o.model_id",
          "AND g.source_namespace=o.object_source_namespace AND g.gene_key=o.object_key",
          "LEFT JOIN %s t ON o.object_kind='transcript' AND t.model_id=o.model_id",
          "AND t.source_namespace=o.object_source_namespace AND t.transcript_key=o.object_key",
          "LEFT JOIN %s p ON o.object_kind='translation' AND p.model_id=o.model_id",
          "AND p.source_namespace=o.object_source_namespace AND p.translation_key=o.object_key",
          "WHERE o.object_xref_key=0 OR length(trim(o.source_internal_id))=0",
          "OR o.object_kind NOT IN ('gene','transcript','translation')",
          "OR os.source_namespace IS NULL OR xs.source_namespace IS NULL OR x.xref_key IS NULL",
          "OR (o.object_kind='gene' AND g.gene_key IS NULL)",
          "OR (o.object_kind='transcript' AND t.transcript_key IS NULL)",
          "OR (o.object_kind='translation' AND p.translation_key IS NULL)"
        ),
        relation_sql("model_object_xref"), relation_sql("model_source"),
        relation_sql("model_source"), relation_sql("model_xref"),
        relation_sql("model_gene"), relation_sql("model_transcript"),
        relation_sql("model_translation")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
    check_query(
      validator, "DVMA3521", "model_object_xref",
      "display-xref marker disagrees with the source gene/transcript display pointer",
      sprintf(
        paste(
          "SELECT o.model_id, o.object_source_namespace, o.xref_source_namespace,",
          "o.object_xref_key, o.object_kind, o.object_key, o.xref_key",
          "FROM %s o LEFT JOIN %s g ON o.object_kind='gene'",
          "AND g.model_id=o.model_id AND g.source_namespace=o.object_source_namespace",
          "AND g.gene_key=o.object_key",
          "LEFT JOIN %s t ON o.object_kind='transcript'",
          "AND t.model_id=o.model_id AND t.source_namespace=o.object_source_namespace",
          "AND t.transcript_key=o.object_key",
          "WHERE (o.object_kind='gene' AND o.is_display_xref",
          "IS DISTINCT FROM (o.xref_source_namespace=o.object_source_namespace",
          "AND g.display_xref_key IS NOT NULL AND g.display_xref_key=o.xref_key))",
          "OR (o.object_kind='transcript' AND o.is_display_xref",
          "IS DISTINCT FROM (o.xref_source_namespace=o.object_source_namespace",
          "AND t.display_xref_key IS NOT NULL AND t.display_xref_key=o.xref_key))",
          "OR (o.object_kind='translation' AND o.is_display_xref)"
        ),
        relation_sql("model_object_xref"), relation_sql("model_gene"),
        relation_sql("model_transcript")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
    check_query(
      validator, "DVMA3522", "model_object_xref",
      "every gene/transcript display pointer must resolve to exactly one matching object_xref",
      sprintf(
        paste(
          "WITH expected AS (SELECT model_id,source_namespace,'gene' AS object_kind,",
          "gene_key AS object_key,display_xref_key AS xref_key FROM %s",
          "WHERE display_xref_key IS NOT NULL UNION ALL",
          "SELECT model_id,source_namespace,'transcript' AS object_kind,",
          "transcript_key AS object_key,display_xref_key AS xref_key FROM %s",
          "WHERE display_xref_key IS NOT NULL)",
          "SELECT e.model_id,e.source_namespace,e.object_kind,e.object_key,e.xref_key,",
          "count(o.object_xref_key) AS matching_rows FROM expected e LEFT JOIN %s o",
          "ON o.model_id=e.model_id AND o.object_source_namespace=e.source_namespace",
          "AND o.xref_source_namespace=e.source_namespace",
          "AND o.object_kind=e.object_kind AND o.object_key=e.object_key",
          "AND o.xref_key=e.xref_key AND o.is_display_xref",
          "GROUP BY e.model_id,e.source_namespace,e.object_kind,e.object_key,e.xref_key",
          "HAVING count(o.object_xref_key)<>1"
        ),
        relation_sql("model_gene"), relation_sql("model_transcript"),
        relation_sql("model_object_xref")
      ),
      "model_id, source_namespace, object_kind, object_key"
    )
  }

  if (relations_ready(validator, c(
      "model_seq_region_synonym", "model_seq_region", "model_external_db"
  ))) {
    check_query(
      validator, "DVMA3530", "model_seq_region_synonym",
      "sequence-region synonym has invalid source identity, owner, or optional external database",
      sprintf(
        paste(
          "SELECT s.model_id, s.source_namespace, s.synonym_key FROM %s s",
          "LEFT JOIN %s r USING (model_id, source_namespace, seq_region_key)",
          "LEFT JOIN %s d ON d.model_id=s.model_id",
          "AND d.source_namespace=s.source_namespace AND d.external_db_key=s.external_db_key",
          "WHERE s.synonym_key=0 OR length(trim(s.source_internal_id))=0",
          "OR length(trim(s.synonym))=0 OR r.seq_region_key IS NULL",
          "OR (s.external_db_key IS NOT NULL AND d.external_db_key IS NULL)"
        ),
        relation_sql("model_seq_region_synonym"), relation_sql("model_seq_region"),
        relation_sql("model_external_db")
      ),
      "model_id, source_namespace, synonym_key"
    )
  }

  if (relations_ready(validator, c("model_xref_identity", "model_object_xref", "model_xref"))) {
    check_query(
      validator, "DVMA3540", "model_xref_identity",
      "xref identity does not resolve to its distinct object_xref source row",
      sprintf(
        paste(
          "SELECT i.model_id, i.object_source_namespace, i.xref_source_namespace,",
          "i.object_xref_key FROM %s i LEFT JOIN %s o",
          "USING (model_id, object_source_namespace, xref_source_namespace, object_xref_key)",
          "WHERE o.object_xref_key IS NULL",
          "OR i.identity_status NOT IN ('unverified','exact','mismatch')",
          "OR (i.identity_status='unverified' AND",
          "(i.verified_sequence_role IS NOT NULL OR i.verified_sequence_sha256 IS NOT NULL))",
          "OR (i.identity_status IN ('exact','mismatch') AND",
          "(i.verified_sequence_role IS NULL OR",
          "NOT regexp_full_match(i.verified_sequence_sha256, '[0-9a-f]{64}')))"
        ),
        relation_sql("model_xref_identity"), relation_sql("model_object_xref")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
    check_query(
      validator, "DVMA3541", "model_object_xref",
      "every retained object_xref requires exactly one explicit identity status",
      sprintf(
        paste(
          "SELECT o.model_id, o.object_source_namespace, o.xref_source_namespace,",
          "o.object_xref_key, count(i.object_xref_key) AS identity_count FROM %s o",
          "LEFT JOIN %s i USING",
          "(model_id, object_source_namespace, xref_source_namespace, object_xref_key)",
          "GROUP BY o.model_id, o.object_source_namespace,",
          "o.xref_source_namespace, o.object_xref_key",
          "HAVING count(i.object_xref_key)<>1"
        ),
        relation_sql("model_object_xref"), relation_sql("model_xref_identity")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
    check_query(
      validator, "DVMA3542", "model_xref_identity",
      "sequence identity comparison requires a versioned accession",
      sprintf(
        paste(
          "SELECT i.model_id, i.object_source_namespace, i.xref_source_namespace,",
          "i.object_xref_key FROM %s i JOIN %s o",
          "USING (model_id, object_source_namespace, xref_source_namespace, object_xref_key)",
          "JOIN %s x ON x.model_id=o.model_id",
          "AND x.source_namespace=o.xref_source_namespace AND x.xref_key=o.xref_key",
          "WHERE i.identity_status IN ('exact','mismatch') AND",
          "(x.accession_version IS NULL OR length(trim(x.accession_version))=0)"
        ),
        relation_sql("model_xref_identity"), relation_sql("model_object_xref"),
        relation_sql("model_xref")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
  }

  if (relations_ready(validator, c(
      "model_xref_identity", "model_object_xref", "model_sequence_state"
  ))) {
    check_query(
      validator, "DVMA3543", "model_xref_identity",
      "sequence identity status disagrees with the named model sequence state",
      sprintf(
        paste(
          "SELECT i.model_id,i.object_source_namespace,i.xref_source_namespace,",
          "i.object_xref_key,i.identity_status,i.verified_sequence_role",
          "FROM %s i JOIN %s o USING",
          "(model_id,object_source_namespace,xref_source_namespace,object_xref_key)",
          "LEFT JOIN %s s ON s.model_id=o.model_id",
          "AND s.source_namespace=o.object_source_namespace",
          "AND s.owner_kind=o.object_kind AND s.owner_key=o.object_key",
          "AND s.role=i.verified_sequence_role AND s.status='present'",
          "WHERE i.identity_status IN ('exact','mismatch') AND",
          "(o.object_kind NOT IN ('transcript','translation')",
          "OR s.sequence_sha256 IS NULL",
          "OR (o.object_kind='transcript' AND i.verified_sequence_role",
          "NOT IN ('raw_spliced','edited_spliced'))",
          "OR (o.object_kind='translation' AND i.verified_sequence_role",
          "NOT IN ('translatable_cds','peptide_unedited','peptide_final'))",
          "OR (i.identity_status='exact' AND",
          "i.verified_sequence_sha256<>s.sequence_sha256)",
          "OR (i.identity_status='mismatch' AND",
          "i.verified_sequence_sha256=s.sequence_sha256))"
        ),
        relation_sql("model_xref_identity"), relation_sql("model_object_xref"),
        relation_sql("model_sequence_state")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
  }
}

check_sequences <- function(validator) {
  if (validator$ready[["model_sequence_blob"]]) {
    check_query(
      validator, "DVMA3600", "model_sequence_blob",
      "sequence blob digest, alphabet, byte count, or exact uppercase bytes are invalid",
      sprintf(
        paste(
          "WITH decoded AS (SELECT *,try(decode(sequence_bytes)) AS sequence_text FROM %s)",
          "SELECT sequence_sha256, alphabet, byte_count,",
          "octet_length(sequence_bytes) AS actual_byte_count FROM decoded WHERE",
          "NOT regexp_full_match(sequence_sha256, '[0-9a-f]{64}')",
          "OR alphabet NOT IN ('dna','peptide') OR byte_count<0",
          "OR byte_count<>octet_length(sequence_bytes)",
          "OR sha256(sequence_bytes)<>sequence_sha256",
          "OR sequence_text IS NULL",
          "OR (alphabet='dna' AND NOT regexp_full_match(sequence_text, '[A-Z]*'))",
          "OR (alphabet='peptide' AND NOT regexp_full_match(sequence_text, '[A-Z*]*'))"
        ),
        relation_sql("model_sequence_blob")
      ),
      "sequence_sha256, alphabet"
    )
  }

  if (relations_ready(validator, c(
      "model_sequence_state", "model_transcript", "model_translation"
  ))) {
    check_query(
      validator, "DVMA3610", "model_sequence_state",
      "sequence-state tagged owner, role, alphabet, status, or key coupling is invalid",
      sprintf(
        paste(
          "SELECT model_id, source_namespace, owner_kind, owner_key, role FROM %s WHERE",
          "owner_kind NOT IN ('transcript','translation') OR owner_key=0",
          "OR status NOT IN ('present','absent','provider_error')",
          "OR NOT ((owner_kind='transcript'",
          "AND role IN ('raw_spliced','edited_spliced') AND alphabet='dna')",
          "OR (owner_kind='translation' AND role='translatable_cds' AND alphabet='dna')",
          "OR (owner_kind='translation'",
          "AND role IN ('peptide_unedited','peptide_final') AND alphabet='peptide'))",
          "OR ((status='present') <> (sequence_sha256 IS NOT NULL))",
          "OR ((status='present') <> (byte_count IS NOT NULL))",
          "OR (byte_count IS NOT NULL AND byte_count<0)"
        ),
        relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, owner_kind, owner_key, role"
    )
    check_query(
      validator, "DVMA3611", "model_sequence_state",
      "sequence state does not resolve to its tagged transcript/translation owner",
      sprintf(
        paste(
          "SELECT s.model_id, s.source_namespace, s.owner_kind, s.owner_key, s.role",
          "FROM %s s LEFT JOIN %s t ON s.owner_kind='transcript'",
          "AND t.model_id=s.model_id AND t.source_namespace=s.source_namespace",
          "AND t.transcript_key=s.owner_key",
          "LEFT JOIN %s x ON s.owner_kind='translation'",
          "AND x.model_id=s.model_id AND x.source_namespace=s.source_namespace",
          "AND x.translation_key=s.owner_key",
          "WHERE (s.owner_kind='transcript' AND t.transcript_key IS NULL)",
          "OR (s.owner_kind='translation' AND x.translation_key IS NULL)"
        ),
        relation_sql("model_sequence_state"), relation_sql("model_transcript"),
        relation_sql("model_translation")
      ),
      "model_id, source_namespace, owner_kind, owner_key, role"
    )
  }

  if (relations_ready(validator, c("model_sequence_state", "model_sequence_blob"))) {
    check_query(
      validator, "DVMA3612", "model_sequence_state",
      "present sequence state does not resolve by composite digest/alphabet or byte count",
      sprintf(
        paste(
          "SELECT s.model_id, s.source_namespace, s.owner_kind, s.owner_key, s.role,",
          "s.sequence_sha256, s.alphabet FROM %s s LEFT JOIN %s b",
          "USING (sequence_sha256, alphabet) WHERE s.status='present'",
          "AND (b.sequence_sha256 IS NULL OR b.byte_count<>s.byte_count)"
        ),
        relation_sql("model_sequence_state"), relation_sql("model_sequence_blob")
      ),
      "model_id, source_namespace, owner_kind, owner_key, role"
    )
  }

  if (relations_ready(validator, c("model_sequence_state", "model_transcript_exon"))) {
    check_query(
      validator, "DVMA3613", "model_sequence_state",
      "raw_spliced byte count disagrees with source-faithful exon/cDNA geometry",
      sprintf(
        paste(
          "WITH lengths AS (SELECT model_id,source_namespace,transcript_key,",
          "max(raw_cdna_end1) AS raw_length FROM %s",
          "GROUP BY model_id,source_namespace,transcript_key)",
          "SELECT s.model_id,s.source_namespace,s.owner_key,s.role,s.byte_count,l.raw_length",
          "FROM %s s JOIN lengths l ON l.model_id=s.model_id",
          "AND l.source_namespace=s.source_namespace AND l.transcript_key=s.owner_key",
          "WHERE s.owner_kind='transcript' AND s.role='raw_spliced'",
          "AND s.status='present' AND s.byte_count<>l.raw_length"
        ),
        relation_sql("model_transcript_exon"), relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, owner_key"
    )
  }

  if (relations_ready(validator, c("model_sequence_state", "model_transcript_edit"))) {
    check_query(
      validator, "DVMA3614", "model_sequence_state",
      "transcript without applied nucleotide edits must have identical raw and edited states",
      sprintf(
        paste(
          "WITH states AS (SELECT model_id,source_namespace,owner_key,",
          "max(sequence_sha256) FILTER (WHERE role='raw_spliced' AND status='present') AS raw_sha,",
          "max(byte_count) FILTER (WHERE role='raw_spliced' AND status='present') AS raw_count,",
          "max(sequence_sha256) FILTER (WHERE role='edited_spliced' AND status='present') AS edited_sha,",
          "max(byte_count) FILTER (WHERE role='edited_spliced' AND status='present') AS edited_count",
          "FROM %s WHERE owner_kind='transcript'",
          "GROUP BY model_id,source_namespace,owner_key),",
          "edits AS (SELECT model_id,source_namespace,transcript_key AS owner_key,count(*) AS n",
          "FROM %s WHERE status='applied' GROUP BY model_id,source_namespace,transcript_key)",
          "SELECT s.model_id,s.source_namespace,s.owner_key,s.raw_sha,s.edited_sha",
          "FROM states s LEFT JOIN edits e USING (model_id,source_namespace,owner_key)",
          "WHERE coalesce(e.n,0)=0 AND s.raw_sha IS NOT NULL AND s.edited_sha IS NOT NULL",
          "AND (s.raw_sha<>s.edited_sha OR s.raw_count<>s.edited_count)"
        ),
        relation_sql("model_sequence_state"), relation_sql("model_transcript_edit")
      ),
      "model_id, source_namespace, owner_key"
    )
  }

  if (relations_ready(validator, c("model_sequence_state", "model_translation"))) {
    check_query(
      validator, "DVMA3615", "model_sequence_state",
      "edited coding bounds exceed edited_spliced sequence state",
      sprintf(
        paste(
          "SELECT x.model_id,x.source_namespace,x.translation_key,",
          "x.edited_cdna_coding_end1,s.byte_count AS edited_length FROM %s x",
          "LEFT JOIN %s s ON s.model_id=x.model_id",
          "AND s.source_namespace=x.source_namespace AND s.owner_kind='transcript'",
          "AND s.owner_key=x.transcript_key AND s.role='edited_spliced'",
          "AND s.status='present' WHERE s.owner_key IS NOT NULL",
          "AND x.edited_cdna_coding_end1>s.byte_count"
        ),
        relation_sql("model_translation"), relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, translation_key"
    )
    check_query(
      validator, "DVMA3616", "model_sequence_state",
      "translatable_cds byte count disagrees with edited coding span plus phase padding",
      sprintf(
        paste(
          "SELECT x.model_id,x.source_namespace,x.translation_key,s.byte_count,",
          "CAST(x.edited_cdna_coding_end1 AS HUGEINT)",
          "-CAST(x.edited_cdna_coding_start1 AS HUGEINT)+1",
          "+CAST(x.start_phase_padding AS HUGEINT) AS expected_count FROM %s x JOIN %s s",
          "ON s.model_id=x.model_id AND s.source_namespace=x.source_namespace",
          "AND s.owner_kind='translation' AND s.owner_key=x.translation_key",
          "AND s.role='translatable_cds' AND s.status='present'",
          "WHERE s.byte_count<x.start_phase_padding",
          "OR s.byte_count-x.start_phase_padding<>",
          "x.edited_cdna_coding_end1-x.edited_cdna_coding_start1+1"
        ),
        relation_sql("model_translation"), relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, translation_key"
    )
    if (validator$ready[["model_sequence_blob"]]) {
      check_query(
        validator, "DVMA3618", "model_sequence_state",
        "translatable_cds bytes disagree with edited transcript bounds and explicit phase padding",
        sprintf(
          paste(
            "SELECT x.model_id,x.source_namespace,x.translation_key FROM %s x",
            "JOIN %s es ON es.model_id=x.model_id",
            "AND es.source_namespace=x.source_namespace AND es.owner_kind='transcript'",
            "AND es.owner_key=x.transcript_key AND es.role='edited_spliced'",
            "AND es.status='present' JOIN %s eb",
            "ON eb.sequence_sha256=es.sequence_sha256 AND eb.alphabet=es.alphabet",
            "JOIN %s cs ON cs.model_id=x.model_id",
            "AND cs.source_namespace=x.source_namespace AND cs.owner_kind='translation'",
            "AND cs.owner_key=x.translation_key AND cs.role='translatable_cds'",
            "AND cs.status='present' JOIN %s cb",
            "ON cb.sequence_sha256=cs.sequence_sha256 AND cb.alphabet=cs.alphabet",
            "WHERE x.start_phase_padding BETWEEN 0 AND 2",
            "AND x.edited_cdna_coding_start1>=1",
            "AND x.edited_cdna_coding_end1>=x.edited_cdna_coding_start1",
            "AND x.edited_cdna_coding_end1<=es.byte_count AND",
            "(try(decode(eb.sequence_bytes)) IS NULL",
            "OR try(decode(cb.sequence_bytes)) IS NULL",
            "OR try(decode(cb.sequence_bytes))<>",
            "repeat('N',x.start_phase_padding)||substr(try(decode(eb.sequence_bytes)),",
            "x.edited_cdna_coding_start1,",
            "x.edited_cdna_coding_end1-x.edited_cdna_coding_start1+1))"
          ),
          relation_sql("model_translation"), relation_sql("model_sequence_state"),
          relation_sql("model_sequence_blob"), relation_sql("model_sequence_state"),
          relation_sql("model_sequence_blob")
        ),
        "model_id, source_namespace, translation_key"
      )
    }
  }

  if (relations_ready(validator, c("model_sequence_state", "model_translation_edit"))) {
    check_query(
      validator, "DVMA3617", "model_sequence_state",
      "translation without applied peptide edits must have identical unedited/final peptide states",
      sprintf(
        paste(
          "WITH states AS (SELECT model_id,source_namespace,owner_key,",
          "max(sequence_sha256) FILTER (WHERE role='peptide_unedited' AND status='present') AS raw_sha,",
          "max(byte_count) FILTER (WHERE role='peptide_unedited' AND status='present') AS raw_count,",
          "max(sequence_sha256) FILTER (WHERE role='peptide_final' AND status='present') AS final_sha,",
          "max(byte_count) FILTER (WHERE role='peptide_final' AND status='present') AS final_count",
          "FROM %s WHERE owner_kind='translation'",
          "GROUP BY model_id,source_namespace,owner_key),",
          "edits AS (SELECT model_id,source_namespace,translation_key AS owner_key,count(*) AS n",
          "FROM %s WHERE status='applied' GROUP BY model_id,source_namespace,translation_key)",
          "SELECT s.model_id,s.source_namespace,s.owner_key,s.raw_sha,s.final_sha",
          "FROM states s LEFT JOIN edits e USING (model_id,source_namespace,owner_key)",
          "WHERE coalesce(e.n,0)=0 AND s.raw_sha IS NOT NULL AND s.final_sha IS NOT NULL",
          "AND (s.raw_sha<>s.final_sha OR s.raw_count<>s.final_count)"
        ),
        relation_sql("model_sequence_state"), relation_sql("model_translation_edit")
      ),
      "model_id, source_namespace, owner_key"
    )
  }

  if (validator$profile != "structural" && relations_ready(validator, c(
      "model_sequence_state", "model_transcript", "model_translation"
  ))) {
    check_query(
      validator, "DVMA3620", "model_sequence_state",
      "every transcript requires present raw_spliced and edited_spliced states",
      sprintf(
        paste(
          "WITH required(role,alphabet) AS",
          "(VALUES ('raw_spliced','dna'),('edited_spliced','dna'))",
          "SELECT t.model_id,t.source_namespace,'transcript' AS owner_kind,",
          "t.transcript_key AS owner_key,r.role FROM %s t CROSS JOIN required r",
          "LEFT JOIN %s s ON s.model_id=t.model_id",
          "AND s.source_namespace=t.source_namespace AND s.owner_kind='transcript'",
          "AND s.owner_key=t.transcript_key AND s.role=r.role",
          "AND s.alphabet=r.alphabet AND s.status='present' WHERE s.owner_key IS NULL"
        ),
        relation_sql("model_transcript"), relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, owner_key, role"
    )
    check_query(
      validator, "DVMA3621", "model_sequence_state",
      "every translation requires present CDS, unedited-peptide, and final-peptide states",
      sprintf(
        paste(
          "WITH required(role,alphabet) AS (VALUES ('translatable_cds','dna'),",
          "('peptide_unedited','peptide'),('peptide_final','peptide'))",
          "SELECT x.model_id,x.source_namespace,'translation' AS owner_kind,",
          "x.translation_key AS owner_key,r.role FROM %s x CROSS JOIN required r",
          "LEFT JOIN %s s ON s.model_id=x.model_id",
          "AND s.source_namespace=x.source_namespace AND s.owner_kind='translation'",
          "AND s.owner_key=x.translation_key AND s.role=r.role",
          "AND s.alphabet=r.alphabet AND s.status='present' WHERE s.owner_key IS NULL"
        ),
        relation_sql("model_translation"), relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, owner_key, role"
    )
    check_query(
      validator, "DVMA3622", "model_sequence_state",
      "model candidates cannot admit absent or provider-error required sequence states",
      sprintf(
        "SELECT model_id,source_namespace,owner_kind,owner_key,role,status FROM %s WHERE status<>'present'",
        relation_sql("model_sequence_state")
      ),
      "model_id, source_namespace, owner_kind, owner_key, role"
    )
  }

  if (relations_ready(validator, c(
      "model_xref_identity", "model_object_xref", "model_sequence_state"
  ))) {
    check_query(
      validator, "DVMA3630", "model_xref_identity",
      "verified xref sequence identity does not match the named model sequence state",
      sprintf(
        paste(
          "SELECT i.model_id,i.object_source_namespace,i.xref_source_namespace,",
          "i.object_xref_key FROM %s i JOIN %s o",
          "USING (model_id,object_source_namespace,xref_source_namespace,object_xref_key)",
          "LEFT JOIN %s s ON s.model_id=i.model_id",
          "AND s.source_namespace=i.object_source_namespace",
          "AND s.owner_kind=o.object_kind AND s.owner_key=o.object_key",
          "AND s.role=i.verified_sequence_role AND s.status='present'",
          "WHERE i.identity_status IN ('exact','mismatch') AND",
          "(NOT ((o.object_kind='transcript' AND",
          "i.verified_sequence_role IN ('raw_spliced','edited_spliced')) OR",
          "(o.object_kind='translation' AND",
          "i.verified_sequence_role IN ('peptide_unedited','peptide_final'))) OR",
          "s.owner_key IS NULL OR (i.identity_status='exact' AND",
          "s.sequence_sha256<>i.verified_sequence_sha256) OR",
          "(i.identity_status='mismatch' AND",
          "s.sequence_sha256=i.verified_sequence_sha256))"
        ),
        relation_sql("model_xref_identity"), relation_sql("model_object_xref"),
        relation_sql("model_sequence_state")
      ),
      "model_id, object_source_namespace, xref_source_namespace, object_xref_key"
    )
  }
}

check_deterministic_source_keys <- function(validator) {
  specs <- list(
    model_seq_region = c("source_namespace", "seq_region_key"),
    model_gene = c("source_namespace", "gene_key"),
    model_transcript = c("source_namespace", "transcript_key"),
    model_exon = c("source_namespace", "exon_key"),
    model_translation = c("source_namespace", "translation_key"),
    model_attribute_type = c("source_namespace", "attrib_type_key"),
    model_external_db = c("source_namespace", "external_db_key"),
    model_xref = c("source_namespace", "xref_key"),
    model_seq_region_synonym = c("source_namespace", "synonym_key")
  )
  for (table in names(specs)) {
    if (!validator$ready[[table]]) next
    namespace <- specs[[table]][[1L]]
    key <- specs[[table]][[2L]]
    check_query(
      validator, "DVMA3700", table,
      "pack-local source key is not a bijection assigned by canonical source identity order",
      sprintf(
        paste(
          "WITH numbered AS (SELECT model_id,%s,source_internal_id,%s,",
          "row_number() OVER (PARTITION BY model_id,%s",
          "ORDER BY source_internal_id) AS expected_key,",
          "count(*) OVER (PARTITION BY model_id,%s,source_internal_id) AS identity_count",
          "FROM %s) SELECT model_id,%s AS source_namespace,source_internal_id,",
          "%s AS actual_key,expected_key,identity_count FROM numbered",
          "WHERE %s<>expected_key OR identity_count<>1"
        ),
        sql_ident(namespace), sql_ident(key), sql_ident(namespace), sql_ident(namespace),
        relation_sql(table), sql_ident(namespace), sql_ident(key), sql_ident(key)
      ),
      "model_id, source_namespace, source_internal_id"
    )
  }

  if (validator$ready[["model_object_xref"]]) {
    check_query(
      validator, "DVMA3700", "model_object_xref",
      "object_xref key is not a bijection assigned by canonical source identity order",
      sprintf(
        paste(
          "WITH numbered AS (SELECT *,row_number() OVER",
          "(PARTITION BY model_id,object_source_namespace,xref_source_namespace",
          "ORDER BY source_internal_id) AS expected_key,",
          "count(*) OVER (PARTITION BY model_id,object_source_namespace,",
          "xref_source_namespace,source_internal_id) AS identity_count FROM %s)",
          "SELECT model_id,object_source_namespace,xref_source_namespace,",
          "source_internal_id,object_xref_key,expected_key,identity_count FROM numbered",
          "WHERE object_xref_key<>expected_key OR identity_count<>1"
        ),
        relation_sql("model_object_xref")
      ),
      "model_id, object_source_namespace, xref_source_namespace, source_internal_id"
    )
  }
}

finish_validation <- function(validator) {
  if (validator$profile == "conformance_gate") {
    add_issue(
      validator, "DVMA2901", "model_manifest",
      paste(
        "canonical model_id/model_build_id authentication is unavailable:",
        "canonical row serialization and golden vectors are not yet specified"
      )
    )
    add_issue(
      validator, "DVMA2902", "external_attestation",
      paste(
        "external physical-pack attestation verification is unavailable:",
        "an enclosing DuckDB file or Parquet bundle cannot hash itself"
      )
    )
  }
  if (length(validator$issues)) {
    keys <- vapply(validator$issues, function(issue) {
      paste(issue$code, issue$relation, issue$message, sep = "\r")
    }, character(1))
    validator$issues <- validator$issues[order(keys, method = "radix")]
  }
  failure_count <- sum(vapply(
    validator$issues, function(issue) as.numeric(issue$count), numeric(1)
  ))
  list(
    contract = CONTRACT,
    failure_count = failure_count,
    identity_authentication = "unavailable",
    issue_count = length(validator$issues),
    issues = validator$issues,
    ok = length(validator$issues) == 0L,
    physical_attestation = "external_required",
    profile = validator$profile,
    relations = validator$row_counts,
    format_version = FORMAT_VERSION,
    validator_version = VALIDATOR_VERSION
  )
}

validate_connection <- function(connection, expected, profile, max_samples, pack_kind) {
  validator <- new_validator(connection, expected, profile, max_samples, pack_kind)
  inspect_shape(validator)
  check_generic_constraints(validator)
  check_manifest_sources_build(validator)
  check_artifacts_and_pack_references(validator)
  check_nonempty_candidate(validator)
  check_geometry(validator)
  check_membership_and_translation(validator)
  check_attributes(validator)
  check_edits_and_annotations(validator)
  check_xrefs(validator)
  check_sequences(validator)
  check_deterministic_source_keys(validator)
  finish_validation(validator)
}

attach_parquet_pack <- function(connection, directory) {
  DBI::dbExecute(connection, sprintf("CREATE SCHEMA %s", sql_ident(SCHEMA_NAME)))
  paths <- sort(list.files(
    directory, pattern = "\\.parquet$", full.names = TRUE, recursive = FALSE
  ), method = "radix")
  for (path in paths) {
    table <- sub("\\.parquet$", "", basename(path))
    DBI::dbExecute(
      connection,
      sprintf(
        "CREATE VIEW %s AS SELECT * FROM read_parquet(%s)",
        relation_sql(table), sql_literal(normalizePath(path, mustWork = TRUE))
      )
    )
  }
}

initialize_database <- function(path, pretty) {
  if (file.exists(path) || dir.exists(path)) {
    tool_error("DVMA0002", sprintf("refusing to overwrite existing path: %s", path))
  }
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    tool_error("DVMA0003", sprintf("parent directory does not exist: %s", parent))
  }
  connection <- connect_duckdb(path)
  success <- FALSE
  on.exit({
    disconnect_duckdb(connection)
    if (!success && file.exists(path)) unlink(path)
  }, add = TRUE)
  DBI::dbExecute(connection, read_schema())
  success <- TRUE
  emit_json(list(
    contract = CONTRACT,
    database = normalizePath(path, mustWork = TRUE),
    format_version = FORMAT_VERSION,
    ok = TRUE,
    operation = "init"
  ), pretty)
  0L
}

validate_pack <- function(path, profile, max_samples, pretty) {
  expected <- expected_schema()
  if (dir.exists(path)) {
    connection <- connect_duckdb()
    pack_kind <- "parquet_directory"
    attach_parquet_pack(connection, path)
  } else if (file.exists(path)) {
    connection <- connect_duckdb(path, read_only = TRUE)
    pack_kind <- "duckdb_database"
  } else {
    tool_error("DVMA0004", sprintf("pack path is neither a file nor directory: %s", path))
  }
  on.exit(disconnect_duckdb(connection), add = TRUE)
  result <- validate_connection(connection, expected, profile, max_samples, pack_kind)
  result$operation <- "validate"
  result$pack <- normalizePath(path, mustWork = TRUE)
  result$pack_kind <- pack_kind
  emit_json(result, pretty)
  if (isTRUE(result$ok)) 0L else 1L
}

usage <- function() {
  paste(
    "usage:",
    "  Rscript scripts/duckvep_model.R init DATABASE [--pretty]",
    "  Rscript scripts/duckvep_model.R validate PACK",
    "      [--profile structural|model_candidate|conformance_gate]",
    "      [--max-samples N] [--pretty]",
    "  Rscript scripts/duckvep_model.R schema",
    sep = "\n"
  )
}

parse_cli <- function(args) {
  if (!length(args)) tool_error("DVMA0005", usage())
  command <- args[[1L]]
  rest <- args[-1L]
  pretty <- any(rest == "--pretty")
  rest <- rest[rest != "--pretty"]

  if (command == "schema") {
    if (length(rest)) tool_error("DVMA0005", usage())
    return(list(command = command, pretty = pretty))
  }
  if (command == "init") {
    if (length(rest) != 1L || startsWith(rest[[1L]], "--")) {
      tool_error("DVMA0005", usage())
    }
    return(list(command = command, path = rest[[1L]], pretty = pretty))
  }
  if (command != "validate") tool_error("DVMA0005", usage())

  profile <- "conformance_gate"
  max_samples <- 5L
  positional <- character()
  index <- 1L
  while (index <= length(rest)) {
    value <- rest[[index]]
    if (startsWith(value, "--profile=")) {
      profile <- sub("^--profile=", "", value)
    } else if (value == "--profile") {
      index <- index + 1L
      if (index > length(rest)) tool_error("DVMA0005", usage())
      profile <- rest[[index]]
    } else if (startsWith(value, "--max-samples=")) {
      max_samples <- suppressWarnings(as.integer(sub("^--max-samples=", "", value)))
    } else if (value == "--max-samples") {
      index <- index + 1L
      if (index > length(rest)) tool_error("DVMA0005", usage())
      max_samples <- suppressWarnings(as.integer(rest[[index]]))
    } else if (startsWith(value, "--")) {
      tool_error("DVMA0005", sprintf("unknown option: %s\n%s", value, usage()))
    } else {
      positional <- c(positional, value)
    }
    index <- index + 1L
  }
  if (!profile %in% PROFILES) {
    tool_error("DVMA0005", sprintf("unknown validation profile: %s", profile))
  }
  if (is.na(max_samples) || max_samples < 0L || max_samples > 100L) {
    tool_error("DVMA0005", "--max-samples must be an integer from 0 through 100")
  }
  if (length(positional) != 1L) tool_error("DVMA0005", usage())
  list(
    command = command, path = positional[[1L]], profile = profile,
    max_samples = max_samples, pretty = pretty
  )
}

main <- function() {
  args <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (args$command == "schema") {
    cat(read_schema(), "\n", sep = "")
    return(0L)
  }
  if (args$command == "init") {
    return(initialize_database(args$path, args$pretty))
  }
  validate_pack(args$path, args$profile, args$max_samples, args$pretty)
}

status <- tryCatch(
  main(),
  duckvep_model_tool_error = function(error) {
    pretty <- "--pretty" %in% commandArgs(trailingOnly = TRUE)
    emit_json(list(
      contract = CONTRACT,
      error = list(code = error$code, message = conditionMessage(error)),
      ok = FALSE,
      operation = if (length(commandArgs(trailingOnly = TRUE)))
        commandArgs(trailingOnly = TRUE)[[1L]] else "unknown"
    ), pretty)
    2L
  },
  error = function(error) {
    pretty <- "--pretty" %in% commandArgs(trailingOnly = TRUE)
    emit_json(list(
      contract = CONTRACT,
      error = list(code = "DVMA0001", message = conditionMessage(error)),
      ok = FALSE,
      operation = if (length(commandArgs(trailingOnly = TRUE)))
        commandArgs(trailingOnly = TRUE)[[1L]] else "unknown"
    ), pretty)
    2L
  }
)
quit(save = "no", status = status, runLast = FALSE)
