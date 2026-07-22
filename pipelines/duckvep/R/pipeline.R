duckvep_targets_bool <- function(value, field) {
  value <- tolower(trimws(as.character(value)))
  if (length(value) != 1L || !(value %in% c("true", "false"))) {
    stop(field, " must be true or false", call. = FALSE)
  }
  identical(value, "true")
}

duckvep_targets_path <- function(path, root) {
  if (!nzchar(path) || identical(path, ":memory:")) return(path)
  path <- path.expand(path)
  absolute <- grepl("^/", path) ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    grepl("^\\\\\\\\", path)
  if (!absolute) path <- file.path(root, path)
  duckvep_targets_canonical_path(path)
}

duckvep_targets_canonical_path <- function(path) {
  if (length(path) != 1L || !nzchar(path)) {
    stop("path must be one non-empty string", call. = FALSE)
  }
  path <- path.expand(path)
  absolute <- grepl("^/", path) ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    grepl("^\\\\\\\\", path)
  if (!absolute) path <- file.path(getwd(), path)
  path <- chartr("\\", "/", path)
  prefix <- if (grepl("^[A-Za-z]:/", path)) {
    substring(path, 1L, 3L)
  } else if (startsWith(path, "//")) {
    parts <- strsplit(substring(path, 3L), "/+", perl = TRUE)[[1L]]
    if (length(parts) < 2L) stop("invalid UNC path: ", path, call. = FALSE)
    paste0("//", parts[[1L]], "/", parts[[2L]])
  } else {
    "/"
  }
  remainder <- substring(path, nchar(prefix) + 1L)
  parts <- strsplit(remainder, "/+", perl = TRUE)[[1L]]
  current <- normalizePath(prefix, winslash = "/", mustWork = TRUE)
  for (part in parts) {
    if (!nzchar(part) || identical(part, ".")) next
    if (identical(part, "..")) {
      current <- dirname(current)
    } else {
      candidate <- file.path(current, part)
      current <- if (file.exists(candidate)) {
        normalizePath(candidate, winslash = "/", mustWork = TRUE)
      } else {
        candidate
      }
    }
  }
  current
}

duckvep_targets_read_campaigns <- function(
    path,
    root,
    micromamba_default,
    vep_prefix_default,
    output_root,
    extension_paths,
    targets_store) {
  value <- utils::read.delim(
    path,
    colClasses = "character",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  allowed <- c(
    "id", "enabled", "corpus", "event_mode", "vcf", "gff",
    "gff_index_policy", "cache_dir",
    "cache_info", "cache_receipt", "assembly", "species", "cache_version", "fasta",
    "database", "model_sql", "model_name", "sample_per_shape",
    "max_allele_length", "split_multiallelic", "stratify_raw_allele_length",
    "seed", "chrom", "distance", "max_sv_size", "regulatory", "hgvs",
    "fork", "vep_buffer_size", "vep_prefix", "micromamba",
    "nmd_plugin_dir", "oracle_environment_receipt", "duckdb_memory_limit",
    "duckdb_threads", "source_url", "source_version", "source_checksum"
  )
  unknown <- setdiff(names(value), allowed)
  if (length(unknown)) {
    stop("unknown campaign field(s): ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  if (!("id" %in% names(value)) || !nrow(value)) {
    stop("campaign manifest needs at least one id", call. = FALSE)
  }
  if (!("enabled" %in% names(value))) value$enabled <- "true"
  enabled <- vapply(
    value$enabled,
    duckvep_targets_bool,
    logical(1L),
    field = "enabled"
  )
  value <- value[enabled, , drop = FALSE]
  if (!nrow(value)) stop("campaign manifest enables no campaigns", call. = FALSE)
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", value$id)) ||
      anyDuplicated(value$id)) {
    stop("campaign ids must be unique path-safe identifiers", call. = FALSE)
  }
  path_fields <- intersect(
    c(
      "vcf", "gff", "cache_dir", "cache_info", "cache_receipt", "fasta", "database",
      "model_sql", "nmd_plugin_dir", "oracle_environment_receipt",
      "vep_prefix"
    ),
    names(value)
  )
  for (field in path_fields) {
    value[[field]] <- vapply(
      value[[field]],
      duckvep_targets_path,
      character(1L),
      root = root
    )
  }
  campaigns <- lapply(seq_len(nrow(value)), function(index) {
    row <- as.list(value[index, , drop = FALSE])
    row <- lapply(row, function(item) item[[1L]])
    row <- row[nzchar(vapply(row, as.character, character(1L)))]
    micromamba <- row$micromamba %||% micromamba_default
    if (!file.exists(micromamba) &&
        !grepl("[/\\\\]", micromamba)) {
      micromamba <- unname(Sys.which(micromamba))
    } else {
      micromamba <- duckvep_targets_path(micromamba, root)
    }
    if (!nzchar(micromamba)) stop("micromamba is unavailable", call. = FALSE)
    row$micromamba <- normalizePath(micromamba, mustWork = TRUE)
    vep_prefix <- row$vep_prefix %||% vep_prefix_default
    row$vep_prefix <- duckvep_targets_path(vep_prefix, root)
    row$gff_index_policy <- row$gff_index_policy %||% "ignore"
    if (!(row$gff_index_policy %in% c("require", "ignore"))) {
      stop(
        "targets campaigns require gff_index_policy=require or ignore",
        call. = FALSE
      )
    }
    if (!nzchar(row$oracle_environment_receipt %||% "")) {
      stop("campaigns must declare oracle_environment_receipt", call. = FALSE)
    }
    required <- c("corpus", "event_mode", "fasta", "model_name")
    missing <- required[!nzchar(vapply(
      required,
      function(field) row[[field]] %||% "",
      character(1L)
    ))]
    if (length(missing)) {
      stop(
        "campaign ", row$id, " is missing field(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    if (!nzchar(row$gff %||% "") && !nzchar(row$cache_dir %||% "")) {
      stop("campaign ", row$id, " needs gff or cache_dir", call. = FALSE)
    }
    database <- row$database %||% ""
    model_paths <- c(
      nzchar(database) && !identical(database, ":memory:"),
      nzchar(row$model_sql %||% "")
    )
    if (sum(model_paths) != 1L) {
      stop(
        "campaign ", row$id, " needs exactly one of database or model_sql",
        call. = FALSE
      )
    }
    row
  })
  duckvep_targets_attach_outputs(
    campaigns,
    output_root,
    path,
    extension_paths,
    targets_store,
    root
  )
}

duckvep_targets_source_revision <- function(root) {
  environment <- new.env(parent = baseenv())
  sys.source(
    file.path(root, "scripts", "duckvep_evidence.R"),
    envir = environment
  )
  environment$duckvep_evidence_assert_checkout(
    root,
    context = "DuckVEP targets workflow"
  )
}

duckvep_targets_build_extension <- function(root, source_revision, receipt_path) {
  environment <- new.env(parent = baseenv())
  sys.source(
    file.path(root, "scripts", "duckvep_evidence.R"),
    envir = environment
  )
  extension <- file.path(
    root,
    "build",
    "release",
    "duckhts.duckdb_extension"
  )
  built <- environment$duckvep_evidence_build_extension(
    root,
    extension,
    source_revision
  )
  receipt <- environment$duckvep_evidence_write_extension_receipt(
    receipt_path,
    source_revision,
    built
  )
  c(extension = built$path, receipt = receipt)
}

duckvep_targets_campaign_inputs <- function(campaign) {
  inputs <- character()
  add_file <- function(path, label) {
    if (!nzchar(path) || identical(path, ":memory:")) return(invisible(NULL))
    if (!file.exists(path) || dir.exists(path)) {
      stop(label, " is not a file: ", path, call. = FALSE)
    }
    inputs <<- c(inputs, normalizePath(path, mustWork = TRUE))
    invisible(NULL)
  }
  add_file(campaign$vcf %||% "", "vcf")
  add_file(campaign$gff %||% "", "gff")
  gff_index <- if (identical(campaign$gff_index_policy, "require")) {
    paste0(campaign$gff %||% "", ".tbi")
  } else {
    ""
  }
  if (nzchar(gff_index)) {
    add_file(gff_index, "GFF index")
  }
  add_file(campaign$fasta %||% "", "fasta")
  if (nzchar(campaign$fasta %||% "")) {
    add_file(paste0(campaign$fasta, ".fai"), "fasta index")
    if (grepl("[.](?:bgz|gz)$", campaign$fasta)) {
      add_file(paste0(campaign$fasta, ".gzi"), "fasta gzi")
    }
  }
  add_file(campaign$database %||% "", "model database")
  add_file(campaign$model_sql %||% "", "model SQL")
  add_file(campaign$cache_info %||% "", "VEP cache info")
  add_file(campaign$cache_receipt %||% "", "VEP cache receipt")
  add_file(
    campaign$oracle_environment_receipt %||% "",
    "VEP environment receipt"
  )
  add_file(campaign$micromamba %||% "", "micromamba executable")
  vep_history <- file.path(
    campaign$vep_prefix %||% "",
    "conda-meta",
    "history"
  )
  add_file(vep_history, "VEP environment history")
  if (nzchar(campaign$nmd_plugin_dir %||% "")) {
    add_file(file.path(campaign$nmd_plugin_dir, "NMD.pm"), "NMD plugin")
  }
  if (nzchar(campaign$cache_dir %||% "") &&
      (!nzchar(campaign$cache_info %||% "") ||
        !nzchar(campaign$cache_receipt %||% ""))) {
    stop(
      "cache-backed campaigns must declare cache_info and cache_receipt",
      call. = FALSE
    )
  }
  if (nzchar(campaign$cache_dir %||% "")) {
    expected_cache_info <- duckvep_evidence_cache_info_path(
      campaign$cache_dir,
      campaign$species %||% "homo_sapiens",
      campaign$cache_version %||% "116",
      campaign$assembly %||% "GRCh38"
    )
    actual_cache_info <- normalizePath(
      campaign$cache_info,
      mustWork = TRUE
    )
    if (!identical(actual_cache_info, expected_cache_info)) {
      stop(
        "cache_info must be the cache root's exact species/release/assembly marker",
        call. = FALSE
      )
    }
  }
  unique(inputs)
}

duckvep_targets_cache_state <- function(campaign) {
  if (!nzchar(campaign$cache_dir %||% "")) return("gff")
  state <- duckvep_evidence_validate_cache_receipt(
    campaign$cache_receipt,
    campaign$cache_dir,
    campaign$species %||% "homo_sapiens",
    campaign$cache_version %||% "116",
    campaign$assembly %||% "GRCh38"
  )
  paste(
    state$receipt_sha256,
    state$inventory_sha256,
    state$entries,
    state$bytes,
    sep = ":"
  )
}

`%||%` <- function(left, right) {
  if (is.null(left) || !length(left)) right else left
}

duckvep_targets_campaign_outputs <- function(campaign, output_root) {
  directory <- file.path(output_root, campaign$id)
  stem <- file.path(directory, campaign$id)
  outputs <- c(
    annotations = paste0(stem, "_annotations.parquet"),
    statistical_conformance = paste0(stem, "_statistical_conformance.csv"),
    methodology_audit = paste0(stem, "_methodology_audit.csv"),
    consequence_pairs = paste0(stem, "_pairs.parquet"),
    nmd_conformance = paste0(stem, "_nmd_conformance.csv"),
    sample_vcf = paste0(stem, "_sample.vcf"),
    log = paste0(stem, ".log")
  )
  if (identical(campaign$event_mode %||% "small", "small")) {
    outputs <- c(
      outputs,
      eligibility = paste0(stem, "_eligibility.csv")
    )
  }
  if (duckvep_targets_bool(campaign$hgvs %||% "false", "hgvs")) {
    outputs <- c(
      outputs,
      hgvs_summary = paste0(stem, "_hgvs_conformance.csv"),
      hgvs_pairs = paste0(stem, "_hgvs_pairs.parquet")
    )
  }
  resolved <- normalizePath(unname(outputs), mustWork = FALSE)
  names(resolved) <- names(outputs)
  resolved
}

duckvep_targets_campaign_args <- function(
    campaign,
    extension,
    extension_receipt,
    outputs,
    published_outputs = outputs) {
  values <- c(
    corpus = "corpus",
    event_mode = "event-mode",
    vcf = "vcf",
    gff = "gff",
    gff_index_policy = "gff-index-policy",
    cache_dir = "cache-dir",
    cache_info = "cache-info",
    cache_receipt = "cache-receipt",
    assembly = "assembly",
    species = "species",
    cache_version = "cache-version",
    fasta = "fasta",
    database = "database",
    model_sql = "model-sql",
    model_name = "model-name",
    sample_per_shape = "sample-per-shape",
    max_allele_length = "max-allele-length",
    seed = "seed",
    chrom = "chrom",
    distance = "distance",
    max_sv_size = "max-sv-size",
    fork = "fork",
    vep_buffer_size = "vep-buffer-size",
    vep_prefix = "vep-prefix",
    micromamba = "micromamba",
    nmd_plugin_dir = "nmd-plugin-dir",
    oracle_environment_receipt = "oracle-environment-receipt",
    duckdb_memory_limit = "duckdb-memory-limit",
    duckdb_threads = "duckdb-threads",
    source_url = "source-url",
    source_version = "source-version",
    source_checksum = "source-checksum"
  )
  flags <- c(
    split_multiallelic = "split-multiallelic",
    stratify_raw_allele_length = "stratify-raw-allele-length",
    regulatory = "regulatory",
    hgvs = "hgvs"
  )
  arguments <- character()
  for (field in names(values)) {
    value <- campaign[[field]] %||% ""
    if (nzchar(value)) {
      arguments <- c(arguments, paste0("--", values[[field]]), value)
    }
  }
  for (field in names(flags)) {
    value <- campaign[[field]] %||% "false"
    if (duckvep_targets_bool(value, field)) {
      arguments <- c(arguments, paste0("--", flags[[field]]))
    }
  }
  c(
    arguments,
    "--extension", extension,
    "--extension-build-receipt", extension_receipt,
    "--annotations-out", outputs[["annotations"]],
    "--annotations-label", published_outputs[["annotations"]],
    "--sample-vcf", outputs[["sample_vcf"]],
    "--keep-sample-vcf",
    if ("eligibility" %in% names(outputs)) {
      c("--eligibility-out", outputs[["eligibility"]])
    },
    if ("hgvs_summary" %in% names(outputs)) {
      c(
        "--hgvs-out", outputs[["hgvs_summary"]],
        "--hgvs-pairs-out", outputs[["hgvs_pairs"]]
      )
    }
  )
}

duckvep_targets_extension_paths <- function(extension_bundle) {
  extension <- extension_bundle[endsWith(
    extension_bundle,
    "duckhts.duckdb_extension"
  )]
  receipt <- extension_bundle[endsWith(extension_bundle, "extension.tsv")]
  if (length(extension) != 1L || length(receipt) != 1L) {
    stop("extension bundle must contain one extension and one receipt", call. = FALSE)
  }
  list(extension = unname(extension), receipt = unname(receipt))
}

duckvep_targets_path_key <- function(paths) {
  keys <- unname(vapply(paths, duckvep_targets_canonical_path, character(1L)))
  if (identical(.Platform$OS.type, "windows")) tolower(keys) else keys
}

duckvep_targets_path_within <- function(path, directory) {
  path <- duckvep_targets_path_key(path)
  directory <- duckvep_targets_path_key(directory)
  identical(path, directory) || startsWith(path, paste0(directory, "/"))
}

duckvep_targets_attach_outputs <- function(
    campaigns,
    output_root,
    campaign_manifest,
    extension_bundle,
    targets_store,
    root) {
  outputs <- stats::setNames(
    lapply(
      campaigns,
      duckvep_targets_campaign_outputs,
      output_root = output_root
    ),
    vapply(campaigns, `[[`, character(1L), "id")
  )
  flat_outputs <- unlist(outputs, use.names = FALSE)
  output_keys <- duckvep_targets_path_key(flat_outputs)
  if (anyDuplicated(output_keys)) {
    stop("campaign outputs collide after canonicalization", call. = FALSE)
  }
  all_inputs <- unique(unlist(
    lapply(campaigns, duckvep_targets_campaign_inputs),
    use.names = FALSE
  ))
  protected <- c(campaign_manifest, all_inputs, extension_bundle)
  overlap <- flat_outputs[
    output_keys %in% duckvep_targets_path_key(protected)
  ]
  if (length(overlap)) {
    stop(
      "campaign output overlaps a manifest, input, or shared build artifact: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }
  campaign_directories <- unique(dirname(flat_outputs))
  enclosed_protected <- protected[vapply(
    protected,
    function(path) any(vapply(
      campaign_directories,
      function(directory) duckvep_targets_path_within(path, directory),
      logical(1L)
    )),
    logical(1L)
  )]
  if (length(enclosed_protected)) {
    stop(
      "campaign directory contains a manifest, input, or shared build artifact: ",
      paste(enclosed_protected, collapse = ", "),
      call. = FALSE
    )
  }
  allowed_repo_outputs <- file.path(
    root,
    ".tmp",
    "duckvep_targets",
    "results"
  )
  if (duckvep_targets_path_within(output_root, root) &&
      !duckvep_targets_path_within(output_root, allowed_repo_outputs)) {
    stop(
      "DUCKVEP_TARGETS_OUTPUT may be outside the repository or below ",
      allowed_repo_outputs,
      call. = FALSE
    )
  }
  if (duckvep_targets_path_within(output_root, targets_store) ||
      duckvep_targets_path_within(targets_store, output_root)) {
    stop("campaign output root and targets store must not overlap", call. = FALSE)
  }
  for (index in seq_along(campaigns)) {
    campaigns[[index]]$outputs <- outputs[[campaigns[[index]]$id]]
  }
  campaigns
}

duckvep_targets_backup_directories <- function(final_directory) {
  parent <- dirname(final_directory)
  prefix <- paste0(".previous-", basename(final_directory), "-")
  candidates <- list.files(
    parent,
    all.files = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )
  candidates[
    startsWith(basename(candidates), prefix) & dir.exists(candidates)
  ]
}

duckvep_targets_recover_campaign_directory <- function(final_directory) {
  backups <- duckvep_targets_backup_directories(final_directory)
  if (length(backups) > 1L) {
    stop(
      "multiple interrupted publication backups require inspection: ",
      paste(backups, collapse = ", "),
      call. = FALSE
    )
  }
  if (!length(backups)) return(invisible(NULL))
  if (dir.exists(final_directory)) {
    unlink(backups, recursive = TRUE, force = TRUE)
    if (dir.exists(backups)) {
      stop("cannot remove completed publication backup: ", backups,
        call. = FALSE
      )
    }
  } else if (!file.rename(backups, final_directory)) {
    stop("cannot restore interrupted campaign output: ", backups,
      call. = FALSE
    )
  }
  invisible(NULL)
}

duckvep_targets_publish_outputs <- function(stage_directory, final_directory) {
  parent <- dirname(final_directory)
  if (!identical(dirname(stage_directory), parent)) {
    stop("stage and final campaign directories must be siblings", call. = FALSE)
  }
  backup_directory <- tempfile(
    paste0(".previous-", basename(final_directory), "-"),
    tmpdir = parent
  )
  had_previous <- dir.exists(final_directory)
  if (had_previous && !file.rename(final_directory, backup_directory)) {
    stop("cannot stage previous campaign directory for replacement",
      call. = FALSE
    )
  }
  if (!file.rename(stage_directory, final_directory)) {
    if (had_previous && !file.rename(backup_directory, final_directory)) {
      stop(
        "cannot publish or restore campaign output; backup retained at ",
        backup_directory,
        call. = FALSE
      )
    }
    stop(
      "cannot publish campaign directory; staged files remain in ",
      stage_directory,
      call. = FALSE
    )
  }
  if (had_previous) {
    unlink(backup_directory, recursive = TRUE, force = TRUE)
    if (dir.exists(backup_directory)) {
      stop("published output but cannot remove prior backup: ", backup_directory,
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

duckvep_targets_run_campaign <- function(
    campaign,
    campaign_inputs,
    cache_state,
    extension_bundle,
    root) {
  if (!requireNamespace("blit", quietly = TRUE) ||
      utils::packageVersion("blit") < "0.2.0.9000") {
    stop("blit >= 0.2.0.9000 is required", call. = FALSE)
  }
  invisible(campaign_inputs)
  if (length(cache_state) != 1L || !nzchar(cache_state)) {
    stop("campaign has no validated cache state", call. = FALSE)
  }
  extension_paths <- duckvep_targets_extension_paths(extension_bundle)
  outputs <- campaign$outputs
  if (is.null(outputs)) stop("campaign has no validated outputs", call. = FALSE)
  directory <- dirname(outputs[[1L]])
  parent <- dirname(directory)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  duckvep_targets_recover_campaign_directory(directory)
  stage_directory <- tempfile(
    paste0(".run-", campaign$id, "-"),
    tmpdir = parent
  )
  dir.create(stage_directory, recursive = FALSE, showWarnings = FALSE)
  staged_outputs <- file.path(stage_directory, basename(outputs))
  names(staged_outputs) <- names(outputs)
  arguments <- duckvep_targets_campaign_args(
    campaign,
    extension_paths$extension,
    extension_paths$receipt,
    staged_outputs,
    outputs
  )
  runner <- file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "corpus_differential.R"
  )
  command <- do.call(
    blit::exec,
    as.list(duckvep_blit_quote(c(
      file.path(R.home("bin"), "Rscript"),
      runner,
      arguments
    )))
  )
  command <- blit::cmd_wd(command, root)
  status <- suppressWarnings(blit::cmd_run(
    command,
    stdout = staged_outputs[["log"]],
    stderr = "2>&1",
    stdin = NULL,
    verbose = FALSE
  ))
  if (status != 0L) {
    stop(
      "DuckVEP campaign ", campaign$id, " failed; retained log: ",
      staged_outputs[["log"]],
      call. = FALSE
    )
  }
  missing <- staged_outputs[!file.exists(staged_outputs)]
  if (length(missing)) {
    stop(
      "DuckVEP campaign ", campaign$id, " omitted output(s): ",
      paste(names(missing), collapse = ", "),
      "; retained run directory: ", stage_directory,
      call. = FALSE
    )
  }
  duckvep_targets_publish_outputs(stage_directory, directory)
  outputs
}
