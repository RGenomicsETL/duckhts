#!/usr/bin/env Rscript
# Seeded exact-component expansion and evidence-status checks against base R.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(optparse) })
op <- OptionParser()
op <- add_option(op, '--extension', default = 'build/release/duckhts.duckdb_extension')
op <- add_option(op, '--trials', type = 'integer', default = 100000L)
op <- add_option(op, '--seed', type = 'integer', default = 173L)
op <- add_option(op, '--out', default = '')
opt <- parse_args(op)
stopifnot(!is.na(opt$trials), opt$trials >= 7L, !is.na(opt$seed))

check_comparisons <- function(x, ids) {
  stopifnot(nrow(x) == length(ids), !anyDuplicated(x$scene), identical(x$scene, ids))
  same_sequence <- (is.na(x$actual_sequence) & is.na(x$expected_sequence)) |
    (!is.na(x$actual_sequence) & !is.na(x$expected_sequence) &
       x$actual_sequence == x$expected_sequence)
  stopifnot(all(same_sequence), identical(x$actual_status, x$expected_status))
  TRUE
}

run <- function() {
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  out <- if (nzchar(opt$out)) opt$out else tempfile(paste0('repeat_sequence_seed', opt$seed, '_'),
    'test/duckvep/conformance/results')
  stopifnot(!dir.exists(out), dir.create(out, recursive = TRUE))
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = 'true')))
  on.exit(dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(dbQuoteString(con, x))
  dbExecute(con, paste('LOAD', q(extension)))
  dbExecute(con, 'SET threads=1')
  set.seed(opt$seed)
  alphabet <- strsplit('ACGTRYSWKMBDHVNacgtryswkmbdhvn', '', fixed = TRUE)[[1L]]
  coverage <- setNames(integer(7L), c('exact', 'summary', 'missing_count', 'missing_unit',
                                    'fractional', 'empty', 'missing_components'))
  controls <- NULL
  capacity_checks <- list()
  for (first in seq.int(1L, opt$trials, by = 1000L)) {
    ids <- seq.int(first, min(first + 999L, opt$trials))
    cases <- data.frame(scene = ids, sequence_exact = TRUE, capacity = 5000,
      has_components = TRUE, expected_sequence = NA_character_, expected_status = 'ok')
    components <- vector('list', length(ids))
    for (i in seq_along(ids)) {
      state <- (ids[i] - 1L) %% length(coverage) + 1L
      coverage[state] <- coverage[state] + 1L
      n <- sample.int(6L, 1L)
      units <- vapply(seq_len(n), function(j)
        paste0(sample(alphabet, sample.int(9L, 1L), replace = TRUE), collapse = ''), '')
      counts <- as.double(sample(c(0:3, 10, 31, 100), n, replace = TRUE))
      if (state == 1L) {
        cases$expected_sequence[i] <- paste0(strrep(units, counts), collapse = '')
        cases$capacity[i] <- nchar(cases$expected_sequence[i]) + sample(0:2, 1L)
      } else if (state == 2L) {
        cases$sequence_exact[i] <- FALSE
        cases$expected_status[i] <- 'summary_only'
        cases$capacity[i] <- 0
        counts[sample.int(n, 1L)] <- 1e300
      } else if (state == 3L) {
        counts[sample.int(n, 1L)] <- NA_real_
        cases$expected_status[i] <- 'incomplete_input'
      } else if (state == 4L) {
        units[sample.int(n, 1L)] <- NA_character_
        cases$expected_status[i] <- 'incomplete_input'
      } else if (state == 5L) {
        at <- sample.int(n, 1L)
        counts[at] <- counts[at] + 0.5
        cases$expected_status[i] <- 'nonintegral_count'
      } else {
        units <- character(); counts <- numeric(); n <- 0L
        if (state == 6L) {
          cases$expected_sequence[i] <- ''
          cases$capacity[i] <- 0
        } else {
          cases$has_components[i] <- FALSE
          cases$expected_status[i] <- 'incomplete_input'
        }
      }
      components[[i]] <- data.frame(scene = rep(ids[i], n), ordinal = seq_len(n),
                                     unit = units, unit_count = counts)
    }
    components <- do.call(rbind, components)
    dbWriteTable(con, 'repeat_cases', cases, overwrite = TRUE)
    dbWriteTable(con, 'repeat_components', components, overwrite = TRUE)
    query <- paste(
      'WITH grouped AS (SELECT scene,list(struct_pack(unit:=unit,count:=unit_count)',
      'ORDER BY ordinal) parts FROM repeat_components GROUP BY scene), evaluated AS (',
      'SELECT c.*,duckvep_repeat_sequence(CASE WHEN c.has_components THEN',
      'coalesce(g.parts,[]::STRUCT(unit VARCHAR,count DOUBLE)[]) ELSE NULL END,',
      'c.sequence_exact,max_sequence_bases:=c.capacity) r',
      'FROM repeat_cases c LEFT JOIN grouped g USING(scene))',
      'SELECT scene,expected_sequence,expected_status,r.sequence AS actual_sequence,',
      'r.status AS actual_status FROM evaluated ORDER BY scene')
    actual <- try(dbGetQuery(con, query), silent = TRUE)
    if (inherits(actual, 'try-error')) {
      saveRDS(list(seed = opt$seed, first = first, cases = cases, components = components,
        error = as.character(actual)), file.path(out, 'counterexample.rds'))
      stop(actual)
    }
    dbWriteTable(con, 'cases_all', cases, append = first != 1L)
    dbWriteTable(con, 'components_all', components, append = first != 1L)
    dbWriteTable(con, 'comparisons', actual, append = first != 1L)
    ok <- try(check_comparisons(actual, ids), silent = TRUE)
    if (inherits(ok, 'try-error')) {
      saveRDS(list(seed = opt$seed, first = first, cases = cases, components = components,
        comparisons = actual, error = as.character(ok)), file.path(out, 'counterexample.rds'))
      stop(ok)
    }
    if (is.null(controls)) {
      duplicate <- actual; duplicate$scene[2L] <- duplicate$scene[1L]
      altered <- actual; altered$actual_sequence[1L] <- 'DELIBERATE_CORRUPTION'
      status <- actual; status$actual_status[2L] <- 'ok'
      controls <- vapply(list(drop = actual[-1L, ], duplicate = duplicate,
        altered_sequence = altered, altered_status = status), function(x)
        inherits(try(check_comparisons(x, ids), silent = TRUE), 'try-error'), TRUE)
      stopifnot(all(controls))
      eligible <- which(actual$expected_status == 'ok' & nchar(actual$expected_sequence) > 0L)
      for (i in head(eligible, 32L)) {
        limit <- nchar(actual$expected_sequence[i]) - 1L
        dbExecute(con, paste('UPDATE repeat_cases SET capacity=', limit, 'WHERE scene=', ids[i]))
        failure <- try(dbGetQuery(con, query), silent = TRUE)
        dbExecute(con, paste('UPDATE repeat_cases SET capacity=', cases$capacity[i],
                             'WHERE scene=', ids[i]))
        passed <- inherits(failure, 'try-error') &&
          grepl('exceeds max_sequence_bases=', as.character(failure), fixed = TRUE)
        capacity_checks[[length(capacity_checks) + 1L]] <- data.frame(scene = ids[i],
          capacity = limit, required = nchar(actual$expected_sequence[i]), passed = passed)
        if (!passed) {
          saveRDS(list(cases = cases, components = components, scene = ids[i],
            capacity = limit, result = failure), file.path(out, 'capacity_counterexample.rds'))
          stop('capacity exhaustion did not report the named limit')
        }
      }
      stopifnot(identical(dbGetQuery(con, query), actual))
    }
  }
  stopifnot(sum(coverage) == opt$trials, all(coverage >= opt$trials %/% 7L),
    dbGetQuery(con, 'SELECT count(*) n FROM comparisons')$n == opt$trials,
    dbGetQuery(con, 'SELECT count(DISTINCT scene) n FROM comparisons')$n == opt$trials)
  for (table in c('cases_all', 'components_all', 'comparisons'))
    dbExecute(con, paste('COPY', table, 'TO', q(file.path(out, paste0(table, '.parquet'))),
                        '(FORMAT PARQUET)'))
  write.csv(do.call(rbind, capacity_checks), file.path(out, 'capacity_checks.csv'), row.names = FALSE)
  sha <- function(path) digest::digest(file = path, algo = 'sha256', serialize = FALSE)
  sources <- c('test/duckvep/conformance/repeat_sequence_differential.R',
               'src/duckvep/duckvep_sql.c')
  for (path in sources) {
    destination <- file.path(out, 'source', path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, destination), sha(path) == sha(destination))
  }
  files <- list.files(out, recursive = TRUE, full.names = TRUE)
  jsonlite::write_json(list(scope = 'ordered exact repeat expansion and explicit evidence statuses',
    seed = opt$seed, trials = opt$trials, threads = 1L, batch_size = 1000L,
    strata = as.list(coverage), corruption_controls = as.list(controls),
    capacity_exhaustion_checks = length(capacity_checks), recovery_exact = TRUE,
    oracle = 'base R strrep over independently generated ordered components',
    R_version = R.version.string, RNG_kind = RNGkind(), duckdb_version = as.character(packageVersion('duckdb')),
    source_revision = system2('git', c('rev-parse', 'HEAD'), stdout = TRUE),
    source_dirty = length(system2('git', c('status', '--porcelain'), stdout = TRUE)) > 0L,
    build_binding = 'diagnostic_unbound', extension_sha256 = sha(extension),
    failed = 0L, skipped = 0L, failures_waived = 0L,
    interpretation = 'Generated component scenes and status strata, not independent biological samples or a VEP parser-conformance claim.',
    sha256 = as.list(setNames(vapply(files, sha, ''), substring(files, nchar(out) + 2L)))),
    file.path(out, 'receipt.json'), pretty = TRUE, auto_unbox = TRUE)
  cat(opt$trials, 'scenes passed; all seven strata and four corruption controls; ', out, '\n')
}
run()
