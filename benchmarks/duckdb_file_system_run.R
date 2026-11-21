# read_bam through each file-access route: htslib's own file access (libcurl
# for http) and the DuckDB file system of the connection (httpfs for http).
# Source from benchmarks/benchmark_duckdb_file_system.Rmd at the repository
# root. Every measurement is one fresh R process.

file_system_routes <- c("htslib", "duckdb")

# htslib tuning that DuckHTS applies to a remote alignment handle, from
# src/include/hts_io_tuning.h (native build). Both routes receive it: the
# profile depends on the path the user gave, not on who reads it.
file_system_tuning <- data.frame(
  profile = c("indexed region or parallel scan", "single-handle streaming scan"),
  hfile_block_mib = c(1, 4),
  bgzf_cache_mib = c(64, 256)
)

# Largest single request the DuckDB route issued while this benchmark was
# designed (httpfs grows its read-ahead to this size). The peak-RSS ceilings
# below budget one such buffer per open handle.
file_system_read_ahead_mib <- 32

# Peak-RSS ceiling in MiB: the idle process plus three times the live I/O state
# of `handles` open alignment handles under one tuning profile.
file_system_rss_ceiling_mib <- function(idle_mib, profile, handles) {
  tuning <- file_system_tuning[file_system_tuning$profile == profile, ]
  stopifnot(nrow(tuning) == 1L)
  idle_mib + 3 * handles *
    (tuning$hfile_block_mib + tuning$bgzf_cache_mib + file_system_read_ahead_mib)
}

file_system_start_server <- function(directory, latency_ms, python = "python3") {
  log <- tempfile("range_server_", fileext = ".tsv")
  file.create(log)
  process <- processx::process$new(
    python, c("scripts/benchmark_range_server.py", directory, log, as.character(latency_ms)),
    stdout = "|", stderr = "|"
  )
  process$poll_io(20000)
  ready <- process$read_output_lines()
  if (length(ready) != 1L || !grepl("^ready [0-9]+$", ready)) {
    process$kill()
    stop("the range server did not start: ", paste(c(ready, process$read_error_lines()), collapse = "\n"))
  }
  list(process = process, log = log, latency_ms = latency_ms,
       base = paste0("http://127.0.0.1:", sub("^ready ", "", ready)))
}

# Requests the server logged since the log was last emptied.
file_system_requests <- function(log) {
  lines <- readLines(log)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  stopifnot(all(lengths(fields) == 5L))
  method <- vapply(fields, `[[`, "", 1L)
  range <- vapply(fields, `[[`, "", 3L)
  written <- as.numeric(vapply(fields, `[[`, "", 5L))
  bounds <- regmatches(range, regexec("^bytes=([0-9]+)-([0-9]*)$", range))
  ranged <- lengths(bounds) == 3L
  last <- vapply(bounds[ranged], `[[`, "", 3L)
  first <- as.numeric(vapply(bounds[ranged], `[[`, "", 2L))
  closed <- nzchar(last)
  closed_bytes <- as.numeric(last[closed]) - first[closed] + 1
  data.frame(
    requests = length(lines),
    get_requests = sum(method == "GET"),
    head_requests = sum(method == "HEAD"),
    open_ended_ranges = sum(!closed),
    closed_ranges = sum(closed),
    median_closed_range_bytes = if (any(closed)) stats::median(closed_bytes) else NA_real_,
    max_closed_range_bytes = if (any(closed)) max(closed_bytes) else NA_real_,
    body_bytes_written = sum(written)
  )
}

# httpfs settings measured on the DuckDB route, as the SET statements that
# select them. Each differs from the defaults in one respect.
file_system_httpfs_settings <- list(
  defaults = character(),
  force_download = "SET force_download = true",
  no_external_file_cache = "SET enable_external_file_cache = false",
  grid_request_sizing = "SET external_file_cache_request_sizing = 'GRID'",
  cache_block_16mib = "SET external_file_cache_remote_max_block_size = 16777216",
  grid_and_cache_block_16mib = c("SET external_file_cache_request_sizing = 'GRID'",
                                 "SET external_file_cache_remote_max_block_size = 16777216"),
  no_keep_alive = "SET http_keep_alive = false",
  no_connection_caching = "SET httpfs_connection_caching = false",
  curl_client = "SET httpfs_client_implementation = 'curl'"
)

# One query in a fresh process. `region = NA` scans the whole file; `idle`
# loads everything and reads nothing, which gives the fixed process overhead.
# `settings` are SET statements run after the route is chosen.
file_system_measure <- function(extension, source, region, route, threads, idle = FALSE,
                                settings = character()) {
  callr::r(function(extension, source, region, route, threads, idle, settings) {
    con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
    quote <- function(text) as.character(DBI::dbQuoteString(con, text))
    DBI::dbExecute(con, paste("LOAD", quote(extension)))
    # Both routes load httpfs for an http source, so that neither timed query
    # pays for loading it.
    if (startsWith(source, "http://")) DBI::dbExecute(con, "LOAD httpfs")
    DBI::dbExecute(con, sprintf("SET threads = %d", threads))
    DBI::dbExecute(con, "SET memory_limit = '2GiB'")
    DBI::dbExecute(con, "SET max_temp_directory_size = '0B'")
    DBI::dbExecute(con, paste("SET duckhts_file_system =", quote(route)))
    for (statement in settings) DBI::dbExecute(con, statement)
    arguments <- if (is.na(region)) quote(source) else paste0(quote(source), ", region := ", quote(region))
    sql <- if (idle) "SELECT 0::BIGINT AS records, '' AS digest" else paste0(
      "SELECT count(*) AS records, CAST(sum(hash(QNAME, FLAG, POS, CIGAR)) AS VARCHAR) AS digest ",
      "FROM read_bam(", arguments, ")"
    )
    seconds <- system.time(result <- DBI::dbGetQuery(con, sql))[["elapsed"]]
    memory <- DBI::dbGetQuery(con, paste(
      "SELECT coalesce(sum(memory_usage_bytes), 0) AS buffer_bytes,",
      "coalesce(sum(temporary_storage_bytes), 0) AS temporary_bytes FROM duckdb_memory()"
    ))
    status <- readLines("/proc/self/status")
    peak <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", status[grepl("^VmHWM:", status)]))
    data.frame(records = as.numeric(result$records), digest = result$digest, seconds = seconds,
               peak_rss_mib = peak / 1024,
               duckdb_buffer_bytes_after = as.numeric(memory$buffer_bytes),
               duckdb_temporary_bytes_after = as.numeric(memory$temporary_bytes))
  }, list(extension, source, region, route, threads, idle, settings))
}

# The 100 kb and 10 Mb regions of `bam` on the DuckDB route under each httpfs
# setting, at each server delay, with the htslib route at its defaults as the
# reference. One thread. Returns the oracle and one row per measurement.
file_system_settings_benchmark <- function(extension, bam, repetitions = 3L,
                                           latencies_ms = c(0, 20), python = "python3") {
  extension <- normalizePath(extension, mustWork = TRUE)
  bam <- normalizePath(bam, mustWork = TRUE)
  stopifnot(file.exists(paste0(bam, ".bai")), nzchar(Sys.which("samtools")), 0 %in% latencies_ms)
  served <- tempfile("served_")
  dir.create(served)
  stopifnot(all(file.symlink(c(bam, paste0(bam, ".bai")), served)))
  servers <- lapply(stats::setNames(latencies_ms, latencies_ms), function(latency) {
    file_system_start_server(served, latency, python)
  })
  on.exit({
    for (server in servers) server$process$kill()
    unlink(served, recursive = TRUE)
  })

  regions <- c("chr20:10000000-10100000", "chr20:10000000-20000000")
  source <- paste0(servers[["0"]]$base, "/", basename(bam))
  idle_mib <- file_system_measure(extension, source, NA, "duckdb", 1L, idle = TRUE)$peak_rss_mib
  handle_ceiling <- file_system_rss_ceiling_mib(idle_mib, "indexed region or parallel scan", 1)
  arms <- c(list(htslib = character()), file_system_httpfs_settings)
  plan <- expand.grid(setting = names(arms), region = regions, latency_ms = latencies_ms,
                      stringsAsFactors = FALSE)
  plan$route <- ifelse(plan$setting == "htslib", "htslib", "duckdb")
  # A forced download holds the whole file.
  plan$rss_ceiling_mib <- handle_ceiling +
    ifelse(plan$setting == "force_download", 3 * file.size(bam) / 2^20, 0)

  oracle <- do.call(rbind, lapply(regions, function(region) {
    cbind(region = region, file_system_oracle(bam, region))
  }))
  runs <- list()
  for (repetition in seq_len(repetitions)) {
    for (row in seq_len(nrow(plan))) {
      arm <- plan[row, ]
      server <- servers[[as.character(arm$latency_ms)]]
      file.create(server$log)
      measured <- file_system_measure(extension, paste0(server$base, "/", basename(bam)), arm$region,
                                      arm$route, 1L, settings = arms[[arm$setting]])
      runs[[length(runs) + 1L]] <- cbind(arm, repetition = repetition, measured,
                                         file_system_requests(server$log), row.names = NULL)
    }
  }
  runs <- do.call(rbind, runs)
  expected <- oracle$records[match(runs$region, oracle$region)]
  stopifnot(all(runs$records == expected),
            all(tapply(runs$digest, runs$region, function(digest) length(unique(digest))) == 1L),
            all(runs$duckdb_temporary_bytes_after == 0),
            all(runs$peak_rss_mib <= runs$rss_ceiling_mib))
  list(oracle = oracle, runs = runs)
}

# Median time and request pattern of each setting, region and delay.
file_system_settings_summary <- function(runs) {
  key <- paste(runs$setting, runs$region, runs$latency_ms, sep = "|")
  rows <- lapply(split(runs, factor(key, levels = unique(key))), function(group) {
    data.frame(
      region = group$region[[1L]], latency_ms = group$latency_ms[[1L]], setting = group$setting[[1L]],
      runs = nrow(group), median_seconds = stats::median(group$seconds),
      min_seconds = min(group$seconds), max_seconds = max(group$seconds),
      max_peak_rss_mib = max(group$peak_rss_mib),
      get_requests = stats::median(group$get_requests),
      max_closed_range_mib = max(group$max_closed_range_bytes) / 2^20,
      body_mib_written = stats::median(group$body_bytes_written) / 2^20
    )
  })
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$region, summary$latency_ms), ]
  rownames(summary) <- NULL
  summary
}

# Record count of a file or region from samtools, the independent oracle.
file_system_oracle <- function(bam, region) {
  arguments <- c("view", "-c", "-@", "4", shQuote(bam), if (!is.na(region)) shQuote(region))
  seconds <- system.time(count <- system2("samtools", arguments, stdout = TRUE))[["elapsed"]]
  stopifnot(length(count) == 1L, grepl("^[0-9]+$", count))
  data.frame(records = as.numeric(count), seconds = seconds)
}

# Runs `plan` (one row per workload) `repetitions` times, alternating routes
# inside each repetition, and returns one row per measurement.
file_system_run_plan <- function(extension, plan, repetitions, servers) {
  runs <- list()
  for (row in seq_len(nrow(plan))) {
    workload <- plan[row, ]
    server <- if (is.na(workload$latency_ms)) NULL else servers[[as.character(workload$latency_ms)]]
    source <- if (is.null(server)) workload$path else paste0(server$base, "/", basename(workload$path))
    for (repetition in seq_len(repetitions)) {
      for (route in file_system_routes) {
        if (!is.null(server)) file.create(server$log)
        measured <- file_system_measure(extension, source, workload$region, route, workload$threads)
        requests <- if (is.null(server)) file_system_requests(nullfile()) else file_system_requests(server$log)
        runs[[length(runs) + 1L]] <- cbind(workload, route = route, repetition = repetition,
                                           measured, requests, row.names = NULL)
      }
    }
  }
  do.call(rbind, runs)
}

# `bams` names the 1x, 2x and 4x alignment files, each with a .bai beside it.
# Returns the idle overhead, the samtools oracle, and every measurement.
file_system_benchmark <- function(extension, bams, region_repetitions = 5L, scan_repetitions = 3L,
                                  latencies_ms = c(0, 20), scan_threads = c(1L, 4L),
                                  python = "python3") {
  extension <- normalizePath(extension, mustWork = TRUE)
  bams <- vapply(bams, normalizePath, "", mustWork = TRUE)
  stopifnot(identical(names(bams), c("1x", "2x", "4x")), all(file.exists(paste0(bams, ".bai"))),
            nzchar(Sys.which("samtools")), 0 %in% latencies_ms)

  served <- tempfile("served_")
  dir.create(served)
  stopifnot(all(file.symlink(c(bams, paste0(bams, ".bai")), served)))
  servers <- lapply(stats::setNames(latencies_ms, latencies_ms), function(latency) {
    file_system_start_server(served, latency, python)
  })
  on.exit({
    for (server in servers) server$process$kill()
    unlink(served, recursive = TRUE)
  })

  start <- 10000000
  widths <- c(100000, 2500000, 5000000, 10000000)
  regions <- sprintf("chr20:%.0f-%.0f", start, start + widths)
  region_plan <- expand.grid(
    workload = "region", scale = "1x", region = regions, latency_ms = latencies_ms, threads = 1L,
    stringsAsFactors = FALSE
  )
  scan_plan <- rbind(
    expand.grid(workload = "scan", scale = names(bams), region = NA_character_, latency_ms = 0,
                threads = scan_threads, stringsAsFactors = FALSE),
    data.frame(workload = "scan", scale = "1x", region = NA_character_, latency_ms = NA_real_,
               threads = 1L)
  )
  region_plan$path <- bams[region_plan$scale]
  scan_plan$path <- bams[scan_plan$scale]

  idle <- do.call(rbind, lapply(file_system_routes, function(route) {
    source <- paste0(servers[["0"]]$base, "/", basename(bams[["1x"]]))
    cbind(route = route, file_system_measure(extension, source, NA, route, max(scan_threads), idle = TRUE))
  }))
  idle_mib <- max(idle$peak_rss_mib)
  ceiling_for <- function(plan) {
    parallel <- plan$workload == "region" | plan$threads > 1L
    ifelse(parallel,
           file_system_rss_ceiling_mib(idle_mib, "indexed region or parallel scan", plan$threads),
           file_system_rss_ceiling_mib(idle_mib, "single-handle streaming scan", 1))
  }
  region_plan$rss_ceiling_mib <- ceiling_for(region_plan)
  scan_plan$rss_ceiling_mib <- ceiling_for(scan_plan)

  targets <- unique(rbind(region_plan, scan_plan)[c("scale", "region", "path")])
  oracle <- do.call(rbind, lapply(seq_len(nrow(targets)), function(row) {
    cbind(targets[row, c("scale", "region")],
          file_system_oracle(targets$path[[row]], targets$region[[row]]), row.names = NULL)
  }))

  runs <- rbind(
    file_system_run_plan(extension, region_plan, region_repetitions, servers),
    file_system_run_plan(extension, scan_plan, scan_repetitions, servers)
  )
  runs$path <- NULL
  list(idle = idle, oracle = oracle, runs = runs)
}

# Stops unless every measurement returned the oracle's record count, one digest
# per input, no spill, and a peak RSS within its ceiling.
file_system_check <- function(result) {
  runs <- result$runs
  key <- function(frame) paste(frame$scale, frame$region)
  expected <- result$oracle$records[match(key(runs), key(result$oracle))]
  stopifnot(!anyNA(expected), all(runs$records == expected),
            all(tapply(runs$digest, key(runs), function(digest) length(unique(digest))) == 1L),
            all(runs$duckdb_temporary_bytes_after == 0),
            all(runs$peak_rss_mib <= runs$rss_ceiling_mib))
  invisible(result)
}

# Median and range of the timings and the request pattern of each workload.
file_system_summary <- function(runs) {
  groups <- c("workload", "scale", "region", "latency_ms", "threads", "route")
  key <- do.call(paste, c(lapply(runs[groups], as.character), sep = "|"))
  rows <- lapply(split(runs, factor(key, levels = unique(key))), function(group) {
    cbind(group[1L, groups], data.frame(
      runs = nrow(group), records = group$records[[1L]],
      median_seconds = stats::median(group$seconds),
      min_seconds = min(group$seconds), max_seconds = max(group$seconds),
      median_peak_rss_mib = stats::median(group$peak_rss_mib),
      max_peak_rss_mib = max(group$peak_rss_mib),
      rss_ceiling_mib = group$rss_ceiling_mib[[1L]],
      get_requests = stats::median(group$get_requests),
      open_ended_ranges = stats::median(group$open_ended_ranges),
      closed_ranges = stats::median(group$closed_ranges),
      max_closed_range_mib = max(group$max_closed_range_bytes) / 2^20,
      body_mib_written = stats::median(group$body_bytes_written) / 2^20
    ), row.names = NULL)
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary
}
