# Counts from complete HTSlib VCF text, independent of either DuckHTS reader.
# This workload has one sample and VCF 4.2 GT spelling. Missing list elements
# occupy slots; an absent list has no child slots. Scalars occupy one per call.
genotype_format_counts <- function(lines) {
  totals <- c(records = length(lines), calls = length(lines), gt_slots = 0, ps_values = 0,
              ad_slots = 0, ad_values = 0, dp_values = 0, gq_values = 0)
  if (!length(lines)) return(totals)
  rows <- strsplit(lines, "\t", fixed = TRUE)
  stopifnot(all(lengths(rows) == 10L))
  formats <- vapply(rows, `[[`, character(1), 9L)
  for (format in unique(formats)) {
    members <- strsplit(format, ":", fixed = TRUE)[[1L]]
    calls <- strsplit(vapply(rows[formats == format], `[[`, character(1), 10L), ":", fixed = TRUE)
    field <- function(name) {
      index <- match(name, members)
      if (is.na(index)) return(rep(NA_character_, length(calls)))
      vapply(calls, function(call) if (length(call) >= index) call[[index]] else ".", character(1))
    }
    gt <- field("GT")
    present <- !is.na(gt)
    stopifnot(all(grepl("^[0-9.]+([|/][0-9.]+)*$", gt[present])))
    totals["gt_slots"] <- totals["gt_slots"] + sum(nchar(gsub("[^|/]", "", gt[present]))) + sum(present)
    ps <- field("PS")
    totals["ps_values"] <- totals["ps_values"] + sum(!is.na(ps) & ps != ".")
    ad <- field("AD")
    values <- strsplit(ad[!is.na(ad)], ",", fixed = TRUE)
    totals["ad_slots"] <- totals["ad_slots"] + sum(lengths(values))
    totals["ad_values"] <- totals["ad_values"] + sum(unlist(values, use.names = FALSE) != ".")
    for (name in c("DP", "GQ")) {
      values <- field(name)
      key <- paste0(tolower(name), "_values")
      totals[key] <- totals[key] + sum(!is.na(values) & values != ".")
    }
  }
  totals
}

genotype_format_oracle <- function(input) {
  output <- tempfile("genotype-format-oracle-", fileext = ".vcf")
  on.exit(unlink(output))
  stopifnot(system2("bcftools", shQuote(c("view", "-H", input)), stdout = output) == 0L)
  connection <- file(output, "r")
  on.exit(close(connection), add = TRUE, after = FALSE)
  counts <- genotype_format_counts(character())
  repeat {
    lines <- readLines(connection, n = 65536L)
    if (!length(lines)) break
    counts <- counts + genotype_format_counts(lines)
  }
  counts
}

# One CTAS in a fresh process. Snapshots/checksums/denominators are outside timing.
genotype_format_run <- function(extension, input, fields, calls_projected, snapshot = "") {
  directory <- tempfile("genotype-format-run-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  database <- file.path(directory, "result.duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")),
                        dbdir = database)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  quote <- function(text) as.character(DBI::dbQuoteString(con, text))
  query <- function(sql) DBI::dbGetQuery(con, sql)
  DBI::dbExecute(con, paste("LOAD", quote(extension)))
  DBI::dbExecute(con, "SET threads=1")
  DBI::dbExecute(con, "SET memory_limit='8GB'")
  selection <- paste(vapply(fields, quote, character(1)), collapse = ",")
  columns <- if (calls_projected) "*" else "record_index, CHROM, POS, ID, REF, ALT"
  sql <- sprintf(paste("CREATE TABLE result AS SELECT %s FROM read_geno(%s,",
    "format_fields := [%s], scan_mode := 'sequential', decompression_threads := 0,",
    "decode_error_policy := 'error')"), columns, quote(input), selection)
  elapsed <- system.time(DBI::dbExecute(con, sql))
  status <- readLines("/proc/self/status")
  rss <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", status[grepl("^VmHWM:", status)]))
  counts <- c(records = query("SELECT count(*) AS n FROM result")$n, calls = 0, gt_slots = 0,
              ps_values = 0, ad_slots = 0, ad_values = 0, dp_values = 0, gq_values = 0)
  ordinals <- query("SELECT min(record_index) AS first, max(record_index) AS last, count(DISTINCT record_index) AS n FROM result")
  stopifnot(ordinals$first == 0, ordinals$last == counts["records"] - 1,
            ordinals$n == counts["records"])
  if (calls_projected) {
    DBI::dbExecute(con, "CREATE VIEW expanded AS SELECT unnest(calls) AS c FROM result")
    common <- query(paste("SELECT count(*) AS calls, sum(len(c.alleles)) AS gt_slots,",
                          "count(c.phase_set) AS ps_values FROM expanded"))
    counts[names(common)] <- unlist(common, use.names = FALSE)
    if ("AD" %in% fields) {
      ad <- query("SELECT coalesce(sum(len(c.format.AD)),0) AS ad_slots, coalesce(sum(list_count(c.format.AD)),0) AS ad_values FROM expanded")
      counts[names(ad)] <- unlist(ad, use.names = FALSE)
    }
    for (field in intersect(c("DP", "GQ"), fields))
      counts[paste0(tolower(field), "_values")] <- query(sprintf(
        "SELECT count(c.format.%s) AS n FROM expanded", field))$n
  }
  checksum <- query("SELECT sum(hash(r)::HUGEINT)::VARCHAR AS hash FROM result r")$hash
  DBI::dbExecute(con, "CHECKPOINT")
  bytes <- file.info(database)$size
  if (nzchar(snapshot)) DBI::dbExecute(con, sprintf(
    "COPY result TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", quote(snapshot)))
  data.frame(as.list(counts), extra_slots = counts["ad_slots"] +
    counts["calls"] * sum(c("DP", "GQ") %in% fields),
    elapsed = unname(elapsed["elapsed"]), cpu = sum(elapsed[c("user.self", "sys.self")]),
    peak_rss_kib = rss, database_bytes = bytes, checksum = checksum, row.names = NULL)
}

genotype_format_difference <- function(con, left, right) {
  sql <- sprintf(paste("SELECT count(*) AS n FROM (((%s) EXCEPT ALL (%s))",
                       "UNION ALL ((%s) EXCEPT ALL (%s)))"), left, right, right, left)
  as.numeric(DBI::dbGetQuery(con, sql)$n)
}

# Read the exact rendered cohort table; require identical registered inputs and
# denominator tables before making a revision comparison.
genotype_hprc_comparison <- function(baseline, current) {
  reports <- lapply(c(baseline, current), readLines)
  identity <- function(lines) lines[grepl("^\\| (VCF|BCF) +\\| geno_hprc", lines)]
  denominators <- function(lines) lines[grepl("^\\| (full|selected|sparse|carriers) +\\|", lines)]
  stopifnot(length(identity(reports[[1]])) == 2L,
            identical(identity(reports[[1]]), identity(reports[[2]])),
            length(denominators(reports[[1]])) == 4L,
            identical(denominators(reports[[1]]), denominators(reports[[2]])))
  rows <- lapply(reports, function(lines) {
    first <- grep("^\\| format +\\| reader +\\| workload +\\| elapsed +\\|", lines)
    stopifnot(length(first) == 1L)
    last <- first + 2L
    while (last <= length(lines) && startsWith(lines[last], "|")) last <- last + 1L
    split <- function(x) trimws(strsplit(substring(x, 2, nchar(x) - 1), "|", fixed=TRUE)[[1L]])
    table <- as.data.frame(do.call(rbind,lapply(lines[seq.int(first+2L,last-1L)],split)))
    names(table) <- split(lines[first])
    stopifnot(nrow(table) == 16L)
    table[c("format","reader","workload","elapsed","peak_rss_kib")]
  })
  result <- merge(rows[[1]],rows[[2]],by=c("format","reader","workload"),suffixes=c("_baseline","_current"))
  stopifnot(nrow(result) == 16L, !anyDuplicated(result[c("format","reader","workload")]))
  for (field in setdiff(names(result),c("format","reader","workload")))
    result[[field]] <- as.numeric(result[[field]])
  stopifnot(!anyNA(result), all(result$elapsed_baseline > 0))
  result$elapsed_change_percent <- 100 * (result$elapsed_current / result$elapsed_baseline - 1)
  result
}
