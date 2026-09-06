# One measured materialization per fresh process. The report supplies registry
# paths and independent denominators; this driver never stages or downloads.
genotype_reader_run <- function(extension, input, reader, workload, selector,
                                sample_names, snapshot, threads = 1L) {
  directory <- tempfile("duckhts-genotype-run-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  database <- file.path(directory, "result.duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")), dbdir = database)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  quote <- function(text) as.character(DBI::dbQuoteString(con, text))
  query <- function(sql) DBI::dbGetQuery(con, sql)
  DBI::dbExecute(con, paste("LOAD", quote(extension)))
  DBI::dbExecute(con, paste("SET threads=", threads))
  DBI::dbExecute(con, "SET memory_limit='8GB'")
  DBI::dbWriteTable(con, "samples", data.frame(sample_index = seq_along(sample_names) - 1L,
                                               sample_name = sample_names))
  arguments <- paste0(quote(input), ", scan_mode:='sequential', decompression_threads:=0",
                      if (nzchar(selector)) paste0(", samples:=", quote(selector)) else "")
  sparse <- workload == "sparse"
  if (reader == "read_geno") {
    source <- paste0("SELECT * FROM read_geno(", arguments,
                     ", non_reference_only:=", if (sparse) "true" else "false", ")")
    normalized <- function(table) paste0(
      "SELECT CHROM, POS, ID, REF, ALT, c.sample_index, c.alleles, c.phase_before, c.phase_set ",
      "FROM (SELECT CHROM,POS,ID,REF,ALT,unnest(calls) AS c FROM ", table, ")")
  } else {
    source <- paste0("SELECT CHROM,POS,ID,REF,ALT,SAMPLE_ID,FORMAT_GT FROM read_bcf(",
                     arguments, ", tidy_format:=true)",
                     if (sparse) " WHERE regexp_matches(FORMAT_GT,'(^|[|/])[1-9][0-9]*($|[|/])')" else "")
    # The registered source is VCF 4.2, GT-only. HTSlib's leading-slot flag is
    # implied by the other separators (or by a known haploid); no PS is invented.
    normalized <- function(table) paste0(
      "SELECT CHROM,POS,ID,REF,ALT,s.sample_index,",
      "list_transform(tokens,x->CASE WHEN x='.' THEN NULL ELSE x::INTEGER END) AS alleles,",
      "list_transform(range(1,len(tokens)+1),i->CASE WHEN i=1 THEN CASE WHEN len(tokens)=1 ",
      "THEN tokens[1]<>'.' ELSE NOT contains(FORMAT_GT,'/') END ELSE separators[i-1]='|' END) AS phase_before,",
      "NULL::BIGINT AS phase_set FROM (SELECT *,regexp_split_to_array(FORMAT_GT,'[|/]') AS tokens,",
      "regexp_extract_all(FORMAT_GT,'[|/]') AS separators FROM ", table, ") q JOIN samples s ON SAMPLE_ID=s.sample_name")
  }
  carriers <- function(calls) paste0(
    "SELECT CHROM,POS,ID,REF,ALT[allele_index] AS ALT,sample_index,slot,allele_index,phase_before,phase_set FROM (",
    "SELECT CHROM,POS,ID,REF,ALT,sample_index,generate_subscripts(alleles,1)-1 AS slot,",
    "unnest(alleles) AS allele_index,unnest(phase_before) AS phase_before,phase_set FROM (", calls,
    ")) WHERE allele_index>0")
  sql <- if (workload == "carriers") carriers(normalized(paste0("(", source, ")"))) else source
  elapsed <- system.time(DBI::dbExecute(con, paste("CREATE TABLE result AS", sql)))
  status <- readLines("/proc/self/status")
  rss <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*", "\\1", status[grepl("^VmHWM:", status)]))
  output_rows <- as.numeric(query("SELECT count(*) AS n FROM result")$n)
  fingerprint <- query("SELECT sum(hash(r)::HUGEINT)::VARCHAR AS hash FROM result r")$hash
  # Storage after checkpoint is an observed materialization denominator, not
  # an estimate of vector capacity or total first-party allocation.
  DBI::dbExecute(con, "CHECKPOINT")
  materialized_bytes <- file.info(database)$size
  DBI::dbExecute(con, paste("CREATE VIEW normalized AS", if (workload == "carriers")
    "SELECT * FROM result" else normalized("result")))
  if (nzchar(snapshot)) DBI::dbExecute(con, paste0(
    "COPY normalized TO ", quote(snapshot), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  counts <- if (workload == "carriers") data.frame(calls = NA_real_, slots = output_rows) else query(
    "SELECT count(*) AS calls,sum(len(alleles)) AS slots FROM normalized")
  data.frame(reader, workload, output_rows, output_calls = as.numeric(counts$calls),
             output_slots = as.numeric(counts$slots), checksum = fingerprint,
             elapsed = unname(elapsed["elapsed"]), cpu = sum(elapsed[c("user.self", "sys.self")]),
             peak_rss_kib = rss, materialized_database_bytes = materialized_bytes)
}
