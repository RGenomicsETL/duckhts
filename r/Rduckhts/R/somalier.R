#' Import an Already Selected Somalier Sites VCF or BCF
#'
#' Convert an existing Somalier-compatible sites file into the canonical typed
#' panel and population-frequency relation used by DuckHTS extraction,
#' relatedness, and contamination functions. REF and ALT are oriented into
#' lexical A/B order and alternate-allele frequency is flipped with the alleles,
#' so `population_b_af` always describes `allele_b`. Exact Somalier v0.3.4 X/Y
#' aliases are excluded, matching its autosomal frequency importer.
#'
#' This function does not select sites from a population VCF. Somalier's
#' `find-sites` algorithm has separate AF/AN, QC, interval-exclusion, and spacing
#' semantics and is not implied by importing its output.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param path One already selected sites VCF/BCF path or URI. `INFO/AF` must be
#'   declared `Number=A,Type=Float` and every retained record must be one
#'   canonical biallelic SNV with one finite AF value.
#' @param assembly Nonempty assembly identifier attached to every site.
#' @param max_sites Positive input-site limit, at most 100,000,000.
#' @param table_name Optional output table. `NULL` returns a data frame.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_import_sites <- function(
  con, path, assembly, max_sites = 1000000,
  table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(path, "path")
  .somalier_scalar_text(assembly, "assembly")
  if (nchar(assembly, type = "bytes") > 1024L) {
    stop("assembly must be at most 1,024 bytes", call. = FALSE)
  }
  max_sites <- .somalier_bounded_whole_number(
    max_sites, "max_sites", 1, 100000000
  )
  query <- sprintf(
    paste0(
      "SELECT * FROM duckhts_somalier_import_sites(%s, %s, max_sites := %s)"
    ),
    sql_quote_string(con, path), sql_quote_string(con, assembly), max_sites
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Extract Panel-Aligned Counts from VCF or BCF
#'
#' Produce the complete sample-by-panel count relation consumed by the
#' Somalier-derived relatedness and contamination functions. `FORMAT/AD` must
#' declare `Number=R,Type=Integer`; A and B slots are matched by exact REF/ALT
#' identity, and `other` sums only the remaining declared-allele slots. Missing
#' sites and unavailable AD remain rows with three NULL counts, distinct from
#' measured zero depth. The panel can be any typed table/view or ordinary
#' Parquet file with the canonical six panel identity columns.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param path One VCF/BCF path or URI.
#' @param panel_table Name of the ordered panel table or view.
#' @param panel_parquet Ordinary panel Parquet path, instead of `panel_table`.
#' @param samples Optional HTSlib sample selector: comma-separated inclusion,
#'   leading `^` exclusion, `"-"` for all, or `""` for none.
#' @param filter_policy Record FILTER policy: `"pass_or_unapplied"` makes named
#'   failures unavailable, `"include_all"` uses their AD, and `"error"` rejects
#'   a selected panel record with a named failure.
#' @param table_name Optional output table. `NULL` returns a data frame.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_vcf_counts <- function(
  con, path, panel_table = NULL, panel_parquet = NULL, samples = NULL,
  filter_policy = c("pass_or_unapplied", "include_all", "error"),
  table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(path, "path")
  if (!is.null(samples)) .somalier_scalar_text(samples, "samples", allow_empty = TRUE)
  filter_policy <- match.arg(filter_policy)

  panel <- .somalier_source(con, panel_table, panel_parquet, "panel")
  if (panel$temporary) {
    on.exit(.somalier_drop_view(con, panel$name), add = TRUE, after = FALSE)
  }
  arguments <- c(
    sql_quote_string(con, path), sql_quote_string(con, panel$name),
    if (!is.null(samples)) paste0("samples := ", sql_quote_string(con, samples)),
    paste0("filter_policy := ", sql_quote_string(con, filter_policy))
  )
  query <- paste0(
    "SELECT * FROM duckhts_somalier_vcf_counts(",
    paste(arguments, collapse = ", "), ")"
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Extract Panel-Aligned Counts from BAM or CRAM
#'
#' Count observed A, B, and other query bases at every site in a typed panel.
#' The panel is prepared once, then each scan worker owns one indexed
#' multi-region scan and independent BAM/CRAM, index, reference, pileup, and
#' overlap state. `worker_count` bounds the number of DuckDB-scheduled panel
#' shards; the connection's thread setting bounds how many can run concurrently.
#' `decompression_threads` separately controls htslib decompression workers per
#' source handle.
#' Valid uncovered sites are measured zero depth; reference or alignment-header
#' mismatches remain rows with NULL counts and a named status. The panel can be
#' a committed table/view or an ordinary Parquet file. Caller-local temporary
#' relations and uncommitted changes are not visible during panel preparation.
#' One retained-connection preparation slot is shared by concurrent calls;
#' nested or concurrent preparation errors and callers may retry.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param source_path One indexed BAM/CRAM path or URI.
#' @param sample_id Nonempty sample identity assigned to the extracted rows.
#' @param reference_path Reference FASTA matching the panel and alignments.
#' @param panel_table Name of a committed ordered panel table or view.
#' @param panel_parquet Ordinary panel Parquet path, instead of `panel_table`.
#' @param index_path Optional explicit BAM/CRAM index path.
#' @param reference_index_path Optional explicit FASTA index path.
#' @param min_mapq Minimum alignment mapping quality.
#' @param min_baseq Minimum observed-base quality. Missing qualities pass only
#'   when this is zero.
#' @param require_flags SAM flag bits that every retained alignment must have.
#' @param exclude_flags SAM flag bits that exclude an alignment.
#' @param overlap_policy Either `"hileup_v0.1.0"` encounter-order mate
#'   suppression or `"none"`.
#' @param decompression_threads Number of htslib decompression worker threads.
#' @param worker_count Number of independently schedulable panel shards, from
#'   1 through 64. Each nonempty shard opens its own reader/reference state.
#' @param max_depth Maximum admitted pileup depth before an explicit error.
#' @param max_overlap_qnames Per-site, per-job capacity for overlap-suppression
#'   names.
#' @param max_sites Maximum panel cardinality.
#' @param max_region_bytes Per-job capacity for the indexed multi-region request.
#' @param remote_block_bytes Remote alignment block size per worker handle.
#' @param remote_cache_bytes Remote alignment cache size per worker handle.
#' @param reference_cache_bytes Remote reference cache size per worker handle.
#' @param table_name Optional output table. `NULL` returns a data frame.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame ordered by `site_index` if `table_name` is `NULL`;
#'   otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_bam_counts <- function(
  con, source_path, sample_id, reference_path,
  panel_table = NULL, panel_parquet = NULL,
  index_path = NULL, reference_index_path = NULL,
  min_mapq = 1, min_baseq = 0, require_flags = 0, exclude_flags = 1796,
  overlap_policy = c("hileup_v0.1.0", "none"), decompression_threads = 0,
  worker_count = 1, max_depth = 100000, max_overlap_qnames = 100000,
  max_sites = 1000000,
  max_region_bytes = 67108864, remote_block_bytes = 1048576,
  remote_cache_bytes = 67108864, reference_cache_bytes = 67108864,
  table_name = NULL, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  .somalier_scalar_text(source_path, "source_path")
  .somalier_scalar_text(sample_id, "sample_id")
  .somalier_scalar_text(reference_path, "reference_path")
  if (nchar(sample_id, type = "bytes") > 1024L) {
    stop("sample_id must be at most 1,024 bytes", call. = FALSE)
  }
  if (!is.null(index_path)) .somalier_scalar_text(index_path, "index_path")
  if (!is.null(reference_index_path)) {
    .somalier_scalar_text(reference_index_path, "reference_index_path")
  }
  overlap_policy <- match.arg(overlap_policy)

  min_mapq <- .somalier_bounded_whole_number(min_mapq, "min_mapq", 0, 255)
  min_baseq <- .somalier_bounded_whole_number(min_baseq, "min_baseq", 0, 255)
  require_flags <- .somalier_bounded_whole_number(
    require_flags, "require_flags", 0, 65535
  )
  exclude_flags <- .somalier_bounded_whole_number(
    exclude_flags, "exclude_flags", 0, 65535
  )
  decompression_threads <- .somalier_bounded_whole_number(
    decompression_threads, "decompression_threads", 0, .Machine$integer.max
  )
  worker_count <- .somalier_bounded_whole_number(
    worker_count, "worker_count", 1, 64
  )
  max_depth <- .somalier_bounded_whole_number(
    max_depth, "max_depth", 1, .Machine$integer.max - 1
  )
  max_overlap_qnames <- .somalier_bounded_whole_number(
    max_overlap_qnames, "max_overlap_qnames", 0, .Machine$integer.max
  )
  max_sites <- .somalier_bounded_whole_number(
    max_sites, "max_sites", 1, 100000000
  )
  max_region_bytes <- .somalier_bounded_whole_number(
    max_region_bytes, "max_region_bytes", 1, .Machine$integer.max
  )
  remote_block_bytes <- .somalier_bounded_whole_number(
    remote_block_bytes, "remote_block_bytes", 0, .Machine$integer.max
  )
  remote_cache_bytes <- .somalier_bounded_whole_number(
    remote_cache_bytes, "remote_cache_bytes", 0, .Machine$integer.max
  )
  reference_cache_bytes <- .somalier_bounded_whole_number(
    reference_cache_bytes, "reference_cache_bytes", 0, .Machine$integer.max
  )

  if (is.null(panel_table) == is.null(panel_parquet)) {
    stop("panel: supply exactly one table name or Parquet path", call. = FALSE)
  }
  if (!is.null(panel_table)) {
    .somalier_validate_name(panel_table, "panel_table")
  } else {
    .somalier_validate_name(panel_parquet, "panel_parquet")
  }
  arguments <- c(
    sql_quote_string(con, source_path),
    if (is.null(panel_table)) "NULL" else sql_quote_string(con, panel_table),
    sql_quote_string(con, sample_id), sql_quote_string(con, reference_path),
    if (!is.null(panel_parquet)) {
      paste0("panel_parquet := ", sql_quote_string(con, panel_parquet))
    },
    if (!is.null(index_path)) {
      paste0("index_path := ", sql_quote_string(con, index_path))
    },
    if (!is.null(reference_index_path)) {
      paste0(
        "reference_index_path := ",
        sql_quote_string(con, reference_index_path)
      )
    },
    paste0("min_mapq := ", min_mapq),
    paste0("min_baseq := ", min_baseq),
    paste0("require_flags := ", require_flags),
    paste0("exclude_flags := ", exclude_flags),
    paste0("overlap_policy := ", sql_quote_string(con, overlap_policy)),
    paste0("decompression_threads := ", decompression_threads),
    paste0("worker_count := ", worker_count),
    paste0("max_depth := ", max_depth),
    paste0("max_overlap_qnames := ", max_overlap_qnames),
    paste0("max_sites := ", max_sites),
    paste0("max_region_bytes := ", max_region_bytes),
    paste0("remote_block_bytes := ", remote_block_bytes),
    paste0("remote_cache_bytes := ", remote_cache_bytes),
    paste0("reference_cache_bytes := ", reference_cache_bytes)
  )
  query <- paste0(
    "SELECT * FROM duckhts_somalier_bam_counts(",
    paste(arguments, collapse = ", "), ") ORDER BY site_index"
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Prepare Somalier-Derived Sample Sketches
#'
#' Build packed, panel-verified relatedness sketches from measured A/B/other
#' count evidence. The panel must contain `assembly`, zero-based `site_index`,
#' `region`, one-based `position`, and uppercase single-base `allele_a` and
#' `allele_b`. Evidence must contain the same site identity columns plus
#' `sample_id` and nullable count columns `a`, `b`, and `other`. All three counts
#' are NULL for unavailable evidence; three measured zeros are not unavailable.
#' The native SQL preparation checks every evidence site's geometry and A/B
#' orientation against the ordered panel before computing its digest. Panel
#' alleles must be distinct uppercase single-base A/C/G/T with lexical A < B;
#' the exact X/Y aliases excluded by Somalier v0.3.4 are rejected. Other
#' contig aliases cannot be classified biologically from the region string.
#' The three-state calculation assumes diploid sites; count evidence alone
#' does not prove sample ploidy.
#'
#' Supply each source as either a table/view name or an ordinary Parquet path.
#' Parquet inputs are exposed through query-scoped temporary views; no private
#' sketch format or user-supplied panel digest is involved. The result has one
#' `sketch` struct per selected sample. Its packed words can be persisted with
#' DuckDB's usual Parquet `COPY` statement.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param evidence_table Name of an ordinary evidence table or view.
#' @param evidence_parquet Path to an evidence Parquet file, instead of
#'   `evidence_table`.
#' @param panel_table Name of the required ordered panel table or view. Panel
#'   assembly and region values are each limited to 1,024 bytes.
#' @param panel_parquet Path to the required ordered panel Parquet file, instead
#'   of `panel_table`.
#' @param table_name Optional output table. `NULL` returns a data frame.
#' @param sample_ids Optional nonempty vector of distinct sample IDs to retain.
#' @param min_depth Minimum A+B count depth for a relatedness genotype call.
#' @param min_het_balance Lower B/(A+B) balance accepted as heterozygous.
#' @param hom_balance_cutoff B/(A+B) balance below which a site is homozygous A;
#'   its upper symmetric limit determines homozygous B.
#' @param max_sites Positive per-sample panel capacity, at most 100,000,000.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_sketches <- function(
  con, evidence_table = NULL, evidence_parquet = NULL,
  panel_table = NULL, panel_parquet = NULL, table_name = NULL,
  sample_ids = NULL, min_depth = 7, min_het_balance = 0.3,
  hom_balance_cutoff = 0.01, max_sites = 1000000, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  min_depth <- .somalier_whole_number(min_depth, "min_depth", 2^53)
  max_sites <- .somalier_whole_number(max_sites, "max_sites", 100000000)
  min_het_balance <- .somalier_fraction(min_het_balance, "min_het_balance")
  hom_balance_cutoff <- .somalier_fraction(hom_balance_cutoff, "hom_balance_cutoff")
  if (hom_balance_cutoff > min_het_balance || min_het_balance > 0.5) {
    stop("require 0 <= hom_balance_cutoff <= min_het_balance <= 0.5", call. = FALSE)
  }
  .somalier_validate_sample_ids(sample_ids)

  evidence <- .somalier_source(con, evidence_table, evidence_parquet, "evidence")
  if (evidence$temporary) {
    on.exit(.somalier_drop_view(con, evidence$name), add = TRUE, after = FALSE)
  }
  panel <- .somalier_source(con, panel_table, panel_parquet, "panel")
  if (panel$temporary) {
    on.exit(.somalier_drop_view(con, panel$name), add = TRUE, after = FALSE)
  }
  evidence$name <- .somalier_select_samples(con, evidence$name, sample_ids)
  if (!is.null(sample_ids)) {
    on.exit(.somalier_drop_view(con, evidence$name), add = TRUE, after = FALSE)
  }

  query <- sprintf(
    paste0("SELECT * FROM duckhts_somalier_prepare_sketches(",
           "%s, %s, %s, %s, %s, max_sites := %s)"),
    sql_quote_string(con, evidence$name), sql_quote_string(con, panel$name),
    .somalier_quote_number(con, min_depth),
    .somalier_quote_number(con, min_het_balance),
    .somalier_quote_number(con, hom_balance_cutoff),
    .somalier_quote_number(con, max_sites)
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Compare Somalier-Derived Sample Sketches
#'
#' Compute the named relatedness and concordance statistics for either every
#' distinct pair in one sketch relation or the ordered pairs in `pairs_table`.
#' The sketch relation must contain one non-NULL `sketch` struct per sample.
#' Selected pairs are an ordinary relation with `sample_a` and `sample_b`
#' columns. Missing or duplicate sample and pair identities error instead of
#' silently dropping or multiplying requested comparisons. The native kernel
#' checks assembly, ordered-panel digest, classification settings, mask shape,
#' and mask contents for each comparison. No SQL row-order guarantee is implied.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param sketches_table Name of a prepared sketch table or view.
#' @param sketches_parquet Path to ordinary Parquet-persisted sketches instead
#'   of `sketches_table`.
#' @param pairs_table Optional name of an ordered-pair table or view; `NULL`
#'   requests all distinct unordered sample pairs.
#' @param table_name Optional output table. `NULL` returns a data frame, which
#'   should be used only for a result small enough to fit in R memory.
#' @param max_sites Positive per-pair panel capacity, at most 100,000,000.
#' @param overwrite Whether an existing output table may be replaced.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_relatedness <- function(
  con, sketches_table = NULL, sketches_parquet = NULL,
  pairs_table = NULL, table_name = NULL, max_sites = 1000000,
  overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  max_sites <- .somalier_whole_number(max_sites, "max_sites", 100000000)
  if (!is.null(pairs_table)) {
    .somalier_validate_name(pairs_table, "pairs_table")
  }
  sketches <- .somalier_source(con, sketches_table, sketches_parquet, "sketches")
  if (sketches$temporary) {
    on.exit(.somalier_drop_view(con, sketches$name), add = TRUE, after = FALSE)
  }
  sketch_relation <- sql_quote_identifier(con, sketches$name)
  sketch_check <- DBI::dbGetQuery(con, sprintf(
    paste0("SELECT count(*) AS rows, count(DISTINCT sketch.sample_id) AS samples, ",
           "count(*) FILTER (WHERE sketch IS NULL OR sketch.sample_id IS NULL ",
           "OR sketch.sample_id = '') AS invalid FROM %s"),
    sketch_relation
  ))
  if (sketch_check$invalid[[1]] != 0 ||
      sketch_check$rows[[1]] != sketch_check$samples[[1]]) {
    stop("sketches must have one non-NULL sketch per distinct sample", call. = FALSE)
  }

  if (is.null(pairs_table)) {
    from_sql <- sprintf(
      "FROM %s a JOIN %s b ON a.sketch.sample_id < b.sketch.sample_id",
      sketch_relation, sketch_relation
    )
  } else {
    pair_relation <- sql_quote_identifier(con, pairs_table)
    pair_check <- DBI::dbGetQuery(con, sprintf(
      paste0("SELECT count(*) AS rows, ",
             "count(DISTINCT struct_pack(sample_a := p.sample_a, ",
             "sample_b := p.sample_b)) AS distinct_pairs, ",
             "count(*) FILTER (WHERE p.sample_a IS NULL OR p.sample_b IS NULL ",
             "OR p.sample_a = p.sample_b) AS invalid, ",
             "count(*) FILTER (WHERE a.sketch IS NULL OR b.sketch IS NULL) AS missing ",
             "FROM %s p LEFT JOIN %s a ON a.sketch.sample_id = p.sample_a ",
             "LEFT JOIN %s b ON b.sketch.sample_id = p.sample_b"),
      pair_relation, sketch_relation, sketch_relation
    ))
    if (pair_check$invalid[[1]] != 0 || pair_check$missing[[1]] != 0 ||
        pair_check$rows[[1]] != pair_check$distinct_pairs[[1]]) {
      stop("pairs_table requires distinct, non-NULL, different existing sample IDs",
           call. = FALSE)
    }
    from_sql <- sprintf(
      paste0("FROM %s p JOIN %s a ON a.sketch.sample_id = p.sample_a ",
             "JOIN %s b ON b.sketch.sample_id = p.sample_b"),
      pair_relation, sketch_relation, sketch_relation
    )
  }
  query <- sprintf(
    "SELECT unnest(duckhts_somalier_relatedness(a.sketch, b.sketch, %s)) %s",
    .somalier_quote_number(con, max_sites), from_sql
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Estimate Per-Sample Contamination with CHARR
#'
#' Apply the Somalier-derived CHARR estimator to measured A/B/other count
#' evidence aligned to an ordered panel and population-B allele frequencies.
#' `frequency_table` contains the panel's six identity columns plus
#' `population_b_af`; it must cover every panel site exactly once. Evidence and
#' frequency identities, coordinates, and A/B orientation are checked against
#' the panel before their digests are derived. A result with no usable
#' homozygous-like evidence has status `no_evidence` and a NULL estimate.
#'
#' Each of the evidence, panel, and frequency inputs is supplied as exactly one
#' named table/view or ordinary Parquet path. Optional sample selection is exact:
#' every requested ID must occur. Results contain the panel and frequency
#' digests, usable-site denominators, numerical status, and all filter settings.
#'
#' @inheritParams rduckhts_somalier_sketches
#' @param min_depth Minimum measured A+B depth for CHARR homozygous-like eligibility.
#' @param frequency_table Name of the required population-frequency table or view.
#' @param frequency_parquet Path to the required population-frequency Parquet
#'   file, instead of `frequency_table`.
#' @param max_depth Maximum supported A+B depth for the binomial eligibility test,
#'   at most 1,000,000.
#' @param hom_minor_rate Expected minor-read rate used to recognize
#'   homozygous-like anchors.
#' @param hom_tail_alpha Binomial upper-tail threshold for homozygous-like
#'   eligibility.
#' @param max_threshold_work Positive cumulative limit on exact binomial
#'   certification steps, at most 100,000,000. Distinct observed depths are
#'   certified once per call and shared across samples.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_charr <- function(
  con, evidence_table = NULL, evidence_parquet = NULL,
  panel_table = NULL, panel_parquet = NULL,
  frequency_table = NULL, frequency_parquet = NULL,
  table_name = NULL, sample_ids = NULL, min_depth = 15,
  max_depth = 1000000, hom_minor_rate = 0.12, hom_tail_alpha = 0.002,
  max_threshold_work = 16000000, max_sites = 1000000, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  min_depth <- .somalier_whole_number(min_depth, "min_depth", 1000000)
  max_depth <- .somalier_whole_number(max_depth, "max_depth", 1000000)
  max_threshold_work <- .somalier_whole_number(
    max_threshold_work, "max_threshold_work", 100000000
  )
  max_sites <- .somalier_whole_number(max_sites, "max_sites", 100000000)
  if (min_depth > max_depth) {
    stop("min_depth must be no larger than max_depth", call. = FALSE)
  }
  hom_minor_rate <- .somalier_open_fraction(
    hom_minor_rate, "hom_minor_rate", upper = 0.5
  )
  hom_tail_alpha <- .somalier_open_fraction(hom_tail_alpha, "hom_tail_alpha")
  .somalier_validate_sample_ids(sample_ids)

  evidence <- .somalier_source(con, evidence_table, evidence_parquet, "evidence")
  if (evidence$temporary) {
    on.exit(.somalier_drop_view(con, evidence$name), add = TRUE, after = FALSE)
  }
  panel <- .somalier_source(con, panel_table, panel_parquet, "panel")
  if (panel$temporary) {
    on.exit(.somalier_drop_view(con, panel$name), add = TRUE, after = FALSE)
  }
  frequency <- .somalier_source(
    con, frequency_table, frequency_parquet, "frequency"
  )
  if (frequency$temporary) {
    on.exit(.somalier_drop_view(con, frequency$name), add = TRUE, after = FALSE)
  }
  evidence$name <- .somalier_select_samples(con, evidence$name, sample_ids)
  if (!is.null(sample_ids)) {
    on.exit(.somalier_drop_view(con, evidence$name), add = TRUE, after = FALSE)
  }

  query <- sprintf(
    paste0("SELECT unnest(contamination) FROM duckhts_somalier_charr(",
           "%s, %s, %s, min_depth := %s, max_depth := %s, ",
           "hom_minor_rate := %s, hom_tail_alpha := %s, ",
           "max_threshold_work := %s, max_sites := %s)"),
    sql_quote_string(con, evidence$name), sql_quote_string(con, panel$name),
    sql_quote_string(con, frequency$name),
    .somalier_quote_number(con, min_depth),
    .somalier_quote_number(con, max_depth),
    .somalier_quote_number(con, hom_minor_rate),
    .somalier_quote_number(con, hom_tail_alpha),
    .somalier_quote_number(con, max_threshold_work),
    .somalier_quote_number(con, max_sites)
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

#' Estimate Directional Contamination Against Matched Anchors
#'
#' Estimate one contamination fraction for each explicitly ordered
#' receiver/anchor pair. The pair relation must contain distinct nonempty
#' `receiver_id` and `anchor_id` columns. Every listed sample must have one
#' count tuple for every ordered panel site. Both samples and population-B
#' frequencies use the panel's A/B orientation; the anchor supplies the expected
#' uncontaminated receiver genotype and is not interpreted as the contaminating
#' donor. Reversing a pair is therefore a different analysis.
#'
#' Sample identifiers and assembly strings are limited to 1,024 bytes,
#' including when prepared profiles are persisted and supplied directly.
#'
#' The returned score omits alpha-independent binomial coefficients. It compares
#' alpha candidates for the same observed receiver/anchor pair and is not an
#' absolute likelihood comparable between pairs. No usable evidence yields
#' status `no_evidence` with NULL alpha and relative log-likelihood.
#'
#' @inheritParams rduckhts_somalier_charr
#' @param min_depth Minimum measured receiver A+B depth for matched-contamination eligibility.
#' @param pairs_table Name of an ordinary ordered-pair table or view.
#' @param pairs_parquet Path to an ordered-pair Parquet file, instead of
#'   `pairs_table`.
#' @param error_rate Count error probability used by the fitted likelihood.
#' @param min_probability Positive probability floor used in logarithms.
#' @param min_prior_frequency Lower population-frequency clamp.
#' @param alpha_min,alpha_max Closed contamination search range, satisfying
#'   `0 <= alpha_min < alpha_max <= 1`.
#' @param grid_step Positive initial grid spacing, at most one.
#' @param refine_tolerance Positive local-refinement tolerance.
#' @param max_evaluations Positive bound on likelihood evaluations.
#' @return A data frame if `table_name` is `NULL`; otherwise invisible `TRUE`.
#' @export
rduckhts_somalier_matched_contamination <- function(
  con, evidence_table = NULL, evidence_parquet = NULL,
  panel_table = NULL, panel_parquet = NULL,
  frequency_table = NULL, frequency_parquet = NULL,
  pairs_table = NULL, pairs_parquet = NULL, table_name = NULL,
  min_depth = 15, max_depth = 1000000, hom_minor_rate = 0.05,
  hom_tail_alpha = 0.001, error_rate = 0.002,
  min_probability = 1e-10, min_prior_frequency = 1e-6,
  alpha_min = 0, alpha_max = 1, grid_step = 0.01,
  refine_tolerance = 1e-10, max_evaluations = 4096,
  max_threshold_work = 16000000, max_sites = 1000000, overwrite = FALSE
) {
  .somalier_validate_output(con, table_name, overwrite)
  min_depth <- .somalier_whole_number(min_depth, "min_depth", 1000000)
  max_depth <- .somalier_whole_number(max_depth, "max_depth", 1000000)
  max_evaluations <- .somalier_whole_number(
    max_evaluations, "max_evaluations", .Machine$integer.max
  )
  max_threshold_work <- .somalier_whole_number(
    max_threshold_work, "max_threshold_work", 100000000
  )
  max_sites <- .somalier_whole_number(max_sites, "max_sites", 100000000)
  if (min_depth > max_depth) {
    stop("min_depth must be no larger than max_depth", call. = FALSE)
  }
  hom_minor_rate <- .somalier_open_fraction(
    hom_minor_rate, "hom_minor_rate", upper = 0.5
  )
  hom_tail_alpha <- .somalier_open_fraction(hom_tail_alpha, "hom_tail_alpha")
  error_rate <- .somalier_half_open_fraction(error_rate, "error_rate", 0.5)
  min_probability <- .somalier_open_fraction(
    min_probability, "min_probability", upper = 0.5
  )
  min_prior_frequency <- .somalier_open_fraction(
    min_prior_frequency, "min_prior_frequency", upper = 0.5
  )
  alpha_min <- .somalier_fraction(alpha_min, "alpha_min")
  alpha_max <- .somalier_fraction(alpha_max, "alpha_max")
  grid_step <- .somalier_open_fraction(grid_step, "grid_step", upper = 1,
                                       include_upper = TRUE)
  refine_tolerance <- .somalier_positive_number(
    refine_tolerance, "refine_tolerance"
  )
  if (alpha_min >= alpha_max) {
    stop("alpha_min must be smaller than alpha_max", call. = FALSE)
  }

  evidence <- .somalier_source(con, evidence_table, evidence_parquet, "evidence")
  if (evidence$temporary) {
    on.exit(.somalier_drop_view(con, evidence$name), add = TRUE, after = FALSE)
  }
  panel <- .somalier_source(con, panel_table, panel_parquet, "panel")
  if (panel$temporary) {
    on.exit(.somalier_drop_view(con, panel$name), add = TRUE, after = FALSE)
  }
  frequency <- .somalier_source(
    con, frequency_table, frequency_parquet, "frequency"
  )
  if (frequency$temporary) {
    on.exit(.somalier_drop_view(con, frequency$name), add = TRUE, after = FALSE)
  }
  pairs <- .somalier_source(con, pairs_table, pairs_parquet, "pairs")
  if (pairs$temporary) {
    on.exit(.somalier_drop_view(con, pairs$name), add = TRUE, after = FALSE)
  }

  query <- sprintf(
    paste0("SELECT unnest(contamination) FROM ",
           "duckhts_somalier_matched_contamination(%s, %s, %s, %s, ",
           "min_depth := %s, max_depth := %s, hom_minor_rate := %s, ",
           "hom_tail_alpha := %s, error_rate := %s, min_probability := %s, ",
           "min_prior_frequency := %s, alpha_min := %s, alpha_max := %s, ",
           "grid_step := %s, refine_tolerance := %s, max_evaluations := %s, ",
           "max_threshold_work := %s, max_sites := %s)"),
    sql_quote_string(con, evidence$name), sql_quote_string(con, panel$name),
    sql_quote_string(con, frequency$name), sql_quote_string(con, pairs$name),
    .somalier_quote_number(con, min_depth),
    .somalier_quote_number(con, max_depth),
    .somalier_quote_number(con, hom_minor_rate),
    .somalier_quote_number(con, hom_tail_alpha),
    .somalier_quote_number(con, error_rate),
    .somalier_quote_number(con, min_probability),
    .somalier_quote_number(con, min_prior_frequency),
    .somalier_quote_number(con, alpha_min),
    .somalier_quote_number(con, alpha_max),
    .somalier_quote_number(con, grid_step),
    .somalier_quote_number(con, refine_tolerance),
    .somalier_quote_number(con, max_evaluations),
    .somalier_quote_number(con, max_threshold_work),
    .somalier_quote_number(con, max_sites)
  )
  .somalier_publish_query(con, query, table_name, overwrite)
}

.somalier_validate_name <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(name, " must be one nonempty table name or Parquet path", call. = FALSE)
  }
  invisible(value)
}

.somalier_scalar_text <- function(value, name, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!allow_empty && !nzchar(value))) {
    qualifier <- if (allow_empty) "one non-missing string" else "one nonempty string"
    stop(name, " must be ", qualifier, call. = FALSE)
  }
  invisible(value)
}

.somalier_validate_output <- function(con, table_name, overwrite) {
  if (!is.null(table_name)) .somalier_validate_name(table_name, "table_name")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  .duckhts_check_table_target(con, table_name, overwrite)
  invisible(TRUE)
}

.somalier_whole_number <- function(value, name, maximum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value > maximum || value != floor(value)) {
    stop(name, " must be one positive, exactly representable whole number <= ",
         format(maximum, scientific = FALSE), call. = FALSE)
  }
  value
}

.somalier_bounded_whole_number <- function(value, name, minimum, maximum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < minimum || value > maximum || value != floor(value)) {
    stop(
      name, " must be one exactly representable whole number in [",
      format(minimum, scientific = FALSE), ", ",
      format(maximum, scientific = FALSE), "]",
      call. = FALSE
    )
  }
  value
}

.somalier_fraction <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value > 1) {
    stop(name, " must be one finite fraction in [0, 1]", call. = FALSE)
  }
  value
}

.somalier_open_fraction <- function(value, name, upper = 1,
                                    include_upper = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(name, " must be one finite number", call. = FALSE)
  }
  invalid_upper <- if (include_upper) value > upper else value >= upper
  if (value <= 0 || invalid_upper) {
    closing <- if (include_upper) "]" else ")"
    stop(name, " must be one finite number in (0, ", upper, closing,
         call. = FALSE)
  }
  value
}

.somalier_half_open_fraction <- function(value, name, upper) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value >= upper) {
    stop(name, " must be one finite number in [0, ", upper, ")", call. = FALSE)
  }
  value
}

.somalier_positive_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(name, " must be one positive finite number", call. = FALSE)
  }
  value
}

.somalier_validate_sample_ids <- function(sample_ids) {
  if (!is.null(sample_ids) &&
      (!is.character(sample_ids) || !length(sample_ids) || anyNA(sample_ids) ||
       any(!nzchar(sample_ids)) || anyDuplicated(sample_ids))) {
    stop("sample_ids must be a nonempty vector of distinct, non-missing IDs", call. = FALSE)
  }
  invisible(TRUE)
}

.somalier_select_samples <- function(con, evidence_name, sample_ids) {
  if (is.null(sample_ids)) return(evidence_name)
  sample_literal <- sql_varchar_list_literal(con, sample_ids, "sample_ids")
  matched <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(DISTINCT sample_id) AS n FROM %s WHERE sample_id IN (SELECT unnest(%s))",
    sql_quote_identifier(con, evidence_name), sample_literal
  ))$n[[1]]
  if (matched != length(sample_ids)) {
    stop("every sample_ids value must occur in the evidence relation", call. = FALSE)
  }
  selected_name <- .somalier_temp_view_name()
  DBI::dbExecute(con, sprintf(
    "CREATE TEMP VIEW %s AS SELECT * FROM %s WHERE sample_id IN (SELECT unnest(%s))",
    sql_quote_identifier(con, selected_name), sql_quote_identifier(con, evidence_name),
    sample_literal
  ))
  selected_name
}

.somalier_quote_number <- function(con, value) {
  as.character(DBI::dbQuoteLiteral(con, value))
}

.somalier_temp_view_name <- function() {
  basename(tempfile(pattern = "rduckhts_somalier_"))
}

.somalier_drop_view <- function(con, name) {
  try(
    DBI::dbExecute(con, paste("DROP VIEW IF EXISTS", sql_quote_identifier(con, name))),
    silent = TRUE
  )
  invisible(TRUE)
}

.somalier_source <- function(con, table, parquet, name) {
  if (is.null(table) == is.null(parquet)) {
    stop(name, ": supply exactly one table name or Parquet path", call. = FALSE)
  }
  if (!is.null(table)) {
    .somalier_validate_name(table, paste0(name, "_table"))
    view_name <- .somalier_temp_view_name()
    DBI::dbExecute(con, sprintf(
      "CREATE TEMP VIEW %s AS SELECT * FROM %s",
      sql_quote_identifier(con, view_name), sql_quote_identifier(con, table)
    ))
    return(list(name = view_name, temporary = TRUE))
  }
  .somalier_validate_name(parquet, paste0(name, "_parquet"))
  view_name <- .somalier_temp_view_name()
  DBI::dbExecute(con, sprintf(
    "CREATE TEMP VIEW %s AS SELECT * FROM read_parquet(%s)",
    sql_quote_identifier(con, view_name), sql_quote_string(con, parquet)
  ))
  list(name = view_name, temporary = TRUE)
}

.somalier_publish_query <- function(con, query, table_name, overwrite) {
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  .duckhts_create_table(con, table_name, query, overwrite)
  invisible(TRUE)
}
