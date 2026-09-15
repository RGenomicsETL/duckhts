somalier_fixture_counts <- function(sample_number) {
  site_index <- 0:16999
  genotype <- (site_index * 7 + sample_number * 13) %% 10
  unavailable <- (site_index + sample_number * 17) %% 97 == 0
  minor <- as.numeric((site_index + sample_number * 11) %% 23 == 0)
  a <- ifelse(genotype <= 3, 30 - minor, ifelse(genotype <= 6, 15, minor))
  b <- 30 - a
  list(a = a, b = b, genotype = genotype, unavailable = unavailable)
}

somalier_charr_fixture_expectation <- function(sample_number) {
  counts <- somalier_fixture_counts(sample_number)
  af <- (10 + ((0:16999) * 17) %% 39) / 100
  # At depth 30, every fixture minor count is 0 or 1. The declared binomial
  # tail admits one minor count; heterozygous and unavailable sites do not fit.
  stopifnot(stats::pbinom(0, 30, 0.12, lower.tail = FALSE) > 0.002)
  hom_a <- !counts$unavailable & counts$genotype <= 3
  hom_b <- !counts$unavailable & counts$genotype >= 7
  contribution <- c(counts$b[hom_a] / (af[hom_a] * 30),
                    counts$a[hom_b] / ((1 - af[hom_b]) * 30))
  list(unavailable_sites = sum(counts$unavailable), usable_sites = length(contribution),
       usable_hom_a = sum(hom_a), usable_hom_b = sum(hom_b),
       estimate = mean(contribution))
}

somalier_matched_fixture_expectation <- function(receiver_number = 2L,
                                                  anchor_number = 1L) {
  receiver <- somalier_fixture_counts(receiver_number)
  anchor <- somalier_fixture_counts(anchor_number)
  af <- (10 + ((0:16999) * 17) %% 39) / 100
  stopifnot(stats::pbinom(0, 30, 0.05, lower.tail = FALSE) > 0.001)
  hom_a <- !anchor$unavailable & anchor$genotype <= 3
  hom_b <- !anchor$unavailable & anchor$genotype >= 7
  usable <- !receiver$unavailable & (hom_a | hom_b)
  a <- receiver$a[usable]
  b <- receiver$b[usable]
  frequency <- af[usable]
  clean_b <- as.numeric(hom_b[usable])
  priors <- cbind((1 - frequency)^2, 2 * frequency * (1 - frequency),
                  frequency^2)
  score <- function(alpha) {
    terms <- matrix(0, nrow = length(a), ncol = 3L)
    for (genotype in 0:2) {
      latent_b <- (1 - alpha) * clean_b + alpha * genotype / 2
      observed_b <- latent_b * 0.998 + (1 - latent_b) * 0.002
      observed_b <- pmax(1e-10, pmin(1 - 1e-10, observed_b))
      terms[, genotype + 1L] <- log(priors[, genotype + 1L]) +
        b * log(observed_b) + a * log1p(-observed_b)
    }
    maximum <- pmax(terms[, 1L], terms[, 2L], terms[, 3L])
    sum(maximum + log(exp(terms[, 1L] - maximum) +
                      exp(terms[, 2L] - maximum) + exp(terms[, 3L] - maximum)))
  }
  grid <- seq(0, 1, by = 0.01)
  grid_scores <- vapply(grid, score, numeric(1L))
  best <- which.max(grid_scores)
  interval <- c(grid[[max(1L, best - 1L)]], grid[[min(length(grid), best + 1L)]])
  refined <- stats::optimize(score, interval = interval, maximum = TRUE, tol = 1e-9)
  candidates <- c(grid[[best]], refined$maximum, 0, 1)
  values <- vapply(candidates, score, numeric(1L))
  optimum <- candidates[[which.max(values)]]
  list(usable_sites = sum(usable), alpha = optimum,
       relative_log_likelihood = max(values))
}

somalier_peak_rss_kib <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  status <- readLines("/proc/self/status", warn = FALSE)
  peak <- status[startsWith(status, "VmHWM:")]
  if (length(peak) != 1L) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", peak))
}

somalier_benchmark_run <- function(case, extension, paths, output_dir, threads = 4L) {
  stopifnot(case %in% c("panel_hash", "sketches", "related_selected", "related_all",
                        "charr", "matched_anchor", "matched_memory", "end_to_end"),
            file.exists(extension), all(file.exists(paths)), threads >= 1L)
  suppressPackageStartupMessages(library(DBI))
  quote_string <- function(value) as.character(dbQuoteString(con, value))
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste("LOAD", dbQuoteString(con, extension)))
  dbExecute(con, sprintf("SET threads=%d", as.integer(threads)))
  for (name in names(paths)) {
    dbExecute(con, sprintf("CREATE TEMP VIEW %s AS SELECT * FROM read_parquet(%s)",
      dbQuoteIdentifier(con, name), quote_string(paths[[name]])))
  }
  run_matched <- function() dbExecute(con, "CREATE OR REPLACE TEMP TABLE matched_anchor AS
    SELECT * FROM duckhts_somalier_matched_contamination(
      'evidence', 'panel', 'frequency', 'pairs', max_sites := 17000)")
  memory_case <- identical(case, "matched_memory")
  baseline_peak_rss_kib <- NA_real_
  measured_peak_rss_kib <- NA_real_
  if (memory_case) {
    baseline_peak_rss_kib <- somalier_peak_rss_kib()
    started <- proc.time()[["elapsed"]]
    run_matched()
    elapsed <- proc.time()[["elapsed"]] - started
    measured_peak_rss_kib <- somalier_peak_rss_kib()
  }
  input <- dbGetQuery(con, "SELECT
    (SELECT count(*) FROM panel)::DOUBLE AS sites,
    (SELECT count(*) FROM frequency)::DOUBLE AS frequency_rows,
    (SELECT count(*) FROM evidence)::DOUBLE AS evidence_rows,
    (SELECT count(DISTINCT sample_id) FROM evidence)::DOUBLE AS samples,
    (SELECT count(*) FROM pairs)::DOUBLE AS selected_pairs,
    (SELECT count(DISTINCT sample_id) FROM (
      SELECT receiver_id AS sample_id FROM pairs UNION ALL
      SELECT anchor_id AS sample_id FROM pairs))::DOUBLE AS selected_samples,
    (SELECT count(DISTINCT receiver_id) FROM pairs)::DOUBLE AS receiver_samples,
    (SELECT count(DISTINCT anchor_id) FROM pairs)::DOUBLE AS anchor_samples,
    (SELECT count(DISTINCT struct_pack(receiver_id := receiver_id,
      anchor_id := anchor_id)) FROM pairs)::DOUBLE AS distinct_ordered_pairs,
    (SELECT count(*) FROM pairs WHERE receiver_id = 'sample002'
      AND anchor_id = 'sample001')::DOUBLE AS retained_pair,
    (SELECT count(*) FROM pairs WHERE receiver_id = anchor_id
      OR receiver_id IS NULL OR anchor_id IS NULL
      OR TRY_CAST(substr(receiver_id, 7) AS INTEGER) IS NULL
      OR TRY_CAST(substr(anchor_id, 7) AS INTEGER) IS NULL
      OR TRY_CAST(substr(receiver_id, 7) AS INTEGER) NOT BETWEEN 1 AND
         (SELECT count(DISTINCT sample_id) FROM evidence)
      OR TRY_CAST(substr(anchor_id, 7) AS INTEGER) NOT BETWEEN 1 AND
         (SELECT count(DISTINCT sample_id) FROM evidence))::DOUBLE AS invalid_pairs,
    (SELECT count(*) FROM evidence
      WHERE a IS NULL AND b IS NULL AND other IS NULL)::DOUBLE
      AS unavailable_rows,
    (SELECT count(*) FROM evidence
      WHERE a IS NOT NULL AND b IS NOT NULL AND other IS NOT NULL)::DOUBLE
      AS measured_rows,
    (SELECT count(*) FROM evidence WHERE a IS NOT NULL AND
      (a + b != 30 OR other != 0 OR
       NOT (a = 15 AND b = 15 OR least(a, b) <= 1)))::DOUBLE
       AS invalid_measured_rows,
    (SELECT count(DISTINCT struct_pack(region := region, pos1 := position,
      allele_a := allele_a, allele_b := allele_b)) FROM panel)::DOUBLE
      AS physical_sites,
    (SELECT count(DISTINCT struct_pack(sample_id := sample_id, site_index := site_index))
      FROM evidence)::DOUBLE AS sample_sites,
    (SELECT count(*) FROM panel WHERE assembly != 'GRCh38'
      OR site_index >= 17000
      OR region != 'chr' || ((site_index % 22) + 1)::VARCHAR
      OR position != 1000000 + (site_index // 22) * 101
      OR allele_a >= allele_b)::DOUBLE AS invalid_panel_rows,
    (SELECT count(*) FROM frequency
      WHERE population_b_af != (10 + (site_index * 17) % 39)::DOUBLE / 100)::DOUBLE
      AS invalid_frequency_rows")
  stopifnot(input$sites == 17000, input$frequency_rows == 17000,
            input$evidence_rows == input$sites * input$samples,
            input$physical_sites == input$sites, input$sample_sites == input$evidence_rows,
            input$unavailable_rows + input$measured_rows == input$evidence_rows,
            input$invalid_measured_rows == 0,
            input$invalid_panel_rows == 0, input$invalid_frequency_rows == 0,
            input$selected_samples >= 2, input$selected_samples <= input$samples,
            input$selected_pairs >= input$selected_samples - 1,
            input$selected_pairs <= input$selected_samples * (input$selected_samples - 1),
            input$distinct_ordered_pairs == input$selected_pairs,
            input$invalid_pairs == 0, input$retained_pair == 1)
  hashes <- function() dbGetQuery(con, "SELECT
    duckhts_somalier_panel_sha256('panel') AS panel_sha256,
    duckhts_somalier_frequency_sha256('frequency', 'panel') AS frequency_sha256")
  prepare_sketches <- function() dbExecute(con, "CREATE OR REPLACE TEMP TABLE sketches AS
    SELECT * FROM duckhts_somalier_prepare_sketches(
      'evidence', 'panel', 7, 0.3, 0.01, max_sites := 17000)")
  selected_relatedness <- function() dbExecute(con, "CREATE OR REPLACE TEMP TABLE related_selected AS
    SELECT unnest(duckhts_somalier_relatedness(a.sketch, b.sketch, 17000))
    FROM pairs p JOIN sketches a ON a.sketch.sample_id = p.sample_a
    JOIN sketches b ON b.sketch.sample_id = p.sample_b")
  all_relatedness <- function() dbExecute(con, "CREATE OR REPLACE TEMP TABLE related_all AS
    SELECT unnest(duckhts_somalier_relatedness(a.sketch, b.sketch, 17000))
    FROM sketches a JOIN sketches b ON a.sketch.sample_id < b.sketch.sample_id")
  run_charr <- function() dbExecute(con, "CREATE OR REPLACE TEMP TABLE charr AS
    SELECT * FROM duckhts_somalier_charr(
      'evidence', 'panel', 'frequency', max_sites := 17000)")
  if (!memory_case) {
    started <- proc.time()[["elapsed"]]
    if (case == "panel_hash") identity <- hashes()
    if (case == "sketches") prepare_sketches()
    if (case == "related_selected") {
      prepare_sketches()
      started <- proc.time()[["elapsed"]]
      selected_relatedness()
    }
    if (case == "related_all") {
      prepare_sketches()
      started <- proc.time()[["elapsed"]]
      all_relatedness()
    }
    if (case == "charr") run_charr()
    if (case == "matched_anchor") run_matched()
  }
  output_names <- c("sketches", "related_selected", "related_all", "charr", "matched_anchor")
  output_files <- character()
  output_rows_by_relation <- stats::setNames(rep(0, length(output_names)), output_names)
  output_bytes_by_relation <- output_rows_by_relation
  if (case == "end_to_end") {
    identity <- hashes()
    prepare_sketches()
    selected_relatedness()
    all_relatedness()
    run_charr()
    run_matched()
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    orders <- c("sketch.sample_id", "sample_a, sample_b", "sample_a, sample_b",
                "contamination.sample_id",
                "contamination.receiver_id, contamination.anchor_id")
    output_files <- file.path(output_dir, paste0(output_names, ".parquet"))
    for (i in seq_along(output_names)) {
      dbExecute(con, sprintf("COPY (SELECT * FROM %s ORDER BY %s) TO %s
        (FORMAT PARQUET, COMPRESSION ZSTD)", output_names[[i]], orders[[i]],
        quote_string(output_files[[i]])))
    }
  }
  if (!memory_case) elapsed <- proc.time()[["elapsed"]] - started
  if (case %in% c("panel_hash", "end_to_end")) {
    stopifnot(nchar(identity$panel_sha256) == 64L,
              nchar(identity$frequency_sha256) == 64L,
              grepl("^[0-9a-f]{64}$", identity$panel_sha256),
              grepl("^[0-9a-f]{64}$", identity$frequency_sha256),
              identity$panel_sha256 != identity$frequency_sha256)
  }
  sketch_result_rows <- 0
  if (case %in% c("sketches", "related_selected", "related_all", "end_to_end")) {
    sketch_check <- dbGetQuery(con, "SELECT count(*)::DOUBLE AS rows,
      count(DISTINCT sketch.sample_id)::DOUBLE AS samples,
      min(sketch.site_count)::DOUBLE AS min_sites, max(sketch.site_count)::DOUBLE AS max_sites,
      min(len(sketch.hom_a))::DOUBLE AS min_words, max(len(sketch.hom_b))::DOUBLE AS max_words,
      count(DISTINCT sketch.panel_sha256)::DOUBLE AS panel_digests FROM sketches")
    stopifnot(sketch_check$rows == input$samples, sketch_check$samples == input$samples,
              sketch_check$min_sites == 17000, sketch_check$max_sites == 17000,
              sketch_check$min_words == 266, sketch_check$max_words == 266,
              sketch_check$panel_digests == 1)
    sketch_result_rows <- sketch_check$rows[[1L]]
  }
  selected_result_rows <- 0
  related_selected_joint_sites <- 0
  if (case %in% c("related_selected", "end_to_end")) {
    selected_check <- dbGetQuery(con, "WITH raw AS (
      SELECT p.sample_a, p.sample_b, e.site_index,
        CASE WHEN e.a IS NULL THEN -1 WHEN e.b::DOUBLE / (e.a + e.b) < 0.01 THEN 0
          WHEN e.b::DOUBLE / (e.a + e.b) > 0.99 THEN 2
          WHEN e.b::DOUBLE / (e.a + e.b) BETWEEN 0.3 AND 0.7 THEN 1 ELSE -1 END AS ga,
        CASE WHEN f.a IS NULL THEN -1 WHEN f.b::DOUBLE / (f.a + f.b) < 0.01 THEN 0
          WHEN f.b::DOUBLE / (f.a + f.b) > 0.99 THEN 2
          WHEN f.b::DOUBLE / (f.a + f.b) BETWEEN 0.3 AND 0.7 THEN 1 ELSE -1 END AS gb
      FROM pairs p
      JOIN evidence e ON e.sample_id = p.sample_a
      JOIN evidence f ON f.sample_id = p.sample_b
        AND f.site_index = e.site_index),
    oracle AS (SELECT sample_a, sample_b, count(*) FILTER (WHERE ga >= 0 AND gb >= 0) AS called,
      count(*) FILTER (WHERE ga >= 0 AND gb >= 0 AND ga != 1 AND gb != 1 AND ga != gb) AS ibs0,
      count(*) FILTER (WHERE ga >= 0 AND gb >= 0 AND ga = gb) AS ibs2
      FROM raw GROUP BY sample_a, sample_b)
    SELECT count(*) AS rows,
      count(DISTINCT struct_pack(sample_a := r.sample_a,
        sample_b := r.sample_b)) AS distinct_pairs,
      sum(r.jointly_called)::DOUBLE AS joint_sites,
      count(*) FILTER (WHERE r.jointly_called = o.called AND
      r.ibs0 = o.ibs0 AND r.ibs2 = o.ibs2) AS exact FROM related_selected r
    JOIN oracle o USING (sample_a, sample_b)")
    stopifnot(selected_check$rows == input$selected_pairs,
              selected_check$distinct_pairs == input$selected_pairs,
              selected_check$exact == input$selected_pairs)
    selected_result_rows <- selected_check$rows[[1L]]
    related_selected_joint_sites <- selected_check$joint_sites[[1L]]
  }
  all_result_rows <- 0
  related_all_joint_sites <- 0
  if (case %in% c("related_all", "end_to_end")) {
    expected_pairs <- input$samples * (input$samples - 1) / 2
    all_check <- dbGetQuery(con, "SELECT count(*)::DOUBLE AS rows,
      count(DISTINCT struct_pack(sample_a := sample_a, sample_b := sample_b))::DOUBLE AS pairs,
      sum(jointly_called)::DOUBLE AS joint_sites,
      count(*) FILTER (WHERE status = 'ok' AND jointly_called > 0)::DOUBLE AS evaluated
      FROM related_all")
    stopifnot(all_check$rows == expected_pairs, all_check$pairs == expected_pairs,
              all_check$evaluated == expected_pairs)
    key_check <- dbGetQuery(con, "WITH expected AS (
      SELECT a.sketch.sample_id AS sample_a, b.sketch.sample_id AS sample_b
      FROM sketches a JOIN sketches b ON a.sketch.sample_id < b.sketch.sample_id)
      SELECT count(*) FILTER (WHERE r.sample_a IS NULL)::DOUBLE AS missing,
        count(*) FILTER (WHERE e.sample_a IS NULL)::DOUBLE AS extra
      FROM expected e FULL OUTER JOIN related_all r
        ON e.sample_a = r.sample_a AND e.sample_b = r.sample_b")
    stopifnot(key_check$missing == 0, key_check$extra == 0)
    all_result_rows <- all_check$rows[[1L]]
    related_all_joint_sites <- all_check$joint_sites[[1L]]
  } else expected_pairs <- 0
  charr_result_rows <- 0
  charr_usable_sites <- 0
  if (case %in% c("charr", "end_to_end")) {
    charr_check <- dbGetQuery(con, "SELECT count(*)::DOUBLE AS rows,
      sum(contamination.site_count)::DOUBLE AS site_sum,
      sum(contamination.observed_sites)::DOUBLE AS observed_sum,
      sum(contamination.unavailable_sites)::DOUBLE AS unavailable_sum,
      sum(contamination.usable_sites)::DOUBLE AS usable_sum,
      count(*) FILTER (WHERE contamination.status = 'ok' AND
        contamination.usable_sites > 0)::DOUBLE AS evaluated FROM charr")
    stopifnot(charr_check$rows == input$samples,
              charr_check$site_sum == input$evidence_rows,
              charr_check$observed_sum == input$evidence_rows,
              charr_check$unavailable_sum == input$unavailable_rows,
              charr_check$evaluated == input$samples)
    charr_integer_oracle <- dbGetQuery(con, "WITH raw AS (
      SELECT e.sample_id, count(*) FILTER (WHERE e.a IS NULL) AS unavailable,
        count(*) FILTER (WHERE e.a IS NOT NULL AND e.a > e.b AND e.b <= 1) AS hom_a,
        count(*) FILTER (WHERE e.a IS NOT NULL AND e.b > e.a AND e.a <= 1) AS hom_b,
        sum(CASE WHEN e.a IS NOT NULL AND e.a > e.b AND e.b <= 1 THEN
          e.b::DOUBLE / (f.population_b_af * 30)
          WHEN e.a IS NOT NULL AND e.b > e.a AND e.a <= 1 THEN
          e.a::DOUBLE / ((1 - f.population_b_af) * 30) ELSE NULL END) /
        count(*) FILTER (WHERE e.a IS NOT NULL AND
          (e.a > e.b AND e.b <= 1 OR e.b > e.a AND e.a <= 1)) AS estimate
      FROM evidence e JOIN frequency f USING (site_index) GROUP BY e.sample_id)
      SELECT count(*) AS rows, count(*) FILTER (WHERE
        c.contamination.unavailable_sites = r.unavailable AND
        c.contamination.usable_hom_a = r.hom_a AND
        c.contamination.usable_hom_b = r.hom_b AND
        c.contamination.usable_sites = r.hom_a + r.hom_b AND
        abs(c.contamination.estimate - r.estimate) < 1e-8) AS exact
      FROM charr c JOIN raw r ON c.contamination.sample_id = r.sample_id")
    stopifnot(charr_integer_oracle$rows == input$samples,
              charr_integer_oracle$exact == input$samples)
    charr_result_rows <- charr_check$rows[[1L]]
    charr_usable_sites <- charr_check$usable_sum[[1L]]
    retained_charr <- dbGetQuery(con, "SELECT contamination.unavailable_sites AS unavailable,
      contamination.usable_sites AS usable,
      contamination.usable_hom_a AS hom_a,
      contamination.usable_hom_b AS hom_b,
      contamination.estimate AS estimate FROM charr
      WHERE contamination.sample_id = 'sample001'")
    expected_charr <- somalier_charr_fixture_expectation(1L)
    stopifnot(nrow(retained_charr) == 1L,
              retained_charr$unavailable == expected_charr$unavailable_sites,
              retained_charr$usable == expected_charr$usable_sites,
              retained_charr$hom_a == expected_charr$usable_hom_a,
              retained_charr$hom_b == expected_charr$usable_hom_b,
              abs(retained_charr$estimate - expected_charr$estimate) < 1e-10)
  }
  matched_result_rows <- 0
  matched_usable_sites <- 0
  matched_evaluations <- 0
  matched_likelihood_contributions <- 0
  if (case %in% c("matched_anchor", "matched_memory", "end_to_end")) {
    matched_check <- dbGetQuery(con, "SELECT count(*)::DOUBLE AS rows,
      count(DISTINCT struct_pack(receiver_id := contamination.receiver_id,
        anchor_id := contamination.anchor_id))::DOUBLE AS pairs,
      sum(contamination.observed_sites)::DOUBLE AS observed_sum,
      sum(contamination.usable_sites)::DOUBLE AS usable_sum,
      sum(contamination.evaluations)::DOUBLE AS evaluation_sum,
      sum(contamination.usable_sites::DOUBLE * contamination.evaluations)::DOUBLE
        AS contribution_sum,
      count(*) FILTER (WHERE contamination.status = 'ok' AND
        contamination.usable_sites > 0 AND contamination.evaluations > 0 AND
        contamination.evaluations <= contamination.max_evaluations AND
        isfinite(contamination.alpha) AND
        isfinite(contamination.relative_log_likelihood))::DOUBLE AS evaluated
      FROM matched_anchor")
    stopifnot(matched_check$rows == input$selected_pairs,
              matched_check$pairs == input$selected_pairs,
              matched_check$observed_sum == input$selected_pairs * input$sites,
              matched_check$evaluated == input$selected_pairs)
    matched_result_rows <- matched_check$rows[[1L]]
    exact_pairs <- dbGetQuery(con, "SELECT count(*)::DOUBLE AS exact FROM matched_anchor m
      JOIN pairs p ON m.contamination.receiver_id = p.receiver_id
        AND m.contamination.anchor_id = p.anchor_id")$exact[[1L]]
    stopifnot(exact_pairs == input$selected_pairs)
    matched_integer_oracle <- dbGetQuery(con, "WITH raw AS (
      SELECT p.receiver_id, p.anchor_id,
        count(*) FILTER (WHERE e.a IS NULL) AS receiver_unavailable,
        count(*) FILTER (WHERE f.a IS NULL) AS anchor_unavailable,
        count(*) FILTER (WHERE e.a IS NOT NULL AND f.a IS NOT NULL AND
          (f.a > f.b AND f.b <= 1 OR f.b > f.a AND f.a <= 1)) AS usable
      FROM pairs p JOIN evidence e ON e.sample_id = p.receiver_id
      JOIN evidence f ON f.sample_id = p.anchor_id AND f.site_index = e.site_index
      GROUP BY p.receiver_id, p.anchor_id)
      SELECT count(*) AS rows, count(*) FILTER (WHERE
        m.contamination.receiver_unavailable_sites = r.receiver_unavailable AND
        m.contamination.anchor_unavailable_sites = r.anchor_unavailable AND
        m.contamination.usable_sites = r.usable) AS exact
      FROM matched_anchor m JOIN raw r ON
        m.contamination.receiver_id = r.receiver_id AND
        m.contamination.anchor_id = r.anchor_id")
    stopifnot(matched_integer_oracle$rows == input$selected_pairs,
              matched_integer_oracle$exact == input$selected_pairs)
    matched_usable_sites <- matched_check$usable_sum[[1L]]
    matched_evaluations <- matched_check$evaluation_sum[[1L]]
    matched_likelihood_contributions <- matched_check$contribution_sum[[1L]]
    for (receiver_number in 2:min(input$selected_samples, 3L)) {
      retained_matched <- dbGetQuery(con, sprintf("SELECT
        contamination.usable_sites AS usable, contamination.alpha AS alpha,
        contamination.relative_log_likelihood AS score FROM matched_anchor
        WHERE contamination.receiver_id = 'sample%03d'
          AND contamination.anchor_id = 'sample001'", receiver_number))
      expected_matched <- somalier_matched_fixture_expectation(receiver_number, 1L)
      stopifnot(nrow(retained_matched) == 1L,
                retained_matched$usable == expected_matched$usable_sites,
                abs(retained_matched$alpha - expected_matched$alpha) < 1e-4,
                abs(retained_matched$score - expected_matched$relative_log_likelihood) < 1e-3)
    }
  }
  if (length(output_files)) {
    expected_output_rows <- c(input$samples, input$selected_pairs, expected_pairs,
                              input$samples, input$selected_pairs)
    observed_output_rows <- vapply(output_files, function(path) dbGetQuery(con,
      sprintf("SELECT count(*)::DOUBLE AS n FROM read_parquet(%s)", quote_string(path)))$n[[1L]],
      numeric(1L))
    stopifnot(all(observed_output_rows == as.numeric(expected_output_rows)))
    output_rows_by_relation[] <- observed_output_rows
    output_bytes_by_relation[] <- file.info(output_files)$size
  }
  fingerprint <- switch(case,
    panel_hash = paste(identity$panel_sha256, identity$frequency_sha256, sep = ":"),
    sketches = dbGetQuery(con, "SELECT hex(bit_xor(hash(sketch))) AS x FROM sketches")$x[[1L]],
    related_selected = dbGetQuery(con, "SELECT hex(bit_xor(hash(sample_a, sample_b,
      jointly_called, ibs0, ibs2, relatedness))) AS x FROM related_selected")$x[[1L]],
    related_all = dbGetQuery(con, "SELECT hex(bit_xor(hash(sample_a, sample_b,
      jointly_called, ibs0, ibs2, relatedness))) AS x FROM related_all")$x[[1L]],
    charr = dbGetQuery(con, "SELECT hex(bit_xor(hash(contamination))) AS x FROM charr")$x[[1L]],
    matched_anchor = dbGetQuery(con, "SELECT hex(bit_xor(hash(contamination))) AS x FROM matched_anchor")$x[[1L]],
    matched_memory = dbGetQuery(con, "SELECT hex(bit_xor(hash(contamination))) AS x FROM matched_anchor")$x[[1L]],
    end_to_end = paste(vapply(output_files, digest::digest, character(1L),
                              file = TRUE, algo = "sha256"), collapse = ":"))
  peak_kib <- if (memory_case) measured_peak_rss_kib else somalier_peak_rss_kib()
  result_rows <- switch(case, panel_hash = nrow(identity), sketches = sketch_result_rows,
    related_selected = selected_result_rows, related_all = all_result_rows,
    charr = charr_result_rows, matched_anchor = matched_result_rows,
    matched_memory = matched_result_rows,
    end_to_end = sketch_result_rows + selected_result_rows + all_result_rows +
      charr_result_rows + matched_result_rows)
  data.frame(case = case, seconds = elapsed, peak_process_rss_kib = peak_kib,
    baseline_peak_rss_kib = baseline_peak_rss_kib,
    measured_peak_rss_kib = measured_peak_rss_kib,
    peak_increase_kib = if (memory_case && is.finite(peak_kib) &&
                            is.finite(baseline_peak_rss_kib))
      max(0, peak_kib - baseline_peak_rss_kib) else NA_real_,
    sites = input$sites, samples = input$samples, evidence_rows = input$evidence_rows,
    unavailable_evidence_rows = input$unavailable_rows,
    measured_evidence_rows = input$measured_rows,
    selected_pairs = input$selected_pairs, selected_samples = input$selected_samples,
    selected_receiver_samples = input$receiver_samples,
    selected_anchor_samples = input$anchor_samples,
    related_selected_pairs = if (case %in% c("related_selected", "end_to_end"))
      input$selected_pairs else 0,
    related_all_pairs = if (case %in% c("related_all", "end_to_end"))
      expected_pairs else 0,
    related_selected_joint_sites = related_selected_joint_sites,
    related_all_joint_sites = related_all_joint_sites,
    charr_samples = if (case %in% c("charr", "end_to_end")) input$samples else 0,
    charr_usable_sites = charr_usable_sites,
    matched_pairs = if (case %in% c("matched_anchor", "matched_memory", "end_to_end"))
      input$selected_pairs else 0,
    matched_profile_samples = if (case %in% c("matched_anchor", "matched_memory",
                                             "end_to_end")) input$selected_samples else 0,
    matched_profile_input_rows = if (case %in% c("matched_anchor", "matched_memory",
                                                 "end_to_end")) input$selected_samples * input$sites
      else 0,
    matched_profile_elements = if (case %in% c("matched_anchor", "matched_memory",
                                              "end_to_end")) input$selected_samples * input$sites * 4
      else 0,
    matched_profile_payload_bytes = if (case %in% c("matched_anchor", "matched_memory",
                                                   "end_to_end")) input$selected_samples * input$sites * 10
      else 0,
    frequency_profile_elements = if (case %in% c("matched_anchor", "matched_memory",
                                                "end_to_end")) input$sites else 0,
    frequency_profile_payload_bytes = if (case %in% c("matched_anchor", "matched_memory",
                                                     "end_to_end")) input$sites * 8 else 0,
    matched_usable_sites = matched_usable_sites,
    matched_evaluations = matched_evaluations,
    matched_panel_site_visits = matched_evaluations * input$sites +
      if (case %in% c("matched_anchor", "matched_memory", "end_to_end"))
        input$selected_pairs * input$sites else 0,
    matched_likelihood_contributions = matched_likelihood_contributions,
    result_rows = result_rows,
    persisted_output_rows = sum(output_rows_by_relation),
    persisted_output_bytes = sum(output_bytes_by_relation),
    sketches_output_rows = output_rows_by_relation[["sketches"]],
    sketches_output_bytes = output_bytes_by_relation[["sketches"]],
    selected_output_rows = output_rows_by_relation[["related_selected"]],
    selected_output_bytes = output_bytes_by_relation[["related_selected"]],
    all_pairs_output_rows = output_rows_by_relation[["related_all"]],
    all_pairs_output_bytes = output_bytes_by_relation[["related_all"]],
    charr_output_rows = output_rows_by_relation[["charr"]],
    charr_output_bytes = output_bytes_by_relation[["charr"]],
    matched_output_rows = output_rows_by_relation[["matched_anchor"]],
    matched_output_bytes = output_bytes_by_relation[["matched_anchor"]],
    fingerprint = fingerprint, stringsAsFactors = FALSE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 9L)
  paths <- stats::setNames(args[3:6], c("panel", "frequency", "evidence", "pairs"))
  result <- somalier_benchmark_run(args[[1L]], args[[2L]], paths, args[[7L]],
                                   as.integer(args[[8L]]))
  utils::write.table(result, args[[9L]], sep = "\t", row.names = FALSE, quote = FALSE)
}
