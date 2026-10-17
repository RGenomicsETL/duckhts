# The planted-truth evaluation: synthetic chr20 children of two unrelated
# founders of one population, with planted autozygous segments, planted
# genotype errors and simulated read counts, decoded by the genotype and
# read-count arms of roh_truth_declaration and scored against the plan.
#
# Run from the repository root after sourcing common.R.

roh_truth_rng <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

# The founders of each population and the draw of A, B and C per child.
roh_truth_children <- function(pedigree_path, populations, children) {
  pedigree <- utils::read.table(pedigree_path, header = TRUE, sep = "", stringsAsFactors = FALSE,
                                colClasses = "character")
  founders <- pedigree[pedigree$FatherID == "0" & pedigree$MotherID == "0", , drop = FALSE]
  rows <- list()
  for (population in populations) {
    pool <- sort(founders$SampleID[founders$Population == population])
    if (length(pool) < 3L * children) {
      stop("population ", population, " has ", length(pool), " founders, fewer than ",
           3L * children, call. = FALSE)
    }
    drawn <- pool[sample.int(length(pool), 3L * children)]
    rows[[population]] <- data.frame(
      child = sprintf("%s_%02d", population, seq_len(children)), population = population,
      founder_a = drawn[seq_len(children)], founder_b = drawn[children + seq_len(children)],
      founder_c = drawn[2L * children + seq_len(children)], stringsAsFactors = FALSE)
  }
  plan <- do.call(rbind, rows)
  rownames(plan) <- NULL
  plan$haplotype_a <- sample.int(2L, nrow(plan), replace = TRUE)
  plan$haplotype_b <- sample.int(2L, nrow(plan), replace = TRUE)
  list(children = plan, founders = founders[founders$Population %in% populations,
                                            c("SampleID", "Population")])
}

# The planted segments of one child, as roh_truth_declaration describes them.
roh_truth_segments <- function(positions, declaration) {
  first <- positions[[1L]]
  last <- positions[[length(positions)]]
  window <- floor((positions - 1) / declaration$window_bases)
  window_sites <- tabulate(window - min(window) + 1L, nbins = max(window) - min(window) + 1L)
  starts <- numeric()
  ends <- numeric()
  for (length_bases in declaration$planted_lengths) {
    placed <- FALSE
    for (draw in seq_len(declaration$max_draws)) {
      start <- first + floor(stats::runif(1L) * (last - length_bases - first + 2))
      end <- start + length_bases - 1
      windows <- seq(floor((start - 1) / declaration$window_bases),
                     floor((end - 1) / declaration$window_bases)) - min(window) + 1L
      sparse <- any(windows < 1L | windows > length(window_sites)) ||
        any(window_sites[windows] < declaration$min_window_sites)
      apart <- all(start > ends + declaration$min_gap_bases | end < starts - declaration$min_gap_bases)
      if (!sparse && apart) {
        placed <- TRUE
        break
      }
    }
    if (!placed) stop("could not place a planted segment", call. = FALSE)
    starts <- c(starts, start)
    ends <- c(ends, end)
  }
  data.frame(segment = seq_along(starts), start = starts, end = ends,
             length = declaration$planted_lengths)
}

# Scan the source once for every founder of the populations at the panel
# sites. Returns the sites with each population's clamped frequency and, for
# each drawn founder, the code 2 x first allele + second allele of its phased
# genotype. Every founder genotype at every site must be phased, biallelic
# and called.
roh_truth_founder_haplotypes <- function(con, source_path, sites_path, founders, drawn, populations,
                                         af_clamp) {
  q <- roh_truth_quote
  DBI::dbExecute(con, sprintf(paste(
    "CREATE OR REPLACE TEMP TABLE panel_sites AS SELECT",
    "(row_number() OVER (ORDER BY pos) - 1)::INTEGER AS site_index, chrom, pos, ref, alt",
    "FROM read_parquet(%s)"), q(sites_path)))
  shape <- DBI::dbGetQuery(con, "SELECT count(*) AS n, count(DISTINCT chrom) AS chroms,
                                        count(DISTINCT pos) AS positions FROM panel_sites")
  if (shape$chroms != 1 || shape$positions != shape$n) {
    stop("the panel must hold one chromosome with one site per position", call. = FALSE)
  }
  column <- function(sample) sprintf("\"FORMAT_GT_%s\"", sample)
  alt_count <- function(sample) {
    sprintf("(substr(%1$s, 1, 1) = '1')::INTEGER + (substr(%1$s, 3, 1) = '1')::INTEGER", column(sample))
  }
  frequencies <- vapply(populations, function(population) {
    members <- founders$SampleID[founders$Population == population]
    sprintf("greatest(%s, least(%s, (%s) / %d.0)) AS \"af_%s\"", format(af_clamp), format(1 - af_clamp),
            paste(vapply(members, alt_count, ""), collapse = " + "), 2L * length(members), population)
  }, "")
  codes <- vapply(drawn, function(sample) {
    sprintf("(2 * (substr(%1$s, 1, 1) = '1')::INTEGER + (substr(%1$s, 3, 1) = '1')::INTEGER)::TINYINT AS \"code_%2$s\"",
            column(sample), sample)
  }, "")
  malformed <- sprintf("NOT regexp_full_match(%s, '[01][|][01]')", vapply(founders$SampleID, column, ""))
  DBI::dbExecute(con, sprintf(paste(
    "CREATE OR REPLACE TEMP TABLE panel_founders AS SELECT s.site_index, s.chrom, s.pos, s.ref, s.alt,",
    "(%s)::INTEGER AS malformed, %s, %s",
    "FROM read_bcf(%s, samples := %s) AS r",
    "JOIN panel_sites AS s ON r.CHROM = s.chrom AND r.POS = s.pos AND r.REF = s.ref",
    "AND len(r.ALT) = 1 AND r.ALT[1] = s.alt"),
    paste(sprintf("(%s)::INTEGER", malformed), collapse = " + "),
    paste(frequencies, collapse = ", "), paste(codes, collapse = ", "),
    q(source_path), q(paste(founders$SampleID, collapse = ","))))
  check <- DBI::dbGetQuery(con, paste(
    "SELECT count(*) AS rows, count(DISTINCT site_index) AS sites,",
    "coalesce(sum(malformed), 0) AS malformed FROM panel_founders"))
  if (check$rows != shape$n || check$sites != shape$n || check$malformed != 0) {
    stop("the source does not give one phased, biallelic, called genotype per founder and panel site: ",
         check$rows, " rows, ", check$sites, " sites, ", check$malformed, " malformed", call. = FALSE)
  }
  DBI::dbGetQuery(con, "SELECT * EXCLUDE (malformed) FROM panel_founders ORDER BY site_index")
}

roh_truth_genotype_strings <- c("0/0", "0/1", "1/1")

# The background stretches of one child: maximal sets of consecutive panel
# sites without a heterozygous true genotype, from their first to their last
# site, of at least `min_bases`.
roh_truth_background <- function(positions, dosage, min_bases) {
  heterozygous <- which(dosage == 1L)
  first <- c(1L, heterozygous + 1L)
  last <- c(heterozygous - 1L, length(positions))
  keep <- first <= last
  first <- first[keep]
  last <- last[keep]
  span <- positions[last] - positions[first] + 1
  long <- span >= min_bases
  data.frame(start = positions[first[long]], end = positions[last[long]])
}

# Whether each position lies in one of the disjoint intervals.
roh_truth_inside <- function(positions, intervals) {
  if (!nrow(intervals)) return(rep(FALSE, length(positions)))
  intervals <- intervals[order(intervals$start), , drop = FALSE]
  index <- findInterval(positions, intervals$start)
  inside <- index > 0L
  inside[inside] <- positions[inside] <= intervals$end[index[inside]]
  inside
}

# Build the children, their truth and their evidence in `work`: a VCF of the
# called genotypes, a Parquet of read counts, a Parquet of the population
# frequencies, the planted segments and the true genotypes (for the check).
roh_truth_build <- function(con, work, declaration, populations, children_per_population) {
  pedigree_path <- roh_truth_check_file(declaration$pedigree_artifact)
  source_path <- roh_truth_check_file(declaration$source_artifact)
  sites_path <- roh_truth_check_file(declaration$sites_artifact)

  roh_truth_rng(declaration$seed)
  drawn <- roh_truth_children(pedigree_path, populations, children_per_population)
  children <- drawn$children
  founder_samples <- unique(c(children$founder_a, children$founder_b, children$founder_c))
  panel <- roh_truth_founder_haplotypes(con, source_path, sites_path, drawn$founders,
                                        founder_samples, populations, declaration$af_clamp)
  positions <- as.numeric(panel$pos)
  sites <- nrow(panel)
  code <- function(sample) as.integer(panel[[paste0("code_", sample)]])
  allele <- function(sample, haplotype) {
    value <- code(sample)
    if (haplotype == 1L) value %/% 2L else value %% 2L
  }

  segments <- list()
  called <- matrix(0L, nrow = sites, ncol = nrow(children))
  truth <- matrix(0L, nrow = sites, ncol = nrow(children))
  realized <- list()
  heterozygous <- list()
  background <- list()
  counts_table <- "simulated_counts"
  DBI::dbExecute(con, sprintf(paste(
    "CREATE OR REPLACE TEMP TABLE %s (sample_id VARCHAR, chrom VARCHAR, pos BIGINT,",
    "ref_count INTEGER, alt_count INTEGER, af DOUBLE)"), counts_table))
  parameters <- declaration$truth
  for (k in seq_len(nrow(children))) {
    child <- children[k, ]
    planted <- roh_truth_segments(positions, declaration)
    segments[[k]] <- data.frame(child = child$child, population = child$population, planted)
    first <- allele(child$founder_a, child$haplotype_a)
    second <- allele(child$founder_b, child$haplotype_b)
    inside <- rep(FALSE, sites)
    for (s in seq_len(nrow(planted))) {
      inside <- inside | (positions >= planted$start[[s]] & positions <= planted$end[[s]])
    }
    second[inside] <- first[inside]
    dosage <- first + second
    truth[, k] <- dosage

    # Genotype evidence: an error reports one of the two other genotypes.
    error <- stats::runif(sites) < declaration$gt_error_rate
    shift <- sample.int(2L, sites, replace = TRUE)
    reported <- ifelse(error, (dosage + shift) %% 3L, dosage)
    called[, k] <- reported

    # Read evidence.
    contaminant <- code(child$founder_c) %/% 2L + code(child$founder_c) %% 2L
    depth <- stats::rpois(sites, declaration$depth)
    artefact <- stats::runif(sites) < parameters[["artefact_weight"]]
    uniform_alt <- floor(stats::runif(sites) * (depth + 1))
    q <- c(parameters[["seq_error"]], parameters[["allele_balance"]], 1 - parameters[["seq_error"]])
    chance <- (1 - parameters[["contamination"]]) * q[dosage + 1L] +
      parameters[["contamination"]] * q[contaminant + 1L]
    binomial_alt <- stats::rbinom(sites, depth, chance)
    alt <- as.integer(ifelse(artefact, uniform_alt, binomial_alt))
    frequency <- panel[[paste0("af_", child$population)]]
    DBI::dbAppendTable(con, counts_table, data.frame(
      sample_id = child$child, chrom = panel$chrom, pos = positions,
      ref_count = as.integer(depth - alt), alt_count = alt, af = frequency,
      stringsAsFactors = FALSE))
    realized[[k]] <- data.frame(
      child = child$child, population = child$population,
      planted_sites = sum(inside), planted_heterozygous = sum(dosage[inside] == 1L),
      genotype_errors = sum(error), artefact_sites = sum(artefact),
      mean_depth = mean(depth), contaminant_heterozygous = mean(contaminant == 1L))
    heterozygous[[k]] <- data.frame(child = child$child, pos = positions[dosage == 1L])
    background[[k]] <- data.frame(child = child$child,
                                  roh_truth_background(positions, dosage, declaration$min_call_bases))
  }
  segments <- do.call(rbind, segments)
  realized <- do.call(rbind, realized)
  heterozygous <- do.call(rbind, heterozygous)
  background <- do.call(rbind, background)

  # The check of the plan, in SQL from the founder codes, independent of the
  # R construction above: every true genotype is the sum of the two drawn
  # haplotypes, with the second replaced by the first inside a segment.
  DBI::dbWriteTable(con, "truth_children", children, overwrite = TRUE)
  DBI::dbWriteTable(con, "truth_segments", segments, overwrite = TRUE)
  DBI::dbWriteTable(con, "truth_dosage", data.frame(
    child = rep(children$child, each = sites), site_index = rep(seq_len(sites) - 1L, nrow(children)),
    dosage = as.vector(truth)), overwrite = TRUE)
  DBI::dbExecute(con, paste(
    "CREATE OR REPLACE TEMP TABLE founder_codes AS SELECT site_index, pos,",
    "substr(founder, 6) AS founder, code FROM (UNPIVOT (SELECT site_index, pos, COLUMNS('^code_')",
    "FROM panel_founders) ON COLUMNS('^code_') INTO NAME founder VALUE code)"))
  mismatches <- DBI::dbGetQuery(con, paste(
    "WITH expected AS (SELECT c.child, a.site_index,",
    "CASE WHEN c.haplotype_a = 1 THEN a.code // 2 ELSE a.code % 2 END AS first_allele,",
    "CASE WHEN c.haplotype_b = 1 THEN b.code // 2 ELSE b.code % 2 END AS second_allele,",
    "EXISTS (SELECT 1 FROM truth_segments s WHERE s.child = c.child AND a.pos BETWEEN s.start AND s.\"end\") AS planted",
    "FROM truth_children c JOIN founder_codes a ON a.founder = c.founder_a",
    "JOIN founder_codes b ON b.founder = c.founder_b AND b.site_index = a.site_index)",
    "SELECT count(*) AS rows,",
    "count(*) FILTER (WHERE d.dosage <> CASE WHEN e.planted THEN 2 * e.first_allele",
    "ELSE e.first_allele + e.second_allele END) AS mismatches,",
    "count(*) FILTER (WHERE e.planted AND d.dosage = 1) AS planted_heterozygous",
    "FROM expected e JOIN truth_dosage d USING (child, site_index)"))
  if (mismatches$rows != sites * nrow(children) || mismatches$mismatches != 0 ||
      mismatches$planted_heterozygous != 0) {
    stop("the children differ from their plan: ", mismatches$mismatches, " mismatches, ",
         mismatches$planted_heterozygous, " heterozygous planted sites", call. = FALSE)
  }

  # The evidence files.
  counts_path <- file.path(work, "counts.parquet")
  DBI::dbExecute(con, sprintf("COPY (SELECT * FROM %s ORDER BY sample_id, pos) TO %s (FORMAT parquet)",
                              counts_table, roh_truth_quote(counts_path)))
  frequency_path <- file.path(work, "frequencies.parquet")
  frequencies <- do.call(rbind, lapply(populations, function(population) {
    data.frame(population = population, chrom = panel$chrom, pos = positions, ref = panel$ref,
               alt = panel$alt, af = panel[[paste0("af_", population)]], stringsAsFactors = FALSE)
  }))
  DBI::dbWriteTable(con, "truth_frequencies", frequencies, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY truth_frequencies TO %s (FORMAT parquet)",
                              roh_truth_quote(frequency_path)))
  vcf_path <- file.path(work, "children.vcf")
  genotypes <- matrix(roh_truth_genotype_strings[called + 1L], nrow = sites)
  body <- do.call(paste, c(list(panel$chrom, format(positions, scientific = FALSE, trim = TRUE), ".",
                                panel$ref, panel$alt, ".", "PASS", ".", "GT"),
                           as.data.frame(genotypes, stringsAsFactors = FALSE), sep = "\t"))
  # The panel is chr20 of GRCh38; the decodes do not read the contig length.
  if (!identical(panel$chrom[[1L]], "chr20")) stop("the panel must be chr20", call. = FALSE)
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr20,length=64444167>",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       children$child), collapse = "\t"),
               body), vcf_path)
  vcf_sha256 <- roh_truth_sha256(vcf_path)
  # BGZF with an index: read_bcf reports a row estimate (records times
  # samples) only for an indexed scan. Without one the planner took the call
  # rows for the small side of the frequency join and built its hash table on
  # them, about 16 MiB per sample here, where the declared model of the
  # genotype decode builds it on the frequency rows.
  status <- system2(RBCFTools::bgzip_path(), c("-f", shQuote(vcf_path)))
  if (status != 0L) stop("bgzip failed on the children VCF", call. = FALSE)
  compressed_path <- paste0(vcf_path, ".gz")
  DBI::dbGetQuery(con, sprintf("SELECT * FROM bcf_index(%s)", roh_truth_quote(compressed_path)))
  list(children = children, segments = segments, realized = realized, sites = sites,
       positions = positions, heterozygous = heterozygous, background = background,
       span_bases = positions[[sites]] - positions[[1L]] + 1,
       counts_path = counts_path, frequency_path = frequency_path, vcf_path = compressed_path,
       counts_sha256 = roh_truth_sha256(counts_path), vcf_sha256 = vcf_sha256)
}

# Disjoint intervals covering the union of two interval sets of one child.
roh_truth_union <- function(a, b) {
  intervals <- rbind(a[c("start", "end")], b[c("start", "end")])
  if (!nrow(intervals)) return(data.frame(start = numeric(), end = numeric()))
  intervals <- intervals[order(intervals$start, intervals$end), , drop = FALSE]
  starts <- numeric()
  ends <- numeric()
  for (k in seq_len(nrow(intervals))) {
    if (length(ends) && intervals$start[[k]] <= ends[[length(ends)]] + 1) {
      ends[[length(ends)]] <- max(ends[[length(ends)]], intervals$end[[k]])
    } else {
      starts <- c(starts, intervals$start[[k]])
      ends <- c(ends, intervals$end[[k]])
    }
  }
  data.frame(start = starts, end = ends)
}

# The true intervals of every child: version 1, the planted segments; version
# 2 (declaration_v2.R), the planted segments merged with the background
# stretches.
roh_truth_sets <- function(children, segments, background) {
  per_child <- function(version) {
    do.call(rbind, lapply(children$child, function(child) {
      planted <- segments[segments$child == child, c("start", "end"), drop = FALSE]
      intervals <- if (version == "v1") {
        planted
      } else {
        roh_truth_union(planted, background[background$child == child, c("start", "end"), drop = FALSE])
      }
      if (nrow(intervals)) data.frame(child = child, intervals)
    }))
  }
  list(v1 = per_child("v1"), v2 = per_child("v2"))
}

# Intervals (start, end; one-based, inclusive) in both of two disjoint sets.
roh_truth_intersect <- function(a, b) {
  if (!nrow(a) || !nrow(b)) return(data.frame(start = numeric(), end = numeric()))
  first <- outer(as.numeric(a$start), as.numeric(b$start), pmax)
  last <- outer(as.numeric(a$end), as.numeric(b$end), pmin)
  keep <- last >= first
  data.frame(start = first[keep], end = last[keep])
}

# Bases of each run class of every arm against the true intervals of a truth
# version (roh_truth_sets), per child.
roh_truth_score <- function(runs, arms, children, truth, min_call_bases) {
  rows <- list()
  for (arm in arms) {
    for (child in children$child) {
      true_set <- truth[truth$child == child, c("start", "end"), drop = FALSE]
      calls <- runs[runs$arm == arm & runs$sample == child, c("start", "end"), drop = FALSE]
      for (run_class in c("all", "long")) {
        kept <- if (run_class == "long") {
          calls[calls$end - calls$start + 1 >= min_call_bases, , drop = FALSE]
        } else {
          calls
        }
        called <- roh_run_bases(kept)
        true_bases <- roh_run_overlap_bases(kept, true_set)
        rows[[length(rows) + 1L]] <- data.frame(
          arm = arm, run_class = run_class, child = child,
          population = children$population[children$child == child],
          truth_bases = roh_run_bases(true_set), runs = nrow(kept), called_bases = called,
          true_bases = true_bases, false_bases = called - true_bases, stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

# The secondary measures of the declaration: the false called bases inside
# and outside the background stretches, and the sites inside calls.
roh_truth_secondary <- function(runs, arms, children, segments, background, heterozygous, positions,
                                min_call_bases) {
  rows <- list()
  for (child in children$child) {
    planted <- segments[segments$child == child, c("start", "end"), drop = FALSE]
    stretches <- background[background$child == child, c("start", "end"), drop = FALSE]
    child_heterozygous <- heterozygous$pos[heterozygous$child == child]
    outside_planted <- positions[!roh_truth_inside(positions, planted)]
    background_bases <- roh_run_bases(stretches) - roh_run_overlap_bases(stretches, planted)
    for (arm in arms) {
      calls <- runs[runs$arm == arm & runs$sample == child, c("start", "end"), drop = FALSE]
      for (run_class in c("all", "long")) {
        kept <- if (run_class == "long") {
          calls[calls$end - calls$start + 1 >= min_call_bases, , drop = FALSE]
        } else {
          calls
        }
        in_stretches <- roh_truth_intersect(kept, stretches)
        false_in_background <- roh_run_bases(in_stretches) -
          roh_run_overlap_bases(in_stretches, planted)
        false_bases <- roh_run_bases(kept) - roh_run_overlap_bases(kept, planted)
        rows[[length(rows) + 1L]] <- data.frame(
          arm = arm, run_class = run_class, child = child,
          population = children$population[children$child == child],
          background_bases = background_bases, false_bases = false_bases,
          false_in_background = false_in_background,
          false_elsewhere = false_bases - false_in_background,
          heterozygous_sites = length(child_heterozygous),
          heterozygous_in_calls = sum(roh_truth_inside(child_heterozygous, kept)),
          background_sites_in_calls = sum(roh_truth_inside(outside_planted, kept)),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

roh_truth_pooled <- function(scores) {
  sensitivity <- sum(scores$true_bases) / sum(scores$truth_bases)
  called <- sum(scores$called_bases)
  c(sensitivity = sensitivity,
    false_discovery = if (called > 0) sum(scores$false_bases) / called else NA_real_)
}

roh_truth_summary <- function(scores, by) {
  groups <- split(scores, scores[by], drop = TRUE)
  result <- do.call(rbind, lapply(groups, function(rows) {
    pooled <- roh_truth_pooled(rows)
    data.frame(rows[1L, by, drop = FALSE], children = nrow(rows), runs = sum(rows$runs),
               truth_bases = sum(rows$truth_bases), called_bases = sum(rows$called_bases),
               true_bases = sum(rows$true_bases), false_bases = sum(rows$false_bases),
               sensitivity = pooled[["sensitivity"]], false_discovery = pooled[["false_discovery"]],
               row.names = NULL)
  }))
  rownames(result) <- NULL
  result
}

# Paired bootstrap of the pooled difference test - reference, resampling
# children within each population.
roh_truth_compare <- function(scores, test, reference, run_class, declaration) {
  pick <- function(arm) {
    rows <- scores[scores$arm == arm & scores$run_class == run_class, , drop = FALSE]
    rows[order(rows$child), , drop = FALSE]
  }
  a <- pick(test)
  b <- pick(reference)
  stopifnot(identical(a$child, b$child))
  difference <- function(index) {
    roh_truth_pooled(a[index, ]) - roh_truth_pooled(b[index, ])
  }
  observed <- difference(seq_len(nrow(a)))
  strata <- split(seq_len(nrow(a)), a$population)
  roh_truth_rng(declaration$seed + 1L)
  draws <- t(vapply(seq_len(declaration$bootstrap_draws), function(draw) {
    index <- unlist(lapply(strata, function(members) members[sample.int(length(members), replace = TRUE)]),
                    use.names = FALSE)
    difference(index)
  }, numeric(2L)))
  tail <- (1 - declaration$bootstrap_level) / 2
  interval <- function(column, probability) {
    stats::quantile(draws[, column], probability, names = FALSE, na.rm = TRUE)
  }
  data.frame(test = test, reference = reference, run_class = run_class,
             metric = c("sensitivity", "false_discovery"),
             test_value = roh_truth_pooled(a), reference_value = roh_truth_pooled(b),
             difference = observed,
             lower = c(interval(1L, tail), interval(2L, tail)),
             upper = c(interval(1L, 1 - tail), interval(2L, 1 - tail)),
             row.names = NULL)
}

# The share of each planted segment that the long runs of each arm cover.
roh_truth_by_length <- function(runs, arms, segments, min_call_bases) {
  rows <- list()
  for (arm in arms) {
    long <- runs[runs$arm == arm & runs$end - runs$start + 1 >= min_call_bases, , drop = FALSE]
    covered <- vapply(seq_len(nrow(segments)), function(s) {
      calls <- long[long$sample == segments$child[[s]], c("start", "end"), drop = FALSE]
      roh_run_overlap_bases(calls, segments[s, c("start", "end")]) / segments$length[[s]]
    }, numeric(1L))
    for (length_bases in sort(unique(segments$length))) {
      share <- covered[segments$length == length_bases]
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm, length = length_bases, segments = length(share), mean_covered = mean(share),
        fully_covered = mean(share >= 0.99), missed = mean(share == 0))
    }
  }
  do.call(rbind, rows)
}

# Where the bases that stable_fraction and gained_bases count lie: in the
# truth (a true call, or a missed true base that a higher weight recovers) or
# outside it (a false call). `truth` is one version of roh_truth_sets;
# `central` holds the runs of counts_fit with both columns; `low` and `high`
# the runs at w - sd and w + sd.
roh_truth_stability <- function(central, low, high, truth, children) {
  per_child <- list()
  for (child in children$child) {
    true_set <- truth[truth$child == child, c("start", "end"), drop = FALSE]
    runs <- central[central$sample == child, , drop = FALSE]
    both <- roh_truth_intersect(low[low$sample == child, c("start", "end"), drop = FALSE],
                                high[high$sample == child, c("start", "end"), drop = FALSE])
    unstable <- 0
    unstable_true <- 0
    run_true <- 0
    for (r in seq_len(nrow(runs))) {
      run <- runs[r, c("start", "end")]
      stable_part <- roh_truth_intersect(run, both)
      located <- runs$length[[r]] - roh_run_bases(stable_part)
      reported <- runs$length[[r]] * (1 - runs$stable_fraction[[r]])
      if (abs(located - reported) >= 0.5) {
        stop("stable_fraction of a run of ", child, " differs from its located bases (",
             reported, " against ", located, ")", call. = FALSE)
      }
      in_truth <- roh_run_overlap_bases(run, true_set)
      unstable <- unstable + located
      unstable_true <- unstable_true + in_truth - roh_run_overlap_bases(stable_part, true_set)
      run_true <- run_true + in_truth
    }
    run_bases <- sum(as.numeric(runs$length))
    child_high <- high[high$sample == child, c("start", "end"), drop = FALSE]
    attached <- child_high[roh_fitted_meets(child_high, runs), , drop = FALSE]
    gained <- roh_fitted_gained_split(child_high, runs)
    if (sum(as.numeric(runs$gained_bases)) != gained[["attached"]]) {
      stop("gained_bases of ", child, " differ from the decode at w + sd", call. = FALSE)
    }
    gained_true <- roh_run_overlap_bases(attached, true_set) -
      roh_run_overlap_bases(roh_truth_intersect(attached, runs[c("start", "end")]), true_set)
    per_child[[child]] <- data.frame(
      child = child, population = children$population[children$child == child],
      run_bases = run_bases, run_false = run_bases - run_true,
      unstable_bases = unstable, unstable_false = unstable - unstable_true,
      stable_bases = run_bases - unstable,
      stable_false = (run_bases - run_true) - (unstable - unstable_true),
      missed_bases = roh_run_bases(true_set) - run_true,
      gained_bases = gained[["attached"]], gained_true = gained_true,
      gained_false = gained[["attached"]] - gained_true, new_run_bases = gained[["isolated"]])
  }
  result <- do.call(rbind, per_child)
  rownames(result) <- NULL
  result
}

roh_truth_run <- function(extension, cli, output_dir, size = c("full", "reduced"),
                          declaration = roh_truth_declaration) {
  size <- match.arg(size)
  populations <- if (size == "full") declaration$populations else declaration$reduced_populations
  children_per_population <- if (size == "full") declaration$children else declaration$reduced_children
  threads <- declaration$threads
  extension <- normalizePath(extension, mustWork = TRUE)
  work <- tempfile("roh-truth-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  con <- roh_truth_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  built <- roh_truth_build(con, work, declaration, populations, children_per_population)
  children <- built$children
  sites <- built$sites
  rows <- nrow(children) * sites
  q <- roh_truth_quote
  denominators <- DBI::dbGetQuery(con, sprintf(paste(
    "WITH r AS (SELECT sample_id, chrom, (pos - 1) // 10000000 AS block_index,",
    "floor(af * 64) AS bin, ref_count + alt_count AS depth, alt_count AS alt FROM read_parquet(%s)",
    "WHERE af > 0 AND af < 1 AND ref_count + alt_count BETWEEN 1 AND 1000)",
    "SELECT (SELECT count(*) FROM (SELECT DISTINCT sample_id, chrom, block_index, bin, depth, alt FROM r)) AS cells,",
    "(SELECT count(*) FROM (SELECT DISTINCT sample_id, depth, alt FROM r)) AS pairs"),
    q(built$counts_path)))

  recorder <- roh_truth_recorder(cli, extension, work)
  overhead_mib <- recorder$overhead(threads)
  ceiling <- roh_truth_ceiling_bytes
  counts_setup <- sprintf("CREATE TEMP VIEW counts AS SELECT * FROM read_parquet(%s);", q(built$counts_path))
  runs_columns <- "sample, chrom, start, \"end\", length, n_markers, quality"
  all_runs <- list()
  keep_runs <- function(arm, path) {
    runs <- roh_truth_read_runs(path)
    if (nrow(runs)) all_runs[[length(all_runs) + 1L]] <<- data.frame(arm = arm, runs)
    runs
  }

  # Genotype arms, one query per population and gt_error.
  for (arm in names(declaration$gt_errors)) {
    for (population in populations) {
      members <- children$child[children$population == population]
      output <- file.path(work, sprintf("%s_%s.csv", arm, population))
      recorder$record(sprintf("%s_%s", arm, population), "genotype_decode",
        sprintf(paste("CREATE TEMP VIEW population_af AS SELECT chrom, pos, ref, alt, af",
                      "FROM read_parquet(%s) WHERE population = %s;"),
                q(built$frequency_path), q(population)),
        sprintf(paste("SELECT %s FROM duckhts_roh_af_table(%s, 'population_af', gt_error := %s,",
                      "samples := %s) ORDER BY sample, chrom, start"), runs_columns, q(built$vcf_path),
                format(declaration$gt_errors[[arm]]), q(paste(members, collapse = ","))),
        threads, output, live_bytes = ceiling$genotype_decode(length(members), sites, threads),
        overhead_mib = overhead_mib)
      keep_runs(arm, output)
    }
  }

  # The read-count arms.
  output <- file.path(work, "counts_default.csv")
  recorder$record("counts_default", "count_decode", counts_setup,
    sprintf("SELECT %s FROM duckhts_roh_counts('counts') ORDER BY sample, chrom, start", runs_columns),
    threads, output, live_bytes = ceiling$count_decode(rows, sites, threads), overhead_mib = overhead_mib)
  keep_runs("counts_default", output)

  fit_path <- file.path(work, "fit.parquet")
  recorder$record("fit", "fit", counts_setup,
    "SELECT * FROM duckhts_count_error_fit('counts') ORDER BY sample", threads, fit_path, "parquet",
    live_bytes = ceiling$fit(denominators$cells, denominators$pairs, nrow(children)),
    overhead_mib = overhead_mib)
  fit <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) ORDER BY sample", q(fit_path)))
  if (!setequal(fit$sample, children$child) || nrow(fit) != nrow(children)) {
    stop("the fit does not have one row per child", call. = FALSE)
  }
  numbers <- c("seq_error", "contamination", "allele_balance", "artefact_weight")
  usable <- fit$sample[stats::complete.cases(fit[numbers])]

  # One decode of duckhts_roh_counts_fitted, with a fit relation of exactly
  # the numbers in `parameters` for the children it names.
  decode_fitted <- function(name, parameters, stability = FALSE) {
    fit_file <- file.path(work, sprintf("fit_%s.parquet", name))
    DBI::dbWriteTable(con, "decode_fit", parameters[c("sample", numbers, "artefact_weight_sd")],
                      overwrite = TRUE)
    DBI::dbExecute(con, sprintf("COPY decode_fit TO %s (FORMAT parquet)", q(fit_file)))
    output <- file.path(work, sprintf("%s.csv", name))
    columns <- if (stability) paste(runs_columns, "stable_fraction, gained_bases", sep = ", ") else runs_columns
    decoded <- parameters$sample
    recorder$record(name, "count_decode", c(counts_setup,
      sprintf("CREATE TEMP VIEW decoded_counts AS SELECT * FROM counts WHERE sample_id IN (%s);",
              paste(q(decoded), collapse = ", ")),
      sprintf("CREATE TEMP VIEW decode_fit AS SELECT * FROM read_parquet(%s);", q(fit_file))),
      sprintf(paste("SELECT %s FROM duckhts_roh_counts_fitted('decoded_counts', 'decode_fit')",
                    "ORDER BY sample, chrom, start"), columns), threads, output,
      live_bytes = ceiling$count_decode(length(decoded) * sites, sites, threads), overhead_mib = overhead_mib)
    output
  }
  truth_parameters <- data.frame(sample = children$child, t(declaration$truth), artefact_weight_sd = 0)
  keep_runs("counts_truth", decode_fitted("counts_truth", truth_parameters))
  fitted_parameters <- fit[fit$sample %in% usable, c("sample", numbers, "artefact_weight_sd")]
  central <- roh_truth_read_runs(decode_fitted("counts_fit", fitted_parameters, stability = TRUE))
  if (nrow(central)) {
    all_runs[[length(all_runs) + 1L]] <- data.frame(
      arm = "counts_fit", central[setdiff(names(central), c("stable_fraction", "gained_bases"))])
  }
  guard <- c(0, 0.05)
  sd <- ifelse(is.na(fitted_parameters$artefact_weight_sd), 0, fitted_parameters$artefact_weight_sd)
  shifted <- function(direction) {
    parameters <- fitted_parameters
    parameters$artefact_weight <- pmin(guard[[2L]], pmax(guard[[1L]],
                                                         parameters$artefact_weight + direction * sd))
    parameters$artefact_weight_sd <- 0
    parameters
  }
  low <- roh_truth_read_runs(decode_fitted("fit_low", shifted(-1)))
  high <- roh_truth_read_runs(decode_fitted("fit_high", shifted(1)))

  observations <- recorder$table()
  roh_truth_gate(observations, file.path(output_dir, "truth_observations.csv"))
  runs <- do.call(rbind, all_runs)
  # Both truths: version 1 (declaration.R) and version 2 (declaration_v2.R),
  # whose background stretches have the length bound of the secondary
  # measures of version 1.
  stopifnot(roh_truth_declaration_v2$min_stretch_bases == declaration$min_call_bases)
  truths <- roh_truth_sets(children, built$segments, built$background)
  versions <- lapply(names(truths), function(version) {
    scores <- roh_truth_score(runs, declaration$arms, children, truths[[version]],
                              declaration$min_call_bases)
    comparisons <- do.call(rbind, lapply(seq_len(nrow(declaration$comparisons)), function(k) {
      row <- declaration$comparisons[k, ]
      do.call(rbind, lapply(c("long", "all"), function(run_class) {
        data.frame(truth = version, criterion = row$criterion,
                   roh_truth_compare(scores, row$test, row$reference, run_class, declaration))
      }))
    }))
    long_rows <- comparisons[comparisons$run_class == "long", ]
    t1 <- long_rows[long_rows$criterion == "T1", ]
    t2 <- long_rows[long_rows$criterion == "T2", ]
    t1_pass <- t1$lower[t1$metric == "sensitivity"] > 0 &&
      t1$upper[t1$metric == "false_discovery"] < declaration$t1_max_fdr_increase
    t2_pass <- all(t2$lower >= -declaration$t2_margin & t2$upper <= declaration$t2_margin)
    list(scores = data.frame(truth = version, scores),
         pooled = data.frame(truth = version, roh_truth_summary(scores, c("arm", "run_class"))),
         by_population = data.frame(truth = version,
                                    roh_truth_summary(scores, c("arm", "run_class", "population"))),
         comparisons = comparisons,
         criteria = data.frame(truth = version, criterion = c("T1", "T2"), pass = c(t1_pass, t2_pass)))
  })
  gather <- function(part) do.call(rbind, lapply(versions, `[[`, part))
  scores <- gather("scores")
  pooled <- gather("pooled")
  by_population <- gather("by_population")
  comparisons <- gather("comparisons")

  fit_with_population <- merge(fit, children[c("child", "population")], by.x = "sample", by.y = "child")
  t3 <- do.call(rbind, lapply(split(fit_with_population, fit_with_population$population), function(rows) {
    medians <- vapply(numbers, function(number) stats::median(rows[[number]], na.rm = TRUE), numeric(1L))
    data.frame(population = rows$population[[1L]], number = numbers,
               simulated = declaration$truth[numbers], median = medians,
               tolerance = declaration$t3_tolerance[numbers],
               pass = abs(medians - declaration$truth[numbers]) <= declaration$t3_tolerance[numbers],
               row.names = NULL)
  }))
  ok_share <- mean(fit$status == "ok")
  t3_pass <- all(t3$pass) && ok_share >= declaration$t3_min_ok
  criteria <- rbind(gather("criteria"), data.frame(truth = "either", criterion = "T3", pass = t3_pass))

  secondary <- roh_truth_secondary(runs, declaration$arms, children, built$segments, built$background,
                                   built$heterozygous, built$positions, declaration$min_call_bases)
  secondary_pooled <- do.call(rbind, lapply(split(secondary, secondary[c("arm", "run_class")], drop = TRUE),
    function(rows) {
      data.frame(arm = rows$arm[[1L]], run_class = rows$run_class[[1L]],
                 background_bases = sum(rows$background_bases), false_bases = sum(rows$false_bases),
                 false_in_background = sum(rows$false_in_background),
                 false_elsewhere = sum(rows$false_elsewhere),
                 heterozygous_sites = sum(rows$heterozygous_sites),
                 heterozygous_in_calls = sum(rows$heterozygous_in_calls),
                 background_sites_in_calls = sum(rows$background_sites_in_calls))
    }))
  rownames(secondary_pooled) <- NULL
  by_length <- roh_truth_by_length(runs, declaration$arms, built$segments, declaration$min_call_bases)
  decoded_children <- children[children$child %in% usable, , drop = FALSE]
  stability <- do.call(rbind, lapply(names(truths), function(version) {
    data.frame(truth = version, roh_truth_stability(central, low, high, truths[[version]], decoded_children))
  }))

  metadata <- data.frame(
    key = c("revision", "src_tree", "extension_sha256", "duckdb_cli", "r_version", "host_load",
            "cpu_affinity", "threads", "size", "populations", "children", "sites", "span_bases",
            "counts_rows", "fit_cells", "fit_pairs", "overhead_mib", "counts_sha256", "vcf_sha256",
            "children_with_usable_fit"),
    value = c(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
              system2("git", c("rev-parse", "HEAD:src"), stdout = TRUE),
              roh_truth_sha256(extension), system2(cli, "--version", stdout = TRUE)[[1L]],
              R.version.string, sub(".*load average: ", "", system2("uptime", stdout = TRUE)),
              sub("^.*: *", "", system2("taskset", c("-pc", Sys.getpid()), stdout = TRUE)),
              threads, size, paste(populations, collapse = ","), nrow(children), sites,
              format(built$span_bases, scientific = FALSE), format(rows, scientific = FALSE),
              format(denominators$cells), format(denominators$pairs), sprintf("%.1f", overhead_mib),
              built$counts_sha256, built$vcf_sha256, length(usable)),
    stringsAsFactors = FALSE)
  write <- function(data, name) utils::write.csv(data, file.path(output_dir, name), row.names = FALSE)
  write(metadata, "truth_metadata.csv")
  write(children, "truth_children.csv")
  write(built$segments, "truth_segments.csv")
  write(truths$v2, "truth_v2_intervals.csv")
  write(built$realized, "truth_realized.csv")
  write(fit, "truth_fit.csv")
  # The two run tables are the large outputs; they are written compressed.
  write_compressed <- function(data, name) {
    connection <- gzfile(file.path(output_dir, name), "w")
    on.exit(close(connection))
    utils::write.csv(data, connection, row.names = FALSE)
  }
  write_compressed(runs, "truth_runs.csv.gz")
  write_compressed(central, "truth_counts_fit_runs.csv.gz")
  write(scores, "truth_child_scores.csv")
  write(pooled, "truth_pooled.csv")
  write(by_population, "truth_by_population.csv")
  write(comparisons, "truth_comparisons.csv")
  write(by_length, "truth_by_length.csv")
  write(built$background, "truth_background_stretches.csv")
  write(secondary, "truth_child_secondary.csv")
  write(secondary_pooled, "truth_secondary.csv")
  write(t3, "truth_fit_recovery.csv")
  write(criteria, "truth_criteria.csv")
  write(stability, "truth_stability.csv")
  write(observations, "truth_observations.csv")
  list(metadata = metadata, children = children, segments = built$segments, realized = built$realized,
       fit = fit, pooled = pooled, by_population = by_population, comparisons = comparisons,
       by_length = by_length, secondary = secondary_pooled, background = built$background,
       t3 = t3, ok_share = ok_share, criteria = criteria,
       stability = stability, observations = observations, span_bases = built$span_bases,
       runs = runs, counts_fit_runs = central)
}
