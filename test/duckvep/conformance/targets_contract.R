#!/usr/bin/env Rscript

# Focused contracts for the optional targets/blit campaign adapter. This does
# not execute VEP or require the targets package; the checked witness pipeline
# supplies that end-to-end test.

root <- normalizePath(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)[[1L]],
  mustWork = TRUE
)
source(file.path(root, "scripts", "duckvep_evidence.R"), local = TRUE)
source(
  file.path(root, "pipelines", "duckvep", "R", "pipeline.R"),
  local = TRUE
)

expect_error <- function(expression, pattern) {
  message <- tryCatch(
    {
      force(expression)
      ""
    },
    error = conditionMessage
  )
  stopifnot(nzchar(message), grepl(pattern, message))
}

temporary_root <- tempfile("duckvep-targets-contract-")
dir.create(temporary_root)
on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)

micromamba <- file.path(temporary_root, "micromamba")
writeLines("fixture", micromamba)
vep_prefix <- file.path(temporary_root, "vep")
dir.create(file.path(vep_prefix, "conda-meta"), recursive = TRUE)
writeLines("fixture", file.path(vep_prefix, "conda-meta", "history"))
extension <- file.path(temporary_root, "duckhts.duckdb_extension")
extension_receipt <- file.path(temporary_root, "extension.tsv")
invisible(lapply(c(extension, extension_receipt), writeLines, text = "fixture"))
manifest <- file.path(temporary_root, "campaigns.tsv")
campaign <- data.frame(
  id = "contract",
  enabled = "true",
  corpus = "contract",
  event_mode = "small",
  gff = file.path(root, "test", "data", "duckvep", "minimal.gff3"),
  gff_index_policy = "ignore",
  fasta = file.path(root, "test", "data", "duckvep", "minimal.fa"),
  database = ":memory:",
  model_sql = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "minimal_model.sql"
  ),
  model_name = "contract",
  hgvs = "true",
  fork = "1",
  duckdb_threads = "1",
  oracle_environment_receipt = file.path(
    root,
    "test",
    "duckvep",
    "upstream",
    "receipts",
    "vep116_2026-07-22.conda-explicit.txt"
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  campaign,
  manifest,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

output_root <- file.path(temporary_root, "results")
targets_store <- file.path(temporary_root, "store")
campaigns <- duckvep_targets_read_campaigns(
  manifest,
  root,
  micromamba,
  vep_prefix,
  output_root,
  c(extension, extension_receipt),
  targets_store
)
stopifnot(
  length(campaigns) == 1L,
  length(campaigns[[1L]]$outputs) == 10L,
  all(file.exists(duckvep_targets_campaign_inputs(campaigns[[1L]])))
)

traversal <- file.path(output_root, "..", basename(targets_store))
expect_error(
  duckvep_targets_read_campaigns(
    manifest,
    root,
    micromamba,
    vep_prefix,
    traversal,
    c(extension, extension_receipt),
    targets_store
  ),
  "must not overlap"
)
if (!identical(.Platform$OS.type, "windows")) {
  stopifnot(identical(
    duckvep_targets_path_key(c("/tmp/A.fa", "/tmp/a.fa")),
    c("/tmp/A.fa", "/tmp/a.fa")
  ))
  dir.create(targets_store)
  symlink_output <- file.path(temporary_root, "results-link")
  stopifnot(file.symlink(targets_store, symlink_output))
  expect_error(
    duckvep_targets_read_campaigns(
      manifest,
      root,
      micromamba,
      vep_prefix,
      symlink_output,
      c(extension, extension_receipt),
      targets_store
    ),
    "must not overlap"
  )
  symlink_fixture <- file.path(temporary_root, "symlink-fixture")
  dir.create(file.path(symlink_fixture, "a"), recursive = TRUE)
  dir.create(file.path(symlink_fixture, "b", "deep"), recursive = TRUE)
  stopifnot(file.symlink(
    file.path(symlink_fixture, "b", "deep"),
    file.path(symlink_fixture, "a", "link")
  ))
  stopifnot(identical(
    duckvep_targets_canonical_path(file.path(
      symlink_fixture,
      "a",
      "link",
      "..",
      "model"
    )),
    file.path(symlink_fixture, "b", "model")
  ))
}

dangerous_output <- file.path(temporary_root, "dangerous-output")
dangerous_manifest <- file.path(dangerous_output, "contract", "campaigns.tsv")
dir.create(dirname(dangerous_manifest), recursive = TRUE)
stopifnot(file.copy(manifest, dangerous_manifest))
expect_error(
  duckvep_targets_read_campaigns(
    dangerous_manifest,
    root,
    micromamba,
    vep_prefix,
    dangerous_output,
    c(extension, extension_receipt),
    targets_store
  ),
  "campaign directory contains"
)

windows_space <- duckvep_blit_quote(
  "C:/Program Files/R/Rscript.exe",
  os_type = "windows"
)
windows_meta <- duckvep_blit_quote("x&y", os_type = "windows")
stopifnot(
  startsWith(windows_space, "^\""),
  endsWith(windows_space, "^\""),
  grepl("^&", windows_meta, fixed = TRUE)
)

cache_root <- file.path(temporary_root, "cache")
cache_info <- file.path(cache_root, "homo_sapiens", "116_GRCh38", "info.txt")
dir.create(dirname(cache_info), recursive = TRUE)
writeLines(c("species\thomo_sapiens", "assembly\tGRCh38"), cache_info)
decoy_info <- file.path(
  cache_root,
  "decoy",
  "homo_sapiens",
  "116_GRCh38",
  "info.txt"
)
dir.create(dirname(decoy_info), recursive = TRUE)
writeLines(c("species\thomo_sapiens", "assembly\tGRCh38"), decoy_info)
stopifnot(
  identical(
    duckvep_evidence_cache_info_path(
      cache_root,
      "homo_sapiens",
      "116",
      "GRCh38"
    ),
    normalizePath(cache_info)
  ),
  !identical(normalizePath(decoy_info), normalizePath(cache_info))
)

lock_records <- c(
  "# platform: linux-64",
  "@EXPLICIT",
  "https://example.org/a.conda",
  "file:///tmp/local.conda"
)
installed_records <- c(
  "List of packages in environment: /tmp/vep",
  "",
  "file:///tmp/local.conda",
  "https://example.org/a.conda"
)
stopifnot(identical(
  duckvep_evidence_explicit_packages(lock_records),
  duckvep_evidence_explicit_packages(installed_records)
))
expect_error(
  duckvep_evidence_explicit_packages(c(installed_records, "not-a-url")),
  "invalid non-URL"
)

final_directory <- dirname(campaigns[[1L]]$outputs[[1L]])
dir.create(dirname(final_directory), recursive = TRUE)
dir.create(final_directory)
writeLines("old", file.path(final_directory, "old"))
stage_directory <- tempfile(".run-contract-", tmpdir = dirname(final_directory))
dir.create(stage_directory)
writeLines("new", file.path(stage_directory, "new"))
duckvep_targets_publish_outputs(stage_directory, final_directory)
stopifnot(
  !file.exists(file.path(final_directory, "old")),
  identical(readLines(file.path(final_directory, "new")), "new")
)

# Simulate termination after moving the old complete directory aside.
backup <- tempfile(
  paste0(".previous-", basename(final_directory), "-"),
  tmpdir = dirname(final_directory)
)
stopifnot(file.rename(final_directory, backup))
duckvep_targets_recover_campaign_directory(final_directory)
stopifnot(dir.exists(final_directory), !dir.exists(backup))

# Simulate termination after publishing the complete replacement but before
# removing the old directory.
backup <- tempfile(
  paste0(".previous-", basename(final_directory), "-"),
  tmpdir = dirname(final_directory)
)
dir.create(backup)
writeLines("old", file.path(backup, "old"))
duckvep_targets_recover_campaign_directory(final_directory)
stopifnot(
  dir.exists(final_directory),
  !dir.exists(backup),
  identical(readLines(file.path(final_directory, "new")), "new")
)

cat("DuckVEP targets/blit contracts: ok\n")
