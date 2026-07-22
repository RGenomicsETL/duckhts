library(targets)

duckvep_root <- normalizePath(
  Sys.getenv("DUCKHTS_ROOT", unset = "."),
  mustWork = TRUE
)
if (!file.exists(file.path(duckvep_root, "ARCHITECTURE.md"))) {
  stop("run the DuckVEP targets workflow from the DuckHTS repository root")
}
source(file.path(duckvep_root, "scripts", "duckvep_evidence.R"))
source(file.path(duckvep_root, "pipelines", "duckvep", "R", "pipeline.R"))

campaign_manifest <- Sys.getenv(
  "DUCKVEP_TARGETS_CAMPAIGNS",
  unset = file.path(duckvep_root, "pipelines", "duckvep", "campaigns.tsv")
)
campaign_manifest <- duckvep_targets_path(campaign_manifest, duckvep_root)
output_root <- Sys.getenv(
  "DUCKVEP_TARGETS_OUTPUT",
  unset = file.path(duckvep_root, ".tmp", "duckvep_targets", "results")
)
output_root <- duckvep_targets_path(output_root, duckvep_root)
targets_store <- Sys.getenv(
  "DUCKVEP_TARGETS_STORE",
  unset = file.path(duckvep_root, ".tmp", "duckvep_targets", "store")
)
targets_store <- duckvep_targets_path(targets_store, duckvep_root)
extension_receipt <- file.path(
  duckvep_root,
  ".tmp",
  "duckvep_targets",
  "extension.tsv"
)
default_micromamba <- Sys.getenv("MICROMAMBA", unset = "micromamba")
default_vep_prefix <- Sys.getenv(
  "VEP_PREFIX",
  unset = "/root/miniconda3/envs/vep"
)

tar_option_set(
  packages = "blit",
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue",
  workspace_on_error = TRUE
)

list(
  tar_target(
    duckvep_source_revision,
    duckvep_targets_source_revision(duckvep_root),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    duckvep_campaign_manifest,
    normalizePath(campaign_manifest, mustWork = TRUE),
    format = "file"
  ),
  tar_target(
    duckvep_runtime_defaults,
    list(
      micromamba = default_micromamba,
      vep_prefix = default_vep_prefix
    )
  ),
  tar_target(
    duckvep_campaigns,
    duckvep_targets_read_campaigns(
      duckvep_campaign_manifest,
      duckvep_root,
      duckvep_runtime_defaults$micromamba,
      duckvep_runtime_defaults$vep_prefix,
      output_root,
      c(
        file.path(
          duckvep_root,
          "build",
          "release",
          "duckhts.duckdb_extension"
        ),
        extension_receipt
      ),
      targets_store
    ),
    iteration = "list"
  ),
  tar_target(
    duckvep_extension_bundle,
    duckvep_targets_build_extension(
      duckvep_root,
      duckvep_source_revision,
      extension_receipt
    ),
    format = "file"
  ),
  tar_target(
    duckvep_campaign,
    duckvep_campaigns,
    pattern = map(duckvep_campaigns),
    iteration = "list"
  ),
  tar_target(
    duckvep_campaign_inputs,
    duckvep_targets_campaign_inputs(duckvep_campaign),
    pattern = map(duckvep_campaign),
    format = "file"
  ),
  tar_target(
    duckvep_campaign_evidence,
    duckvep_targets_run_campaign(
      duckvep_campaign,
      duckvep_campaign_inputs,
      duckvep_extension_bundle,
      duckvep_root
    ),
    pattern = map(duckvep_campaign, duckvep_campaign_inputs),
    format = "file"
  )
)
