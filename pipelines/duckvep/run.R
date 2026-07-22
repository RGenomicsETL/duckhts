#!/usr/bin/env Rscript

root <- normalizePath(
  Sys.getenv("DUCKHTS_ROOT", unset = "."),
  mustWork = TRUE
)
script <- file.path(root, "pipelines", "duckvep", "_targets.R")
store <- Sys.getenv(
  "DUCKVEP_TARGETS_STORE",
  unset = file.path(root, ".tmp", "duckvep_targets", "store")
)

targets::tar_make(script = script, store = store)
errored <- targets::tar_errored(store = store)
if (length(errored)) {
  stop(
    "DuckVEP targets failed: ",
    paste(errored, collapse = ", "),
    call. = FALSE
  )
}
