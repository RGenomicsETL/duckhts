# Stage the synthetic chr20 children and their planted truth from the registry.
source("benchmarks/roh_synthetic_stage.R")
staged <- stage_roh_synthetic_from_registry("roh_synthetic_chr20_children_bcf",
                                            "roh_synthetic_chr20_truth")
cat(sprintf("records=%.0f;samples=%.0f;records_sha256=%s\nrows=%.0f;sha256=%s\n",
            staged$records, staged$samples, staged$records_sha256,
            staged$truth_rows, staged$truth_sha256))
