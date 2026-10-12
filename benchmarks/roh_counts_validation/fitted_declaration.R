# The declared design of the fitted read-count emission check: decodes,
# perturbations, criteria and bounds. benchmark_roh_fitted_emission.Rmd states
# the same values in prose and reads them from here. They are fixed before any
# decode with a fitted allele balance or artefact weight is run.
#
# The inputs, run classes and overlap metric are those of
# roh_validation_declaration (declaration.R).
roh_fitted_declaration <- list(
  # The genotype decodes that the count decodes are compared with
  # (duckhts_roh_af_table, gt_error). The criterion uses the first.
  gt_errors = c(30, 100),
  # The fit: duckhts_count_error_fit with its default options. A sample whose
  # status is not 'ok' is not decoded with its fit.
  fit_status = "ok",
  # Count decodes, all through duckhts_roh_counts_fitted. Each names the source
  # of its four numbers: "default" is the scalar default of
  # duckhts_roh_counts, "fit" is the sample's fitted value.
  decodes = data.frame(
    decode = c("a", "b", "c"),
    label = c("scalar defaults", "fitted error and contamination", "full fit"),
    seq_error = c("default", "fit", "fit"),
    contamination = c("default", "fit", "fit"),
    allele_balance = c("default", "default", "fit"),
    artefact_weight = c("default", "default", "fit"),
    stringsAsFactors = FALSE),
  defaults = c(seq_error = 0.001, contamination = 0, allele_balance = 0.5,
               artefact_weight = 0),
  # Criterion: on every sample, decode c is at least as close to the genotype
  # decode at this gt_error as decodes a and b, by the absolute difference of
  # bases in runs (all runs), which is the absolute FROH difference times the
  # span.
  criterion_gt_error = 30,
  criterion_against = c("a", "b"),
  # Perturbations of decode c, one number at a time. "add" adds the value,
  # "scale" multiplies by it, and "add_sd" adds the value times the sample's
  # contamination_sd, or times sd_fallback when contamination_sd is NULL.
  perturbations = data.frame(
    parameter = c(rep("allele_balance", 4L), rep("artefact_weight", 2L),
                  rep("seq_error", 2L), rep("contamination", 2L)),
    operation = c(rep("add", 4L), rep("scale", 4L), rep("add_sd", 2L)),
    value = c(-0.01, -0.005, 0.005, 0.01, 0.5, 2, 0.5, 2, -1, 1),
    stringsAsFactors = FALSE),
  sd_fallback = 0.001,
  # A perturbed value outside the range duckhts_roh_counts accepts is set to
  # the end of that range and the row is marked. Halving or doubling a
  # seq_error in (0, 0.25) stays in (0, 0.5), so only the ends listed can act.
  guard_ranges = list(seq_error = c(0, 0.5), contamination = c(0, 1),
                      allele_balance = c(0.35, 0.65), artefact_weight = c(0, 0.05)),
  # Sensitivity bound, for every sample and perturbation, against decode c on
  # all runs: the absolute FROH change is below max_froh_change, and the bases
  # in the runs of exactly one of the two decodes, divided by the bases in the
  # runs of decode c, are below max_changed_share. When decode c has no run in
  # a sample, the share is 0 if the perturbed decode has none either and fails
  # otherwise.
  max_froh_change = 0.005,
  max_changed_share = 0.02,
  # Stability, declared after the first render and before any decode with an
  # artefact_weight_sd. Decode `decode` passes each sample's fitted
  # artefact_weight_sd to duckhts_roh_counts_fitted (a NULL, from a fit with
  # fewer than three blocks, counts as 0 there), so its runs carry
  # stable_fraction; every other count decode, perturbed ones included,
  # passes 0. Reported per sample, with no verdict:
  #   - artefact_weight_sd, and w - sd and w + sd (kept in [0, 0.05]) beside
  #     the declared x 0.5 and x 2 of the artefact_weight perturbations;
  #   - the stable_fraction of the runs of `decode`: the runs, the runs at 1,
  #     the quantiles below, and the stable bases (run length times
  #     stable_fraction, summed) over the run bases;
  #   - for each perturbation of `perturbed_parameter`: the lost bases, those
  #     in runs of `decode` that the perturbed decode does not call, split by
  #     the stable_fraction of the run that holds them (1, or below 1); and the
  #     gained bases, called by the perturbed decode outside any run of
  #     `decode`, which no stable_fraction covers. The FROH change is the
  #     gained minus the lost bases over the span, and the changed bases of
  #     the sensitivity rows are the gained plus the lost bases. The report
  #     gives the share of the lost bases, and of the changed bases, that lies
  #     in runs with stable_fraction below 1.
  stability = list(decode = "c", perturbed_parameter = "artefact_weight",
                   quantiles = c(0, 0.25, 0.5, 0.75, 1)),
  # DuckDB threads of every decode and of the fit.
  threads = 2L)
