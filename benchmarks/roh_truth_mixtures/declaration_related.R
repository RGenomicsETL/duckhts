# Read-level parent/child mixture declaration. Commit before running this arm.
# Receiver HG00403 and donor HG00405 are parent and child in the registered
# 1000 Genomes pedigree (SH001, CHS): HG00405's parents are HG00403/HG00404.
# The 2504-sample release has the two parents; the 698-related-sample release
# adds HG00405 (ERR3988761), sequenced at 30x by NYGC. Its ENA HTTPS CRAM is
# range-readable; the corresponding 2504 S3 data key is absent.
#
# The existing unrelated arm, declaration, inputs and criteria stay separate.
# Use its alpha grid, read-name hash seed (mates together), counting settings,
# frequencies and decoder. Include alpha 0 as the receiver-only control.
#
# M3: at each nonzero alpha, contamination_relative recovers alpha within
# max(0.005, 0.25 * alpha). At alpha >= 0.05 the relative log-likelihood must
# exceed the unrelated log-likelihood. The earlier count-level result separated
# the hypotheses at 5%, not at 2%; at 1%/2% no likelihood sign is required.
# Report those signs rather than treating a small preference as identification.
#
# M4: the unrelated model's contamination recovers alpha / 2 within
# max(0.003, 0.25 * alpha / 2), testing the catalog's parent/child interpretation.
# These tolerances allow read sampling and nuisance fits, not a claim of a
# confidence interval. Report the realised donor share and fit status too:
# equal sampling fractions need not give equal shares at unequal source depths.
#
# R1: |FROH erosion with the relative-fit term| <= |erosion without the term|
# + 0.002 at each nonzero alpha. Erosion is measured from alpha 0's with-term
# decode. Supply contamination_relative in the decoder's contamination field;
# retain the other parameters supplied by duckhts_count_error_fit, whose public
# output exposes the unrelated fit's nuisance estimates. This plug-in decode
# uses the existing population-frequency emission, not a parent/child emission.
#
# Apply the original input identity, no-spill and overhead + 3 x live-state
# memory gates. This fixed chr20 arm is not a 1x/2x/4x scaling experiment.
roh_related_declaration <- list(
  receiver = "HG00403",
  donor = "HG00405",
  slice_artifact = "roh_mixture_hg00405_chr20_cram",
  alphas = c(0.01, 0.02, 0.05, 0.10),
  m3_absolute_floor = 0.005,
  m3_relative_tolerance = 0.25,
  likelihood_min_alpha = 0.05,
  m4_absolute_floor = 0.003,
  m4_relative_tolerance = 0.25,
  erosion_tolerance = 0.002)
