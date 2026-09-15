import math
import tables
import strutils

type allele_count = object
  nref: uint32
  nalt: uint32
  nother: uint32

include "upstream_charr_helpers.nim"
include "upstream_pair_helpers.nim"

proc emit(metric, input: string, value: float64) =
  echo metric, "\t", input, "\t", formatFloat(value, ffScientific, 17)

echo "metric\tinput\tvalue"
for depth in [100, 200, 6000]:
  emit("charr_threshold", $depth,
    max_hom_minor_reads(depth, 0.12, 0.002).float64)

let deep_counts = @[allele_count(nref: 4000, nalt: 2000, nother: 0)]
let deep_result = estimate_charr(deep_counts, @[0.5'f32], 15, 0.12, 0.002)
emit("charr_usable", "4000/2000/0", deep_result.n_sites_usable.float64)
emit("charr_estimate", "4000/2000/0", deep_result.contamination)

let control = estimate_charr(@[allele_count(nref: 199, nalt: 1, nother: 0)],
  @[0.25'f32], 15, 0.12, 0.002)
emit("charr_estimate", "199/1/0", control.contamination)

let priors = compute_log_gt_priors(@[0.5'f32, 0.5'f32])
let sites = @[
  PairSiteObservation(receiver_nref: 160, receiver_nalt: 40,
    anchor_gt: 0, log_contam_priors: priors[0]),
  PairSiteObservation(receiver_nref: 40, receiver_nalt: 160,
    anchor_gt: 2, log_contam_priors: priors[1])
]
let pair_result = estimate_pair_contamination_from_sites(sites, 0.002)
emit("pair_usable", "two_sites", pair_result.n_sites_usable.float64)
emit("pair_alpha", "upstream_search", pair_result.contamination)
emit("pair_log_likelihood", "upstream_search",
  pair_total_log_likelihood(sites, pair_result.contamination, 0.002))
for alpha in [0.0, 0.2, 0.4, 0.5, 1.0]:
  emit("pair_log_likelihood", $alpha,
    pair_total_log_likelihood(sites, alpha, 0.002))
