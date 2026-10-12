# Version 2 of the truth of the planted-truth evaluation. Version 1
# (declaration.R) is kept as it is and its verdicts are still reported.
#
# Why. Version 1 took the background of a child of two unrelated founders to
# hold no runs, so that every called base outside a planted segment was
# false. The full run, with the memory gate recording only (not a render),
# showed that premise false. The long runs that `counts_fit` called outside
# the planted segments held 0.0038 true heterozygous sites per panel site,
# against about 0.30 for the children's background as a whole, and 142 of
# their 270 Mb lay in stretches of at least 500 kb without a single true
# heterozygous site. `gt30` and `counts_truth` showed the same. Those bases are
# homozygous in the constructed child; a decoder that calls them is right.
#
# Truth, version 2: the constructed child's true genotypes, the genotypes of
# the drawn founder haplotypes after planting and before the planted genotype
# errors and the read simulation. A true base of a child is a base of
#   - a planted segment (start to end, as in version 1), or
#   - a background stretch: a maximal set of consecutive panel sites at which
#     the child's true genotype is not heterozygous, from the position of its
#     first site to the position of its last site, one-based and inclusive,
#     whose length (last - first + 1) is at least min_stretch_bases.
# The two sets are merged into disjoint intervals. Consecutive means
# consecutive in the panel's position order; a stretch ends at the site
# before a true heterozygous site and starts at the site after one, or at the
# first or last panel site of the chromosome.
#
# Everything else is as in version 1: the same children, arms, run classes,
# metrics (sensitivity: true called bases over true bases; false discovery:
# called bases outside the truth over called bases), bootstrap, comparisons
# and criteria T1 (the same 0.01 bound on the false-discovery increase) and
# T2 (the same margin), now against this truth. T3 does not depend on the
# truth and is reported once. The stability measures are reported against
# both truths. Both versions' verdicts are reported side by side.
roh_truth_declaration_v2 <- list(
  min_stretch_bases = 5e5)
