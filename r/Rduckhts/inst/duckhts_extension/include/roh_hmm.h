/* roh_hmm.h -- runs-of-homozygosity hidden Markov model, ported from bcftools
   vcfroh.c and HMM.c (two states, autozygous and Hardy-Weinberg; Viterbi path
   and forward-backward quality).

   Copyright (c) 2014-2026 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

   The engine has no DuckDB dependency. Coordinates: sites are pushed as
   zero-based positions (pos0); segments are reported one-based inclusive
   (pos1), as the bcftools RG lines are.
 */
#ifndef DUCKHTS_ROH_HMM_H
#define DUCKHTS_ROH_HMM_H

#include <stddef.h>
#include <stdint.h>

#define DUCKHTS_ROH_PHRED_TABLE 256

typedef struct {
    int32_t start_pos1;
    int32_t end_pos1;
    int32_t n_markers;
    double quality; /* mean forward-backward phred score of the run's sites */
} duckhts_roh_segment_t;

/* Model parameters for one run. Borrowed pointers must stay valid until the
 * last duckhts_roh_decode() of the run. */
typedef struct {
    double hw_to_az;           /* P(AZ|HW) per base pair, bcftools --hw-to-az */
    double az_to_hw;           /* P(HW|AZ) per base pair, bcftools --az-to-hw */
    double rec_rate;           /* > 0: constant rate (-M); scales a map when present */
    const int32_t *map_pos0;   /* genetic map nodes, zero-based, ascending; or NULL */
    const double *map_rate;    /* map nodes in Morgans (cM * 0.01) */
    size_t map_n;
} duckhts_roh_params_t;

typedef struct duckhts_roh duckhts_roh_t;

duckhts_roh_t *duckhts_roh_create(void);
void duckhts_roh_destroy(duckhts_roh_t *roh);

/* Starts a run for a list of at most max_sites sites. Returns 0, or -1 when
 * storage cannot be allocated. */
int duckhts_roh_begin(duckhts_roh_t *roh, const duckhts_roh_params_t *params, size_t max_sites);

/* P(D|genotype) for the three genotypes from phred-scaled likelihoods
 * (values above 255 count as 255). Returns 0 when the site is unusable: the
 * three values are equal or their sum is zero (bcftools skips such sites). */
int duckhts_roh_pdg_from_pl(const duckhts_roh_t *roh, int32_t pl_rr, int32_t pl_ra,
                            int32_t pl_aa, double pdg[3]);

/* Whether three phred-scaled likelihoods are usable evidence: none is negative
 * and they are not all equal. duckhts_roh_pdg_from_pl() applies the same test. */
int duckhts_roh_pl_usable(int32_t pl_rr, int32_t pl_ra, int32_t pl_aa);

/* P(D|genotype) for likelihoods already found usable and capped at 255, as a
 * packed site list stores them. The values are those duckhts_roh_pdg_from_pl()
 * gives for the uncapped likelihoods. */
void duckhts_roh_pdg_from_capped_pl(const duckhts_roh_t *roh, const uint8_t pl[3], double pdg[3]);

/* P(D|genotype) for a called dosage 0, 1 or 2 with error phred gt_error
 * (bcftools -G). Returns 0 for an unusable dosage. */
int duckhts_roh_pdg_from_gt(double gt_error, int dosage, double pdg[3]);

/* The read model of one sample: the four numbers duckhts_count_error_fit()
 * reports under these names. */
typedef struct {
    double seq_error;       /* a read shows the other allele at a homozygous site */
    double contamination;   /* share of reads from a second genome of af's population */
    double allele_balance;  /* share of a heterozygote's reads that show the counted allele */
    double artefact_weight; /* chance that a site's counts are uniform on 0..depth */
} duckhts_roh_count_model_t;

/* P(D|genotype) from allele read counts at a biallelic site. This is a DuckHTS
 * extension, not part of bcftools roh. Reads are independent; each shows the
 * counted allele (the one af refers to) with probability
 *   (1 - contamination) * q_g + contamination * c,
 * where q_g is seq_error, allele_balance and 1 - seq_error for the genotypes
 * with zero, one and two copies, and c is the mean of q_h over a contaminating
 * genotype h in Hardy-Weinberg proportions at af. With chance artefact_weight
 * the counted reads are instead uniform on 0..depth, whatever the genotype:
 *   P(D|g) = (1 - artefact_weight) * Binomial(counted; depth, p_g)
 *            + artefact_weight / (depth + 1).
 * With artefact_weight 0 the binomial coefficient is common to the three
 * genotypes and omitted. Relative likelihoods are floored at 10^(-25.5), the
 * PL path's cap of 255, so a site whose frequency is 1 cannot zero both
 * emissions. The caller guarantees counts >= 0, 0 < seq_error < 1/2,
 * 0 <= contamination < 1, 0 < allele_balance < 1, 0 <= artefact_weight < 1
 * and 0 < af <= 1. Returns 0 when the site has no reads. */
int duckhts_roh_pdg_from_counts(int32_t other_count, int32_t counted_count,
                                const duckhts_roh_count_model_t *model, double af, double pdg[3]);

/* Appends a usable site: the pdg values are normalised here, as bcftools
 * does, then the emission probabilities for the allele frequency af follow. */
void duckhts_roh_push(duckhts_roh_t *roh, int32_t pos0, const double pdg[3], double af);

size_t duckhts_roh_site_count(const duckhts_roh_t *roh);

/* Viterbi and forward-backward over the pushed sites. Returns 0, or -1 when
 * storage cannot be allocated. */
int duckhts_roh_decode(duckhts_roh_t *roh);

/* Segments of the last decode, valid until the next begin(). */
const duckhts_roh_segment_t *duckhts_roh_segments(const duckhts_roh_t *roh, size_t *count);

#endif
