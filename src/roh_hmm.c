/* roh_hmm.c -- runs-of-homozygosity HMM, ported from bcftools vcfroh.c and HMM.c.

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

   DuckHTS changes: the two-state model is fixed, buffers are owned by one
   workspace object, and the output is an array of segments rather than RG
   text lines. The arithmetic, including the order of the floating-point
   operations, follows the pinned bcftools source so that Viterbi paths and
   qualities reproduce `bcftools roh`. duckhts_roh_pdg_from_counts() is a
   DuckHTS addition with no bcftools counterpart.
 */
#include "roh_hmm.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define STATE_HW 0 /* Hardy-Weinberg */
#define STATE_AZ 1 /* autozygous */
#define NSTATES 2
#define NTPROB 10000 /* precomputed transition matrices, bcftools hmm_init(..., 10000) */
#define MAT(matrix, i, j) (matrix)[NSTATES * (i) + (j)] /* P(i|j): transition j -> i */

struct duckhts_roh {
    double pl2p[DUCKHTS_ROH_PHRED_TABLE];

    duckhts_roh_params_t params;
    int igenmap; /* cursor into the genetic map, carried across the three passes */

    double tprob_arr[NTPROB * NSTATES * NSTATES]; /* tprob_arr[i] = T^(i+1) */
    double tprob_t2az, tprob_t2hw;
    int tprob_ready;

    size_t nsites, msites;
    uint32_t *sites; /* zero-based positions */
    double *eprob;   /* [NSTATES * nsites] */
    uint8_t *vpath;  /* [NSTATES * nsites] */
    double *fwd;     /* [NSTATES * (nsites + 1)] */

    size_t nseg, mseg;
    duckhts_roh_segment_t *segments;
};

duckhts_roh_t *duckhts_roh_create(void) {
    duckhts_roh_t *roh = calloc(1, sizeof(*roh));
    if (roh == NULL) return NULL;
    for (int i = 0; i < DUCKHTS_ROH_PHRED_TABLE; i++) roh->pl2p[i] = pow(10., -i / 10.);
    return roh;
}

void duckhts_roh_destroy(duckhts_roh_t *roh) {
    if (roh == NULL) return;
    free(roh->sites);
    free(roh->eprob);
    free(roh->vpath);
    free(roh->fwd);
    free(roh->segments);
    free(roh);
}

static void multiply_matrix(const double *a, const double *b, double *dst) {
    double out[NSTATES * NSTATES];
    for (int i = 0; i < NSTATES; i++) {
        for (int j = 0; j < NSTATES; j++) {
            double val = 0;
            for (int k = 0; k < NSTATES; k++) val += MAT(a, i, k) * MAT(b, k, j);
            MAT(out, i, j) = val;
        }
    }
    memcpy(dst, out, sizeof(out));
}

static void init_tprob(duckhts_roh_t *roh) {
    double t2az = roh->params.hw_to_az;
    double t2hw = roh->params.az_to_hw;
    if (roh->tprob_ready && roh->tprob_t2az == t2az && roh->tprob_t2hw == t2hw) return;

    double *first = roh->tprob_arr;
    MAT(first, STATE_HW, STATE_HW) = 1 - t2az;
    MAT(first, STATE_HW, STATE_AZ) = t2hw;
    MAT(first, STATE_AZ, STATE_HW) = t2az;
    MAT(first, STATE_AZ, STATE_AZ) = 1 - t2hw;
    for (int i = 1; i < NTPROB; i++) {
        multiply_matrix(first, roh->tprob_arr + (size_t)(i - 1) * NSTATES * NSTATES,
                        roh->tprob_arr + (size_t)i * NSTATES * NSTATES);
    }
    roh->tprob_t2az = t2az;
    roh->tprob_t2hw = t2hw;
    roh->tprob_ready = 1;
}

/* Genetic distance in Morgans between two zero-based positions, interpolated
 * from the map nodes exactly as bcftools get_genmap_rate() does, including
 * its use of the slope between the bracketing nodes and its cursor. */
static double get_genmap_rate(duckhts_roh_t *roh, int start, int end) {
    const int32_t *pos = roh->params.map_pos0;
    const double *rate = roh->params.map_rate;
    int ngenmap = (int)roh->params.map_n;

    int i = roh->igenmap;
    if (pos[i] > start) {
        while (i > 0 && pos[i] > start) i--;
    } else {
        while (i + 1 < ngenmap && pos[i + 1] < start) i++;
    }
    int j = i;
    while (j + 1 < ngenmap && pos[j] < end) j++;
    if (i == j) {
        roh->igenmap = i;
        return 0;
    }

    if (start < pos[i]) start = pos[i];
    if (end > pos[j]) end = pos[j];
    double value = (rate[j] - rate[i]) / (pos[j] - pos[i]) * (end - start);
    roh->igenmap = j;
    return value;
}

static void rescale_tprob(double *tprob, double ci) {
    if (ci > 1) ci = 1;
    MAT(tprob, STATE_HW, STATE_AZ) *= ci;
    MAT(tprob, STATE_AZ, STATE_HW) *= ci;
    MAT(tprob, STATE_AZ, STATE_AZ) = 1 - MAT(tprob, STATE_HW, STATE_AZ);
    MAT(tprob, STATE_HW, STATE_HW) = 1 - MAT(tprob, STATE_AZ, STATE_HW);
}

/* Transition matrix between two sites; pos_diff and the hook follow
 * hmm_run_viterbi(): T^(distance) from the precomputed powers, then the map
 * or constant-rate hook rescales the crossover terms. */
static void set_tprob(duckhts_roh_t *roh, uint32_t prev_pos, uint32_t pos, double *curr) {
    int pos_diff = pos == prev_pos ? 0 : (int)(pos - prev_pos - 1);
    int n = pos_diff % NTPROB;
    memcpy(curr, roh->tprob_arr + (size_t)n * NSTATES * NSTATES, sizeof(double) * NSTATES * NSTATES);
    int nblocks = pos_diff / NTPROB;
    const double *block = roh->tprob_arr + (size_t)(NTPROB - 1) * NSTATES * NSTATES;
    for (int i = 0; i < nblocks; i++) multiply_matrix(block, curr, curr);

    if (roh->params.map_pos0 != NULL) {
        double ci = get_genmap_rate(roh, (int)prev_pos, (int)pos);
        if (roh->params.rec_rate != 0) ci *= roh->params.rec_rate;
        rescale_tprob(curr, ci);
    } else if (roh->params.rec_rate > 0) {
        rescale_tprob(curr, (pos - prev_pos) * roh->params.rec_rate);
    }
}

int duckhts_roh_begin(duckhts_roh_t *roh, const duckhts_roh_params_t *params, size_t max_sites) {
    roh->params = *params;
    roh->igenmap = 0;
    roh->nsites = 0;
    roh->nseg = 0;
    init_tprob(roh);

    if (max_sites > roh->msites) {
        if (max_sites > SIZE_MAX / (sizeof(double) * NSTATES) - NSTATES) return -1;
        uint32_t *sites = realloc(roh->sites, max_sites * sizeof(*sites));
        if (sites == NULL) return -1;
        roh->sites = sites;
        double *eprob = realloc(roh->eprob, max_sites * NSTATES * sizeof(*eprob));
        if (eprob == NULL) return -1;
        roh->eprob = eprob;
        uint8_t *vpath = realloc(roh->vpath, max_sites * NSTATES * sizeof(*vpath));
        if (vpath == NULL) return -1;
        roh->vpath = vpath;
        double *fwd = realloc(roh->fwd, (max_sites + 1) * NSTATES * sizeof(*fwd));
        if (fwd == NULL) return -1;
        roh->fwd = fwd;
        roh->msites = max_sites;
    }
    return 0;
}

int duckhts_roh_pdg_from_pl(const duckhts_roh_t *roh, int32_t pl_rr, int32_t pl_ra,
                            int32_t pl_aa, double pdg[3]) {
    if (pl_rr < 0 || pl_ra < 0 || pl_aa < 0) return 0;
    if (pl_rr == pl_ra && pl_rr == pl_aa) return 0;
    int max = DUCKHTS_ROH_PHRED_TABLE - 1;
    pdg[0] = roh->pl2p[pl_rr < max ? pl_rr : max];
    pdg[1] = roh->pl2p[pl_ra < max ? pl_ra : max];
    pdg[2] = roh->pl2p[pl_aa < max ? pl_aa : max];
    return pdg[0] + pdg[1] + pdg[2] != 0;
}

int duckhts_roh_pdg_from_gt(double gt_error, int dosage, double pdg[3]) {
    double unseen_pl = pow(10, -gt_error / 10.);
    if (dosage == 1) {
        pdg[0] = pdg[2] = unseen_pl;
        pdg[1] = 1 - 2 * unseen_pl;
    } else if (dosage == 0) {
        pdg[0] = 1 - unseen_pl - unseen_pl * unseen_pl;
        pdg[1] = unseen_pl;
        pdg[2] = unseen_pl * unseen_pl;
    } else if (dosage == 2) {
        pdg[0] = unseen_pl * unseen_pl;
        pdg[1] = unseen_pl;
        pdg[2] = 1 - unseen_pl - unseen_pl * unseen_pl;
    } else {
        return 0;
    }
    return pdg[0] + pdg[1] + pdg[2] != 0;
}

/* DuckHTS extension: read-count emissions with sequencing error and
 * contamination (see roh_hmm.h). Log likelihoods relative to the largest keep
 * deep sites from underflowing; push() normalises the three values. */
int duckhts_roh_pdg_from_counts(int32_t other_count, int32_t counted_count, double seq_error,
                                double contamination, double af, double pdg[3]) {
    if ((int64_t)other_count + counted_count == 0) return 0;
    double contaminant = af * (1 - seq_error) + (1 - af) * seq_error;
    double q[3] = {seq_error, 0.5, 1 - seq_error};
    double log_likelihood[3];
    double best = -INFINITY;
    for (int g = 0; g < 3; g++) {
        double p = (1 - contamination) * q[g] + contamination * contaminant;
        log_likelihood[g] = counted_count * log(p) + other_count * log1p(-p);
        if (log_likelihood[g] > best) best = log_likelihood[g];
    }
    for (int g = 0; g < 3; g++) pdg[g] = exp(log_likelihood[g] - best);
    return 1;
}

void duckhts_roh_push(duckhts_roh_t *roh, int32_t pos0, const double pdg_in[3], double af) {
    double pdg[3];
    double sum = pdg_in[0] + pdg_in[1] + pdg_in[2];
    for (int j = 0; j < 3; j++) pdg[j] = pdg_in[j] / sum;

    double *eprob = &roh->eprob[NSTATES * roh->nsites];
    eprob[STATE_AZ] = pdg[0] * (1 - af) + pdg[2] * af;
    eprob[STATE_HW] = pdg[0] * (1 - af) * (1 - af) + 2 * pdg[1] * (1 - af) * af + pdg[2] * af * af;
    roh->sites[roh->nsites] = (uint32_t)pos0;
    roh->nsites++;
}

size_t duckhts_roh_site_count(const duckhts_roh_t *roh) {
    return roh->nsites;
}

static double phred_score(double prob) {
    if (prob == 0) return 99;
    prob = -4.3429 * log(prob);
    return prob > 99 ? 99 : prob;
}

static void run_viterbi(duckhts_roh_t *roh) {
    size_t n = roh->nsites;
    double vprob_a[NSTATES] = {1. / NSTATES, 1. / NSTATES};
    double vprob_b[NSTATES];
    double *vprob = vprob_a, *vprob_tmp = vprob_b;
    double curr[NSTATES * NSTATES];
    uint32_t prev_pos = roh->sites[0];

    for (size_t i = 0; i < n; i++) {
        uint8_t *vpath = &roh->vpath[i * NSTATES];
        const double *eprob = &roh->eprob[i * NSTATES];

        set_tprob(roh, prev_pos, roh->sites[i], curr);
        prev_pos = roh->sites[i];

        double vnorm = 0;
        for (int j = 0; j < NSTATES; j++) {
            double vmax = 0;
            int k_vmax = 0;
            for (int k = 0; k < NSTATES; k++) {
                double pval = vprob[k] * MAT(curr, j, k);
                if (vmax < pval) {
                    vmax = pval;
                    k_vmax = k;
                }
            }
            vpath[j] = (uint8_t)k_vmax;
            vprob_tmp[j] = vmax * eprob[j];
            vnorm += vprob_tmp[j];
        }
        for (int j = 0; j < NSTATES; j++) vprob_tmp[j] /= vnorm;
        double *swap = vprob;
        vprob = vprob_tmp;
        vprob_tmp = swap;
    }

    int iptr = 0;
    for (int i = 1; i < NSTATES; i++) {
        if (vprob[iptr] < vprob[i]) iptr = i;
    }

    /* Trace back, reusing vpath[i * NSTATES] for the chosen state. */
    for (size_t i = n; i-- > 0;) {
        int iptr_prev = roh->vpath[i * NSTATES + (size_t)iptr];
        roh->vpath[i * NSTATES] = (uint8_t)iptr;
        iptr = iptr_prev;
    }
}

/* Posterior state probabilities land in fwd[NSTATES * (i + 1) + state]. */
static void run_fwd_bwd(duckhts_roh_t *roh) {
    size_t n = roh->nsites;
    double curr[NSTATES * NSTATES];
    double bwd_a[NSTATES] = {1, 1};
    double bwd_b[NSTATES];
    double *bwd = bwd_a, *bwd_tmp = bwd_b;

    roh->fwd[0] = roh->fwd[1] = 1. / NSTATES;
    uint32_t prev_pos = roh->sites[0];
    for (size_t i = 0; i < n; i++) {
        double *fwd_prev = &roh->fwd[i * NSTATES];
        double *fwd = &roh->fwd[(i + 1) * NSTATES];
        const double *eprob = &roh->eprob[i * NSTATES];

        set_tprob(roh, prev_pos, roh->sites[i], curr);
        prev_pos = roh->sites[i];

        double norm = 0;
        for (int j = 0; j < NSTATES; j++) {
            double pval = 0;
            for (int k = 0; k < NSTATES; k++) pval += fwd_prev[k] * MAT(curr, j, k);
            fwd[j] = pval * eprob[j];
            norm += fwd[j];
        }
        for (int j = 0; j < NSTATES; j++) fwd[j] /= norm;
    }

    prev_pos = roh->sites[n - 1];
    for (size_t i = 0; i < n; i++) {
        double *fwd = &roh->fwd[(n - i) * NSTATES];
        const double *eprob = &roh->eprob[(n - i - 1) * NSTATES];

        set_tprob(roh, roh->sites[n - i - 1], prev_pos, curr);
        prev_pos = roh->sites[n - i - 1];

        double bwd_norm = 0;
        for (int j = 0; j < NSTATES; j++) {
            double pval = 0;
            for (int k = 0; k < NSTATES; k++) pval += bwd[k] * eprob[k] * MAT(curr, k, j);
            bwd_tmp[j] = pval;
            bwd_norm += pval;
        }
        double norm = 0;
        for (int j = 0; j < NSTATES; j++) {
            bwd_tmp[j] /= bwd_norm;
            fwd[j] *= bwd[j];
            norm += fwd[j];
        }
        for (int j = 0; j < NSTATES; j++) fwd[j] /= norm;
        double *swap = bwd_tmp;
        bwd_tmp = bwd;
        bwd = swap;
    }
}

static int append_segment(duckhts_roh_t *roh, uint32_t beg, uint32_t end, uint32_t nqual,
                          double qual) {
    if (roh->nseg == roh->mseg) {
        size_t capacity = roh->mseg ? roh->mseg * 2 : 16;
        if (capacity > SIZE_MAX / sizeof(*roh->segments)) return -1;
        duckhts_roh_segment_t *grown = realloc(roh->segments, capacity * sizeof(*grown));
        if (grown == NULL) return -1;
        roh->segments = grown;
        roh->mseg = capacity;
    }
    duckhts_roh_segment_t *segment = &roh->segments[roh->nseg++];
    segment->start_pos1 = (int32_t)(beg + 1);
    segment->end_pos1 = (int32_t)(end + 1);
    segment->n_markers = (int32_t)nqual;
    segment->quality = qual / nqual;
    return 0;
}

int duckhts_roh_decode(duckhts_roh_t *roh) {
    roh->nseg = 0;
    if (roh->nsites == 0) return 0;

    roh->igenmap = 0;
    run_viterbi(roh);
    run_fwd_bwd(roh);

    int rg_state = 0;
    uint32_t beg = 0, end = 0, nqual = 0;
    double qual_sum = 0;
    for (size_t i = 0; i < roh->nsites; i++) {
        int state = roh->vpath[i * NSTATES] == STATE_AZ ? 1 : 0;
        double qual = phred_score(1.0 - roh->fwd[NSTATES * (i + 1) + (size_t)state]);
        if (state != rg_state) {
            if (!state) {
                if (append_segment(roh, beg, end, nqual, qual_sum) != 0) return -1;
                rg_state = 0;
            } else {
                rg_state = 1;
                beg = end = roh->sites[i];
                qual_sum = qual;
                nqual = 1;
            }
        } else if (state) {
            nqual++;
            qual_sum += qual;
            end = roh->sites[i];
        }
    }
    if (rg_state && append_segment(roh, beg, end, nqual, qual_sum) != 0) return -1;
    return 0;
}

const duckhts_roh_segment_t *duckhts_roh_segments(const duckhts_roh_t *roh, size_t *count) {
    *count = roh->nseg;
    return roh->segments;
}
