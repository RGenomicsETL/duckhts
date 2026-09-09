#include "duckvep_property.h"

/* ===================================================================== *
 * annotate_tile FUSION: the composed engine == the sum of its tested parts.
 *
 * annotate_tile fuses the sweep (candidate pairs) with the classifier (region
 * mask) and the SO mapping. The strongest hermetic check: for every generated
 * scene, the rows annotate_tile emits must equal the INDEPENDENT composition —
 * the same sweep set, each pair mapped to its structural SO terms by a from-spec
 * map written here (not the kernel's mapper). kprop_scene transcripts carry no
 * exons and no CDS, so they exercise the upstream/downstream/intron branches at
 * scale + the sweep<->builder<->filter wiring; the deterministic scene below
 * pins the exon/UTR/CDS/splice branches exactly.
 * ===================================================================== */

static enum theft_trial_res prop_annotate_matches_composition(struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    duckvep_exon_model_t exons;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    static duckvep_consequence_t rows[KPROP_MAX_PAIRS];
    static uint64_t sweep_buf[KPROP_MAX_PAIRS];
    duckvep_result_builder_t rb;
    uint32_t active[KPROP_MAX_TX];
    uint32_t candidates[KPROP_MAX_TX];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    size_t n_sweep, i;
    const uint32_t HALO = KPROP_SWEEP_HALO; /* fixed >0 so options.halo is unambiguous */
    (void)t;

    memset(&exons, 0, sizeof exons); /* no exons -> non_coding_transcript (placement) / up / down */
    memset(&err, 0, sizeof err);

    /* Scene transcripts are sorted by (chrom_id, start1) by construction, so the
     * model must open. */
    if (duckvep_model_open(&s->tx, &exons, NULL, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    memset(&init, 0, sizeof init);
    /* up == down == halo: every sweep-admitted pair is within the directional
     * reach, so the filter is inert, and options.halo (now clamped to
     * max(halo,up,down)) stays exactly HALO. */
    init.upstream_dist = HALO;
    init.downstream_dist = HALO;
    init.halo = HALO;
    if (duckvep_options_open(&init, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }

    /* Independent composition: same sweep, then a from-spec structural map. */
    col.buf = sweep_buf; col.n = 0u; col.cap = KPROP_MAX_PAIRS;
    n_sweep = duckvep_sweep_candidates(&s->v, &s->tx, HALO,
                                       active, KPROP_MAX_TX,
                                       candidates, KPROP_MAX_TX,
                                       pair_sink, &col, &st);
    if (st != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    /* Every in-window pair maps to a non-empty mask (non_coding_transcript/up/down) and the
     * filter is inert, so annotate emits exactly the sweep set in sweep order. */
    if (duckvep_result_builder_count(&rb) != n_sweep) { res = THEFT_TRIAL_FAIL; goto done; }

    for (i = 0u; i < n_sweep; i++) {
        uint32_t vi = (uint32_t)(col.buf[i] >> 32);
        uint32_t ti = (uint32_t)(col.buf[i] & 0xffffffffu);
        uint32_t start = s->v.pos1[vi];
        uint32_t end = s->v.variant_kind[vi] == (uint8_t)DUCKVEP_KIND_SV
                         ? s->v.end1[vi]
                         : start;
        uint32_t ts = s->tx.start1[ti];
        uint32_t te = s->tx.end1[ti];
        int fwd = s->tx.strand[ti] >= 0;
        uint32_t exp_region;
        uint64_t exp_mask;
        const duckvep_consequence_t *row = &rows[i];

        if (end < ts)        exp_region = fwd ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                                              : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        else if (start > te) exp_region = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                                              : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        else                 exp_region = (uint32_t)DUCKVEP_REGION_INTRON;

        if (exp_region == (uint32_t)DUCKVEP_REGION_UPSTREAM)
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE);
        else if (exp_region == (uint32_t)DUCKVEP_REGION_DOWNSTREAM)
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE);
        else
            /* zero-exon non-coding transcript: no exons -> no introns -> within_intron is
             * FALSE, so NO intron_variant; the in-span position is non_coding_transcript_variant
             * (placement). This is VEP-correct (the old INTRON co-emission was an over-call). */
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT);

        if (row->variant_idx != vi || row->tx_idx != ti) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->region_mask != exp_region) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->consequence_mask != exp_mask) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->impact != (uint8_t)DUCKVEP_IMPACT_MODIFIER) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->cdna_pos != -1 || row->cds_pos != -1 || row->protein_pos != -1) {
            res = THEFT_TRIAL_FAIL; goto done;
        }
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_matches_composition_for_any_scene(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile == sweep + classify + structural-SO composition";
    cfg.prop1 = prop_annotate_matches_composition;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

int consequence_rows_equal(const duckvep_consequence_t *a,
                                  const duckvep_consequence_t *b) {
    return a->variant_idx == b->variant_idx &&
           a->tx_idx == b->tx_idx &&
           a->gene_idx == b->gene_idx &&
           a->consequence_mask == b->consequence_mask &&
           a->region_mask == b->region_mask &&
           a->flags == b->flags &&
           a->impact == b->impact &&
           a->sequence_status == b->sequence_status &&
           a->nmd_prediction == b->nmd_prediction &&
           a->nmd_escape_reasons == b->nmd_escape_reasons &&
           a->cdna_pos == b->cdna_pos &&
           a->cds_pos == b->cds_pos &&
           a->protein_pos == b->protein_pos &&
           a->aa_ref == b->aa_ref &&
           a->aa_alt == b->aa_alt &&
           a->overlap_object_kind == b->overlap_object_kind;
}

static int collect_annotation_cursor(
    const duckvep_model_t          *model,
    const duckvep_variant_batch_t  *variants,
    const duckvep_options_t        *options,
    duckvep_workspace_t            *workspace,
    size_t                          chunk_capacity,
    duckvep_consequence_t          *rows,
    size_t                          row_capacity,
    size_t                         *row_count) {

    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_consequence_t chunk[9];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    size_t count = 0u;
    int ok = 1;

    memset(&error, 0, sizeof error);
    if (row_count != NULL) *row_count = 0u;
    if (chunk_capacity == 0u || chunk_capacity > 9u || row_count == NULL ||
        (row_capacity != 0u && rows == NULL)) {
        return 0;
    }
    if (duckvep_annotate_cursor_open(
            model, variants, options, workspace, &cursor,
            &error) != DUCKVEP_OK) {
        return 0;
    }
    while (!duckvep_annotate_cursor_done(cursor)) {
        duckvep_status_t status;
        size_t i;

        duckvep_result_builder_init(&builder, chunk, chunk_capacity);
        status = duckvep_annotate_cursor_fill(cursor, &builder, &error);
        if (status != DUCKVEP_OK && status != DUCKVEP_ERR_RESULT_FULL) {
            ok = 0;
            break;
        }
        if (status == DUCKVEP_ERR_RESULT_FULL &&
            duckvep_result_builder_count(&builder) == 0u) {
            ok = 0;
            break;
        }
        for (i = 0u; i < duckvep_result_builder_count(&builder); i++) {
            if (count >= row_capacity) {
                ok = 0;
                break;
            }
            rows[count++] = chunk[i];
        }
        if (!ok) break;
    }
    duckvep_annotate_cursor_close(cursor);
    if (ok) *row_count = count;
    return ok;
}

static struct {
    uint64_t far_directional_snv;
    uint64_t simple_point_snv;
    uint64_t generalized_pair;
    uint64_t nmd_transcript_rows;
    uint64_t mirna_transcripts;
    uint64_t coding_transcripts;
    uint64_t cursor_splits;
} g_annotation_shortcut_cov;

/* Sorted-stream shortcuts are optimizations, not a second consequence
 * authority.  This property keeps the candidate sweep identical and compares
 * every public consequence field and row order against an internal mode that
 * disables both direct SNV emitters and both sorted exon cursors.  The same
 * comparison is then repeated through arbitrarily small resumable output
 * buffers.  Transcript modes cover non-coding, pre-labelled NMD, miRNA, and
 * topology-only coding models on both strands; SV records in the generated
 * scene remain useful controls and always take the generalized path. */
static enum theft_trial_res prop_annotation_shortcuts_match_generalized(
    struct theft *t,
    void         *arg1) {

    const struct kprop_scene *scene = (const struct kprop_scene *)arg1;
    duckvep_transcript_model_t transcripts = scene->tx;
    duckvep_exon_model_t exons;
    uint64_t flags[KPROP_MAX_TX];
    uint32_t exon_offset[KPROP_MAX_TX];
    uint16_t exon_count[KPROP_MAX_TX];
    uint32_t exon_start1[KPROP_MAX_TX];
    uint32_t exon_end1[KPROP_MAX_TX];
    uint32_t cds_start1[KPROP_MAX_TX];
    uint32_t cds_end1[KPROP_MAX_TX];
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *fast_workspace = NULL;
    duckvep_workspace_t *generalized_workspace = NULL;
    duckvep_options_init_t options_init;
    duckvep_error_t error;
    duckvep_consequence_t fast_rows[KPROP_MAX_PAIRS];
    duckvep_consequence_t generalized_rows[KPROP_MAX_PAIRS];
    duckvep_consequence_t fast_cursor_rows[KPROP_MAX_PAIRS];
    duckvep_consequence_t generalized_cursor_rows[KPROP_MAX_PAIRS];
    duckvep_result_builder_t builder;
    const duckvep_workspace_delta_route_stats_t *fast_stats;
    const duckvep_workspace_delta_route_stats_t *generalized_stats;
    size_t fast_count;
    size_t generalized_count;
    size_t fast_cursor_count;
    size_t generalized_cursor_count;
    size_t chunk_capacity;
    size_t i;
    enum theft_trial_res result = THEFT_TRIAL_PASS;
    (void)t;

    memset(&exons, 0, sizeof exons);
    memset(&options_init, 0, sizeof options_init);
    memset(&error, 0, sizeof error);
    for (i = 0u; i < transcripts.transcript_count; i++) {
        uint32_t mode = (uint32_t)((scene->tstart[i] + i) % 5u);

        flags[i] = 0u;
        exon_offset[i] = (uint32_t)i;
        exon_count[i] = 1u;
        exon_start1[i] = scene->tstart[i];
        exon_end1[i] = scene->tend[i];
        cds_start1[i] = 0u;
        cds_end1[i] = 0u;
        if (mode == 1u) {
            flags[i] = (uint64_t)DUCKVEP_TX_BIOTYPE_NMD;
        } else if (mode == 2u) {
            flags[i] = (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA;
            g_annotation_shortcut_cov.mirna_transcripts++;
        } else if (mode >= 3u) {
            flags[i] = (uint64_t)DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
            if (mode == 4u) flags[i] |= (uint64_t)DUCKVEP_TX_BIOTYPE_NMD;
            cds_start1[i] = scene->tstart[i];
            cds_end1[i] = scene->tend[i];
            g_annotation_shortcut_cov.coding_transcripts++;
        }
    }
    transcripts.flags = flags;
    transcripts.exon_offset = exon_offset;
    transcripts.exon_count = exon_count;
    transcripts.cds_start1 = cds_start1;
    transcripts.cds_end1 = cds_end1;
    exons.start1 = exon_start1;
    exons.end1 = exon_end1;
    exons.exon_count = transcripts.transcript_count;

    /* The scene generator deliberately samples zero, both sides of VEP's
     * 5,000-base default, 10 kb, and 65,535 bases.  Use that distance in the
     * optimized-vs-generalized comparison so a shortcut cannot accidentally
     * bake in the historical SQL default while still passing this oracle. */
    options_init.upstream_dist = scene->halo;
    options_init.downstream_dist = scene->halo;
    options_init.halo = scene->halo;
    options_init.distances_are_explicit = 1u;
    if (duckvep_model_open(
            &transcripts, &exons, NULL, NULL, &model, &error) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(&options_init, &options, &error) != DUCKVEP_OK ||
        duckvep_workspace_open(model, &fast_workspace, &error) != DUCKVEP_OK ||
        duckvep_workspace_open(
            model, &generalized_workspace, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    duckvep_workspace_force_generalized_annotation(
        generalized_workspace, 1);
    duckvep_workspace_delta_route_stats_reset(fast_workspace);
    duckvep_workspace_delta_route_stats_reset(generalized_workspace);

    duckvep_result_builder_init(&builder, fast_rows, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(
            model, &scene->v, options, fast_workspace,
            &builder, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    fast_count = duckvep_result_builder_count(&builder);
    duckvep_result_builder_init(
        &builder, generalized_rows, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(
            model, &scene->v, options, generalized_workspace,
            &builder, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    generalized_count = duckvep_result_builder_count(&builder);
    if (fast_count != generalized_count) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    for (i = 0u; i < fast_count; i++) {
        if (!consequence_rows_equal(&fast_rows[i], &generalized_rows[i])) {
            result = THEFT_TRIAL_FAIL;
            goto done;
        }
        if ((fast_rows[i].consequence_mask &
             DUCKVEP_SO(DUCKVEP_SO_NMD_TRANSCRIPT)) != 0u) {
            g_annotation_shortcut_cov.nmd_transcript_rows++;
        }
    }

    fast_stats = duckvep_workspace_delta_route_stats(fast_workspace);
    generalized_stats =
        duckvep_workspace_delta_route_stats(generalized_workspace);
    if (fast_stats == NULL || generalized_stats == NULL ||
        generalized_stats->far_directional_snv != 0u ||
        generalized_stats->simple_point_snv != 0u) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    g_annotation_shortcut_cov.far_directional_snv +=
        fast_stats->far_directional_snv;
    g_annotation_shortcut_cov.simple_point_snv +=
        fast_stats->simple_point_snv;
    g_annotation_shortcut_cov.generalized_pair +=
        generalized_stats->generalized_pair;

    chunk_capacity = 1u +
        (scene->v.count + scene->tx.transcript_count + scene->halo) % 9u;
    if (!collect_annotation_cursor(
            model, &scene->v, options, fast_workspace, chunk_capacity,
            fast_cursor_rows, KPROP_MAX_PAIRS, &fast_cursor_count) ||
        !collect_annotation_cursor(
            model, &scene->v, options, generalized_workspace, chunk_capacity,
            generalized_cursor_rows, KPROP_MAX_PAIRS,
            &generalized_cursor_count) ||
        fast_cursor_count != fast_count ||
        generalized_cursor_count != generalized_count) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    for (i = 0u; i < fast_count; i++) {
        if (!consequence_rows_equal(&fast_rows[i], &fast_cursor_rows[i]) ||
            !consequence_rows_equal(
                &generalized_rows[i], &generalized_cursor_rows[i])) {
            result = THEFT_TRIAL_FAIL;
            goto done;
        }
    }
    g_annotation_shortcut_cov.cursor_splits++;

done:
    duckvep_workspace_close(generalized_workspace);
    duckvep_workspace_close(fast_workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    return result;
}

TEST annotation_shortcuts_match_generalized_for_any_sorted_scene(void) {
    struct theft_run_config config;

    memset(&config, 0, sizeof config);
    memset(&g_annotation_shortcut_cov, 0,
           sizeof g_annotation_shortcut_cov);
    config.name =
        "optimized sorted annotation == forced generalized full rows";
    config.prop1 = prop_annotation_shortcuts_match_generalized;
    config.type_info[0] = &kprop_scene_info;
    config.trials =
        kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    config.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&config));
    ASSERT(g_annotation_shortcut_cov.far_directional_snv > 0u);
    ASSERT(g_annotation_shortcut_cov.simple_point_snv > 0u);
    ASSERT(g_annotation_shortcut_cov.generalized_pair > 0u);
    ASSERT(g_annotation_shortcut_cov.nmd_transcript_rows > 0u);
    ASSERT(g_annotation_shortcut_cov.mirna_transcripts > 0u);
    ASSERT(g_annotation_shortcut_cov.coding_transcripts > 0u);
    ASSERT(g_annotation_shortcut_cov.cursor_splits > 0u);
    fprintf(stderr,
            "[annotation-shortcut coverage] far=%llu simple=%llu "
            "generalized=%llu nmd_rows=%llu mirna_tx=%llu coding_tx=%llu "
            "cursor_splits=%llu\n",
            (unsigned long long)g_annotation_shortcut_cov.far_directional_snv,
            (unsigned long long)g_annotation_shortcut_cov.simple_point_snv,
            (unsigned long long)g_annotation_shortcut_cov.generalized_pair,
            (unsigned long long)g_annotation_shortcut_cov.nmd_transcript_rows,
            (unsigned long long)g_annotation_shortcut_cov.mirna_transcripts,
            (unsigned long long)g_annotation_shortcut_cov.coding_transcripts,
            (unsigned long long)g_annotation_shortcut_cov.cursor_splits);
    PASS();
}

TEST annotate_cursor_resumes_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2]   = {200u, 200u};
    static const int8_t   strand[2] = {1, 1};
    static const uint64_t flags[2]  = {0u, 0u};
    static const uint32_t zero32[2] = {0u, 0u};
    static const uint16_t zero16[2] = {0u, 0u};
    uint16_t vchrom[1] = {0u};
    uint32_t vpos[1] = {150u};
    uint32_t vend[1] = {150u};
    uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rowbuf[1];
    duckvep_consequence_t got[2];
    duckvep_result_builder_t rb;
    size_t got_n = 0u;

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex); memset(&err, 0, sizeof err);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind; v.count = 1u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(model, &v, opts, ws, &cur, &err));

    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        duckvep_result_builder_init(&rb, rowbuf, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        ASSERT(st == DUCKVEP_OK || st == DUCKVEP_ERR_RESULT_FULL);
        if (duckvep_result_builder_count(&rb) == 1u) {
            ASSERT(got_n < 2u);
            got[got_n++] = rowbuf[0];
        }
    }
    ASSERT_EQ(2u, got_n);
    ASSERT_EQ(0u, got[0].variant_idx); ASSERT_EQ(0u, got[0].tx_idx);
    ASSERT_EQ(0u, got[1].variant_idx); ASSERT_EQ(1u, got[1].tx_idx);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT), got[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT), got[1].consequence_mask);

    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}


int kprop_context_snapshot(const duckvep_coding_context_t *context,
                                struct kprop_context_snapshot *out) {
    memset(out, 0, sizeof *out);
    out->metadata = *context;
    out->metadata.ref_cds = out->metadata.alt_cds = NULL;
    out->metadata.ref_peptide = out->metadata.alt_peptide = NULL;
    out->metadata.single_edit_alt = out->metadata.local_alt_cds = NULL;
    out->metadata.local_ref_peptide = out->metadata.local_alt_peptide = NULL;
    out->metadata.ref_peptide_edit_position1 = NULL;
    out->metadata.ref_peptide_edit_alt = NULL;
    out->metadata.pre_cds_bases = out->metadata.post_cds_bases = NULL;
    out->window_present = duckvep_coding_context_peptide_window_open(context, &out->window);
    for (int alt = 0; alt <= 1; alt++) {
        size_t cds_len = alt ? context->alt_cds_len : context->ref_cds_len;
        size_t protein_len = alt ? context->alt_peptide_len : context->ref_peptide_len;
        size_t local_len = !out->window_present ? 0u
            : alt ? out->window.alt_length : out->window.ref_length;
        if (cds_len > sizeof out->cds[alt] || protein_len > sizeof out->protein[alt] ||
            local_len > sizeof out->local_peptide[alt]) return 0;
        for (size_t i = 0u; i < cds_len; i++)
            out->cds[alt][i] = (uint8_t)duckvep_coding_context_cds_base(context, alt, i);
        for (size_t i = 0u; i < protein_len; i++)
            out->protein[alt][i] = duckvep_coding_context_peptide_base(context, alt, i);
        for (size_t i = 0u; i < local_len; i++)
            out->local_peptide[alt][i] = duckvep_coding_context_peptide_window_base(
                context, &out->window, alt, i);
    }
    return 1;
}

struct annotation_observer_capture {
    const duckvep_variant_batch_t *expected_batch;
    duckvep_consequence_t rows[6];
    duckvep_event_t events[6];
    duckvep_sequence_delta_t deltas[6];
    uint8_t trace_present[6];
    uint8_t transcript_edit_present[6];
    uint8_t transcript_edit_status[6];
    uint8_t coding_context_status[6];
    uint8_t delta_present[6];
    uint8_t coding_context_valid[6];
    uint32_t transcript_edit_first_cdna1[6];
    uint32_t transcript_edit_ref_length[6];
    uint32_t transcript_edit_alt_length[6];
    struct kprop_context_snapshot coding[6];
    size_t count;
};

static int annotation_observer_capture_row(
    void                             *context,
    const duckvep_variant_batch_t    *variants,
    const duckvep_consequence_t      *row,
    const duckvep_pair_facts_t       *facts) {

    struct annotation_observer_capture *capture =
        (struct annotation_observer_capture *)context;
    size_t index;

    if (capture == NULL || variants == NULL ||
        capture->expected_batch == NULL ||
        variants->count != capture->expected_batch->count ||
        variants->pos1 != capture->expected_batch->pos1 || row == NULL ||
        capture->count >= sizeof capture->rows / sizeof capture->rows[0]) {
        return 0;
    }
    index = capture->count++;
    capture->rows[index] = *row;
    capture->trace_present[index] = facts != NULL ? 1u : 0u;
    if (facts != NULL) {
        if (facts->event == NULL) return 0;
        capture->events[index] = *facts->event;
        capture->transcript_edit_status[index] =
            facts->transcript_edit_status;
        capture->coding_context_status[index] =
            facts->coding_context_status;
        capture->delta_present[index] = facts->delta != NULL ? 1u : 0u;
        if (facts->delta) capture->deltas[index] = *facts->delta;
        capture->coding_context_valid[index] = facts->coding_context_valid;
        if (facts->coding_context_valid && (!facts->coding_context ||
            !kprop_context_snapshot(facts->coding_context, &capture->coding[index]))) return 0;
        if (facts->transcript_edit != NULL) {
            capture->transcript_edit_present[index] = 1u;
            capture->transcript_edit_first_cdna1[index] =
                facts->transcript_edit->first.cdna_anchor1;
            capture->transcript_edit_ref_length[index] =
                facts->transcript_edit->ref_length;
            capture->transcript_edit_alt_length[index] =
                facts->transcript_edit->alt_length;
        }
    } else {
        memset(&capture->events[index], 0, sizeof capture->events[index]);
    }
    return 1;
}

TEST observed_single_pair_matches_cursor_without_sweep_storage(void) {
    static const uint8_t unpadded[] = "ATGGCTGCTGCTGCTGCTGCTGCTGCTTAA";
    static const uint8_t flanks[] = "AAAAAAAAA";
    uint32_t random = UINT32_C(173);
    size_t compared = 0u;
    for (unsigned layout = 0u; layout < 2u; layout++)
    for (unsigned phase = 0u; phase <= 2u; phase++)
    for (int strand = -1; strand <= 1; strand += 2) {
        struct kprop_proj_scene s = {0};
        s.chrom = 0u; s.tstart = s.cds_s = 11u; s.tend = s.cds_e = 40u;
        s.strand = (int8_t)strand; s.excnt = 1u;
        s.es[0] = 11u; s.ee[0] = 40u; s.cs[0] = 1u; s.ce[0] = 30u;
        s.phase[0] = s.end_phase[0] = (int8_t)phase;
        if (layout) {
            s.tstart = 1u; s.tend = 84u; s.cds_e = 79u; s.excnt = 3u;
            s.es[0] = 1u; s.ee[0] = 4u; s.ce[0] = 4u;
            s.es[1] = 11u; s.ee[1] = 21u; s.cs[1] = 5u; s.ce[1] = 15u;
            s.es[2] = 61u; s.ee[2] = 84u; s.cs[2] = 16u; s.ce[2] = 39u;
            s.phase[0] = s.end_phase[0] = -1;
            s.phase[1] = (int8_t)phase;
            s.end_phase[1] = s.phase[2] = (int8_t)((phase + 11u) % 3u);
            s.end_phase[2] = (int8_t)phase;
            if (strand < 0) {
                s.tstart = 16u; s.tend = 99u; s.cds_s = 21u; s.cds_e = 89u;
                for (unsigned e = 0u; e < 3u; e++) {
                    uint32_t start = s.es[e];
                    s.es[e] = 100u - s.ee[e]; s.ee[e] = 100u - start;
                }
            }
        }
        s.flags = DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
        kprop_proj_scene_finish(&s);
        uint8_t cds[sizeof unpadded + 2u];
        memset(cds, 'N', phase);
        memcpy(cds + phase, unpadded, sizeof unpadded - 1u);
        uint64_t offset = 0u, post_offset = 4u;
        uint32_t length = sizeof unpadded - 1u + phase;
        uint32_t pre_length = layout ? 4u : 0u, post_length = layout ? 5u : 0u;
        uint8_t table = DUCKVEP_CODON_TABLE_STANDARD;
        duckvep_sequence_pool_t seq = {.cds_bytes = cds, .cds_bytes_len = length,
            .cds_offset = &offset, .cds_length = &length, .codon_table = &table, .transcript_count = 1u};
        if (layout) {
            seq.flank_bytes = flanks; seq.flank_bytes_len = sizeof flanks - 1u;
            seq.pre_cds_offset = &offset; seq.pre_cds_length = &pre_length;
            seq.post_cds_offset = &post_offset; seq.post_cds_length = &post_length;
            seq.flanks_complete = 1u;
        }
        duckvep_model_t *model = NULL;
        duckvep_workspace_t *workspace = NULL;
        duckvep_options_t *options = NULL;
        duckvep_error_t error = {0};
        duckvep_options_init_t init = {.distances_are_explicit = 1u,
            .compatibility_profile = DUCKVEP_COMPAT_VEP_116};
        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &seq, NULL, &model, &error));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &workspace, &error));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &options, &error));
        uint8_t genomic[30];
        for (size_t i = 0u; i < sizeof genomic; i++) {
            uint8_t base = unpadded[strand > 0 ? i : sizeof genomic - 1u - i];
            genomic[i] = strand > 0 ? base : base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
        }
        duckvep_haplotype_edit_t edits[32];
        uint8_t alt_cds[80], ref_peptide[32], alt_peptide[32];
        duckvep_delta_scratch_t scratch = {edits, 32u, alt_cds, sizeof alt_cds,
            ref_peptide, sizeof ref_peptide, alt_peptide, sizeof alt_peptide};
        for (uint16_t rl = 1u; rl <= 4u; rl++) for (uint16_t al = 1u; al <= 4u; al++)
            for (unsigned draw = 0u; draw < 64u; draw++) {
                random = random * UINT32_C(1664525) + UINT32_C(1013904223);
                uint32_t pos = 11u + random % (31u - rl), end = pos + rl - 1u;
                uint8_t alleles[8];
                memcpy(alleles, genomic + pos - 11u, rl);
                if (layout) {
                    /* Select every wholly exonic REF placement across both coding
                     * exons; their junction splits a codon, not a source record. */
                    uint32_t slot = random % (32u - 2u * rl);
                    uint32_t first = slot < 12u - rl ? slot : slot + rl - 1u;
                    uint32_t a, b;
                    ASSERT(proj_brute_cdna_to_genomic(&s, first + 5u, &a, NULL));
                    ASSERT(proj_brute_cdna_to_genomic(&s, first + rl + 4u, &b, NULL));
                    pos = a < b ? a : b; end = pos + rl - 1u;
                    for (unsigned j = 0u; j < rl; j++) {
                        uint8_t base = unpadded[first + (strand > 0 ? j : rl - 1u - j)];
                        alleles[j] = strand > 0 ? base
                            : base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
                    }
                }
                for (unsigned i = 0u; i < al; i++) {
                    random = random * UINT32_C(1664525) + UINT32_C(1013904223);
                    alleles[rl + i] = (uint8_t)"ACGT"[(random >> 16u) & 3u];
                }
                if (rl == al && !memcmp(alleles, alleles + rl, rl))
                    alleles[rl] = alleles[rl] == 'A' ? 'C' : 'A';
                duckvep_event_t event = {0};
                ASSERT(duckvep_event_prepare_small(pos, alleles, rl, alleles + rl, al, &event));
                event.chrom_id = s.chrom;
                uint32_t ro = 0u, ao = rl;
                duckvep_variant_batch_t variant = {.chrom_id = &s.chrom, .pos1 = &pos, .end1 = &end,
                    .ref_offset = &ro, .alt_offset = &ao, .ref_length = &rl, .alt_length = &al,
                    .allele_bytes = alleles, .allele_bytes_len = (size_t)rl + al,
                    .variant_kind = &event.kind, .count = 1u};
                struct annotation_observer_capture cursor_capture = {.expected_batch = &variant};
                struct annotation_observer_capture direct_capture = {.expected_batch = &variant};
                duckvep_annotate_cursor_t *cursor = NULL;
                duckvep_consequence_t row;
                duckvep_result_builder_t result;
                duckvep_result_builder_init(&result, &row, 1u);
                ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(
                    model, &variant, options, workspace, &cursor, &error));
                duckvep_annotate_cursor_set_observer(cursor, annotation_observer_capture_row, &cursor_capture);
                duckvep_status_t filled = duckvep_annotate_cursor_fill(cursor, &result, &error);
                ASSERT(filled == DUCKVEP_OK || filled == DUCKVEP_ERR_RESULT_FULL);
                ASSERT(result.count <= 1u);
                ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_pair_observed(model, &variant, &event, 0u,
                    &scratch, NULL, annotation_observer_capture_row, &direct_capture, &error));
                ASSERT_MEM_EQ(cursor_capture.rows, direct_capture.rows, sizeof cursor_capture.rows);
                ASSERT_MEM_EQ(cursor_capture.events, direct_capture.events, sizeof cursor_capture.events);
                ASSERT_MEM_EQ(cursor_capture.deltas, direct_capture.deltas, sizeof cursor_capture.deltas);
                ASSERT_MEM_EQ(&cursor_capture, &direct_capture, sizeof cursor_capture);
                duckvep_haplotype_edit_t projected;
                ASSERT_EQ(DUCKVEP_CDS_EDIT_OK, duckvep_variant_cds_edit_build(
                    &s.tx, &s.ex, &seq, &variant, 0u, 0u, (int8_t)strand, &projected));
                duckvep_haplotype_edit_t unchanged = projected;
                struct annotation_observer_capture reused = {.expected_batch = &variant};
                ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_pair_observed(model, &variant, &event, 0u,
                    &scratch, &projected, annotation_observer_capture_row, &reused, &error));
                ASSERT_MEM_EQ(&direct_capture, &reused, sizeof reused);
                ASSERT_MEM_EQ(&unchanged, &projected, sizeof projected);
                duckvep_result_builder_reset(&result);
                ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_fill(cursor, &result, &error));
                ASSERT_EQ(0u, result.count);
                duckvep_annotate_cursor_close(cursor);
                compared++;
                ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_annotate_pair_observed(model, &variant, &event, 1u,
                    &scratch, NULL, annotation_observer_capture_row, &direct_capture, &error));
                if (draw == 0u) {
                    struct annotation_observer_capture rejected = {.expected_batch = &variant};
                    for (unsigned fault = 0u; fault < 8u; fault++) {
                        duckvep_event_t bad = event;
                        if (fault == 0u) bad.chrom_id++;
                        if (fault == 1u) bad.raw_start1++;
                        if (fault == 2u) bad.raw_end1++;
                        if (fault == 3u) bad.ref_diff_length = UINT16_MAX;
                        if (fault == 4u) bad.alt_diff_length = UINT16_MAX;
                        if (fault == 5u) bad.feature_allele_offset = UINT16_MAX;
                        if (fault == 6u) bad.anchor_ref_offset = UINT16_MAX;
                        if (fault == 7u) bad.has_mate = 1u;
                        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_annotate_pair_observed(
                            model, &variant, &bad, 0u, &scratch, NULL,
                            annotation_observer_capture_row, &rejected, &error));
                    }
                    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_annotate_pair_observed(
                        model, &variant, NULL, 0u, &scratch, NULL,
                        annotation_observer_capture_row, &rejected, &error));
                    variant.allele_bytes_len--;
                    ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE, duckvep_annotate_pair_observed(
                        model, &variant, &event, 0u, &scratch, NULL,
                        annotation_observer_capture_row, &rejected, &error));
                    variant.allele_bytes_len++;
                    for (unsigned fault = 0u; fault < 6u; fault++) {
                        duckvep_haplotype_edit_t bad = projected;
                        if (fault == 0u) bad.cds_start = 0u;
                        if (fault == 1u) bad.variant_strand = -1;
                        if (fault == 2u) bad.ref_len++;
                        if (fault == 3u) bad.alt_len++;
                        if (fault == 4u) bad.alt = bad.alt ? NULL : alleles;
                        if (fault == 5u) bad.cds_start = UINT32_MAX;
                        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_annotate_pair_observed(
                            model, &variant, &event, 0u, &scratch, &bad,
                            annotation_observer_capture_row, &rejected, &error));
                    }
                    ASSERT_EQ(0u, rejected.count);
                }
            }
        duckvep_options_close(options);
        duckvep_workspace_close(workspace);
        duckvep_model_close(model);
    }
    ASSERT_EQ(12288u, compared);
    PASS();
}

TEST annotation_pair_facts_share_projection_and_coding_state(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {205u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {
        DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING
    };
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {2u};
    static const uint32_t cds_start[1] = {100u};
    static const uint32_t cds_end[1] = {108u};
    static const uint32_t exon_start[2] = {100u, 200u};
    static const uint32_t exon_end[2] = {108u, 205u};
    static const uint32_t cdna_start[2] = {1u, 10u};
    static const uint32_t cdna_end[2] = {9u, 15u};
    static const int8_t phase[2] = {0, -1};
    static const uint8_t cds_bytes[9] = {
        'A','T','G','A','A', 'A','T','A','A'
    };
    static const uint64_t cds_offset[1] = {0u};
    static const uint32_t cds_length[1] = {9u};
    static const uint8_t codon_table[1] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3] = {95u, 103u, 150u};
    static const uint8_t vkind[3] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t alleles[6] = {'A','T', 'A','G', 'A','T'};
    static const uint32_t ref_offset[3] = {0u, 2u, 4u};
    static const uint32_t alt_offset[3] = {1u, 3u, 5u};
    static const uint16_t allele_length[3] = {1u, 1u, 1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t sequences;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_error_t error;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t builder;
    struct annotation_observer_capture capture;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&sequences, 0, sizeof sequences);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    memset(&capture, 0, sizeof capture);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 2u;
    sequences.cds_bytes = cds_bytes;
    sequences.cds_bytes_len = sizeof cds_bytes;
    sequences.cds_offset = cds_offset; sequences.cds_length = cds_length;
    sequences.codon_table = codon_table; sequences.transcript_count = 1u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vpos;
    variants.variant_kind = vkind; variants.ref_offset = ref_offset;
    variants.alt_offset = alt_offset; variants.ref_length = allele_length;
    variants.alt_length = allele_length; variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles; variants.count = 3u;
    capture.expected_batch = &variants;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, &sequences, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(
        model, &variants, options, workspace, &cursor, &error));
    duckvep_annotate_cursor_set_observer(
        cursor, annotation_observer_capture_row, &capture);
    duckvep_result_builder_init(&builder, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_cursor_fill(cursor, &builder, &error));
    ASSERT_EQ(3u, builder.count);
    ASSERT_EQ(3u, capture.count);

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT,
              capture.transcript_edit_status[0]);
    ASSERT_EQ(0u, capture.transcript_edit_present[0]);
    ASSERT_EQ(0u, capture.delta_present[0]);
    ASSERT_EQ(0u, capture.coding_context_valid[0]);

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              capture.transcript_edit_status[1]);
    ASSERT_EQ(1u, capture.transcript_edit_present[1]);
    ASSERT_EQ(4u, capture.transcript_edit_first_cdna1[1]);
    ASSERT_EQ(1u, capture.transcript_edit_ref_length[1]);
    ASSERT_EQ(1u, capture.transcript_edit_alt_length[1]);
    ASSERT_EQ(1u, capture.delta_present[1]);
    ASSERT_EQ_FMT((uint8_t)DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  capture.coding_context_status[1], "%u");
    ASSERT_EQ(1u, capture.coding_context_valid[1]);
    ASSERT_EQ(4, rows[1].cds_pos);

    ASSERT_EQ(DUCKVEP_TRANSCRIPT_EDIT_OK,
              capture.transcript_edit_status[2]);
    ASSERT_EQ(1u, capture.transcript_edit_present[2]);
    ASSERT_EQ(0u, capture.delta_present[2]);
    ASSERT_EQ(0u, capture.coding_context_valid[2]);
    ASSERT_EQ_FMT((uint32_t)DUCKVEP_REGION_INTRON,
                  rows[2].region_mask, "%u");

    duckvep_annotate_cursor_close(cursor);
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* The cumulative HGVS adapter observes each consequence synchronously because
 * transcript traces borrow worker scratch. A full one-row output chunk must
 * pause only after that row and its callback are complete, then resume at the
 * next transcript/regulatory/motif object without replaying or dropping an
 * observation. */
TEST annotation_observer_resumes_across_transcript_and_interval_rows(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {150u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t zero32[1] = {0u};
    static const uint16_t zero16[1] = {0u};
    static const uint16_t fchrom[2] = {0u, 0u};
    static const uint32_t fstart[2] = {100u, 100u};
    static const uint32_t fend[2] = {150u, 150u};
    static const uint8_t fkind[2] = {
        (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
        (uint8_t)DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE
    };
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {110u, 111u};
    static const uint32_t vend[2] = {110u, 111u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_interval_feature_model_t features;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_error_t error;
    duckvep_consequence_t row;
    duckvep_result_builder_t builder;
    struct annotation_observer_capture capture;
    size_t emitted = 0u;
    size_t fill_count = 0u;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&features, 0, sizeof features);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    memset(&capture, 0, sizeof capture);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 1u;
    features.chrom_id = fchrom; features.start1 = fstart; features.end1 = fend;
    features.kind = fkind; features.feature_count = 2u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind; variants.count = 2u;
    capture.expected_batch = &variants;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(
        &tx, &exons, NULL, &features, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &workspace, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(
        model, &variants, options, workspace, &cursor, &error));
    duckvep_annotate_cursor_set_observer(
        cursor, annotation_observer_capture_row, &capture);

    while (!duckvep_annotate_cursor_done(cursor)) {
        duckvep_status_t status;
        size_t before = capture.count;

        duckvep_result_builder_init(&builder, &row, 1u);
        status = duckvep_annotate_cursor_fill(cursor, &builder, &error);
        ASSERT(status == DUCKVEP_OK || status == DUCKVEP_ERR_RESULT_FULL);
        ASSERT_EQ(duckvep_result_builder_count(&builder), capture.count - before);
        if (duckvep_result_builder_count(&builder) == 1u) {
            ASSERT(emitted < 6u);
            ASSERT(consequence_rows_equal(&row, &capture.rows[emitted]));
            emitted++;
        }
        fill_count++;
        ASSERT(fill_count <= 7u);
    }

    ASSERT_EQ(6u, emitted);
    ASSERT_EQ(6u, capture.count);
    ASSERT_EQ(7u, fill_count);
    for (emitted = 0u; emitted < 6u; emitted++) {
        uint32_t expected_variant = (uint32_t)(emitted / 3u);
        uint8_t expected_object = (uint8_t)(emitted % 3u);

        ASSERT_EQ(expected_variant, capture.rows[emitted].variant_idx);
        ASSERT_EQ(expected_object, capture.rows[emitted].overlap_object_kind);
        if (expected_object == (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT) {
            ASSERT_EQ(1u, capture.trace_present[emitted]);
            ASSERT_EQ(vpos[expected_variant], capture.events[emitted].raw_start1);
            ASSERT_EQ((uint8_t)DUCKVEP_KIND_SNV, capture.events[emitted].kind);
        } else {
            ASSERT_EQ(0u, capture.trace_present[emitted]);
            ASSERT_EQ((uint32_t)(expected_object - 1u),
                      capture.rows[emitted].interval_feature_idx);
        }
    }

    duckvep_annotate_cursor_close(cursor);
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

TEST sorted_point_cursor_survives_tiles_and_resets_on_rewind(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {310u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {3u};
    static const uint32_t zero[1] = {0u};
    static const uint32_t es[3] = {100u, 200u, 300u};
    static const uint32_t ee[3] = {120u, 220u, 310u};
    uint16_t vchrom[3] = {0u, 0u, 0u};
    uint32_t vpos[3];
    uint32_t vend[3];
    uint8_t vkind[3] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_model_t *other_model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[3];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = zero; tx.cds_end1 = zero;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 3u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend;
    v.variant_kind = vkind;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &ex, NULL, NULL, &other_model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));

    vpos[0] = vend[0] = 105u;
    vpos[1] = vend[1] = 150u;
    v.count = 2u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(other_model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(65u, err.where_code);
    duckvep_result_builder_reset(&rb);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[1].region_mask);

    vpos[0] = vend[0] = 205u;
    vpos[1] = vend[1] = 250u;
    vpos[2] = vend[2] = 305u;
    v.count = 3u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[1].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[2].region_mask);

    /* A new run may rewind. The workspace detects that boundary and invalidates
     * its per-transcript ranks lazily, rather than carrying exon 3 into exon 1. */
    vpos[0] = vend[0] = 110u;
    v.count = 1u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(other_model);
    duckvep_model_close(model);
    PASS();
}

TEST padded_snv_rewind_uses_vep_feature_span(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {220u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t zero[1] = {0u};
    static const uint32_t es[2] = {100u, 200u};
    static const uint32_t ee[2] = {120u, 220u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {100u, 101u};
    static const uint32_t vend[2] = {210u, 105u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint32_t roff[2] = {0u, 222u};
    static const uint32_t aoff[2] = {111u, 227u};
    static const uint16_t rlen[2] = {111u, 5u};
    static const uint16_t alen[2] = {111u, 5u};
    uint8_t alleles[232];
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(alleles, 'A', sizeof alleles);
    alleles[aoff[0] + alen[0] - 1u] = 'C'; /* effect position 210 */
    alleles[aoff[1] + alen[1] - 1u] = 'C'; /* effect position 105 */
    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = zero; tx.cds_end1 = zero;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend;
    v.variant_kind = vkind; v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen;
    v.allele_bytes = alleles; v.allele_bytes_len = sizeof alleles;
    v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_EXON | DUCKVEP_REGION_INTRON),
              rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[1].region_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST span_cursor_resets_after_nonmonotone_tile_skips_transcript(void) {
    enum {
        LONG_PREFIX = 60000,
        LONG_REF_LENGTH = LONG_PREFIX + 2,
        LONG_ALT_LENGTH = LONG_PREFIX + 1,
        FIRST_ALLELE_BYTES = LONG_REF_LENGTH + LONG_ALT_LENGTH + 4
    };
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {60000u};
    static const uint32_t tend[1] = {61100u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t zero[1] = {0u};
    static const uint32_t es[2] = {60000u, 61000u};
    static const uint32_t ee[2] = {60100u, 61100u};
    static uint8_t first_alleles[FIRST_ALLELE_BYTES];
    static const uint8_t second_alleles[8] = {
        'A', 'A', 'C', 'C', 'A', 'A', 'C', 'C'
    };
    static const uint16_t chrom[2] = {0u, 0u};
    static const uint32_t first_pos[2] = {1000u, 10000u};
    static const uint32_t first_end[2] = {61001u, 10001u};
    static const uint8_t first_kind[2] = {
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint32_t first_ref_offset[2] = {
        0u, LONG_REF_LENGTH + LONG_ALT_LENGTH
    };
    static const uint32_t first_alt_offset[2] = {
        LONG_REF_LENGTH, LONG_REF_LENGTH + LONG_ALT_LENGTH + 2u
    };
    static const uint16_t first_ref_length[2] = {
        LONG_REF_LENGTH, 2u
    };
    static const uint16_t first_alt_length[2] = {
        LONG_ALT_LENGTH, 2u
    };
    static const uint32_t second_pos[2] = {56000u, 60050u};
    static const uint32_t second_end[2] = {56001u, 60051u};
    static const uint8_t second_kind[2] = {
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint32_t second_ref_offset[2] = {0u, 4u};
    static const uint32_t second_alt_offset[2] = {2u, 6u};
    static const uint16_t second_length[2] = {2u, 2u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t first;
    duckvep_variant_batch_t second;
    duckvep_options_init_t init;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *reused = NULL;
    duckvep_workspace_t *fresh = NULL;
    duckvep_consequence_t first_rows[2];
    duckvep_consequence_t reused_rows[2];
    duckvep_consequence_t fresh_rows[2];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(first_alleles, 'A', sizeof first_alleles);
    first_alleles[LONG_PREFIX] = 'G';
    first_alleles[LONG_PREFIX + 1u] = 'C';
    first_alleles[LONG_REF_LENGTH + LONG_PREFIX] = 'T';
    first_alleles[first_ref_offset[1]] = 'A';
    first_alleles[first_ref_offset[1] + 1u] = 'A';
    first_alleles[first_alt_offset[1]] = 'C';
    first_alleles[first_alt_offset[1] + 1u] = 'C';

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&first, 0, sizeof first);
    memset(&second, 0, sizeof second);
    memset(&init, 0, sizeof init);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = zero; tx.cds_end1 = zero;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;

    first.chrom_id = chrom; first.pos1 = first_pos; first.end1 = first_end;
    first.variant_kind = first_kind; first.ref_offset = first_ref_offset;
    first.alt_offset = first_alt_offset; first.ref_length = first_ref_length;
    first.alt_length = first_alt_length; first.allele_bytes = first_alleles;
    first.allele_bytes_len = sizeof first_alleles; first.count = 2u;

    second.chrom_id = chrom; second.pos1 = second_pos; second.end1 = second_end;
    second.variant_kind = second_kind; second.ref_offset = second_ref_offset;
    second.alt_offset = second_alt_offset; second.ref_length = second_length;
    second.alt_length = second_length; second.allele_bytes = second_alleles;
    second.allele_bytes_len = sizeof second_alleles; second.count = 2u;

    init.distances_are_explicit = 1u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &reused, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &fresh, &err));

    duckvep_result_builder_init(&rb, first_rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &first, opts, reused, &rb, &err));
    ASSERT_EQ(1u, rb.count);

    duckvep_result_builder_init(&rb, reused_rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &second, opts, reused, &rb, &err));
    ASSERT_EQ(1u, rb.count);
    duckvep_result_builder_init(&rb, fresh_rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &second, opts, fresh, &rb, &err));
    ASSERT_EQ(1u, rb.count);
    ASSERT(consequence_rows_equal(&reused_rows[0], &fresh_rows[0]));
    ASSERT_EQ(1u, reused_rows[0].variant_idx);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, reused_rows[0].region_mask);

    duckvep_workspace_close(fresh);
    duckvep_workspace_close(reused);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* The resumable cursor is a transport primitive: arbitrary output-buffer splits
 * must equal one full annotate_tile call over the same tile. This property uses
 * the same random zero-exon scenes as the composition oracle and varies the chunk
 * cap from the generated dimensions, so it proves cursor state across variant,
 * active-set, and output-buffer boundaries without restating the sweep. */
static enum theft_trial_res prop_annotate_cursor_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    duckvep_exon_model_t exons;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    static duckvep_consequence_t full[KPROP_MAX_PAIRS];
    static duckvep_consequence_t got[KPROP_MAX_PAIRS];
    duckvep_consequence_t chunk[9];
    duckvep_result_builder_t rb;
    duckvep_status_t st;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    size_t full_n, got_n = 0u, i;
    size_t chunk_cap;
    const uint32_t HALO = KPROP_SWEEP_HALO;
    (void)t;

    memset(&exons, 0, sizeof exons);
    memset(&err, 0, sizeof err);
    memset(&init, 0, sizeof init);
    init.upstream_dist = HALO;
    init.downstream_dist = HALO;
    init.halo = HALO;
    chunk_cap = 1u + ((s->v.count + s->tx.transcript_count + (size_t)s->halo) % 9u);

    if (duckvep_model_open(&s->tx, &exons, NULL, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(&init, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, full, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    full_n = duckvep_result_builder_count(&rb);

    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_result_builder_init(&rb, chunk, chunk_cap);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL && duckvep_result_builder_count(&rb) == 0u) {
            res = THEFT_TRIAL_FAIL; goto done;
        }
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (got_n >= KPROP_MAX_PAIRS) { res = THEFT_TRIAL_FAIL; goto done; }
            got[got_n++] = chunk[i];
        }
    }

    if (got_n != full_n) { res = THEFT_TRIAL_FAIL; goto done; }
    for (i = 0u; i < full_n; i++) {
        if (!consequence_rows_equal(&got[i], &full[i])) { res = THEFT_TRIAL_FAIL; goto done; }
    }

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor output splits == one annotate_tile";
    cfg.prop1 = prop_annotate_cursor_matches_tile;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

TEST coding_delta_so_and_independent_class_gates_are_separate(void) {
    static const struct { size_t offset; uint64_t pre, so; } fields[] = {
#define FIELD(member, pre, so) {offsetof(duckvep_sequence_delta_t, member), DUCKVEP_PRE(pre), DUCKVEP_SO(so)}
        FIELD(synonymous, DUCKVEP_PRE_SYNONYMOUS, DUCKVEP_SO_SYNONYMOUS),
        FIELD(missense, DUCKVEP_PRE_MISSENSE, DUCKVEP_SO_MISSENSE),
        FIELD(stop_gained, DUCKVEP_PRE_STOP_GAINED, DUCKVEP_SO_STOP_GAINED),
        FIELD(stop_lost, DUCKVEP_PRE_STOP_LOST, DUCKVEP_SO_STOP_LOST),
        FIELD(stop_retained, DUCKVEP_PRE_STOP_RETAINED, DUCKVEP_SO_STOP_RETAINED),
        FIELD(start_lost, DUCKVEP_PRE_START_LOST, DUCKVEP_SO_START_LOST),
        FIELD(start_retained, DUCKVEP_PRE_START_RETAINED, DUCKVEP_SO_START_RETAINED),
        FIELD(frameshift, DUCKVEP_PRE_FRAMESHIFT, DUCKVEP_SO_FRAMESHIFT),
        FIELD(inframe_deletion, DUCKVEP_PRE_INFRAME_DELETION, DUCKVEP_SO_INFRAME_DELETION),
        FIELD(inframe_insertion, DUCKVEP_PRE_INFRAME_INSERTION, DUCKVEP_SO_INFRAME_INSERTION),
        FIELD(protein_altering, DUCKVEP_PRE_PROTEIN_ALTERING, DUCKVEP_SO_PROTEIN_ALTERING),
        FIELD(coding_unknown, DUCKVEP_PRE_CODING_UNKNOWN, DUCKVEP_SO_CODING_SEQUENCE),
        FIELD(partial_codon, DUCKVEP_PRE_PARTIAL_CODON, DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON)
#undef FIELD
    };
    const uint64_t classes[] = {DUCKVEP_PRE(DUCKVEP_PRE_SNP),
        DUCKVEP_PRE(DUCKVEP_PRE_INSERTION), DUCKVEP_PRE(DUCKVEP_PRE_DELETION)};
    ASSERT_EQ(0u, duckvep_effect_eval_coding_delta(NULL));
    for (uint32_t bits = 0u; bits < (1u << 13u); bits++) for (unsigned valid = 0u; valid < 2u; valid++) {
        duckvep_sequence_delta_t delta = {0};
        uint64_t expected_so = 0u;
        delta.valid = (uint8_t)valid;
        for (unsigned i = 0u; i < 13u; i++) if (bits & (1u << i)) {
            *((uint8_t *)&delta + fields[i].offset) = 1u;
            expected_so |= fields[i].so;
        }
        ASSERT_EQ(valid ? expected_so : 0u, duckvep_effect_eval_coding_delta(&delta));
        for (unsigned shape = 0u; shape < 8u; shape++) {
            duckvep_effect_ctx_t ctx = {0};
            ctx.pre_bits = DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM);
            for (unsigned i = 0u; i < 3u; i++) if (shape & (1u << i)) ctx.pre_bits |= classes[i];
            uint64_t expected = ctx.pre_bits;
            if (valid) {
                expected |= DUCKVEP_PRE(DUCKVEP_PRE_DELTA);
                for (unsigned i = 0u; i < 13u; i++) if (bits & (1u << i)) {
                    uint64_t pre = fields[i].pre;
                    if (pre == DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE) && (shape & 6u)) continue;
                    if (pre == DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT) && (shape & 1u)) continue;
                    if (pre == DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION) && !(shape & 2u)) continue;
                    if (pre == DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION) && !(shape & 4u)) continue;
                    expected |= pre;
                }
            }
            duckvep_effect_ctx_apply_delta(&ctx, &delta);
            ASSERT_EQ(expected, ctx.pre_bits);
        }
    }
    PASS();
}

TEST event_length_delta_pre_bits_follow_trimmed_alleles(void) {
    duckvep_effect_ctx_t ctx;
    duckvep_event_t event;
    duckvep_sequence_delta_t delta;

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 1u;
    event.alt_diff_length = 4u;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    /* VEP classifies ordinary alleles from the complete uploaded feature
     * lengths, before transcript projection or common-affix trimming. A
     * feature spanning an intron can therefore remain SNP-class even when its
     * outer mapped CDS replacement changes frame. */
    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    memset(&delta, 0, sizeof delta);
    event.kind = (uint8_t)DUCKVEP_KIND_MNV;
    event.ref_diff_length = 1u;
    event.alt_diff_length = 4u;
    event.feature_length_relation =
        (uint8_t)DUCKVEP_FEATURE_LENGTH_EQUAL;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SNP)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);
    delta.valid = 1u;
    delta.missense = 1u;
    delta.frameshift = 1u;
    delta.inframe_insertion = 1u;
    delta.inframe_deletion = 1u;
    duckvep_effect_ctx_apply_delta(&ctx, &delta);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    memset(&delta, 0, sizeof delta);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 4u;
    event.alt_diff_length = 4u;
    event.feature_length_relation =
        (uint8_t)DUCKVEP_FEATURE_LENGTH_INCREASE;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    delta.valid = 1u;
    delta.missense = 1u;
    delta.inframe_insertion = 1u;
    duckvep_effect_ctx_apply_delta(&ctx, &delta);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION)) != 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    memset(&delta, 0, sizeof delta);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 4u;
    event.alt_diff_length = 4u;
    event.feature_length_relation =
        (uint8_t)DUCKVEP_FEATURE_LENGTH_DECREASE;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u);
    delta.valid = 1u;
    delta.missense = 1u;
    delta.inframe_deletion = 1u;
    duckvep_effect_ctx_apply_delta(&ctx, &delta);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION)) != 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 4u;
    event.alt_diff_length = 1u;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_MNV;
    event.ref_diff_length = 2u;
    event.alt_diff_length = 2u;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);

    PASS();
}

/* VariationEffect::feature_ablation is not structural-only. A literal allele
 * whose normalized deletion fact contains the complete transcript must enter
 * the same tier-1 transcript_ablation rule as a symbolic deletion. */
TEST ordinary_complete_feature_deletion_is_ablation_known_scene(void) {
    duckvep_effect_ctx_t ctx;
    duckvep_event_t event;

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    ctx.region_state.complete_overlap_feature = 1u;
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 7u;
    event.alt_diff_length = 1u;

    duckvep_effect_ctx_apply_event(NULL, &ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_ABLATION)) != 0u);
    duckvep_effect_ctx_finalize(&ctx);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_ABLATION),
              duckvep_effect_eval(ctx.pre_bits));
    PASS();
}

/* BaseVariationFeatureOverlapAllele::coding_unknown returns false for complete
 * transcript overlap before peptide predicates are considered. An explicit
 * unknown fact from a sequence delta must therefore also be suppressed. */
TEST complete_feature_overlap_suppresses_explicit_coding_unknown_known_scene(void) {
    duckvep_effect_ctx_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&ctx, 0, sizeof ctx);
    memset(&delta, 0, sizeof delta);
    ctx.region_state.complete_overlap_feature = 1u;
    ctx.pre_bits = DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                   DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                   DUCKVEP_PRE(DUCKVEP_PRE_UTR5) |
                   DUCKVEP_PRE(DUCKVEP_PRE_UTR3);
    delta.valid = 1u;
    delta.coding_unknown = 1u;

    duckvep_effect_ctx_apply_delta(&ctx, &delta);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN)) != 0u);
    duckvep_effect_ctx_finalize(&ctx);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN)) == 0u);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
              DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
              duckvep_effect_eval(ctx.pre_bits));
    PASS();
}

/* VEP runs its four-comparison UTR overlap predicate on empty endpoint
 * intervals. For an equal-length uploaded span containing a one-exon coding
 * transcript whose CDS equals the transcript, both strands consequently emit
 * both UTR terms. Complete overlap also makes coding_unknown false. */
TEST annotate_complete_equal_length_span_keeps_empty_endpoint_utrs_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1008u, 2008u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2000u};
    static const uint32_t tcds_e[2] = {1008u, 2008u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1008u, 2008u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {9u, 9u};
    static const int8_t ephase[2] = {0, 0};
    static const uint16_t vchrom[2] = {0u, 1u};
    static const uint32_t vpos[2] = {999u, 1999u};
    static const uint32_t vend[2] = {1009u, 2009u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t abytes[] =
        "AAAAAAAAAAA" "AAAAACAAAAA"
        "CCCCCCCCCCC" "CCCCCGCCCCC";
    static const uint32_t roff[2] = {0u, 22u};
    static const uint32_t aoff[2] = {11u, 33u};
    static const uint16_t alen[2] = {11u, 11u};
    static const uint64_t expected =
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR);
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s;
    tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind; variants.allele_bytes = abytes;
    variants.allele_bytes_len = sizeof abytes - 1u;
    variants.ref_offset = roff; variants.alt_offset = aoff;
    variants.ref_length = alen; variants.alt_length = alen;
    variants.count = 2u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&builder, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &builder, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&builder));
    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected, rows[i].consequence_mask);
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

struct kprop_complete_overlap {
    int8_t strand;
    uint16_t transcript_length;
    uint16_t left_padding;
    uint16_t right_padding;
};

static enum theft_alloc_res kprop_complete_overlap_alloc(
    struct theft *t, void *env, void **instance) {
    struct kprop_complete_overlap *s;
    uint8_t long_side;
    (void)env;

    s = (struct kprop_complete_overlap *)calloc(1u, sizeof *s);
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->strand = theft_random_bits(t, 1) == 0u ? 1 : -1;
    s->transcript_length = (uint16_t)(3u * (1u + kprop_bounded(t, 100u)));
    s->left_padding = (uint16_t)(1u + kprop_bounded(t, 64u));
    s->right_padding = (uint16_t)kprop_bounded(t, 65u);
    if (kprop_bounded(t, 8u) == 0u) {
        long_side = (uint8_t)theft_random_bits(t, 1);
        if (long_side == 0u) {
            s->left_padding = (uint16_t)(5001u + kprop_bounded(t, 1000u));
        } else {
            s->right_padding = (uint16_t)(5001u + kprop_bounded(t, 1000u));
        }
    }
    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_complete_overlap_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static struct theft_type_info kprop_complete_overlap_info = {
    .alloc = kprop_complete_overlap_alloc,
    .free = kprop_complete_overlap_free,
};

static struct {
    uint32_t forward;
    uint32_t reverse;
    uint32_t right_endpoint;
    uint32_t over_5000;
} g_complete_overlap_cov;

static enum theft_trial_res prop_complete_literal_spans_match_vep_source_semantics(
    struct theft *t, void *arg1) {
    const struct kprop_complete_overlap *s =
        (const struct kprop_complete_overlap *)arg1;
    const uint32_t transcript_start = 10000u;
    const uint32_t transcript_end =
        transcript_start + (uint32_t)s->transcript_length - 1u;
    const uint32_t raw_start = transcript_start - (uint32_t)s->left_padding;
    const size_t raw_length = (size_t)s->left_padding +
        (size_t)s->transcript_length + (size_t)s->right_padding;
    const uint32_t raw_end = raw_start + (uint32_t)raw_length - 1u;
    uint16_t tchrom[1] = {0u};
    uint32_t tstart[1] = {transcript_start};
    uint32_t tend[1] = {transcript_end};
    int8_t tstrand[1] = {s->strand};
    uint64_t tflags[1] = {0u};
    uint32_t texoff[1] = {0u};
    uint16_t texcnt[1] = {1u};
    uint32_t tcds_s[1] = {transcript_start};
    uint32_t tcds_e[1] = {transcript_end};
    uint32_t estart[1] = {transcript_start};
    uint32_t eend[1] = {transcript_end};
    uint32_t ecdna_s[1] = {1u};
    uint32_t ecdna_e[1] = {(uint32_t)s->transcript_length};
    int8_t ephase[1] = {0};
    uint16_t vchrom[2] = {0u, 0u};
    uint32_t vpos[2] = {raw_start, raw_start};
    uint32_t vend[2] = {raw_end, raw_end};
    uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_DEL
    };
    uint32_t roff[2];
    uint32_t aoff[2];
    uint16_t rlen[2];
    uint16_t alen[2];
    uint8_t *alleles = NULL;
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    enum theft_trial_res result = THEFT_TRIAL_PASS;
    uint64_t expected_span_mask;
    (void)t;

    alleles = (uint8_t *)malloc(raw_length * 3u + 1u);
    if (alleles == NULL) return THEFT_TRIAL_ERROR;
    memset(alleles, 'A', raw_length * 3u + 1u);
    /* Keep the equal-length event a genuine MNV after shared-edge trimming,
     * while independently varying uploaded padding outside the transcript. */
    alleles[raw_length + (size_t)s->left_padding] = (uint8_t)'C';
    alleles[raw_length + (size_t)s->left_padding +
            (size_t)s->transcript_length - 1u] = (uint8_t)'G';
    roff[0] = 0u;
    aoff[0] = (uint32_t)raw_length;
    roff[1] = (uint32_t)(raw_length * 2u);
    aoff[1] = (uint32_t)(raw_length * 3u);
    rlen[0] = rlen[1] = (uint16_t)raw_length;
    alen[0] = (uint16_t)raw_length;
    alen[1] = 1u;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s;
    tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind; variants.allele_bytes = alleles;
    variants.allele_bytes_len = raw_length * 3u + 1u;
    variants.ref_offset = roff; variants.alt_offset = aoff;
    variants.ref_length = rlen; variants.alt_length = alen;
    variants.count = 2u;

    /* The uploaded feature always extends left of the transcript, satisfying
     * VEP's inverted pre-CDS interval comparison. It satisfies the inverted
     * post-CDS interval only when it also extends right of the transcript. */
    expected_span_mask = s->strand > 0
        ? DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR)
        : DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR);
    if (s->right_padding > 0u) {
        expected_span_mask |= s->strand > 0
            ? DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR)
            : DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR);
    }

    if (duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (duckvep_options_open(NULL, &options, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (duckvep_workspace_open(model, &workspace, &error) != DUCKVEP_OK) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    duckvep_result_builder_init(&builder, rows, 2u);
    if (duckvep_annotate_tile(model, &variants, options, workspace,
                              &builder, &error) != DUCKVEP_OK ||
        duckvep_result_builder_count(&builder) != 2u ||
        rows[0].variant_idx != 0u ||
        rows[0].consequence_mask != expected_span_mask ||
        rows[1].variant_idx != 1u ||
        rows[1].consequence_mask !=
            DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_ABLATION)) {
        result = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (s->strand > 0) g_complete_overlap_cov.forward++;
    else g_complete_overlap_cov.reverse++;
    if (s->right_padding == 0u) g_complete_overlap_cov.right_endpoint++;
    if (raw_length > 5000u) g_complete_overlap_cov.over_5000++;

done:
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    free(alleles);
    return result;
}

TEST complete_literal_spans_match_vep_source_semantics_for_any_scene(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "complete literal spans == VEP complete-overlap source semantics";
    cfg.prop1 = prop_complete_literal_spans_match_vep_source_semantics;
    cfg.type_info[0] = &kprop_complete_overlap_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_complete_overlap_cov, 0, sizeof g_complete_overlap_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_complete_overlap_cov.forward > 0u);
    ASSERT(g_complete_overlap_cov.reverse > 0u);
    ASSERT(g_complete_overlap_cov.right_endpoint > 0u);
    ASSERT(g_complete_overlap_cov.over_5000 > 0u);
    fprintf(stderr,
            "[complete-overlap coverage] forward=%u reverse=%u "
            "right_endpoint=%u over_5000=%u\n",
            g_complete_overlap_cov.forward, g_complete_overlap_cov.reverse,
            g_complete_overlap_cov.right_endpoint,
            g_complete_overlap_cov.over_5000);
    PASS();
}

/* Deterministic SO-mapping anchor: a 2-exon coding transcript on chr0 (exons
 * [1000,1300] & [1700,2000], CDS [1100,1900]) and a single-exon non-coding
 * transcript on chr1 ([1000,1300]). Each variant lands in exactly one structural
 * bucket; the expected consequence_mask + impact are hand-computed from the SO
 * spec, so this cannot pass vacuously and pins every branch of the mapping. */
TEST annotate_structural_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 1000u};
    static const uint32_t tend[2]    = {2000u, 1300u};
    static const int8_t   tstrand[2] = {1, 1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 2u};
    static const uint16_t texcnt[2]  = {2u, 1u};
    static const uint32_t tcds_s[2]  = {1100u, 0u};
    static const uint32_t tcds_e[2]  = {1900u, 0u};
    static const uint32_t estart[3]  = {1000u, 1700u, 1000u};
    static const uint32_t eend[3]    = {1300u, 2000u, 1300u};

    static const uint16_t vchrom[9] = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 1u};
    static const uint32_t vpos[9]   = {900u, 1050u, 1200u, 1301u, 1500u, 1750u, 1950u, 2100u, 1150u};
    static const uint8_t  vkind[9]  = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};

    static const uint64_t exp_mask[9] = {
        DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE),                                       /* 900  up   */
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),                                         /* 1050 5'utr*/
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),                                     /* 1200 cds  */
        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),                                       /* 1301 = 1st intronic base (essential donor); splice_donor ALONE (intron_variant suppressed at the essential site, VEP-faithful), HIGH impact flows through annotate_tile */
        DUCKVEP_SO(DUCKVEP_SO_INTRON),                                              /* 1500 intr */
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),                                     /* 1750 cds  */
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),                                         /* 1950 3'utr*/
        DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),                                     /* 2100 down */
        DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON)                          /* 1150 nc exon -> exon term only */
    };
    static const uint8_t exp_impact[9] = {
        DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER,
        DUCKVEP_IMPACT_HIGH,     DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER,
        DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[16];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 3u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 9u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err)); /* defaults */
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 16u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(9u, duckvep_result_builder_count(&rb));

    for (i = 0u; i < 9u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_impact[i], rows[i].impact);
        ASSERT_EQ(-1, rows[i].cds_pos);
        if (i == 2u || i == 5u) {
            ASSERT((rows[i].flags &
                    (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0u);
            ASSERT_EQ(DUCKVEP_SEQUENCE_MISSING, rows[i].sequence_status);
        } else {
            ASSERT_EQ(0u, rows[i].flags &
                         (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
            ASSERT_EQ(DUCKVEP_SEQUENCE_NOT_APPLICABLE,
                      rows[i].sequence_status);
        }
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_padded_small_variants_use_differing_region_topology(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {100u};
    static const uint32_t tend[1]    = {300u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {120u};
    static const uint32_t tcds_e[1]  = {220u};
    static const uint32_t estart[2]  = {100u, 200u};
    static const uint32_t eend[2]    = {150u, 300u};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {150u, 199u, 220u};
    static const uint32_t vend[3]   = {151u, 200u, 221u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_DEL
    };
    static const uint8_t  bytes[9] = {
        'C','G', 'C', /* changed base 151 = first intronic donor base */
        'G','T', 'G', /* changed base 200 = exon start, not acceptor dinucleotide */
        'T','A', 'T'  /* changed base 221 = 3'UTR, not CDS stop anchor */
    };
    static const uint32_t roff[3] = {0u, 3u, 6u};
    static const uint32_t aoff[3] = {2u, 5u, 8u};
    static const uint16_t rlen[3] = {2u, 2u, 2u};
    static const uint16_t alen[3] = {1u, 1u, 1u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR)
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cds_pos);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* End-to-end structural/CNV anchor. Copy direction is explicit: an undirected
 * CNV does not become a gain or loss merely from its interval. Tier-1 complete
 * feature events suppress the ordinary placement terms; contained events retain
 * both their structural predicate and the applicable topology consequence. */
TEST annotate_sv_cnv_known_scene(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {300u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t cds_s[1] = {120u};
    static const uint32_t cds_e[1] = {280u};
    static const uint32_t es[2] = {100u, 250u};
    static const uint32_t ee[2] = {150u, 300u};

    static const uint16_t vchrom[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vstart[6] = {50u, 60u, 125u, 135u, 180u, 260u};
    static const uint32_t vend[6] = {350u, 340u, 130u, 140u, 220u, 265u};
    static const uint8_t vkind[6] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV,
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[6] = {
        DUCKVEP_SV_CNV, DUCKVEP_SV_DELETION, DUCKVEP_SV_CNV,
        DUCKVEP_SV_CNV, DUCKVEP_SV_INVERSION, DUCKVEP_SV_CNV
    };
    static const uint8_t copy_change[6] = {
        DUCKVEP_COPY_CHANGE_GAIN, DUCKVEP_COPY_CHANGE_LOSS,
        DUCKVEP_COPY_CHANGE_GAIN, DUCKVEP_COPY_CHANGE_LOSS,
        DUCKVEP_COPY_CHANGE_UNKNOWN, DUCKVEP_COPY_CHANGE_NEUTRAL
    };
    static const uint64_t expected[6] = {
        DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_AMPLIFICATION),
        DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_ABLATION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_ELONGATION) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
            DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION),
        DUCKVEP_SO(DUCKVEP_SO_INTRON),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    duckvep_error_t err;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.sv_type = sv_type; v.copy_change = copy_change;
    v.count = 6u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(6u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 6u; i++) {
        ASSERT_EQ((uint32_t)i, rows[i].variant_idx);
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(i < 4u ? (uint8_t)DUCKVEP_IMPACT_HIGH
                         : (uint8_t)DUCKVEP_IMPACT_MODIFIER,
                  rows[i].impact);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(-1, rows[i].protein_pos);
    }

    /* Direct kernel callers bypass the adapter tile, so the borrowed-view
     * validator must also reject contradictory structural metadata. */
    {
        static const uint8_t bad_type[1] = {DUCKVEP_SV_DELETION};
        static const uint8_t bad_copy[1] = {DUCKVEP_COPY_CHANGE_GAIN};
        duckvep_variant_batch_t bad = v;
        bad.count = 1u;
        bad.sv_type = bad_type;
        bad.copy_change = bad_copy;
        duckvep_result_builder_reset(&rb);
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
                  duckvep_annotate_tile(model, &bad, opts, ws, &rb, &err));
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_breakend_pairs_keep_local_topology_and_mate_truncation(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {100u, 400u};
    static const uint32_t tend[2] = {300u, 500u};
    static const int8_t strand[2] = {1, 1};
    static const uint64_t flags[2] = {0u, 0u};
    static const uint32_t exoff[2] = {0u, 2u};
    static const uint16_t excnt[2] = {2u, 1u};
    static const uint32_t cds_s[2] = {120u, 0u};
    static const uint32_t cds_e[2] = {280u, 0u};
    static const uint32_t es[3] = {100u, 250u, 400u};
    static const uint32_t ee[3] = {150u, 300u, 500u};
    static const uint16_t vchrom[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vpos[6] = {99u, 150u, 150u, 150u, 179u, 300u};
    static const uint32_t vend[6] = {99u, 150u, 150u, 150u, 179u, 300u};
    static const uint16_t mate_chrom[6] = {0u, 0u, 0u, 0u, 0u, 1u};
    static const uint32_t mate_pos[6] = {180u, 301u, 5301u, 250u, 301u, 450u};
    static const uint8_t kind[6] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV,
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[6] = {
        DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND,
        DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND
    };
    static const uint8_t copy_change[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t pair_variant[7] = {0u, 1u, 2u, 3u, 4u, 5u, 5u};
    static const uint32_t pair_tx[7] = {0u, 0u, 0u, 0u, 0u, 0u, 1u};
    /* At the first intron base, structural predicates have no ordinary term.
     * A close extragenic mate defaults to intergenic before the allele union.
     * A distant mate, an intragenic mate, or a shared intron term must not. */
    static const uint64_t expected[7] = {
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
            DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) | DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) | DUCKVEP_SO(DUCKVEP_SO_INTRON),
        DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION)
    };
    static const uint32_t expected_region[7] = {
        DUCKVEP_REGION_UTR,
        DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE,
        DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE,
        DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE,
        DUCKVEP_REGION_INTRON,
        DUCKVEP_REGION_DOWNSTREAM,
        0u
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t variants;
    duckvep_candidate_pairs_t pairs;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[7];
    duckvep_result_builder_t results;
    duckvep_error_t error;
    size_t row;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&variants, 0, sizeof variants);
    memset(&pairs, 0, sizeof pairs);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = cds_s; tx.cds_end1 = cds_e;
    tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 3u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.mate_chrom_id = mate_chrom; variants.mate_pos1 = mate_pos;
    variants.variant_kind = kind; variants.sv_type = sv_type;
    variants.copy_change = copy_change; variants.count = 6u;
    pairs.variant_idx = pair_variant; pairs.tx_idx = pair_tx; pairs.count = 7u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&results, rows, 7u);
    ASSERT_EQ(DUCKVEP_ERR_UNSUPPORTED,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &results, &error));
    duckvep_result_builder_reset(&results);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_pairs(model, &variants, &pairs, options,
                                     workspace, &results, &error));
    ASSERT_EQ(7u, duckvep_result_builder_count(&results));
    for (row = 0u; row < 7u; row++) {
        ASSERT_EQ(pair_variant[row], rows[row].variant_idx);
        ASSERT_EQ(pair_tx[row], rows[row].tx_idx);
        ASSERT_EQ(expected[row], rows[row].consequence_mask);
        ASSERT_EQ(expected_region[row], rows[row].region_mask);
    }
    {
        static const uint32_t duplicate_variant[2] = {0u, 0u};
        static const uint32_t duplicate_tx[2] = {0u, 0u};
        duckvep_candidate_pairs_t bad_pairs = pairs;

        bad_pairs.variant_idx = duplicate_variant;
        bad_pairs.tx_idx = duplicate_tx;
        bad_pairs.count = 2u;
        duckvep_result_builder_reset(&results);
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
                  duckvep_annotate_pairs(model, &variants, &bad_pairs,
                                         options, workspace, &results, &error));
    }
    {
        duckvep_variant_batch_t missing_mate = variants;
        duckvep_candidate_pairs_t one_pair = pairs;

        missing_mate.mate_chrom_id = NULL;
        one_pair.count = 1u;
        duckvep_result_builder_reset(&results);
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
                  duckvep_annotate_pairs(model, &missing_mate, &one_pair,
                                         options, workspace, &results, &error));
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* RegFeat::annotate_InputBuffer evaluates the local VEP-shifted BND point and
 * every raw mate point. Candidate discovery may conservatively include a
 * feature at the same numeric position on a different contig; the kernel must
 * reject it and deduplicate a feature hit by both endpoints. A local hit keeps
 * the ordinary regulatory/motif term; a mate-only hit takes VEP's exceptional
 * feature_truncation chromosome-breakpoint branch. If that mate-only object's
 * local point is outside but no farther than VEP's 5000-base structural-feature
 * admission distance, the same object also receives intergenic_variant. If
 * both points hit the same object, the local term wins. */
TEST annotate_breakend_interval_feature_pairs_use_both_endpoints(void) {
    static const uint16_t feature_chrom[4] = {0u, 0u, 1u, 1u};
    static const uint32_t feature_start[4] = {100u, 200u, 100u, 400u};
    static const uint32_t feature_end[4] = {120u, 220u, 120u, 450u};
    static const uint8_t feature_kind[4] = {
        DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
        DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
        DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
        DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE
    };
    static const uint16_t chrom[3] = {0u, 0u, 0u};
    static const uint32_t pos[3] = {99u, 130u, 199u};
    static const uint32_t end[3] = {99u, 130u, 199u};
    static const uint16_t mate_chrom[3] = {1u, 0u, 0u};
    static const uint32_t mate_pos[3] = {430u, 110u, 210u};
    static const uint8_t kind[3] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[3] = {
        DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND
    };
    static const uint8_t copy_change[3] = {
        DUCKVEP_COPY_CHANGE_UNKNOWN, DUCKVEP_COPY_CHANGE_UNKNOWN,
        DUCKVEP_COPY_CHANGE_UNKNOWN
    };
    static const uint32_t pair_variant[5] = {0u, 0u, 0u, 1u, 2u};
    static const uint32_t pair_feature[5] = {0u, 2u, 3u, 0u, 1u};
    static const uint32_t expected_variant[4] = {0u, 0u, 1u, 2u};
    static const uint32_t expected_feature[4] = {0u, 3u, 0u, 1u};
    static const uint64_t expected_mask[4] = {
        DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
            DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
        DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE)
    };
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_interval_feature_model_t features;
    duckvep_variant_batch_t variants;
    duckvep_interval_feature_pairs_t pairs;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[5];
    duckvep_result_builder_t results;
    duckvep_error_t error;
    size_t row;

    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&features, 0, sizeof features);
    memset(&variants, 0, sizeof variants);
    memset(&pairs, 0, sizeof pairs);
    memset(&error, 0, sizeof error);
    features.chrom_id = feature_chrom;
    features.start1 = feature_start;
    features.end1 = feature_end;
    features.kind = feature_kind;
    features.feature_count = 4u;
    variants.chrom_id = chrom;
    variants.pos1 = pos;
    variants.end1 = end;
    variants.mate_chrom_id = mate_chrom;
    variants.mate_pos1 = mate_pos;
    variants.variant_kind = kind;
    variants.sv_type = sv_type;
    variants.copy_change = copy_change;
    variants.count = 3u;
    pairs.variant_idx = pair_variant;
    pairs.feature_idx = pair_feature;
    pairs.count = 5u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(
        &transcripts, &exons, NULL, &features, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&results, rows, 5u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_interval_feature_pairs(
        model, &variants, &pairs, options, workspace, &results, &error));
    ASSERT_EQ(4u, duckvep_result_builder_count(&results));
    for (row = 0u; row < 4u; row++) {
        ASSERT_EQ(expected_variant[row], rows[row].variant_idx);
        ASSERT_EQ(expected_feature[row], rows[row].interval_feature_idx);
        ASSERT_EQ(expected_mask[row], rows[row].consequence_mask);
        ASSERT_EQ(UINT32_MAX, rows[row].gene_idx);
    }
    {
        static const uint32_t duplicate_variant[2] = {0u, 0u};
        static const uint32_t duplicate_feature[2] = {0u, 0u};
        duckvep_interval_feature_pairs_t bad = pairs;

        bad.variant_idx = duplicate_variant;
        bad.feature_idx = duplicate_feature;
        bad.count = 2u;
        duckvep_result_builder_reset(&results);
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
                  duckvep_annotate_interval_feature_pairs(
                      model, &variants, &bad, options, workspace, &results,
                      &error));
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* VEP exposes a configurable transcript-direction distance, whose default is
 * 5000, and separately hardcodes StructuralVariationOverlap allele admission
 * to MAX_DISTANCE_FROM_TRANSCRIPT (also 5000). A wider caller window must not
 * widen that structural-allele cap, but ordinary predicates on a mate-created
 * allele still inspect the local feature and may emit a directional term past
 * the fixed cap. Exercise both transcripts and interval features exactly at
 * 5000 and one base beyond under deliberately 10000- and zero-base caller
 * windows. */
TEST annotate_breakend_fixed_admission_is_not_configurable_window(void) {
    static const uint16_t transcript_chrom[1] = {0u};
    static const uint32_t transcript_start[1] = {100u};
    static const uint32_t transcript_end[1] = {250u};
    static const int8_t transcript_strand[1] = {1};
    static const uint64_t transcript_flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {2u};
    static const uint32_t cds_start[1] = {0u};
    static const uint32_t cds_end[1] = {0u};
    static const uint32_t exon_start[2] = {100u, 200u};
    static const uint32_t exon_end[2] = {150u, 250u};
    static const uint16_t feature_chrom[1] = {0u};
    static const uint32_t feature_start[1] = {1000u};
    static const uint32_t feature_end[1] = {1020u};
    static const uint8_t feature_kind[1] = {
        DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION
    };
    /* The kernel shifts raw local BND POS by one. The resulting local points
     * are transcript_end + {5000,5001} and feature_end + {5000,5001}. */
    static const uint16_t chrom[4] = {0u, 0u, 0u, 0u};
    static const uint32_t pos[4] = {5249u, 5250u, 6019u, 6020u};
    static const uint32_t end[4] = {5249u, 5250u, 6019u, 6020u};
    static const uint16_t mate_chrom[4] = {0u, 0u, 0u, 0u};
    static const uint32_t mate_pos[4] = {124u, 124u, 1010u, 1010u};
    static const uint8_t kind[4] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV,
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[4] = {
        DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND,
        DUCKVEP_SV_BREAKEND, DUCKVEP_SV_BREAKEND
    };
    static const uint8_t copy_change[4] = {
        DUCKVEP_COPY_CHANGE_UNKNOWN, DUCKVEP_COPY_CHANGE_UNKNOWN,
        DUCKVEP_COPY_CHANGE_UNKNOWN, DUCKVEP_COPY_CHANGE_UNKNOWN
    };
    static const uint32_t transcript_variant[2] = {0u, 1u};
    static const uint32_t transcript_index[2] = {0u, 0u};
    static const uint32_t feature_variant[2] = {2u, 3u};
    static const uint32_t feature_index[2] = {0u, 0u};
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_interval_feature_model_t features;
    duckvep_variant_batch_t variants;
    duckvep_candidate_pairs_t transcript_pairs;
    duckvep_interval_feature_pairs_t feature_pairs;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t results;
    duckvep_options_init_t init;
    duckvep_error_t error;

    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&features, 0, sizeof features);
    memset(&variants, 0, sizeof variants);
    memset(&transcript_pairs, 0, sizeof transcript_pairs);
    memset(&feature_pairs, 0, sizeof feature_pairs);
    memset(&init, 0, sizeof init);
    memset(&error, 0, sizeof error);
    transcripts.chrom_id = transcript_chrom;
    transcripts.start1 = transcript_start;
    transcripts.end1 = transcript_end;
    transcripts.strand = transcript_strand;
    transcripts.flags = transcript_flags;
    transcripts.exon_offset = exon_offset;
    transcripts.exon_count = exon_count;
    transcripts.cds_start1 = cds_start;
    transcripts.cds_end1 = cds_end;
    transcripts.transcript_count = 1u;
    exons.start1 = exon_start;
    exons.end1 = exon_end;
    exons.exon_count = 2u;
    features.chrom_id = feature_chrom;
    features.start1 = feature_start;
    features.end1 = feature_end;
    features.kind = feature_kind;
    features.feature_count = 1u;
    variants.chrom_id = chrom;
    variants.pos1 = pos;
    variants.end1 = end;
    variants.mate_chrom_id = mate_chrom;
    variants.mate_pos1 = mate_pos;
    variants.variant_kind = kind;
    variants.sv_type = sv_type;
    variants.copy_change = copy_change;
    variants.count = 4u;
    transcript_pairs.variant_idx = transcript_variant;
    transcript_pairs.tx_idx = transcript_index;
    transcript_pairs.count = 2u;
    feature_pairs.variant_idx = feature_variant;
    feature_pairs.feature_idx = feature_index;
    feature_pairs.count = 2u;
    init.distances_are_explicit = 1u;
    init.upstream_dist = 10000u;
    init.downstream_dist = 10000u;
    init.halo = 10000u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(
        &transcripts, &exons, NULL, &features, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));

    duckvep_result_builder_init(&results, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_pairs(model, &variants, &transcript_pairs,
                                     options, workspace, &results, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
              DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),
              rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_REGION_DOWNSTREAM, rows[0].region_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
              DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),
              rows[1].consequence_mask);
    ASSERT_EQ(DUCKVEP_REGION_DOWNSTREAM, rows[1].region_mask);

    duckvep_result_builder_reset(&results);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_interval_feature_pairs(
        model, &variants, &feature_pairs, options, workspace, &results,
        &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
              DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
              rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
              rows[1].consequence_mask);

    /* Narrowing the caller-managed directional window to zero must not remove
     * an endpoint that StructuralVariationOverlap admitted through its fixed
     * 5000-base rule. The local overlap allele now has no directional
     * predicate and therefore contributes VEP's default intergenic term; the
     * mate overlap allele independently contributes feature_truncation. */
    duckvep_options_close(options);
    options = NULL;
    init.upstream_dist = 0u;
    init.downstream_dist = 0u;
    init.halo = 0u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &options, &error));
    duckvep_result_builder_reset(&results);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_pairs(model, &variants, &transcript_pairs,
                                     options, workspace, &results, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
              DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
              rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
              rows[1].consequence_mask);

    /* The interval-feature endpoint rule is the same fixed constant and does
     * not consult the transcript directional option. */
    duckvep_result_builder_reset(&results);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_interval_feature_pairs(
        model, &variants, &feature_pairs, options, workspace, &results,
        &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
              DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
              rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION),
              rows[1].consequence_mask);

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* TranscriptStructuralVariation clamps an insertion at an internal exon
 * entrance to the exonic flank for within_cdna, while the generic exon
 * predicate retains the reversed P+1,P interval. The exact VEP result is thus
 * feature elongation plus non-coding transcript, not exon or intron placement. */
TEST annotate_sv_insertion_at_noncoding_exon_entrance(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {250u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t cds_s[1] = {0u};
    static const uint32_t cds_e[1] = {0u};
    static const uint32_t es[2] = {100u, 200u};
    static const uint32_t ee[2] = {150u, 250u};
    static const uint32_t ecs[2] = {1u, 52u};
    static const uint32_t ece[2] = {51u, 102u};
    static const int8_t phase[2] = {-1, -1};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vstart[1] = {199u};
    static const uint32_t vend[1] = {199u};
    static const uint8_t vkind[1] = {DUCKVEP_KIND_SV};
    static const uint8_t sv_type[1] = {DUCKVEP_SV_INSERTION};
    static const uint8_t copy_change[1] = {DUCKVEP_COPY_CHANGE_UNKNOWN};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = cds_s; tx.cds_end1 = cds_e;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = ecs; ex.cdna_end1 = ece;
    ex.phase = phase; ex.end_phase = phase; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.sv_type = sv_type; v.copy_change = copy_change;
    v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FEATURE_ELONGATION) |
              DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
              row.consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Complete neutral/undirected structural spans use VEP's transcript-level
 * fallback predicates. They must not leak the point-era coding_sequence or
 * non_coding_transcript_exon fallback merely because the span crosses CDS/exon. */
TEST annotate_complete_neutral_sv_uses_transcript_fallbacks(void) {
    static const uint16_t tchrom[3] = {0u, 1u, 2u};
    static const uint32_t tstart[3] = {100u, 100u, 100u};
    static const uint32_t tend[3] = {300u, 300u, 300u};
    static const int8_t strand[3] = {1, 1, 1};
    static const uint64_t flags[3] = {
        DUCKVEP_TX_HAS_TRANSLATION | DUCKVEP_TX_BIOTYPE_PROTEIN_CODING,
        0u,
        DUCKVEP_TX_HAS_TRANSLATION
    };
    static const uint32_t exoff[3] = {0u, 1u, 2u};
    static const uint16_t excnt[3] = {1u, 1u, 1u};
    static const uint32_t cds_s[3] = {100u, 0u, 100u};
    static const uint32_t cds_e[3] = {300u, 0u, 300u};
    static const uint32_t es[3] = {100u, 100u, 100u};
    static const uint32_t ee[3] = {300u, 300u, 300u};
    static const uint16_t vchrom[3] = {0u, 1u, 2u};
    static const uint32_t vstart[3] = {50u, 50u, 50u};
    static const uint32_t vend[3] = {350u, 350u, 350u};
    static const uint8_t vkind[3] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[3] = {
        DUCKVEP_SV_CNV, DUCKVEP_SV_CNV, DUCKVEP_SV_CNV
    };
    static const uint8_t copy_change[3] = {
        DUCKVEP_COPY_CHANGE_NEUTRAL, DUCKVEP_COPY_CHANGE_NEUTRAL,
        DUCKVEP_COPY_CHANGE_NEUTRAL
    };
    static const uint64_t expected[3] = {
        DUCKVEP_SO(DUCKVEP_SO_CODING_TRANSCRIPT),
        DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
        DUCKVEP_SO(DUCKVEP_SO_INTERGENIC)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[3];
    duckvep_result_builder_t rb;
    duckvep_error_t err;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 3u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 3u;
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.sv_type = sv_type; v.copy_change = copy_change;
    v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ((uint32_t)i, rows[i].variant_idx);
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_IMPACT_MODIFIER, rows[i].impact);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Regression: the emitted row->region_mask is reconciled to AGREE with emission.
 * duckvep_region_mask (the coarse pre-classifier) sets REGION_SPLICE at EVERY exon
 * boundary including the transcript's outer 5'/3' ends, which are NOT real splice
 * sites — so it over-calls. annotate_pair must emit the authoritative splice fact
 * (splice_classify), never the coarse second opinion. This pins: (a) a position in
 * the coarse splice reach of the OUTER 5' end where region_mask must have NO SPLICE
 * bit and emit only 5_prime_utr; (b) a real donor site where SPLICE is retained. */
TEST annotate_region_mask_truthful_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2000u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1100u};
    static const uint32_t tcds_e[1]  = {1900u};
    static const uint32_t estart[2]  = {1000u, 1700u};
    static const uint32_t eend[2]    = {1300u, 2000u};

    /* 1002 = 3rd base of the FIRST exon: inside the coarse exonic splice reach of the
     * transcript's outer 5' end (1000), but no intron there -> not a splice site.
     * 1301 = first intronic base of intron [1301,1699] -> a real essential donor. */
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1002u, 1301u};
    static const uint8_t  vkind[2]  = {0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    /* Document the coarse over-call we are suppressing: the pre-classifier DOES flag
     * SPLICE at 1002 (outer-end reach). The emitted row below must not. */
    ASSERT((duckvep_region_mask(&tx, &exons, 0, 1002u, 3u, 8u) &
            (uint32_t)DUCKVEP_REGION_SPLICE) != 0u);

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));

    /* (a) outer-end over-call suppressed: pure 5'UTR, region_mask has NO SPLICE bit. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR), rows[0].consequence_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UTR, rows[0].region_mask);
    ASSERT_EQ(0u, rows[0].region_mask & (uint32_t)DUCKVEP_REGION_SPLICE);

    /* (b) real essential donor: emitted as splice_donor ALONE because the dinucleotide is a
     * splice site, not within_intron. The region_mask still records intron placement plus
     * splice proximity; it is a structural diagnostic, not the emitted consequence set. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR), rows[1].consequence_mask);
    ASSERT((rows[1].region_mask & (uint32_t)DUCKVEP_REGION_SPLICE) != 0u);
    ASSERT((rows[1].region_mask & (uint32_t)DUCKVEP_REGION_INTRON) != 0u);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Reverse-strand SO-mapping anchor: same coding transcript geometry as above but
 * on the '-' strand. Strand must flip upstream<->downstream and 5'<->3' UTR
 * relative to genomic coordinates; CDS and intron stay strand-independent. */
TEST annotate_reverse_strand_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2000u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1100u};
    static const uint32_t tcds_e[1]  = {1900u};
    static const uint32_t estart[2]  = {1700u, 1000u};
    static const uint32_t eend[2]    = {2000u, 1300u};

    static const uint16_t vchrom[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vpos[6]   = {900u, 1050u, 1200u, 1500u, 1950u, 2100u};
    static const uint8_t  vkind[6]  = {0u, 0u, 0u, 0u, 0u, 0u};

    static const uint64_t exp_mask[6] = {
        DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),    /* 900  genomic-upstream -> downstream on '-' */
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),        /* 1050 before CDS start -> 3' on '-'         */
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),    /* 1200 cds                                    */
        DUCKVEP_SO(DUCKVEP_SO_INTRON),             /* 1500 intron                                 */
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),        /* 1950 after CDS end -> 5' on '-'             */
        DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE)       /* 2100 genomic-downstream -> upstream on '-'  */
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 6u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(6u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 6u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Codon-bucket fusion: a single-exon coding transcript on '+' strand with a CDS
 * sequence pool. annotate_tile must refine the generic coding_sequence_variant
 * into the specific codon consequence via the projection + coding-SNV kernels.
 * CDS = ATG AAA TTT TGG CCC (M K F W P). Expected, hand-computed from the genetic
 * code (the oracle): a synonymous, a missense, and a stop_gained SNV, with the
 * right cds_pos / protein_pos / aa. */
TEST annotate_codon_snv_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};

    /* genomic == transcript (+ strand); cds_bytes[0] is CDS position 1. */
    static const uint8_t cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'T','T','T', 'T','G','G', 'C','C','C'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* 3 SNVs: cds6 A>G (AAA->AAG, K->K syn), cds8 T>A (TTT->TAT, F->Y mis),
     * cds11 G>A (TGG->TAG, W->* stop_gained). */
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1005u, 1007u, 1010u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u}; /* SNV */
    static const uint8_t  abytes[6] = {'A','G', 'T','A', 'G','A'};
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)
    };
    static const int32_t exp_cds[3]     = {6, 8, 11};
    static const int32_t exp_protein[3] = {2, 3, 4};
    static const char    exp_aa_ref[3]  = {'K', 'F', 'W'};
    static const char    exp_aa_alt[3]  = {'K', 'Y', '*'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));

    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
        ASSERT_EQ(DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Reverse-strand codon fusion: single-exon '-' transcript over genomic
 * [1000,1014]; the transcript CDS is the reverse-complement of the genome, so
 * cdna = 1015 - genomic_pos. Transcript CDS = ATG AAA TTT TGG CCC (M K F W P), so
 * the GENOMIC bases are the complement of the transcript bases. Same three edits
 * as the '+' scene, expressed in genomic orientation; annotate_tile must reverse-
 * complement them (via coding_snv_from_cds) and land the identical aa/cds/protein.
 * Variants are listed in ascending genomic order (the kernel's sort precondition);
 * each is hand-derived genomic->transcript via complement (audited before run):
 *   1004 C>T genomic -> transcript G>A at cds11 -> TGG->TAG = W->* stop_gained
 *   1007 A>T genomic -> transcript T>A at cds8  -> TTT->TAT = F->Y missense
 *   1009 T>C genomic -> transcript A>G at cds6  -> AAA->AAG = K->K synonymous
 */
TEST annotate_codon_reverse_strand_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};  /* cdna 1 at the high genomic end (1014) */
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};

    /* TRANSCRIPT-oriented CDS (what the kernel translates). */
    static const uint8_t cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'T','T','T', 'T','G','G', 'C','C','C'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* GENOMIC ref/alt (complement of the transcript bases). */
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1004u, 1007u, 1009u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[6] = {'C','T', 'A','T', 'T','C'}; /* sg: C>T, mis: A>T, syn: T>C */
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)
    };
    static const int32_t exp_cds[3]     = {11, 8, 6};
    static const int32_t exp_protein[3] = {4, 3, 2};
    static const char    exp_aa_ref[3]  = {'W', 'F', 'K'};
    static const char    exp_aa_alt[3]  = {'*', 'Y', 'K'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Ref-mismatch fallback (BLACK-BOX behavior, not a branch proof): a coding SNV
 * whose stated genomic REF does not match the reference base must NOT be
 * force-classified into a codon term — the row stays the generic
 * coding_sequence_variant with no codon coordinates (cds_pos == -1). (The
 * REF_MISMATCH status itself is asserted directly by the coding_snv_from_cds unit
 * tests; here we only pin annotate_tile's observable fallback.) */
TEST annotate_codon_ref_mismatch_falls_back(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* pos 1003 = cds 4 (ref 'A'); we LIE and say REF is 'C'. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1003u};
    static const uint8_t  vkind[1]  = {0u};
    static const uint8_t  abytes[2] = {'C', 'G'}; /* wrong ref C, alt G */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {1u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cds_pos);      /* not refined */
    ASSERT_EQ(0, (int)rows[0].aa_ref);
    ASSERT_EQ(DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              rows[0].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* stop_retained anchor: a synonymous change AT the stop codon (stop -> stop) is
 * stop_retained_variant, not synonymous_variant. CDS = ATG AAA TAA (M K *), with a
 * 3' UTR after it (exon to 1014, CDS to 1008) so the terminal stop is not within
 * splice reach of the exon end. SNV at cds8 A>G: TAA -> TGA (still a stop). */
TEST annotate_codon_stop_retained_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u}; /* CDS ends before the exon (3' UTR follows) */
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1007u}; /* cds 8, 2nd base of the stop codon */
    static const uint8_t  vkind[1]  = {0u};
    static const uint8_t  abytes[2] = {'A', 'G'}; /* TAA -> TGA */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {1u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED), rows[0].consequence_mask);
    ASSERT_EQ(8, rows[0].cds_pos);
    ASSERT_EQ(3, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* start_lost anchor: VEP's SNV start_lost is not just "missense at protein 1".
 * It requires overlap with the annotated start codon plus translation_start==1
 * and compares the complete reference and alternate peptide alleles. Ensembl's
 * initial_met edit therefore changes the reference used by the predicate without
 * rewriting the stored CDS. */
TEST annotate_codon_start_lost_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {900u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {900u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {115u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'T','A','C', 'A','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint32_t peptide_edit_offset[2] = {0u, 1u};
    static const uint32_t peptide_edit_position1[1] = {1u};
    static const uint8_t  peptide_edit_alt[1] = {(uint8_t)'M'};
    static const uint64_t nf_flags[1] = {(uint64_t)DUCKVEP_TX_CDS_START_NF};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1000u, 1002u, 1002u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[6] = {'T','A', 'C','A', 'C','T'};
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)
    };
    static const char exp_alt[3] = {'N', '*', 'Y'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(1, rows[i].protein_pos);
        ASSERT_EQ((uint8_t)'Y', rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_alt[i], rows[i].aa_alt);
    }

    duckvep_workspace_close(ws); ws = NULL;
    duckvep_options_close(opts); opts = NULL;
    duckvep_model_close(model); model = NULL;

    /* Ensembl initial_met changes the reference peptide used by TVA::peptide,
     * without rewriting the CDS. The raw TAC->TAT change is no longer
     * synonymous Y->Y: the annotated reference residue is M, so VEP reports
     * start_lost and suppresses the ordinary missense label. */
    seq.peptide_edit_offset = peptide_edit_offset;
    seq.peptide_edit_position1 = peptide_edit_position1;
    seq.peptide_edit_alt = peptide_edit_alt;
    seq.peptide_edit_count = 1u;
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
              DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              rows[1].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), rows[2].consequence_mask);
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ((uint8_t)'M', rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_alt[i], rows[i].aa_alt);
    }
    duckvep_workspace_close(ws); ws = NULL;
    duckvep_options_close(opts); opts = NULL;
    duckvep_model_close(model); model = NULL;

    /* VEP _overlaps_start_codon returns false for cds_start_NF transcripts, so the
     * same first-codon SNV falls back to ordinary missense rather than start_lost. */
    tx.flags = nf_flags;
    v.count = 1u;
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Multi-exon codon stitching: a 2-exon '+' transcript (exon0 [1000,1010] cdna
 * 1-11, exon1 [2000,2010] cdna 12-22), CDS cdna 1..21 = ATG AAA AAA CGT AAA AAA AAA
 * (M K K R K K K). codon4 (CGT) STRADDLES the junction: cds10,11 in exon0
 * (genomic 1009,1010), cds12 in exon1 (genomic 2000). Three uploaded alleles:
 *   - genomic 1004 (cds5, mid-exon0, no splice): AAA->ATA = K->I missense.
 *   - genomic 1009 CGA>CTA spans the exon/intron boundary but differs only at
 *     genomic 1010. VEP maps the whole equal-length feature, cannot form peptide
 *     alleles, and emits coding_sequence_variant plus splice_region_variant.
 *   - genomic 2000 (cds12, the split base, exon1 ACCEPTOR): CGT->CGC = R->R
 *     synonymous AND splice_region (a real internal boundary). Correctly emitting
 *     synonymous here proves the codon was assembled across the junction: a wrong
 *     contiguous-genomic read of 2000..2002 would be TAA (a stop), not Arg, so an
 *     R->R synonymous result is only reachable via the spliced CDS (1009,1010,2000). */
TEST annotate_codon_multi_exon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2010u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {2009u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1010u, 2010u};
    static const uint32_t ecdna_s[2] = {1u, 12u};
    static const uint32_t ecdna_e[2] = {11u, 22u};
    static const int8_t   ephase[2]  = {0, 0};
    /* spliced transcript CDS (21 nt). */
    static const uint8_t cds_bytes[21] = {
        'A','T','G', 'A','A','A', 'A','A','A', 'C','G','T',
        'A','A','A', 'A','A','A', 'A','A','A'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {21u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1004u, 1009u, 2000u};
    static const uint32_t vend[3]   = {1004u, 1011u, 2000u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[10] = {
        'A','T', 'C','G','A', 'C','T','A', 'T','C'
    }; /* 1004 A>T, 1009 CGA>CTA, 2000 T>C */
    static const uint32_t roff[3]   = {0u, 2u, 8u};
    static const uint32_t aoff[3]   = {1u, 5u, 9u};
    static const uint16_t rlen[3]   = {1u, 3u, 1u};
    static const uint16_t alen[3]   = {1u, 3u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION)
    };
    static const int32_t exp_cds[3]     = {5, -1, 12};
    static const int32_t exp_protein[3] = {2, -1, 4};
    static const char    exp_aa_ref[3]  = {'K', '\0', 'R'};
    static const char    exp_aa_alt[3]  = {'I', '\0', 'R'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 21u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
    }
    ASSERT_EQ(DUCKVEP_SEQUENCE_NOT_APPLICABLE, rows[1].sequence_status);
    ASSERT_EQ(0u, rows[1].flags &
                   (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* A 95-base literal REF/ALT feature with exonic endpoints and one internal
 * intron remains an ordinary VEP VariationFeature. The independent-event
 * compatibility path replaces the outer mapped CDS range with the complete
 * ALT, while the raw equal-length predicate suppresses frame/in-frame labels.
 * This pins the composed annotation route, not only its projection helper. */
TEST annotate_long_literal_internal_intron_uses_uploaded_feature_context(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {208u};
    static const int8_t tstrand[1] = {1};
    static const uint32_t texoff[1] = {0u};
    static const uint16_t texcnt[1] = {2u};
    static const uint32_t tcds_s[1] = {100u};
    static const uint32_t tcds_e[1] = {208u};
    static const uint32_t estart[2] = {100u, 200u};
    static const uint32_t eend[2] = {108u, 208u};
    static const uint32_t ecdna_s[2] = {1u, 10u};
    static const uint32_t ecdna_e[2] = {9u, 18u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[18] = {
        'A','T','G', 'A','A','A', 'C','C','C',
        'G','G','G', 'T','T','T', 'T','A','A'
    };
    static const uint64_t cds_offset[1] = {0u};
    static const uint32_t cds_length[1] = {18u};
    static const uint8_t codon_table[1] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1] = {107u};
    static const uint32_t vend[1] = {201u};
    static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint32_t ref_offset[1] = {0u};
    static const uint32_t alt_offset[1] = {95u};
    static const uint16_t allele_length[1] = {95u};
    uint8_t allele_bytes[190];
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    const duckvep_workspace_delta_route_stats_t *stats;
    uint64_t forbidden;

    memset(allele_bytes, (int)'C', 95u);
    memset(allele_bytes + 95u, (int)'A', 95u);
    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = k_zero_flags;
    tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e;
    tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_offset; seq.cds_length = cds_length;
    seq.codon_table = codon_table; seq.transcript_count = 1u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind;
    variants.allele_bytes = allele_bytes;
    variants.allele_bytes_len = sizeof allele_bytes;
    variants.ref_offset = ref_offset; variants.alt_offset = alt_offset;
    variants.ref_length = allele_length; variants.alt_length = allele_length;
    variants.count = 1u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_workspace_delta_route_stats_reset(workspace);
    duckvep_result_builder_init(&builder, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &builder, &error));
    ASSERT_EQ(1u, duckvep_result_builder_count(&builder));
    stats = duckvep_workspace_delta_route_stats(workspace);
    ASSERT(stats != NULL);
    ASSERT_EQ(UINT64_C(1), stats->uploaded_feature_context);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, row.sequence_status);
    forbidden = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
                DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION);
    ASSERT_EQ(UINT64_C(0), row.consequence_mask & forbidden);
    ASSERT(row.consequence_mask != 0u);

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* Same-codon MNV reference: a multi-base substitution contained within one codon uses the
 * same codon-change oracle as SNVs. The direct reference also accepts a narrow two-codon
 * missense case; production MNVs always use the shared CodingContext interpreter. */
TEST annotate_codon_mnv_same_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* 2bp MNV at genomic 1003-1004 (cds 4-5), AA>GG. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1003u};
    static const uint32_t vend[1]   = {1004u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[4] = {'A','A', 'G','G'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t rlen[1]   = {2u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(1u, stats->substitution_context);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_start_lost_route_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1000u};
    static const uint32_t vend[1]   = {1002u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[6] = {'A','T','G', 'G','C','A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {3u};
    static const uint16_t rlen[1]   = {3u};
    static const uint16_t alen[1]   = {3u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), rows[0].consequence_mask);
    ASSERT_EQ(1, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'M', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'A', rows[0].aa_alt);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(1u, stats->substitution_context);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_cursor_mnv_start_lost_route_matches_tile_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1000u};
    static const uint32_t vend[1]   = {1002u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[6] = {'A','T','G', 'G','C','A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {3u};
    static const uint16_t rlen[1]   = {3u};
    static const uint16_t alen[1]   = {3u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t cursor_n = 0u;
    int saw_full = 0;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats); memset(&cursor_stats, 0, sizeof cursor_stats);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(model, &v, opts, ws, &cur, &err));
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        ASSERT(st == DUCKVEP_OK || st == DUCKVEP_ERR_RESULT_FULL);
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            ASSERT(cursor_n < 4u);
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    cursor_stats = *stats;

    ASSERT(saw_full);
    ASSERT_EQ(1u, cursor_n);
    ASSERT(consequence_rows_equal(&tile_rows[0], &cursor_rows[0]));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), cursor_rows[0].consequence_mask);
    ASSERT_EQ(1u, tile_stats.substitution_context);
    ASSERT_EQ(tile_stats.substitution_context, cursor_stats.substitution_context);

    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_len3_and_cross_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* Variant 0 replaces all of codon 2 (AAA>GGG = K->G). Variant 1 spans codons
     * 2/3 (cds 6-7) but touches the terminal codon in this 9nt CDS, so it stays the
     * generic CDS bucket until terminal-boundary edit-set support. */
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1003u, 1005u};
    static const uint32_t vend[2]   = {1005u, 1006u};
    static const uint8_t  vkind[2]  = {(uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[10] = {'A','A','A', 'G','G','G', 'A','T', 'G','C'};
    static const uint32_t roff[2]   = {0u, 6u};
    static const uint32_t aoff[2]   = {3u, 8u};
    static const uint16_t rlen[2]   = {3u, 2u};
    static const uint16_t alen[2]   = {3u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);
    /* Cross-codon MNV (codon 2 synonymous K->K, codon 3 missense F->L): the generalized
     * window classifier resolves the whole window as a coarse missense_variant, where the
     * old two-codon slice bailed to coding_sequence_variant. Both variants now route through
     * the authoritative interpreter (accepted), so there is no fallback. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[1].consequence_mask);
    ASSERT_EQ(-1, rows[1].cdna_pos);
    ASSERT_EQ(-1, rows[1].cds_pos);
    ASSERT_EQ(-1, rows[1].protein_pos);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(2u, stats->substitution_context);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_cross_codon_missense_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'T','T','A', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'T','T','A', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[5] = {0u, 0u, 0u, 0u, 1u};
    static const uint32_t vpos[5]   = {1002u, 1005u, 1005u, 1011u, 2008u};
    static const uint32_t vend[5]   = {1003u, 1006u, 1006u, 1012u, 2009u};
    static const uint8_t  vkind[5]  = {
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t  abytes[20] = {
        'G','A', 'T','T',  /* +: ATG/AAA -> ATT/TAA = start_lost (M->I) & stop_gained (K->*) */
        'A','T', 'G','G',  /* +: AAA/TTA -> AAG/GTA = one syn, one missense => missense */
        'A','T', 'G','C',  /* +: AAA/TTA -> AAG/CTA = both synonymous => synonymous */
        'G','T', 'A','A',  /* +: GGG/TTT -> GGA/ATT = one syn, one missense => missense */
        'A','T', 'C','C'   /* -: genomic AT>CC -> transcript AT>GG => missense */
    };
    static const uint32_t roff[5] = {0u, 4u, 8u, 12u, 16u};
    static const uint32_t aoff[5] = {2u, 6u, 10u, 14u, 18u};
    static const uint16_t rlen[5] = {2u, 2u, 2u, 2u, 2u};
    static const uint16_t alen[5] = {2u, 2u, 2u, 2u, 2u};
    /* The generalized window classifier resolves every two-codon window (including the
     * start&stop composite and both-synonymous), so none fall back to coding_sequence. */
    static const uint64_t exp_mask[5] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE)
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 5u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(5u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 5u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(-1, rows[i].protein_pos);
        ASSERT_EQ((uint8_t)0u, rows[i].aa_ref);
        ASSERT_EQ((uint8_t)0u, rows[i].aa_alt);
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(5u, stats->substitution_context);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_reverse_strand_same_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* Genomic 1004-1005 is transcript CDS 5-4 on the reverse strand. TT>CC genomic
     * reverse-complements to AA>GG in transcript codon 2: AAA -> GGA = K->G. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1004u};
    static const uint32_t vend[1]   = {1005u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[4] = {'T','T', 'C','C'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t rlen[1]   = {2u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* VEP's _before_coding/_after_coding predicates special-case the reversed
 * insertion interval P+1,P at a CDS endpoint. The mapped flank can be coding
 * while the topology point selected for the insertion is in the adjacent
 * intron; within_cdna is still true and the insertion carries the appropriate
 * UTR term. Pin the 5' boundary on both transcript strands. */
TEST annotate_insertion_at_cds_utr_boundary_both_strands_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1028u, 2034u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 2u};
    static const uint16_t texcnt[2] = {2u, 2u};
    static const uint32_t tcds_s[2] = {1020u, 2000u};
    static const uint32_t tcds_e[2] = {1028u, 2008u};
    static const uint32_t estart[4] = {1000u, 1020u, 2026u, 2000u};
    static const uint32_t eend[4] = {1002u, 1028u, 2034u, 2008u};
    static const uint32_t ecdna_s[4] = {1u, 4u, 1u, 10u};
    static const uint32_t ecdna_e[4] = {3u, 12u, 9u, 18u};
    static const int8_t ephase[4] = {-1, 0, -1, 0};

    static const uint16_t vchrom[2] = {0u, 1u};
    static const uint32_t vpos[2] = {1019u, 2008u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t alleles[6] = {'A', 'A', 'G', 'T', 'T', 'G'};
    static const uint32_t ref_offset[2] = {0u, 3u};
    static const uint32_t alt_offset[2] = {1u, 4u};
    static const uint16_t ref_length[2] = {1u, 1u};
    static const uint16_t alt_length[2] = {2u, 2u};
    const uint64_t expected =
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
        DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION);
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s;
    tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 4u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vpos;
    variants.variant_kind = vkind; variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles;
    variants.ref_offset = ref_offset; variants.alt_offset = alt_offset;
    variants.ref_length = ref_length; variants.alt_length = alt_length;
    variants.count = 2u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&builder, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &builder, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&builder));
    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        if (rows[i].consequence_mask != expected) {
            fprintf(stderr,
                    "\n[CDS/UTR insertion boundary] case=%zu mask=%llu "
                    "region=%u sequence_status=%u flags=%u\n",
                    i, (unsigned long long)rows[i].consequence_mask,
                    rows[i].region_mask, rows[i].sequence_status,
                    rows[i].flags);
        }
        ASSERT_EQ(expected, rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE,
                  rows[i].sequence_status);
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* Equal-length features can cover one CDS endpoint and one base outside the
 * transcript. At the 5' end, VEP's mapper Gap prevents a start peptide and the
 * inverted empty UTR still yields generic coding plus 5' UTR. At the 3' end,
 * a one- or two-base prepared terminal codon remains classifiable and adds the
 * incomplete-terminal term. Pin both outer boundaries on both strands. */
TEST annotate_equal_length_outer_transcript_gaps_both_strands_known_scene(void) {
    static const uint16_t tchrom[4] = {0u, 1u, 2u, 3u};
    static const uint32_t tstart[4] = {1000u, 2000u, 3000u, 4000u};
    static const uint32_t tend[4] = {1008u, 2008u, 3007u, 4007u};
    static const int8_t tstrand[4] = {1, -1, 1, -1};
    static const uint64_t tflags[4] = {0u, 0u, 0u, 0u};
    static const uint32_t texoff[4] = {0u, 1u, 2u, 3u};
    static const uint16_t texcnt[4] = {1u, 1u, 1u, 1u};
    static const uint32_t tcds_s[4] = {1000u, 2000u, 3000u, 4000u};
    static const uint32_t tcds_e[4] = {1008u, 2008u, 3007u, 4007u};
    static const uint32_t estart[4] = {1000u, 2000u, 3000u, 4000u};
    static const uint32_t eend[4] = {1008u, 2008u, 3007u, 4007u};
    static const uint32_t ecdna_s[4] = {1u, 1u, 1u, 1u};
    static const uint32_t ecdna_e[4] = {9u, 9u, 8u, 8u};
    static const int8_t ephase[4] = {0, 0, 0, 0};
    static const uint8_t cds_bytes[34] = {
        'A','T','G', 'A','A','A', 'T','A','A',
        'A','T','G', 'A','A','A', 'T','A','A',
        'A','T','G', 'A','A','A', 'G','G',
        'A','T','G', 'A','A','A', 'G','G'
    };
    static const uint64_t cds_offset[4] = {0u, 9u, 18u, 26u};
    static const uint32_t cds_length[4] = {9u, 9u, 8u, 8u};
    static const uint8_t codon_table[4] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint16_t vchrom[4] = {0u, 1u, 2u, 3u};
    static const uint32_t vpos[4] = {999u, 2008u, 3007u, 3999u};
    static const uint32_t vend[4] = {1000u, 2009u, 3008u, 4000u};
    static const uint8_t vkind[4] = {
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t alleles[16] = {
        'C','A', 'G','T',
        'T','C', 'A','G',
        'G','A', 'T','C',
        'A','C', 'G','T'
    };
    static const uint32_t ref_offset[4] = {0u, 4u, 8u, 12u};
    static const uint32_t alt_offset[4] = {2u, 6u, 10u, 14u};
    static const uint16_t allele_length[4] = {2u, 2u, 2u, 2u};
    static const uint64_t expected[4] = {
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s;
    tx.cds_end1 = tcds_e; tx.transcript_count = 4u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 4u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_offset; seq.cds_length = cds_length;
    seq.codon_table = codon_table; seq.transcript_count = 4u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind; variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles;
    variants.ref_offset = ref_offset; variants.alt_offset = alt_offset;
    variants.ref_length = allele_length; variants.alt_length = allele_length;
    variants.count = 4u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&builder, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &builder, &error));
    ASSERT_EQ(4u, duckvep_result_builder_count(&builder));
    for (i = 0u; i < 4u; i++) {
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(i < 2u ? (uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE
                         : (uint8_t)DUCKVEP_SEQUENCE_RESOLVED,
                  rows[i].sequence_status);
        ASSERT_EQ(0u, rows[i].flags &
                      (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

/* VEP keeps an equal-length uploaded feature intact for CDS mapping even when
 * semantic trimming leaves one changed coding base. If the uploaded span also
 * reaches the transcript's 3' UTR, one CDS endpoint is a mapper gap and VEP
 * cannot form peptide alleles: coding_unknown wins. The adjacent CDS-only
 * controls prove this is a mapping boundary, not blanket stop-loss suppression. */
TEST annotate_equal_length_cds_to_utr3_mapping_gap_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1011u, 2011u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2003u};
    static const uint32_t tcds_e[2] = {1008u, 2011u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1011u, 2011u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {12u, 12u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[18] = {
        'A','T','G', 'A','A','A', 'T','A','A',
        'A','T','G', 'A','A','A', 'T','A','A'
    };
    static const uint64_t cds_off[2] = {0u, 9u};
    static const uint32_t cds_len[2] = {9u, 9u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint8_t post_cds[6] = {'A','A','A', 'A','A','A'};
    static const uint64_t pre_cds_off[2] = {0u, 3u};
    static const uint32_t pre_cds_len[2] = {0u, 0u};
    static const uint64_t post_cds_off[2] = {0u, 3u};
    static const uint32_t post_cds_len[2] = {3u, 3u};

    static const uint16_t vchrom[4] = {0u, 0u, 1u, 1u};
    static const uint32_t vpos[4] = {1007u, 1008u, 2002u, 2003u};
    static const uint32_t vend[4] = {1008u, 1009u, 2003u, 2004u};
    static const uint8_t vkind[4] = {
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t abytes[16] = {
        'A','A', 'C','A', /* + CDS-only: TAA -> TCA, stop_lost */
        'A','A', 'T','A', /* + CDS/3'UTR feature: coding_unknown */
        'T','T', 'T','G', /* - 3'UTR/CDS feature: coding_unknown */
        'T','T', 'G','T'  /* - CDS-only: TAA -> TAC, stop_lost */
    };
    static const uint32_t roff[4] = {0u, 4u, 8u, 12u};
    static const uint32_t aoff[4] = {2u, 6u, 10u, 14u};
    static const uint16_t rlen[4] = {2u, 2u, 2u, 2u};
    static const uint16_t alen[4] = {2u, 2u, 2u, 2u};
    static const uint64_t expected[4] = {
        DUCKVEP_SO(DUCKVEP_SO_STOP_LOST),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_STOP_LOST)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.flank_bytes = post_cds; seq.flank_bytes_len = sizeof post_cds;
    seq.pre_cds_offset = pre_cds_off; seq.pre_cds_length = pre_cds_len;
    seq.post_cds_offset = post_cds_off; seq.post_cds_length = post_cds_len;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(4u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 4u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        if (i == 1u || i == 2u) {
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE,
                      rows[i].sequence_status);
            ASSERT_EQ(0u, rows[i].flags &
                          (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
        }
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* VEP maps the first coding piece of a CDS/3'-UTR boundary feature even when
 * the other piece is a mapper gap. A change in the terminal one- or two-base
 * codon retains the incomplete-terminal-codon predicate, while ALT translation
 * can borrow the 3-prime UTR and independently create a stop. The + transcript has no
 * CDS_END_NF flag, matching Ensembl's mitochondrial partial-stop models; the
 * prepared CDS length is the authority used by VEP's predicate. */
TEST annotate_partial_codon_cds_to_utr3_mapping_gap_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1013u, 2013u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {
        0u,
        (uint64_t)DUCKVEP_TX_CDS_END_NF
    };
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2003u};
    static const uint32_t tcds_e[2] = {1010u, 2013u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1013u, 2013u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {14u, 14u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[22] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G',
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G'
    };
    static const uint64_t cds_off[2] = {0u, 11u};
    static const uint32_t cds_len[2] = {11u, 11u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint8_t post_cds[6] = {'A','A','A', 'A','A','A'};
    static const uint64_t pre_cds_off[2] = {0u, 0u};
    static const uint32_t pre_cds_len[2] = {0u, 0u};
    static const uint64_t post_cds_off[2] = {0u, 3u};
    static const uint32_t post_cds_len[2] = {3u, 3u};

    static const uint16_t vchrom[8] = {
        0u, 0u, 0u, 0u, 1u, 1u, 1u, 1u
    };
    static const uint32_t vpos[8] = {
        1007u, 1008u, 1009u, 1010u,
        2002u, 2003u, 2004u, 2004u
    };
    static const uint32_t vend[8] = {
        1009u, 1009u, 1010u, 1011u,
        2003u, 2004u, 2005u, 2006u
    };
    static const uint8_t vkind[8] = {
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t abytes[36] = {
        'C','C','G', 'A','A','A', /* + table 2: PX/Q* after borrowing UTR A */
        'C','G', 'T','A', /* + table 2: PX/P*, retained full peptide plus stop */
        'G','G', 'A','A', /* +: wholly inside the two-base partial codon */
        'G','A', 'T','A', /* +: terminal CDS G changes; trailing A is 3' UTR */
        'T','C', 'T','A', /* -: genomic C>A is transcript G>T */
        'C','C', 'A','A', /* -: wholly inside the two-base partial codon */
        'C','G', 'T','A', /* -: complete codon P->P, then partial X */
        'C','G','G', 'T','T','T' /* -: complete codon P->Q, then partial X */
    };
    static const uint32_t roff[8] = {
        0u, 6u, 10u, 14u, 18u, 22u, 26u, 30u
    };
    static const uint32_t aoff[8] = {
        3u, 8u, 12u, 16u, 20u, 24u, 28u, 33u
    };
    static const uint16_t allele_len[8] = {
        3u, 2u, 2u, 2u, 2u, 2u, 2u, 3u
    };
    static const uint64_t expected[8] = {
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_MISSENSE)
    };
    /* The native scalar slot is absent for a two-residue window. Rich
     * presentation keeps the independently projected protein start/end. */
    static const int32_t expected_protein[8] = {-1, -1, 4, 4, 4, 4, -1, -1};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.flank_bytes = post_cds; seq.flank_bytes_len = sizeof post_cds;
    seq.pre_cds_offset = pre_cds_off; seq.pre_cds_length = pre_cds_len;
    seq.post_cds_offset = post_cds_off; seq.post_cds_length = post_cds_len;
    seq.flanks_complete = 1u; seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = allele_len; v.alt_length = allele_len; v.count = 8u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(8u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 8u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(expected_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_padded_small_variant_delta_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1002u, 1002u, 1005u};
    static const uint32_t vend[3]   = {1005u, 1006u, 1006u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t  abytes[24] = {
        /* Complete feature spans codons 1-2 while only CDS5 changes. VEP sees
         * MK>MR, so start_lost and start_retained_variant are both true. */
        'G','A','A','A', 'G','A','G','A',
        'G','A','A','A','C', 'G','C',       /* GAAAC>GC deletes codon2 AAA */
        'A','C', 'A','G','C','C','C'        /* AC>AGCCC inserts GCC after codon2 */
    };
    static const uint32_t roff[3] = {0u, 8u, 15u};
    static const uint32_t aoff[3] = {4u, 13u, 17u};
    static const uint16_t rlen[3] = {4u, 5u, 2u};
    static const uint16_t alen[3] = {4u, 2u, 5u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION),
        DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)
    };
    static const int32_t exp_protein[3] = {-1, 2, 3};
    static const int32_t exp_cdna[3] = {-1, -1, -1};
    static const int32_t exp_cds[3] = {-1, -1, -1};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cdna[i], rows[i].cdna_pos);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
    }
    ASSERT_EQ((uint8_t)0u, rows[0].aa_ref);
    ASSERT_EQ((uint8_t)0u, rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_equal_length_feature_window_both_strands_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1014u, 2014u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2000u};
    static const uint32_t tcds_e[2] = {1014u, 2014u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[30] = {
        'A','T','G', 'G','A','A', 'C','C','C', 'T','G','G', 'T','A','A',
        'A','T','G', 'G','A','A', 'C','C','C', 'T','G','G', 'T','A','A'
    };
    static const uint64_t cds_off[2] = {0u, 15u};
    static const uint32_t cds_len[2] = {15u, 15u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    /* Each two-base feature differs at one base. Its paired one-base record has
     * the same semantic CDS edit but does not select the adjacent start/stop codon. */
    static const uint16_t vchrom[12] = {
        0u, 0u, 0u, 0u, 0u, 0u, 1u, 1u, 1u, 1u, 1u, 1u
    };
    static const uint32_t vpos[12] = {
        1002u, 1003u, 1011u, 1011u, 1011u, 1011u,
        2002u, 2002u, 2003u, 2003u, 2011u, 2011u
    };
    static const uint32_t vend[12] = {
        1003u, 1003u, 1011u, 1011u, 1012u, 1012u,
        2003u, 2003u, 2003u, 2003u, 2012u, 2011u
    };
    static const uint8_t vkind[12] = {
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t abytes[36] = {
        'G','G', 'G','A', 'G','A', 'G','A', 'G','T',
        'G','T', 'A','T', 'G','T', 'T','T',
        'A','C', 'A','A', 'A','C', 'A','T', 'C','A', 'C','T',
        'C','C', 'T','C', 'C','T'
    };
    static const uint32_t roff[12] = {
        0u, 4u, 6u, 8u, 10u, 14u, 18u, 22u, 26u, 28u, 30u, 34u
    };
    static const uint32_t aoff[12] = {
        2u, 5u, 7u, 9u, 12u, 16u, 20u, 24u, 27u, 29u, 32u, 35u
    };
    static const uint16_t rlen[12] = {
        2u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 1u, 1u, 2u, 1u
    };
    static const uint16_t alen[12] = {
        2u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 1u, 1u, 2u, 1u
    };
    static const uint64_t expected[12] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE)
    };
    static const int32_t expected_protein[12] = {
        -1, 2, 4, 4, -1, -1, -1, -1, 4, 4, -1, 2
    };
    static const uint8_t expected_ref_aa[12] = {
        0u, 'E', 'W', 'W', 0u, 0u, 0u, 0u, 'W', 'W', 0u, 'E'
    };
    static const uint8_t expected_alt_aa[12] = {
        0u, 'K', '*', 'C', 0u, 0u, 0u, 0u, 'C', '*', 0u, 'K'
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[12];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 12u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 12u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(12u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 12u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(expected_protein[i], rows[i].protein_pos);
        ASSERT_EQ(expected_ref_aa[i], rows[i].aa_ref);
        ASSERT_EQ(expected_alt_aa[i], rows[i].aa_alt);
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(6u, stats->substitution_context);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* VEP maps a same-length feature crossing the 5' UTR/CDS boundary to the
 * start codon only. The long cases also change the next coding codon but remain
 * start_retained_variant when ATG survives; they must not gain missense. */
TEST annotate_equal_length_utr5_start_boundary_both_strands_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1011u, 2011u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1003u, 2000u};
    static const uint32_t tcds_e[2] = {1011u, 2008u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1011u, 2011u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {12u, 12u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[18] = {
        'A','T','G', 'G','A','A', 'T','A','A',
        'A','T','G', 'G','A','A', 'T','A','A'
    };
    static const uint64_t cds_off[2] = {0u, 9u};
    static const uint32_t cds_len[2] = {9u, 9u};
    static const uint8_t flank_bytes[] = "ACGAGC";
    static const uint64_t pre_off[2] = {0u, 3u};
    static const uint64_t post_off[2] = {3u, 6u};
    static const uint32_t pre_len[2] = {3u, 3u};
    static const uint32_t post_len[2] = {0u, 0u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[7] = {0u, 0u, 0u, 0u, 1u, 1u, 1u};
    static const uint32_t vpos[7] = {
        1001u, 1002u, 1002u, 1002u, 2004u, 2008u, 2008u
    };
    static const uint32_t vend[7] = {
        1003u, 1003u, 1007u, 1005u, 2009u, 2009u, 2010u
    };
    static const uint8_t vkind[7] = {
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t abytes[52] = {
        'C','G','A', 'G','A','A',
        'G','A', 'A','C',
        'G','A','T','G','G','A', 'T','A','T','G','A','T',
        'G','A','T','G', 'G','T','A','A',
        'T','C','C','A','T','G', 'G','G','C','A','T','A',
        'T','G', 'C','A',
        'T','G','C', 'T','A','A'
    };
    static const uint32_t roff[7] = {0u, 6u, 10u, 22u, 30u, 42u, 46u};
    static const uint32_t aoff[7] = {3u, 8u, 16u, 26u, 36u, 44u, 49u};
    static const uint16_t rlen[7] = {3u, 2u, 6u, 4u, 6u, 2u, 3u};
    static const uint16_t alen[7] = {3u, 2u, 6u, 4u, 6u, 2u, 3u};
    static const uint64_t expected[7] = {
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[7];
    duckvep_consequence_t cursor_rows[7];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t cursor_count = 0u;
    int saw_full = 0;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    seq.flank_bytes = flank_bytes; seq.flank_bytes_len = sizeof flank_bytes - 1u;
    seq.pre_cds_offset = pre_off; seq.pre_cds_length = pre_len;
    seq.post_cds_offset = post_off; seq.post_cds_length = post_len;
    seq.flanks_complete = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 7u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 7u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(7u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 7u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
        ASSERT_EQ(0u, rows[i].flags &
                      (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(7u, stats->boundary_context);

    duckvep_workspace_delta_route_stats_reset(ws);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_cursor_open(model, &v, opts, ws, &cursor, &err));
    while (!duckvep_annotate_cursor_done(cursor)) {
        duckvep_status_t status;
        duckvep_result_builder_init(&rb, chunk, 1u);
        status = duckvep_annotate_cursor_fill(cursor, &rb, &err);
        ASSERT(status == DUCKVEP_OK || status == DUCKVEP_ERR_RESULT_FULL);
        if (status == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        ASSERT(cursor_count + duckvep_result_builder_count(&rb) <= 7u);
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            cursor_rows[cursor_count++] = chunk[i];
        }
    }
    ASSERT(saw_full);
    ASSERT_EQ(7u, cursor_count);
    for (i = 0u; i < 7u; i++) {
        ASSERT(consequence_rows_equal(&rows[i], &cursor_rows[i]));
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(7u, stats->boundary_context);

    duckvep_annotate_cursor_close(cursor);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    /* Withheld UTR is missing evidence, not permission to reinterpret the
     * same uploaded features as clipped CDS-only substitutions. */
    seq.flanks_complete = 0u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 7u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(7u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 7u; i++) {
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK,
                  rows[i].sequence_status);
        ASSERT(rows[i].flags & (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
                  DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[i].consequence_mask);
    }
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* These source-derived VEP-116 witnesses exercise the allocation-free UTR/CDS
 * boundary path on both strands. The first and fourth intentionally carry VEP's
 * observable start_lost + start_retained_variant combination. The fifth pins
 * VEP's transcript-associated default when every consequence predicate is false. */
TEST annotate_length_changing_cds_boundaries_both_strands_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {100u, 300u};
    static const uint32_t tend[2] = {201u, 401u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {120u, 310u};
    static const uint32_t tcds_e[2] = {191u, 381u};
    static const uint32_t estart[2] = {100u, 300u};
    static const uint32_t eend[2] = {201u, 401u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {102u, 102u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[] =
        "ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA"
        "ATGGTACGTACGTACGTACGTACGTACGTACTACGTACGTACGTACGTACGTACGTACGTACGTACTGGTAA";
    static const uint64_t cds_off[2] = {0u, 72u};
    static const uint32_t cds_len[2] = {72u, 72u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint8_t flank_bytes[] =
        "TACGTACGTACGTACGTACGACGTACGTAC"
        "TACGTACGTACGTACGTACGACGTACGTAC";
    static const uint64_t pre_off[2] = {0u, 30u};
    static const uint32_t pre_len[2] = {20u, 20u};
    static const uint64_t post_off[2] = {20u, 50u};
    static const uint32_t post_len[2] = {10u, 10u};

    static const uint16_t vchrom[5] = {0u, 0u, 0u, 1u, 1u};
    static const uint32_t vpos[5] = {117u, 184u, 190u, 309u, 380u};
    static const uint32_t vend[5] = {120u, 191u, 193u, 310u, 383u};
    static const uint8_t vkind[5] = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_DEL
    };
    static const uint8_t allele_bytes[] = {
        'A','C','G','A', 'A',
        'A','A','A','C', 'A',
        'T','T', 'G',
        'A','T','C','G', 'A',
        'A','C','T','G','G','T','A','A',
        'A','A','C','C','G','G','T','T','G','A','C','A','C','T','A','T','T','A',
        'C','T','C','A','T','A','C','C','A','A','T','G','G','G','T','G','C'
    };
    static const uint32_t ref_off[5] = {0u, 18u, 5u, 10u, 13u};
    static const uint32_t alt_off[5] = {4u, 26u, 9u, 12u, 17u};
    static const uint16_t ref_len[5] = {4u, 8u, 4u, 2u, 4u};
    static const uint16_t alt_len[5] = {1u, 35u, 1u, 1u, 1u};
    static const uint64_t expected[5] = {
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_STOP_LOST),
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) |
            DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[5];
    duckvep_consequence_t cursor_rows[5];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    duckvep_event_t first_event;
    duckvep_region_state_t first_region;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t cursor_count = 0u;
    int saw_full = 0;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&variants, 0, sizeof variants);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e;
    tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend; exons.cdna_start1 = ecdna_s;
    exons.cdna_end1 = ecdna_e; exons.phase = ephase;
    exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes - 1u;
    seq.cds_offset = cds_off; seq.cds_length = cds_len;
    seq.codon_table = cds_tab; seq.transcript_count = 2u;
    seq.flank_bytes = flank_bytes; seq.flank_bytes_len = sizeof flank_bytes - 1u;
    seq.pre_cds_offset = pre_off; seq.pre_cds_length = pre_len;
    seq.post_cds_offset = post_off; seq.post_cds_length = post_len;
    seq.flanks_complete = 1u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
    variants.variant_kind = vkind; variants.allele_bytes = allele_bytes;
    variants.allele_bytes_len = sizeof allele_bytes;
    variants.ref_offset = ref_off; variants.alt_offset = alt_off;
    variants.ref_length = ref_len; variants.alt_length = alt_len;
    variants.count = 5u;

    duckvep_event_load(&variants, 0u, &first_event);
    ASSERT_EQ(118u, first_event.feature_start1);
    ASSERT_EQ(120u, first_event.feature_end1);
    first_region = duckvep_region_classify_span(
        &tx, &exons, 0u, first_event.feature_start1,
        first_event.feature_end1, 3u, 8u);
    ASSERT(first_region.overlaps_utr5);
    ASSERT(first_region.overlaps_cds);

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 5u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, opts, ws, &rb, &err));
    ASSERT_EQ(5u, duckvep_result_builder_count(&rb));
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    for (i = 0u; i < 5u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ_FMT((uint8_t)DUCKVEP_SEQUENCE_RESOLVED,
                      rows[i].sequence_status, "%u");
        ASSERT_EQ_FMT(expected[i], rows[i].consequence_mask, "%" PRIu64);
    }
    ASSERT_EQ_FMT(UINT64_C(4), stats->boundary_context, "%" PRIu64);
    duckvep_workspace_delta_route_stats_reset(ws);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(
        model, &variants, opts, ws, &cursor, &err));
    while (!duckvep_annotate_cursor_done(cursor)) {
        duckvep_status_t status;

        duckvep_result_builder_init(&rb, chunk, 1u);
        status = duckvep_annotate_cursor_fill(cursor, &rb, &err);
        ASSERT(status == DUCKVEP_OK || status == DUCKVEP_ERR_RESULT_FULL);
        if (status == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        ASSERT(cursor_count + duckvep_result_builder_count(&rb) <= 5u);
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            cursor_rows[cursor_count++] = chunk[i];
        }
    }
    ASSERT(saw_full);
    ASSERT_EQ(5u, cursor_count);
    for (i = 0u; i < 5u; i++) {
        ASSERT(consequence_rows_equal(&rows[i], &cursor_rows[i]));
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(4u, stats->boundary_context);

    duckvep_annotate_cursor_close(cursor);
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    cursor = NULL; ws = NULL; model = NULL;

    /* Bytes alone are insufficient: only a receipt-backed complete-flank model
     * may resolve this VEP string predicate. */
    seq.flanks_complete = 0u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 5u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, opts, ws, &rb, &err));
    ASSERT_EQ(5u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 5u; i++) {
        if (i == 1u) {
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED,
                      rows[i].sequence_status);
            ASSERT_EQ(0u, rows[i].flags &
                          (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
            ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTERGENIC),
                      rows[i].consequence_mask);
        } else {
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_FLANK,
                      rows[i].sequence_status);
            ASSERT(rows[i].flags &
                   (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
        }
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* VEP-116 mapper-endpoint witnesses discovered by seed 71. A length-changing
 * feature with one transcript endpoint in a mapper Gap cannot form ordinary
 * codon/peptide alleles. VEP emits coding_unknown for the mapped CDS portion,
 * except that its independent genomic terminal-codon test turns a missing
 * sequence-edit result into stop_retained_variant. Mirror both states on the
 * negative strand to pin transcript-oriented endpoint and stop selection. */
TEST length_changing_mapper_endpoint_gap_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2] = {250u, 250u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 2u};
    static const uint16_t texcnt[2] = {2u, 2u};
    static const uint32_t tcds_s[2] = {120u, 120u};
    static const uint32_t tcds_e[2] = {240u, 240u};
    static const uint32_t estart[4] = {100u, 200u, 200u, 100u};
    static const uint32_t eend[4] = {150u, 250u, 250u, 150u};
    static const uint32_t ecdna_s[4] = {1u, 52u, 1u, 52u};
    static const uint32_t ecdna_e[4] = {51u, 102u, 51u, 102u};
    static const int8_t ephase[4] = {0, 0, 0, 0};
    static const uint8_t coding_gap_alleles[] =
        "CGTACGTACGTACGTACGATGGTACGTACGTACGTACGTACGTACGTACG"
        "C";
    static const uint8_t stop_gap_alleles[] =
        "CGTACGTACGTACGTACGTACGTACGTACGTACTGGTAAACGTACGTACG"
        "C";
    static const uint32_t raw_pos1[4] = {102u, 202u, 198u, 98u};
    static const uint32_t feature_start1[4] = {103u, 203u, 199u, 99u};
    static const uint32_t feature_end1[4] = {151u, 251u, 247u, 147u};
    static const size_t tx_idx[4] = {0u, 0u, 1u, 1u};
    static const uint8_t expect_stop_retained[4] = {0u, 1u, 0u, 1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    size_t i;

    ASSERT_EQ(51u, sizeof coding_gap_alleles - 1u);
    ASSERT_EQ(51u, sizeof stop_gap_alleles - 1u);
    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = tflags; tx.exon_offset = texoff;
    tx.exon_count = texcnt; tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e;
    tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend; exons.cdna_start1 = ecdna_s;
    exons.cdna_end1 = ecdna_e; exons.phase = ephase;
    exons.end_phase = ephase; exons.exon_count = 4u;

    for (i = 0u; i < 4u; i++) {
        uint16_t vchrom[1] = {tchrom[tx_idx[i]]};
        uint32_t vpos[1] = {raw_pos1[i]};
        uint32_t vend[1] = {raw_pos1[i] + 49u};
        uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_DEL};
        uint32_t roff[1] = {0u};
        uint32_t aoff[1] = {50u};
        uint16_t rlen[1] = {50u};
        uint16_t alen[1] = {1u};
        duckvep_variant_batch_t variants;
        duckvep_event_t event;
        duckvep_sequence_delta_route_t route;
        duckvep_sequence_delta_t delta;
        const uint8_t *alleles = (i % 2u) == 0u
            ? coding_gap_alleles : stop_gap_alleles;

        memset(&variants, 0, sizeof variants);
        variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vend;
        variants.variant_kind = vkind; variants.allele_bytes = alleles;
        variants.allele_bytes_len = 51u;
        variants.ref_offset = roff; variants.alt_offset = aoff;
        variants.ref_length = rlen; variants.alt_length = alen;
        variants.count = 1u;
        duckvep_event_load(&variants, 0u, &event);
        ASSERT_EQ(feature_start1[i], event.feature_start1);
        ASSERT_EQ(feature_end1[i], event.feature_end1);

        duckvep_sequence_delta_fill_for_annotation_trace(
            DUCKVEP_KIND_DEL, &tx, &exons, NULL, &variants, 0u,
            tx_idx[i], event.feature_start1, tstrand[tx_idx[i]], NULL,
            &event, UINT32_MAX, UINT32_MAX, &route, &delta);
        ASSERT_EQ(DUCKVEP_DELTA_ROUTE_BOUNDARY_CONTEXT, route);
        ASSERT(delta.valid);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED,
                  delta.sequence_status);
        ASSERT_EQ(expect_stop_retained[i], delta.stop_retained);
        ASSERT_EQ((uint8_t)!expect_stop_retained[i], delta.coding_unknown);
        ASSERT(!delta.stop_lost && !delta.frameshift &&
               !delta.inframe_deletion && !delta.inframe_insertion);
    }
    PASS();
}

TEST annotate_codon_indel_frameshift_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[7] = {0u, 0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vpos[7]   = {1000u, 1002u, 1003u, 1004u, 1006u, 1006u, 1006u};
    static const uint32_t vend[7]   = {1000u, 1005u, 1003u, 1005u, 1006u, 1007u, 1006u};
    static const uint8_t  vkind[7]  = {(uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[23] = {'A', 'A','T', 'G','A','A','A', 'G', 'A', 'A','T', 'A','A', 'A', 'A', 'A','T', 'A','A', 'A', 'T', 'T','A'};
    static const uint32_t roff[7]   = {0u, 3u, 8u, 11u, 14u, 17u, 20u};
    static const uint32_t aoff[7]   = {1u, 7u, 9u, 13u, 15u, 19u, 21u};
    static const uint16_t rlen[7]   = {1u, 4u, 1u, 2u, 1u, 2u, 1u};
    static const uint16_t alen[7]   = {2u, 1u, 2u, 1u, 2u, 1u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 7u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(7u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                  DUCKVEP_SO(DUCKVEP_SO_START_LOST),
              rows[0].consequence_mask);
    ASSERT((rows[0].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT) != 0u);
    ASSERT((rows[0].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_START_LOST) != 0u);
    ASSERT((rows[0].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID) != 0u);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION), rows[1].consequence_mask);
    ASSERT_EQ(0u, rows[1].flags &
        (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT);
    ASSERT_EQ(-1, rows[1].cds_pos);
    ASSERT_EQ(2, rows[1].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), rows[2].consequence_mask);
    ASSERT((rows[2].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT) != 0u);
    ASSERT_EQ(-1, rows[2].cds_pos);
    ASSERT_EQ(2, rows[2].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), rows[3].consequence_mask);
    ASSERT((rows[3].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT) != 0u);
    ASSERT_EQ(-1, rows[3].cds_pos);
    ASSERT_EQ(2, rows[3].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[4].consequence_mask); /* wrong INS REF */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[5].consequence_mask); /* wrong DEL REF */
    /* Final-codon +1 INS with the ATG start intact: the general CodingContext resolves
     * the frameshift the direct body-only restriction rejected at the terminal codon.
     * VEP's coordinate-only endpoint reconstruction also calls the non-stop TTT
     * endpoint altered, so stop_lost coexists. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                  DUCKVEP_SO(DUCKVEP_SO_STOP_LOST),
              rows[6].consequence_mask);
    ASSERT((rows[6].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT) != 0u);
    ASSERT((rows[6].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_STOP_LOST) != 0u);
    ASSERT_EQ(-1, rows[6].cds_pos);
    ASSERT_EQ(3, rows[6].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_terminal_stop_frame_change_uses_transcript_tail(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {1000u};
    static const uint32_t tend[1] = {1011u};
    static const int8_t tstrand[1] = {1};
    static const uint32_t texoff[1] = {0u};
    static const uint16_t texcnt[1] = {1u};
    static const uint32_t tcds_s[1] = {1000u};
    static const uint32_t tcds_e[1] = {1008u};
    static const uint32_t estart[1] = {1000u};
    static const uint32_t eend[1] = {1011u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {12u};
    static const int8_t ephase[1] = {0};
    static const uint8_t cds_bytes[9] = {
        'A','T','G', 'A','A','A', 'T','A','A'
    };
    static const uint8_t post_cds[3] = {'A','C','G'};
    static const uint64_t pre_cds_off[1] = {0u};
    static const uint32_t pre_cds_len[1] = {0u};
    static const uint64_t post_cds_off[1] = {0u};
    static const uint32_t post_cds_len[1] = {3u};
    static const uint64_t cds_off[1] = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t cds_tab[1] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {1006u, 1006u};
    static const uint32_t vend[2] = {1007u, 1006u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t abytes[6] = {
        'T','A', 'T',
        'T', 'T','T'
    };
    static const uint32_t roff[2] = {0u, 3u};
    static const uint32_t aoff[2] = {2u, 4u};
    static const uint16_t rlen[2] = {2u, 1u};
    static const uint16_t alen[2] = {1u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.flank_bytes = post_cds; seq.flank_bytes_len = sizeof post_cds;
    seq.pre_cds_offset = pre_cds_off; seq.pre_cds_length = pre_cds_len;
    seq.post_cds_offset = post_cds_off; seq.post_cds_length = post_cds_len;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED), rows[0].consequence_mask);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[0].sequence_status);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_alt);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_LOST), rows[1].consequence_mask);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[1].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);

    seq.flank_bytes = NULL; seq.flank_bytes_len = 0u;
    seq.pre_cds_offset = NULL; seq.pre_cds_length = NULL;
    seq.post_cds_offset = NULL; seq.post_cds_length = NULL;
    model = NULL; opts = NULL; ws = NULL;
    memset(&err, 0, sizeof err);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    v.count = 1u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);
    ASSERT((rows[0].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0u);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL,
              rows[0].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_delins_frameshift_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[4] = {0u, 0u, 0u, 1u};
    static const uint32_t vpos[4]   = {1003u, 1004u, 1004u, 2009u};
    static const uint32_t vend[4]   = {1005u, 1005u, 1005u, 2010u};
    static const uint8_t  vkind[4]  = {
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL
    };
    static const uint8_t  abytes[18] = {
        'A','A','A', 'A','G',          /* + padded: raw AAA>AG trims to AA>G */
        'A','A', 'G',                  /* +: cds5-6 AA>G, net -1 frameshift */
        'A','A', 'G','G','G','G','G',  /* +: net +3 delins, edit-set backlog */
        'T','T', 'C'                   /* -: genomic TT>C == transcript AA>G */
    };
    static const uint32_t roff[4] = {0u, 5u, 8u, 15u};
    static const uint32_t aoff[4] = {3u, 7u, 10u, 17u};
    static const uint16_t rlen[4] = {3u, 2u, 2u, 2u};
    static const uint16_t alen[4] = {2u, 1u, 5u, 1u};
    static const uint64_t exp_mask[4] = {
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT),
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)
    };
    static const int32_t exp_protein[4] = {2, 2, 2, 2};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[5];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 5u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(4u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 4u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_inframe_insertion_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1005u, 1005u};
    static const uint32_t vend[2]   = {1005u, 1005u};
    static const uint8_t  vkind[2]  = {(uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[10] = {'A', 'A','G','C','C', 'A', 'A','T','A','A'};
    static const uint32_t roff[2]   = {0u, 5u};
    static const uint32_t aoff[2]   = {1u, 6u};
    static const uint16_t rlen[2]   = {1u, 1u};
    static const uint16_t alen[2]   = {4u, 4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(3, rows[0].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
                  DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              rows[1].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_inframe_insertion_reverse_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1008u}; /* transcript CDS 7 on the '-' strand */
    static const uint32_t vend[1]   = {1008u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[5] = {'G', 'G','G','G','C'}; /* ref G, inserted transcript GCC */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
    ASSERT_EQ(3, rows[0].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_protein_altering_insertion_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[3] = {0u, 0u, 1u};
    static const uint32_t vpos[3]   = {1003u, 1003u, 2010u};
    static const uint32_t vend[3]   = {1003u, 1003u, 2010u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t  abytes[15] = {
        'A', 'A','G','C','C',  /* + strand: A + transcript GCC -> protein_altering */
        'A', 'A','C','C','T',  /* + strand: protein-altering junction + TAA stop */
        'T', 'T','G','G','C'   /* - strand: transcript GCC after CDS 4 */
    };
    static const uint32_t roff[3] = {0u, 5u, 10u};
    static const uint32_t aoff[3] = {1u, 6u, 11u};
    static const uint16_t rlen[3] = {1u, 1u, 1u};
    static const uint16_t alen[3] = {4u, 4u, 4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[6];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 6u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[0].consequence_mask);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING) |
                  DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              rows[1].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[2].consequence_mask);
    ASSERT_EQ(2, rows[2].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_indel_reverse_ref_mismatch_falls_back(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1004u};
    static const uint32_t vend[1]   = {1004u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[3] = {'A', 'A','T'}; /* wrong genomic REF at 1004; actual is T */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_unsorted_variant_batch(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {160u, 150u};
    static const uint8_t vkind[2]   = {0u, 0u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_missing_alleles_for_nonpoint_small_variant(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {151u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_DEL};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_accepts_right_anchored_pure_insertion(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t bytes[3]   = {'A', 'T', 'A'}; /* A>TA inserts before POS */
    static const uint32_t roff[1]   = {0u};
    static const uint16_t rlen[1]   = {1u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t alen[1]   = {2u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.ref_length = rlen; v.alt_offset = aoff; v.alt_length = alen;
    v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_kind_allele_shape_mismatch(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t bytes[3]   = {'T', 'A', 'T'};
    static const uint32_t roff[1]   = {0u};
    static const uint16_t rlen[1]   = {2u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t alen[1]   = {1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.ref_length = rlen; v.alt_offset = aoff; v.alt_length = alen;
    v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_allele_slice_outside_pool(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint64_t cds_off[1] = {0u};
    static const uint32_t cds_len[1] = {0u};
    static const uint8_t codon_table[1] = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t alleles[1] = {'A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u}; /* one past the pool */
    static const uint16_t one[1]    = {1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = codon_table;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.ref_offset = roff; v.ref_length = one; v.alt_offset = aoff; v.alt_length = one;
    v.allele_bytes = alleles; v.allele_bytes_len = sizeof alleles; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

/* A full result builder reports DUCKVEP_ERR_RESULT_FULL and never truncates a
 * partial row: two intronic pairs, capacity one. */
TEST annotate_result_full_is_reported(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t   tstrand[1]= {1};
    static const uint32_t tcds[1]   = {0u};
    static const uint32_t teoff[1]  = {0u};
    static const uint16_t tecnt[1]  = {0u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {150u, 160u};
    static const uint8_t  vkind[2]  = {0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[1];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = teoff; tx.exon_count = tecnt;
    tx.cds_start1 = tcds; tx.cds_end1 = tcds; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 1u);
    ASSERT_EQ(DUCKVEP_ERR_RESULT_FULL,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb)); /* exactly the first row, no partial */

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_result_count_past_capacity(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u}, tend[1] = {200u}, zero32[1] = {0u};
    static const uint16_t zero16[1] = {0u};
    static const int8_t tstrand[1] = {1};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1] = {150u};
    static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    rb.count = 2u;
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

/* Distances are uint32 coordinates. A signed narrowing would turn this 4-billion
 * base upstream distance negative and bypass the 3-billion-base directional cut. */
TEST annotate_directional_distance_uses_u32(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {UINT32_C(4000000001)};
    static const uint32_t tend[1] = {UINT32_C(4000000010)};
    static const uint32_t zero32[1] = {0u};
    static const uint16_t zero16[1] = {0u};
    static const int8_t tstrand[1] = {1};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1] = {1u};
    static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&init, 0, sizeof init); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 1u;
    init.upstream_dist = UINT32_C(3000000000);
    init.downstream_dist = 1u;
    init.halo = UINT32_MAX;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, rb.count);

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

TEST options_reject_unknown_compatibility_profile(void) {
    duckvep_options_init_t init;
    duckvep_options_t *options = NULL;
    duckvep_error_t error;

    memset(&init, 0, sizeof init);
    memset(&error, 0, sizeof error);
    init.compatibility_profile = UINT8_MAX;
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_options_open(&init, &options, &error));
    ASSERT_EQ(NULL, options);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, error.status);
    PASS();
}

/* The symmetric sweep halo over-admits when up/downstream distances differ; the
 * directional filter must drop pairs beyond the per-direction window. halo 5000,
 * up = down = 100: variants 200bp away are dropped, 50bp away are kept. */
TEST annotate_directional_distance_filter(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {1000u};
    static const uint32_t tend[1]   = {2000u};
    static const int8_t   tstrand[1]= {1};
    static const uint32_t tcds[1]   = {0u};
    static const uint32_t teoff[1]  = {0u};
    static const uint16_t tecnt[1]  = {0u};
    static const uint16_t vchrom[4] = {0u, 0u, 0u, 0u};
    static const uint32_t vpos[4]   = {800u, 950u, 2050u, 2300u}; /* up200, up50, down50, down300 */
    static const uint8_t  vkind[4]  = {0u, 0u, 0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = teoff; tx.exon_count = tecnt;
    tx.cds_start1 = tcds; tx.cds_end1 = tcds; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 4u;

    memset(&init, 0, sizeof init);
    init.upstream_dist = 100u; init.downstream_dist = 100u; init.halo = 5000u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));

    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ_FMT(1u, rows[0].variant_idx, "%u"); /* up50 kept */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE), rows[0].consequence_mask);
    ASSERT_EQ_FMT(2u, rows[1].variant_idx, "%u"); /* down50 kept */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE), rows[1].consequence_mask);

    duckvep_options_close(opts);
    opts = NULL;
    memset(&init, 0, sizeof init);
    init.distances_are_explicit = 1u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    duckvep_result_builder_reset(&rb);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_honors_splice_region_options(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {300u};
    static const int8_t tstrand[1] = {1};
    static const uint32_t zero32[1] = {0u};
    static const uint16_t two16[1] = {2u};
    static const uint32_t estart[2] = {100u, 250u};
    static const uint32_t eend[2] = {150u, 300u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {146u, 159u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_options_init_t init;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t results;
    duckvep_error_t error;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants); memset(&init, 0, sizeof init);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = k_zero_flags;
    tx.exon_offset = zero32; tx.exon_count = two16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32;
    tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vpos;
    variants.variant_kind = vkind; variants.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &workspace, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    duckvep_result_builder_init(&results, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &results, &error));
    ASSERT_EQ(2u, results.count);
    ASSERT_EQ(0u, rows[0].consequence_mask &
                  DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION));
    ASSERT_EQ(0u, rows[1].consequence_mask &
                  DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION));
    duckvep_options_close(options); options = NULL;

    init.splice_region_exonic = 5u;
    init.splice_region_intronic = 12u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &options, &error));
    duckvep_result_builder_init(&results, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &results, &error));
    ASSERT_EQ(2u, results.count);
    ASSERT(rows[0].consequence_mask & DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION));
    ASSERT(rows[1].consequence_mask & DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION));
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_EXON | DUCKVEP_REGION_SPLICE),
              rows[0].region_mask);
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE),
              rows[1].region_mask);

    duckvep_options_close(options);
    duckvep_workspace_close(workspace);
    duckvep_model_close(model);
    PASS();
}

TEST differing_regions_honor_wider_splice_windows(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {300u};
    static const int8_t tstrand[1] = {1};
    static const uint32_t zero32[1] = {0u};
    static const uint16_t two16[1] = {2u};
    static const uint32_t estart[2] = {100u, 250u};
    static const uint32_t eend[2] = {150u, 300u};
    static const uint8_t ref[2] = {'A', 'A'};
    static const uint8_t alt[2] = {'C', 'A'};
    static const uint8_t alt_short[1] = {'C'};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_splice_state_t splice;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = tstrand; tx.flags = k_zero_flags;
    tx.exon_offset = zero32; tx.exon_count = two16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32;
    tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;

    /* Exon base five and intron base twelve are outside the default cache
     * gate but inside the requested generic splice-region windows. */
    splice = duckvep_splice_classify_differing_regions_with_windows(
        &tx, &exons, 0u, 146u, ref, 2u, alt, 2u, 5u, 12u);
    ASSERT(splice.splice_region);
    splice = duckvep_splice_classify_differing_regions_with_windows(
        &tx, &exons, 0u, 162u, ref, 2u, alt, 2u, 5u, 12u);
    ASSERT(splice.splice_region);

    splice = duckvep_splice_classify_differing_regions(
        &tx, &exons, 0u, 146u, ref, 2u, alt, 2u);
    ASSERT(!splice.splice_region);
    splice = duckvep_splice_classify_differing_regions(
        &tx, &exons, 0u, 162u, ref, 2u, alt, 2u);
    ASSERT(!splice.splice_region);

    splice = duckvep_splice_classify_differing_regions_with_windows(
        &tx, &exons, 0u, 146u, ref, 2u, alt_short, 1u, 5u, 12u);
    ASSERT(splice.splice_region);
    splice = duckvep_splice_classify_differing_regions(
        &tx, &exons, 0u, 146u, ref, 2u, alt_short, 1u);
    ASSERT(!splice.splice_region);
    PASS();
}

/* model_open validates once and rejects malformed models with STABLE where_codes
 * (golden anchors: 13 exon-range, 14 cds-range, 15 unsorted — see duckvep_kernel.c). */
TEST model_open_rejects_invalid_models(void) {
    duckvep_model_t *model = NULL;
    duckvep_error_t err;
    duckvep_exon_model_t exons;

#if SIZE_MAX > UINT32_MAX
    /* tx_idx is uint32 throughout the ABI; reject an oversized borrowed view
     * before touching any of its one-element sentinel columns. */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t z32[1] = {0u};
        static const int8_t st[1] = {1};
        static const uint16_t z16[1] = {0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = z32; tx.end1 = z32; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = z32; tx.exon_count = z16;
        tx.cds_start1 = z32; tx.cds_end1 = z32;
        tx.transcript_count = (size_t)UINT32_MAX + 1u;
        ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE,
                  duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
        ASSERT_EQ(18u, err.where_code);
        ASSERT(model == NULL);
    }
#endif

    /* exon slice out of range */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {5u};
        static const uint32_t cz[1] = {0u};
        static const uint32_t es[2] = {100u, 150u}, ee[2] = {120u, 200u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 1u;
        exons.start1 = es; exons.end1 = ee; exons.exon_count = 2u; /* slice [0,5) exceeds 2 */
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
        ASSERT_EQ(13u, err.where_code);
        ASSERT_EQ(NULL, model);
    }
    /* cds outside the transcript span */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {0u};
        static const uint32_t cs[1] = {50u}, ce[1] = {150u}; /* 50 < start1 100 */
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cs; tx.cds_end1 = ce;
        tx.transcript_count = 1u;
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
        ASSERT_EQ(14u, err.where_code);
    }
    /* transcripts not sorted by (chrom_id, start1) */
    {
        static const uint16_t c[2] = {0u, 0u};
        static const uint32_t s1[2] = {200u, 100u}; /* descending start on same chrom */
        static const uint32_t e1[2] = {300u, 150u};
        static const int8_t st[2] = {1, 1};
        static const uint32_t eoff[2] = {0u, 0u}; static const uint16_t ecnt[2] = {0u, 0u};
        static const uint32_t cz[2] = {0u, 0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 2u;
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
        ASSERT_EQ(15u, err.where_code);
    }
    /* NULL flags with transcript_count>0 is rejected (distilled tx_flags are mandatory;
     * the kernel will read them for biotype terms — pin the invariant now). */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {0u};
        static const uint32_t cz[1] = {0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = NULL; /* the invariant under test */
        tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 1u;
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_model_open(&tx, &exons, NULL, NULL, &model, &err));
        ASSERT_EQ(11u, err.where_code); /* DVW_MODEL_NULL_VIEW */
        ASSERT(model == NULL);
    }
    PASS();
}
