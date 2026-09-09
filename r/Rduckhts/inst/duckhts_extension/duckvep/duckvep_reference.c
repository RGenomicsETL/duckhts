#include "duckvep_reference.h"
#include "duckvep_model.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

#define DUCKVEP_REFERENCE_READ_AHEAD 65536u

int duckvep_reference_reader_init(duckvep_reference_reader_t *reader,
    const duckvep_owned_model_t *model, char *bases, size_t capacity,
    char *error, size_t error_size) {
    if (!reader || reader->fai || !model || (capacity && !bases)) {
        duckvep_sql_set_error(error, error_size, "DuckVEP: invalid reference reader initialization");
        return 0;
    }
    reader->model = model; reader->bases = bases; reader->capacity = capacity;
    reader->length = 0u;
    if (!model->reference_fasta_path) return 1;
    if (!duckvep_model_reference_identity_matches(model)) {
        duckvep_sql_set_error(error, error_size,
            "DuckVEP: model reference FASTA or index changed after model load");
        return 0;
    }
    reader->fai = fai_load3_format(model->reference_fasta_open_path,
        model->reference_fai_open_path, model->reference_gzi_open_path, 0, FAI_FASTA);
    if (!reader->fai) {
        duckvep_sql_set_error(error, error_size, "DuckVEP: could not open the model reference FASTA/index");
        return 0;
    }
    if (!duckvep_model_reference_identity_matches(model)) {
        fai_destroy(reader->fai); reader->fai = NULL;
        duckvep_sql_set_error(error, error_size,
            "DuckVEP: model reference FASTA or index changed while a worker was opening it");
        return 0;
    }
    return 1;
}

int duckvep_reference_reader_windows(duckvep_reference_reader_t *reader,
    const duckvep_event_t *event, int *available,
    duckvep_hgvs_reference_window_t *shift, duckvep_hgvs_reference_window_t *lookup,
    char *error, size_t error_size) {
    if (available) *available = 0;
    if (shift) memset(shift, 0, sizeof *shift);
    if (lookup) memset(lookup, 0, sizeof *lookup);
    if (!reader || !reader->model || !event || !available || !shift || !lookup) {
        duckvep_sql_set_error(error, error_size, "DuckVEP: invalid reference-window state");
        return 0;
    }
    const duckvep_owned_model_t *model = reader->model;
    if (!model->reference_fasta_path) return 1;
    if (!reader->fai || !reader->bases || !reader->capacity) {
        duckvep_sql_set_error(error, error_size, "DuckVEP: reference workspace is not initialized");
        return 0;
    }
    size_t begin = 0u, end = model->known_seq_region_count;
    while (begin < end) {
        size_t middle = begin + (end - begin) / 2u;
        if (model->known_seq_regions[middle] < event->chrom_id) begin = middle + 1u;
        else end = middle;
    }
    if (begin == model->known_seq_region_count || model->known_seq_regions[begin] != event->chrom_id ||
        !model->sequence_names || !model->sequence_names[begin] || !model->sequence_lengths) {
        duckvep_sql_set_error(error, error_size,
            "DuckVEP: sequence-region name is absent from the reference-enabled model");
        return 0;
    }
    uint32_t shift_start1, shift_end1, fetch_start1, fetch_end1;
    uint32_t sequence_length = model->sequence_lengths[begin];
    if (duckvep_hgvs_genomic_search_interval(event, sequence_length, &shift_start1, &shift_end1) !=
            DUCKVEP_HGVS_OK ||
        duckvep_hgvs_reference_fetch_interval(event, sequence_length, &fetch_start1, &fetch_end1) !=
            DUCKVEP_HGVS_OK || shift_start1 < fetch_start1 || shift_end1 > fetch_end1) {
        duckvep_sql_set_error(error, error_size,
            "DuckVEP: semantic edit has no valid bounded VEP reference interval");
        return 0;
    }
    uint64_t cached_end1 = (uint64_t)reader->start1 + reader->length;
    if (!reader->length || reader->chrom_id != event->chrom_id ||
        fetch_start1 < reader->start1 || fetch_end1 >= cached_end1) {
        if (!duckvep_model_reference_identity_matches(model)) {
            duckvep_sql_set_error(error, error_size,
                "DuckVEP: pinned reference FASTA or index changed after model load");
            return 0;
        }
        const char *name = model->sequence_names[begin];
        uint32_t cache_end1 = sequence_length - fetch_end1 < DUCKVEP_REFERENCE_READ_AHEAD
            ? sequence_length : fetch_end1 + DUCKVEP_REFERENCE_READ_AHEAD;
        hts_pos_t fetched_length = -1;
        size_t required = 0u;
        int status = faidx_fetch_seq64_into(reader->fai, name, fetch_start1 - 1,
            cache_end1 - 1, NULL, 0u, &fetched_length, &required);
        if (status != -1 || errno != ENOSPC) {
            duckvep_sql_set_error(error, error_size, "DuckVEP: invalid reference fetch capacity");
            return 0;
        }
        if (required > reader->capacity && cache_end1 != fetch_end1) {
            cache_end1 = fetch_end1;
            status = faidx_fetch_seq64_into(reader->fai, name, fetch_start1 - 1,
                cache_end1 - 1, NULL, 0u, &fetched_length, &required);
            if (status != -1 || errno != ENOSPC) {
                duckvep_sql_set_error(error, error_size, "DuckVEP: invalid reference fetch capacity");
                return 0;
            }
        }
        if (required > reader->capacity) {
            snprintf(error, error_size,
                "DuckVEP: reference workspace bytes=%zu, required=%zu at %s:%u-%u in %s",
                reader->capacity, required, name, fetch_start1, fetch_end1, model->reference_fasta_path);
            return 0;
        }
        /* Invalidate before overwriting: a short read must never expose a
         * partially overwritten window as a successful cache hit. */
        reader->length = 0u;
        if (faidx_fetch_seq64_into(reader->fai, name, fetch_start1 - 1, cache_end1 - 1,
                reader->bases, reader->capacity, &fetched_length, &required) != 0 ||
            fetched_length < 0 || (uint64_t)fetched_length != (uint64_t)cache_end1 - fetch_start1 + 1u) {
            duckvep_sql_set_error(error, error_size,
                "DuckVEP: reference FASTA fetch did not return the requested interval");
            return 0;
        }
        if (!duckvep_model_reference_identity_matches(model)) {
            duckvep_sql_set_error(error, error_size,
                "DuckVEP: pinned reference FASTA or index changed during fetch");
            return 0;
        }
        reader->length = (size_t)fetched_length;
        reader->start1 = fetch_start1; reader->chrom_id = event->chrom_id;
    }
    *lookup = (duckvep_hgvs_reference_window_t){(const uint8_t *)reader->bases,
        reader->length, reader->start1, reader->chrom_id};
    *shift = (duckvep_hgvs_reference_window_t){(const uint8_t *)reader->bases + (shift_start1 - reader->start1),
        (size_t)((uint64_t)shift_end1 - shift_start1 + 1u), shift_start1, reader->chrom_id};
    *available = 1;
    return 1;
}
