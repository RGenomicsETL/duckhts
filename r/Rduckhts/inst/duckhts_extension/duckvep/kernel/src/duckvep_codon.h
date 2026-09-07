/*
 * duckvep_codon.h — codon translation + coding-change classification (INTERNAL).
 *
 * The CDS bucket of the kernel (the ~0.5% of candidates inside coding sequence)
 * lands here: a reference codon and the alt codon (3 ASCII bases each, A/C/G/T,
 * already strand-oriented and frame-aligned by the coding projection) are
 * translated and classified.
 *
 * CODON TABLE is a parameter, not a constant. The model compiler reads the
 * Ensembl sequence-region `codon_table` attribute and passes its NCBI table id;
 * the kernel never guesses from a chromosome name. The supported set mirrors
 * BioPerl 1.7.8, which is the translator used by VEP 116.
 *
 * Translation uses base order T=0,C=1,A=2,G=3 with codon index b1*16+b2*4+b3,
 * indexing a canonical amino-acid string per table — so the table itself is the
 * oracle (standard: ATG->M, TGG->W, TAA/TAG/TGA->*; mito: TGA->W, ATA->M,
 * AGA/AGG->*).
 */
#ifndef DUCKVEP_CODON_H
#define DUCKVEP_CODON_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* NCBI translation-table id. Named constants cover the common human cases;
 * duckvep_codon_table_supported() is authoritative for the complete set. */
typedef enum duckvep_codon_table {
    DUCKVEP_CODON_TABLE_STANDARD  = 1, /* nuclear, transl_table=1 */
    DUCKVEP_CODON_TABLE_VERT_MITO = 2  /* vertebrate mitochondrial, transl_table=2 */
} duckvep_codon_table_t;

int duckvep_codon_table_supported(duckvep_codon_table_t table);

/* BioPerl 1.7.8 start predicate: true when ANY A/C/G/T expansion is a start.
 * Three borrowed A/C/G/T/U/N bytes, case-insensitive; invalid input/table is false.
 * This predicate is distinct from ordinary translation of the same codon. */
int duckvep_codon_is_start(const uint8_t *codon3, duckvep_codon_table_t table);

/* Borrow the immutable 64-amino-acid translation table for a supported NCBI
 * table id. The caller retains no ownership; NULL means unsupported. This is
 * the bulk-translation path after bases and table id have been validated. */
const char *duckvep_codon_table_amino_acids(duckvep_codon_table_t table);

typedef enum duckvep_translation_status {
    DUCKVEP_TRANSLATION_OK = 0,
    DUCKVEP_TRANSLATION_INVALID_ARG,
    DUCKVEP_TRANSLATION_BUFFER_TOO_SMALL,
    DUCKVEP_TRANSLATION_INVALID_BASE
} duckvep_translation_status_t;

typedef struct duckvep_translation {
    size_t length; /* Every complete codon, including residues after internal stops. */
    size_t first_stop_position1; /* Zero means no stop; otherwise also the visible protein length. */
    uint8_t unambiguous; /* All CDS bases, including a trailing partial codon, are A/C/G/T/U. */
} duckvep_translation_t;

typedef enum duckvep_translation_ambiguity {
    DUCKVEP_TRANSLATION_N_UNKNOWN = 0, /* Conservative independent coding predicates. */
    DUCKVEP_TRANSLATION_N_CONSENSUS   /* BioPerl: resolve all A/C/G/T expansions of N. */
} duckvep_translation_ambiguity_t;

/* Translate once into distinct caller storage of at least cds_length/3 + 1
 * bytes. Full peptide output is NUL-terminated; callers select the prefix
 * through first_stop_position1 when they need the stop-truncated protein.
 * A/C/G/T/U and N are case-insensitive. N_UNKNOWN emits X for any N-containing
 * codon; N_CONSENSUS emits the common amino acid of its expansions, B for D/N,
 * Z for E/Q, and X otherwise. The unambiguous input fact remains false even when
 * the amino acid is resolved. Validate
 * every input base, including partial codons and sequence beyond the first stop.
 * No allocation. Result is zero on failure; bytes may be partial on invalid
 * input. Result storage must not overlap either byte span. */
duckvep_translation_status_t duckvep_translate_cds(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    duckvep_translation_ambiguity_t ambiguity,
    uint8_t *peptide, size_t peptide_capacity, duckvep_translation_t *result);

/* Return the first raw-CDS stop as a one-based peptide position, or zero when
 * no complete codon is a stop. N-containing codons translate to X. The return
 * value is false only for an invalid base, unsupported table, or NULL output. */
int duckvep_cds_first_stop_position1(
    const uint8_t          *cds,
    size_t                  cds_length,
    duckvep_codon_table_t   table,
    uint32_t               *position1_out);

/* Coding-change classification bits. */
enum {
    DUCKVEP_CODON_SYNONYMOUS  = 1u << 0,
    DUCKVEP_CODON_MISSENSE    = 1u << 1,
    DUCKVEP_CODON_STOP_GAINED = 1u << 2,
    DUCKVEP_CODON_STOP_LOST   = 1u << 3,
    DUCKVEP_CODON_INVALID     = 1u << 4
};

typedef struct duckvep_codon_result {
    char     aa_ref;  /* 1-letter AA, '*' = stop, 'X' = invalid input */
    char     aa_alt;
    uint32_t change;  /* one DUCKVEP_CODON_* bit */
} duckvep_codon_result_t;

/* Translate a 3-base ASCII codon to a 1-letter amino acid under `table`. Returns
 * 'X' if any base is not A/C/G/T/U (case-insensitive), or if the table id is
 * not supported. */
char duckvep_translate_codon(const char *codon3, duckvep_codon_table_t table);

/* Classify a reference->alt codon change under `table`. start-codon context
 * (start_lost) is decided by the caller, which knows whether this is the
 * translation start. */
duckvep_codon_result_t duckvep_codon_change(const char *ref3, const char *alt3,
                                            duckvep_codon_table_t table);

/* Equivalent classifier for transcript-oriented codons already normalized to
 * uppercase A/C/G/T/N by the validated model/edit path. N remains invalid; the
 * caller must not pass any other byte. This avoids repeating six general ASCII
 * normalization operations in the coding-SNV hot loop. */
duckvep_codon_result_t duckvep_codon_change_prepared(
    const char *ref3,
    const char *alt3,
    duckvep_codon_table_t table);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_CODON_H */
