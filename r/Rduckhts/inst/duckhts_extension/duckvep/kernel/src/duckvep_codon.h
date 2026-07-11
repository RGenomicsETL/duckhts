/*
 * duckvep_codon.h — codon translation + coding-change classification (INTERNAL).
 *
 * The CDS bucket of the kernel (the ~0.5% of candidates inside coding sequence)
 * lands here: a reference codon and the alt codon (3 ASCII bases each, A/C/G/T,
 * already strand-oriented and frame-aligned by the coding projection) are
 * translated and classified.
 *
 * CODON TABLE is a parameter, not a constant: nuclear transcripts use the
 * standard genetic code (NCBI transl_table=1); mitochondrial transcripts use the
 * vertebrate mitochondrial code (transl_table=2), where TGA->Trp, ATA->Met, and
 * AGA/AGG->Stop. The kernel is string-free — it never sees a chromosome name.
 * Choosing the table per transcript (i.e. detecting the mitochondrion across the
 * chrM / chrMT / MT / NC_012920 naming zoo) is the ADAPTER's job at model-prepare
 * time; it passes the resolved duckvep_codon_table_t in. [[chr-prefix normalization]]
 *
 * Translation uses base order T=0,C=1,A=2,G=3 with codon index b1*16+b2*4+b3,
 * indexing a canonical amino-acid string per table — so the table itself is the
 * oracle (standard: ATG->M, TGG->W, TAA/TAG/TGA->*; mito: TGA->W, ATA->M,
 * AGA/AGG->*).
 */
#ifndef DUCKVEP_CODON_H
#define DUCKVEP_CODON_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* NCBI translation table ids (the two duckvep needs today). */
typedef enum duckvep_codon_table {
    DUCKVEP_CODON_TABLE_STANDARD  = 1, /* nuclear, transl_table=1 */
    DUCKVEP_CODON_TABLE_VERT_MITO = 2  /* vertebrate mitochondrial, transl_table=2 */
} duckvep_codon_table_t;

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
 * 'X' if any base is not A/C/G/T/U (case-insensitive). An unknown table falls
 * back to the standard code. */
char duckvep_translate_codon(const char *codon3, duckvep_codon_table_t table);

/* Classify a reference->alt codon change under `table`. start-codon context
 * (start_lost) is decided by the caller, which knows whether this is the
 * translation start. */
duckvep_codon_result_t duckvep_codon_change(const char *ref3, const char *alt3,
                                            duckvep_codon_table_t table);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_CODON_H */
