/* Fail each wrapper allocation while producing non-applicable norm results. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

static int allocation_count;
static int fail_allocation;

static void *
fault_malloc(size_t size)
{
    allocation_count++;
    return allocation_count == fail_allocation ? NULL : malloc(size);
}

static void *
fault_calloc(size_t count, size_t size)
{
    allocation_count++;
    return allocation_count == fail_allocation ? NULL : calloc(count, size);
}

#define malloc fault_malloc
#define calloc fault_calloc
#include "../../src/bcftools_norm_udf.c"
#undef malloc
#undef calloc
DUCKDB_EXTENSION_GLOBAL

void
duckhts_reference_uppercase(char *sequence)
{
    for (; *sequence != '\0'; sequence++)
        *sequence = (char)toupper((unsigned char)*sequence);
}

char *
duckhts_reference_fetch(duckhts_reference_entry_t *entry, const char *chrom,
    hts_pos_t start, hts_pos_t end, duckhts_reference_end_policy_t policy)
{
    (void)entry;
    (void)start;
    (void)end;
    (void)policy;
    if (strcmp(chrom, "missing") == 0)
        return NULL;
    char *sequence = malloc(2u);
    if (sequence == NULL)
        exit(2);
    memcpy(sequence, "A", 2u);
    return sequence;
}

typedef struct {
    const char *chrom;
    const char *ref;
    char *alts[2];
    int n_alt;
    const char *status;
} passthrough_case_t;

static int
check_case(const passthrough_case_t *test)
{
    norm_result_t result;
    char *error = NULL;
    fail_allocation = 0;
    allocation_count = 0;
    int ok = normalize_variant_row(NULL, test->chrom, 1, test->ref,
        (char **)test->alts, test->n_alt, 0, 0, 0, 0, &result, &error);
    int allocations = allocation_count;
    if (!ok || error != NULL || result.applicable ||
        strcmp(result.status, test->status) != 0 ||
        !same_nullable_string(result.ref, test->ref) ||
        !compare_alt_arrays(result.alts, result.n_alt, (char **)test->alts, test->n_alt)) {
        fprintf(stderr, "invalid control for %s\n", test->status);
        free_norm_result(&result);
        free(error);
        return 1;
    }
    free_norm_result(&result);
    int failures = 0;
    for (int allocation = 1; allocation <= allocations; allocation++) {
        allocation_count = 0;
        fail_allocation = allocation;
        error = NULL;
        ok = normalize_variant_row(NULL, test->chrom, 1, test->ref,
            (char **)test->alts, test->n_alt, 0, 0, 0, 0, &result, &error);
        if (ok || error == NULL || strcmp(error, "bcftools_norm_row: out of memory") != 0) {
            fprintf(stderr, "%s allocation %d/%d: success=%d error=%s\n",
                test->status, allocation, allocations, ok, error == NULL ? "<none>" : error);
            failures++;
        }
        free_norm_result(&result);
        free(error);
    }
    return failures;
}

int
main(void)
{
    const passthrough_case_t cases[] = {
        {NULL, "A", {"T", "G"}, 2, "NullInput"},
        {"chr1", "A", {NULL, "G"}, 2, "NullInput"},
        {"chr1", "A", {"A[CHR2:2[", "G"}, 2, "Breakend"},
        {"chr1", "A", {"*", "G"}, 2, "SpanningDeletion"},
        {"chr1", "Z", {"T", "G"}, 2, "InvalidRefAllele"},
        {"missing", "A", {"T", "G"}, 2, "MissingContig"},
        {"chr1", "T", {"C", "G"}, 2, "RefMismatch"},
        {"chr1", "A", {"Z", "G"}, 2, "InvalidAltAllele"},
        {"chr1", "A", {"A", "G"}, 2, "DuplicateAllele"},
        {"chr1", "A", {"<DUP>", NULL}, 1, "SymbolicInsufficientMetadata"},
        {"chr1", "A", {"<INS>", NULL}, 1, "UnsupportedSymbolicAllele"},
        {"chr1", "A", {"A[CHR2:2[", "<NON_REF>"}, 2, "Breakend"}
    };
    int failures = 0;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
        failures += check_case(&cases[i]);
    return failures == 0 ? 0 : 1;
}
