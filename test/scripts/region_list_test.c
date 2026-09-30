/* Exercise the production helper without DuckDB. Only its single allocation is
 * intercepted; HTSlib remains the real name/coordinate parser. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned allocation_calls;
static int fail_allocation;
static void *region_test_malloc(size_t bytes) {
    allocation_calls++;
    return fail_allocation ? NULL : malloc(bytes);
}
#define malloc region_test_malloc
#include "../../src/region_list.c"
#undef malloc

#define CHECK(test) do { if (!(test)) { \
    fprintf(stderr, "region list: line %d: %s\n", __LINE__, #test); return 0; \
} } while (0)

static int name2id(void *header, const char *name) {
    (void)header;
    const char *names[] = {"chr1", "chr2", "chr,part", "chr1:100-200"};
    for (int i = 0; i < 4; i++) if (strcmp(name, names[i]) == 0) return i;
    return strcmp(name, "bad_header") == 0 ? -2 : -1;
}

static int examples(void) {
    char **items = NULL;
    unsigned int count = 0;
    char error[256];
    CHECK(duckhts_region_list_parse(NULL, &items, &count, error, sizeof(error)));
    CHECK(items == NULL && count == 0);
    CHECK(duckhts_region_list_parse("", &items, &count, error, sizeof(error)));
    CHECK(items == NULL && count == 0);
    const char *empty[] = {",", " , , ", "\t", "chr1,", ",chr1", "chr1, ,chr2"};
    for (unsigned i = 0; i < sizeof(empty) / sizeof(*empty); i++) {
        allocation_calls = 0;
        CHECK(!duckhts_region_list_parse(empty[i], &items, &count, error, sizeof(error)));
        CHECK(items == NULL && count == 0 && allocation_calls == 0);
        CHECK(strstr(error, "empty item") != NULL);
    }
    CHECK(duckhts_region_list_parse(" {chr,part}:1-5 ,chr1,chr1 ",
                                    &items, &count, error, sizeof(error)));
    CHECK(count == 3 && strcmp(items[0], "{chr,part}:1-5") == 0);
    CHECK(strcmp(items[1], "chr1") == 0 && strcmp(items[2], "chr1") == 0);
    CHECK(duckhts_region_list_validate(items, count, name2id, NULL, error, sizeof(error)));
    free(items);
    const char *valid[] = {"chr1:1-10", "chr1:100-", "chr1:-10", "absent:1-10",
                           "{chr1:100-200}", ".", "*"};
    const char *invalid[] = {"chr1:10-2", "chr1:abc", "{chr1", "chr1:100-200", "bad_header"};
    for (unsigned i = 0; i < sizeof(valid) / sizeof(*valid); i++) {
        CHECK(duckhts_region_list_parse(valid[i], &items, &count, error, sizeof(error)));
        CHECK(duckhts_region_list_validate(items, count, name2id, NULL, error, sizeof(error)));
        free(items);
    }
    for (unsigned i = 0; i < sizeof(invalid) / sizeof(*invalid); i++) {
        CHECK(duckhts_region_list_parse(invalid[i], &items, &count, error, sizeof(error)));
        CHECK(!duckhts_region_list_validate(items, count, name2id, NULL, error, sizeof(error)));
        CHECK(strstr(error, "invalid item") != NULL);
        free(items);
    }
    return 1;
}

static uint32_t random_u32(uint32_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    return *state;
}

static int properties(void) {
    uint32_t seed = 190;
    for (unsigned trial = 0; trial < 10000; trial++) {
        char input[8192] = {0}, expected[64][64], error[256];
        char **items = NULL;
        unsigned int count = 0;
        unsigned int n = 1 + random_u32(&seed) % 64;
        size_t used = 0;
        for (unsigned int i = 0; i < n; i++) {
            unsigned int pos = 1 + random_u32(&seed) % 100000;
            snprintf(expected[i], sizeof(expected[i]), "chr%u:%u-%u",
                     1 + random_u32(&seed) % 2, pos, pos + 5);
            int bytes = snprintf(input + used, sizeof(input) - used,
                                 "%s \t%s \n", i ? "," : "", expected[i]);
            CHECK(bytes > 0 && (size_t)bytes < sizeof(input) - used);
            used += (size_t)bytes;
        }
        allocation_calls = 0;
        fail_allocation = 1;
        CHECK(!duckhts_region_list_parse(input, &items, &count, error, sizeof(error)));
        CHECK(allocation_calls == 1 && items == NULL && count == 0);
        CHECK(strstr(error, "out of memory") != NULL);
        fail_allocation = 0;
        allocation_calls = 0;
        CHECK(duckhts_region_list_parse(input, &items, &count, error, sizeof(error)));
        CHECK(allocation_calls == 1 && count == n);
        for (unsigned int i = 0; i < n; i++) CHECK(strcmp(items[i], expected[i]) == 0);
        CHECK(duckhts_region_list_validate(items, count, name2id, NULL, error, sizeof(error)));
        free(items);
        CHECK(used + 2 <= sizeof(input));
        input[used++] = ',';
        input[used] = '\0';
        allocation_calls = 0;
        CHECK(!duckhts_region_list_parse(input, &items, &count, error, sizeof(error)));
        CHECK(items == NULL && count == 0 && allocation_calls == 0);
    }
    return 1;
}

static int plan_add(duckhts_interval_plan_t *plan, const char *chrom, int64_t start, int64_t end,
                    size_t item) {
    char error[256];
    return duckhts_interval_plan_add(plan, chrom, start, end, item, error, sizeof(error));
}

static int intervals(void) {
    char error[256];
    duckhts_interval_plan_t plan = {0};
    /* Empty plans are valid. */
    CHECK(duckhts_interval_plan_finish(&plan, error, sizeof(error)));
    CHECK(plan.finished && plan.count == 0);
    duckhts_interval_plan_destroy(&plan);
    duckhts_interval_plan_destroy(&plan);
    /* Sort by (name, start); merge overlap and adjacency; never broaden a gap. */
    CHECK(plan_add(&plan, "b", 10, 20, 1) && plan_add(&plan, "a", 30, 40, 2) &&
          plan_add(&plan, "a", 5, 10, 3) && plan_add(&plan, "a", 10, 12, 4) &&
          plan_add(&plan, "a", 11, 13, 5) && plan_add(&plan, "a", 14, 15, 6) &&
          plan_add(&plan, "b", 10, 20, 7) && plan_add(&plan, "a:1", 0, 1, 8) &&
          plan_add(&plan, "a", 0, INT64_MAX, 9));
    CHECK(duckhts_interval_plan_finish(&plan, error, sizeof(error)));
    CHECK(plan.count == 3 && plan.added == 9);
    CHECK(strcmp(plan.items[0].chrom, "a") == 0 && plan.items[0].start == 0 &&
          plan.items[0].end == INT64_MAX);
    CHECK(strcmp(plan.items[1].chrom, "a:1") == 0 && plan.items[1].start == 0 && plan.items[1].end == 1);
    CHECK(strcmp(plan.items[2].chrom, "b") == 0 && plan.items[2].start == 10 && plan.items[2].end == 20);
    duckhts_interval_plan_destroy(&plan);
    CHECK(plan_add(&plan, "a", 5, 10, 1) && plan_add(&plan, "a", 30, 40, 2) &&
          plan_add(&plan, "a", 10, 12, 3) && plan_add(&plan, "a", 14, 15, 4) &&
          plan_add(&plan, "a", 11, 13, 5));
    CHECK(duckhts_interval_plan_finish(&plan, error, sizeof(error)));
    CHECK(plan.count == 3);
    CHECK(plan.items[0].start == 5 && plan.items[0].end == 13);
    CHECK(plan.items[1].start == 14 && plan.items[1].end == 15);
    CHECK(plan.items[2].start == 30 && plan.items[2].end == 40);
    duckhts_interval_plan_destroy(&plan);
    /* Invalid intervals are explicit and leave the plan unchanged. */
    CHECK(!duckhts_interval_plan_add(&plan, NULL, 0, 1, 3, error, sizeof(error)) &&
          strstr(error, "interval 3") && strstr(error, "NULL chrom"));
    CHECK(!duckhts_interval_plan_add(&plan, "a", -1, 1, 4, error, sizeof(error)) &&
          strstr(error, "negative start"));
    CHECK(!duckhts_interval_plan_add(&plan, "a", 5, 5, 5, error, sizeof(error)) &&
          strstr(error, "end <= start"));
    CHECK(!duckhts_interval_plan_add(&plan, "a", 6, 5, 6, error, sizeof(error)) &&
          strstr(error, "end <= start"));
    CHECK(plan.count == 0);
    /* Both caps are checked before storage, with overflow-safe arithmetic. */
    plan.added = DUCKHTS_INTERVAL_PLAN_MAX_INTERVALS;
    CHECK(!duckhts_interval_plan_add(&plan, "a", 0, 1, 7, error, sizeof(error)) &&
          strstr(error, "cap is 1000000"));
    plan.added = 0;
    plan.payload = DUCKHTS_INTERVAL_PLAN_MAX_BYTES - 1;
    CHECK(!duckhts_interval_plan_add(&plan, "a", 0, 1, 8, error, sizeof(error)) &&
          strstr(error, "128 MiB"));
    plan.payload = UINT64_MAX;
    CHECK(!duckhts_interval_plan_add(&plan, "a", 0, 1, 9, error, sizeof(error)) &&
          strstr(error, "128 MiB"));
    plan.payload = 0;
    CHECK(plan.count == 0);
    duckhts_interval_plan_destroy(&plan);
    /* A finished plan is immutable. */
    CHECK(duckhts_interval_plan_finish(&plan, error, sizeof(error)));
    CHECK(!duckhts_interval_plan_add(&plan, "a", 0, 1, 1, error, sizeof(error)));
    duckhts_interval_plan_destroy(&plan);
    return 1;
}

int main(void) {
    if (!examples() || !properties() || !intervals()) return 1;
    puts("region list: 10000 trials, seed=190; allocation failure/recovery, exact item order and typed interval plans: OK");
    return 0;
}
