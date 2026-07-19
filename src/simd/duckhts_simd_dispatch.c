#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_ATOMICS__) && \
    !defined(__EMSCRIPTEN__) && !defined(DUCKDB_WASM_EXTENSION)
#include <stdatomic.h>
#define DUCKHTS_SIMD_HAVE_C11_ATOMICS 1
#else
#define DUCKHTS_SIMD_HAVE_C11_ATOMICS 0
#endif

#if !DUCKHTS_SIMD_HAVE_C11_ATOMICS && (defined(__GNUC__) || defined(__clang__))
#define DUCKHTS_SIMD_HAVE_GNU_ATOMICS 1
#else
#define DUCKHTS_SIMD_HAVE_GNU_ATOMICS 0
#endif

#include "duckhts_simd_internal.h"

enum {
    DUCKHTS_SIMD_INIT_UNINITIALIZED = 0,
    DUCKHTS_SIMD_INIT_IN_PROGRESS = 1,
    DUCKHTS_SIMD_INIT_DONE = 2,
    DUCKHTS_SIMD_PRIORITY_UNSET = 0x7fffffff
};

static const char *duckhts_simd_dispatch_mode_name = "capability-mask-dispatch-table";

static duckhts_simd_dispatch_table_t duckhts_simd_table_auto;
static duckhts_simd_dispatch_table_t duckhts_simd_table_scalar;
static duckhts_simd_dispatch_table_t duckhts_simd_table_avx2;
static duckhts_simd_dispatch_table_t duckhts_simd_table_avx512;
static duckhts_simd_dispatch_table_t duckhts_simd_table_neon;
static duckhts_simd_dispatch_table_t duckhts_simd_table_wasm_simd128;
static const duckhts_simd_dispatch_table_t *duckhts_simd_resolve_backend_table(const char *backend,
                                                                                  char *err,
                                                                                  size_t err_len);

#if DUCKHTS_SIMD_HAVE_C11_ATOMICS
static _Atomic(const duckhts_simd_dispatch_table_t *) duckhts_simd_current_table;
static _Atomic int duckhts_simd_initialized;
static const duckhts_simd_dispatch_table_t *duckhts_simd_load_table(void) {
    return atomic_load_explicit(&duckhts_simd_current_table, memory_order_acquire);
}
static void duckhts_simd_store_table(const duckhts_simd_dispatch_table_t *table) {
    atomic_store_explicit(&duckhts_simd_current_table, table, memory_order_release);
}
static int duckhts_simd_load_init_state(void) {
    return atomic_load_explicit(&duckhts_simd_initialized, memory_order_acquire);
}
static int duckhts_simd_is_initialized(void) {
    return duckhts_simd_load_init_state() == DUCKHTS_SIMD_INIT_DONE;
}
static int duckhts_simd_try_begin_init(void) {
    int expected = DUCKHTS_SIMD_INIT_UNINITIALIZED;
    return atomic_compare_exchange_strong_explicit(&duckhts_simd_initialized, &expected,
                                                   DUCKHTS_SIMD_INIT_IN_PROGRESS,
                                                   memory_order_acq_rel, memory_order_acquire);
}
static void duckhts_simd_finish_init(void) {
    atomic_store_explicit(&duckhts_simd_initialized, DUCKHTS_SIMD_INIT_DONE, memory_order_release);
}
static void duckhts_simd_wait_initialized(void) {
    while (duckhts_simd_load_init_state() != DUCKHTS_SIMD_INIT_DONE) {
    }
}
#elif DUCKHTS_SIMD_HAVE_GNU_ATOMICS
static const duckhts_simd_dispatch_table_t *duckhts_simd_current_table;
static int duckhts_simd_initialized;
static const duckhts_simd_dispatch_table_t *duckhts_simd_load_table(void) {
    return __atomic_load_n(&duckhts_simd_current_table, __ATOMIC_ACQUIRE);
}
static void duckhts_simd_store_table(const duckhts_simd_dispatch_table_t *table) {
    __atomic_store_n(&duckhts_simd_current_table, table, __ATOMIC_RELEASE);
}
static int duckhts_simd_load_init_state(void) {
    return __atomic_load_n(&duckhts_simd_initialized, __ATOMIC_ACQUIRE);
}
static int duckhts_simd_is_initialized(void) {
    return duckhts_simd_load_init_state() == DUCKHTS_SIMD_INIT_DONE;
}
static int duckhts_simd_try_begin_init(void) {
    int expected = DUCKHTS_SIMD_INIT_UNINITIALIZED;
    return __atomic_compare_exchange_n(&duckhts_simd_initialized, &expected,
                                       DUCKHTS_SIMD_INIT_IN_PROGRESS, 0,
                                       __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE);
}
static void duckhts_simd_finish_init(void) {
    __atomic_store_n(&duckhts_simd_initialized, DUCKHTS_SIMD_INIT_DONE, __ATOMIC_RELEASE);
}
static void duckhts_simd_wait_initialized(void) {
    while (duckhts_simd_load_init_state() != DUCKHTS_SIMD_INIT_DONE) {
    }
}
#else
/* MSVC is not a supported DuckHTS build path.  This fallback exists only for
 * strictly single-threaded C targets without C11 or GNU atomics. */
static const duckhts_simd_dispatch_table_t *duckhts_simd_current_table;
static int duckhts_simd_initialized;
static const duckhts_simd_dispatch_table_t *duckhts_simd_load_table(void) { return duckhts_simd_current_table; }
static void duckhts_simd_store_table(const duckhts_simd_dispatch_table_t *table) { duckhts_simd_current_table = table; }
static int duckhts_simd_is_initialized(void) { return duckhts_simd_initialized == DUCKHTS_SIMD_INIT_DONE; }
static int duckhts_simd_try_begin_init(void) {
    if (duckhts_simd_initialized != DUCKHTS_SIMD_INIT_UNINITIALIZED) return 0;
    duckhts_simd_initialized = DUCKHTS_SIMD_INIT_IN_PROGRESS;
    return 1;
}
static void duckhts_simd_finish_init(void) { duckhts_simd_initialized = DUCKHTS_SIMD_INIT_DONE; }
static void duckhts_simd_wait_initialized(void) {
    while (duckhts_simd_initialized != DUCKHTS_SIMD_INIT_DONE) {
    }
}
#endif

typedef int (*duckhts_simd_status_fn)(void);

typedef struct duckhts_simd_backend_entry {
    const char *name;
    duckhts_simd_cap_t cap;
    int selectable;
    duckhts_simd_status_fn compiled;
    duckhts_simd_status_fn cpu_supported;
    int priority;
    duckhts_simd_dispatch_table_t *table;
} duckhts_simd_backend_entry_t;

static int duckhts_simd_true(void) { return 1; }
static int duckhts_simd_false(void) { return 0; }

static int duckhts_backend_name_equal(const char *a, const char *b) {
    return a && b && strcmp(a, b) == 0;
}

static int duckhts_cpu_has_sse2(void) {
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
    __builtin_cpu_init();
    return __builtin_cpu_supports("sse2") != 0;
#else
    return 0;
#endif
}

static int duckhts_cpu_has_sse41(void) {
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
    __builtin_cpu_init();
    return __builtin_cpu_supports("sse4.1") != 0;
#else
    return 0;
#endif
}

/* Lower priority values are preferred by auto-dispatch.  Placeholder rows
 * (selectable = 0) document known CPU features before DuckHTS has a backend
 * implementation for them. */
static const duckhts_simd_backend_entry_t duckhts_simd_backend_entries[] = {
    {"scalar",       DUCKHTS_SIMD_CAP_SCALAR,       1, duckhts_simd_true,                  duckhts_simd_true,                        1000, &duckhts_simd_table_scalar},
    {"sse2",         DUCKHTS_SIMD_CAP_SSE2,         0, duckhts_simd_false,                 duckhts_cpu_has_sse2,                     40,   (duckhts_simd_dispatch_table_t *)0},
    {"sse41",        DUCKHTS_SIMD_CAP_SSE41,        0, duckhts_simd_false,                 duckhts_cpu_has_sse41,                    30,   (duckhts_simd_dispatch_table_t *)0},
    {"avx2",         DUCKHTS_SIMD_CAP_AVX2,         1, duckhts_simd_avx2_compiled,         duckhts_simd_avx2_cpu_supported,         20,   &duckhts_simd_table_avx2},
    {"avx512",       DUCKHTS_SIMD_CAP_AVX512,       1, duckhts_simd_avx512_compiled,       duckhts_simd_avx512_cpu_supported,       10,   &duckhts_simd_table_avx512},
    {"neon",         DUCKHTS_SIMD_CAP_NEON,         1, duckhts_simd_neon_compiled,         duckhts_simd_neon_cpu_supported,         50,   &duckhts_simd_table_neon},
    {"wasm_simd128", DUCKHTS_SIMD_CAP_WASM_SIMD128, 1, duckhts_simd_wasm_simd128_compiled, duckhts_simd_wasm_simd128_cpu_supported, 60,   &duckhts_simd_table_wasm_simd128}
};

static size_t duckhts_simd_backend_count(void) {
    return sizeof(duckhts_simd_backend_entries) / sizeof(duckhts_simd_backend_entries[0]);
}

const char *duckhts_simd_kernel_name(duckhts_kernel_kind_t kernel) {
    static const char *names[] = {
#define DUCKHTS_SIMD_KERNEL(id, field, sig, name) [DUCKHTS_KERNEL_##id] = name,
#include "duckhts_simd_kernels.def"
#undef DUCKHTS_SIMD_KERNEL
    };
    if ((unsigned)kernel >= (unsigned)DUCKHTS_KERNEL_COUNT) return "unknown";
    return names[kernel] ? names[kernel] : "unknown";
}

const char *duckhts_simd_cap_name(duckhts_simd_cap_t cap) {
    switch (cap) {
    case DUCKHTS_SIMD_CAP_SCALAR: return "scalar";
    case DUCKHTS_SIMD_CAP_SSE2: return "sse2";
    case DUCKHTS_SIMD_CAP_SSE41: return "sse41";
    case DUCKHTS_SIMD_CAP_AVX2: return "avx2";
    case DUCKHTS_SIMD_CAP_AVX512: return "avx512";
    case DUCKHTS_SIMD_CAP_NEON: return "neon";
    case DUCKHTS_SIMD_CAP_WASM_SIMD128: return "wasm_simd128";
    default: return "none";
    }
}

static const duckhts_simd_backend_entry_t *duckhts_simd_find_backend(const char *backend) {
    if (!backend) return (const duckhts_simd_backend_entry_t *)0;
    for (size_t i = 0; i < duckhts_simd_backend_count(); i++) {
        if (duckhts_backend_name_equal(backend, duckhts_simd_backend_entries[i].name)) {
            return &duckhts_simd_backend_entries[i];
        }
    }
    return (const duckhts_simd_backend_entry_t *)0;
}

static const char *duckhts_simd_backend_canonical(const char *backend) {
    const duckhts_simd_backend_entry_t *entry;
    if (!backend || duckhts_backend_name_equal(backend, "auto")) return "auto";
    entry = duckhts_simd_find_backend(backend);
    return entry ? entry->name : (const char *)0;
}

static int duckhts_simd_entry_compiled(const duckhts_simd_backend_entry_t *entry) {
    return entry && entry->compiled && entry->compiled() != 0;
}

static int duckhts_simd_entry_cpu_supported(const duckhts_simd_backend_entry_t *entry) {
    return entry && entry->cpu_supported && entry->cpu_supported() != 0;
}

static int duckhts_simd_entry_available(const duckhts_simd_backend_entry_t *entry) {
    return duckhts_simd_entry_compiled(entry) && duckhts_simd_entry_cpu_supported(entry);
}

static duckhts_simd_cap_t duckhts_simd_available_caps(void) {
    duckhts_simd_cap_t caps = DUCKHTS_SIMD_CAP_SCALAR;
    for (size_t i = 0; i < duckhts_simd_backend_count(); i++) {
        const duckhts_simd_backend_entry_t *entry = &duckhts_simd_backend_entries[i];
        if (entry->selectable && duckhts_simd_entry_available(entry)) caps |= entry->cap;
    }
    return caps;
}

static void duckhts_simd_builder_init(duckhts_simd_builder_t *builder,
                                      duckhts_simd_dispatch_table_t *table,
                                      const char *requested_backend,
                                      duckhts_simd_cap_t viable_caps) {
    memset(table, 0, sizeof(*table));
    table->requested_backend = requested_backend;
    table->viable_caps = viable_caps | DUCKHTS_SIMD_CAP_SCALAR;
    for (int i = 0; i < DUCKHTS_KERNEL_COUNT; i++) {
        table->slots[i].kernel = duckhts_simd_kernel_name((duckhts_kernel_kind_t)i);
        table->slots[i].backend = "unresolved";
        table->slots[i].cap = DUCKHTS_SIMD_CAP_NONE;
        table->slots[i].priority = DUCKHTS_SIMD_PRIORITY_UNSET;
        table->slots[i].scalar_fallback = 0;
        builder->best_priority[i] = DUCKHTS_SIMD_PRIORITY_UNSET;
    }
    builder->table = table;
}

void duckhts_simd_builder_consider_base_counts(duckhts_simd_builder_t *builder,
                                               duckhts_simd_cap_t cap,
                                               const char *backend,
                                               int priority,
                                               duckhts_base_counts_fn fn) {
    duckhts_simd_dispatch_table_t *table;
    duckhts_simd_kernel_slot_t *slot;
    int kernel = DUCKHTS_KERNEL_SEQ_BASE_COUNTS;

    if (!builder || !builder->table || !fn) return;
    table = builder->table;
    if ((table->viable_caps & cap) == 0) return;
    if (priority >= builder->best_priority[kernel]) return;

    table->base_counts = fn;
    slot = &table->slots[kernel];
    slot->kernel = duckhts_simd_kernel_name((duckhts_kernel_kind_t)kernel);
    slot->backend = backend;
    slot->cap = cap;
    slot->priority = priority;
    slot->scalar_fallback = (cap == DUCKHTS_SIMD_CAP_SCALAR &&
                             (table->viable_caps & ~DUCKHTS_SIMD_CAP_SCALAR) != 0);
    builder->best_priority[kernel] = priority;
}

void duckhts_simd_builder_consider_bam_nt16_counts(duckhts_simd_builder_t *builder,
                                                   duckhts_simd_cap_t cap,
                                                   const char *backend,
                                                   int priority,
                                                   duckhts_bam_nt16_counts_fn fn) {
    duckhts_simd_dispatch_table_t *table;
    duckhts_simd_kernel_slot_t *slot;
    int kernel = DUCKHTS_KERNEL_BAM_NT16_COUNTS;

    if (!builder || !builder->table || !fn) return;
    table = builder->table;
    if ((table->viable_caps & cap) == 0) return;
    if (priority >= builder->best_priority[kernel]) return;

    table->bam_nt16_counts = fn;
    slot = &table->slots[kernel];
    slot->kernel = duckhts_simd_kernel_name((duckhts_kernel_kind_t)kernel);
    slot->backend = backend;
    slot->cap = cap;
    slot->priority = priority;
    slot->scalar_fallback = (cap == DUCKHTS_SIMD_CAP_SCALAR &&
                             (table->viable_caps & ~DUCKHTS_SIMD_CAP_SCALAR) != 0);
    builder->best_priority[kernel] = priority;
}

void duckhts_simd_builder_consider_nt16_gc_counts(duckhts_simd_builder_t *builder,
                                                  duckhts_simd_cap_t cap,
                                                  const char *backend,
                                                  int priority,
                                                  duckhts_nt16_gc_counts_fn fn) {
    duckhts_simd_dispatch_table_t *table;
    duckhts_simd_kernel_slot_t *slot;
    int kernel = DUCKHTS_KERNEL_NT16_GC_COUNTS;

    if (!builder || !builder->table || !fn) return;
    table = builder->table;
    if ((table->viable_caps & cap) == 0) return;
    if (priority >= builder->best_priority[kernel]) return;

    table->nt16_gc_counts = fn;
    slot = &table->slots[kernel];
    slot->kernel = duckhts_simd_kernel_name((duckhts_kernel_kind_t)kernel);
    slot->backend = backend;
    slot->cap = cap;
    slot->priority = priority;
    slot->scalar_fallback = (cap == DUCKHTS_SIMD_CAP_SCALAR &&
                             (table->viable_caps & ~DUCKHTS_SIMD_CAP_SCALAR) != 0);
    builder->best_priority[kernel] = priority;
}

void duckhts_simd_builder_consider_fastq_qc(duckhts_simd_builder_t *builder,
                                            duckhts_simd_cap_t cap,
                                            const char *backend,
                                            int priority,
                                            duckhts_fastq_qc_fn fn) {
    duckhts_simd_dispatch_table_t *table;
    duckhts_simd_kernel_slot_t *slot;
    int kernel = DUCKHTS_KERNEL_FASTQ_QC;

    if (!builder || !builder->table || !fn) return;
    table = builder->table;
    if ((table->viable_caps & cap) == 0) return;
    if (priority >= builder->best_priority[kernel]) return;

    table->fastq_qc = fn;
    slot = &table->slots[kernel];
    slot->kernel = duckhts_simd_kernel_name((duckhts_kernel_kind_t)kernel);
    slot->backend = backend;
    slot->cap = cap;
    slot->priority = priority;
    slot->scalar_fallback = (cap == DUCKHTS_SIMD_CAP_SCALAR &&
                             (table->viable_caps & ~DUCKHTS_SIMD_CAP_SCALAR) != 0);
    builder->best_priority[kernel] = priority;
}

static const char *duckhts_simd_selected_label(const duckhts_simd_dispatch_table_t *table) {
    const char *label = (const char *)0;

    if (!table) return "scalar";
    if (!duckhts_backend_name_equal(table->requested_backend, "auto")) return table->requested_backend;

    for (int i = 0; i < DUCKHTS_KERNEL_COUNT; i++) {
        const char *backend = table->slots[i].backend;
        if (!backend || duckhts_backend_name_equal(backend, "unresolved")) continue;
        if (!label) {
            label = backend;
            continue;
        }
        if (!duckhts_backend_name_equal(label, backend)) return "mixed";
    }
    return label ? label : "scalar";
}

static int duckhts_simd_table_resolved(const duckhts_simd_dispatch_table_t *table) {
    if (!table) return 0;
    for (int i = 0; i < DUCKHTS_KERNEL_COUNT; i++) {
        if (table->slots[i].cap == DUCKHTS_SIMD_CAP_NONE ||
            !table->slots[i].backend ||
            duckhts_backend_name_equal(table->slots[i].backend, "unresolved")) {
            return 0;
        }
    }
#define DUCKHTS_SIMD_KERNEL(id, field, sig, name) if (!table->field) return 0;
#include "duckhts_simd_kernels.def"
#undef DUCKHTS_SIMD_KERNEL
    return 1;
}

static void duckhts_simd_require_table_resolved(duckhts_simd_dispatch_table_t *table) {
    const char *requested = table && table->requested_backend ? table->requested_backend : "unknown";
    duckhts_simd_builder_t builder;
    if (!table || duckhts_simd_table_resolved(table)) return;

    fprintf(stderr,
            "duckhts SIMD dispatch table for request '%s' is missing entries; falling back to scalar backend\n",
            requested);

    // Defensive fallback: rebuild the whole table from scalar registrations so
    // future kernels cannot retain NULL function pointers after a partial miss.
    duckhts_simd_builder_init(&builder, table, requested, DUCKHTS_SIMD_CAP_SCALAR);
    duckhts_simd_scalar_register(&builder);
    table->selected_backend = "scalar";
}

static void duckhts_simd_build_table(duckhts_simd_dispatch_table_t *table,
                                     const char *requested_backend,
                                     duckhts_simd_cap_t viable_caps) {
    duckhts_simd_builder_t builder;
    duckhts_simd_builder_init(&builder, table, requested_backend, viable_caps);
    duckhts_simd_scalar_register(&builder);
    duckhts_simd_avx512_register(&builder);
    duckhts_simd_avx2_register(&builder);
    duckhts_simd_neon_register(&builder);
    duckhts_simd_wasm_simd128_register(&builder);
    table->selected_backend = duckhts_simd_selected_label(table);
    duckhts_simd_require_table_resolved(table);
}

static void duckhts_simd_build_all_tables(void) {
    duckhts_simd_build_table(&duckhts_simd_table_auto, "auto", duckhts_simd_available_caps());

    for (size_t i = 0; i < duckhts_simd_backend_count(); i++) {
        const duckhts_simd_backend_entry_t *entry = &duckhts_simd_backend_entries[i];
        duckhts_simd_cap_t caps;
        if (!entry->table) continue;
        caps = DUCKHTS_SIMD_CAP_SCALAR;
        if (entry->selectable && duckhts_simd_entry_available(entry)) caps |= entry->cap;
        duckhts_simd_build_table(entry->table, entry->name, caps);
    }
}

void duckhts_simd_init(void) {
    if (duckhts_simd_is_initialized()) return;
    if (!duckhts_simd_try_begin_init()) {
        duckhts_simd_wait_initialized();
        return;
    }
    duckhts_simd_build_all_tables();
    duckhts_simd_store_table(&duckhts_simd_table_auto);
    duckhts_simd_finish_init();
}

static const duckhts_simd_dispatch_table_t *duckhts_simd_current(void) {
    const duckhts_simd_dispatch_table_t *table;
    if (!duckhts_simd_is_initialized()) duckhts_simd_init();
    table = duckhts_simd_load_table();
    return table ? table : &duckhts_simd_table_scalar;
}

int duckhts_simd_backend_compiled(const char *backend) {
    return duckhts_simd_entry_compiled(duckhts_simd_find_backend(backend));
}

int duckhts_simd_backend_cpu_supported(const char *backend) {
    return duckhts_simd_entry_cpu_supported(duckhts_simd_find_backend(backend));
}

int duckhts_simd_backend_available(const char *backend) {
    return duckhts_simd_entry_available(duckhts_simd_find_backend(backend));
}

static const duckhts_simd_dispatch_table_t *duckhts_simd_resolve_backend_table(const char *backend,
                                                                                  char *err,
                                                                                  size_t err_len) {
    const duckhts_simd_backend_entry_t *entry = (const duckhts_simd_backend_entry_t *)0;
    const duckhts_simd_dispatch_table_t *table;
    const char *canonical;

    if (!duckhts_simd_is_initialized()) duckhts_simd_init();

    canonical = duckhts_simd_backend_canonical(backend);
    if (!canonical) {
        if (err && err_len) {
            snprintf(err, err_len,
                     "unknown SIMD backend '%s' (expected auto, scalar, or a compiled SIMD backend)",
                     backend ? backend : "");
        }
        return (const duckhts_simd_dispatch_table_t *)0;
    }

    if (duckhts_backend_name_equal(canonical, "auto")) {
        table = &duckhts_simd_table_auto;
    }
    else {
        entry = duckhts_simd_find_backend(canonical);
        if (entry && !entry->selectable) {
            if (err && err_len) {
                snprintf(err, err_len,
                         "SIMD backend '%s' is recognized but has no selectable implementation in this build",
                         canonical);
            }
            return (const duckhts_simd_dispatch_table_t *)0;
        }
        if (entry && !duckhts_simd_entry_compiled(entry)) {
            if (err && err_len) {
                snprintf(err, err_len, "SIMD backend '%s' was not compiled into this build", canonical);
            }
            return (const duckhts_simd_dispatch_table_t *)0;
        }
        if (entry && !duckhts_simd_entry_cpu_supported(entry)) {
            if (err && err_len) {
                snprintf(err, err_len, "SIMD backend '%s' is not supported by this CPU/runtime", canonical);
            }
            return (const duckhts_simd_dispatch_table_t *)0;
        }
        table = entry ? entry->table : (const duckhts_simd_dispatch_table_t *)0;
    }

    if (!table) {
        if (err && err_len) snprintf(err, err_len, "SIMD backend '%s' is not available in this process", canonical);
        return (const duckhts_simd_dispatch_table_t *)0;
    }

    return table;
}

int duckhts_simd_set_backend(const char *backend, char *err, size_t err_len) {
    const duckhts_simd_dispatch_table_t *table = duckhts_simd_resolve_backend_table(backend, err, err_len);
    if (!table) return 0;
    duckhts_simd_store_table(table);
    return 1;
}

const char *duckhts_simd_selected_backend(void) {
    return duckhts_simd_current()->selected_backend;
}

const char *duckhts_simd_requested_backend(void) {
    return duckhts_simd_current()->requested_backend;
}

const duckhts_simd_dispatch_table_t *duckhts_simd_dispatch_snapshot(void) {
    return duckhts_simd_current();
}

void duckhts_simd_base_counts_with_table(const duckhts_simd_dispatch_table_t *table,
                                         const char *seq,
                                         size_t len,
                                         duckhts_simd_base_counts_t *out) {
    if (!table) table = duckhts_simd_current();
    table->base_counts(seq, len, out);
}

void duckhts_simd_base_counts(const char *seq, size_t len,
                              duckhts_simd_base_counts_t *out) {
    duckhts_simd_base_counts_with_table(duckhts_simd_dispatch_snapshot(), seq, len, out);
}

void duckhts_simd_bam_nt16_counts_with_table(const duckhts_simd_dispatch_table_t *table,
                                             const uint8_t *packed_seq,
                                             int32_t n_bases,
                                             duckhts_simd_bam_nt16_counts_t *out) {
    if (!table) table = duckhts_simd_current();
    table->bam_nt16_counts(packed_seq, n_bases, out);
}

void duckhts_simd_bam_nt16_counts(const uint8_t *packed_seq, int32_t n_bases,
                                  duckhts_simd_bam_nt16_counts_t *out) {
    duckhts_simd_bam_nt16_counts_with_table(duckhts_simd_dispatch_snapshot(), packed_seq, n_bases, out);
}

void duckhts_simd_nt16_gc_counts_with_table(const duckhts_simd_dispatch_table_t *table,
                                            const uint8_t *codes, size_t n,
                                            duckhts_simd_base_counts_t *out) {
    if (!table) table = duckhts_simd_current();
    table->nt16_gc_counts(codes, n, out);
}

void duckhts_simd_fastq_qc_read_with_table(const duckhts_simd_dispatch_table_t *table,
                                           const char *sequence,
                                           const char *quality,
                                           size_t len,
                                           duckhts_simd_fastq_cycle_t *cycles,
                                           duckhts_simd_fastq_read_t *out) {
    if (!table) table = duckhts_simd_current();
    table->fastq_qc(sequence, quality, len, cycles, out);
}

void duckhts_simd_fastq_qc_read(const char *sequence, const char *quality,
                                size_t len,
                                duckhts_simd_fastq_cycle_t *cycles,
                                duckhts_simd_fastq_read_t *out) {
    duckhts_simd_fastq_qc_read_with_table(duckhts_simd_dispatch_snapshot(),
                                          sequence, quality, len, cycles, out);
}

static inline void set_null_at(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
}

static inline int row_is_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    if (!validity) return 1;
    return duckdb_validity_row_is_valid(validity, row);
}

static char ascii_lower(char c) {
    if (c >= 'A' && c <= 'Z') return (char)(c - 'A' + 'a');
    return c;
}

static int ascii_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

static int normalize_backend_name(const char *src, idx_t src_len, char *dst, size_t dst_len) {
    idx_t start = 0;
    idx_t end = src_len;
    size_t out_len;

    while (start < end && ascii_space(src[start])) start++;
    while (end > start && ascii_space(src[end - 1])) end--;

    out_len = (size_t)(end - start);
    if (out_len == 0 || out_len >= dst_len) return 0;

    for (size_t i = 0; i < out_len; i++) {
        dst[i] = ascii_lower(src[start + (idx_t)i]);
    }
    dst[out_len] = '\0';
    return 1;
}

static int backend_arg_at(duckdb_vector vector, idx_t row, char *dst, size_t dst_len) {
    duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vector);
    duckdb_string_t *value = &data[row];
    return normalize_backend_name(duckdb_string_t_data(value), duckdb_string_t_length(*value), dst, dst_len);
}

typedef struct duckhts_simd_set_backend_bind_data {
    char backend[32];
} duckhts_simd_set_backend_bind_data_t;

typedef struct duckhts_simd_set_backend_init_data {
    int emitted;
} duckhts_simd_set_backend_init_data_t;

static void duckhts_simd_set_backend_bind(duckdb_bind_info info) {
    duckdb_value value = (duckdb_value)0;
    char *raw = (char *)0;
    char backend[32];
    char err[160];
    duckhts_simd_set_backend_bind_data_t *bind;
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);

    if (duckdb_bind_get_parameter_count(info) != 1) {
        duckdb_bind_set_error(info, "duckhts_simd_set_backend: expected exactly one backend argument");
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }

    value = duckdb_bind_get_parameter(info, 0);
    if (!value || duckdb_is_null_value(value)) {
        duckdb_bind_set_error(info, "duckhts_simd_set_backend: backend must be a non-NULL constant string");
        if (value) duckdb_destroy_value(&value);
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }

    raw = duckdb_get_varchar(value);
    duckdb_destroy_value(&value);
    if (!raw) {
        duckdb_bind_set_error(info, "duckhts_simd_set_backend: backend must be a non-NULL constant string");
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }
    if (!normalize_backend_name(raw, (idx_t)strlen(raw), backend, sizeof(backend))) {
        duckdb_free(raw);
        duckdb_bind_set_error(info, "duckhts_simd_set_backend: backend must be a non-empty short string");
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }
    duckdb_free(raw);

    if (!duckhts_simd_resolve_backend_table(backend, err, sizeof(err))) {
        char msg[208];
        snprintf(msg, sizeof(msg), "duckhts_simd_set_backend: %s", err);
        duckdb_bind_set_error(info, msg);
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }

    bind = (duckhts_simd_set_backend_bind_data_t *)duckdb_malloc(sizeof(*bind));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_simd_set_backend: out of memory");
        duckdb_destroy_logical_type(&varchar_type);
        return;
    }
    memcpy(bind->backend, backend, sizeof(bind->backend));
    duckdb_bind_add_result_column(info, "backend", varchar_type);
    duckdb_bind_set_cardinality(info, 1, true);
    duckdb_bind_set_bind_data(info, bind, duckdb_free);
    duckdb_destroy_logical_type(&varchar_type);
}

static void duckhts_simd_backend_scalar(duckdb_function_info info,
                                        duckdb_data_chunk input,
                                        duckdb_vector output) {
    (void)info;
    idx_t row_count = duckdb_data_chunk_get_size(input);
    const char *backend = duckhts_simd_selected_backend();
    idx_t len = (idx_t)strlen(backend);
    for (idx_t row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element_len(output, row, backend, len);
    }
}

static void duckhts_simd_requested_backend_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    (void)info;
    idx_t row_count = duckdb_data_chunk_get_size(input);
    const char *backend = duckhts_simd_requested_backend();
    idx_t len = (idx_t)strlen(backend);
    for (idx_t row = 0; row < row_count; row++) {
        duckdb_vector_assign_string_element_len(output, row, backend, len);
    }
}

static void duckhts_simd_backend_status_scalar(duckdb_function_info info,
                                                duckdb_data_chunk input,
                                                duckdb_vector output,
                                                int (*status_fn)(const char *)) {
    (void)info;
    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(input, 0);
    bool *out_data = (bool *)duckdb_vector_get_data(output);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (idx_t row = 0; row < row_count; row++) {
        char backend[32];
        if (!row_is_valid(backend_vec, row)) {
            set_null_at(output, row);
            continue;
        }
        if (!backend_arg_at(backend_vec, row, backend, sizeof(backend))) {
            out_data[row] = false;
            continue;
        }
        out_data[row] = status_fn(backend) != 0;
    }
}

static void duckhts_simd_backend_compiled_scalar(duckdb_function_info info,
                                                 duckdb_data_chunk input,
                                                 duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_compiled);
}

static void duckhts_simd_backend_cpu_supported_scalar(duckdb_function_info info,
                                                      duckdb_data_chunk input,
                                                      duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_cpu_supported);
}

static void duckhts_simd_backend_available_scalar(duckdb_function_info info,
                                                  duckdb_data_chunk input,
                                                  duckdb_vector output) {
    duckhts_simd_backend_status_scalar(info, input, output, duckhts_simd_backend_available);
}

static void duckhts_simd_set_backend_init(duckdb_init_info info) {
    duckhts_simd_set_backend_init_data_t *init =
        (duckhts_simd_set_backend_init_data_t *)duckdb_malloc(sizeof(*init));
    if (!init) {
        duckdb_init_set_error(info, "duckhts_simd_set_backend: failed to allocate init state");
        return;
    }
    init->emitted = 0;
    duckdb_init_set_init_data(info, init, duckdb_free);
}

static void duckhts_simd_set_backend_scan(duckdb_function_info info, duckdb_data_chunk output) {
    const duckhts_simd_set_backend_bind_data_t *bind =
        (const duckhts_simd_set_backend_bind_data_t *)duckdb_function_get_bind_data(info);
    duckhts_simd_set_backend_init_data_t *init =
        (duckhts_simd_set_backend_init_data_t *)duckdb_function_get_init_data(info);
    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(output, 0);
    char err[160];
    const char *selected;

    if (!init || init->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }
    if (!bind) {
        duckdb_function_set_error(info, "duckhts_simd_set_backend: missing backend bind data");
        duckdb_data_chunk_set_size(output, 0);
        return;
    }
    if (!duckhts_simd_set_backend(bind->backend, err, sizeof(err))) {
        char msg[208];
        snprintf(msg, sizeof(msg), "duckhts_simd_set_backend: %s", err);
        duckdb_function_set_error(info, msg);
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    selected = duckhts_simd_selected_backend();
    duckdb_vector_assign_string_element(backend_vec, 0, selected);
    init->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

static void register_simd_noarg_string_function(duckdb_connection connection,
                                                const char *name,
                                                duckdb_scalar_function_t fn_ptr) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_scalar_function_set_name(fn, name);
    duckdb_scalar_function_set_return_type(fn, varchar_type);
    duckdb_scalar_function_set_volatile(fn);
    duckdb_scalar_function_set_function(fn, fn_ptr);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_scalar_function(&fn);
}

static void register_simd_bool_backend_function(duckdb_connection connection,
                                                const char *name,
                                                duckdb_scalar_function_t fn_ptr) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_scalar_function_set_name(fn, name);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, bool_type);
    duckdb_scalar_function_set_function(fn, fn_ptr);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_scalar_function(&fn);
}

typedef struct duckhts_simd_scan_init {
    idx_t offset;
    const duckhts_simd_dispatch_table_t *table;
} duckhts_simd_scan_init_t;

static void destroy_simd_scan_init(void *ptr) {
    duckdb_free(ptr);
}

static void duckhts_simd_scan_init(duckdb_init_info info) {
    duckhts_simd_scan_init_t *init = (duckhts_simd_scan_init_t *)duckdb_malloc(sizeof(*init));
    if (!init) {
        duckdb_init_set_error(info, "duckhts SIMD table function: failed to allocate init state");
        return;
    }
    init->offset = 0;
    init->table = duckhts_simd_current();
    duckdb_init_set_init_data(info, init, destroy_simd_scan_init);
}

static int duckhts_simd_table_uses_backend(const duckhts_simd_dispatch_table_t *table, const char *backend) {
    if (!table || !backend) return 0;
    for (int i = 0; i < DUCKHTS_KERNEL_COUNT; i++) {
        if (duckhts_backend_name_equal(table->slots[i].backend, backend)) return 1;
    }
    return 0;
}

static void duckhts_simd_info_bind(duckdb_bind_info info) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_bind_add_result_column(info, "backend", varchar_type);
    duckdb_bind_add_result_column(info, "selectable", bool_type);
    duckdb_bind_add_result_column(info, "compiled", bool_type);
    duckdb_bind_add_result_column(info, "cpu_supported", bool_type);
    duckdb_bind_add_result_column(info, "available", bool_type);
    duckdb_bind_add_result_column(info, "selected", bool_type);
    duckdb_bind_add_result_column(info, "requested", bool_type);
    duckdb_bind_add_result_column(info, "dispatch_mode", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
}

static void duckhts_simd_info_scan(duckdb_function_info info, duckdb_data_chunk output) {
    duckhts_simd_scan_init_t *init = (duckhts_simd_scan_init_t *)duckdb_function_get_init_data(info);
    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;
    idx_t n_backends = (idx_t)duckhts_simd_backend_count();
    const duckhts_simd_dispatch_table_t *table;
    const char *requested;

    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(output, 0);
    bool *selectable = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 1));
    bool *compiled = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2));
    bool *cpu_supported = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 3));
    bool *available = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 4));
    bool *selected_out = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 5));
    bool *requested_out = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 6));
    duckdb_vector dispatch_mode_vec = duckdb_data_chunk_get_vector(output, 7);

    if (!init) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }
    table = init->table ? init->table : duckhts_simd_current();
    requested = table->requested_backend;

    while (row_count < vector_size && init->offset < n_backends) {
        const duckhts_simd_backend_entry_t *entry = &duckhts_simd_backend_entries[init->offset];
        const char *backend = entry->name;
        int is_selectable = entry->selectable != 0;
        int is_compiled = duckhts_simd_entry_compiled(entry) != 0;
        int is_cpu_supported = duckhts_simd_entry_cpu_supported(entry) != 0;

        duckdb_vector_assign_string_element(backend_vec, row_count, backend);
        selectable[row_count] = is_selectable;
        compiled[row_count] = is_compiled;
        cpu_supported[row_count] = is_cpu_supported;
        available[row_count] = duckhts_simd_entry_available(entry) != 0;
        selected_out[row_count] = duckhts_simd_table_uses_backend(table, backend) != 0;
        requested_out[row_count] = duckhts_backend_name_equal(backend, requested) != 0;
        duckdb_vector_assign_string_element(dispatch_mode_vec, row_count, duckhts_simd_dispatch_mode_name);

        row_count++;
        init->offset++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

static void duckhts_simd_kernel_info_bind(duckdb_bind_info info) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);

    duckdb_bind_add_result_column(info, "kernel", varchar_type);
    duckdb_bind_add_result_column(info, "selected_backend", varchar_type);
    duckdb_bind_add_result_column(info, "selected_capability", varchar_type);
    duckdb_bind_add_result_column(info, "requested_backend", varchar_type);
    duckdb_bind_add_result_column(info, "scalar_fallback", bool_type);
    duckdb_bind_add_result_column(info, "priority", int_type);
    duckdb_bind_add_result_column(info, "dispatch_mode", varchar_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&int_type);
}

static void duckhts_simd_kernel_info_scan(duckdb_function_info info, duckdb_data_chunk output) {
    duckhts_simd_scan_init_t *init = (duckhts_simd_scan_init_t *)duckdb_function_get_init_data(info);
    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;
    const duckhts_simd_dispatch_table_t *table;

    duckdb_vector kernel_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector backend_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector cap_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector requested_vec = duckdb_data_chunk_get_vector(output, 3);
    bool *scalar_fallback = (bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 4));
    int32_t *priority = (int32_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 5));
    duckdb_vector dispatch_mode_vec = duckdb_data_chunk_get_vector(output, 6);

    if (!init) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }
    table = init->table ? init->table : duckhts_simd_current();

    while (row_count < vector_size && init->offset < (idx_t)DUCKHTS_KERNEL_COUNT) {
        const duckhts_simd_kernel_slot_t *slot = &table->slots[init->offset];

        duckdb_vector_assign_string_element(kernel_vec, row_count, slot->kernel);
        duckdb_vector_assign_string_element(backend_vec, row_count, slot->backend);
        duckdb_vector_assign_string_element(cap_vec, row_count, duckhts_simd_cap_name(slot->cap));
        duckdb_vector_assign_string_element(requested_vec, row_count, table->requested_backend);
        scalar_fallback[row_count] = slot->scalar_fallback != 0;
        priority[row_count] = slot->priority == DUCKHTS_SIMD_PRIORITY_UNSET ? -1 : (int32_t)slot->priority;
        duckdb_vector_assign_string_element(dispatch_mode_vec, row_count, duckhts_simd_dispatch_mode_name);

        row_count++;
        init->offset++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

static void register_simd_info_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_simd_info");
    duckdb_table_function_set_bind(tf, duckhts_simd_info_bind);
    duckdb_table_function_set_init(tf, duckhts_simd_scan_init);
    duckdb_table_function_set_function(tf, duckhts_simd_info_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

static void register_simd_kernel_info_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_simd_kernel_info");
    duckdb_table_function_set_bind(tf, duckhts_simd_kernel_info_bind);
    duckdb_table_function_set_init(tf, duckhts_simd_scan_init);
    duckdb_table_function_set_function(tf, duckhts_simd_kernel_info_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

static void register_simd_set_backend_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_set_name(tf, "duckhts_simd_set_backend");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_set_bind(tf, duckhts_simd_set_backend_bind);
    duckdb_table_function_set_init(tf, duckhts_simd_set_backend_init);
    duckdb_table_function_set_function(tf, duckhts_simd_set_backend_scan);
    duckdb_register_table_function(connection, tf);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_table_function(&tf);
}

void register_duckhts_simd_functions(duckdb_connection connection) {
    register_simd_noarg_string_function(connection, "duckhts_simd_backend", duckhts_simd_backend_scalar);
    register_simd_noarg_string_function(connection, "duckhts_simd_requested_backend", duckhts_simd_requested_backend_scalar);
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_compiled", duckhts_simd_backend_compiled_scalar);
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_cpu_supported", duckhts_simd_backend_cpu_supported_scalar);
    register_simd_bool_backend_function(connection, "duckhts_simd_backend_available", duckhts_simd_backend_available_scalar);
    register_simd_info_function(connection);
    register_simd_kernel_info_function(connection);
    register_simd_set_backend_function(connection);
}
