#ifndef DUCKHTS_ATTRIBUTES_H
#define DUCKHTS_ATTRIBUTES_H

/* Mark a printf-style function so the compiler checks its callers' formats.
 * MinGW's GCC treats the printf archetype as Microsoft's format rules, which reject
 * %lld/%llu; use its C99-conforming stdio rules instead. Clang (Windows ARM64) has no
 * gnu_printf archetype and accepts %llu as printf. */
#if defined(__MINGW32__) && !defined(__clang__)
#define DUCKHTS_PRINTF(format_index, first_arg) \
    __attribute__((format(gnu_printf, format_index, first_arg)))
#else
#define DUCKHTS_PRINTF(format_index, first_arg) \
    __attribute__((format(printf, format_index, first_arg)))
#endif

#endif
