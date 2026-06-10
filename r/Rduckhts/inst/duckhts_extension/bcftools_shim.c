#define DUCKHTS_BCFTOOLS_SHIM_IMPLEMENTATION
#include "bcftools_shim.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>

typedef struct {
    jmp_buf env_stack[16];
    int depth;
} duckhts_filter_try_t;

#if defined(_MSC_VER)
#define DUCKHTS_THREAD_LOCAL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#define DUCKHTS_THREAD_LOCAL __thread
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define DUCKHTS_THREAD_LOCAL _Thread_local
#else
#define DUCKHTS_THREAD_LOCAL
#endif

static DUCKHTS_THREAD_LOCAL duckhts_filter_try_t g_filter_try = {0};
static DUCKHTS_THREAD_LOCAL char g_filter_last_error[1024] = {0};

char *bcftools_version(void) {
    return "bcftools-shim";
}

void duckhts_bcftools_error(const char *format, ...) {
    va_list ap;
    g_filter_last_error[0] = '\0';
    va_start(ap, format);
    vsnprintf(g_filter_last_error, sizeof(g_filter_last_error), format, ap);
    va_end(ap);
    if (g_filter_try.depth > 0) {
        longjmp(g_filter_try.env_stack[g_filter_try.depth - 1], 1);
    }
    fputs(g_filter_last_error, stderr);
    fputc('\n', stderr);
}

void duckhts_bcftools_error_errno(const char *format, ...) {
    va_list ap;
    size_t n;
    g_filter_last_error[0] = '\0';
    va_start(ap, format);
    vsnprintf(g_filter_last_error, sizeof(g_filter_last_error), format, ap);
    va_end(ap);
    if (errno != 0) {
        n = strlen(g_filter_last_error);
        if (n + 3 < sizeof(g_filter_last_error)) {
            snprintf(g_filter_last_error + n, sizeof(g_filter_last_error) - n, ": %s", strerror(errno));
        }
    }
    if (g_filter_try.depth > 0) {
        longjmp(g_filter_try.env_stack[g_filter_try.depth - 1], 1);
    }
    fputs(g_filter_last_error, stderr);
    fputc('\n', stderr);
}

int duckhts_filter_try_begin_impl(void) {
    g_filter_last_error[0] = '\0';
    if (g_filter_try.depth >= (int)(sizeof(g_filter_try.env_stack) / sizeof(g_filter_try.env_stack[0]))) {
        snprintf(g_filter_last_error, sizeof(g_filter_last_error), "bcftools_score: filter try depth exceeded");
        // Overflow did not push an env, so the caller must not run its
        // longjmp/try_end recovery against a frame it never owns. Unwind to the
        // nearest active env when one exists; only signal a plain error return
        // when there is no enclosing try to recover into.
        if (g_filter_try.depth > 0) {
            longjmp(g_filter_try.env_stack[g_filter_try.depth - 1], 1);
        }
        return 1;
    }
    g_filter_try.depth++;
    return 0;
}

jmp_buf *duckhts_filter_try_current_env(void) {
    if (g_filter_try.depth <= 0) return NULL;
    return &g_filter_try.env_stack[g_filter_try.depth - 1];
}

void duckhts_filter_try_end(void) {
    if (g_filter_try.depth > 0) g_filter_try.depth--;
}

const char *duckhts_filter_last_error(void) {
    return g_filter_last_error;
}
