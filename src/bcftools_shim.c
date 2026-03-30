#include "bcftools_shim.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>

typedef struct {
    jmp_buf env;
    int active;
} duckhts_filter_try_t;

static duckhts_filter_try_t g_filter_try = {0};
static char g_filter_last_error[1024] = {0};

char *bcftools_version(void) {
    return "bcftools-shim";
}

void duckhts_bcftools_error(const char *format, ...) {
    va_list ap;
    g_filter_last_error[0] = '\0';
    va_start(ap, format);
    vsnprintf(g_filter_last_error, sizeof(g_filter_last_error), format, ap);
    va_end(ap);
    if (g_filter_try.active) {
        longjmp(g_filter_try.env, 1);
    }
    fputs(g_filter_last_error, stderr);
    abort();
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
    if (g_filter_try.active) {
        longjmp(g_filter_try.env, 1);
    }
    fputs(g_filter_last_error, stderr);
    fputc('\n', stderr);
    abort();
}

int duckhts_filter_try_begin(void) {
    g_filter_try.active = 1;
    g_filter_last_error[0] = '\0';
    return setjmp(g_filter_try.env);
}

void duckhts_filter_try_end(void) {
    g_filter_try.active = 0;
}

const char *duckhts_filter_last_error(void) {
    return g_filter_last_error;
}
