#include "bcftools_shim.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

static void check_true(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "CHECK_FAILED: %s\n", message);
        failures++;
    }
}

static void probe_guarded_error(void) {
    int jumped = duckhts_filter_try_begin();
    if (jumped == 0) {
        duckhts_bcftools_error("guarded %s", "boom");
        check_true(0, "guarded duckhts_bcftools_error returned instead of longjmp");
        duckhts_filter_try_end();
        return;
    }

    check_true(strcmp(duckhts_filter_last_error(), "guarded boom") == 0,
               "guarded error text was not preserved");
    duckhts_filter_try_end();
    printf("GUARDED_LONGJMP=1\n");
}

static void probe_guarded_errno_error(void) {
    int jumped = duckhts_filter_try_begin();
    if (jumped == 0) {
        errno = ENOENT;
        duckhts_bcftools_error_errno("guarded errno");
        check_true(0, "guarded duckhts_bcftools_error_errno returned instead of longjmp");
        duckhts_filter_try_end();
        return;
    }

    check_true(strstr(duckhts_filter_last_error(), "guarded errno") != NULL,
               "guarded errno prefix was not preserved");
    check_true(strstr(duckhts_filter_last_error(), "No such file") != NULL,
               "guarded errno strerror text was not appended");
    duckhts_filter_try_end();
    printf("GUARDED_ERRNO_LONGJMP=1\n");
}

static void probe_unguarded_error_return(void) {
    int caller_state = 1;

    duckhts_bcftools_error("unguarded boom");
    caller_state = 2;

    check_true(caller_state == 2,
               "unguarded duckhts_bcftools_error did not return to caller");
    check_true(strcmp(duckhts_filter_last_error(), "unguarded boom") == 0,
               "unguarded error text was not preserved");
    printf("UNGUARDED_RETURNED=1\n");
    printf("UNGUARDED_CALLER_STATE=%d\n", caller_state);
}

static void probe_unguarded_errno_error_return(void) {
    int caller_state = 1;

    errno = ENOENT;
    duckhts_bcftools_error_errno("unguarded errno");
    caller_state = 2;

    check_true(caller_state == 2,
               "unguarded duckhts_bcftools_error_errno did not return to caller");
    check_true(strstr(duckhts_filter_last_error(), "unguarded errno") != NULL,
               "unguarded errno prefix was not preserved");
    check_true(strstr(duckhts_filter_last_error(), "No such file") != NULL,
               "unguarded errno strerror text was not appended");
    printf("UNGUARDED_ERRNO_RETURNED=1\n");
    printf("UNGUARDED_ERRNO_CALLER_STATE=%d\n", caller_state);
}

int main(void) {
    probe_guarded_error();
    probe_guarded_errno_error();
    probe_unguarded_error_return();
    probe_unguarded_errno_error_return();

    if (failures != 0) {
        fprintf(stderr, "FAILURES=%d\n", failures);
        return 1;
    }

    printf("BCFTOOLS_SHIM_ERROR_PATH_PROBE_OK=1\n");
    return 0;
}
