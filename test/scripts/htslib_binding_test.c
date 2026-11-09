/*
 * The host process must not be able to replace the htslib that DuckHTS bundles.
 *
 * This program stands for a host that carries its own htslib: it defines
 * hopen(), which the dynamic loader finds before the extension's copy. The
 * bundled hts_open() (hts.c) calls hopen() (hfile.c). If that call goes
 * through the dynamic symbol table, it lands here.
 */
#include <errno.h>
#include <stdio.h>

#include "htslib/hfile.h"
#include "htslib/hts.h"

static int host_hopen_calls;

hFILE *hopen(const char *filename, const char *mode, ...) {
    (void)filename;
    (void)mode;
    host_hopen_calls++;
    errno = ENOENT;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: htslib_binding_test MISSING_FILE\n");
        return 2;
    }
    htsFile *file = hts_open(argv[1], "r");
    if (file != NULL) {
        fprintf(stderr, "hts_open opened a file that must not exist: %s\n", argv[1]);
        hts_close(file);
        return 2;
    }
    if (host_hopen_calls != 0) {
        fprintf(stderr, "bundled hts_open called the host's hopen %d time(s)\n", host_hopen_calls);
        return 1;
    }
    puts("htslib_binding_test: OK");
    return 0;
}
