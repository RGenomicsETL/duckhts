/*
 * Stands for a host process that carries its own htslib.
 *
 * Preloaded into the process that loads DuckHTS, this library defines hopen()
 * and seq_nt16_str, which the dynamic loader finds before the copies bundled
 * in the extension. The bundled hts_open() (hts.c) calls hopen() (hfile.c),
 * and read_fastq decodes each base through seq_nt16_str. Either reference
 * lands here if it goes through the dynamic symbol table.
 */
#include <errno.h>
#include <stddef.h>

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

/* Every 4-bit base code decodes to X; htslib's table is "=ACMGRSVTWYHKDBN". */
const char seq_nt16_str[] = "XXXXXXXXXXXXXXXX";

int htslib_binding_host_hopen_calls(void) {
    return host_hopen_calls;
}
