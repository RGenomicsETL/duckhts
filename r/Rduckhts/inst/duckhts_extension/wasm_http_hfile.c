/*
 * wasm_http_hfile.c -- browser-native http/https hFILE backend for webR.
 *
 * Written against vendored htslib internals. Re-check hFILE_scheme_handler
 * layout and browser acceptance tests on every htslib vendor bump.
 *
 * Transport: synchronous XMLHttpRequest, which is allowed in Web Workers.
 * webR already relies on synchronous XHR in its own worker helpers.
 *
 * This backend does not use libcurl or SOCKFS socket emulation. It is compiled
 * only for Emscripten builds and registered from the DuckDB extension entry
 * point before the first hopen("http://...") attempt.
 *
 * Browser constraints still apply: same-origin URLs work, remote URLs require
 * permissive CORS headers, and S3/GCS remain out of scope for the wasm build.
 *
 * Seek semantics: SEEK_SET and SEEK_CUR are pure C; SEEK_END issues one HEAD
 * request to get Content-Length and caches the result per handle.
 *
 * Write support is intentionally absent. All remote hFILE access is read-only.
 */

#ifdef __EMSCRIPTEN__

#define HTS_BUILDING_LIBRARY  /* enables HTSLIB_EXPORT visibility in headers */

#include "hfile_internal.h"   /* hFILE_backend, hFILE_scheme_handler, hfile_init/destroy */
#include "htslib/hts_defs.h"

#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <emscripten.h>

/* -------------------------------------------------------------------------
 * Our hFILE subclass.  The base member must be first so that (hFILE *) casts
 * to/from (hFILE_wasm_http *) are well-defined.
 * ---------------------------------------------------------------------- */
typedef struct {
    hFILE base;          /* MUST be first */
    char *url;           /* owned, null-terminated */
    off_t http_offset;   /* next byte position for HTTP range requests */
    off_t file_size;     /* Content-Length from HEAD; -1 if not yet fetched */
} hFILE_wasm_http;

/* -------------------------------------------------------------------------
 * wasm_xhr_range_read -- synchronous range GET via XMLHttpRequest.
 *
 * Returns: bytes read (>= 0), 0 at EOF (HTTP 416), or -1 on error.
 *
 * Parameters exposed to EM_ASM:
 *   $0  url_ptr  -- char *, pointer into Emscripten heap
 *   $1  from     -- double, byte offset (double is precise up to 2^53 = 9 PB)
 *   $2  len      -- int,    number of bytes requested
 *   $3  buf      -- uint8_t *, destination in Emscripten heap
 *
 * The JS writes into HEAPU8 at the buf pointer and returns the byte count.
 * ---------------------------------------------------------------------- */
static ssize_t wasm_http_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_wasm_http *fp = (hFILE_wasm_http *)fpv;
    int got;

    if (nbytes == 0) return 0;

    got = EM_ASM_INT({
        var url  = UTF8ToString($0);
        var from = $1;          /* double: start offset */
        var len  = $2;          /* int: bytes requested */
        var buf  = $3;          /* pointer: destination in HEAPU8 */

        var xhr = new XMLHttpRequest();
        xhr.open("GET", url, false);   /* false = synchronous; allowed in Workers */
        xhr.setRequestHeader("Range", "bytes=" + from + "-" + (from + len - 1));
        xhr.responseType = "arraybuffer";
        try {
            xhr.send(null);
        } catch (e) {
            return -1;
        }
        /* 416 = Range Not Satisfiable: caller is past end-of-file */
        if (xhr.status === 416) return 0;
        /* 200 (server ignored Range) or 206 (partial content) are both OK */
        if (xhr.status !== 200 && xhr.status !== 206) return -1;

        var data = new Uint8Array(xhr.response || new ArrayBuffer(0));
        var start = 0;
        if (xhr.status === 200 && from > 0) {
            if (from >= data.length) return 0;
            start = from;
        }
        var available = data.length - start;
        var n = (available < len) ? available : len;
        HEAPU8.set(data.subarray(start, start + n), buf);
        return n;
    }, fp->url, (double)fp->http_offset, (int)nbytes, (uint8_t *)buffer);

    if (got < 0) { errno = EIO; return (ssize_t)-1; }
    fp->http_offset += (off_t)got;
    return (ssize_t)got;
}

/* write: read-only backend */
static ssize_t wasm_http_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    (void)fpv; (void)buffer; (void)nbytes;
    errno = EROFS;
    return (ssize_t)-1;
}

/* -------------------------------------------------------------------------
 * wasm_http_seek -- update the logical read position.
 *
 * SEEK_SET / SEEK_CUR: pure arithmetic, no XHR.
 * SEEK_END: issues one synchronous HEAD to get Content-Length (cached).
 * ---------------------------------------------------------------------- */
static off_t wasm_http_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_wasm_http *fp = (hFILE_wasm_http *)fpv;
    off_t new_offset;

    switch (whence) {
    case SEEK_SET:
        new_offset = offset;
        break;
    case SEEK_CUR:
        new_offset = fp->http_offset + offset;
        break;
    case SEEK_END:
        if (fp->file_size < 0) {
            double sz = EM_ASM_DOUBLE({
                var url = UTF8ToString($0);
                var xhr = new XMLHttpRequest();
                xhr.open("HEAD", url, false);
                try {
                    xhr.send(null);
                } catch (e) {
                    return -1.0;
                }
                if (xhr.status !== 200 && xhr.status !== 206) return -1.0;
                var cl = xhr.getResponseHeader("Content-Length");
                return (cl !== null) ? parseFloat(cl) : -1.0;
            }, fp->url);
            if (sz < 0.0) { errno = ESPIPE; return (off_t)-1; }
            fp->file_size = (off_t)sz;
        }
        new_offset = fp->file_size + offset;
        break;
    default:
        errno = EINVAL;
        return (off_t)-1;
    }

    if (new_offset < 0) { errno = EINVAL; return (off_t)-1; }
    fp->http_offset = new_offset;
    return new_offset;
}

/* flush: no-op for a read-only stream */
static int wasm_http_flush(hFILE *fpv)
{
    (void)fpv;
    return 0;
}

/* close: free url; htslib frees the hFILE struct via hfile_destroy */
static int wasm_http_close(hFILE *fpv)
{
    hFILE_wasm_http *fp = (hFILE_wasm_http *)fpv;
    free(fp->url);
    fp->url = NULL;
    return 0;
}

static const struct hFILE_backend wasm_http_backend = {
    wasm_http_read,
    wasm_http_write,
    wasm_http_seek,
    wasm_http_flush,
    wasm_http_close,
};

/* -------------------------------------------------------------------------
 * Scheme handler open(): allocate hFILE_wasm_http and initialise it.
 * ---------------------------------------------------------------------- */
static hFILE *wasm_http_open(const char *url, const char *mode)
{
    hFILE_wasm_http *fp;
    int access_mode = hfile_oflags(mode) & O_ACCMODE;

    /* Only read mode is supported. */
    if (access_mode != O_RDONLY) {
        errno = EROFS;
        return NULL;
    }

    fp = (hFILE_wasm_http *)hfile_init(sizeof(hFILE_wasm_http), mode, 0);
    if (fp == NULL) return NULL;

    fp->url = strdup(url);
    if (fp->url == NULL) {
        hfile_destroy(&fp->base);
        errno = ENOMEM;
        return NULL;
    }
    fp->http_offset = 0;
    fp->file_size   = -1;
    fp->base.backend = &wasm_http_backend;
    return &fp->base;
}

/* htslib uses the priority field modulo 1000. Any positive value is fine. */
static const struct hFILE_scheme_handler wasm_http_handler = {
    wasm_http_open,
    hfile_always_remote,
    "wasm-xhr",
    50,
};

void register_wasm_http_hfile_backend(void)
{
    static int registered = 0;
    const char *known_schemes[1];
    int nschemes = 0;

    if (registered) return;

    /* Force htslib's lazy scheme-table initialisation before registration. */
    if (hfile_list_schemes(NULL, known_schemes, &nschemes) < 0) return;

    hfile_add_scheme_handler("http",  &wasm_http_handler);
    hfile_add_scheme_handler("https", &wasm_http_handler);
    registered = 1;
}

#endif /* __EMSCRIPTEN__ */
