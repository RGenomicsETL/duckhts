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
#include "htslib/hts_log.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <emscripten.h>

#define WASM_HTTP_LARGE_FULL_DOWNLOAD_BYTES (64.0 * 1024.0 * 1024.0)

/* -------------------------------------------------------------------------
 * Our hFILE subclass.  The base member must be first so that (hFILE *) casts
 * to/from (hFILE_wasm_http *) are well-defined.
 * ---------------------------------------------------------------------- */
typedef struct {
    hFILE base;          /* MUST be first */
    char *url;           /* owned, null-terminated */
    off_t http_offset;   /* next byte position for HTTP range requests */
    off_t file_size;     /* Content-Length from HEAD; -1 if not yet fetched */
    int warned_no_range; /* whether we've warned that Range is ignored */
    int warned_large_full_download; /* whether we've warned about large 200 payload */
} hFILE_wasm_http;

/*
 * Optional browser-side configuration object:
 *   Module.duckhtsWasmHttpConfig = {
 *     headers: {"Authorization": "Bearer ...", "X-Foo": "bar"},
 *     allowHosts: ["example.org", ".ebi.ac.uk"],
 *     enforceHostAllowlist: false,
 *     withCredentials: false,
 *     allowInsecureAuth: false
 *   }
 *
 * Security model:
 * - Custom headers are applied only when allowHosts matches URL hostname.
 * - Optional hard host allowlist: enforceHostAllowlist=true blocks requests to
 *   hosts that do not match allowHosts.
 * - Authorization is stripped on non-HTTPS URLs unless allowInsecureAuth=true.
 */

static double wasm_http_discover_size(const char *url, uintptr_t cache_key)
{
    return EM_ASM_DOUBLE({
        var url = UTF8ToString($0);
        var key = String($1 >>> 0);
        var cache = Module.duckhtsWasmHttpFullObjectCache;
        var cfg = Module.duckhtsWasmHttpConfig || null;
        var xhr = new XMLHttpRequest();
        var cl = null;
        var cr = null;
        var totalStr = null;
        var slash = -1;

        function hostMatchesAllowlist(targetUrl, allowHosts) {
            var host = "";
            var allow = allowHosts;
            var i;
            if (!allow) return false;
            if (typeof allow === "string") {
                allow = allow.split(",");
            }
            try {
                host = new URL(targetUrl).hostname.toLowerCase();
            } catch (e) {
                return false;
            }
            if (!Array.isArray(allow)) return false;
            for (i = 0; i < allow.length; i++) {
                var raw = allow[i];
                if (raw === null || raw === undefined) continue;
                var rule = String(raw).trim().toLowerCase();
                if (!rule) continue;
                if (rule.charAt(0) === ".") {
                    if (host.length > rule.length && host.endsWith(rule)) return true;
                } else {
                    if (host === rule) return true;
                }
            }
            return false;
        }

        function applyConfiguredHeaders(xhrObj, targetUrl) {
            var headers;
            var keyName;
            var isHttps = (typeof targetUrl === "string") && targetUrl.toLowerCase().startsWith("https://");
            var allowInsecureAuth = !!(cfg && cfg.allowInsecureAuth);
            if (!cfg || !cfg.headers) return;
            if (!hostMatchesAllowlist(targetUrl, cfg.allowHosts)) return;
            headers = cfg.headers;
            for (keyName in headers) {
                if (!Object.prototype.hasOwnProperty.call(headers, keyName)) continue;
                if (headers[keyName] === null || headers[keyName] === undefined) continue;
                if (!allowInsecureAuth && keyName.toLowerCase() === "authorization" && !isHttps) continue;
                try {
                    xhrObj.setRequestHeader(String(keyName), String(headers[keyName]));
                } catch (e) {
                }
            }
            if (cfg.withCredentials === true) xhrObj.withCredentials = true;
        }

        function shouldAllowRequest(targetUrl) {
            if (!cfg || cfg.enforceHostAllowlist !== true) return true;
            return hostMatchesAllowlist(targetUrl, cfg.allowHosts);
        }

        if (!shouldAllowRequest(url)) {
            return -1.0;
        }

        if (!cache) {
            cache = Object.create(null);
            Module.duckhtsWasmHttpFullObjectCache = cache;
        }

        if (cache[key]) {
            return cache[key].length;
        }

        xhr.open("HEAD", url, false);
        applyConfiguredHeaders(xhr, url);
        try {
            xhr.send(null);
        } catch (e) {
            xhr = null;
        }

        if (xhr && (xhr.status === 200 || xhr.status === 206)) {
            cl = xhr.getResponseHeader("Content-Length");
            if (cl !== null) {
                var clv = parseFloat(cl);
                if (!isNaN(clv) && clv >= 0) return clv;
            }
        }

        xhr = new XMLHttpRequest();
        xhr.open("GET", url, false);
        xhr.setRequestHeader("Range", "bytes=0-0");
        applyConfiguredHeaders(xhr, url);
        xhr.responseType = "arraybuffer";
        try {
            xhr.send(null);
        } catch (e) {
            return -1.0;
        }

        if (xhr.status !== 200 && xhr.status !== 206) return -1.0;

        cr = xhr.getResponseHeader("Content-Range");
        if (cr !== null) {
            slash = cr.lastIndexOf("/");
            if (slash >= 0) {
                totalStr = cr.substring(slash + 1).trim();
                if (totalStr !== "*") {
                    var total = parseFloat(totalStr);
                    if (!isNaN(total) && total >= 0) return total;
                }
            }
        }

        cl = xhr.getResponseHeader("Content-Length");
        if (cl !== null && xhr.status === 200) {
            var cl2 = parseFloat(cl);
            if (!isNaN(cl2) && cl2 >= 0) return cl2;
        }

        if (xhr.status === 200) {
            var data = new Uint8Array(xhr.response || new ArrayBuffer(0));
            cache[key] = data;
            return data.length;
        }

        return -1.0;
    }, url, (unsigned int)cache_key);
}

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
    double discovered_size = -1.0;
    double probed_size;
    double full_response_bytes = 0.0;
    int request_len;
    int got;
    int range_ignored = 0;

    if (nbytes == 0) return 0;

    if (fp->file_size < 0 && fp->http_offset > 0) {
        probed_size = wasm_http_discover_size(fp->url, (uintptr_t)fp);
        if (probed_size >= 0.0) fp->file_size = (off_t)probed_size;
    }

    if (fp->file_size >= 0 && fp->http_offset >= fp->file_size) return 0;

    request_len = (nbytes > (size_t)INT_MAX) ? INT_MAX : (int)nbytes;

    got = EM_ASM_INT({
        var url  = UTF8ToString($0);
        var from = $1;          /* double: start offset */
        var len  = $2;          /* int: bytes requested */
        var buf  = $3;          /* pointer: destination in HEAPU8 */
        var outp = $4;          /* pointer: discovered size (double) */
        var knownSize = $5;     /* double: cached size, -1 if unknown */
        var outRangeIgnored = $6; /* pointer: int flag */
        var outFullBytes = $7;    /* pointer: full response bytes (double) */
        var key = String($8 >>> 0); /* per-handle cache key */
        var end = from + len - 1;
        var cl = null;
        var cr = null;
        var totalStr = null;
        var slash = -1;
        var cache = Module.duckhtsWasmHttpFullObjectCache;
        var cached = null;
        var cfg = Module.duckhtsWasmHttpConfig || null;

        function hostMatchesAllowlist(targetUrl, allowHosts) {
            var host = "";
            var allow = allowHosts;
            var i;
            if (!allow) return false;
            if (typeof allow === "string") {
                allow = allow.split(",");
            }
            try {
                host = new URL(targetUrl).hostname.toLowerCase();
            } catch (e) {
                return false;
            }
            if (!Array.isArray(allow)) return false;
            for (i = 0; i < allow.length; i++) {
                var raw = allow[i];
                if (raw === null || raw === undefined) continue;
                var rule = String(raw).trim().toLowerCase();
                if (!rule) continue;
                if (rule.charAt(0) === ".") {
                    if (host.length > rule.length && host.endsWith(rule)) return true;
                } else {
                    if (host === rule) return true;
                }
            }
            return false;
        }

        function applyConfiguredHeaders(xhrObj, targetUrl) {
            var headers;
            var keyName;
            var isHttps = (typeof targetUrl === "string") && targetUrl.toLowerCase().startsWith("https://");
            var allowInsecureAuth = !!(cfg && cfg.allowInsecureAuth);
            if (!cfg || !cfg.headers) return;
            if (!hostMatchesAllowlist(targetUrl, cfg.allowHosts)) return;
            headers = cfg.headers;
            for (keyName in headers) {
                if (!Object.prototype.hasOwnProperty.call(headers, keyName)) continue;
                if (headers[keyName] === null || headers[keyName] === undefined) continue;
                if (!allowInsecureAuth && keyName.toLowerCase() === "authorization" && !isHttps) continue;
                try {
                    xhrObj.setRequestHeader(String(keyName), String(headers[keyName]));
                } catch (e) {
                }
            }
            if (cfg.withCredentials === true) xhrObj.withCredentials = true;
        }

        function shouldAllowRequest(targetUrl) {
            if (!cfg || cfg.enforceHostAllowlist !== true) return true;
            return hostMatchesAllowlist(targetUrl, cfg.allowHosts);
        }

        HEAPF64[outp >> 3] = -1.0;
        HEAP32[outRangeIgnored >> 2] = 0;
        HEAPF64[outFullBytes >> 3] = 0.0;

        if (!shouldAllowRequest(url)) {
            return -2;
        }

        if (!cache) {
            cache = Object.create(null);
            Module.duckhtsWasmHttpFullObjectCache = cache;
        }

        cached = cache[key] || null;
        if (cached) {
            HEAPF64[outp >> 3] = cached.length;
            HEAP32[outRangeIgnored >> 2] = 1;
            HEAPF64[outFullBytes >> 3] = cached.length;
            if (from >= cached.length) return 0;
            var cn = ((cached.length - from) < len) ? (cached.length - from) : len;
            HEAPU8.set(cached.subarray(from, from + cn), buf);
            return cn;
        }

        if (knownSize >= 0 && from >= knownSize) return 0;
        if (knownSize >= 0 && end >= knownSize) end = knownSize - 1;
        if (end < from) return 0;

        var xhr = new XMLHttpRequest();
        xhr.open("GET", url, false);   /* false = synchronous; allowed in Workers */
        xhr.setRequestHeader("Range", "bytes=" + from + "-" + end);
        applyConfiguredHeaders(xhr, url);
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

        if (xhr.status === 206) {
            cr = xhr.getResponseHeader("Content-Range");
            if (cr !== null) {
                slash = cr.lastIndexOf("/");
                if (slash >= 0) {
                    totalStr = cr.substring(slash + 1).trim();
                    if (totalStr !== "*") {
                        var total = parseFloat(totalStr);
                        if (!isNaN(total) && total >= 0) HEAPF64[outp >> 3] = total;
                    }
                }
            }
        } else {
            cl = xhr.getResponseHeader("Content-Length");
            if (from === 0 && cl !== null) {
                var clv = parseFloat(cl);
                if (!isNaN(clv) && clv >= 0) HEAPF64[outp >> 3] = clv;
            }
            HEAP32[outRangeIgnored >> 2] = 1;
        }

        var data = new Uint8Array(xhr.response || new ArrayBuffer(0));
        if (xhr.status === 200) {
            HEAPF64[outp >> 3] = data.length;
            HEAPF64[outFullBytes >> 3] = data.length;
            cache[key] = data;
        }
        var start = 0;
        if (xhr.status === 200 && from > 0) {
            if (from >= data.length) return 0;
            start = from;
        }
        var available = data.length - start;
        var n = (available < len) ? available : len;
        HEAPU8.set(data.subarray(start, start + n), buf);
        return n;
    }, fp->url, (double)fp->http_offset, request_len, (uint8_t *)buffer,
       &discovered_size, (double)fp->file_size, &range_ignored, &full_response_bytes,
       (unsigned int)(uintptr_t)fp);

    if (got < 0) {
        errno = (got == -2) ? EACCES : EIO;
        return (ssize_t)-1;
    }
    if (range_ignored && !fp->warned_no_range) {
        hts_log_warning(
            "wasm_http_hfile: server ignored HTTP Range for %s; "
            "the browser backend will fetch full-object responses and slice bytes locally",
            fp->url);
        fp->warned_no_range = 1;
    }
    if (range_ignored && full_response_bytes >= WASM_HTTP_LARGE_FULL_DOWNLOAD_BYTES
        && !fp->warned_large_full_download) {
        hts_log_warning(
            "wasm_http_hfile: large fallback download %.2f MiB for %s "
            "(offset=%lld request=%d); consider a Range-capable endpoint to avoid repeated full downloads",
            full_response_bytes / (1024.0 * 1024.0), fp->url,
            (long long)fp->http_offset, request_len);
        fp->warned_large_full_download = 1;
    }
    if (discovered_size >= 0.0) fp->file_size = (off_t)discovered_size;
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
            double sz = wasm_http_discover_size(fp->url, (uintptr_t)fp);
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
    EM_ASM({
        var key = String($0 >>> 0);
        var cache = Module.duckhtsWasmHttpFullObjectCache;
        if (cache && cache[key]) {
            delete cache[key];
        }
    }, (unsigned int)(uintptr_t)fp);
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
    fp->warned_no_range = 0;
    fp->warned_large_full_download = 0;
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

#else

/* Keep non-Emscripten builds non-empty under -Wpedantic. */
extern int duckhts_wasm_http_hfile_non_emscripten_anchor;

#endif /* __EMSCRIPTEN__ */
