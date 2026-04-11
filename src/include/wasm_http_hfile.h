/*
 * wasm_http_hfile.h -- browser-native http/https hFILE backend declaration.
 *
 * Only the registration call is needed by callers; the backend is self-contained
 * in wasm_http_hfile.c and compiled away entirely on non-Emscripten targets.
 */

#ifndef WASM_HTTP_HFILE_H
#define WASM_HTTP_HFILE_H

/*
 * register_wasm_http_hfile_backend() - register browser-native http and https
 * handlers with htslib's scheme table.
 *
 * On __EMSCRIPTEN__ builds this installs the synchronous XHR-backed handlers.
 * On all other builds it is a no-op inline stub.
 */

#ifdef __EMSCRIPTEN__
void register_wasm_http_hfile_backend(void);
#else
static inline void register_wasm_http_hfile_backend(void) {}
#endif

#endif /* WASM_HTTP_HFILE_H */
