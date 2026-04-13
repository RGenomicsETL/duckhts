#ifndef DUCKHTS_WASM_SOCKET_COMPAT_H
#define DUCKHTS_WASM_SOCKET_COMPAT_H

#ifdef __EMSCRIPTEN__

#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef DUCKHTS_WASM_DUCKDB_RUNTIME
/*
 * duckdb-wasm's MAIN_MODULE exports legalized wrappers for several libc
 * symbols that use 64-bit arguments/returns, plus native i64 versions under
 * orig$*. Bind the community-extension side-module path to the native
 * symbols so its ABI matches the duckdb-wasm host. webR does not need this
 * remap, so it is gated separately from __EMSCRIPTEN__.
 */
extern off_t duckhts_wasm_orig_lseek(int fd, off_t offset, int whence) __asm__("orig$lseek");
extern int duckhts_wasm_orig_ftruncate(int fd, off_t length) __asm__("orig$ftruncate");
extern long long duckhts_wasm_orig_strtoll(const char *nptr, char **endptr, int base) __asm__("orig$strtoll");
extern unsigned long long duckhts_wasm_orig_strtoull(const char *nptr, char **endptr, int base) __asm__("orig$strtoull");
extern long long duckhts_wasm_orig_atoll(const char *nptr) __asm__("orig$atoll");
extern time_t duckhts_wasm_orig_time(time_t *timer) __asm__("orig$time");

#define lseek duckhts_wasm_orig_lseek
#define ftruncate duckhts_wasm_orig_ftruncate
#define strtoll duckhts_wasm_orig_strtoll
#define strtoull duckhts_wasm_orig_strtoull
#define atoll duckhts_wasm_orig_atoll
#define time duckhts_wasm_orig_time
#endif

static inline ssize_t duckhts_wasm_recv(int fd, void *buf, size_t len, int flags) {
    (void) flags;
    return read(fd, buf, len);
}

static inline ssize_t duckhts_wasm_send(int fd, const void *buf, size_t len, int flags) {
    (void) flags;
    return write(fd, buf, len);
}

static inline int duckhts_wasm_closesocket(int fd) {
    return close(fd);
}

#define recv duckhts_wasm_recv
#define send duckhts_wasm_send
#define closesocket duckhts_wasm_closesocket

#endif

#endif
