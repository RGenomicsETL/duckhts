#ifndef DUCKHTS_WASM_SOCKET_COMPAT_H
#define DUCKHTS_WASM_SOCKET_COMPAT_H

#ifdef __EMSCRIPTEN__

#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

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
