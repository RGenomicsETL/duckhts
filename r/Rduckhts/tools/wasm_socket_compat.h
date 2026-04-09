#ifndef RDUCKHTS_WASM_SOCKET_COMPAT_H
#define RDUCKHTS_WASM_SOCKET_COMPAT_H

#ifdef __EMSCRIPTEN__

#include <sys/types.h>
#include <unistd.h>

static inline ssize_t rduckhts_wasm_recv(int fd, void *buf, size_t len, int flags) {
    (void) flags;
    return read(fd, buf, len);
}

static inline ssize_t rduckhts_wasm_send(int fd, const void *buf, size_t len, int flags) {
    (void) flags;
    return write(fd, buf, len);
}

static inline int rduckhts_wasm_closesocket(int fd) {
    return close(fd);
}

#define recv rduckhts_wasm_recv
#define send rduckhts_wasm_send
#define closesocket rduckhts_wasm_closesocket

#endif

#endif
