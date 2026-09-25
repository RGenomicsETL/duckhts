#ifndef DUCKHTS_OUTPUT_FILE_H
#define DUCKHTS_OUTPUT_FILE_H

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#ifdef _WIN32
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

typedef struct {
    int opened;
    int identity_valid;
    unsigned long long device;
    unsigned long long inode;
} duckhts_output_owner_t;

/* New files use 0666 subject to umask, like fopen("wb") and htslib. */
static inline int duckhts_output_open(const char *path, int overwrite, duckhts_output_owner_t *owner) {
    int flags = O_WRONLY | O_CREAT | (overwrite ? O_TRUNC : O_EXCL);
    int fd;
    owner->opened = 0;
    owner->identity_valid = 0;
#ifdef _WIN32
    fd = _open(path, flags | _O_BINARY, _S_IREAD | _S_IWRITE);
    if (fd >= 0) {
        BY_HANDLE_FILE_INFORMATION info;
        HANDLE handle = (HANDLE)_get_osfhandle(fd);
        if (GetFileInformationByHandle(handle, &info)) {
            owner->device = (unsigned long long)info.dwVolumeSerialNumber;
            owner->inode = ((unsigned long long)info.nFileIndexHigh << 32) | info.nFileIndexLow;
            owner->identity_valid = 1;
        }
    }
#else
    fd = open(path, flags, 0666);
    if (fd >= 0) {
        struct stat st;
        if (fstat(fd, &st) == 0) {
            owner->device = (unsigned long long)st.st_dev;
            owner->inode = (unsigned long long)st.st_ino;
            owner->identity_valid = 1;
        }
    }
#endif
    if (fd >= 0) owner->opened = 1;
    return fd;
}

static inline void duckhts_output_close_fd(int fd) {
#ifdef _WIN32
    _close(fd);
#else
    close(fd);
#endif
}

static inline FILE *duckhts_output_fdopen(int fd) {
#ifdef _WIN32
    return _fdopen(fd, "wb");
#else
    return fdopen(fd, "wb");
#endif
}

/* Compare the opened file's identity with the current destination before cleanup. */
static inline void duckhts_output_cleanup(const char *path, duckhts_output_owner_t owner) {
    if (!owner.opened || !owner.identity_valid || !path) return;
#ifdef _WIN32
    BY_HANDLE_FILE_INFORMATION info;
    HANDLE handle = CreateFileA(path, FILE_READ_ATTRIBUTES,
                                FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                                NULL, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, NULL);
    if (handle == INVALID_HANDLE_VALUE) return;
    if (GetFileInformationByHandle(handle, &info) &&
        owner.device == (unsigned long long)info.dwVolumeSerialNumber &&
        owner.inode == (((unsigned long long)info.nFileIndexHigh << 32) | info.nFileIndexLow)) {
        CloseHandle(handle);
        remove(path);
        return;
    }
    CloseHandle(handle);
#else
    struct stat st;
    if (stat(path, &st) == 0 && owner.device == (unsigned long long)st.st_dev &&
        owner.inode == (unsigned long long)st.st_ino) remove(path);
#endif
}

#endif
