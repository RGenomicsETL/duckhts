#ifndef DUCKHTS_OUTPUT_FILE_H
#define DUCKHTS_OUTPUT_FILE_H

#include <errno.h>
#include <fcntl.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#include <process.h>
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
    char *temporary_path;
} duckhts_output_owner_t;

/* New files use 0666 subject to umask, like fopen("wb") and htslib.
 * Overwrites are staged in the destination directory until finish succeeds. */
static inline int duckhts_output_open(const char *path, int overwrite, duckhts_output_owner_t *owner) {
    static atomic_ulong sequence = 0;
    int flags = O_WRONLY | O_CREAT | O_EXCL;
    int fd = -1;
    unsigned long pid;
    memset(owner, 0, sizeof(*owner));
#ifdef _WIN32
    flags |= _O_BINARY;
    pid = (unsigned long)_getpid();
#else
    pid = (unsigned long)getpid();
#ifdef O_NOFOLLOW
    flags |= O_NOFOLLOW;
#endif
#endif
    if (overwrite) {
        size_t len = strlen(path);
        if (len > SIZE_MAX - 80) {
            errno = ENAMETOOLONG;
            return -1;
        }
        owner->temporary_path = malloc(len + 80);
        if (!owner->temporary_path) {
            errno = ENOMEM;
            return -1;
        }
    }
    for (int attempt = 0; attempt < 128; attempt++) {
        const char *destination = path;
        if (overwrite) {
            unsigned long serial = atomic_fetch_add(&sequence, 1);
            sprintf(owner->temporary_path, "%s.duckhts-%lu-%lu.tmp", path, pid, serial);
            destination = owner->temporary_path;
        }
#ifdef _WIN32
        fd = _open(destination, flags, _S_IREAD | _S_IWRITE);
#else
        fd = open(destination, flags, 0666);
#endif
        if (fd >= 0 || !overwrite || errno != EEXIST) break;
    }
    if (fd < 0) {
        int saved_errno = errno;
        free(owner->temporary_path);
        owner->temporary_path = NULL;
        errno = saved_errno;
        return -1;
    }
    owner->opened = 1;
#ifdef _WIN32
    BY_HANDLE_FILE_INFORMATION info;
    HANDLE handle = (HANDLE)_get_osfhandle(fd);
    if (GetFileInformationByHandle(handle, &info)) {
        owner->device = (unsigned long long)info.dwVolumeSerialNumber;
        owner->inode = ((unsigned long long)info.nFileIndexHigh << 32) | info.nFileIndexLow;
        owner->identity_valid = 1;
    }
#else
    struct stat st;
    if (fstat(fd, &st) == 0) {
        owner->device = (unsigned long long)st.st_dev;
        owner->inode = (unsigned long long)st.st_ino;
        owner->identity_valid = 1;
    }
#endif
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

static inline const char *duckhts_output_path(const char *path, const duckhts_output_owner_t *owner) {
    return owner->temporary_path ? owner->temporary_path : path;
}

/* Compare the opened file's identity with the current path before cleanup. */
static inline void duckhts_output_cleanup(const char *path, duckhts_output_owner_t *owner) {
    path = duckhts_output_path(path, owner);
    if (owner->opened && owner->identity_valid && path) {
#ifdef _WIN32
        BY_HANDLE_FILE_INFORMATION info;
        HANDLE handle = CreateFileA(path, FILE_READ_ATTRIBUTES,
                                    FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                                    NULL, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OPEN_REPARSE_POINT, NULL);
        if (handle != INVALID_HANDLE_VALUE) {
            int matches = GetFileInformationByHandle(handle, &info) &&
                !(info.dwFileAttributes & (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)) &&
                owner->device == (unsigned long long)info.dwVolumeSerialNumber &&
                owner->inode == (((unsigned long long)info.nFileIndexHigh << 32) | info.nFileIndexLow);
            CloseHandle(handle);
            if (matches) remove(path);
        }
#else
        struct stat st;
        if (lstat(path, &st) == 0 && S_ISREG(st.st_mode) &&
            owner->device == (unsigned long long)st.st_dev &&
            owner->inode == (unsigned long long)st.st_ino) unlink(path);
#endif
    }
    free(owner->temporary_path);
    memset(owner, 0, sizeof(*owner));
}

/* Call only after closing and checking the writer. A failed publish also removes
 * the staged file; it never unlinks the destination being replaced. */
static inline int duckhts_output_finish(const char *path, duckhts_output_owner_t *owner, int success) {
    if (success && owner->opened && owner->temporary_path) {
#ifdef _WIN32
        if (!MoveFileExA(owner->temporary_path, path, MOVEFILE_REPLACE_EXISTING)) {
            errno = EIO;
#else
        if (rename(owner->temporary_path, path) != 0) {
#endif
            int saved_errno = errno;
            duckhts_output_cleanup(path, owner);
            errno = saved_errno;
            return -1;
        }
    }
    if (!success) {
        duckhts_output_cleanup(path, owner);
    } else {
        free(owner->temporary_path);
        memset(owner, 0, sizeof(*owner));
    }
    return 0;
}

#endif
