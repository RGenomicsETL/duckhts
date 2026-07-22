/* Read-only libBigWig transport implemented with htslib hFILE. */

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hfile.h>

#include "hts_io_tuning.h"
#include "bigWig.h"

size_t GLOBAL_DEFAULTBUFFERSIZE;

static hFILE *
duckhts_bigwig_hfile(URL_t *url)
{
    return (hFILE *)url->x.fp;
}

static int
duckhts_bigwig_fill(URL_t *url, size_t position)
{
    hFILE *file;
    off_t current;
    size_t filled;

    if (url == NULL || url->x.fp == NULL || url->memBuf == NULL ||
        url->bufSize == 0u || (uintmax_t)position > INT64_MAX)
        return 0;
    file = duckhts_bigwig_hfile(url);
    current = htell(file);
    if (current < 0 || (uintmax_t)current != (uintmax_t)position) {
        if (hseek(file, (off_t)position, SEEK_SET) < 0)
            return 0;
    }
    filled = 0u;
    while (filled < url->bufSize) {
        ssize_t count = hread(file, (char *)url->memBuf + filled,
            url->bufSize - filled);

        if (count < 0)
            return 0;
        if (count == 0)
            break;
        filled += (size_t)count;
    }
    url->filePos = position;
    url->bufPos = 0u;
    url->bufLen = filled;
    return filled != 0u;
}

size_t
urlRead(URL_t *url, void *buffer, size_t size)
{
    size_t copied;

    if (url == NULL || url->x.fp == NULL)
        return 0u;
    if (url->type == BWG_FILE) {
        ssize_t count = hread(duckhts_bigwig_hfile(url), buffer, size);
        return count < 0 ? 0u : (size_t)count;
    }
    copied = 0u;
    while (copied < size) {
        size_t available, take, next_position;

        if (url->bufPos >= url->bufLen) {
            if (url->filePos == SIZE_MAX)
                next_position = 0u;
            else {
                if (url->bufLen > SIZE_MAX - url->filePos)
                    return copied;
                next_position = url->filePos + url->bufLen;
            }
            if (!duckhts_bigwig_fill(url, next_position))
                return copied;
        }
        available = url->bufLen - url->bufPos;
        take = size - copied < available ? size - copied : available;
        memcpy((char *)buffer + copied,
            (char *)url->memBuf + url->bufPos, take);
        url->bufPos += take;
        copied += take;
    }
    return copied;
}

CURLcode
urlSeek(URL_t *url, size_t position)
{
    off_t offset;

    if (url == NULL || url->x.fp == NULL)
        return CURLE_FAILED_INIT;
    if (url->type != BWG_FILE && url->filePos != SIZE_MAX &&
        position >= url->filePos &&
        position - url->filePos <= url->bufLen) {
        url->bufPos = position - url->filePos;
        return CURLE_OK;
    }
    offset = (off_t)position;
    if (offset < 0 || (uintmax_t)offset != (uintmax_t)position)
        return CURLE_FAILED_INIT;
    if (url->type != BWG_FILE) {
        url->filePos = position;
        url->bufPos = 0u;
        url->bufLen = 0u;
        return CURLE_OK;
    }
    return hseek(duckhts_bigwig_hfile(url), offset, SEEK_SET) < 0
        ? CURLE_FAILED_INIT : CURLE_OK;
}

int64_t
urlTell(URL_t *url)
{
    off_t offset;

    if (url == NULL || url->x.fp == NULL)
        return -1;
    if (url->type != BWG_FILE && url->filePos != SIZE_MAX) {
        if ((uintmax_t)url->filePos > INT64_MAX ||
            url->bufPos > (size_t)INT64_MAX - url->filePos)
            return -1;
        return (int64_t)(url->filePos + url->bufPos);
    }
    offset = htell(duckhts_bigwig_hfile(url));
    if (offset < 0 || (uintmax_t)offset > INT64_MAX)
        return -1;
    return (int64_t)offset;
}

URL_t *
urlOpen(const char *path, CURLcode (*callback)(CURL *), const char *mode)
{
    URL_t *url;
    hFILE *file;

    (void)callback;
    if (path == NULL || (mode != NULL && strchr(mode, 'w') != NULL)) {
        errno = EINVAL;
        return NULL;
    }
    file = hopen(path, "r");
    if (file == NULL)
        return NULL;
    duckhts_apply_remote_hfile_tuning(file, path,
        DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    url = calloc(1u, sizeof(*url));
    if (url == NULL) {
        hclose_abruptly(file);
        return NULL;
    }
    url->x.fp = (FILE *)file;
    url->fname = path;
    url->type = hisremote(path) ? BWG_HTTPS : BWG_FILE;
    url->filePos = SIZE_MAX;
    if (url->type != BWG_FILE) {
        url->bufSize = GLOBAL_DEFAULTBUFFERSIZE;
        url->memBuf = malloc(url->bufSize);
        if (url->bufSize == 0u || url->memBuf == NULL) {
            hclose_abruptly(file);
            free(url->memBuf);
            free(url);
            return NULL;
        }
    }
    return url;
}

void
urlClose(URL_t *url)
{
    int close_status;

    if (url == NULL)
        return;
    close_status = url->x.fp == NULL ? 0 :
        hclose(duckhts_bigwig_hfile(url));
    (void)close_status;
    free(url->memBuf);
    free(url);
}

/* DuckHTS exposes a reader only.  This satisfies upstream bwClose() without
 * pulling the writer implementation and rejects accidental write finalization. */
int
bwFinalize(bigWigFile_t *file)
{
    return file != NULL && file->isWrite ? 1 : 0;
}
