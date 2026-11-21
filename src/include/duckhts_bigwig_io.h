#ifndef DUCKHTS_BIGWIG_IO_H
#define DUCKHTS_BIGWIG_IO_H

#include "bigWig.h"
#include "hts_file_access.h"

/*
 * bwOpen(path, NULL, "r") with the file read through access
 * (hts_file_access.h). path is borrowed by the returned file and must outlive
 * it. Returns NULL when the file does not open or is not a bigWig.
 */
bigWigFile_t *duckhts_bigwig_open(const duckhts_file_access_t *access, const char *path);

#endif
