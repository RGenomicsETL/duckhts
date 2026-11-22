#ifndef DUCKHTS_HTS_FILE_ACCESS_H
#define DUCKHTS_HTS_FILE_ACCESS_H

/*
 * File access for every DuckHTS reader.
 *
 * A file access is the pair of settings a connection reads its inputs with:
 *
 *   duckhts_file_system        'htslib' (the default) keeps htslib's own file
 *                              access; 'duckdb' reads through the DuckDB file
 *                              system of the connection, with its secrets,
 *                              settings and registered file systems.
 *   duckhts_remote_block_size  bytes htslib asks its file backend for in each
 *                              read of a remote file; 0 (the default) keeps
 *                              the size DuckHTS chooses for the scan.
 *
 * A table function captures it once, at bind, and keeps it in its bind data;
 * a scalar function captures it in the state of each worker. Every open of
 * an input and of its index then goes through the functions below,
 * which take the paths as the user wrote them. With 'htslib' each function is
 * the htslib call it is named after. With 'duckdb' it opens the same paths
 * through the connection (duckdb_fs_hfile.h) and finds an index by opening
 * the names htslib would look for beside the file; it never reads or writes
 * the working directory and never builds a missing index.
 *
 * A NULL file access is htslib's own file access with the automatic block
 * size. Writers and index builders use it: the DuckDB route is read-only.
 *
 * A file access is immutable and may be used from any thread. It must be
 * destroyed, and the files it opened through DuckDB closed, before its
 * connection closes.
 */

#include <stdbool.h>

#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>

#include "duckdb_extension.h"
#include "hts_io_tuning.h"

#define DUCKHTS_FILE_SYSTEM_OPTION "duckhts_file_system"
#define DUCKHTS_REMOTE_BLOCK_SIZE_OPTION "duckhts_remote_block_size"

/*
 * Bounds of a remote block size in bytes. htslib opens every file with a
 * 128 KiB block and cannot shrink a block that already holds data, so a
 * smaller value would have no effect. Each open handle allocates one block.
 */
#define DUCKHTS_REMOTE_BLOCK_SIZE_MIN (128 * 1024)
#define DUCKHTS_REMOTE_BLOCK_SIZE_MAX (256 * 1024 * 1024)

typedef struct duckhts_file_access duckhts_file_access_t;

/* Registers the two settings. Returns false when DuckDB rejects one. */
bool duckhts_register_file_access_options(duckdb_connection connection);

/*
 * Captures the file access of the connection that binds a table function
 * (duckhts_file_access_bind) or that initialises the local state of a scalar
 * function (duckhts_file_access_scalar_init). On an invalid setting each sets
 * the bind or init error, prefixed with function_name, and returns false with
 * *access NULL. The caller releases *access with duckhts_file_access_destroy.
 */
bool duckhts_file_access_bind(duckdb_bind_info info, const char *function_name,
                              duckhts_file_access_t **access);
bool duckhts_file_access_scalar_init(duckdb_init_info info, const char *function_name,
                                     duckhts_file_access_t **access);

void duckhts_file_access_destroy(duckhts_file_access_t *access);

/*
 * True when access reads through the DuckDB file system. A cache that keeps
 * htslib handles beyond one query may hold only handles for which this is
 * false: the others belong to the connection that opened them.
 */
bool duckhts_file_access_uses_duckdb(const duckhts_file_access_t *access);

/*
 * Whether path is read over a network: htslib's hisremote(path) with
 * 'htslib'; with 'duckdb', whether path starts with a URL scheme other than
 * file, because the DuckDB file system has no such query. Remote tuning
 * applies only to remote paths.
 */
bool duckhts_file_access_is_remote(const duckhts_file_access_t *access, const char *path);

/*
 * hts_open(path, "r"), also for sam_open, bcf_open and vcf_open, then the
 * remote tuning of profile and the block size of the file access.
 */
htsFile *duckhts_file_access_open(const duckhts_file_access_t *access, const char *path,
                                  duckhts_hts_io_profile_t profile);

/*
 * bcf_sr_add_reader(readers, path), for a synced reader that allows a missing
 * index (BCF_SR_ALLOW_NO_IDX). Returns the htslib status.
 */
int duckhts_file_access_sr_add_reader(const duckhts_file_access_t *access, bcf_srs_t *readers,
                                      const char *path);

/* hopen(path, "r"), then the remote tuning of profile and the block size. */
hFILE *duckhts_file_access_hopen(const duckhts_file_access_t *access, const char *path,
                                 duckhts_hts_io_profile_t profile);

/*
 * hts_readlist(path, 1, n): the non-empty lines of a plain or BGZF-compressed
 * text file, or NULL when it cannot be read. The caller frees each string and
 * the array with free(). With 'duckdb' the file is read through the route, and
 * a failure writes DuckDB's reason into reason as " (<reason>)"; reason is
 * empty otherwise. The signature matches filter_readlist_fn, with access as
 * its data.
 */
char **duckhts_file_access_readlist(void *access, const char *path, int *n,
                                    char *reason, size_t reason_size);

/*
 * sam_index_load3, bcf_index_load3 and tbx_index_load3. With 'duckdb' and no
 * index_path the candidates are, in htslib's order, .csi then .bai (BAM, SAM)
 * or .tbi (tabix), each appended to the path and then in place of its
 * extension; .crai for CRAM; .csi alone for BCF. HTS_IDX_SAVE_REMOTE has no
 * effect with 'duckdb'.
 */
hts_idx_t *duckhts_file_access_sam_index_load(const duckhts_file_access_t *access, htsFile *fp,
                                              const char *path, const char *index_path, int flags);
hts_idx_t *duckhts_file_access_bcf_index_load(const duckhts_file_access_t *access,
                                              const char *path, const char *index_path, int flags);
tbx_t *duckhts_file_access_tbx_index_load(const duckhts_file_access_t *access,
                                          const char *path, const char *index_path, int flags);

/*
 * fai_load3_format, then the remote cache tuning of profile. A NULL fai_path
 * or gzi_path is the path with .fai or .gzi appended. FAI_CREATE has no
 * effect with 'duckdb': a missing index is an error.
 */
faidx_t *duckhts_file_access_fai_load(const duckhts_file_access_t *access, const char *path,
                                      const char *fai_path, const char *gzi_path, int flags,
                                      enum fai_format_options format,
                                      duckhts_hts_io_profile_t profile);

/*
 * The reference of an open alignment file. locator is a FASTA path or
 * "fasta" HTS_IDX_DELIM "fai" (duckhts_reference_locator). The first is
 * hts_set_opt(fp, CRAM_OPT_REFERENCE, locator); the second is
 * hts_set_fai_filename(fp, locator). Both return the htslib status. With
 * 'duckdb' the FASTA, its .fai and its .gzi are read through the connection;
 * a reference htslib finds by REF_PATH, REF_CACHE or the header's UR field
 * stays with htslib.
 */
int duckhts_file_access_set_cram_reference(const duckhts_file_access_t *access, htsFile *fp,
                                           const char *locator);
int duckhts_file_access_set_fai_filename(const duckhts_file_access_t *access, htsFile *fp,
                                         const char *locator);

/*
 * DuckDB's reason for the most recent failure on path, as " (reason)" ready
 * to append to an error message, or "" with 'htslib' and when the most recent
 * failure of the connection's route concerns another path. The text belongs
 * to the calling thread and is valid until its next call.
 */
const char *duckhts_file_access_reason(const duckhts_file_access_t *access, const char *path);

/*
 * duckhts_file_access_reason() for a FASTA that duckhts_file_access_fai_load
 * could not open. The failure may concern the FASTA, its index (fai_path, or
 * <fasta>.fai when NULL) or its BGZF index (gzi_path, or <fasta>.gzi when
 * NULL); the reason text names the file it concerns. Same lifetime as
 * duckhts_file_access_reason().
 */
const char *duckhts_file_access_fasta_reason(const duckhts_file_access_t *access, const char *fasta_path,
                                             const char *fai_path, const char *gzi_path);

#endif
