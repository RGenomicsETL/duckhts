#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 8 ]]; then
    echo "Usage: $0 <bam_path> [index_path] [bin_width] [chrom] [output_prefix] [wise_mapq] [nipter_mapq] [nipter_exclude_duplicate_flag]" >&2
    echo "Example: $0 HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam" >&2
    exit 1
fi

bam_path=$1
index_path=${2:-"${bam_path}.bai"}
bin_width=${3:-500}
chrom_filter=${4:-}
output_prefix=${5:-}
wise_mapq=${6:-1}
nipter_mapq=${7:-1}
nipter_exclude_duplicate_flag=${8:-0}
extension_path=${DUCKHTS_EXTENSION_PATH:-./build/release/duckhts.duckdb_extension}

if [[ ! -f "$bam_path" ]]; then
    echo "BAM not found: $bam_path" >&2
    exit 1
fi

if [[ ! -f "$index_path" ]]; then
    echo "Index not found: $index_path" >&2
    exit 1
fi

if [[ ! -f "$extension_path" ]]; then
    echo "Extension not found: $extension_path" >&2
    exit 1
fi

if [[ -z "$output_prefix" ]]; then
    bam_base=$(basename "$bam_path")
    bam_base=${bam_base%.bam}
    output_prefix="${bam_base}.bins.${bin_width}"
    if [[ -n "$chrom_filter" ]]; then
        safe_chrom=$(printf '%s' "$chrom_filter" | tr -c '[:alnum:]_.-' '_')
        output_prefix="${output_prefix}.${safe_chrom}"
    fi
fi

combined_bed_path="${output_prefix}.combined.bed"
combined_bed_gz_path="${combined_bed_path}.gz"
combined_bed_tbi_path="${combined_bed_gz_path}.tbi"

for path in "$combined_bed_path" "$combined_bed_gz_path" "$combined_bed_tbi_path"; do
    if [[ -e "$path" ]]; then
        echo "Refusing to overwrite existing output: $path" >&2
        exit 1
    fi
done

normalized_chrom_expr="CASE WHEN lower(RNAME) LIKE 'chr%' THEN substr(RNAME, 4) ELSE RNAME END"
canonical_where="normalized_chrom IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')"
normalized_filter=${chrom_filter#chr}

filter_clause="TRUE"
if [[ -n "$chrom_filter" ]]; then
    filter_clause="normalized_chrom = '${normalized_filter}'"
fi

nipter_duplicate_clause="TRUE"
if [[ "$nipter_exclude_duplicate_flag" != "0" ]]; then
    nipter_duplicate_clause="NOT is_duplicate(FLAG)"
fi

read -r -d '' sql_export <<SQL || true
LOAD '${extension_path}';

COPY (
    WITH base_reads AS (
        SELECT
            ${normalized_chrom_expr} AS normalized_chrom,
            POS,
            FLAG,
            MAPQ
        FROM read_bam('${bam_path}', index_path := '${index_path}')
        WHERE RNAME IS NOT NULL
          AND POS IS NOT NULL
          AND POS > 0
    ),
    bins_all AS (
        SELECT
            normalized_chrom AS chrom,
            (((POS - 1) // ${bin_width}) * ${bin_width})::BIGINT AS start,
            ((((POS - 1) // ${bin_width}) + 1) * ${bin_width})::BIGINT AS "end",
            COUNT(*) AS count_total_no_filter
        FROM base_reads
        WHERE ${filter_clause}
        GROUP BY normalized_chrom, ((POS - 1) // ${bin_width})
    ),
    bins_wise AS (
        SELECT
            normalized_chrom AS chrom,
            (((POS - 1) // ${bin_width}) * ${bin_width})::BIGINT AS start,
            COUNT(*) AS count_total_wisecondorx
        FROM base_reads
        WHERE ${canonical_where}
          AND ${filter_clause}
          AND MAPQ >= ${wise_mapq}
          AND is_proper_pair(FLAG)
        GROUP BY normalized_chrom, ((POS - 1) // ${bin_width})
    ),
    bins_nipter AS (
        SELECT
            normalized_chrom AS chrom,
            (((POS - 1) // ${bin_width}) * ${bin_width})::BIGINT AS start,
            COUNT(*) AS count_total_nipter,
            SUM(CASE WHEN is_forward_aligned(FLAG) THEN 1 ELSE 0 END) AS count_fwd_nipter,
            SUM(CASE WHEN is_reverse_complemented(FLAG) AND NOT is_unmapped(FLAG) THEN 1 ELSE 0 END) AS count_rev_nipter
        FROM base_reads
        WHERE ${filter_clause}
          AND MAPQ >= ${nipter_mapq}
          AND NOT is_unmapped(FLAG)
          AND ${nipter_duplicate_clause}
        GROUP BY normalized_chrom, ((POS - 1) // ${bin_width})
    ),
    combined AS (
        SELECT
            a.chrom AS chrom,
            a.start AS start,
            a."end" AS "end",
            a.count_total_no_filter AS count_total_no_filter,
            COALESCE(w.count_total_wisecondorx, 0) AS count_total_wisecondorx,
            COALESCE(n.count_total_nipter, 0) AS count_total_nipter,
            COALESCE(n.count_fwd_nipter, 0) AS count_fwd_nipter,
            COALESCE(n.count_rev_nipter, 0) AS count_rev_nipter
        FROM bins_all a
        LEFT JOIN bins_wise w
          ON a.chrom = w.chrom AND a.start = w.start
        LEFT JOIN bins_nipter n
          ON a.chrom = n.chrom AND a.start = n.start
    ),
    metadata_lines AS (
        SELECT 0 AS ord, printf('#duckhts_benchmark_bin_counts_sql') AS line
        UNION ALL SELECT 1, printf('#bam=%s', '${bam_path}')
        UNION ALL SELECT 2, printf('#index=%s', '${index_path}')
        UNION ALL SELECT 3, printf('#bin_width=%s', '${bin_width}')
        UNION ALL SELECT 4, printf('#chrom_filter_raw=%s', '${chrom_filter:-ALL}')
        UNION ALL SELECT 5, printf('#chrom_filter_normalized=%s', '${normalized_filter:-ALL}')
        UNION ALL SELECT 6, printf('#chrom_normalization=strip_leading_chr_prefix')
        UNION ALL SELECT 7, printf('#wisecondorx_profile=canonical_contigs_only,mapq>=%s,require_proper_pair=TRUE,no_adjacent_duplicate_suppression', '${wise_mapq}')
        UNION ALL SELECT 8, printf('#nipter_profile=mapq>=%s,exclude_unmapped=TRUE,exclude_duplicate_flag=%s,proper_pair_not_required,strand_split=TRUE', '${nipter_mapq}', '${nipter_exclude_duplicate_flag}')
        UNION ALL SELECT 9, printf('#columns=chrom,start,end,count_total_no_filter,count_total_wisecondorx,count_total_nipter,count_fwd_nipter,count_rev_nipter')
        UNION ALL SELECT 10, '#chrom	start	end	count_total_no_filter	count_total_wisecondorx	count_total_nipter	count_fwd_nipter	count_rev_nipter'
    ),
    data_lines AS (
        SELECT
            11 + row_number() OVER (ORDER BY chrom, start) AS ord,
            chrom || '	' ||
            CAST(start AS VARCHAR) || '	' ||
            CAST("end" AS VARCHAR) || '	' ||
            CAST(count_total_no_filter AS VARCHAR) || '	' ||
            CAST(count_total_wisecondorx AS VARCHAR) || '	' ||
            CAST(count_total_nipter AS VARCHAR) || '	' ||
            CAST(count_fwd_nipter AS VARCHAR) || '	' ||
            CAST(count_rev_nipter AS VARCHAR) AS line
        FROM combined
    )
    SELECT line
    FROM (
        SELECT ord, line FROM metadata_lines
        UNION ALL
        SELECT ord, line FROM data_lines
    ) lines
    ORDER BY ord
) TO '${combined_bed_path}' (HEADER FALSE, QUOTE '', ESCAPE '');

WITH base_reads AS (
    SELECT
        ${normalized_chrom_expr} AS normalized_chrom,
        POS,
        FLAG,
        MAPQ
    FROM read_bam('${bam_path}', index_path := '${index_path}')
    WHERE RNAME IS NOT NULL
      AND POS IS NOT NULL
      AND POS > 0
)
SELECT
    COUNT(*) FILTER (WHERE ${filter_clause}) AS reads_total_no_filter,
    COUNT(*) FILTER (
        WHERE ${canonical_where}
          AND ${filter_clause}
          AND MAPQ >= ${wise_mapq}
          AND is_proper_pair(FLAG)
    ) AS reads_total_wisecondorx,
    COUNT(*) FILTER (
        WHERE ${filter_clause}
          AND MAPQ >= ${nipter_mapq}
          AND NOT is_unmapped(FLAG)
          AND ${nipter_duplicate_clause}
    ) AS reads_total_nipter,
    SUM(CASE
        WHEN ${filter_clause}
         AND MAPQ >= ${nipter_mapq}
         AND ${nipter_duplicate_clause}
         AND is_forward_aligned(FLAG) THEN 1 ELSE 0
    END) AS reads_fwd_nipter,
    SUM(CASE
        WHEN ${filter_clause}
         AND MAPQ >= ${nipter_mapq}
         AND ${nipter_duplicate_clause}
         AND is_reverse_complemented(FLAG)
         AND NOT is_unmapped(FLAG) THEN 1 ELSE 0
    END) AS reads_rev_nipter
FROM base_reads;
SQL

read -r -d '' sql_index <<SQL || true
LOAD '${extension_path}';
SELECT * FROM bgzip('${combined_bed_path}', output_path := '${combined_bed_gz_path}', keep := TRUE, overwrite := FALSE);
SELECT * FROM tabix_index('${combined_bed_gz_path}', preset := 'bed', comment_char := '#');
SQL

echo "Benchmarking SQL binning baseline"
echo "  bam:         $bam_path"
echo "  index:       $index_path"
echo "  bin_width:   $bin_width"
echo "  wise_mapq:   $wise_mapq"
echo "  nipter_mapq: $nipter_mapq"
echo "  nipter_exclude_duplicate_flag: $nipter_exclude_duplicate_flag"
if [[ -n "$chrom_filter" ]]; then
    echo "  chrom:       ${normalized_filter} (normalized from ${chrom_filter})"
fi
echo "  output:      $combined_bed_path"

/usr/bin/time -f 'elapsed=%E user=%U sys=%S maxrss_kb=%M' \
    duckdb -unsigned -c "$sql_export"

duckdb -unsigned -c "$sql_index"
