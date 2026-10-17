/* Weighted summaries from Fulcrum Genomics Riker (MIT), v0.4.1/40d5e326:
 * counter.rs:110-115 mode ties; 157-175 population moments; 186-240
 * interpolated median and doubled integer MAD; commands/isize.rs:245-290
 * orientation filtering, moment trim and shared histogram cutoff.
 * Histograms stay sparse: SQL never expands read_pairs into individual rows. */
#include "duckhts_registration.h"

bool register_duckhts_insert_size_sql(duckhts_registration_t *registration) {
    /* riker 0.4.1 isize.rs:197-218 filters/counts;
     * sam/pair_orientation.rs:30-59 compares inclusive alignment ends. */
    if (!duckhts_register_sql(registration,
        "CREATE OR REPLACE MACRO duckhts_bam_insert_size_counts_sql("
        "path, reference := NULL, include_duplicates := false, "
        "decompression_threads := 0, scan_mode := 'auto', index_path := NULL) AS TABLE "
        "WITH projected AS ("
        " SELECT FLAG, RNAME, RNEXT, POS::BIGINT AS pos, PNEXT::BIGINT AS pnext, "
        " TLEN::BIGINT AS tlen, CIGAR FROM read_bam(path, reference := reference, "
        " cigar_representation := 'binary', decompression_threads := decompression_threads, "
        " scan_mode := scan_mode, index_path := index_path)"
        "), eligible AS ("
        " SELECT *, is_reverse_complemented(FLAG) AS rev, "
        " is_next_segment_reverse_complemented(FLAG) AS mate_rev "
        " FROM projected WHERE is_paired(FLAG) AND NOT is_unmapped(FLAG) "
        " AND NOT is_next_segment_unmapped(FLAG) AND NOT is_first_segment(FLAG) "
        " AND NOT is_secondary(FLAG) AND NOT is_supplementary(FLAG) "
        " AND NOT is_qc_fail(FLAG) "
        " AND (include_duplicates OR NOT is_duplicate(FLAG)) "
        " AND tlen <> 0 AND (RNEXT = '=' OR RNEXT = RNAME)"
        "), oriented AS ("
        " SELECT CASE WHEN rev = mate_rev THEN 'TANDEM' "
        " WHEN rev THEN CASE WHEN pnext < pos + "
        " greatest(coalesce(cigar_reference_length(CIGAR), 0), 1) - 1 "
        " THEN 'FR' ELSE 'RF' END "
        " ELSE CASE WHEN pos < pos + tlen THEN 'FR' ELSE 'RF' END END AS pair_orientation, "
        " abs(tlen) AS insert_size FROM eligible"
        ") SELECT pair_orientation, insert_size, count(*)::UBIGINT AS read_pairs "
        "FROM oriented GROUP BY pair_orientation, insert_size")) {
        return false;
    }
    return duckhts_register_sql(registration,
        "CREATE OR REPLACE MACRO duckhts_insert_size_metrics(counts_table, "
        "deviations := 10.0, min_frac := 0.05, sample_name := NULL) AS TABLE ("
        "WITH options AS (SELECT CASE WHEN deviations IS NULL OR NOT isfinite(deviations) OR deviations < 0 "
        "OR min_frac IS NULL OR NOT isfinite(min_frac) OR min_frac < 0 "
        "THEN error('deviations and min_frac must be finite and nonnegative') ELSE 1 END AS ok), "
        "checked AS (SELECT CASE WHEN pair_orientation IS NULL OR pair_orientation NOT IN ('FR','RF','TANDEM') "
        "OR insert_size IS NULL OR insert_size < 1 OR insert_size > 4294967295 OR insert_size != trunc(insert_size) "
        "OR read_pairs IS NULL OR read_pairs < 1 OR read_pairs != trunc(read_pairs) "
        "THEN error('counts require FR/RF/TANDEM, whole insert_size 1..UINT32_MAX and positive whole read_pairs') "
        "ELSE pair_orientation END AS o, insert_size::BIGINT AS x, read_pairs::HUGEINT AS c FROM query_table(counts_table)), "
        "hist AS (SELECT o,x,sum(c) AS c FROM checked GROUP BY o,x), "
        "ranked AS (SELECT *,sum(c) OVER (PARTITION BY o ORDER BY x) AS r, "
        "sum(c) OVER (PARTITION BY o) AS n FROM hist), "
        "med AS (SELECT o,max(n) AS n,(min(x) FILTER (WHERE r > (n-1)//2)::DOUBLE + "
        "min(x) FILTER (WHERE r > n//2)::DOUBLE)/2 AS median FROM ranked GROUP BY o), "
        "dev AS (SELECT h.o,least(4294967295,round(abs(h.x-m.median)*2))::BIGINT AS d,sum(h.c) AS c "
        "FROM hist h JOIN med m USING(o) GROUP BY h.o,d), "
        "dev_ranked AS (SELECT *,sum(c) OVER (PARTITION BY o ORDER BY d) AS r FROM dev), "
        "mad AS (SELECT m.o,m.n,m.median,(min(d) FILTER (WHERE r > (m.n-1)//2)::DOUBLE + "
        "min(d) FILTER (WHERE r > m.n//2)::DOUBLE)/4 AS mad FROM med m JOIN dev_ranked USING(o) "
        "GROUP BY m.o,m.n,m.median), "
        "stats AS (SELECT m.o,m.n,m.median,m.mad,min(x) AS lo,max(x) AS hi, "
        "first(x ORDER BY c DESC,x DESC) AS mode, "
        "sum(c) FILTER (WHERE x <= median+deviations*mad) AS trimmed_n, "
        "sum(x::DOUBLE*c::DOUBLE ORDER BY x) FILTER (WHERE x <= median+deviations*mad) AS sx, "
        "sum(x::DOUBLE*x::DOUBLE*c::DOUBLE ORDER BY x) FILTER (WHERE x <= median+deviations*mad) AS sx2 "
        "FROM hist JOIN mad m USING(o) GROUP BY m.o,m.n,m.median,m.mad), "
        "selected AS (SELECT * FROM stats WHERE n::DOUBLE >= (SELECT sum(n)::DOUBLE FROM med)*min_frac), "
        "result AS (SELECT sample_name::VARCHAR AS sample,o AS pair_orientation,n AS read_pairs, "
        "median AS median_insert_size,mad AS median_absolute_deviation, "
        "coalesce(sx/trimmed_n::DOUBLE,0) AS mean_insert_size, "
        "sqrt(greatest(0,coalesce(sx2/trimmed_n::DOUBLE-pow(sx/trimmed_n::DOUBLE,2),0))) AS standard_deviation, "
        "mode AS mode_insert_size,lo AS min_insert_size,hi AS max_insert_size, "
        "floor(least(4294967295,max(median+deviations*mad) OVER ()))::BIGINT AS histogram_cutoff FROM selected) "
        "SELECT * FROM result WHERE (SELECT ok FROM options)=1 "
        "UNION ALL SELECT NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL "
        "WHERE (SELECT ok FROM options)=0)"
    );
}
