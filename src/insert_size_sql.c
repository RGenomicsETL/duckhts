/* Weighted summaries from Fulcrum Genomics Riker (MIT), v0.4.1/40d5e326:
 * counter.rs:110-115 mode ties; 157-175 population moments; 186-240
 * interpolated median and doubled integer MAD; commands/isize.rs:245-290
 * orientation filtering, moment trim and shared histogram cutoff.
 * Histograms stay sparse: SQL never expands read_pairs into individual rows. */
#include "duckhts_registration.h"

bool register_duckhts_insert_size_sql(duckhts_registration_t *registration) {
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
        "floor(max(median+deviations*mad) OVER ())::BIGINT AS histogram_cutoff FROM selected) "
        "SELECT * FROM result WHERE (SELECT ok FROM options)=1 "
        "UNION ALL SELECT NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL "
        "WHERE (SELECT ok FROM options)=0)"
    );
}
