-- Independent projected-read_bam oracle: Fulcrum Genomics riker v0.4.1
-- (MIT), tag 40d5e326. Picard is context only, not the oracle.
-- Caller sets insert_path, insert_reference and insert_duplicates; DuckDB's
-- threads setting controls the projected reader's processing parallelism.
WITH projected AS (
    SELECT flag, rname, pos,
           pos + greatest(coalesce(cigar_reference_length(cigar), 0), 1) - 1 AS end_pos,
           rnext, pnext, tlen::BIGINT AS tlen
    FROM read_bam(getvariable('insert_path'),
                  reference := getvariable('insert_reference'))
), eligible AS (
    -- commands/isize.rs:197-206: no MAPQ or proper-pair predicate; reject
    -- first segments rather than requiring the second-segment flag.
    SELECT * FROM projected
    WHERE (flag & 1) != 0 AND (flag & (4 | 8 | 64 | 256 | 512 | 2048)) = 0
      AND (getvariable('insert_duplicates') OR (flag & 1024) = 0)
      AND tlen != 0 AND rname = rnext
), oriented AS (
    -- sam/pair_orientation.rs:30-59: equal strands TANDEM; reverse end
    -- is inclusive end_pos; forward mate reverse end is POS + TLEN.
    SELECT CASE
        WHEN ((flag & 16) != 0) = ((flag & 32) != 0) THEN 'TANDEM'
        WHEN (flag & 16) != 0 THEN CASE WHEN pnext < end_pos THEN 'FR' ELSE 'RF' END
        ELSE CASE WHEN pos < pos + tlen THEN 'FR' ELSE 'RF' END
    END AS pair_orientation, abs(tlen) AS insert_size
    FROM eligible
)
-- commands/isize.rs:211-218: count unsigned absolute TLEN by orientation.
SELECT pair_orientation, insert_size, count(*)::UBIGINT AS read_pairs
FROM oriented GROUP BY pair_orientation, insert_size;
