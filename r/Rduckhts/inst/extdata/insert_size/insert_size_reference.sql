-- Fulcrum Genomics riker 0.4.1 (MIT), tag 40d5e326, is the oracle.
-- src/commands/isize.rs:197-218: non-first, primary, mapped pair, no QC-fail,
-- optional duplicates, same reference and nonzero TLEN; absolute TLEN counts.
-- src/sam/pair_orientation.rs:30-59: equal strands TANDEM; compare inclusive
-- ends, with forward mate end POS + TLEN (not minus one); ties are RF.
-- Missing CIGAR spans one base (noodles alignment_end, called at :49).
WITH projected AS (
  SELECT FLAG, RNAME, RNEXT, POS::BIGINT AS pos, PNEXT::BIGINT AS pnext,
         TLEN::BIGINT AS tlen, CIGAR
  FROM read_bam(getvariable('insert_path'), reference := getvariable('insert_reference'),
                cigar_representation := 'binary')
), eligible AS (
  SELECT *, is_reverse_complemented(FLAG) AS rev,
         is_next_segment_reverse_complemented(FLAG) AS mate_rev
  FROM projected
  WHERE is_paired(FLAG) AND NOT is_unmapped(FLAG)
    AND NOT is_next_segment_unmapped(FLAG) AND NOT is_first_segment(FLAG)
    AND NOT is_secondary(FLAG) AND NOT is_supplementary(FLAG)
    AND NOT is_qc_fail(FLAG)
    AND (getvariable('insert_duplicates') OR NOT is_duplicate(FLAG))
    AND tlen <> 0 AND (RNEXT = '=' OR RNEXT = RNAME)
), oriented AS (
  SELECT CASE WHEN rev = mate_rev THEN 'TANDEM'
              WHEN rev THEN CASE WHEN pnext < pos + greatest(coalesce(cigar_reference_length(CIGAR), 0), 1) - 1
                                 THEN 'FR' ELSE 'RF' END
              ELSE CASE WHEN pos < pos + tlen THEN 'FR' ELSE 'RF' END
         END AS pair_orientation,
         abs(tlen) AS insert_size
  FROM eligible
)
SELECT pair_orientation, insert_size, count(*)::UBIGINT AS read_pairs
FROM oriented
GROUP BY pair_orientation, insert_size;
