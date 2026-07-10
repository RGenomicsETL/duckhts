-- Apply after test/duckvep_model_valid.sql.
UPDATE duckvep_model.model_transcript_exon
SET exon_rank = 3
WHERE transcript_key = 1 AND exon_key = 2;
