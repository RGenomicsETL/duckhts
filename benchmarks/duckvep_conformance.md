DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 `--gff` output. Unresolved DuckVEP rows remain
in the denominator; they are not discarded as unsupported cases. The CSV
is append-only by source revision, corpus, and resident model.

## Latest tested revision

| revision | corpus                    | model        | oracle    |  pairs | exact        | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:--------------------------|:-------------|:----------|-------:|:-------------|-----------:|-----------------------:|:------------------------|
| 24bb1714 | state_exploration_seed_29 | differential | VEP 116.0 | 100242 | 85238/100242 |      28109 |                   1453 | 2.12%                   |
| 24bb1714 | witnesses                 | differential | VEP 116.0 |    242 | 238/242      |         11 |                      0 | 1.58%                   |

Detailed consequence tables below use the largest corpus at this tested
revision.

## History

| run_date   | source_revision | corpus                    | model        |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:--------------------------|:-------------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses                 | differential |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29 | differential | 100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                   |
| 2026-07-13 | 24bb1714        | witnesses                 | differential |    242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                   |
| 2026-07-13 | 34b37ca1        | witnesses                 | differential |    242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                   |
| 2026-07-13 | 87f03a2a        | witnesses                 | differential |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | eb212de3        | witnesses                 | differential |    242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                   |

## Official Ensembl release corpus in Parquet

The official release consequence VCF is already BGZF-compressed. This
table measures its complete typed DuckHTS reader projection and the
narrower `VE` plus CSQ projection used by the bulk oracle lane. It is a
storage comparison, not a claim that the Parquet projection can
reproduce the original VCF byte-for-byte.

| revision | release | assembly | chromosome | projection  | columns | records    | ALT_alleles | CSQ_entries | source_MiB | parquet_MiB | parquet_of_source | elapsed_seconds | records_per_second |
|:---------|--------:|:---------|:-----------|:------------|--------:|:-----------|:------------|:------------|:-----------|:------------|:------------------|:----------------|:-------------------|
| 55c55238 |     116 | GRCh38   | 22         | full_typed  |      51 | 14,920,904 | 17,767,586  | 30,199,106  | 265.6      | 219.8       | 82.7%             | 55.2            | 270,179            |
| 55c55238 |     116 | GRCh38   | 22         | consequence |      14 | 14,920,904 | 17,767,586  | 30,199,106  | 265.6      | 155.9       | 58.7%             | 38.1            | 391,872            |

The ledger records the official source URL, SHA-256 of every input and
output, DuckHTS and DuckDB versions, compression, row-group size, thread
count, machine, and exact byte sizes.

## Randomized pure-C properties

The property ledger is separate from the VEP differential. It records
successful runs of each randomized oracle, including the exact seed and
duplicate count. A failed suite does not append rows.

| run_date   | source_revision | seed               | randomized_targets | trials     | passed     | failed | duplicates | suite_tests | suite_assertions | suite_elapsed_seconds | compiler                                   |
|:-----------|:----------------|:-------------------|-------------------:|:-----------|:-----------|-------:|-----------:|------------:|:-----------------|----------------------:|:-------------------------------------------|
| 2026-07-11 | 8cc22218        | 0xd0c0ffee12345678 |                 39 | 3,800,500  | 3,800,500  |      0 |          0 |         133 | 189,981          |                15.400 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-13 | 24bb1714        | 0x000000000000001d |                 40 | 39,000,500 | 39,000,500 |      0 |          0 |         142 | 1,873,864        |               187.865 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-13 | 34b37ca1        | 0xd0c0ffee12345678 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         140 | 190,041          |                18.113 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-13 | 87f03a2a        | 0xd0c0ffee12345678 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         139 | 190,024          |                18.541 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-13 | eb212de3        | 0xd0c0ffee12345678 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         141 | 190,064          |                23.425 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

| target                                                                 | trials  | passed  | failed | skipped | duplicates |
|:-----------------------------------------------------------------------|:--------|:--------|:-------|:--------|:-----------|
| allele-normalized sweep candidates == independent trim oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor cross-codon MNV route == tile                          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor DEL route == tile under output splits                  | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor INS route == tile under output splits                  | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor output splits == one annotate_tile                     | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor padded SNV == tile under output splits                 | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile == sweep + classify + structural-SO composition          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon refinement == coding-SNV kernel oracle             | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon-aligned in-frame deletion == CDS-position oracle   | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon-boundary in-frame insertion == CDS-position oracle | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile non-boundary in-frame insertion == peptide-window oracle | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile rejects NULL model without reading the batch             | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile same-codon MNV == codon oracle                           | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile simple frameshift indel == CDS-position oracle           | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile start_lost SNV == start-codon oracle                     | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile two-codon body MNV missense == codon-window oracle       | 100,000 | 100,000 | 0      | 0       | 0          |
| cgranges-seeded first event + sweep == brute-force candidates          | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context == direct CDS splice + full peptide oracles             | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta == single-codon oracle                            | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame deletion == edit-origin oracle           | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame delins == edit-origin oracle             | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame insertion == edit-origin oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| codon change classification consistent with translation                | 100,000 | 100,000 | 0      | 0       | 0          |
| coordinate projection == brute-force transcript-order base walk        | 100,000 | 100,000 | 0      | 0       | 0          |
| event differing-region normalization == independent trim oracle        | 100,000 | 100,000 | 0      | 0       | 0          |
| haplotype blocks preserve every frame and same-codon interaction       | 100,000 | 100,000 | 0      | 0       | 0          |
| multi-edit CDS haplotype apply == left-to-right rebuild oracle         | 100,000 | 100,000 | 0      | 0       | 0          |
| region mask structural invariants                                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta annotation wrapper MNV == direct shape                  | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch INDEL == in-frame delins oracle                 | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch MNV == single-codon oracle                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch two-codon MNV window == codon-window oracle     | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence-backed SNV codon edit == codon-slice edit oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
| sorted point cursor classifier == exhaustive exon/gap scans            | 100,000 | 100,000 | 0      | 0       | 0          |
| sweep candidate set == brute-force candidate set                       | 100,000 | 100,000 | 0      | 0       | 0          |
| tile_controller_preserves_sorted_stream                                | 500     | 500     | 0      | 0       | 0          |
| variant CDS edit builder == direct CDS splice oracle                   | 100,000 | 100,000 | 0      | 0       | 0          |
| variant CDS edit-set builder == single-edit splice oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
| variant CDS edit-set builder splits MNV diff islands                   | 100,000 | 100,000 | 0      | 0       | 0          |
| variant coding context == direct CDS splice + full peptide oracles     | 100,000 | 100,000 | 0      | 0       | 0          |

## Individual Sequence Ontology terms

For each transcript pair, this compares the union of terms emitted by
either engine. A missing or extra term is therefore visible under its
own SO name. Rows must not be summed across terms because one pair can
carry several terms.

| consequence_class                   | impact   |     n | exact_agree | unresolved | resolved_n | resolved_agree | resolved_discordant | term_mismatch | engine_extra | engine_missing | resolved_error_upper_95 |
|:------------------------------------|:---------|------:|------------:|-----------:|-----------:|---------------:|--------------------:|--------------:|-------------:|---------------:|:------------------------|
| inframe_insertion                   | MODERATE |  4694 |        4124 |         65 |       4629 |           4124 |                 505 |           570 |          249 |            321 | 11.84%                  |
| stop_retained_variant               | LOW      |  2630 |          38 |       2296 |        334 |             38 |                 296 |          2592 |          172 |           2420 | 91.82%                  |
| coding_sequence_variant             | MODIFIER | 38860 |       25201 |      28109 |      10751 |          10514 |                 237 |         13659 |        13422 |            237 | 2.50%                   |
| start_retained_variant              | LOW      |   343 |          91 |         89 |        254 |             91 |                 163 |           252 |            0 |            252 | 70.07%                  |
| intron_variant                      | MODIFIER | 38610 |       38354 |      10525 |      28085 |          27940 |                 145 |           256 |           41 |            215 | 0.61%                   |
| missense_variant                    | MODERATE |  4392 |        4252 |          0 |       4392 |           4252 |                 140 |           140 |          126 |             14 | 3.75%                   |
| splice_region_variant               | LOW      |  8192 |        7996 |        709 |       7483 |           7344 |                 139 |           196 |           67 |            129 | 2.19%                   |
| protein_altering_variant            | MODERATE |  2584 |         941 |       1510 |       1074 |            941 |                 133 |          1643 |            0 |           1643 | 14.50%                  |
| splice_acceptor_variant             | HIGH     | 11687 |       11584 |       5286 |       6401 |           6303 |                  98 |           103 |          100 |              3 | 1.86%                   |
| splice_polypyrimidine_tract_variant | LOW      | 10426 |       10247 |        104 |      10322 |          10247 |                  75 |           179 |           72 |            107 | 0.91%                   |
| start_lost                          | HIGH     |  8244 |        1625 |       6546 |       1698 |           1625 |                  73 |          6619 |            0 |           6619 | 5.38%                   |
| stop_gained                         | HIGH     |  8707 |        7986 |        661 |       8046 |           7986 |                  60 |           721 |           60 |            661 | 0.96%                   |
| stop_lost                           | HIGH     |  3735 |        1016 |       2668 |       1067 |           1016 |                  51 |          2719 |           51 |           2668 | 6.24%                   |
| splice_donor_5th_base_variant       | LOW      | 11101 |       10978 |       4927 |       6174 |           6126 |                  48 |           123 |          105 |             18 | 1.03%                   |
| 5_prime_UTR_variant                 | MODIFIER | 15185 |       15077 |       7287 |       7898 |           7852 |                  46 |           108 |          108 |              0 | 0.78%                   |
| 3_prime_UTR_variant                 | MODIFIER | 13992 |       13888 |       8255 |       5737 |           5692 |                  45 |           104 |          104 |              0 | 1.05%                   |
| splice_donor_variant                | HIGH     | 12706 |       12614 |       5976 |       6730 |           6685 |                  45 |            92 |           90 |              2 | 0.89%                   |
| inframe_deletion                    | MODERATE |   452 |         371 |         43 |        409 |            371 |                  38 |            81 |           38 |             43 | 12.53%                  |
| splice_donor_region_variant         | LOW      |  4765 |        4655 |       1047 |       3718 |           3693 |                  25 |           110 |           54 |             56 | 0.99%                   |
| downstream_gene_variant             | MODIFIER |   409 |         391 |          0 |        409 |            391 |                  18 |            18 |            0 |             18 | 6.87%                   |
| frameshift_variant                  | HIGH     | 16257 |       16039 |        218 |      16039 |          16039 |                   0 |           218 |            0 |            218 | 0.02%                   |
| intergenic_variant                  | MODIFIER |    15 |           0 |         15 |          0 |              0 |                   0 |            15 |            0 |             15 | not estimable           |
| synonymous_variant                  | LOW      |   103 |         103 |          0 |        103 |            103 |                   0 |             0 |            0 |              0 | 3.52%                   |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |     n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 54268 |       43042 |      20947 |      33321 |                 720 | 79.31%     | 2.32%                   |
| LOW      | 14990 |       12762 |       2237 |      12753 |                 212 | 85.14%     | 1.90%                   |
| MODERATE |  9423 |        8090 |        947 |       8476 |                 386 | 85.85%     | 5.02%                   |
| MODIFIER | 21561 |       21344 |       3978 |      17583 |                 135 | 98.99%     | 0.91%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
