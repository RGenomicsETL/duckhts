DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 `--gff` output. Unresolved DuckVEP rows remain
in the denominator; they are not discarded as unsupported cases. The CSV
is append-only by source revision, corpus, and resident model.

## Latest tested revision

| revision | corpus                     | model        | oracle    |  pairs | exact        | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:---------------------------|:-------------|:----------|-------:|:-------------|-----------:|-----------------------:|:------------------------|
| 2ab08e2f | state_exploration_seed_197 | differential | VEP 116.0 | 100248 | 88021/100248 |      22598 |                    550 | 0.77%                   |
| 2ab08e2f | state_exploration_seed_211 | differential | VEP 116.0 | 100250 | 87815/100250 |      22801 |                    578 | 0.81%                   |
| 2ab08e2f | state_exploration_seed_71  | differential | VEP 116.0 | 100242 | 87916/100242 |      22502 |                    546 | 0.76%                   |
| 2ab08e2f | witnesses                  | differential | VEP 116.0 |    258 | 254/258      |         10 |                      0 | 1.48%                   |

Detailed consequence tables below use the largest corpus at this tested
revision.

## History

| run_date   | source_revision | corpus                     | model        |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:---------------------------|:-------------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses                  | differential |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29  | differential | 100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                   |
| 2026-07-13 | 24bb1714        | witnesses                  | differential |    242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                   |
| 2026-07-13 | 34b37ca1        | witnesses                  | differential |    242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                   |
| 2026-07-13 | 87f03a2a        | witnesses                  | differential |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | defc9a1c        | state_exploration_seed_113 | differential | 100246 |       85598 |      28442 |      71804 |                1086 | 85.39%     | 1.60%                   |
| 2026-07-13 | defc9a1c        | state_exploration_seed_71  | differential | 100242 |       85646 |      27946 |      72296 |                1103 | 85.44%     | 1.62%                   |
| 2026-07-13 | defc9a1c        | witnesses                  | differential |    246 |         242 |         11 |        235 |                   0 | 98.37%     | 1.56%                   |
| 2026-07-13 | eb212de3        | witnesses                  | differential |    242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_197 | differential | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_211 | differential | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_71  | differential | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 2ab08e2f        | witnesses                  | differential |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |

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
| 2026-07-13 | defc9a1c        | 0x0000000000000071 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         142 | 189,962          |                18.793 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-13 | eb212de3        | 0xd0c0ffee12345678 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         141 | 190,064          |                23.425 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-14 | 2ab08e2f        | 0x0000000000000139 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         144 | 201,583          |                18.940 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

| target                                                                 | trials  | passed  | failed | skipped | duplicates |
|:-----------------------------------------------------------------------|:--------|:--------|:-------|:--------|:-----------|
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
| coding context delins shape == local-edge oracle                       | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta == single-codon oracle                            | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame deletion == edit-origin oracle           | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame insertion == edit-origin oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| codon change classification consistent with translation                | 100,000 | 100,000 | 0      | 0       | 0          |
| coordinate projection == brute-force transcript-order base walk        | 100,000 | 100,000 | 0      | 0       | 0          |
| event differing-region normalization == independent trim oracle        | 100,000 | 100,000 | 0      | 0       | 0          |
| haplotype blocks preserve every frame and same-codon interaction       | 100,000 | 100,000 | 0      | 0       | 0          |
| multi-edit CDS haplotype apply == left-to-right rebuild oracle         | 100,000 | 100,000 | 0      | 0       | 0          |
| region mask structural invariants                                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta annotation wrapper MNV == direct shape                  | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch INDEL == local delins-shape oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
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
| VEP feature-span sweep candidates == independent parser oracle         | 100,000 | 100,000 | 0      | 0       | 0          |

## Individual Sequence Ontology terms

For each transcript pair, this compares the union of terms emitted by
either engine. A missing or extra term is therefore visible under its
own SO name. Rows must not be summed across terms because one pair can
carry several terms.

| consequence_class                   | impact   |     n | exact_agree | unresolved | resolved_n | resolved_agree | resolved_discordant | term_mismatch | engine_extra | engine_missing | resolved_error_upper_95 |
|:------------------------------------|:---------|------:|------------:|-----------:|-----------:|---------------:|--------------------:|--------------:|-------------:|---------------:|:------------------------|
| splice_region_variant               | LOW      |  8120 |        7916 |        676 |       7444 |           7287 |                 157 |           204 |           78 |            126 | 2.46%                   |
| splice_acceptor_variant             | HIGH     | 11580 |       11456 |       5351 |       6229 |           6119 |                 110 |           124 |          118 |              6 | 2.12%                   |
| start_retained_variant              | LOW      |   333 |         168 |         75 |        258 |            168 |                  90 |           165 |            0 |            165 | 41.04%                  |
| splice_polypyrimidine_tract_variant | LOW      | 10384 |       10189 |        121 |      10263 |          10189 |                  74 |           195 |           70 |            125 | 0.90%                   |
| 5_prime_UTR_variant                 | MODIFIER | 15125 |       15014 |       7188 |       7937 |           7867 |                  70 |           111 |          111 |              0 | 1.11%                   |
| splice_donor_5th_base_variant       | LOW      | 11330 |       11207 |       4964 |       6366 |           6309 |                  57 |           123 |          104 |             19 | 1.16%                   |
| splice_donor_variant                | HIGH     | 12888 |       12821 |       5950 |       6938 |           6883 |                  55 |            67 |           66 |              1 | 1.03%                   |
| 3_prime_UTR_variant                 | MODIFIER | 14023 |       13925 |       4533 |       9490 |           9435 |                  55 |            98 |           98 |              0 | 0.75%                   |
| intron_variant                      | MODIFIER | 38632 |       38540 |      10673 |      27959 |          27913 |                  46 |            92 |           86 |              6 | 0.22%                   |
| splice_donor_region_variant         | LOW      |  4672 |        4547 |       1025 |       3647 |           3613 |                  34 |           125 |           65 |             60 | 1.30%                   |
| downstream_gene_variant             | MODIFIER |   438 |         423 |          0 |        438 |            423 |                  15 |            15 |            0 |             15 | 5.59%                   |
| stop_gained                         | HIGH     |  8717 |        8612 |         93 |       8624 |           8612 |                  12 |           105 |           12 |             93 | 0.24%                   |
| protein_altering_variant            | MODERATE |  2728 |        2589 |        129 |       2599 |           2589 |                  10 |           139 |            0 |            139 | 0.71%                   |
| intergenic_variant                  | MODIFIER |    14 |           0 |          5 |          9 |              0 |                   9 |            14 |            0 |             14 | 100.00%                 |
| inframe_insertion                   | MODERATE |  4216 |        4207 |          0 |       4216 |           4207 |                   9 |             9 |            0 |              9 | 0.40%                   |
| start_lost                          | HIGH     |  8253 |        1892 |       6357 |       1896 |           1892 |                   4 |          6361 |            0 |           6361 | 0.54%                   |
| stop_lost                           | HIGH     |  3730 |        1079 |       2649 |       1081 |           1079 |                   2 |          2651 |            2 |           2649 | 0.67%                   |
| stop_retained_variant               | LOW      |  2428 |         140 |       2286 |        142 |            140 |                   2 |          2288 |            0 |           2288 | 5.00%                   |
| missense_variant                    | MODERATE |  4201 |        4200 |          0 |       4201 |           4200 |                   1 |             1 |            1 |              0 | 0.13%                   |
| coding_sequence_variant             | MODIFIER | 37366 |       25650 |      22801 |      14565 |          14565 |                   0 |         11716 |        11716 |              0 | 0.03%                   |
| frameshift_variant                  | HIGH     | 16455 |       16256 |        199 |      16256 |          16256 |                   0 |           199 |            0 |            199 | 0.02%                   |
| inframe_deletion                    | MODERATE |   408 |         367 |         41 |        367 |            367 |                   0 |            41 |            0 |             41 | 1.00%                   |
| synonymous_variant                  | LOW      |   113 |         113 |          0 |        113 |            113 |                   0 |             0 |            0 |              0 | 3.21%                   |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |     n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 54546 |       44489 |      20328 |      34218 |                 306 | 81.56%     | 1.00%                   |
| LOW      | 14822 |       12637 |       2251 |      12571 |                 177 | 85.26%     | 1.63%                   |
| MODERATE |  9345 |        9300 |         12 |       9333 |                  33 | 99.52%     | 0.50%                   |
| MODIFIER | 21537 |       21389 |        210 |      21327 |                  62 | 99.31%     | 0.37%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
