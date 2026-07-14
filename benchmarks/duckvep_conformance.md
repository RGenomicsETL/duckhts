DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 `--gff` output. Unresolved DuckVEP rows remain
in the denominator; they are not discarded as unsupported cases. The CSV
is append-only by source revision, corpus, and resident model.

## Latest tested revision

| revision | corpus              | model                  | oracle    |  pairs | exact         | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:--------------------|:-----------------------|:----------|-------:|:--------------|-----------:|-----------------------:|:------------------------|
| 5e5bc1e2 | clinvar_chr21_seed1 | ensembl116_grch38_core | VEP 116.0 | 126320 | 126320/126320 |          0 |                      0 | 0.00%                   |

Detailed consequence tables below use the largest corpus at this tested
revision.

## History

| run_date   | source_revision | corpus                     | model                  |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:---------------------------|:-----------------------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses                  | differential           |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29  | differential           | 100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                   |
| 2026-07-13 | 24bb1714        | witnesses                  | differential           |    242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                   |
| 2026-07-13 | 34b37ca1        | witnesses                  | differential           |    242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                   |
| 2026-07-13 | 87f03a2a        | witnesses                  | differential           |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | defc9a1c        | state_exploration_seed_113 | differential           | 100246 |       85598 |      28442 |      71804 |                1086 | 85.39%     | 1.60%                   |
| 2026-07-13 | defc9a1c        | state_exploration_seed_71  | differential           | 100242 |       85646 |      27946 |      72296 |                1103 | 85.44%     | 1.62%                   |
| 2026-07-13 | defc9a1c        | witnesses                  | differential           |    246 |         242 |         11 |        235 |                   0 | 98.37%     | 1.56%                   |
| 2026-07-13 | eb212de3        | witnesses                  | differential           |    242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_197 | differential           | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_211 | differential           | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_71  | differential           | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 2ab08e2f        | witnesses                  | differential           |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_113 | differential           | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_197 | differential           | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_211 | differential           | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_71  | differential           | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | witnesses                  | differential           |    268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                   |
| 2026-07-14 | 5e5bc1e2        | clinvar_chr21_seed1        | ensembl116_grch38_core | 126320 |      126320 |          0 |     126320 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_197 | differential           | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_211 | differential           | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_71  | differential           | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 8b2a2dbc        | witnesses                  | differential           |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_197 | differential           | 100248 |       91159 |      19473 |      80775 |                 537 | 90.93%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_211 | differential           | 100250 |       90951 |      19676 |      80574 |                 567 | 90.72%     | 0.76%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_71  | differential           | 100242 |       90932 |      19500 |      80742 |                 532 | 90.71%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | witnesses                  | differential           |    262 |         258 |         10 |        252 |                   0 | 98.47%     | 1.45%                   |

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
| 2026-07-14 | 3c427df4        | 0x0000000000000139 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         149 | 202,513          |                19.027 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-14 | 8b2a2dbc        | 0x0000000000000139 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         145 | 201,595          |                21.217 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-14 | fe6f0634        | 0x0000000000000139 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         146 | 201,660          |                18.570 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

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
| intron_variant                      | MODIFIER | 35998 |       35998 |          0 |      35998 |          35998 |                   0 |             0 |            0 |              0 | 0.01%                   |
| frameshift_variant                  | HIGH     | 25230 |       25230 |          0 |      25230 |          25230 |                   0 |             0 |            0 |              0 | 0.01%                   |
| NMD_transcript_variant              | MODIFIER | 21153 |       21153 |          0 |      21153 |          21153 |                   0 |             0 |            0 |              0 | 0.02%                   |
| downstream_gene_variant             | MODIFIER | 17123 |       17123 |          0 |      17123 |          17123 |                   0 |             0 |            0 |              0 | 0.02%                   |
| missense_variant                    | MODERATE | 10334 |       10334 |          0 |      10334 |          10334 |                   0 |             0 |            0 |              0 | 0.04%                   |
| 3_prime_UTR_variant                 | MODIFIER |  9804 |        9804 |          0 |       9804 |           9804 |                   0 |             0 |            0 |              0 | 0.04%                   |
| splice_polypyrimidine_tract_variant | LOW      |  8422 |        8422 |          0 |       8422 |           8422 |                   0 |             0 |            0 |              0 | 0.04%                   |
| upstream_gene_variant               | MODIFIER |  8197 |        8197 |          0 |       8197 |           8197 |                   0 |             0 |            0 |              0 | 0.04%                   |
| splice_region_variant               | LOW      |  6191 |        6191 |          0 |       6191 |           6191 |                   0 |             0 |            0 |              0 | 0.06%                   |
| inframe_deletion                    | MODERATE |  6052 |        6052 |          0 |       6052 |           6052 |                   0 |             0 |            0 |              0 | 0.06%                   |
| non_coding_transcript_exon_variant  | MODIFIER |  3848 |        3848 |          0 |       3848 |           3848 |                   0 |             0 |            0 |              0 | 0.10%                   |
| synonymous_variant                  | LOW      |  3655 |        3655 |          0 |       3655 |           3655 |                   0 |             0 |            0 |              0 | 0.10%                   |
| non_coding_transcript_variant       | MODIFIER |  3109 |        3109 |          0 |       3109 |           3109 |                   0 |             0 |            0 |              0 | 0.12%                   |
| inframe_insertion                   | MODERATE |  2444 |        2444 |          0 |       2444 |           2444 |                   0 |             0 |            0 |              0 | 0.15%                   |
| splice_acceptor_variant             | HIGH     |  2249 |        2249 |          0 |       2249 |           2249 |                   0 |             0 |            0 |              0 | 0.16%                   |
| splice_donor_variant                | HIGH     |  2175 |        2175 |          0 |       2175 |           2175 |                   0 |             0 |            0 |              0 | 0.17%                   |
| coding_sequence_variant             | MODIFIER |  1248 |        1248 |          0 |       1248 |           1248 |                   0 |             0 |            0 |              0 | 0.30%                   |
| splice_donor_5th_base_variant       | LOW      |  1176 |        1176 |          0 |       1176 |           1176 |                   0 |             0 |            0 |              0 | 0.31%                   |
| splice_donor_region_variant         | LOW      |  1161 |        1161 |          0 |       1161 |           1161 |                   0 |             0 |            0 |              0 | 0.32%                   |
| stop_gained                         | HIGH     |   808 |         808 |          0 |        808 |            808 |                   0 |             0 |            0 |              0 | 0.46%                   |
| 5_prime_UTR_variant                 | MODIFIER |   709 |         709 |          0 |        709 |            709 |                   0 |             0 |            0 |              0 | 0.52%                   |
| protein_altering_variant            | MODERATE |   207 |         207 |          0 |        207 |            207 |                   0 |             0 |            0 |              0 | 1.77%                   |
| stop_retained_variant               | LOW      |    93 |          93 |          0 |         93 |             93 |                   0 |             0 |            0 |              0 | 3.89%                   |
| start_lost                          | HIGH     |    73 |          73 |          0 |         73 |             73 |                   0 |             0 |            0 |              0 | 4.93%                   |
| start_retained_variant              | LOW      |    69 |          69 |          0 |         69 |             69 |                   0 |             0 |            0 |              0 | 5.21%                   |
| stop_lost                           | HIGH     |    37 |          37 |          0 |         37 |             37 |                   0 |             0 |            0 |              0 | 9.49%                   |
| mature_miRNA_variant                | MODIFIER |     2 |           2 |          0 |          2 |              2 |                   0 |             0 |            0 |              0 | 84.19%                  |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |     n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 30169 |       30169 |          0 |      30169 |                   0 | 100.00%    | 0.01%                   |
| LOW      | 14689 |       14689 |          0 |      14689 |                   0 | 100.00%    | 0.03%                   |
| MODERATE | 18986 |       18986 |          0 |      18986 |                   0 | 100.00%    | 0.02%                   |
| MODIFIER | 62476 |       62476 |          0 |      62476 |                   0 | 100.00%    | 0.01%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
