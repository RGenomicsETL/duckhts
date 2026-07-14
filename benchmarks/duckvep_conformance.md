DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 `--gff` output. Unresolved DuckVEP rows remain
in the denominator; they are not discarded as unsupported cases. The CSV
is append-only by source revision, corpus, and resident model.

## Latest tested revision

| revision | corpus                     | model        | oracle    |  pairs | exact         | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:---------------------------|:-------------|:----------|-------:|:--------------|-----------:|-----------------------:|:------------------------|
| 3c427df4 | state_exploration_seed_113 | differential | VEP 116.0 | 100268 | 100268/100268 |          0 |                      0 | 0.00%                   |
| 3c427df4 | state_exploration_seed_197 | differential | VEP 116.0 | 100268 | 100268/100268 |          0 |                      0 | 0.00%                   |
| 3c427df4 | state_exploration_seed_211 | differential | VEP 116.0 | 100268 | 100268/100268 |          0 |                      0 | 0.00%                   |
| 3c427df4 | state_exploration_seed_71  | differential | VEP 116.0 | 100268 | 100268/100268 |          0 |                      0 | 0.00%                   |
| 3c427df4 | witnesses                  | differential | VEP 116.0 |    268 | 268/268       |          0 |                      0 | 1.37%                   |

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
| 2026-07-14 | 3c427df4        | state_exploration_seed_113 | differential | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_197 | differential | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_211 | differential | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_71  | differential | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | witnesses                  | differential |    268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_197 | differential | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_211 | differential | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_71  | differential | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 8b2a2dbc        | witnesses                  | differential |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_197 | differential | 100248 |       91159 |      19473 |      80775 |                 537 | 90.93%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_211 | differential | 100250 |       90951 |      19676 |      80574 |                 567 | 90.72%     | 0.76%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_71  | differential | 100242 |       90932 |      19500 |      80742 |                 532 | 90.71%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | witnesses                  | differential |    262 |         258 |         10 |        252 |                   0 | 98.47%     | 1.45%                   |

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
| intron_variant                      | MODIFIER | 38757 |       38757 |          0 |      38757 |          38757 |                   0 |             0 |            0 |              0 | 0.01%                   |
| coding_sequence_variant             | MODIFIER | 25833 |       25833 |          0 |      25833 |          25833 |                   0 |             0 |            0 |              0 | 0.01%                   |
| frameshift_variant                  | HIGH     | 16210 |       16210 |          0 |      16210 |          16210 |                   0 |             0 |            0 |              0 | 0.02%                   |
| 5_prime_UTR_variant                 | MODIFIER | 14904 |       14904 |          0 |      14904 |          14904 |                   0 |             0 |            0 |              0 | 0.02%                   |
| 3_prime_UTR_variant                 | MODIFIER | 13953 |       13953 |          0 |      13953 |          13953 |                   0 |             0 |            0 |              0 | 0.03%                   |
| splice_donor_variant                | HIGH     | 12821 |       12821 |          0 |      12821 |          12821 |                   0 |             0 |            0 |              0 | 0.03%                   |
| splice_acceptor_variant             | HIGH     | 11618 |       11618 |          0 |      11618 |          11618 |                   0 |             0 |            0 |              0 | 0.03%                   |
| splice_donor_5th_base_variant       | LOW      | 11077 |       11077 |          0 |      11077 |          11077 |                   0 |             0 |            0 |              0 | 0.03%                   |
| splice_polypyrimidine_tract_variant | LOW      | 10415 |       10415 |          0 |      10415 |          10415 |                   0 |             0 |            0 |              0 | 0.04%                   |
| stop_gained                         | HIGH     |  8551 |        8551 |          0 |       8551 |           8551 |                   0 |             0 |            0 |              0 | 0.04%                   |
| start_lost                          | HIGH     |  8145 |        8145 |          0 |       8145 |           8145 |                   0 |             0 |            0 |              0 | 0.05%                   |
| splice_region_variant               | LOW      |  8071 |        8071 |          0 |       8071 |           8071 |                   0 |             0 |            0 |              0 | 0.05%                   |
| splice_donor_region_variant         | LOW      |  4713 |        4713 |          0 |       4713 |           4713 |                   0 |             0 |            0 |              0 | 0.08%                   |
| missense_variant                    | MODERATE |  4354 |        4354 |          0 |       4354 |           4354 |                   0 |             0 |            0 |              0 | 0.08%                   |
| inframe_insertion                   | MODERATE |  4319 |        4319 |          0 |       4319 |           4319 |                   0 |             0 |            0 |              0 | 0.09%                   |
| stop_lost                           | HIGH     |  3755 |        3755 |          0 |       3755 |           3755 |                   0 |             0 |            0 |              0 | 0.10%                   |
| protein_altering_variant            | MODERATE |  2633 |        2633 |          0 |       2633 |           2633 |                   0 |             0 |            0 |              0 | 0.14%                   |
| stop_retained_variant               | LOW      |  2469 |        2469 |          0 |       2469 |           2469 |                   0 |             0 |            0 |              0 | 0.15%                   |
| inframe_deletion                    | MODERATE |   420 |         420 |          0 |        420 |            420 |                   0 |             0 |            0 |              0 | 0.87%                   |
| downstream_gene_variant             | MODIFIER |   396 |         396 |          0 |        396 |            396 |                   0 |             0 |            0 |              0 | 0.93%                   |
| start_retained_variant              | LOW      |   322 |         322 |          0 |        322 |            322 |                   0 |             0 |            0 |              0 | 1.14%                   |
| synonymous_variant                  | LOW      |   106 |         106 |          0 |        106 |            106 |                   0 |             0 |            0 |              0 | 3.42%                   |
| intergenic_variant                  | MODIFIER |    18 |          18 |          0 |         18 |             18 |                   0 |             0 |            0 |              0 | 18.53%                  |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |     n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 54344 |       54344 |          0 |      54344 |                   0 | 100.00%    | 0.01%                   |
| LOW      | 14877 |       14877 |          0 |      14877 |                   0 | 100.00%    | 0.02%                   |
| MODERATE |  9521 |        9521 |          0 |       9521 |                   0 | 100.00%    | 0.04%                   |
| MODIFIER | 21526 |       21526 |          0 |      21526 |                   0 | 100.00%    | 0.02%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
