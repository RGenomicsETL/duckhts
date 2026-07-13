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
| defc9a1c | state_exploration_seed_113 | differential | VEP 116.0 | 100246 | 85598/100246 |      28442 |                   1086 | 1.60%                   |
| defc9a1c | state_exploration_seed_71  | differential | VEP 116.0 | 100242 | 85646/100242 |      27946 |                   1103 | 1.62%                   |
| defc9a1c | witnesses                  | differential | VEP 116.0 |    246 | 242/246      |         11 |                      0 | 1.56%                   |

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
| inframe_insertion                   | MODERATE |  4546 |        4264 |         41 |       4505 |           4264 |                 241 |           282 |          228 |             54 | 6.05%                   |
| splice_region_variant               | LOW      |  8154 |        7937 |        764 |       7390 |           7236 |                 154 |           217 |           84 |            133 | 2.44%                   |
| start_retained_variant              | LOW      |   320 |          93 |         85 |        235 |             93 |                 142 |           227 |            0 |            227 | 66.72%                  |
| missense_variant                    | MODERATE |  4465 |        4331 |          0 |       4465 |           4331 |                 134 |           134 |          111 |             23 | 3.54%                   |
| protein_altering_variant            | MODERATE |  2631 |         912 |       1590 |       1041 |            912 |                 129 |          1719 |            0 |           1719 | 14.55%                  |
| intron_variant                      | MODIFIER | 38790 |       38550 |      10759 |      28031 |          27913 |                 118 |           240 |           40 |            200 | 0.50%                   |
| splice_acceptor_variant             | HIGH     | 11728 |       11609 |       5359 |       6369 |           6259 |                 110 |           119 |          114 |              5 | 2.08%                   |
| coding_sequence_variant             | MODIFIER | 39275 |       25750 |      28442 |      10833 |          10755 |                  78 |         13525 |        13447 |             78 | 0.90%                   |
| splice_polypyrimidine_tract_variant | LOW      | 10467 |       10303 |        105 |      10362 |          10303 |                  59 |           164 |           53 |            111 | 0.73%                   |
| splice_donor_variant                | HIGH     | 12923 |       12818 |       6067 |       6856 |           6799 |                  57 |           105 |          104 |              1 | 1.08%                   |
| stop_retained_variant               | LOW      |  2471 |         122 |       2295 |        176 |            122 |                  54 |          2349 |            3 |           2346 | 38.06%                  |
| 5_prime_UTR_variant                 | MODIFIER | 15004 |       14898 |       7186 |       7818 |           7765 |                  53 |           106 |          106 |              0 | 0.89%                   |
| start_lost                          | HIGH     |  8138 |        1632 |       6453 |       1685 |           1632 |                  53 |          6506 |            0 |           6506 | 4.09%                   |
| stop_lost                           | HIGH     |  3808 |        1004 |       2751 |       1057 |           1004 |                  53 |          2804 |           53 |           2751 | 6.51%                   |
| splice_donor_5th_base_variant       | LOW      | 11179 |       11056 |       5062 |       6117 |           6066 |                  51 |           123 |          104 |             19 | 1.09%                   |
| 3_prime_UTR_variant                 | MODIFIER | 14055 |       13952 |       8451 |       5604 |           5559 |                  45 |           103 |          103 |              0 | 1.07%                   |
| inframe_deletion                    | MODERATE |   458 |         377 |         43 |        415 |            377 |                  38 |            81 |           38 |             43 | 12.35%                  |
| stop_gained                         | HIGH     |  8583 |        7909 |        638 |       7945 |           7909 |                  36 |           674 |           36 |            638 | 0.63%                   |
| splice_donor_region_variant         | LOW      |  4771 |        4661 |       1018 |       3753 |           3726 |                  27 |           110 |           58 |             52 | 1.05%                   |
| downstream_gene_variant             | MODIFIER |   396 |         378 |          0 |        396 |            378 |                  18 |            18 |            0 |             18 | 7.09%                   |
| frameshift_variant                  | HIGH     | 16206 |       15986 |        220 |      15986 |          15986 |                   0 |           220 |            0 |            220 | 0.02%                   |
| intergenic_variant                  | MODIFIER |    18 |           0 |         18 |          0 |              0 |                   0 |            18 |            0 |             18 | not estimable           |
| synonymous_variant                  | LOW      |   106 |         106 |          0 |        106 |            106 |                   0 |             0 |            0 |              0 | 3.42%                   |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |     n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 54327 |       43311 |      21069 |      33258 |                 573 | 79.72%     | 1.87%                   |
| LOW      | 14874 |       12629 |       2274 |      12600 |                 221 | 84.91%     | 2.00%                   |
| MODERATE |  9520 |        8358 |       1000 |       8520 |                 162 | 87.79%     | 2.21%                   |
| MODIFIER | 21525 |       21300 |       4099 |      17426 |                 130 | 98.95%     | 0.89%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
