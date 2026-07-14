DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 executable, using either its declared indexed
cache or a staged GFF oracle. Unresolved DuckVEP rows remain in the
denominator; they are not discarded as unsupported cases. The CSV is
append-only by source revision, corpus, and resident model. Independent
frozen distributions and seeds are kept separate so a fix cannot improve
its own hand-picked witnesses and hide a regression elsewhere.

## Latest tested revision

| revision | corpus                                | model                    | assembly       | species               | oracle_source | oracle    |  pairs | exact         | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:--------------------------------------|:-------------------------|:---------------|:----------------------|:--------------|:----------|-------:|:--------------|-----------:|-----------------------:|:------------------------|
| b204dd49 | final_clinvar_coding_seed113          | final-coding             | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 287859 | 287829/287859 |          4 |                     27 | 0.01%                   |
| b204dd49 | final_clinvar_crosschrom_seed17       | final-clinvar            | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 316399 | 316388/316399 |          2 |                     10 | 0.01%                   |
| b204dd49 | final_dbsnp157_windows_seed29         | final-dbsnp              | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  73620 | 73620/73620   |          0 |                      0 | 0.01%                   |
| b204dd49 | final_giab_grch38_seed71              | final-giab               | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  54905 | 54905/54905   |          0 |                      0 | 0.01%                   |
| b204dd49 | final_grch37_cache_seed37             | final-grch37             | GRCh37         | homo_sapiens          | cache         | VEP 116.0 | 486468 | 486332/486468 |        102 |                     80 | 0.02%                   |
| b204dd49 | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63 | GCA000002765v3 | plasmodium_falciparum | cache         | VEP 116.0 |  40734 | 40730/40734   |         24 |                      4 | 0.03%                   |

Detailed consequence tables below use the largest corpus at this tested
revision.

## Prepared model receipts

| revision | species               | release | assembly       | regions | transcripts | coding_backed | exons     | mature_miRNA_segments | peptide_edits | codon_tables     | model_sha256 |
|:---------|:----------------------|:--------|:---------------|:--------|:------------|:--------------|:----------|:----------------------|:--------------|:-----------------|:-------------|
| 8498b92a | homo_sapiens          | 116     | GRCh38         | 194     | 644,427     | 369,631       | 5,068,416 | 2,806                 | 389           | 1:369618;2:13    | 9c3404101b2b |
| 8498b92a | homo_sapiens          | 116     | GRCh37         | 84      | 195,379     | 94,610        | 1,186,433 | 3,788                 | 129           | 1:94597;2:13     | 25459e62e50d |
| 8498b92a | plasmodium_falciparum | 63/116  | GCA000002765v3 | 16      | 5,791       | 5,389         | 15,097    | 0                     | 4             | 1:5356;4:3;11:30 | c011cdd4deab |

These are complete model-build receipts, not counts inferred from a
differential. The ledger retains the full source-manifest, reference,
and model SHA-256 values, the exact VEP transcript filter, every count
above, CDS/flank base totals, and the external artifact name. The
Plasmodium row is an Ensembl Genomes release-63 cache paired with the
VEP/core-116 executable libraries, which is why both release numbers are
recorded.

## History

| run_date   | source_revision | corpus                                | model                    |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:--------------------------------------|:-------------------------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses                             | differential             |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29             | differential             | 100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                   |
| 2026-07-13 | 24bb1714        | witnesses                             | differential             |    242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                   |
| 2026-07-13 | 34b37ca1        | witnesses                             | differential             |    242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                   |
| 2026-07-13 | 87f03a2a        | witnesses                             | differential             |    242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | defc9a1c        | state_exploration_seed_113            | differential             | 100246 |       85598 |      28442 |      71804 |                1086 | 85.39%     | 1.60%                   |
| 2026-07-13 | defc9a1c        | state_exploration_seed_71             | differential             | 100242 |       85646 |      27946 |      72296 |                1103 | 85.44%     | 1.62%                   |
| 2026-07-13 | defc9a1c        | witnesses                             | differential             |    246 |         242 |         11 |        235 |                   0 | 98.37%     | 1.56%                   |
| 2026-07-13 | eb212de3        | witnesses                             | differential             |    242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_197            | differential             | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_211            | differential             | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_71             | differential             | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 2ab08e2f        | witnesses                             | differential             |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_113            | differential             | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_197            | differential             | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_211            | differential             | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_71             | differential             | 100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | witnesses                             | differential             |    268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                   |
| 2026-07-14 | 5e5bc1e2        | clinvar_chr21_seed1                   | ensembl116_grch38_core   | 126320 |      126320 |          0 |     126320 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 8498b92a        | final_clinvar_coding_seed113          | final-coding             | 287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_clinvar_crosschrom_seed17       | final-clinvar            | 316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_dbsnp157_windows_seed29         | final-dbsnp              |  73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_giab_grch38_seed71              | final-giab               |  54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_grch37_cache_seed37             | final-grch37             | 486468 |      482665 |        102 |     486366 |                3747 | 99.22%     | 0.80%                   |
| 2026-07-14 | 8498b92a        | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63 |  40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_197            | differential             | 100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_211            | differential             | 100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_71             | differential             | 100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 8b2a2dbc        | witnesses                             | differential             |    258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | b204dd49        | final_clinvar_coding_seed113          | final-coding             | 287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                   |
| 2026-07-14 | b204dd49        | final_clinvar_crosschrom_seed17       | final-clinvar            | 316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_dbsnp157_windows_seed29         | final-dbsnp              |  73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_giab_grch38_seed71              | final-giab               |  54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_grch37_cache_seed37             | final-grch37             | 486468 |      486332 |        102 |     486366 |                  80 | 99.97%     | 0.02%                   |
| 2026-07-14 | b204dd49        | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63 |  40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_197            | differential             | 100248 |       91159 |      19473 |      80775 |                 537 | 90.93%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_211            | differential             | 100250 |       90951 |      19676 |      80574 |                 567 | 90.72%     | 0.76%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_71             | differential             | 100242 |       90932 |      19500 |      80742 |                 532 | 90.71%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | witnesses                             | differential             |    262 |         258 |         10 |        252 |                   0 | 98.47%     | 1.45%                   |

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
| 2026-07-14 | b204dd49        | 0xd0c0ffee12345678 |                 40 | 3,900,500  | 3,900,500  |      0 |          0 |         154 | 204,654          |                20.979 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
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

| consequence_class                   | impact   |      n | exact_agree | unresolved | resolved_n | resolved_agree | resolved_discordant | term_mismatch | engine_extra | engine_missing | resolved_error_upper_95 |
|:------------------------------------|:---------|-------:|------------:|-----------:|-----------:|---------------:|--------------------:|--------------:|-------------:|---------------:|:------------------------|
| incomplete_terminal_codon_variant   | LOW      |     80 |           0 |         39 |         41 |              0 |                  41 |            80 |            0 |             80 | 100.00%                 |
| coding_sequence_variant             | MODIFIER |   4738 |        4701 |        102 |       4636 |           4607 |                  29 |            37 |            8 |             29 | 0.90%                   |
| frameshift_variant                  | HIGH     |  27248 |       27225 |          8 |      27240 |          27225 |                  15 |            23 |           15 |              8 | 0.09%                   |
| synonymous_variant                  | LOW      |   1908 |        1894 |          0 |       1908 |           1894 |                  14 |            14 |           14 |              0 | 1.23%                   |
| splice_polypyrimidine_tract_variant | LOW      |   2461 |        2449 |          0 |       2461 |           2449 |                  12 |            12 |           12 |              0 | 0.85%                   |
| inframe_insertion                   | MODERATE |      8 |           0 |          0 |          8 |              0 |                   8 |             8 |            0 |              8 | 100.00%                 |
| downstream_gene_variant             | MODIFIER |  93074 |       93070 |          0 |      93074 |          93070 |                   4 |             4 |            4 |              0 | 0.01%                   |
| start_lost                          | HIGH     |   1541 |        1537 |          0 |       1541 |           1537 |                   4 |             4 |            2 |              2 | 0.66%                   |
| missense_variant                    | MODERATE |  23313 |       23303 |          8 |      23305 |          23303 |                   2 |            10 |            2 |              8 | 0.03%                   |
| splice_acceptor_variant             | HIGH     |   8568 |        8566 |          0 |       8568 |           8566 |                   2 |             2 |            2 |              0 | 0.08%                   |
| stop_lost                           | HIGH     |    576 |         574 |          0 |        576 |            574 |                   2 |             2 |            0 |              2 | 1.25%                   |
| 5_prime_UTR_variant                 | MODIFIER |  10550 |       10549 |         36 |      10514 |          10513 |                   1 |             1 |            0 |              1 | 0.05%                   |
| splice_donor_region_variant         | LOW      |    160 |         159 |          0 |        160 |            159 |                   1 |             1 |            1 |              0 | 3.43%                   |
| splice_donor_5th_base_variant       | LOW      |     58 |          57 |          0 |         58 |             57 |                   1 |             1 |            1 |              0 | 9.24%                   |
| NMD_transcript_variant              | MODIFIER |  24184 |       24184 |         20 |      24164 |          24164 |                   0 |             0 |            0 |              0 | 0.02%                   |
| splice_region_variant               | LOW      |  46539 |       46539 |          1 |      46538 |          46538 |                   0 |             0 |            0 |              0 | 0.01%                   |
| stop_gained                         | HIGH     |   1257 |        1256 |          1 |       1256 |           1256 |                   0 |             1 |            0 |              1 | 0.29%                   |
| intron_variant                      | MODIFIER | 140741 |      140741 |          0 |     140741 |         140741 |                   0 |             0 |            0 |              0 | 0.00%                   |
| upstream_gene_variant               | MODIFIER |  80210 |       80210 |          0 |      80210 |          80210 |                   0 |             0 |            0 |              0 | 0.00%                   |
| non_coding_transcript_variant       | MODIFIER |  46096 |       46096 |          0 |      46096 |          46096 |                   0 |             0 |            0 |              0 | 0.01%                   |
| non_coding_transcript_exon_variant  | MODIFIER |  43896 |       43896 |          0 |      43896 |          43896 |                   0 |             0 |            0 |              0 | 0.01%                   |
| mature_miRNA_variant                | MODIFIER |  37942 |       37942 |          0 |      37942 |          37942 |                   0 |             0 |            0 |              0 | 0.01%                   |
| 3_prime_UTR_variant                 | MODIFIER |   9312 |        9312 |          0 |       9312 |           9312 |                   0 |             0 |            0 |              0 | 0.04%                   |
| splice_donor_variant                | HIGH     |   8082 |        8082 |          0 |       8082 |           8082 |                   0 |             0 |            0 |              0 | 0.05%                   |
| start_retained_variant              | LOW      |    274 |         274 |          0 |        274 |            274 |                   0 |             0 |            0 |              0 | 1.34%                   |
| stop_retained_variant               | LOW      |     38 |          38 |          0 |         38 |             38 |                   0 |             0 |            0 |              0 | 9.25%                   |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact            |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:------------------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| (no_vep_emission) |      4 |           0 |          0 |          4 |                   4 | 0.00%      | 100.00%                 |
| HIGH              |  46645 |       46632 |          9 |      46636 |                   4 | 99.97%     | 0.02%                   |
| LOW               |  23354 |       23279 |         39 |      23315 |                  36 | 99.68%     | 0.21%                   |
| MODERATE          |  23319 |       23303 |          8 |      23311 |                   8 | 99.93%     | 0.07%                   |
| MODIFIER          | 393146 |      393118 |         46 |     393100 |                  28 | 99.99%     | 0.01%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
