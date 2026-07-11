DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact transcript-level consequence agreement with
the real Ensembl VEP 116 `--gff` output. Unresolved DuckVEP rows remain
in the denominator; they are not discarded as unsupported cases. The CSV
is append-only by source revision, corpus, and resident model.

## Latest run

| revision | corpus    | model        | oracle    | pairs | exact   | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:----------|:-------------|:----------|------:|:--------|-----------:|-----------------------:|:------------------------|
| 8cc22218 | witnesses | differential | VEP 116.0 |   242 | 203/242 |         33 |                     15 | 11.56%                  |

## History

| run_date   | source_revision | corpus    | model        |   n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:----------|:-------------|----:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses | differential | 242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |

## Randomized pure-C properties

The property ledger is separate from the VEP differential. It records
successful runs of each randomized oracle, including the exact seed and
duplicate count. A failed suite does not append rows.

| run_date   | source_revision | seed               | randomized_targets | trials    | passed    | failed | duplicates | suite_tests | suite_assertions | suite_elapsed_seconds | compiler                                   |
|:-----------|:----------------|:-------------------|-------------------:|:----------|:----------|-------:|-----------:|------------:|:-----------------|----------------------:|:-------------------------------------------|
| 2026-07-11 | 8cc22218        | 0xd0c0ffee12345678 |                 39 | 3,800,500 | 3,800,500 |      0 |          0 |         133 | 189,981          |                  15.4 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

| target                                                                 | trials  | passed  | failed | skipped | duplicates |
|:-----------------------------------------------------------------------|:--------|:--------|:-------|:--------|:-----------|
| allele-normalized sweep candidates == independent trim oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor cross-codon MNV route == tile                          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor DEL route == tile under output splits                  | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor INS route == tile under output splits                  | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor MNV route == tile under output splits                  | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor output splits == one annotate_tile                     | 100,000 | 100,000 | 0      | 0       | 0          |
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
| sequence delta annotation wrapper MNV == legacy                        | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch INDEL == in-frame delins oracle                 | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch MNV == single-codon oracle                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch two-codon MNV window == codon-window oracle     | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence-backed SNV codon edit == codon-slice edit oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
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

| consequence_class                   | impact   |   n | exact_agree | unresolved | resolved_n | resolved_agree | resolved_discordant | term_mismatch | engine_extra | engine_missing | resolved_error_upper_95 |
|:------------------------------------|:---------|----:|------------:|-----------:|-----------:|---------------:|--------------------:|--------------:|-------------:|---------------:|:------------------------|
| frameshift_variant                  | HIGH     |  54 |          31 |         12 |         42 |             31 |                  11 |            23 |           10 |             13 | 42.04%                  |
| stop_lost                           | HIGH     |  21 |           9 |          4 |         17 |              9 |                   8 |            12 |            0 |             12 | 72.19%                  |
| 3_prime_UTR_variant                 | MODIFIER |  29 |          25 |          6 |         23 |             20 |                   3 |             4 |            2 |              2 | 33.59%                  |
| stop_retained_variant               | LOW      |   6 |           2 |          2 |          4 |              2 |                   2 |             4 |            0 |              4 | 93.24%                  |
| downstream_gene_variant             | MODIFIER |   5 |           3 |          0 |          5 |              3 |                   2 |             2 |            0 |              2 | 85.34%                  |
| inframe_insertion                   | MODERATE |   8 |           5 |          2 |          6 |              5 |                   1 |             3 |            1 |              2 | 64.12%                  |
| protein_altering_variant            | MODERATE |   2 |           1 |          0 |          2 |              1 |                   1 |             1 |            0 |              1 | 98.74%                  |
| coding_sequence_variant             | MODIFIER |  33 |           9 |         33 |          0 |              0 |                   0 |            24 |           24 |              0 | not estimable           |
| start_lost                          | HIGH     |  27 |          12 |         15 |         12 |             12 |                   0 |            15 |            0 |             15 | 26.46%                  |
| splice_donor_variant                | HIGH     |  24 |          24 |          4 |         20 |             20 |                   0 |             0 |            0 |              0 | 16.84%                  |
| splice_acceptor_variant             | HIGH     |  17 |          17 |          4 |         13 |             13 |                   0 |             0 |            0 |              0 | 24.71%                  |
| inframe_deletion                    | MODERATE |   8 |           6 |          2 |          6 |              6 |                   0 |             2 |            0 |              2 | 45.93%                  |
| splice_region_variant               | LOW      |  55 |          55 |          1 |         54 |             54 |                   0 |             0 |            0 |              0 | 6.60%                   |
| stop_gained                         | HIGH     |   2 |           1 |          1 |          1 |              1 |                   0 |             1 |            0 |              1 | 97.50%                  |
| start_retained_variant              | LOW      |   1 |           0 |          1 |          0 |              0 |                   0 |             1 |            0 |              1 | not estimable           |
| intron_variant                      | MODIFIER |  32 |          32 |          0 |         32 |             32 |                   0 |             0 |            0 |              0 | 10.89%                  |
| 5_prime_UTR_variant                 | MODIFIER |  22 |          22 |          0 |         22 |             22 |                   0 |             0 |            0 |              0 | 15.44%                  |
| synonymous_variant                  | LOW      |  13 |          13 |          0 |         13 |             13 |                   0 |             0 |            0 |              0 | 24.71%                  |
| splice_donor_region_variant         | LOW      |  12 |          12 |          0 |         12 |             12 |                   0 |             0 |            0 |              0 | 26.46%                  |
| splice_polypyrimidine_tract_variant | LOW      |  11 |          11 |          0 |         11 |             11 |                   0 |             0 |            0 |              0 | 28.49%                  |
| missense_variant                    | MODERATE |  10 |          10 |          0 |         10 |             10 |                   0 |             0 |            0 |              0 | 30.85%                  |
| splice_donor_5th_base_variant       | LOW      |   7 |           7 |          0 |          7 |              7 |                   0 |             0 |            0 |              0 | 40.96%                  |

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once.

| impact   |   n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:---------|----:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| HIGH     | 123 |          93 |         29 |         94 |                   9 | 75.61%     | 17.40%                  |
| LOW      |  45 |          42 |          1 |         44 |                   2 | 93.33%     | 15.47%                  |
| MODERATE |  24 |          22 |          1 |         23 |                   1 | 91.67%     | 21.95%                  |
| MODIFIER |  50 |          46 |          2 |         48 |                   3 | 92.00%     | 17.20%                  |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.
