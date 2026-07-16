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

## Latest tested revision per corpus

| revision | corpus                          | model                      | assembly       | species               | oracle_source | oracle    |  pairs | exact         | unresolved | resolved_disagreements | resolved_error_upper_95 |
|:---------|:--------------------------------|:---------------------------|:---------------|:----------------------|:--------------|:----------|-------:|:--------------|-----------:|-----------------------:|:------------------------|
| b204dd49 | GRCh38 dbSNP                    | final-dbsnp                | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  73620 | 73620/73620   |          0 |                      0 | 0.01%                   |
| b204dd49 | GRCh38 GIAB                     | final-giab                 | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  54905 | 54905/54905   |          0 |                      0 | 0.01%                   |
| 22803a40 | GRCh38 ClinVar coding           | final-coding               | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 287836 | 287836/287836 |          0 |                      0 | 0.00%                   |
| 22803a40 | GRCh38 ClinVar cross-chromosome | final-clinvar              | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 316397 | 316397/316397 |          0 |                      0 | 0.00%                   |
| 7dd90ce8 | GRCh37                          | final-grch37               | GRCh37         | homo_sapiens          | cache         | VEP 116.0 | 486464 | 486464/486464 |          0 |                      0 | 0.00%                   |
| 7dd90ce8 | P. falciparum                   | plasmodium-falciparum-63   | GCA000002765v3 | plasmodium_falciparum | cache         | VEP 116.0 |  40732 | 40732/40732   |          0 |                      0 | 0.01%                   |
| 360619ed | GRCh38 paired BND               | grch38_breakend_multichrom | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  91428 | 91428/91428   |          0 |                      0 | 0.00%                   |

Each row is the newest tested ancestor of the current source for that
named corpus. Expensive corpora do not inherit evidence from a later run
of another corpus, and a newly tested corpus does not hide older
still-applicable evidence. The SO and impact tables keep the same runs
separate so the largest corpus cannot hide a smaller species- or
assembly-specific frontier.

## Paired-breakend differential

| revision | generated_events | transcript_pairs | exact         | unresolved | extra | missing | resolved_error_upper_95 |
|:---------|:-----------------|:-----------------|:--------------|-----------:|------:|--------:|:------------------------|
| 360619ed | 1,004            | 91,428           | 91,428/91,428 |          0 |     0 |       0 | 0.00%                   |

The generated seed-31 corpus spans chromosomes 1, 2, 7, 21, and X;
intra- and interchromosomal mates; all four VCF bracket orientations;
and transcript, exon, intron, CDS, and directional-flank endpoint
states. The comparison is the union of consequences produced by both
breakend endpoints for each transcript, which is VEP 116’s
transcript-level paired-breakend contract.

VEP 116’s buffered BND path inserts mate coordinates into a
chromosome-blind interval tree. A multichromosome batch can therefore
omit valid transcript pairs even when the input chromosomes are
contiguous and position-sorted. The oracle command uses
`--buffer_size 1` in one Perl process so every event is evaluated
independently. This is oracle isolation, not a DuckVEP compatibility
rule; the ledger records `breakend_buffer_size=1` and the artifact hash.

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

| run_date   | source_revision | corpus                                | model                      |       n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:-----------|:----------------|:--------------------------------------|:---------------------------|--------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| 2026-07-11 | 8cc22218        | witnesses                             | differential               |     242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29             | differential               |  100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                   |
| 2026-07-13 | 24bb1714        | witnesses                             | differential               |     242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                   |
| 2026-07-13 | 34b37ca1        | witnesses                             | differential               |     242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                   |
| 2026-07-13 | 87f03a2a        | witnesses                             | differential               |     242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                  |
| 2026-07-13 | defc9a1c        | state_exploration_seed_113            | differential               |  100246 |       85598 |      28442 |      71804 |                1086 | 85.39%     | 1.60%                   |
| 2026-07-13 | defc9a1c        | state_exploration_seed_71             | differential               |  100242 |       85646 |      27946 |      72296 |                1103 | 85.44%     | 1.62%                   |
| 2026-07-13 | defc9a1c        | witnesses                             | differential               |     246 |         242 |         11 |        235 |                   0 | 98.37%     | 1.56%                   |
| 2026-07-13 | eb212de3        | witnesses                             | differential               |     242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_197            | differential               |  100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_211            | differential               |  100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_71             | differential               |  100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 2ab08e2f        | witnesses                             | differential               |     258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_113            | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_197            | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_211            | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | state_exploration_seed_71             | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 3c427df4        | witnesses                             | differential               |     268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                   |
| 2026-07-14 | 5e5bc1e2        | clinvar_chr21_seed1                   | ensembl116_grch38_core     |  126320 |      126320 |          0 |     126320 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-14 | 8498b92a        | final_clinvar_coding_seed113          | final-coding               |  287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_clinvar_crosschrom_seed17       | final-clinvar              |  316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_dbsnp157_windows_seed29         | final-dbsnp                |   73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_giab_grch38_seed71              | final-giab                 |   54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | 8498b92a        | final_grch37_cache_seed37             | final-grch37               |  486468 |      482665 |        102 |     486366 |                3747 | 99.22%     | 0.80%                   |
| 2026-07-14 | 8498b92a        | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63   |   40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_197            | differential               |  100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_211            | differential               |  100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                   |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_71             | differential               |  100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                   |
| 2026-07-14 | 8b2a2dbc        | witnesses                             | differential               |     258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                   |
| 2026-07-14 | b204dd49        | final_clinvar_coding_seed113          | final-coding               |  287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                   |
| 2026-07-14 | b204dd49        | final_clinvar_crosschrom_seed17       | final-clinvar              |  316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_dbsnp157_windows_seed29         | final-dbsnp                |   73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_giab_grch38_seed71              | final-giab                 |   54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-14 | b204dd49        | final_grch37_cache_seed37             | final-grch37               |  486468 |      486332 |        102 |     486366 |                  80 | 99.97%     | 0.02%                   |
| 2026-07-14 | b204dd49        | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63   |   40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_197            | differential               |  100248 |       91159 |      19473 |      80775 |                 537 | 90.93%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_211            | differential               |  100250 |       90951 |      19676 |      80574 |                 567 | 90.72%     | 0.76%                   |
| 2026-07-14 | fe6f0634        | state_exploration_seed_71             | differential               |  100242 |       90932 |      19500 |      80742 |                 532 | 90.71%     | 0.72%                   |
| 2026-07-14 | fe6f0634        | witnesses                             | differential               |     262 |         258 |         10 |        252 |                   0 | 98.47%     | 1.45%                   |
| 2026-07-15 | 7dd90ce8        | final_grch37_cache_seed37             | final-grch37               |  486464 |      486464 |          0 |     486464 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-15 | 7dd90ce8        | plasmodium-falciparum-vep63-seed11663 | plasmodium-falciparum-63   |   40732 |       40732 |          0 |      40732 |                   0 | 100.00%    | 0.01%                   |
| 2026-07-15 | c361346f        | nmd_clinvar_chr21                     | ensembl116-grch38-final    | 1331664 |     1331664 |          0 |    1331664 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 22803a40        | final_clinvar_coding_seed113          | final-coding               |  287836 |      287836 |          0 |     287836 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 22803a40        | final_clinvar_crosschrom_seed17       | final-clinvar              |  316397 |      316397 |          0 |     316397 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr1_seed29              | differential               |  124896 |      124896 |          0 |     124896 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr11_seed307            | differential               |  528847 |      528847 |          0 |     528847 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr17_seed97             | differential               |  120821 |      120821 |          0 |     120821 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr2_seed211             | differential               |  484044 |      484044 |          0 |     484044 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr22_seed401            | differential               |  547182 |      547182 |          0 |     547182 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr6_seed71              | differential               |  110704 |      110704 |          0 |     110704 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_chrX_seed113             | differential               |   98072 |       98072 |          0 |      98072 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | generated_sv_seed17                   | differential               |  126345 |      126345 |          0 |     126345 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 24a5cf2a        | state_exploration_seed_20260716       | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                   |
| 2026-07-16 | 360619ed        | breakend_multichrom_seed31_isolated   | grch38_breakend_multichrom |   91428 |       91428 |          0 |      91428 |                   0 | 100.00%    | 0.00%                   |

## Randomized executable-VEP state exploration

This is the anti-overfitting lane against the VEP executable, not an
internal property test. Each seed contains the fixed predicate witnesses
plus 100,000 unique alleles. Three quarters of random positions are
within six bases of splice, exon, translation-start, and
translation-stop boundaries; one quarter is uniform across the
transcript. SNVs, MNVs, insertions, deletions, and delins are sampled
with equal probability, with differing alleles up to 49 bases.

| revision | seed     | pairs   | exact   | unresolved | resolved_disagreements | exact_error_upper_95_ppm |
|:---------|:---------|:--------|:--------|:-----------|:-----------------------|:-------------------------|
| 24a5cf2a | 20260716 | 100,268 | 100,268 | 0          | 0                      | 36.8                     |
| 24a5cf2a | combined | 100,268 | 100,268 | 0          | 0                      | 36.8                     |

The same campaign covered the following SO terms. Counts are term
memberships, not distinct transcript pairs, because one pair may carry
several terms and the four seeds deliberately retain the same fixed
witnesses.

| consequence_class                   | impact   | seeds_observed |     n | unresolved | term_mismatch | engine_extra | engine_missing |
|:------------------------------------|:---------|---------------:|------:|-----------:|--------------:|-------------:|---------------:|
| intron_variant                      | MODIFIER |              1 | 27021 |          0 |             0 |            0 |              0 |
| frameshift_variant                  | HIGH     |              1 | 19827 |          0 |             0 |            0 |              0 |
| 5_prime_UTR_variant                 | MODIFIER |              1 | 13501 |          0 |             0 |            0 |              0 |
| 3_prime_UTR_variant                 | MODIFIER |              1 | 12213 |          0 |             0 |            0 |              0 |
| missense_variant                    | MODERATE |              1 | 11147 |          0 |             0 |            0 |              0 |
| splice_region_variant               | LOW      |              1 | 10302 |          0 |             0 |            0 |              0 |
| splice_polypyrimidine_tract_variant | LOW      |              1 |  9985 |          0 |             0 |            0 |              0 |
| coding_sequence_variant             | MODIFIER |              1 |  7825 |          0 |             0 |            0 |              0 |
| start_lost                          | HIGH     |              1 |  4832 |          0 |             0 |            0 |              0 |
| stop_gained                         | HIGH     |              1 |  4592 |          0 |             0 |            0 |              0 |
| splice_acceptor_variant             | HIGH     |              1 |  4534 |          0 |             0 |            0 |              0 |
| splice_donor_variant                | HIGH     |              1 |  4310 |          0 |             0 |            0 |              0 |
| inframe_insertion                   | MODERATE |              1 |  4156 |          0 |             0 |            0 |              0 |
| protein_altering_variant            | MODERATE |              1 |  3384 |          0 |             0 |            0 |              0 |
| splice_donor_region_variant         | LOW      |              1 |  3189 |          0 |             0 |            0 |              0 |
| splice_donor_5th_base_variant       | LOW      |              1 |  3166 |          0 |             0 |            0 |              0 |
| stop_lost                           | HIGH     |              1 |  2945 |          0 |             0 |            0 |              0 |
| downstream_gene_variant             | MODIFIER |              1 |   709 |          0 |             0 |            0 |              0 |
| inframe_deletion                    | MODERATE |              1 |   424 |          0 |             0 |            0 |              0 |
| start_retained_variant              | LOW      |              1 |   410 |          0 |             0 |            0 |              0 |
| stop_retained_variant               | LOW      |              1 |   381 |          0 |             0 |            0 |              0 |
| synonymous_variant                  | LOW      |              1 |   269 |          0 |             0 |            0 |              0 |
| intergenic_variant                  | MODIFIER |              1 |    25 |          0 |             0 |            0 |              0 |

This distribution deliberately stresses local allele and boundary states
on one engineered transcript. It does not replace the indexed-cache
corpora, which add real transcript density, imported flags, exceptional
peptide edits, codon tables, assemblies, and species. The revision is
shown explicitly because this expensive campaign is not silently
attributed to later code.

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
| 2026-07-16 | 22803a40        | 0x0000000020260716 |                 43 | 4,200,500  | 4,200,500  |      0 |          0 |         170 | 205,585          |                25.162 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-16 | 360619ed        | 0x0000000020260716 |                 43 | 4,200,500  | 4,200,500  |      0 |          0 |         171 | 205,610          |                24.408 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

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
| phased SNV set == equivalent MNV coding facts                          | 100,000 | 100,000 | 0      | 0       | 0          |
| region mask structural invariants                                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta annotation wrapper MNV == direct shape                  | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta exon hint == unhinted projection                        | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch INDEL == local delins-shape oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch MNV == single-codon oracle                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch two-codon MNV window == codon-window oracle     | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence-backed SNV codon edit == codon-slice edit oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
| sorted point cursor classifier == exhaustive exon/gap scans            | 100,000 | 100,000 | 0      | 0       | 0          |
| sorted span cursor classifier == exhaustive exon/gap scans             | 100,000 | 100,000 | 0      | 0       | 0          |
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

| corpus                          | observed_terms | terms_with_mismatch | term_mismatches | engine_extra | engine_missing | terms_with_unresolved | unresolved_term_memberships |
|:--------------------------------|---------------:|--------------------:|----------------:|-------------:|---------------:|----------------------:|----------------------------:|
| GRCh38 dbSNP                    |             22 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 GIAB                     |             19 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 ClinVar coding           |             27 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 ClinVar cross-chromosome |             28 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh37                          |             26 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| P. falciparum                   |             20 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 paired BND               |             14 |                   0 |               0 |            0 |              0 |                     0 |                           0 |

`term_mismatches` counts an SO term that is missing or extra on a
transcript pair. `unresolved_term_memberships` is reported separately:
an unresolved pair can still carry the exact VEP term set, and a
multi-term pair appears once under each term.

No SO-term mismatch or unresolved membership remains in the latest
declared runs.

Terms absent from this frontier table were exact and resolved everywhere
they were observed. The complete zero and nonzero strata remain in the
CSV ledger.

## VEP impact classes

This table uses full consequence sets, so each transcript pair is
counted once within each corpus.

| corpus                          | impact   |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | resolved_error_upper_95 |
|:--------------------------------|:---------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:------------------------|
| GRCh38 paired BND               | HIGH     |  69654 |       69654 |          0 |      69654 |                   0 | 100.00%    | 0.01%                   |
| GRCh38 ClinVar coding           | HIGH     |  99103 |       99103 |          0 |      99103 |                   0 | 100.00%    | 0.00%                   |
| GRCh38 ClinVar cross-chromosome | HIGH     |  79663 |       79663 |          0 |      79663 |                   0 | 100.00%    | 0.00%                   |
| GRCh38 dbSNP                    | HIGH     |    118 |         118 |          0 |        118 |                   0 | 100.00%    | 3.08%                   |
| GRCh38 GIAB                     | HIGH     |      4 |           4 |          0 |          4 |                   0 | 100.00%    | 60.24%                  |
| GRCh37                          | HIGH     |  46645 |       46645 |          0 |      46645 |                   0 | 100.00%    | 0.01%                   |
| P. falciparum                   | HIGH     |   4309 |        4309 |          0 |       4309 |                   0 | 100.00%    | 0.09%                   |
| GRCh38 ClinVar coding           | LOW      |   5952 |        5952 |          0 |       5952 |                   0 | 100.00%    | 0.06%                   |
| GRCh38 ClinVar cross-chromosome | LOW      |  28776 |       28776 |          0 |      28776 |                   0 | 100.00%    | 0.01%                   |
| GRCh38 dbSNP                    | LOW      |    349 |         349 |          0 |        349 |                   0 | 100.00%    | 1.05%                   |
| GRCh38 GIAB                     | LOW      |    260 |         260 |          0 |        260 |                   0 | 100.00%    | 1.41%                   |
| GRCh37                          | LOW      |  23354 |       23354 |          0 |      23354 |                   0 | 100.00%    | 0.02%                   |
| P. falciparum                   | LOW      |    210 |         210 |          0 |        210 |                   0 | 100.00%    | 1.74%                   |
| GRCh38 ClinVar coding           | MODERATE |  59108 |       59108 |          0 |      59108 |                   0 | 100.00%    | 0.01%                   |
| GRCh38 ClinVar cross-chromosome | MODERATE |  43286 |       43286 |          0 |      43286 |                   0 | 100.00%    | 0.01%                   |
| GRCh38 dbSNP                    | MODERATE |    125 |         125 |          0 |        125 |                   0 | 100.00%    | 2.91%                   |
| GRCh38 GIAB                     | MODERATE |     31 |          31 |          0 |         31 |                   0 | 100.00%    | 11.22%                  |
| GRCh37                          | MODERATE |  23319 |       23319 |          0 |      23319 |                   0 | 100.00%    | 0.02%                   |
| P. falciparum                   | MODERATE |   1937 |        1937 |          0 |       1937 |                   0 | 100.00%    | 0.19%                   |
| GRCh38 paired BND               | MODIFIER |  21774 |       21774 |          0 |      21774 |                   0 | 100.00%    | 0.02%                   |
| GRCh38 ClinVar coding           | MODIFIER | 123673 |      123673 |          0 |     123673 |                   0 | 100.00%    | 0.00%                   |
| GRCh38 ClinVar cross-chromosome | MODIFIER | 164672 |      164672 |          0 |     164672 |                   0 | 100.00%    | 0.00%                   |
| GRCh38 dbSNP                    | MODIFIER |  73028 |       73028 |          0 |      73028 |                   0 | 100.00%    | 0.01%                   |
| GRCh38 GIAB                     | MODIFIER |  54610 |       54610 |          0 |      54610 |                   0 | 100.00%    | 0.01%                   |
| GRCh37                          | MODIFIER | 393146 |      393146 |          0 |     393146 |                   0 | 100.00%    | 0.00%                   |
| P. falciparum                   | MODIFIER |  34276 |       34276 |          0 |      34276 |                   0 | 100.00%    | 0.01%                   |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.

## Variant-induced NMD

This is a separate executable differential against the pinned VEP
Plugins release/116 `NMD.pm`. It compares `triggering`, `escaping`, and
`unresolved` for every eligible transcript pair; it does not infer NMD
from the core `NMD_transcript_variant` biotype consequence.

| revision | corpus            | model                   | exact       | mismatches | not comparable | VEP unresolved | DuckVEP unresolved | resolved error upper 95% |
|:---------|:------------------|:------------------------|:------------|-----------:|---------------:|---------------:|-------------------:|:-------------------------|
| c361346f | nmd_clinvar_chr21 | ensembl116-grch38-final | 68554/68554 |          0 |              0 |          29416 |              29416 | 0.01%                    |

The ledger keeps the prediction confusion matrix rather than only the
total:

| revision | corpus            | VEP_prediction | DuckVEP_prediction |     n |
|:---------|:------------------|:---------------|:-------------------|------:|
| c361346f | nmd_clinvar_chr21 | escaping       | escaping           |  6937 |
| c361346f | nmd_clinvar_chr21 | triggering     | triggering         | 32201 |
| c361346f | nmd_clinvar_chr21 | unresolved     | unresolved         | 29416 |

VEP projects the complete uploaded `VariationFeature` for the plugin’s
CDS and exon-position rules. DuckVEP therefore retains both geometries:
minimized edit coordinates drive consequence and sequence changes, while
the original feature endpoints drive NMD. A padded and a minimal allele
can encode the same sequence edit but cross the plugin’s inclusive
positional threshold differently.
