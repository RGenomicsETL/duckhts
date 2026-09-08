DuckVEP conformance against Ensembl VEP 116
================

<!-- duckvep_conformance.md is generated from duckvep_conformance.Rmd. -->

This report records exact consequence agreement with the real Ensembl
VEP 116 executable at the object VEP annotates: a transcript,
RegulatoryFeature, or MotifFeature. It uses either VEP’s declared
indexed cache or a staged GFF oracle. Unresolved DuckVEP rows remain in
the denominator; they are not discarded as unsupported cases. The CSV is
append-only by source revision, corpus, and resident model. Independent
frozen distributions and seeds are kept separate so a fix cannot improve
its own hand-picked witnesses and hide a regression elsewhere.

Official Ensembl Variation release VCFs provide a separate product-audit
lane. Their indexed `VE` relation can be compared in ordinary CI without
starting Perl VEP; `CSQ` is a lossy presentation of those stored rows.
This is not an executable-VEP oracle: in release 116, `X/Y:276322 G>A`
is published as `intergenic_variant`, while cache-mode VEP with
`--distance 0` emits three path-specific `5_prime_UTR_variant` rows on
each chromosome. Pinned release shards are therefore useful lineage
evidence, while the executable/cache combination remains the semantic
compatibility authority. Matching full models belong in external
versioned artifacts, not in git or the network-free extension build
step.

## Declared conformance closure

The independent-event consequence engine is closed as a semantic
implementation campaign for the declared model and event surfaces:
admitted Ensembl transcript, mature-miRNA, RegulatoryFeature, and
MotifFeature objects; independent literal small alleles; exact typed
DEL, DUP, tandem-DUP, INV, INS, and CNV events; structural tandem
repeats (`STR`); paired breakends; supported BioPerl codon tables and
exceptional Ensembl peptide edits; and the separately declared VEP
NMD-plugin result. DEL/DUP/tandem-DUP/INV/INS/CNV and BND have generated
executable-VEP differentials. Structural `STR` has source-derived
VEP-116 semantics plus fixed SQL/R and randomized C coverage; raw repeat
reconstruction is a separate input-preparation operation. The evidence
spans GRCh38, GRCh37, and *P. falciparum*, executable witnesses,
indexed-cache corpora, generated state exploration, sanitizer runs, and
pure-C oracle properties.

“Closed” means future consequence changes are routine engineering behind
these regression gates. VEP 116 parses `CIPOS`/`CIEND` into inner/outer
structural coordinates, but its registered consequence predicates use
nominal `POS`/`END`; DuckVEP therefore annotates that nominal span while
the surrounding relation preserves the uncertainty metadata. The
checked-in 12-record GRCh38 confidence witness records this directly:
nominal and `IMPRECISE;CIPOS;CIEND` forms of CNV, DEL, DUP, tandem DUP,
INV, and INS produced 466/466 exact transcript pairs, and both engines
had equal nominal/imprecise consequence multisets for all six event-kind
pairs. VEP can also expand a bounded `<CNV:TR>` from `RN`, `RUS`, and
`RUC` or `RB` into a literal allele before consequence calculation.
Implementing that lossless expansion and mapping VEP’s finite supported
symbolic vocabulary into the typed event API are narrower input-
preparation tasks, not missing consequence predicates. VEP itself
rejects unrecognised types such as CPX, so this closure does not promise
arbitrary symbolic parsing. Untested species/releases and phased
multi-record haplotypes remain outside the closure. Haplotype grouping
and combined consequence attribution are the next semantic vertical and
require their own executable oracle and performance campaign. Any newly
observed fixed-event mismatch reopens this contract rather than being
relabelled as unsupported.

## Evidence units and statistical interpretation

The release evidence has three distinct units. A complete-corpus audit
establishes exactness only for every event in that named corpus. The
generated state-exploration campaign deliberately over-samples rare
splice, coding, strand, length, and structural states to discover
defects and demonstrate transition coverage; it is not a sample of
deployment prevalence. A population error-rate statement would require a
separately specified probability sample whose primary outcome is whether
the complete output multiset for each independently sampled input event
differs from VEP.

The append-only history retains its original `upper95` field for schema
compatibility. Tables below label that value as a *descriptive
independent-pair* Clopper–Pearson calculation. Transcript and other
annotated-object pairs are clustered within input events, transcripts,
and genes, and the targeted generators are intentionally not deployment
samples. These calculations therefore are not general engine error-rate
bounds; separate stratum values are also not simultaneous confidence
intervals. The release gate itself is deterministic and stricter: any
discordance, unresolved state, extra emission, or missing emission fails
the audited run.

## Latest tested revision per corpus

| revision | corpus                            | model                      | assembly       | species               | oracle_source | oracle    |  pairs | exact         | unresolved | resolved_disagreements | descriptive_independent_pair_upper_95 |
|:---------|:----------------------------------|:---------------------------|:---------------|:----------------------|:--------------|:----------|-------:|:--------------|-----------:|-----------------------:|:--------------------------------------|
| b204dd49 | GRCh38 dbSNP                      | final-dbsnp                | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  73620 | 73620/73620   |          0 |                      0 | 0.01%                                 |
| b204dd49 | GRCh38 GIAB                       | final-giab                 | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  54905 | 54905/54905   |          0 |                      0 | 0.01%                                 |
| 22803a40 | GRCh38 ClinVar coding             | final-coding               | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 287836 | 287836/287836 |          0 |                      0 | 0.00%                                 |
| 22803a40 | GRCh38 ClinVar cross-chromosome   | final-clinvar              | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 316397 | 316397/316397 |          0 |                      0 | 0.00%                                 |
| 7dd90ce8 | GRCh37                            | final-grch37               | GRCh37         | homo_sapiens          | cache         | VEP 116.0 | 486464 | 486464/486464 |          0 |                      0 | 0.00%                                 |
| 7dd90ce8 | P. falciparum                     | plasmodium-falciparum-63   | GCA000002765v3 | plasmodium_falciparum | cache         | VEP 116.0 |  40732 | 40732/40732   |          0 |                      0 | 0.01%                                 |
| 360619ed | GRCh38 paired BND                 | grch38_breakend_multichrom | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  91428 | 91428/91428   |          0 |                      0 | 0.00%                                 |
| 96b4cd45 | GRCh38 GIAB + core regulation     | differential               | GRCh38         | homo_sapiens          | cache         | VEP 116.0 |  14955 | 14955/14955   |          0 |                      0 | 0.02%                                 |
| 96b4cd45 | GRCh38 exact SV + core regulation | differential               | GRCh38         | homo_sapiens          | cache         | VEP 116.0 | 120224 | 120224/120224 |          0 |                      0 | 0.00%                                 |

Each row is the newest tested ancestor of the current source for that
named corpus. Expensive corpora do not inherit evidence from a later run
of another corpus, and a newly tested corpus does not hide older
still-applicable evidence. The SO and impact tables keep the same runs
separate so the largest corpus cannot hide a smaller species- or
assembly-specific frontier.

## Independent-event HGVS differential

| revision | corpus                          | model        | metric | exact           | match | both_absent | discordant |
|:---------|:--------------------------------|:-------------|:-------|:----------------|------:|------------:|-----------:|
| 6ce2ddd8 | clinvar_chr21_hgvs_seed113      | differential | HGVSC  | 56,998/56,998   | 44871 |       12127 |          0 |
| 6ce2ddd8 | clinvar_chr21_hgvs_seed113      | differential | HGVSP  | 56,998/56,998   | 20782 |       36216 |          0 |
| 6ce2ddd8 | hgvs_terminal_multiplication    | differential | HGVSC  | 4/4             |     1 |           3 |          0 |
| 6ce2ddd8 | hgvs_terminal_multiplication    | differential | HGVSP  | 4/4             |     0 |           4 |          0 |
| 6ce2ddd8 | state_exploration_seed_16180339 | differential | HGVSC  | 100,268/100,268 | 99169 |        1099 |          0 |
| 6ce2ddd8 | state_exploration_seed_16180339 | differential | HGVSP  | 100,268/100,268 | 31403 |       68865 |          0 |
| 6ce2ddd8 | state_exploration_seed_27182818 | differential | HGVSC  | 100,268/100,268 | 99135 |        1133 |          0 |
| 6ce2ddd8 | state_exploration_seed_27182818 | differential | HGVSP  | 100,268/100,268 | 31227 |       69041 |          0 |
| 6ce2ddd8 | state_exploration_seed_31415927 | differential | HGVSC  | 100,268/100,268 | 99146 |        1122 |          0 |
| 6ce2ddd8 | state_exploration_seed_31415927 | differential | HGVSP  | 100,268/100,268 | 31021 |       69247 |          0 |

| revision | corpus                          | extension_build               | extension    | model_kind | model        | reference    | reference_index | source_vcf   | input_vcf    | pair_artifact |
|:---------|:--------------------------------|:------------------------------|:-------------|:-----------|:-------------|:-------------|:----------------|:-------------|:-------------|:--------------|
| 6ce2ddd8 | clinvar_chr21_hgvs_seed113      | htslib_distclean_make_release | 3213f0a209bf | duckdb     | 8a59b14eed5c | 1e74081a49ce | 0998f61682f4    | 7ecec9a75071 | 7ecec9a75071 | f4df0ad05234  |
| 6ce2ddd8 | hgvs_terminal_multiplication    | htslib_distclean_make_release | 3213f0a209bf | sql        | b21fbeac2c28 | 01d1f0252130 | 154cbe440869    | bfa15d2786f3 | c4182bf1b769 | eb37a23b382c  |
| 6ce2ddd8 | state_exploration_seed_16180339 | htslib_distclean_make_release | 3213f0a209bf | sql        | b21fbeac2c28 | 01d1f0252130 | 154cbe440869    | 2d8315a4926a | 53150698e457 | ce7688057f5b  |
| 6ce2ddd8 | state_exploration_seed_27182818 | htslib_distclean_make_release | 3213f0a209bf | sql        | b21fbeac2c28 | 01d1f0252130 | 154cbe440869    | ec5a793adc6b | 8486a6b4c05e | 2aba7e180609  |
| 6ce2ddd8 | state_exploration_seed_31415927 | htslib_distclean_make_release | 3213f0a209bf | sql        | b21fbeac2c28 | 01d1f0252130 | 154cbe440869    | 1c5cbf73b5f6 | beab52a9d117 | 107ea9953774  |

This is exact string agreement for independent transcript events with
VEP 116 invoked using `--hgvs`. A comparison is exact when both engines
emit the same string or both omit that HGVS field. Unresolved, missing,
extra, and unequal strings remain discordant; none is removed from the
denominator. The checked ledger accepts only a pair artifact produced
from the current clean source revision by a vendored-htslib distclean
followed by an in-tree release build. The table retains complete SHA-256
receipts for the extension, model, FASTA and index, source VCF, exact
sampled VCF passed to VEP, and pair-level Parquet; shortened digests are
rendered above. Historical HGVS rows recorded before build receipts were
introduced remain in the append-only CSV but are not presented as
checked evidence.

The `b7c7237ee686` ClinVar HGVS and NMD runs use a freshly acquired
`vep116_grch38_cache_chr21`. Its complete 27,644,657,162-byte source
archive was verified against registry SHA-256
`014c7dd9bb5ad06665866d62eb80f31ca761197bcb9b59280300676f996e600d`
before publication. The cache contains 98 files / 326,767,268 bytes;
every retained file matches the preserved earlier cache byte-for-byte.
Earlier HTTP-identity-only acquisition receipts remain historical
observations, not content-checksum evidence. The refreshed HGVS run
retains the same 1,864-variant input and all 56,998 transcript pairs,
with no unresolved, missing, extra or discordant pairs.

The `6ce2ddd85df7` rerun again retains that exact source/input VCF and
all 56,998 pairs, with exact consequences and HGVSc/HGVSp agreement. It
resolves the model, reference and checksum-verified chromosome-21 cache
through the artifact registry. Its physical model-file receipt differs
from the older run; this is not a claim of byte-identical model files.
The full GIAB conformance campaign has not been refreshed to this
revision. A separate chromosome-21 attempt found no eligible model joins
because the raw callset uses `chr21` and the model uses `21`; it stopped
before VEP, was not counted as a pass, and did not replace the full GIAB
gate.

## Paired-breakend differential

| revision | generated_events | transcript_pairs | exact         | unresolved | extra | missing | descriptive_independent_pair_upper_95 |
|:---------|:-----------------|:-----------------|:--------------|-----------:|------:|--------:|:--------------------------------------|
| 360619ed | 1,004            | 91,428           | 91,428/91,428 |          0 |     0 |       0 | 0.00%                                 |

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

## Core regulation and motif differential

| workload            | consequence_class               | memberships | exact       |
|:--------------------|:--------------------------------|:------------|:------------|
| generated exact SVs | TFBS_ablation                   | 766         | 766/ 766    |
| generated exact SVs | TFBS_amplification              | 1,532       | 1,532/1,532 |
| GIAB chromosome 21  | TF_binding_site_variant         | 2           | 2/ 2        |
| generated exact SVs | TF_binding_site_variant         | 4,354       | 4,354/4,354 |
| generated exact SVs | regulatory_region_ablation      | 358         | 358/ 358    |
| generated exact SVs | regulatory_region_amplification | 716         | 716/ 716    |
| GIAB chromosome 21  | regulatory_region_variant       | 54          | 54/ 54      |
| generated exact SVs | regulatory_region_variant       | 2,689       | 2,689/2,689 |

The GIAB run checks ordinary alleles against transcript and core funcgen
objects. The generated structural run deliberately crosses, contains,
exactly matches, and partially overlaps RegulatoryFeature and
MotifFeature intervals under DEL, DUP, TDUP, INV, INS, and CNV
operations. Structural `STR` is covered separately by the source-derived
VEP-116 rule, fixed SQL/R adapter tests, and randomized C oracles; this
generated executable-VEP run does not reconstruct raw repeat metadata.
The resident model contains only VEP-admitted core funcgen objects: VEP
116 removes `epigenetically_modified_region` rows before constructing
RegulatoryFeature overlap objects, so DuckVEP excludes them during
deterministic model preparation rather than filtering output after
candidate traversal.

## Prepared model receipts

| revision | species               | release | assembly       | regions | transcripts | coding_backed | exons     | mature_miRNA_segments | peptide_edits | regulatory_regions | motif_features | codon_tables     | model_sha256 |
|:---------|:----------------------|:--------|:---------------|:--------|:------------|:--------------|:----------|:----------------------|:--------------|:-------------------|:---------------|:-----------------|:-------------|
| 96b4cd45 | homo_sapiens          | 116     | GRCh38         | 194     | 644,427     | 369,631       | 5,068,416 | 2,806                 | 389           | 380,818            | 1,002,762      | 1:369618;2:13    | 296bc9063356 |
| 8498b92a | homo_sapiens          | 116     | GRCh37         | 84      | 195,379     | 94,610        | 1,186,433 | 3,788                 | 129           | 0                  | 0              | 1:94597;2:13     | 25459e62e50d |
| 8498b92a | plasmodium_falciparum | 63/116  | GCA000002765v3 | 16      | 5,791       | 5,389         | 15,097    | 0                     | 4             | 0                  | 0              | 1:5356;4:3;11:30 | c011cdd4deab |

These are complete model-build receipts, not counts inferred from a
differential. The ledger retains the full source-manifest, reference,
and model SHA-256 values, the exact VEP transcript filter, every count
above, CDS/flank base totals, and the external artifact name. The
Plasmodium row is an Ensembl Genomes release-63 cache paired with the
VEP/core-116 executable libraries, which is why both release numbers are
recorded.

## History

| run_date   | source_revision | corpus                                                   | model                      |       n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | descriptive_independent_pair_upper_95 |
|:-----------|:----------------|:---------------------------------------------------------|:---------------------------|--------:|------------:|-----------:|-----------:|--------------------:|:-----------|:--------------------------------------|
| 2026-07-11 | 8cc22218        | witnesses                                                | differential               |     242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                                |
| 2026-07-13 | 24bb1714        | state_exploration_seed_29                                | differential               |  100242 |       85238 |      28109 |      72133 |                1453 | 85.03%     | 2.12%                                 |
| 2026-07-13 | 24bb1714        | witnesses                                                | differential               |     242 |         238 |         11 |        231 |                   0 | 98.35%     | 1.58%                                 |
| 2026-07-13 | 34b37ca1        | witnesses                                                | differential               |     242 |         209 |         32 |        210 |                  10 | 86.36%     | 8.58%                                 |
| 2026-07-13 | 87f03a2a        | witnesses                                                | differential               |     242 |         203 |         33 |        209 |                  15 | 83.88%     | 11.56%                                |
| 2026-07-13 | defc9a1c        | state_exploration_seed_113                               | differential               |  100246 |       85598 |      28442 |      71804 |                1086 | 85.39%     | 1.60%                                 |
| 2026-07-13 | defc9a1c        | state_exploration_seed_71                                | differential               |  100242 |       85646 |      27946 |      72296 |                1103 | 85.44%     | 1.62%                                 |
| 2026-07-13 | defc9a1c        | witnesses                                                | differential               |     246 |         242 |         11 |        235 |                   0 | 98.37%     | 1.56%                                 |
| 2026-07-13 | eb212de3        | witnesses                                                | differential               |     242 |         219 |         32 |        210 |                   0 | 90.50%     | 1.74%                                 |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_197                               | differential               |  100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                                 |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_211                               | differential               |  100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                                 |
| 2026-07-14 | 2ab08e2f        | state_exploration_seed_71                                | differential               |  100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                                 |
| 2026-07-14 | 2ab08e2f        | witnesses                                                | differential               |     258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                                 |
| 2026-07-14 | 3c427df4        | state_exploration_seed_113                               | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-14 | 3c427df4        | state_exploration_seed_197                               | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-14 | 3c427df4        | state_exploration_seed_211                               | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-14 | 3c427df4        | state_exploration_seed_71                                | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-14 | 3c427df4        | witnesses                                                | differential               |     268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                                 |
| 2026-07-14 | 5e5bc1e2        | clinvar_chr21_seed1                                      | ensembl116_grch38_core     |  126320 |      126320 |          0 |     126320 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-14 | 8498b92a        | final_clinvar_coding_seed113                             | final-coding               |  287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                                 |
| 2026-07-14 | 8498b92a        | final_clinvar_crosschrom_seed17                          | final-clinvar              |  316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                                 |
| 2026-07-14 | 8498b92a        | final_dbsnp157_windows_seed29                            | final-dbsnp                |   73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-14 | 8498b92a        | final_giab_grch38_seed71                                 | final-giab                 |   54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-14 | 8498b92a        | final_grch37_cache_seed37                                | final-grch37               |  486468 |      482665 |        102 |     486366 |                3747 | 99.22%     | 0.80%                                 |
| 2026-07-14 | 8498b92a        | plasmodium-falciparum-vep63-seed11663                    | plasmodium-falciparum-63   |   40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                                 |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_197                               | differential               |  100248 |       88021 |      22598 |      77650 |                 550 | 87.80%     | 0.77%                                 |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_211                               | differential               |  100250 |       87815 |      22801 |      77449 |                 578 | 87.60%     | 0.81%                                 |
| 2026-07-14 | 8b2a2dbc        | state_exploration_seed_71                                | differential               |  100242 |       87916 |      22502 |      77740 |                 546 | 87.70%     | 0.76%                                 |
| 2026-07-14 | 8b2a2dbc        | witnesses                                                | differential               |     258 |         254 |         10 |        248 |                   0 | 98.45%     | 1.48%                                 |
| 2026-07-14 | b204dd49        | final_clinvar_coding_seed113                             | final-coding               |  287859 |      287829 |          4 |     287855 |                  27 | 99.99%     | 0.01%                                 |
| 2026-07-14 | b204dd49        | final_clinvar_crosschrom_seed17                          | final-clinvar              |  316399 |      316388 |          2 |     316397 |                  10 | 100.00%    | 0.01%                                 |
| 2026-07-14 | b204dd49        | final_dbsnp157_windows_seed29                            | final-dbsnp                |   73620 |       73620 |          0 |      73620 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-14 | b204dd49        | final_giab_grch38_seed71                                 | final-giab                 |   54905 |       54905 |          0 |      54905 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-14 | b204dd49        | final_grch37_cache_seed37                                | final-grch37               |  486468 |      486332 |        102 |     486366 |                  80 | 99.97%     | 0.02%                                 |
| 2026-07-14 | b204dd49        | plasmodium-falciparum-vep63-seed11663                    | plasmodium-falciparum-63   |   40734 |       40730 |         24 |      40710 |                   4 | 99.99%     | 0.03%                                 |
| 2026-07-14 | fe6f0634        | state_exploration_seed_197                               | differential               |  100248 |       91159 |      19473 |      80775 |                 537 | 90.93%     | 0.72%                                 |
| 2026-07-14 | fe6f0634        | state_exploration_seed_211                               | differential               |  100250 |       90951 |      19676 |      80574 |                 567 | 90.72%     | 0.76%                                 |
| 2026-07-14 | fe6f0634        | state_exploration_seed_71                                | differential               |  100242 |       90932 |      19500 |      80742 |                 532 | 90.71%     | 0.72%                                 |
| 2026-07-14 | fe6f0634        | witnesses                                                | differential               |     262 |         258 |         10 |        252 |                   0 | 98.47%     | 1.45%                                 |
| 2026-07-15 | 7dd90ce8        | final_grch37_cache_seed37                                | final-grch37               |  486464 |      486464 |          0 |     486464 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-15 | 7dd90ce8        | plasmodium-falciparum-vep63-seed11663                    | plasmodium-falciparum-63   |   40732 |       40732 |          0 |      40732 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-15 | c361346f        | nmd_clinvar_chr21                                        | ensembl116-grch38-final    | 1331664 |     1331664 |          0 |    1331664 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 22803a40        | final_clinvar_coding_seed113                             | final-coding               |  287836 |      287836 |          0 |     287836 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 22803a40        | final_clinvar_crosschrom_seed17                          | final-clinvar              |  316397 |      316397 |          0 |     316397 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr1_seed29                                 | differential               |  124896 |      124896 |          0 |     124896 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr11_seed307                               | differential               |  528847 |      528847 |          0 |     528847 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr17_seed97                                | differential               |  120821 |      120821 |          0 |     120821 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr2_seed211                                | differential               |  484044 |      484044 |          0 |     484044 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr22_seed401                               | differential               |  547182 |      547182 |          0 |     547182 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chr6_seed71                                 | differential               |  110704 |      110704 |          0 |     110704 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_chrX_seed113                                | differential               |   98072 |       98072 |          0 |      98072 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | generated_sv_seed17                                      | differential               |  126345 |      126345 |          0 |     126345 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 24a5cf2a        | state_exploration_seed_20260716                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 360619ed        | breakend_multichrom_seed31_isolated                      | grch38_breakend_multichrom |   91428 |       91428 |          0 |      91428 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-16 | 96b4cd45        | regulation_giab_chr21_seed1                              | differential               |   14955 |       14955 |          0 |      14955 |                   0 | 100.00%    | 0.02%                                 |
| 2026-07-16 | 96b4cd45        | regulation_sv_chr21_seed17                               | differential               |  120224 |      120224 |          0 |     120224 |                   0 | 100.00%    | 0.00%                                 |
| 2026-07-19 | e7c3623d        | breakend_regulation_chr21_22_seed20260719_distance_0     | breakend_distance_0        |   22380 |       22380 |          0 |      22380 |                   0 | 100.00%    | 0.02%                                 |
| 2026-07-19 | e7c3623d        | breakend_regulation_chr21_22_seed20260719_distance_10000 | breakend_distance_10000    |   29304 |       29304 |          0 |      29304 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-19 | e7c3623d        | breakend_regulation_chr21_22_seed20260719_distance_137   | breakend_distance_137      |   24970 |       24970 |          0 |      24970 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-19 | e7c3623d        | breakend_regulation_chr21_22_seed20260719_distance_5000  | breakend_distance_5000     |   29304 |       29304 |          0 |      29304 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-19 | f97101e1        | breakend_regulation_chr21_22_seed20260719_distance_0     | breakend_distance_0        |   22380 |       22380 |          0 |      22380 |                   0 | 100.00%    | 0.02%                                 |
| 2026-07-19 | f97101e1        | breakend_regulation_chr21_22_seed20260719_distance_10000 | breakend_distance_10000    |   29304 |       29304 |          0 |      29304 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-19 | f97101e1        | breakend_regulation_chr21_22_seed20260719_distance_137   | breakend_distance_137      |   24970 |       24970 |          0 |      24970 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-19 | f97101e1        | breakend_regulation_chr21_22_seed20260719_distance_5000  | breakend_distance_5000     |   29304 |       29304 |          0 |      29304 |                   0 | 100.00%    | 0.01%                                 |
| 2026-07-20 | e25c1513        | sv_confidence_grch38                                     | differential               |     466 |         466 |          0 |        466 |                   0 | 100.00%    | 0.79%                                 |
| 2026-07-20 | e25c1513        | witnesses                                                | differential               |     268 |         268 |          0 |        268 |                   0 | 100.00%    | 1.37%                                 |
| 2026-07-22 | 05620047        | state_exploration_seed_31415927                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-06 | a84ff150        | clinvar_chr21_hgvs_seed113                               | differential               |   56998 |       56998 |          0 |      56998 |                   0 | 100.00%    | 0.01%                                 |
| 2026-09-06 | a84ff150        | nmd_clinvar_chr21                                        | ensembl116-grch38-final    | 1353288 |     1353288 |          0 |    1353288 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-06 | b7c7237e        | clinvar_chr21_hgvs_seed113                               | differential               |   56998 |       56998 |          0 |      56998 |                   0 | 100.00%    | 0.01%                                 |
| 2026-09-06 | b7c7237e        | nmd_clinvar_chr21                                        | ensembl116-grch38-final    | 1353288 |     1353288 |          0 |    1353288 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 15417633        | state_exploration_seed_31415927                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | cc1993fd        | state_exploration_seed_31415927                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 7d40756a        | state_exploration_seed_16180339                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 7d40756a        | state_exploration_seed_27182818                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 6ce2ddd8        | clinvar_chr21_hgvs_seed113                               | differential               |   56998 |       56998 |          0 |      56998 |                   0 | 100.00%    | 0.01%                                 |
| 2026-09-07 | 6ce2ddd8        | hgvs_terminal_multiplication                             | differential               |       4 |           4 |          0 |          4 |                   0 | 100.00%    | 60.24%                                |
| 2026-09-07 | 6ce2ddd8        | state_exploration_seed_16180339                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 6ce2ddd8        | state_exploration_seed_27182818                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |
| 2026-09-07 | 6ce2ddd8        | state_exploration_seed_31415927                          | differential               |  100268 |      100268 |          0 |     100268 |                   0 | 100.00%    | 0.00%                                 |

## Randomized executable-VEP state exploration

This is the anti-overfitting lane against the VEP executable, not an
internal property test. Each seed contains the fixed predicate witnesses
plus 100,000 unique alleles. Three quarters of random positions are
within six bases of splice sites, exon endpoints, and translation starts
and stops; one quarter is uniform across the transcript. SNVs, MNVs,
insertions, deletions, and delins are sampled with equal probability,
with differing alleles up to 49 bases.

| revision | seed     | pairs   | exact   | unresolved | resolved_disagreements | descriptive_independent_pair_upper_95_ppm |
|:---------|:---------|:--------|:--------|:-----------|:-----------------------|:------------------------------------------|
| 6ce2ddd8 | 16180339 | 100,268 | 100,268 | 0          | 0                      | 36.8                                      |
| 6ce2ddd8 | 27182818 | 100,268 | 100,268 | 0          | 0                      | 36.8                                      |
| 6ce2ddd8 | 31415927 | 100,268 | 100,268 | 0          | 0                      | 36.8                                      |
| 6ce2ddd8 | combined | 300,804 | 300,804 | 0          | 0                      | 12.3                                      |

The combined denominator counts pair comparisons across seed runs, not
distinct alleles: the 268 fixed witnesses are deliberately shared and
random draws may overlap. Neither the original generator nor its
acceptance rules were changed.

### Fresh-seed counterexample retained

Seed 27182818 at `7d40756adc75` matched every consequence pair but
emitted one extra HGVSc: `chrDuck:250 CGT>CCC`, transcript `DUCK1-201`,
yielded `c.*10[3]` where VEP emitted no HGVSc. The pure-C properties
passed on that revision too. VEP skips transcript allele clipping only
for two-copy duplication; larger multiplications must undergo clipping
and coordinate projection. The fix removes the early repeat-formatting
path instead of changing the oracle or excluding the event. The original
failing HGVS rows remain in the append-only ledger.

| revision |  pairs | match | both_absent | discordant |
|:---------|-------:|------:|------------:|-----------:|
| 7d40756a | 100268 | 99135 |        1132 |          1 |
| 6ce2ddd8 | 100268 | 99135 |        1133 |          0 |

The same three frozen 100,268-pair corpora pass consequences and HGVS at
`6ce2ddd85df7`. A separate four-record executable witness retains the
discovered event, its two-copy positive control and larger-copy
variants. Native regression tests additionally enumerate 120
transcript-end projection cases: four bases, both strands and copy
counts 2 through 16. Those finite cases exercise the shared edit/HGVS
path; they are not 120 additional executable-VEP comparisons or proof of
whole-haplotype SO/HGVS semantics.

The same campaign covered the following SO terms. Counts are term
memberships, not distinct transcript pairs, because one pair may carry
several terms and seed runs deliberately retain the same fixed
witnesses.

| consequence_class                   | impact   | seeds_observed |      n | unresolved | term_mismatch | engine_extra | engine_missing |
|:------------------------------------|:---------|---------------:|-------:|-----------:|--------------:|-------------:|---------------:|
| intron_variant                      | MODIFIER |              3 | 116484 |          0 |             0 |            0 |              0 |
| coding_sequence_variant             | MODIFIER |              3 |  76854 |          0 |             0 |            0 |              0 |
| frameshift_variant                  | HIGH     |              3 |  49155 |          0 |             0 |            0 |              0 |
| 5_prime_UTR_variant                 | MODIFIER |              3 |  44738 |          0 |             0 |            0 |              0 |
| 3_prime_UTR_variant                 | MODIFIER |              3 |  41181 |          0 |             0 |            0 |              0 |
| splice_donor_variant                | HIGH     |              3 |  38750 |          0 |             0 |            0 |              0 |
| splice_acceptor_variant             | HIGH     |              3 |  34670 |          0 |             0 |            0 |              0 |
| splice_donor_5th_base_variant       | LOW      |              3 |  33786 |          0 |             0 |            0 |              0 |
| splice_polypyrimidine_tract_variant | LOW      |              3 |  31101 |          0 |             0 |            0 |              0 |
| stop_gained                         | HIGH     |              3 |  25890 |          0 |             0 |            0 |              0 |
| start_lost                          | HIGH     |              3 |  24495 |          0 |             0 |            0 |              0 |
| splice_region_variant               | LOW      |              3 |  24208 |          0 |             0 |            0 |              0 |
| splice_donor_region_variant         | LOW      |              3 |  14109 |          0 |             0 |            0 |              0 |
| missense_variant                    | MODERATE |              3 |  13035 |          0 |             0 |            0 |              0 |
| inframe_insertion                   | MODERATE |              3 |  12993 |          0 |             0 |            0 |              0 |
| stop_lost                           | HIGH     |              3 |  10983 |          0 |             0 |            0 |              0 |
| protein_altering_variant            | MODERATE |              3 |   8157 |          0 |             0 |            0 |              0 |
| stop_retained_variant               | LOW      |              3 |   7289 |          0 |             0 |            0 |              0 |
| inframe_deletion                    | MODERATE |              3 |   1237 |          0 |             0 |            0 |              0 |
| downstream_gene_variant             | MODIFIER |              3 |   1227 |          0 |             0 |            0 |              0 |
| start_retained_variant              | LOW      |              3 |   1051 |          0 |             0 |            0 |              0 |
| synonymous_variant                  | LOW      |              3 |    325 |          0 |             0 |            0 |              0 |
| intergenic_variant                  | MODIFIER |              3 |     52 |          0 |             0 |            0 |              0 |

This distribution deliberately stresses local alleles and positions near
exon, splice-site, and CDS endpoints on one engineered transcript. It
does not replace the indexed-cache corpora, which add real transcript
density, imported flags, exceptional peptide edits, codon tables,
assemblies, and species. The revision is shown explicitly because this
expensive campaign is not silently attributed to later code.

## Official Ensembl release corpus in Parquet

The official release consequence VCF is already BGZF-compressed. This
table measures its complete typed DuckHTS reader projection and the
narrower `VE` plus CSQ projection used by the bulk oracle lane. It is a
storage comparison, not a claim that the Parquet projection can
reproduce the original VCF byte-for-byte.

The consequence projection is also the natural CI payload: retain
deterministic shards with source URL, Ensembl release/species/assembly,
artifact digest, row cardinalities, and CSQ schema. A scheduled
full-release job may pair the complete projection with a published
receipt-hashed DuckDB model. Neither form broadens the supported
consequence contract; it only makes a large known-variant regression
cheap to replay.

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
| 2026-07-19 | 3feb3bf         | 0x0000000001352770 |                 45 | 4,400,500  | 4,400,500  |      0 |          0 |         180 | 206,342          |                29.515 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-19 | 5778e2b         | 0x000000000135276f |                 44 | 4,300,500  | 4,300,500  |      0 |          0 |         176 | 204,772          |                27.995 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-19 | f97101e         | 0x000000000135276f |                 44 | 4,300,500  | 4,300,500  |      0 |          0 |         176 | 204,781          |                28.512 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-20 | 0714235a        | 0x0000000001352770 |                 49 | 4,800,500  | 4,800,500  |      0 |          0 |         204 | 206,671          |                27.745 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-20 | 7dae50cd        | 0x0000000001352770 |                 49 | 4,800,500  | 4,800,500  |      0 |          0 |         206 | 206,710          |                27.654 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-20 | e25c1513        | 0x0000000001352770 |                 51 | 5,000,500  | 5,000,500  |      0 |          0 |         209 | 208,879          |                40.954 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-22 | 05620047        | 0x0000000001df5e77 |                 51 | 5,000,500  | 5,000,500  |      0 |          0 |         212 | 209,576          |                50.925 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-07-22 | 6eebf9b0        | 0x6a09e667f3bcc909 |                 52 | 5,100,500  | 5,100,500  |      0 |          0 |         214 | 211,624          |                41.131 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 15417633        | 0x0000000001df5e77 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         253 | 27,529,678       |                45.867 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 6ce2ddd8        | 0x0000000000f6e473 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         254 | 27,530,551       |                43.970 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 6ce2ddd8        | 0x00000000019ec6e2 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         254 | 27,529,725       |                54.879 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 6ce2ddd8        | 0x0000000001df5e77 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         254 | 27,530,051       |                43.830 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 7d40756a        | 0x0000000000f6e473 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         253 | 27,530,178       |                43.737 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | 7d40756a        | 0x00000000019ec6e2 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         253 | 27,529,352       |                43.699 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-07 | cc1993fd        | 0x0000000001df5e77 |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         251 | 27,322,306       |                44.148 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | 47eaa7b1        | 0x00000000000000ad |                 56 | 5,600,000  | 5,600,000  |      0 |          0 |         265 | 27,739,311       |                46.712 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | 47eaa7b1        | 0x000000000135282a |                 56 | 5,600,000  | 5,600,000  |      0 |          0 |         265 | 27,739,627       |                46.705 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | 6bcff335        | 0x00000000000000ad |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         262 | 27,738,016       |                47.299 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | 6bcff335        | 0x000000000135282a |                 55 | 5,500,000  | 5,500,000  |      0 |          0 |         262 | 27,738,332       |                47.452 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | ae131a9d        | 0x00000000000000ad |                 56 | 5,600,000  | 5,600,000  |      0 |          0 |         269 | 34,092,011       |                46.210 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| 2026-09-08 | ae131a9d        | 0x000000000135282a |                 56 | 5,600,000  | 5,600,000  |      0 |          0 |         269 | 34,092,327       |                45.861 | cc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |

| target                                                                          | trials  | passed  | failed | skipped | duplicates |
|:--------------------------------------------------------------------------------|:--------|:--------|:-------|:--------|:-----------|
| annotate cursor cross-codon MNV route == tile                                   | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor DEL route == tile under output splits                           | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor INS route == tile under output splits                           | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor output splits == one annotate_tile                              | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate cursor padded SNV == tile under output splits                          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile == sweep + classify + structural-SO composition                   | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon refinement == coding-SNV kernel oracle                      | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon-aligned in-frame deletion == CDS-position oracle            | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile codon-boundary in-frame insertion == CDS-position oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile non-boundary in-frame insertion == peptide-window oracle          | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile rejects NULL model without reading the batch                      | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile same-codon MNV == codon oracle                                    | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile simple frameshift indel == CDS-position oracle                    | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile start_lost SNV == start-codon oracle                              | 100,000 | 100,000 | 0      | 0       | 0          |
| annotate_tile two-codon body MNV missense == codon-window oracle                | 100,000 | 100,000 | 0      | 0       | 0          |
| breakend_parser_recovers_constructed_components                                 | 100,000 | 100,000 | 0      | 0       | 0          |
| cgranges-seeded first event + sweep == brute-force candidates                   | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context == direct CDS splice + full peptide oracles                      | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delins shape == local-edge oracle                                | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta == single-codon oracle                                     | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame deletion == edit-origin oracle                    | 100,000 | 100,000 | 0      | 0       | 0          |
| coding context delta in-frame insertion == edit-origin oracle                   | 100,000 | 100,000 | 0      | 0       | 0          |
| codon change classification consistent with translation                         | 100,000 | 100,000 | 0      | 0       | 0          |
| complete literal spans == VEP complete-overlap source semantics                 | 100,000 | 100,000 | 0      | 0       | 0          |
| coordinate projection == brute-force transcript-order base walk                 | 100,000 | 100,000 | 0      | 0       | 0          |
| event differing-region normalization == independent trim oracle                 | 100,000 | 100,000 | 0      | 0       | 0          |
| haplotype block spans reconstruct the independently replayed CDS                | 100,000 | 100,000 | 0      | 0       | 0          |
| haplotype blocks preserve every frame and same-codon interaction                | 100,000 | 100,000 | 0      | 0       | 0          |
| HGVS genomic 3-prime shift == independent reference byte-walk                   | 100,000 | 100,000 | 0      | 0       | 0          |
| HGVSp fact replay == independently translated edited CDS                        | 100,000 | 100,000 | 0      | 0       | 0          |
| HGVSp frameshift fact == independently extended translation                     | 100,000 | 100,000 | 0      | 0       | 0          |
| multi-edit CDS haplotype apply == left-to-right rebuild oracle                  | 100,000 | 100,000 | 0      | 0       | 0          |
| optimized sorted annotation == forced generalized full rows                     | 100,000 | 100,000 | 0      | 0       | 0          |
| ordered source replacements == literal replay, net spans and applied provenance | 100,000 | 100,000 | 0      | 0       | 0          |
| owned haplotype replay == dense genomic edits in coexisting models              | 100,000 | 100,000 | 0      | 0       | 0          |
| phased SNV set == equivalent MNV coding facts                                   | 100,000 | 100,000 | 0      | 0       | 0          |
| region mask structural invariants                                               | 100,000 | 100,000 | 0      | 0       | 0          |
| regulation sweep/BND pairs == independent feature oracles                       | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta annotation wrapper MNV == direct shape                           | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta exon hint == unhinted projection                                 | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch INDEL == local delins-shape oracle                       | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch MNV == single-codon oracle                               | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence delta scratch two-codon MNV window == codon-window oracle              | 100,000 | 100,000 | 0      | 0       | 0          |
| sequence-backed SNV codon edit == codon-slice edit oracle                       | 100,000 | 100,000 | 0      | 0       | 0          |
| simple indel route == generalized CodingContext                                 | 100,000 | 100,000 | 0      | 0       | 0          |
| sorted point cursor classifier == exhaustive exon/gap scans                     | 100,000 | 100,000 | 0      | 0       | 0          |
| sorted span cursor classifier == exhaustive exon/gap scans                      | 100,000 | 100,000 | 0      | 0       | 0          |
| sparse carrier paths == dense event matrix across input batches                 | 100,000 | 100,000 | 0      | 0       | 0          |
| sweep candidate set == brute-force candidate set                                | 100,000 | 100,000 | 0      | 0       | 0          |
| terminal partial-codon insertion == codon-rounded VEP translation oracle        | 100,000 | 100,000 | 0      | 0       | 0          |
| transcript coordinate == brute-force exon/intron walk                           | 100,000 | 100,000 | 0      | 0       | 0          |
| variant CDS edit builder == direct CDS splice oracle                            | 100,000 | 100,000 | 0      | 0       | 0          |
| variant CDS edit-set builder == single-edit splice oracle                       | 100,000 | 100,000 | 0      | 0       | 0          |
| variant CDS edit-set builder splits MNV diff islands                            | 100,000 | 100,000 | 0      | 0       | 0          |
| variant coding context == direct CDS splice + full peptide oracles              | 100,000 | 100,000 | 0      | 0       | 0          |
| VEP feature-span sweep candidates == independent parser oracle                  | 100,000 | 100,000 | 0      | 0       | 0          |

Passing the requested number of trials is necessary but does not prove
that a generator visited the states named by its contract. Randomized
properties therefore emit distribution counters, and the recorder stores
each counter as a separate numeric row. The table below is the latest
complete run’s state distribution; the long-form CSV remains the
machine-readable authority. Zero is evidence too: it identifies a state
that the declared seed did not exercise and must not be hidden by the
suite-level pass count.

| randomized distribution             | observed states                                                                                                                                                                                                                                                                                    |
|:------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| allele sweep coverage               | del= 812,847; indel= 811,257; ins= 810,886; interbase= 901,094; mnv= 811,941; prefix= 2,841,350; suffix= 2,247,130; tail= 2,683,865                                                                                                                                                                |
| annotation-shortcut coverage        | coding_tx= 641,166; cursor_splits= 100,000; far= 4,940,445; generalized=12,283,823; mirna_tx= 319,849; nmd_rows= 2,281,720; simple= 1,709,247                                                                                                                                                      |
| cds-edit-builder coverage           | body= 41,654; del= 19,940; fwd= 49,955; indel= 19,988; ins= 20,031; mnv= 20,094; rev= 50,045; snv= 19,947; start= 29,239; stop= 29,107                                                                                                                                                             |
| cds-edit-set coverage               | body= 41,654; cap0= 100,000; del= 19,940; fwd= 49,955; indel= 19,988; ins= 20,031; mnv= 20,094; rev= 50,045; snv= 19,947; start= 29,239; stop= 29,107                                                                                                                                              |
| cds-edit-set-mnv coverage           | body= 33,177; capfail= 100,000; fwd= 50,079; multi= 100,000; rev= 49,921; start= 33,617; stop= 33,206                                                                                                                                                                                              |
| coding-context coverage             | capfail= 300,000; del= 19,940; fwd= 49,955; indel= 19,988; ins= 20,031; mnv= 20,094; pep_diff= 86,666; pep_same= 13,334; rev= 50,045; snv= 19,947                                                                                                                                                  |
| codon coverage                      | mis= 68,020; stop_gained= 3,909; stop_lost= 3,974; stop_retained= 673; syn= 23,424                                                                                                                                                                                                                 |
| complete-overlap coverage           | forward= 49,807; over_5000= 12,415; reverse= 50,193; right_endpoint= 1,465                                                                                                                                                                                                                         |
| context-delins-shape coverage       | forward= 50,151; inframe= 49,984; lengthen= 49,864; protein_altering= 50,016; reverse= 49,849; shorten= 50,136                                                                                                                                                                                     |
| context-delta coverage              | fwd= 50,015; mis= 20,049; rev= 49,985; stop_gained= 19,948; stop_lost= 20,007; stop_retained= 19,902; syn= 20,094                                                                                                                                                                                  |
| context-inframe-deletion coverage   | forward= 50,130; reverse= 49,870                                                                                                                                                                                                                                                                   |
| context-inframe-insertion coverage  | forward= 50,050; reverse= 49,950                                                                                                                                                                                                                                                                   |
| cross-mnv coverage                  | fwd= 50,132; len2= 50,187; len3= 49,813; missense= 50,032; rev= 49,868; stop_gained= 24,898; synonymous= 25,070                                                                                                                                                                                    |
| cursor-cross-route coverage         | context= 100,000; fwd= 50,132; len2= 50,187; len3= 49,813; rev= 49,868                                                                                                                                                                                                                             |
| cursor-del-route coverage           | forward= 50,130; full= 100,000; reverse= 49,870                                                                                                                                                                                                                                                    |
| cursor-ins-route coverage           | forward= 50,050; full= 100,000; reverse= 49,950                                                                                                                                                                                                                                                    |
| cursor-route coverage               | full= 100,000; fwd= 50,015; mis= 20,049; rev= 49,985; stop_gained= 19,948; stop_lost= 20,007; stop_retained= 19,902; syn= 20,094                                                                                                                                                                   |
| delta-cross-scratch coverage        | fwd= 50,132; len2= 50,187; len3= 49,813; missense= 50,032; rev= 49,868; stop_gained= 24,898; synonymous= 25,070                                                                                                                                                                                    |
| delta-exon-hint coverage            | del= 19,940; fwd= 49,955; indel= 19,988; ins= 20,031; mnv= 20,094; rev= 50,045; snv= 19,947                                                                                                                                                                                                        |
| delta-scratch coverage              | capfail= 100,000; fwd= 50,015; mis= 20,049; rev= 49,985; stop_gained= 19,948; stop_lost= 20,007; stop_retained= 19,902; syn= 20,094                                                                                                                                                                |
| delta-scratch-indel coverage        | forward= 50,151; lengthen= 49,864; reverse= 49,849; shorten= 50,136                                                                                                                                                                                                                                |
| delta-wrapper coverage              | fwd= 50,015; mis= 20,049; rev= 49,985; stop_gained= 19,948; stop_lost= 20,007; stop_retained= 19,902; syn= 20,094                                                                                                                                                                                  |
| event normalization coverage        | del= 25,028; indel= 24,718; ins= 24,954; interbase= 27,720; prefix= 81,125; prefix0_interbase= 7,019; sub= 25,300; suffix= 77,968                                                                                                                                                                  |
| frameshift coverage                 | -1= 7,294; -2= 7,155; +1= 7,530; +2= 7,131; del= 29,432; delins= 29,110; ins= 41,458; reverse= 47,774; stop_gained= 1,562; terminal_cil_protein_altering= 9; terminal_cil_retained= 47; terminal_endpoint= 16,522; terminal_missing_tail= 4,182; terminal_nonstop= 12,362; terminal_reverse= 8,203 |
| frameshift length-oracle coverage   | frameshift= 62,587; inframe_len= 12,446; stop_gained= 1,893                                                                                                                                                                                                                                        |
| haplotype-MNV equivalence coverage  | body= 33,177; fwd= 50,079; one_codon= 16,102; rev= 49,921; several_codons= 83,898; start= 33,617; stop= 33,206                                                                                                                                                                                     |
| HGVS shift coverage                 | at_vep_limit= 0; composed= 96,683; del= 48,656; dup= 37,932; fwd= 49,955; ins= 51,189; nonlocal_ref_replay= 1,542; protein= 90,059; rev= 49,890; rotated= 20,818; terminal_duplication= 3,162                                                                                                      |
| HGVSp frameshift coverage           | del= 24,911; delins= 29,221; eligible= 91,789; equal_stop= 19; fs= 87,367; fwd= 45,781; immediate_stop= 4,403; ins= 37,657; non_fs= 8,211; rev= 46,008; shortened= 0; ter_known= 17,365; ter_unknown= 70,002                                                                                       |
| HGVSp replay coverage               | del= 3,941; delins= 3,844; dup= 388; equal= 8,346; fwd= 17,190; ins= 1,582; replayed= 34,988; rev= 17,798; special= 63,926; sub= 16,887; terminal_not_applicable= 1,086; vep_position_zero= 1,324; vep_stop_equal= 0                                                                               |
| inframe_deletion coverage           | forward= 50,130; reverse= 49,870                                                                                                                                                                                                                                                                   |
| inframe_insertion coverage          | forward= 50,050; reverse= 49,950                                                                                                                                                                                                                                                                   |
| mnv coverage                        | len2= 50,110; len3= 49,890                                                                                                                                                                                                                                                                         |
| non-boundary insertion coverage     | forward= 50,132; inframe_insertion= 49,985; protein_altering= 50,015; reverse= 49,868                                                                                                                                                                                                              |
| ordered-replacement coverage        | clipped= 85,602; empty= 6,365; forward= 49,777; merged_sources= 135,000; noop= 102,730; ref_slot= 220,687; reverse= 50,223; tied= 519,063                                                                                                                                                          |
| simple-indel equivalence coverage   | del= 11,394; delins= 5,652; fallback= 38,829; fast= 21,130; frameshift= 18,801; fwd= 10,366; inframe_del= 1,847; inframe_ins= 482; ins= 4,084; rev= 10,764                                                                                                                                         |
| start-codon coverage                | co_stop_gained= 4,129; co_synonymous= 24,245; lost_and_retained= 1,596; start_lost= 100,000; start_retained= 1,596; synonymous= 24,245                                                                                                                                                             |
| terminal-partial-insertion coverage | after_tail_rejected= 41,724; length_mod0= 33,507; length_mod1= 33,297; length_mod2= 33,196; mitochondrial= 50,038; nonstop= 76,796; reverse_orientation= 50,001; same_orientation= 49,999; site_first= 41,607; site_internal= 16,669; standard= 49,962; stop= 11,190; tail1= 49,970; tail2= 50,030 |
| variant-coding-context coverage     | capfail= 400,000; del= 19,940; fwd= 49,955; indel= 19,988; ins= 20,031; mnv= 20,094; pep_diff= 86,666; pep_same= 13,334; rev= 50,045; snv= 19,947                                                                                                                                                  |

The run observed all 263 required nonzero counters. The other 3 counters
have named fixed witnesses in the coverage manifest; their absence from
a random draw is not counted as statistical coverage. These counters
describe the declared generators, not an exhaustive enumeration of
biological configurations.

| revision | seed               | targets |  trials | required_counters_observed | minimum_required_counter_hits | fixed_witness_counters_not_hit |
|:---------|:-------------------|--------:|--------:|---------------------------:|------------------------------:|-------------------------------:|
| ae131a9d | 0x00000000000000ad |      56 | 5600000 |                        263 |                             6 |                              3 |
| ae131a9d | 0x000000000135282a |      56 | 5600000 |                        263 |                             9 |                              3 |

Together these seeds executed 11,200,000 property trials on the shown
revision. Millions of passing trials do not make a counter with
single-digit observations densely explored, and marginal counters do not
establish coverage of their cross-products. The fresh-seed HGVS failure
above is direct evidence of this limit. Dedicated rare-state strata and
retained counterexamples complement broad draws; they do not justify a
population error-rate claim. Full phased SO/HGVS, broader structural
composition and stale real-corpus campaigns still need their own
current-revision evidence.

## Phased replay with noncoding contributors

| source_revision |     seed | policy        | input_records | input_calls | input_allele_slots | observed_carriers | provenance_memberships |
|:----------------|---------:|:--------------|--------------:|------------:|-------------------:|------------------:|-----------------------:|
| eb83f6ff        |      173 | strict        |          3764 |       11292 |              22584 |              6000 |                  14292 |
| eb83f6ff        |      173 | vep116_compat |          3764 |       11292 |              22584 |              6000 |                  14292 |
| eb83f6ff        | 20260906 | strict        |          3802 |       11406 |              22812 |              6000 |                  14406 |
| eb83f6ff        | 20260906 | vep116_compat |          3802 |       11406 |              22812 |              6000 |                  14406 |
| d1c591b7        |      173 | strict        |          3764 |       11292 |              22584 |              6000 |                  14292 |
| d1c591b7        |      173 | vep116_compat |          3764 |       11292 |              22584 |              6000 |                  14292 |
| d1c591b7        | 20260906 | strict        |          3802 |       11406 |              22812 |              6000 |                  14406 |
| d1c591b7        | 20260906 | vep116_compat |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 8f9987e3        |      173 | strict        |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 8f9987e3        |      173 | vep116_compat |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 8f9987e3        | 20260906 | strict        |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 8f9987e3        | 20260906 | vep116_compat |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 7f4a4e28        |      173 | strict        |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 7f4a4e28        |      173 | vep116_compat |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 7f4a4e28        | 20260906 | strict        |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 7f4a4e28        | 20260906 | vep116_compat |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 9b1b8a0d        |      173 | strict        |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 9b1b8a0d        |      173 | vep116_compat |          3764 |       11292 |              22584 |              6000 |                  14292 |
| 9b1b8a0d        | 20260906 | strict        |          3802 |       11406 |              22812 |              6000 |                  14406 |
| 9b1b8a0d        | 20260906 | vep116_compat |          3802 |       11406 |              22812 |              6000 |                  14406 |

Sources eb83f6ff1d03d05a3c9f8135c8ef355b7f431ee7,
d1c591b76f8a9a07036736ac0666a004eb58e0eb,
8f9987e3826018cfa73c155eebac7a956b6dc024,
7f4a4e28bff31a13f14ee1cab25049408740ef65,
9b1b8a0d4e12c99f8f08d9c013d840abf5eb7ec1 were built from clean
checkouts, including an HTSlib clean rebuild. The ledger retains the
extension hash, input/run receipt hashes and pinned VEP/variation
revisions. DuckDB used four threads; these are correctness counts, not
timing or memory measurements.

Each unchanged 1,000-transcript corpus first passes its original public
replay checks: 6,000 complete lanes, 4,000 occupied carriers, 3,000
output leaves, 22,000 oracle comparisons and 4,000 first-stop/frame
comparisons per policy. The supplemental corpus adds one homozygous
intronic SNV per transcript while preserving every original VCF record.
Running unmodified Haplosaurus on those inputs leaves its complete
observations unchanged. DuckHTS must retain the intronic contributors,
including on previously implicit reference lanes, with unchanged literal
CDS/protein and an `outside_cds` contributor status.

Both policies pass all 120,000 carrier comparisons and 286,980
provenance memberships. The same biological lanes are counted separately
under each policy and revision; these are not independent statistical
trials. Five deliberately corrupted outputs per policy/seed are
rejected. Fixed SQL/R tests additionally cover UTRs, insertions, missing
calls and coding-overlapping projection failures; a native two-strand
span enumeration supplies 6,774 assertions for the fix.

This certifies the declared literal-replay cases, not altered splicing,
combined SO/HGVS, broad phase compatibility or exhaustive rare
configurations. The [conformance
driver](../test/duckvep/conformance/haplotype_sql_differential.R) keeps
this augmentation opt-in and does not replace the original corpus or
verifier. The [phased replay benchmark](duckvep_haplotypes.md) records
sorted native and public SQL execution separately, with workspace and
process memory. The current SQL benchmark includes local coding-block
SO; these Haplosaurus comparisons do not certify those masks.
Whole-haplotype SO/HGVS remains unfinished.

## Raw genotype compatibility audit

| ploidy | cases | disagreements | oracle_lanes | native_lanes | native_unavailable_carriers |
|-------:|------:|--------------:|-------------:|-------------:|----------------------------:|
|      1 |    12 |            12 |           24 |           12 |                           3 |
|      2 |    96 |            66 |          192 |          192 |                          84 |
|      3 |   768 |           768 |         1536 |         2304 |                        1332 |
|      4 |  6144 |          6144 |        12288 |        24576 |                       16800 |

Source ae131a9d08f8cc3e8b38c034337740db0da8dbd4 records a **failing
decoded-call/raw-parser comparison**: 6990 disagreements in 7020
profiles. It is not a population error rate or a replacement for the
passing literal-sequence corpus. This lane uses public
`input_mode := 'alt_events'`; the separate raw source-record input is
checked below.

The [R
driver](../test/duckvep/conformance/haplotype_phase_differential.R)
enumerates every GT over `0`, `1`, `2`, `.` at ploidies 1–4, every
intervening separator pattern, and absent, `/`, or `|` leading prefixes.
No profile is sampled or excluded. Each has one multiallelic site plus a
homozygous second site in a different phase set, using the registered
180-base reference. There are 14,040 source records/genotype calls,
21,060 source ALT events/candidate calls, and 54,168 input allele slots.
Every REF is checked before execution. These are correctness
denominators, not timing measurements.

The pinned, unmodified Haplosaurus runner and decoded `vep116_compat`
executor consume the same VCF/GFF/FASTA. Comparisons retain complete
CDS/protein multisets, source-record contributors and carrier counts.
Eighteen ordinary called diploid profiles agree; four deliberate
sequence/protein/provenance/duplicate corruptions are rejected. All
other disagreements remain failures, including missing-input NULL
sequences and the difference between explicit source ploidy and
Haplosaurus’s file-input diploid fallback. The command exits nonzero
after writing full observations, comparisons and source-bound receipts.

There are 1471 groups in which distinct raw GT spellings have
**identical HTSlib alleles and phase flags but different Haplosaurus
outputs**. A witness is `0|1` versus `|0|1`: both decode to alleles
`[0,1]`, phase flags `[true,true]`. VEP-116’s parser retains a leading
empty split field as REF before the container consumes its allele slots,
changing the result. The same model, remaining call and phase sets are
used on both sides.

Thus raw-parser compatibility cannot be recovered from typed calls
alone. This does not justify changing HTSlib-faithful genotype decoding
or treating an unknown call as biologically known. Exact raw-input
emulation needs retained source GT and source-record allele context,
with upstream conditional sequence explicitly distinguished from
strict-phase evidence. Whole-haplotype SO/HGVS and typed structural
composition remain separate unfinished requirements.

The constant-space native raw-GT parser separately records **0
disagreements across 14040 source calls**. An optional observer sidecar
reads actual Haplosaurus genotype objects and its file-profile ploidy;
it does not override parsing or sequence construction. The comparison
checks retained/omitted calls, parsed slot counts, the two consumed
allele ordinals, source ploidy and missingness, with seven rejected
field corruptions. The standard output on the 1,000-transcript seed-173
corpus matches its recorded oracle output. Parser code, bridge, compiler
identity, binary and observations are hashed in the same clean-build
receipt.

The native source-record stream records **0 disagreements across 7020
complete CDS/protein, carrier-count and physical-edit provenance
comparisons**, and **0 disagreements across 28080 per-lane source-record
observations**. The latter checks selected allele ordinals, record and
aggregate evidence, and sequence availability against the upstream
object sidecar; four deliberate field corruptions are rejected.

Raw replay keeps full source REF spans and actual ALT ordinals. An
undefined consumed slot follows the pinned runner’s full-REF deletion,
with `conditional` sequence status and evidence distinct from a called
ALT. Missing REF observations retain provenance without inventing
physical edits. Reference validation and coding projection failures
still prevent sequence availability.

This native result covers the enumerated two-site, single-exon grammar,
not overlapping replacements, splicing or combined SO/HGVS. Public
decoded-call comparison objects, including all 6990 failures, match the
preceding audit exactly. The passing native lane does not waive those
failures or make conditional sequence biologically known.

Public `input_mode := 'source_records'` with
`phase_policy := 'vep116_compat'` records **0 disagreements across 7020
complete CDS/protein multiset, carrier-count and physical-edit
provenance comparisons**, plus **0 disagreements across 28080 per-lane
record observations**. It consumes original GT text and complete ALT
lists from the same VCF fixture. The verifier additionally checks exact
global record IDs, regions, positions and source REF/ALT bytes; repeated
local site labels cannot substitute for record identity.

The public lane exercises input materialization, sorting, transcript
closing, DuckDB vector transitions and nested output. Fixed SQL/R
controls reject duplicate calls, NULL keys, inconsistent record/GT
identities, invalid ALT lists and GTs, and ploidy-limit exhaustion. This
establishes the declared finite raw-input contract, not overlapping
replacements, broader phase conformance, whole-haplotype SO/HGVS or
biological certainty for conditional sequence.

## Rare source-record geometry

|     seed | cohort           | profiles | full failures | available sequence differences | profiles with unavailable sequence | carrier-total failures |
|---------:|:-----------------|---------:|--------------:|-------------------------------:|-----------------------------------:|-----------------------:|
|      173 | fixed_and_random |      656 |             0 |                              0 |                                  0 |                      0 |
| 20260906 | fixed_and_random |      656 |             0 |                              0 |                                  0 |                      0 |
|      173 | rare             |    36096 |             0 |                              0 |                                  0 |                      0 |
| 20260906 | rare             |    36096 |             0 |                              0 |                                  0 |                      0 |

Source ae131a9d08f8cc3e8b38c034337740db0da8dbd4 runs the public
`source_records`/`vep116_compat` path against pinned, unmodified
Haplosaurus. Across the two seeds there are **73,504 profiles**, 220,512
source records and 147,008 oracle file lanes; DuckHTS returns 147,008
carrier memberships. The full comparison retains CDS/protein multisets,
counts and applied-record identity sets within each equal-sequence
group. It does not certify physical-edit multiplicity or per-lane
sequence association. Carrier-key uniqueness and totals are checked
separately. Available includes explicitly conditional replay; a NULL
sequence remains a disagreement, not a waiver. All failure columns count
profiles; 0 profiles have unavailable sequence. Availability does not
imply agreement. The receipt ledger retains every run, generator,
coverage and artifact hash.

| source_revision                          | profiles | failures | available_sequence_failures | profiles_with_unavailable |
|:-----------------------------------------|---------:|---------:|----------------------------:|--------------------------:|
| 47eaa7b1460a6b29049cdb50928956a055d6abb9 |    73504 |    15953 |                       14694 |                         0 |
| ae131a9d08f8cc3e8b38c034337740db0da8dbd4 |    73504 |        0 |                           0 |                         0 |

|     seed | still passing | resolved | regressed | still failing |
|---------:|--------------:|---------:|----------:|--------------:|
|      173 |         28807 |     7945 |         0 |             0 |
| 20260906 |         28744 |     8008 |         0 |             0 |

The identical-input comparison resolves **15,953 retained failures**,
with **0 regressions** and **0 remaining failures**. These transitions
apply only to the declared generated inputs; they do not establish
absence of regression outside this campaign. Every earlier failing
comparison remains in its revision-labelled artifact.

Each seed includes 144 fixed profiles and 512 general-random overlapping
pairs. The rare lane requires **32 draws in each of 1,128 cells**: 12
geometries × two strands × 47 valid GT-pattern/source-ploidy
combinations. Ploidies are 1, 2, 4, 8, 16 and 64. Patterns include mixed
and leading separators, first/last missing slots, all-missing calls and
ALT calls beyond the two consumed file lanes. Positions, replacement
lengths and bases are seeded; actual GT strings and all cell quotas are
checked. A pair shares its GT-pattern/ploidy class; this does not
enumerate every pairing of different classes. Each profile has two
tested records and a separate homozygous anchor. Both policies’ existing
decoded-input corpora and the finite raw-GT audit remain independent
evidence lanes.

All 36,752 generated inputs per seed match source
47eaa7b1460a6b29049cdb50928956a055d6abb9 exactly, and each complete
oracle output is byte-identical. The comparison rules and all preceding
comparison artifacts are preserved. Four deliberate output corruptions
are rejected; the 24 fixed disjoint/adjacent controls per seed pass.
Fixed profiles repeat across seeds, and all profiles use the same
registered 180-base, single-exon reference and standard codon table.
These are quota-controlled generated configurations, not independent
biological observations, population error rates, strict-phase
certification, combined SO/HGVS or structural-event conformance. Four
DuckDB threads are used; these counts are not performance measurements.

### Full-span replacement witness

At source `6bcff3353d7c9a2609996aae9828f9da3fb5f238`, the fixed
positive-strand `mnv_retained_middle` cis profile has these records:

| Source | Genomic position | REF     | ALT     | GT   |
|--------|-----------------:|---------|---------|------|
| a      |               40 | TGCTGCT | AGCTGCA | 1\|0 |
| b      |               43 | T       | A       | 1\|0 |
| anchor |              165 | C       | G       | 1\|1 |

On the changed lane, CDS bases 25–40 are `GCTGCAGCTGCAGCTG` in
Haplosaurus. DuckHTS source `6bcff3353d7c9a2609996aae9828f9da3fb5f238`
returns `GCTGCAGCAGCAGCTG`; source
`47eaa7b1460a6b29049cdb50928956a055d6abb9` matches Haplosaurus’s
complete sequence and applied-record provenance. The pinned
[`_mutate_sequences`
implementation](https://github.com/Ensembl/ensembl-variation/blob/2fb834b987ede3824e200197a838ce11e91aeb4b/modules/Bio/EnsEMBL/Variation/TranscriptHaplotypeContainer.pm#L1160)
uses the [Haplo parser’s raw
spans](https://github.com/Ensembl/ensembl-vep/blob/57ea5c52340acc1f156267f810ad162e26597082/modules/Bio/EnsEMBL/VEP/Haplo/Parser/VCF.pm#L198)
and applies full REF spans in descending mapped-CDS order: record a
restores the retained middle base after b changes it. Applied-record
provenance retains b even when its sequence change is overwritten.
Ordered raw replacements use bounded native storage and preserve net
component spans; local SO and displaced-frame facts are unavailable for
that path. Source-buffer ordering and duplicate selection are exercised
by the quota-controlled campaigns in this report. Multi-transcript and
cross-exon mapping require separate conformance. These observations are
not evidence that VEP is wrong.

### Paired source-context configurations

|     seed | profiles | input_records | observed_carriers | failures | oracle_context_changed | observed_context_changed |
|---------:|---------:|--------------:|------------------:|---------:|-----------------------:|-------------------------:|
|      173 |    19664 |        345840 |             39328 |        0 |                   5376 |                     5376 |
| 20260906 |    19664 |        345840 |             39328 |        0 |                   5376 |                     5376 |

Source ae131a9d08f8cc3e8b38c034337740db0da8dbd4 includes **38,016 paired
context profiles** from 1152 seeded edit templates. The table also
includes 656 fixed/general-random controls per seed; their inputs and
complete oracle observations match the corresponding profiles in the
rare-GT campaign. The ledger retains both the one-draw smoke runs and
the 8-draw runs.

Every one of 2,376 geometry × GT-pattern × strand × placement ×
neutral-count cells receives 8 draws per seed. Each edit template is
reused in 33 contexts: reference-only records occur before, between or
after the tested edits, with counts including the neighbourhoods of the
pinned interval tree’s root changes. The generator verifies unchanged
edit geometry, alleles and GTs across those contexts, the complete
upstream source buffer, and the absence of neutral records from retained
genotype objects.

In **10,752 profiles**, the oracle’s complete
sequence/count/applied-record-set observation differs from the
zero-neutral-record case; DuckHTS records 10,752 such changes. The full
per-profile comparison, not equality of those two totals, determines
agreement. These are paired observations, not independent biological
samples or population error-rate estimates. The single-exon, one-sample
reference and comparison limitations above also apply here; rare GT
classes and neutral contexts are separate campaigns, not their complete
cross-product.

### Shared-transcript rare configurations

|     seed | profiles | source_records | observed_carriers | failures | count_failures | input_provenance_failures | mapping_failures |
|---------:|---------:|---------------:|------------------:|---------:|---------------:|--------------------------:|-----------------:|
|      173 |    18432 |         166224 |            110592 |     2304 |              0 |                         0 |                0 |
| 20260906 |    18432 |         166224 |            110592 |     2304 |              0 |                         0 |                0 |

Source 2c5bcb21ba6f054639b2031f8ca465be3b8b27d5 has **36,864 transcript
cases**, including 33,792 quota-controlled cases and the complete
3,072-case diagnostic baseline. Every region has three diploid samples,
one full-exon transcript and an overlapping two-exon transcript (exons
11–70 and 101–190). Source coordinates and spanning lengths vary;
transcript geometry and the registered 180-base reference are fixed.

Each seed supplies 32 draws in every geometry × cohort-GT-pattern ×
strand × neutral-count cell: four record geometries, three sample
configurations, two strands and eleven neutral source counts through 36.
This crosses shared-transcript mapping with source-context pressure; it
is not a distribution of independent biological samples, arbitrary
transcript models or an estimate of population error rates.

All **4,608 disagreements** are exon-spanning records in the two-exon
model. DuckHTS withholds the sequence, while pinned Haplosaurus filters
the unmapped record and replays the other edits. These are failing
compatibility cases, not exclusions or evidence that VEP is wrong. Every
baseline input column, complete oracle JSON line and sequence verdict
agrees with the retained diagnostic; its 384 disagreements remain in the
denominator.

Full CDS/protein groups include counts per sample; exchanging samples is
a rejected corruption even when pooled counts agree. Applied-source
identity sets are checked within equal-sequence groups. Native
contributor checks independently retain source IDs, regions, positions,
REF/ALT interpretations and occupied sample/file-lane keys. Upstream
checks retain the complete source buffer, genotype multiplicities and
constructor-owned mappings, including exon-repeated and unselected
duplicate sources. All thirteen corruption controls pass for each seed.
Source sets do not certify physical-edit multiplicity or
source-to-sample association within a shared upstream sequence group.
Full phase/PS behavior, compound SO/HGVS and structural composition
remain outside this campaign’s tested scope.

## Individual Sequence Ontology terms

For each transcript pair, this compares the union of terms emitted by
either engine. A missing or extra term is therefore visible under its
own SO name. Rows must not be summed across terms because one pair can
carry several terms.

| corpus                            | observed_terms | terms_with_mismatch | term_mismatches | engine_extra | engine_missing | terms_with_unresolved | unresolved_term_memberships |
|:----------------------------------|---------------:|--------------------:|----------------:|-------------:|---------------:|----------------------:|----------------------------:|
| GRCh38 dbSNP                      |             22 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 GIAB                       |             19 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 ClinVar coding             |             27 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 ClinVar cross-chromosome   |             28 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh37                            |             26 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| P. falciparum                     |             20 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 paired BND                 |             14 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 GIAB + core regulation     |             15 |                   0 |               0 |            0 |              0 |                     0 |                           0 |
| GRCh38 exact SV + core regulation |             27 |                   0 |               0 |            0 |              0 |                     0 |                           0 |

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

| corpus                            | impact   |      n | exact_agree | unresolved | resolved_n | resolved_discordant | exact_rate | descriptive_independent_pair_upper_95 |
|:----------------------------------|:---------|-------:|------------:|-----------:|-----------:|--------------------:|:-----------|:--------------------------------------|
| GRCh38 paired BND                 | HIGH     |  69654 |       69654 |          0 |      69654 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 ClinVar coding             | HIGH     |  99103 |       99103 |          0 |      99103 |                   0 | 100.00%    | 0.00%                                 |
| GRCh38 ClinVar cross-chromosome   | HIGH     |  79663 |       79663 |          0 |      79663 |                   0 | 100.00%    | 0.00%                                 |
| GRCh38 dbSNP                      | HIGH     |    118 |         118 |          0 |        118 |                   0 | 100.00%    | 3.08%                                 |
| GRCh38 GIAB                       | HIGH     |      4 |           4 |          0 |          4 |                   0 | 100.00%    | 60.24%                                |
| GRCh37                            | HIGH     |  46645 |       46645 |          0 |      46645 |                   0 | 100.00%    | 0.01%                                 |
| P. falciparum                     | HIGH     |   4309 |        4309 |          0 |       4309 |                   0 | 100.00%    | 0.09%                                 |
| GRCh38 GIAB + core regulation     | HIGH     |      3 |           3 |          0 |          3 |                   0 | 100.00%    | 70.76%                                |
| GRCh38 exact SV + core regulation | HIGH     |  40135 |       40135 |          0 |      40135 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 ClinVar coding             | LOW      |   5952 |        5952 |          0 |       5952 |                   0 | 100.00%    | 0.06%                                 |
| GRCh38 ClinVar cross-chromosome   | LOW      |  28776 |       28776 |          0 |      28776 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 dbSNP                      | LOW      |    349 |         349 |          0 |        349 |                   0 | 100.00%    | 1.05%                                 |
| GRCh38 GIAB                       | LOW      |    260 |         260 |          0 |        260 |                   0 | 100.00%    | 1.41%                                 |
| GRCh37                            | LOW      |  23354 |       23354 |          0 |      23354 |                   0 | 100.00%    | 0.02%                                 |
| P. falciparum                     | LOW      |    210 |         210 |          0 |        210 |                   0 | 100.00%    | 1.74%                                 |
| GRCh38 GIAB + core regulation     | LOW      |     98 |          98 |          0 |         98 |                   0 | 100.00%    | 3.69%                                 |
| GRCh38 exact SV + core regulation | LOW      |    821 |         821 |          0 |        821 |                   0 | 100.00%    | 0.45%                                 |
| GRCh38 ClinVar coding             | MODERATE |  59108 |       59108 |          0 |      59108 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 ClinVar cross-chromosome   | MODERATE |  43286 |       43286 |          0 |      43286 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 dbSNP                      | MODERATE |    125 |         125 |          0 |        125 |                   0 | 100.00%    | 2.91%                                 |
| GRCh38 GIAB                       | MODERATE |     31 |          31 |          0 |         31 |                   0 | 100.00%    | 11.22%                                |
| GRCh37                            | MODERATE |  23319 |       23319 |          0 |      23319 |                   0 | 100.00%    | 0.02%                                 |
| P. falciparum                     | MODERATE |   1937 |        1937 |          0 |       1937 |                   0 | 100.00%    | 0.19%                                 |
| GRCh38 GIAB + core regulation     | MODERATE |     41 |          41 |          0 |         41 |                   0 | 100.00%    | 8.60%                                 |
| GRCh38 exact SV + core regulation | MODERATE |    766 |         766 |          0 |        766 |                   0 | 100.00%    | 0.48%                                 |
| GRCh38 paired BND                 | MODIFIER |  21774 |       21774 |          0 |      21774 |                   0 | 100.00%    | 0.02%                                 |
| GRCh38 ClinVar coding             | MODIFIER | 123673 |      123673 |          0 |     123673 |                   0 | 100.00%    | 0.00%                                 |
| GRCh38 ClinVar cross-chromosome   | MODIFIER | 164672 |      164672 |          0 |     164672 |                   0 | 100.00%    | 0.00%                                 |
| GRCh38 dbSNP                      | MODIFIER |  73028 |       73028 |          0 |      73028 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 GIAB                       | MODIFIER |  54610 |       54610 |          0 |      54610 |                   0 | 100.00%    | 0.01%                                 |
| GRCh37                            | MODIFIER | 393146 |      393146 |          0 |     393146 |                   0 | 100.00%    | 0.00%                                 |
| P. falciparum                     | MODIFIER |  34276 |       34276 |          0 |      34276 |                   0 | 100.00%    | 0.01%                                 |
| GRCh38 GIAB + core regulation     | MODIFIER |  14813 |       14813 |          0 |      14813 |                   0 | 100.00%    | 0.02%                                 |
| GRCh38 exact SV + core regulation | MODIFIER |  78502 |       78502 |          0 |      78502 |                   0 | 100.00%    | 0.00%                                 |

The source artifact hash and exact Ensembl core/variation build remain
in `test/duckvep/conformance/data/conformance_history.csv` for audit and
reruns.

## Variant-induced NMD

This is a separate executable differential against the pinned VEP
Plugins release/116 `NMD.pm`. It compares `triggering`, `escaping`, and
`unresolved` for every eligible transcript pair; it does not infer NMD
from the core `NMD_transcript_variant` biotype consequence.

The `b7c7237ee686` run uses registered `variantkey_clinvar_20260706` on
chromosome 21: 49,937 source records yield 49,781 eligible alleles under
the existing 50-base limit, with no duplicate eligible alleles removed.
The oracle uses registered `vep116_grch38_cache_chr21`, containing every
chromosome-21 cache file and root metadata; the native model is
`duckvep_ensembl116_model`. Its 1,353,288 native pairs all match. The
70,521 NMD classifications include 29,954 unresolved by both engines,
not resolved predictions. These reproduce the `a84ff1500149` counts with
the checksum-verified cache. The older run with 1,331,664 native pairs
and 68,554 NMD classifications has different denominators; comparison
with that run is not an identical-workload comparison.

| revision | corpus            | model                   | exact       | mismatches | not comparable | VEP unresolved | DuckVEP unresolved | descriptive independent-pair upper 95% |
|:---------|:------------------|:------------------------|:------------|-----------:|---------------:|---------------:|-------------------:|:---------------------------------------|
| b7c7237e | nmd_clinvar_chr21 | ensembl116-grch38-final | 70521/70521 |          0 |              0 |          29954 |              29954 | 0.01%                                  |

The ledger keeps the prediction confusion matrix rather than only the
total:

| revision | corpus            | VEP_prediction | DuckVEP_prediction |     n |
|:---------|:------------------|:---------------|:-------------------|------:|
| b7c7237e | nmd_clinvar_chr21 | escaping       | escaping           |  7175 |
| b7c7237e | nmd_clinvar_chr21 | triggering     | triggering         | 33392 |
| b7c7237e | nmd_clinvar_chr21 | unresolved     | unresolved         | 29954 |

VEP projects the complete uploaded `VariationFeature` for the plugin’s
CDS and exon-position rules. DuckVEP therefore retains both geometries:
minimized edit coordinates drive consequence and sequence changes, while
the original feature endpoints drive NMD. A padded and a minimal allele
can encode the same sequence edit but cross the plugin’s inclusive
positional threshold differently.
