ROH against planted truth and read mixtures
================

Two follow-ups of <https://github.com/RGenomicsETL/duckhts/issues/349> test
the ROH decodes where the answer is known:

- **Planted truth.** Synthetic chr20 children carry autozygous segments of
  known lengths, genotype errors at a known rate and read counts simulated
  with known depth, read error, contamination, allele balance and artefact
  share. How much of the planted segments does each decode find, and how much
  does it call outside them? Does `duckhts_count_error_fit` recover the
  simulated numbers, and do the runs it decodes match the runs of the true
  numbers? Where do the bases that `stable_fraction` and `gained_bases`
  count lie?
- **Read mixtures.** Reads of a second 30× sample are mixed into a receiver’s
  chr20 reads at 1, 2, 5 and 10%. Does the fitted `contamination` follow the
  share, and does the contamination term of the read-count decode keep the
  receiver’s runs?

An earlier planted design
(<https://github.com/RGenomicsETL/duckhts/issues/329>, comment “Synthetic
planted truth”) used real children as the background: their own runs had no
truth, and the planted segments had no genotype error, so every arm found
them. Here the background is built from two unrelated founders and the
evidence carries errors and artefacts.

The scope is chr20. The report makes no timing or scaling claim; it measures
memory only to gate it.

## Inputs

From the `r/duckhtsbench` registry:

- the phased 1000 Genomes chr20 panel (`roh_ancestry_chr20_source`) and its
  pedigree (`roh_ancestry_pedigree`);
- the 108,757 chr20 evaluation sites (`roh_ancestry_chr20_af_sites`);
- chr20 CRAM slices of the 30× CRAMs of HG00403 and NA18507
  (`roh_mixture_hg00403_chr20_cram`, `roh_mixture_na18507_chr20_cram`), cut
  from the registered public CRAMs by the staging function of the BAM
  mismatch slices (`benchmarks/bam_mismatch_stage.R`), with the GRCh38
  reference (`liftover_grch38_fasta`);
- the registered read counts of the same two samples
  (`roh_counts_chr20_hg00403`, `roh_counts_chr20_na18507`), whose settings
  the mixture counts reuse.

Before anything is decoded, the render compares each input with the identity
its registry row records. A file that differs stops the render.

## Design, declared before any arm ran

The values are in `benchmarks/roh_truth_mixtures/declaration.R`, committed on
its own before the driver decoded anything. The secondary measures of the
planted truth were added in a second commit, after a reduced driver check and
before the full run; they add reported numbers and change no criterion. A
second version of the truth (`declaration_v2.R`, its own commit) followed the
first full run, which stopped at the memory gate; it is described after the
first version, whose verdicts the report keeps. This section states the design
in words.

### Planted truth

**Children.** Five populations (YRI, CEU, CHS, GIH, PEL), 20 children each.
The founders of a population are its pedigree rows without parents. For each
child, three founders are drawn without replacement: A and B give one
haplotype each, drawn at random, and C is the second genome of the read
counts. No two children share a founder. The frequencies are the ALT allele
frequencies among all founders of the child’s population, clamped to
\[0.001, 0.999\] as the staged read counts clamp theirs.

**Planted segments.** Each child gets one segment of 5, 3, 2 and 1 Mb. A start
is uniform over the positions that keep the segment between the first and the
last site; a draw is rejected when the segment overlaps a 100 kb window with
fewer than 20 panel sites, or lies within 1 Mb of another segment of the
child. Inside a segment, B’s haplotype is replaced by A’s, so the child is
homozygous there. The render checks in SQL, from the founder haplotypes and
independent of the R code that builds the children, that every true genotype
is what the plan says.

**Evidence.**

- Genotypes: at each site, with chance 1%, the reported genotype is one of
  the two other genotypes, each with chance 1/2. The children are written as
  a BGZF VCF with an index, so that `read_bcf` reports a row estimate and the
  frequency join is built on the frequency rows, as the declared memory model
  of the genotype decode assumes.
- Read counts: depth is Poisson with mean 30. With chance 0.3% a site’s alt
  count is uniform on 0 to depth (an artefactual site). Otherwise it is
  binomial with chance (1 − c) q(g) + c q(g_C), where g is the child’s and
  g_C founder C’s ALT dosage, c = 0.02 and q = 0.0005, 0.48 and 0.9995 for 0,
  1 and 2 copies. So `seq_error` is 0.0005, `contamination` 0.02,
  `allele_balance` 0.48 and `artefact_weight` 0.003.

**Arms.**

| Arm              | Decode                                                                                                                                 |
|:-----------------|:---------------------------------------------------------------------------------------------------------------------------------------|
| `gt30`, `gt100`  | `duckhts_roh_af_table` on the reported genotypes with the population frequencies, `gt_error := 30` and `100`, one query per population |
| `counts_default` | `duckhts_roh_counts` with its scalar defaults                                                                                          |
| `counts_truth`   | `duckhts_roh_counts_fitted` with the four simulated numbers                                                                            |
| `counts_fit`     | `duckhts_roh_counts_fitted` with each child’s row of `duckhts_count_error_fit` (default options), `artefact_weight_sd` included        |

A child whose fit has a NULL among the four numbers is not decoded by
`counts_fit` and counts as a child with no runs. Two more decodes, at the
fitted `artefact_weight` w − and + `artefact_weight_sd` (kept in \[0, 0.05\]),
locate the bases that `stable_fraction` and `gained_bases` count.

**Metrics.** Per child and arm: a called base is in a run of the arm; it is
true when it is in a planted segment and false otherwise. Sensitivity is true
called bases over planted bases; false discovery is false called bases over
called bases. Pooled values sum over children before the ratio. The criteria
use runs of at least 500 kb (“long”); all runs are reported too. Intervals are
95% percentile intervals of 2,000 paired bootstrap draws of the pooled
difference, resampling children within each population.

**Criteria.**

- **T1**, the artefact and fit terms help: on long runs, the interval of
  `counts_fit` − `counts_default` lies above 0 for sensitivity and below
  0.01 for false discovery.
- **T2**, the fit decodes as well as the truth: both intervals of
  `counts_fit` − `counts_truth` lie in \[−0.01, 0.01\].
- **T3**, the fit recovers the simulated numbers: in every population, the
  median fitted value lies within 0.0002 of `seq_error`, 0.003 of
  `contamination`, 0.003 of `allele_balance` and 0.001 of `artefact_weight`,
  and at least 90% of the children have status `ok`.

`gt30` − `gt100` and `counts_fit` − `gt30` are reported with intervals and no
verdict.

**Secondary measures**, declared after the reduced driver check (CEU, four
children; not evidence) and before the full run. That check showed every
arm, `gt30` included, calling long runs outside the planted segments, so the
background of two unrelated founders still holds homozygous stretches. Per
arm, without a verdict:

- the background stretches: maximal sets of consecutive panel sites at which
  the child’s true genotype is not heterozygous, from their first to their
  last site, of at least 500 kb, and their bases outside the planted
  segments;
- the false called bases inside a background stretch and elsewhere;
- the true heterozygous sites inside called runs, and the panel sites inside
  called runs outside the planted segments.

**Stability**, without a verdict: of the bases of the runs of `counts_fit`
that `stable_fraction` marks as unstable, how many are false; of the bases
that `gained_bases` counts, how many are true (bases that `counts_fit`
missed) and how many are false.

**Truth, version 2.** The first full run, with the memory gate recording
only, showed the premise of version 1 false: the long runs that `counts_fit`
called outside the planted segments held 0.0038 true heterozygous sites per
panel site, against about 0.30 for the children’s background, and 142 of
their 270 Mb lay in stretches of at least 500 kb without one true
heterozygous site. Such bases are homozygous in the constructed child. In
version 2 a true base is a base of a planted segment, or of a background
stretch: a maximal set of consecutive panel sites at which the child’s true
genotype (after planting, before the genotype errors and the read simulation)
is not heterozygous, from the position of its first to that of its last site,
of at least 500 kb (last − first + 1). The two sets are merged. Arms, run
classes, metrics, bootstrap, comparisons and the bounds of T1 and T2 are those
of version 1; T3 does not depend on the truth. Both versions’ verdicts are
reported.

### Unrelated read mixtures

The receiver is HG00403 and the donor NA18507: two staged slices of
different populations. This arm reports `contamination_relative` under its
hypothesis; the parent/child arm below tests it with a related donor.

A template (read pair) of the receiver is kept with chance 1 − α and one of
the donor with chance α, by `samtools view --subsample` with
`--subsample-seed 0`. samtools keeps a record when a hash of its read name is
below the share (`sam_view.c` lines 177-179 of samtools 1.23 in `.sync/`; the
1.24 binary bundled by RBCFTools has the same rule and, with an explicit
seed, the same seed handling), so mates stay together. The two sets are
merged into one indexed CRAM, counted and deleted before the next α. α is
0.01, 0.02, 0.05 and 0.10; α = 0 is the receiver slice alone.

The counts use the query that staged the registered counts
(`roh_counts_query_sql` in `benchmarks/roh_counts_stage.R`) with their
settings: `duckhts_somalier_bam_counts`, MAPQ ≥ 20, base quality ≥ 13,
excluded flags 1796, mate-overlap suppression `hileup_v0.1.0`, the panel’s
`INFO_AF` clamped to \[0.001, 0.999\]. The counts of each slice alone must equal
its registered counts, or the render stops.

`duckhts_count_error_fit` (default options) fits α = 0 and every mixture in
one query. Each is decoded by `duckhts_roh_counts_fitted` with its fit
(“with term”) and with `contamination` 0 and the other three numbers fitted
(“without term”). FROH is run bases over the span of the panel sites, all
runs; erosion is a decode’s FROH minus that of α = 0 with the term.

- **M1**: the fitted `contamination` increases strictly with α over 0 and
  the four shares.
- **M2**: at every α, \|erosion with term\| ≤ \|erosion without term\| + 0.002.

### Related second genome

The registered 1000 Genomes pedigree identifies HG00403 and HG00404 as
parents of HG00405 (family SH001, CHS). HG00403 receives reads from its child
HG00405. The parents have CRAMs at the existing 2504-sample S3 locator family;
HG00405’s corresponding S3 key is absent. The child is in the
[698-related-sample 30x release](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/README_2504_plus_additional_698_related_samples.txt),
with a range-readable ENA CRAM (ERR3988761) and index. The registry names that
public source and the chr20 slice derived by the same staging function as
the unrelated arm. Its independently registered counts must also match.

`benchmarks/roh_truth_mixtures/declaration_related.R` declares the criteria
before this arm runs. The alpha grid, read-name hash seed (mates together),
counting settings, frequencies and decoder are the same as the unrelated arm.
Nominal alpha is a sampling fraction, not necessarily the realised donor
share at unequal source depths; both are reported.

- **M3**: `contamination_relative` is within max(0.005, 0.25 × α) of α at
  every nonzero α. The relative log-likelihood exceeds the unrelated one at
  α = 0.05 and 0.10. At 0.01/0.02 the likelihood sign is descriptive, not a
  gate: the earlier count-level result separated at 5%, not at 2%.
- **M4**: unrelated `contamination` is within max(0.003, 0.25 × α/2) of α/2
  at every nonzero α, as expected for a parent/child donor under that model.
- **R1**: at every α, \|erosion with the relative-fit term\| ≤ \|erosion without
  it\| + 0.002. Erosion is relative to α = 0 with the relative-fit term.

The with-term decode substitutes `contamination_relative` for `contamination`
in the fit table. Other nuisance estimates remain those exposed by the public
fit, from its unrelated hypothesis. This is a plug-in population-frequency
emission, not a relationship-aware ROH emission. The tolerance floors and
25% relative errors allow read sampling and nuisance fitting; they are not
confidence intervals. No unrelated criterion is changed.

### Gates

The render stops when an input differs from its registered identity, when a
child differs from its plan, when a slice’s counts differ from its registered
counts, when a mixture would leave less than 8 GiB free, when the bases that
`stable_fraction` or `gained_bases` count differ from those the two extra
decodes locate, or when a DuckDB query spills or passes its ceiling. Each
query runs in a fresh DuckDB CLI process under `/usr/bin/time`, with the JSON
profiler and `max_temp_directory_size` zero. The ceiling is the overhead (the
median peak RSS of three empty runs at the query’s threads) plus 3 × its live
state:

- count decode (`benchmark_roh_scaling.md`): rows × 24 + threads × sites ×
  36.25 bytes;
- genotype decode (`duckhts_roh_af_table`, the same report): samples ×
  sites × 16 + sites × 128 + threads × sites × 36.25 bytes;
- fit (`benchmark_count_error_fit.md`): cells × 312 + distinct (depth, alt
  count) pairs × 88 + samples × 1,001 × 24 bytes, with the cells and pairs
  counted untimed;
- read counting: workers × (`max_region_bytes` + `reference_cache_bytes` +
  `remote_cache_bytes` + `remote_block_bytes` + the chr20 reference length)
  - sites × 512 bytes, from the per-job bounds of the catalog entry of
    `duckhts_somalier_bam_counts` at its defaults, the contig being the CRAM
    decoder’s copy of the reference.

Decodes and fits run at two threads; counting at four threads and four
workers.

`ROH_TRUTH_SIZE=reduced` checks the driver: CEU only with four children, and
the mixture at α = 0.05 only, written to a temporary directory instead of
`benchmarks/results/roh-truth-mixtures/`. Every gate still runs. A reduced run
is not evidence. The default is `full`. The bundled samtools needs RBCFTools
1.24-1.1.0 or later on the library path.

``` sh
R_LIBS=/tmp/RBCFTools-release-lib taskset -c 12-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_roh_truth_and_mixtures.Rmd")'
# Driver check only; writes the page outside benchmarks/:
ROH_TRUTH_SIZE=reduced taskset -c 0-3 Rscript -e \
  'rmarkdown::render("benchmarks/benchmark_roh_truth_and_mixtures.Rmd", output_dir = tempdir())'
```

## Run

| Item                     | Value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| revision                 | 4d09d23ef2e59b10f394e976d382ed8055674817                         |
| src_tree                 | e724eb4f772457c455d8dd4610c4d7555c0e8876                         |
| extension_sha256         | fd5bff50d6d33fadbc692e472e2e9059c119014fbb29730d757894a20da2a5f5 |
| duckdb_cli               | v1.5.1 (Variegata) 7dbb2e646f                                    |
| r_version                | R version 4.6.0 (2026-04-24)                                     |
| host_load                | 4.59, 6.31, 5.63                                                 |
| cpu_affinity             | 12-15                                                            |
| threads                  | 2                                                                |
| size                     | full                                                             |
| populations              | YRI,CEU,CHS,GIH,PEL                                              |
| children                 | 100                                                              |
| sites                    | 108757                                                           |
| span_bases               | 64251060                                                         |
| counts_rows              | 10875700                                                         |
| fit_cells                | 4228633                                                          |
| fit_pairs                | 90885                                                            |
| overhead_mib             | 40.2                                                             |
| counts_sha256            | ac858137459611400981e25fcf04e71efb49301bccfd3f40fa816f2ab73004e8 |
| vcf_sha256               | 0c9a26a117056d3071c1dbde9581eb9fd1f302891e89db77cbcfb136b3226cae |
| children_with_usable_fit | 100                                                              |
| receiver                 | HG00403                                                          |
| donor                    | NA18507                                                          |
| alphas                   | 0.01,0.02,0.05,0.1                                               |
| size                     | full                                                             |
| sites                    | 108757                                                           |
| span_bases               | 64251060                                                         |
| samtools                 | 1.24                                                             |
| rbcftools                | 1.24.1.1.0                                                       |
| fit_cells                | 252355                                                           |
| fit_pairs                | 6079                                                             |
| overhead_decode_mib      | 40.2                                                             |
| overhead_count_mib       | 40.2                                                             |

## Planted truth

### The children

| Children | Planted segments | Planted bases per child | Planted sites | Heterozygous planted sites | Genotype error share | Artefact site share | Mean depth |
|---------:|-----------------:|------------------------:|--------------:|---------------------------:|---------------------:|--------------------:|-----------:|
|      100 |              400 |              11,000,000 |     2,075,744 |                          0 |               0.0100 |              0.0030 |      30.00 |

### Fit against the simulated numbers (T3)

The median over each population’s children of each fitted number.

| Population | Number          | Simulated | Median fit | Tolerance | Verdict |
|:-----------|:----------------|----------:|-----------:|----------:|:--------|
| CEU        | seq_error       |    0.0005 | 0.00055319 |    0.0002 | pass    |
| CEU        | contamination   |      0.02 |   0.020086 |     0.003 | pass    |
| CEU        | allele_balance  |      0.48 |    0.48008 |     0.003 | pass    |
| CEU        | artefact_weight |     0.003 |  0.0027723 |     0.001 | pass    |
| CHS        | seq_error       |    0.0005 | 0.00048802 |    0.0002 | pass    |
| CHS        | contamination   |      0.02 |    0.01987 |     0.003 | pass    |
| CHS        | allele_balance  |      0.48 |    0.47985 |     0.003 | pass    |
| CHS        | artefact_weight |     0.003 |  0.0029022 |     0.001 | pass    |
| GIH        | seq_error       |    0.0005 | 0.00050704 |    0.0002 | pass    |
| GIH        | contamination   |      0.02 |   0.020395 |     0.003 | pass    |
| GIH        | allele_balance  |      0.48 |     0.4801 |     0.003 | pass    |
| GIH        | artefact_weight |     0.003 |  0.0032384 |     0.001 | pass    |
| PEL        | seq_error       |    0.0005 | 0.00040518 |    0.0002 | pass    |
| PEL        | contamination   |      0.02 |   0.020436 |     0.003 | pass    |
| PEL        | allele_balance  |      0.48 |    0.47997 |     0.003 | pass    |
| PEL        | artefact_weight |     0.003 |  0.0031118 |     0.001 | pass    |
| YRI        | seq_error       |    0.0005 | 0.00050249 |    0.0002 | pass    |
| YRI        | contamination   |      0.02 |   0.020122 |     0.003 | pass    |
| YRI        | allele_balance  |      0.48 |    0.47983 |     0.003 | pass    |
| YRI        | artefact_weight |     0.003 |   0.002777 |     0.001 | pass    |

Status: ok 100 of 100 children (`ok` share 1.00, at least 0.9 required).

### Sensitivity and false discovery

Pooled over all children, against both truths: version 1 (planted segments)
and version 2 (planted segments and background stretches).

| Runs     | Arm            |  Called bases | Sensitivity v1 | False discovery v1 | Sensitivity v2 | False discovery v2 |
|:---------|:---------------|--------------:|---------------:|-------------------:|---------------:|-------------------:|
| ≥ 500 kb | gt30           | 1,367,240,554 |         0.9991 |             0.1962 |         0.9025 |             0.0931 |
| ≥ 500 kb | gt100          |   563,024,502 |         0.4653 |             0.0909 |         0.4088 |             0.0023 |
| ≥ 500 kb | counts_default |   585,774,552 |         0.4781 |             0.1022 |         0.4259 |             0.0010 |
| ≥ 500 kb | counts_truth   | 1,367,913,243 |         0.9993 |             0.1964 |         0.9010 |             0.0951 |
| ≥ 500 kb | counts_fit     | 1,369,713,321 |         0.9993 |             0.1974 |         0.9038 |             0.0934 |
| all      | gt30           | 2,389,992,995 |         0.9991 |             0.5402 |         0.9057 |             0.4793 |
| all      | gt100          | 1,990,346,025 |         0.9478 |             0.4762 |         0.8424 |             0.4185 |
| all      | counts_default | 2,084,548,593 |         0.9788 |             0.4835 |         0.8692 |             0.4271 |
| all      | counts_truth   | 2,410,255,565 |         0.9993 |             0.5439 |         0.9042 |             0.4846 |
| all      | counts_fit     | 2,411,594,067 |         0.9993 |             0.5442 |         0.9070 |             0.4833 |

### Comparisons

Pooled difference, first arm minus second, with its 95% bootstrap interval,
against each truth.

| Criterion | Runs     | Comparison                  | Metric          |                    Version 1 |                    Version 2 |
|:----------|:---------|:----------------------------|:----------------|-----------------------------:|-----------------------------:|
| T1        | ≥ 500 kb | counts_fit − counts_default | false_discovery | +0.0952 \[+0.0747, +0.1166\] | +0.0924 \[+0.0812, +0.1041\] |
| T1        | ≥ 500 kb | counts_fit − counts_default | sensitivity     | +0.5213 \[+0.4989, +0.5439\] | +0.4779 \[+0.4530, +0.5016\] |
| T2        | ≥ 500 kb | counts_fit − counts_truth   | false_discovery | +0.0011 \[-0.0028, +0.0069\] | -0.0017 \[-0.0035, +0.0000\] |
| T2        | ≥ 500 kb | counts_fit − counts_truth   | sensitivity     | -0.0000 \[-0.0000, +0.0000\] | +0.0028 \[-0.0000, +0.0086\] |
|           | ≥ 500 kb | gt30 − gt100                | false_discovery | +0.1052 \[+0.0818, +0.1294\] | +0.0907 \[+0.0772, +0.1054\] |
|           | ≥ 500 kb | gt30 − gt100                | sensitivity     | +0.5338 \[+0.5156, +0.5532\] | +0.4937 \[+0.4748, +0.5138\] |
|           | ≥ 500 kb | counts_fit − gt30           | false_discovery | +0.0013 \[-0.0065, +0.0076\] | +0.0003 \[-0.0086, +0.0075\] |
|           | ≥ 500 kb | counts_fit − gt30           | sensitivity     | +0.0002 \[+0.0001, +0.0004\] | +0.0013 \[-0.0001, +0.0031\] |
| T1        | all      | counts_fit − counts_default | false_discovery | +0.0607 \[+0.0568, +0.0650\] | +0.0562 \[+0.0527, +0.0597\] |
| T1        | all      | counts_fit − counts_default | sensitivity     | +0.0205 \[+0.0189, +0.0222\] | +0.0378 \[+0.0291, +0.0480\] |
| T2        | all      | counts_fit − counts_truth   | false_discovery | +0.0003 \[-0.0011, +0.0021\] | -0.0013 \[-0.0030, +0.0000\] |
| T2        | all      | counts_fit − counts_truth   | sensitivity     | -0.0000 \[-0.0000, +0.0000\] | +0.0028 \[-0.0001, +0.0085\] |
|           | all      | gt30 − gt100                | false_discovery | +0.0640 \[+0.0602, +0.0679\] | +0.0609 \[+0.0569, +0.0651\] |
|           | all      | gt30 − gt100                | sensitivity     | +0.0513 \[+0.0482, +0.0547\] | +0.0633 \[+0.0556, +0.0719\] |
|           | all      | counts_fit − gt30           | false_discovery | +0.0040 \[+0.0017, +0.0060\] | +0.0039 \[+0.0013, +0.0063\] |
|           | all      | counts_fit − gt30           | sensitivity     | +0.0002 \[+0.0001, +0.0004\] | +0.0013 \[+0.0002, +0.0029\] |

| Criterion | Version 1 | Version 2 |
|:----------|:----------|:----------|
| T1        | fail      | fail      |
| T2        | pass      | pass      |
| T3        | pass      | pass      |

T3 compares the fit with the simulated numbers and does not depend on the
truth.

### By population

Runs of at least 500 kb.

| Population | Arm            | Sensitivity v1 | False discovery v1 | Sensitivity v2 | False discovery v2 |
|:-----------|:---------------|---------------:|-------------------:|---------------:|-------------------:|
| YRI        | gt30           |         0.9989 |             0.0713 |         0.8051 |             0.0435 |
| YRI        | gt100          |         0.4408 |             0.0180 |         0.3513 |             0.0000 |
| YRI        | counts_default |         0.4918 |             0.0196 |         0.3926 |             0.0000 |
| YRI        | counts_truth   |         0.9992 |             0.0717 |         0.8053 |             0.0439 |
| YRI        | counts_fit     |         0.9992 |             0.0711 |         0.8053 |             0.0434 |
| CEU        | gt30           |         0.9992 |             0.2083 |         0.9132 |             0.1034 |
| CEU        | gt100          |         0.4571 |             0.0796 |         0.3998 |             0.0025 |
| CEU        | counts_default |         0.4848 |             0.0925 |         0.4311 |             0.0000 |
| CEU        | counts_truth   |         0.9993 |             0.2113 |         0.9133 |             0.1067 |
| CEU        | counts_fit     |         0.9993 |             0.2093 |         0.9134 |             0.1044 |
| CHS        | gt30           |         0.9989 |             0.2168 |         0.9321 |             0.0955 |
| CHS        | gt100          |         0.4726 |             0.1507 |         0.4468 |             0.0061 |
| CHS        | counts_default |         0.4688 |             0.1489 |         0.4450 |             0.0000 |
| CHS        | counts_truth   |         0.9994 |             0.2259 |         0.9361 |             0.1027 |
| CHS        | counts_fit     |         0.9994 |             0.2260 |         0.9361 |             0.1027 |
| GIH        | gt30           |         0.9992 |             0.1931 |         0.9333 |             0.1211 |
| GIH        | gt100          |         0.4854 |             0.0651 |         0.4449 |             0.0009 |
| GIH        | counts_default |         0.4822 |             0.0554 |         0.4378 |             0.0000 |
| GIH        | counts_truth   |         0.9994 |             0.1823 |         0.9339 |             0.1090 |
| GIH        | counts_fit     |         0.9994 |             0.1799 |         0.9339 |             0.1065 |
| PEL        | gt30           |         0.9993 |             0.2670 |         0.9319 |             0.0948 |
| PEL        | gt100          |         0.4707 |             0.1253 |         0.4057 |             0.0016 |
| PEL        | counts_default |         0.4629 |             0.1815 |         0.4251 |             0.0046 |
| PEL        | counts_truth   |         0.9994 |             0.2658 |         0.9200 |             0.1050 |
| PEL        | counts_fit     |         0.9994 |             0.2740 |         0.9334 |             0.1021 |

### By planted length

The share of each planted segment that the runs of at least 500 kb cover,
averaged over the segments of one length; “complete” counts segments covered
to at least 99%, “missed” those not covered at all.

| Arm            | Length (Mb) | Segments | Mean covered | Complete | Missed |
|:---------------|------------:|---------:|-------------:|---------:|-------:|
| gt30           |           1 |      100 |        0.998 |     0.99 |   0.00 |
| gt30           |           2 |      100 |        0.998 |     0.97 |   0.00 |
| gt30           |           3 |      100 |        0.999 |     1.00 |   0.00 |
| gt30           |           5 |      100 |        0.999 |     1.00 |   0.00 |
| gt100          |           1 |      100 |        0.387 |     0.07 |   0.44 |
| gt100          |           2 |      100 |        0.426 |     0.02 |   0.21 |
| gt100          |           3 |      100 |        0.485 |     0.01 |   0.03 |
| gt100          |           5 |      100 |        0.485 |     0.00 |   0.00 |
| counts_default |           1 |      100 |        0.399 |     0.06 |   0.42 |
| counts_default |           2 |      100 |        0.441 |     0.01 |   0.18 |
| counts_default |           3 |      100 |        0.501 |     0.00 |   0.03 |
| counts_default |           5 |      100 |        0.495 |     0.00 |   0.00 |
| counts_truth   |           1 |      100 |        0.998 |     0.99 |   0.00 |
| counts_truth   |           2 |      100 |        0.999 |     0.99 |   0.00 |
| counts_truth   |           3 |      100 |        0.999 |     0.99 |   0.00 |
| counts_truth   |           5 |      100 |        1.000 |     1.00 |   0.00 |
| counts_fit     |           1 |      100 |        0.998 |     0.99 |   0.00 |
| counts_fit     |           2 |      100 |        0.999 |     0.99 |   0.00 |
| counts_fit     |           3 |      100 |        0.999 |     0.99 |   0.00 |
| counts_fit     |           5 |      100 |        1.000 |     1.00 |   0.00 |

### The background (secondary)

Against version 1 of the truth. Background stretches have no heterozygous
true genotype over at least 500 kb. Bases inside one are homozygous in the
child; version 1 counts them false because they are outside the planted
segments, and version 2 counts them true. The heterozygosity columns compare
the true heterozygous sites per panel site inside an arm’s calls outside the
planted segments with the child’s background as a whole.

| Runs     | Arm            |   False bases | In background stretches |     Elsewhere | Heterozygous sites in calls | Background sites in calls | Heterozygosity in calls |
|:---------|:---------------|--------------:|------------------------:|--------------:|----------------------------:|--------------------------:|------------------------:|
| ≥ 500 kb | gt30           |   268,227,680 |             140,989,178 |   127,238,502 |                       1,027 |                   294,649 |                  0.0035 |
| ≥ 500 kb | gt100          |    51,204,163 |              49,881,986 |     1,322,177 |                          14 |                    40,267 |                  0.0003 |
| ≥ 500 kb | counts_default |    59,872,664 |              59,294,470 |       578,194 |                           1 |                    53,423 |                  0.0000 |
| ≥ 500 kb | counts_truth   |   268,633,915 |             138,599,693 |   130,034,222 |                       1,169 |                   307,002 |                  0.0038 |
| ≥ 500 kb | counts_fit     |   270,435,601 |             142,497,400 |   127,938,201 |                       1,164 |                   304,814 |                  0.0038 |
| all      | gt30           | 1,290,980,121 |             145,349,028 | 1,145,631,093 |                      11,306 |                 2,279,359 |                  0.0050 |
| all      | gt100          |   947,715,734 |             114,789,282 |   832,926,452 |                         887 |                 1,744,658 |                  0.0005 |
| all      | counts_default | 1,007,826,062 |             117,528,173 |   890,297,889 |                          24 |                 1,854,819 |                  0.0000 |
| all      | counts_truth   | 1,310,976,237 |             143,022,922 | 1,167,953,315 |                      12,010 |                 2,330,494 |                  0.0052 |
| all      | counts_fit     | 1,312,316,347 |             146,851,205 | 1,165,465,142 |                      11,869 |                 2,326,778 |                  0.0051 |

The background stretches hold 273,892,716 bases outside the
planted segments over 100 children
(2,738,927 per child). The
heterozygosity of the background as a whole is
0.2996
heterozygous sites per panel site.

### Stability and gained bases against truth

Runs of `counts_fit`, all lengths. Unstable bases are the bases of a run that
`stable_fraction` does not count (not called at both w − sd and w + sd);
gained bases are those `gained_bases` counts. A false base is outside the
truth; a true gained base is one that `counts_fit` missed.

| Truth |     Run bases | False share of run bases | Unstable bases | False share of unstable | False share of stable | Missed true bases | Gained bases | Gained, true | Gained, false | Bases of new runs at w + sd |
|:------|--------------:|-------------------------:|---------------:|------------------------:|----------------------:|------------------:|-------------:|-------------:|--------------:|----------------------------:|
| v1    | 2,411,594,067 |                   0.5442 |     55,445,968 |                  0.9815 |                0.5339 |           722,280 |   11,035,870 |        8,742 |    11,027,128 |                   5,670,684 |
| v2    | 2,411,594,067 |                   0.4833 |     55,445,968 |                  0.8548 |                0.4745 |       127,763,791 |   11,035,870 |      267,062 |    10,768,808 |                   5,670,684 |

## Read mixtures

### Identities

| Slice   | Rows    | Counts SHA-256                                                   | Equals registered |
|:--------|:--------|:-----------------------------------------------------------------|:------------------|
| HG00403 | 108,757 | c0698ddc47beccecb3bf68177450ad6507e5edab9d318df4bc3680b7b6adc6da | yes               |
| NA18507 | 108,757 | 5c9e70d6a891f4a8bce5d5158f8b4706cb43b97a35fb37e4ab996394a563ba0e | yes               |

### Mixtures

|    α | Receiver alignments | Donor alignments | Donor share of alignments |
|-----:|--------------------:|-----------------:|--------------------------:|
| 0.01 |          16,647,490 |          176,122 |                    0.0105 |
| 0.02 |          16,479,939 |          349,659 |                    0.0208 |
| 0.05 |          15,975,052 |          874,405 |                    0.0519 |
| 0.10 |          15,135,131 |        1,746,446 |                    0.1035 |

### Fit

|    α | contamination | contamination_sd | contamination / α | contamination_relative | relative − unrelated log-likelihood | seq_error | allele_balance | artefact_weight | status |
|-----:|--------------:|-----------------:|------------------:|-----------------------:|------------------------------------:|----------:|---------------:|----------------:|:-------|
| 0.00 |     0.0003186 |        6.584e-05 |                   |              0.0006385 |                                 0.1 |  0.000129 |          0.498 |        0.002955 | ok     |
| 0.01 |       0.01204 |        0.0003317 |              1.20 |                0.02238 |                               -42.5 |  0.000135 |          0.498 |        0.003001 | ok     |
| 0.02 |       0.02379 |        0.0005304 |              1.19 |                 0.0416 |                              -164.4 | 0.0001513 |          0.498 |        0.002938 | ok     |
| 0.05 |       0.05801 |         0.001076 |              1.16 |                0.09246 |                              -874.5 | 0.0002852 |          0.498 |        0.003461 | ok     |
| 0.10 |        0.1146 |         0.001382 |              1.15 |                 0.1631 |                             -2808.7 | 0.0004022 |          0.498 |        0.003015 | ok     |

### FROH erosion

| Decode       |    α | Runs |   FROH | Erosion | Jaccard |
|:-------------|-----:|-----:|-------:|--------:|--------:|
| with_term    | 0.00 |  189 | 0.3446 | +0.0000 |   1.000 |
| with_term    | 0.01 |  189 | 0.3448 | +0.0003 |   0.999 |
| with_term    | 0.02 |  190 | 0.3439 | -0.0007 |   0.993 |
| with_term    | 0.05 |  189 | 0.3464 | +0.0019 |   0.992 |
| with_term    | 0.10 |  185 | 0.3421 | -0.0024 |   0.980 |
| without_term | 0.00 |  189 | 0.3446 | +0.0000 |   1.000 |
| without_term | 0.01 |  185 | 0.3424 | -0.0021 |   0.994 |
| without_term | 0.02 |  182 | 0.3370 | -0.0076 |   0.978 |
| without_term | 0.05 |  148 | 0.3087 | -0.0358 |   0.896 |
| without_term | 0.10 |   87 | 0.2134 | -0.1311 |   0.619 |

|    α | Erosion with term | Erosion without term | M2   |
|-----:|------------------:|---------------------:|:-----|
| 0.01 |           +0.0003 |              -0.0021 | pass |
| 0.02 |           -0.0007 |              -0.0076 | pass |
| 0.05 |           +0.0019 |              -0.0358 | pass |
| 0.10 |           -0.0024 |              -0.1311 | pass |

| Criterion | Verdict |
|:----------|:--------|
| M1        | pass    |
| M2        | pass    |

## Related-donor read mixtures

Receiver HG00403, child donor HG00405. Alpha 0 is the receiver-only control.

| Slice   | Rows    | Counts SHA-256                                                   | Equals registered |
|:--------|:--------|:-----------------------------------------------------------------|:------------------|
| HG00403 | 108,757 | c0698ddc47beccecb3bf68177450ad6507e5edab9d318df4bc3680b7b6adc6da | yes               |
| HG00405 | 108,757 | 993e802bf4a03df1c838035f7874edb089ff3458fa42347caf709722ba4abfaa | yes               |

|    α | Receiver alignments | Donor alignments | Donor share of alignments |
|-----:|:--------------------|:-----------------|:--------------------------|
| 0.01 | 16,647,490          | 169,581          | 0.0101                    |
| 0.02 | 16,479,939          | 338,044          | 0.0201                    |
| 0.05 | 15,975,052          | 846,591          | 0.0503                    |
| 0.10 | 15,135,131          | 1,696,197        | 0.1008                    |

### Relative and unrelated fits

|    α | contamination_relative | Relative error | M3 tolerance | M3 recovery | contamination | α/2   | M4 tolerance | M4 recovery | status |
|-----:|:-----------------------|:---------------|:-------------|:------------|:--------------|:------|:-------------|:------------|:-------|
| 0.00 | 0.0006385              | 0.0006385      | 0.005        | control     | 0.0003186     | 0     | 0.003        | control     | ok     |
| 0.01 | 0.00812                | -0.00188       | 0.005        | pass        | 0.004083      | 0.005 | 0.003        | pass        | ok     |
| 0.02 | 0.01558                | -0.004416      | 0.005        | pass        | 0.007824      | 0.01  | 0.003        | pass        | ok     |
| 0.05 | 0.0403                 | -0.009696      | 0.0125       | pass        | 0.02026       | 0.025 | 0.00625      | pass        | ok     |
| 0.10 | 0.08458                | -0.01542       | 0.025        | pass        | 0.04297       | 0.05  | 0.0125       | pass        | ok     |

|    α | Relative log-likelihood | Unrelated log-likelihood | Relative − unrelated | Likelihood gate |
|-----:|:------------------------|:-------------------------|:---------------------|:----------------|
| 0.00 | -170943.3               | -170943.4                | +0.076               | descriptive     |
| 0.01 | -180685.9               | -180686.6                | +0.701               | descriptive     |
| 0.02 | -187397.3               | -187405.6                | +8.321               | descriptive     |
| 0.05 | -203016.3               | -203089.3                | +72.964              | pass            |
| 0.10 | -219205.4               | -219605.2                | +399.819             | pass            |

### FROH with the relative-fit term and without

| Decode       |    α | Runs | FROH   | Erosion | Jaccard |
|:-------------|-----:|-----:|:-------|:--------|:--------|
| with_term    | 0.00 |  189 | 0.3446 | +0.0000 | 1.000   |
| with_term    | 0.01 |  189 | 0.3449 | +0.0003 | 0.999   |
| with_term    | 0.02 |  189 | 0.3423 | -0.0023 | 0.991   |
| with_term    | 0.05 |  189 | 0.3287 | -0.0158 | 0.953   |
| with_term    | 0.10 |  189 | 0.3237 | -0.0208 | 0.938   |
| without_term | 0.00 |  189 | 0.3446 | +0.0000 | 1.000   |
| without_term | 0.01 |  188 | 0.3438 | -0.0008 | 0.997   |
| without_term | 0.02 |  186 | 0.3398 | -0.0048 | 0.986   |
| without_term | 0.05 |  178 | 0.3225 | -0.0221 | 0.936   |
| without_term | 0.10 |  154 | 0.2587 | -0.0858 | 0.751   |

|    α | Erosion with relative term | Erosion without term | R1   |
|-----:|:---------------------------|:---------------------|:-----|
| 0.01 | +0.0003                    | -0.0008              | pass |
| 0.02 | -0.0023                    | -0.0048              | pass |
| 0.05 | -0.0158                    | -0.0221              | pass |
| 0.10 | -0.0208                    | -0.0858              | pass |

| Criterion | Verdict |
|:----------|:--------|
| M3        | pass    |
| M4        | pass    |
| R1        | pass    |

## Gates

Every DuckDB query, with its latency, peak RSS, DuckDB’s peak buffer memory,
its ceiling and its temporary bytes. The render stops on any spill or any
peak RSS above its ceiling, and it reached this table.

| Part     | Query             | Kind            | Threads | Seconds | Peak RSS MiB | Peak buffer MiB | Ceiling MiB | Temporary bytes |
|:---------|:------------------|:----------------|--------:|--------:|-------------:|----------------:|------------:|----------------:|
| truth    | empty_1_threads_2 | empty           |       2 |   0.001 |         40.2 |             0.2 |             |               0 |
| truth    | empty_2_threads_2 | empty           |       2 |   0.005 |         40.3 |             0.3 |             |               0 |
| truth    | empty_3_threads_2 | empty           |       2 |   0.001 |         40.1 |             0.2 |             |               0 |
| truth    | gt30_YRI          | genotype_decode |       2 |   1.316 |        119.0 |           124.7 |       202.2 |               0 |
| truth    | gt30_CEU          | genotype_decode |       2 |   1.357 |        119.6 |            92.7 |       202.2 |               0 |
| truth    | gt30_CHS          | genotype_decode |       2 |   1.401 |        120.5 |           124.7 |       202.2 |               0 |
| truth    | gt30_GIH          | genotype_decode |       2 |   1.347 |        111.7 |           124.7 |       202.2 |               0 |
| truth    | gt30_PEL          | genotype_decode |       2 |   1.351 |        117.1 |           124.7 |       202.2 |               0 |
| truth    | gt100_YRI         | genotype_decode |       2 |   1.360 |        121.4 |           124.6 |       202.2 |               0 |
| truth    | gt100_CEU         | genotype_decode |       2 |   1.415 |        128.9 |            92.7 |       202.2 |               0 |
| truth    | gt100_CHS         | genotype_decode |       2 |   1.403 |        118.1 |           124.6 |       202.2 |               0 |
| truth    | gt100_GIH         | genotype_decode |       2 |   1.392 |        116.2 |           124.7 |       202.2 |               0 |
| truth    | gt100_PEL         | genotype_decode |       2 |   1.392 |        122.8 |           124.7 |       202.2 |               0 |
| truth    | counts_default    | count_decode    |       2 |   1.633 |        464.1 |           572.4 |       809.6 |               0 |
| truth    | fit               | fit             |       2 | 176.537 |        632.2 |           986.5 |      3844.6 |               0 |
| truth    | counts_truth      | count_decode    |       2 |   5.771 |        470.4 |           607.1 |       809.6 |               0 |
| truth    | counts_fit        | count_decode    |       2 |   5.126 |        480.6 |           575.4 |       809.6 |               0 |
| truth    | fit_low           | count_decode    |       2 |   5.156 |        472.5 |           607.1 |       809.6 |               0 |
| truth    | fit_high          | count_decode    |       2 |   5.101 |        460.2 |           575.1 |       809.6 |               0 |
| mixtures | empty_1_threads_2 | empty           |       2 |   0.001 |         40.2 |             0.3 |             |               0 |
| mixtures | empty_2_threads_2 | empty           |       2 |   0.025 |         40.3 |             0.3 |             |               0 |
| mixtures | empty_3_threads_2 | empty           |       2 |   0.002 |         40.2 |             0.3 |             |               0 |
| mixtures | empty_1_threads_4 | empty           |       4 |   0.002 |         40.1 |             0.3 |             |               0 |
| mixtures | empty_2_threads_4 | empty           |       4 |   0.001 |         40.2 |             0.2 |             |               0 |
| mixtures | empty_3_threads_4 | empty           |       4 |   0.001 |         40.4 |             0.3 |             |               0 |
| mixtures | count_HG00403     | read_counting   |       4 |   5.390 |        162.9 |            43.8 |      3253.0 |               0 |
| mixtures | count_NA18507     | read_counting   |       4 |   5.689 |        169.6 |            43.8 |      3253.0 |               0 |
| mixtures | count_alpha_0.01  | read_counting   |       4 |   5.106 |        164.7 |            43.8 |      3253.0 |               0 |
| mixtures | count_alpha_0.02  | read_counting   |       4 |   5.110 |        168.0 |            43.8 |      3253.0 |               0 |
| mixtures | count_alpha_0.05  | read_counting   |       4 |   5.159 |        171.0 |            43.8 |      3253.0 |               0 |
| mixtures | count_alpha_0.10  | read_counting   |       4 |   5.216 |        170.1 |            43.8 |      3253.0 |               0 |
| mixtures | mixture_fit       | fit             |       2 |   9.625 |         87.3 |            69.3 |       267.3 |               0 |
| mixtures | with_term         | count_decode    |       2 |   0.304 |         72.6 |            55.4 |       100.1 |               0 |
| mixtures | without_term      | count_decode    |       2 |   0.379 |         79.9 |            39.3 |       100.1 |               0 |
| related  | empty_1_threads_2 | empty           |       2 |   0.001 |         40.3 |             0.3 |             |               0 |
| related  | empty_2_threads_2 | empty           |       2 |   0.004 |         40.3 |             0.3 |             |               0 |
| related  | empty_3_threads_2 | empty           |       2 |   0.002 |         40.1 |             0.3 |             |               0 |
| related  | empty_1_threads_4 | empty           |       4 |   0.003 |         40.5 |             0.3 |             |               0 |
| related  | empty_2_threads_4 | empty           |       4 |   0.002 |         40.2 |             0.3 |             |               0 |
| related  | empty_3_threads_4 | empty           |       4 |   0.001 |         40.4 |             0.2 |             |               0 |
| related  | count_HG00403     | read_counting   |       4 |   5.227 |        168.8 |            43.8 |      3253.2 |               0 |
| related  | count_HG00405     | read_counting   |       4 |   5.145 |        162.2 |            43.8 |      3253.2 |               0 |
| related  | count_alpha_0.01  | read_counting   |       4 |   5.098 |        168.5 |            43.8 |      3253.2 |               0 |
| related  | count_alpha_0.02  | read_counting   |       4 |   5.115 |        165.5 |            43.8 |      3253.2 |               0 |
| related  | count_alpha_0.05  | read_counting   |       4 |   5.054 |        159.8 |            43.8 |      3253.2 |               0 |
| related  | count_alpha_0.10  | read_counting   |       4 |   5.142 |        169.6 |            43.8 |      3253.2 |               0 |
| related  | mixture_fit       | fit             |       2 |   7.751 |         85.9 |            64.9 |       249.2 |               0 |
| related  | with_term         | count_decode    |       2 |   0.296 |         73.2 |            55.3 |       100.2 |               0 |
| related  | without_term      | count_decode    |       2 |   0.295 |         77.1 |            55.3 |       100.2 |               0 |

Every observation is in `benchmarks/results/roh-truth-mixtures/`: the children, the planted segments, the fits, every run (gzip-compressed), the per-child scores and the queries.

## Limitations

- One chromosome. The children are synthetic: real haplotypes, simulated
  errors and reads. The read model of the simulation is the model of the fit
  and of the decode, apart from the depth distribution and the second genome
  being a real founder; real reads have errors that depend on base quality,
  mapping and sequence context.
- The population frequencies include the founders of the children, two of
  about 200 haplotypes per site.
- Version 1 counts calls outside the planted segments as false, but the
  background of two unrelated founders is not free of homozygous stretches.
  Version 2 counts stretches of at least 500 kb without a true heterozygous
  site as true; a call through a longer stretch broken by one or two
  heterozygous sites, or through a shorter stretch, still counts as false
  there. The founders are the pedigree’s founders; cryptic relatedness
  between them is not checked.
- One error rate, one depth and one set of read parameters. Other rates, a
  low depth and a contaminant of another population are not tested here.
- Two mixture pairs: one unrelated pair of two populations and one
  parent/child pair. The receiver’s frequencies are the panel’s pooled
  frequencies, not its population’s, so a donor of another population biases
  the fitted contamination upward, as the count-level mixture of
  `benchmark_count_error_fit.md` shows. The relative fit’s term is tested as a
  plug-in to the existing emission; a relationship-aware decode is not tested.
