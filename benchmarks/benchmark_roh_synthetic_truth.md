ROH arms against a planted truth
================

## Question

`benchmark_roh_ancestry.md` asks whether ancestry-tuned allele
frequencies (`duckhts_roh_ancestry`) call fewer unsupported runs of
homozygosity (ROH) than pooled or single-population frequencies
(`duckhts_roh_af_table`). Its truth is a heuristic with 23 intervals in
18 children, so its criterion compares means that differ by one to three
calls (<https://github.com/RGenomicsETL/duckhts/issues/329>).

This report asks the same question against a truth that is known
exactly. Autozygous segments are planted into the 377 trio children of
chr20, and each arm is scored on the planted bases.

## Evaluation protocol

This section was written and committed before any arm ran on the
synthetic children. Nothing in it was changed after the arm results were
known.

### Planting

The source is registry artifact `roh_ancestry_chr20_children_bcf`: 377
children, 849,143 biallelic SNV records, phased genotypes. Every child
receives three planted segments, one each of 5 Mb, 2 Mb and 1 Mb, drawn
in that order. That is 1,131 segments and 8 Mb per child.

- A start is uniform over the positions that keep the segment between
  the first and the last record of the BCF. The generator is R’s
  Mersenne-Twister with rejection sampling, seed `20261004`. The
  children are drawn in BCF sample order.
- A draw is rejected when the segment overlaps a 100 kb window that
  holds fewer than 100 records, or when it lies within 1 Mb of a segment
  already planted in the same child. The record count of a window comes
  from the source BCF only. This rule keeps segments out of the windows
  next to the centromere and out of the last window.
- Inside a planted segment, haplotype 1 is copied onto haplotype 2:
  every genotype becomes homozygous for the allele of haplotype 1.
  Outside, the genotypes are the source’s. Sites, INFO and sample names
  are copied. INFO is not recomputed.

The planted intervals are the truth. They depend on the seed, the sample
order and the record positions. They do not depend on any ROH call.

### Call sets

An arm’s calls are the segments its macro returns, with the arguments of
`benchmarks/roh_ancestry/run_chr20_arm.R`. The primary call set holds
the calls of at least 500 kb. That is half the shortest planted length,
so a call that trims the ends of a 1 Mb segment still counts. The same
measures are also reported for all calls and for the calls of at least 1
Mb. They are not part of the criterion.

### Measures

Each measure is a ratio of base counts, summed over the children of one
population.

- **Sensitivity** is the planted bases that calls cover, divided by the
  planted bases. It is also reported for each planted length.
- **False discovery** is the unsupported called bases divided by the
  called bases.

A called base outside the planted segments is not always false, because
the source children have real ROH. Each called base therefore falls into
exactly one of four categories, tested in this order:

1.  *planted*: inside a planted segment of the child;
2.  *source low heterozygosity*: in a 100 kb window where the source
    child has at least one called genotype and at most one heterozygous
    call. This is the window rule of the heuristic truth of
    `benchmark_roh_ancestry.md`, without its 1 Mb run length. The
    windows are registry artifact `roh_ancestry_chr20_truth_windows`;
3.  *no record*: in a 100 kb window without any record;
4.  *unsupported*: every other base.

Only the fourth category counts as false. The fraction of called bases
outside the planted segments, whatever their category, is reported next
to it.

### Comparison

The arms are ancestry-tuned, pooled (`INFO/AF`) and single-population
frequencies. The single-population arm uses AFR for ACB, ASW, YRI and
ESN, AMR for CLM, MXL, PEL and PUR, EUR for CEU and EAS for CHS.

Within one population, the ancestry-tuned arm is compared with one other
arm by two differences, ancestry-tuned minus comparator: the sensitivity
difference and the false-discovery difference. Each difference has a 95%
percentile interval from a paired bootstrap over the children of the
population: 2,000 resamples, seed `20261005`, the same resampled
children for both arms.

- Ancestry-tuned is **better** when the interval of the false-discovery
  difference lies below zero and the lower end of the
  sensitivity-difference interval is above -0.01.
- Ancestry-tuned is **not worse** when the upper end of the
  false-discovery-difference interval is below +0.01 and the lower end
  of the sensitivity-difference interval is above -0.01.

The criterion holds when both of these are true:

- in each admixed population (ACB, ASW, CLM, MXL, PEL, PUR),
  ancestry-tuned is better than the pooled arm and better than the
  single-population arm;
- in each control population (YRI, ESN, CEU, CHS), ancestry-tuned is not
  worse than the single-population arm.

A population where an arm has no call in the primary call set has no
false-discovery ratio. Its check is not evaluable and does not count as
a pass. The comparison of the control populations with the pooled arm is
reported and is not part of the criterion.

### Assumption on the ancestry proportions

The ancestry-tuned arm uses the proportions of registry artifact
`roh_ancestry_chr20_q`. They were estimated from the source children.
Planting changes genotypes only inside the segments, 8 Mb of 64.4 Mb per
child, and replaces a heterozygous genotype by a homozygous one for an
allele the child carries. The evaluation assumes that the proportions of
a synthetic child equal those of its source child. They were not
estimated again.

## Data identity

| id                                | locator                                                                                                            | transform                                                            | supplier_identity                                                                                                                                                                                                                                           |
|:----------------------------------|:-------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| roh_ancestry_chr20_children_bcf   | artifact:roh_ancestry_chr20_source; artifact:roh_ancestry_pedigree                                                 | stage_roh_children_from_registry                                     | records=849143; samples=377; records_sha256=a10defdaae8adf4e429baacf43934f51aeb8bd5643b0a304c3f6ce0d4bee673e; bcftools=1.23.1-70-g6dbd8fef; sha256=a53f83cebd7de6fb4661c14ecd93e5346022ff453392095a085e0b7cda116541                                         |
| roh_ancestry_chr20_truth_windows  | artifact:roh_ancestry_chr20_children_bcf                                                                           | read_bcf_count_called_and_heterozygous_GT_floor_POS_minus_1_by_100kb | sha256=77cc5ec52327aa030273d86542e0c1abe47249bae44efb806d5fca478386c83f; children=377; window_rows=243165                                                                                                                                                   |
| roh_ancestry_chr20_reference_long | artifact:ancestry_reference_grch38_parquet                                                                         | filter_chr20_unpivot_21_groups                                       | rows=2332218; sha256=05497926e65b59fe1eeaf2106e1525dfa3c55b0e67117ea99fee4d5cdcfe3bc3                                                                                                                                                                       |
| roh_ancestry_chr20_af_sites       | artifact:roh_ancestry_chr20_children_bcf; artifact:ancestry_reference_grch38_parquet                               | matched_orientation_nonpalindromic_INFO_AF_arrays_first_value        | sites=108757; all_five_AF_fields_nonnull; sha256=fd08cae8ca68d7092a4e5175c382869ac25cb7e1e8276650b8c4988a78e203fd                                                                                                                                           |
| roh_ancestry_chr20_q              | artifact:roh_ancestry_chr20_children_bcf; artifact:ancestry_reference_grch38_parquet; artifact:ancestry_correction | Rduckhts_ancestry_proportions_dosage_min_cor_0.4                     | children=377; sites=108757; all_status=ok; sha256=7449621e7092f311c1e2853431936937a6a33dfc476f5a9f2348a744cfb64448                                                                                                                                          |
| roh_synthetic_chr20_children_bcf  | artifact:roh_ancestry_chr20_children_bcf                                                                           | stage_roh_synthetic_from_registry                                    | seed=20261004; lengths_bp=5000000,2000000,1000000; min_gap_bp=1000000; min_window_records=100; window_bp=100000; records=849143; samples=377; records_sha256=6d8f90531b85db60c551f9149141c0edd8b4530c05a7fed94d84f920db4ae227; bcftools=1.23.1-70-g6dbd8fef |
| roh_synthetic_chr20_truth         | artifact:roh_ancestry_chr20_children_bcf                                                                           | stage_roh_synthetic_from_registry                                    | rows=1131; sha256=ad15f618eb1f80c06761bf7574cb74ce7d36e388405359528f31cc25090567c7                                                                                                                                                                          |

`benchmarks/roh_synthetic_stage.R` stages the synthetic children and the
truth from the registry and refuses to publish them when their record
count, sample count, record digest or truth digest differs from the
registered identity. The record digest is the SHA-256 of the record
lines that `bcftools view -H` prints; the truth digest is the SHA-256 of
the truth CSV. `make test-roh-synthetic-staging` tests the staging
without network access.

## Planted truth

| Population | children | segments |  planted_bp |   records | source_heterozygous |
|:-----------|---------:|---------:|------------:|----------:|--------------------:|
| ACB        |       20 |       60 | 160,000,000 | 2,121,500 |             163,658 |
| ASW        |       13 |       39 | 104,000,000 | 1,449,108 |             115,128 |
| CLM        |       35 |      105 | 280,000,000 | 3,715,362 |             236,678 |
| MXL        |       32 |       96 | 256,000,000 | 3,451,829 |             207,253 |
| PEL        |       35 |      105 | 280,000,000 | 3,723,457 |             189,695 |
| PUR        |       35 |      105 | 280,000,000 | 3,775,681 |             248,458 |
| YRI        |       56 |      168 | 448,000,000 | 5,964,906 |             479,350 |
| ESN        |       43 |      129 | 344,000,000 | 4,520,830 |             353,380 |
| CEU        |       57 |      171 | 456,000,000 | 6,136,355 |             364,822 |
| CHS        |       51 |      153 | 408,000,000 | 5,490,166 |             307,742 |

The truth holds 1,131 segments in 377 children and 3,016,000,000 planted
bases. The segments cover 40,349,194 genotypes, of which 2,666,164 were
heterozygous in the source and are homozygous in the synthetic children.
The full truth table is
[`planted_truth.csv`](results/roh-synthetic/planted_truth.csv).

`benchmarks/roh_synthetic/verify_planting.R` compares the staged BCF
with its source through `read_bcf`, apart from the staging code. All
320,126,911 genotypes of each BCF are in one group per child and
segment, or in the child’s group outside the segments.

| check                                                              | groups |   genotypes | failures |
|:-------------------------------------------------------------------|-------:|------------:|---------:|
| planted segments hold the truth’s record counts                    |  1,131 |  40,349,194 |        0 |
| planted segments replaced the truth’s heterozygous counts          |  1,131 |  40,349,194 |        0 |
| synthetic planted segments hold no heterozygous genotype           |  1,131 |  40,349,194 |        0 |
| synthetic planted segments are homozygous for haplotype 1’s allele |  1,131 |  40,349,194 |        0 |
| genotypes outside planted segments are unchanged                   |    377 | 279,777,717 |        0 |

## Results

Each arm ran once, in a fresh R process with four DuckDB threads and a
12 GB memory limit, on the extension built from the `src` tree of
revision `d0af7e98`, which this branch does not change. The columns
`5 Mb`, `2 Mb` and `1 Mb` are the sensitivity for the planted segments
of that length. The per-child values are in
[`child_metrics.csv`](results/roh-synthetic/child_metrics.csv).

### Primary call set: calls of at least 500 kb

| Population | Arm               | Children | Sensitivity | 5 Mb   | 2 Mb   | 1 Mb   | False discovery | Outside planted | Called bp / child | Unsupported bp / child | Low-heterozygosity bp / child | No-record bp / child |
|:-----------|:------------------|---------:|:------------|:-------|:-------|:-------|:----------------|:----------------|:------------------|:-----------------------|:------------------------------|:---------------------|
| ACB        | Pooled            |       20 | 99.96%      | 99.98% | 99.93% | 99.93% | 7.83%           | 12.25%          | 9,113,579         | 714,004                | 387,503                       | 15,000               |
| ACB        | Single population |       20 | 99.95%      | 99.98% | 99.89% | 99.93% | 6.66%           | 10.33%          | 8,917,409         | 593,573                | 327,503                       | 0                    |
| ACB        | Ancestry tuned    |       20 | 99.95%      | 99.98% | 99.89% | 99.92% | 6.60%           | 10.27%          | 8,911,586         | 587,805                | 327,503                       | 0                    |
| ASW        | Pooled            |       13 | 99.90%      | 99.99% | 99.70% | 99.86% | 11.65%          | 17.23%          | 9,655,525         | 1,124,874              | 492,308                       | 46,154               |
| ASW        | Single population |       13 | 99.84%      | 99.97% | 99.50% | 99.85% | 8.90%           | 13.33%          | 9,214,905         | 820,345                | 384,615                       | 23,077               |
| ASW        | Ancestry tuned    |       13 | 99.89%      | 99.98% | 99.70% | 99.85% | 11.57%          | 17.15%          | 9,645,662         | 1,115,617              | 492,308                       | 46,154               |
| CLM        | Pooled            |       35 | 99.90%      | 99.97% | 99.75% | 99.85% | 21.96%          | 26.60%          | 10,887,746        | 2,390,507              | 496,734                       | 8,571                |
| CLM        | Single population |       35 | 99.89%      | 99.97% | 99.72% | 99.86% | 22.05%          | 26.71%          | 10,904,438        | 2,404,845              | 499,591                       | 8,571                |
| CLM        | Ancestry tuned    |       35 | 99.90%      | 99.97% | 99.74% | 99.86% | 21.86%          | 26.50%          | 10,873,812        | 2,376,494              | 496,734                       | 8,571                |
| MXL        | Pooled            |       32 | 99.94%      | 99.98% | 99.93% | 99.77% | 23.89%          | 28.47%          | 11,176,978        | 2,670,443              | 501,986                       | 9,375                |
| MXL        | Single population |       32 | 99.93%      | 99.98% | 99.89% | 99.77% | 23.83%          | 28.36%          | 11,159,151        | 2,659,633              | 495,736                       | 9,375                |
| MXL        | Ancestry tuned    |       32 | 99.93%      | 99.98% | 99.89% | 99.77% | 23.75%          | 28.33%          | 11,154,715        | 2,648,886              | 501,986                       | 9,375                |
| PEL        | Pooled            |       35 | 99.60%      | 99.97% | 99.93% | 97.08% | 31.48%          | 40.52%          | 13,396,164        | 4,217,202              | 1,211,225                     | 0                    |
| PEL        | Single population |       35 | 99.60%      | 99.97% | 99.92% | 97.09% | 31.18%          | 40.13%          | 13,309,583        | 4,150,415              | 1,191,225                     | 0                    |
| PEL        | Ancestry tuned    |       35 | 99.60%      | 99.97% | 99.92% | 97.09% | 31.16%          | 40.15%          | 13,314,131        | 4,149,319              | 1,196,940                     | 0                    |
| PUR        | Pooled            |       35 | 99.90%      | 99.96% | 99.91% | 99.57% | 15.20%          | 20.41%          | 10,041,063        | 1,526,396              | 522,972                       | 0                    |
| PUR        | Single population |       35 | 99.89%      | 99.96% | 99.89% | 99.56% | 14.27%          | 19.53%          | 9,931,678         | 1,417,173              | 522,972                       | 0                    |
| PUR        | Ancestry tuned    |       35 | 99.87%      | 99.96% | 99.80% | 99.59% | 15.02%          | 20.24%          | 10,017,066        | 1,504,102              | 522,972                       | 0                    |
| YRI        | Pooled            |       56 | 99.73%      | 99.97% | 99.96% | 98.05% | 10.73%          | 17.00%          | 9,611,681         | 1,030,981              | 575,865                       | 26,786               |
| YRI        | Single population |       56 | 99.72%      | 99.97% | 99.93% | 98.03% | 9.64%           | 15.49%          | 9,440,024         | 909,973                | 536,580                       | 16,071               |
| YRI        | Ancestry tuned    |       56 | 99.72%      | 99.97% | 99.93% | 98.03% | 9.54%           | 15.41%          | 9,430,035         | 900,060                | 536,580                       | 16,071               |
| ESN        | Pooled            |       43 | 99.50%      | 99.87% | 99.55% | 97.56% | 12.07%          | 17.92%          | 9,698,474         | 1,170,140              | 533,360                       | 34,884               |
| ESN        | Single population |       43 | 99.24%      | 99.67% | 98.98% | 97.56% | 10.52%          | 16.06%          | 9,457,493         | 994,572                | 496,150                       | 27,907               |
| ESN        | Ancestry tuned    |       43 | 99.24%      | 99.67% | 98.98% | 97.56% | 10.49%          | 16.03%          | 9,454,832         | 991,930                | 496,150                       | 27,907               |
| CEU        | Pooled            |       57 | 99.58%      | 99.97% | 99.35% | 98.12% | 21.65%          | 25.06%          | 10,629,845        | 2,301,702              | 345,833                       | 15,789               |
| CEU        | Single population |       57 | 99.58%      | 99.97% | 99.33% | 98.12% | 21.59%          | 25.03%          | 10,625,856        | 2,294,523              | 349,342                       | 15,789               |
| CEU        | Ancestry tuned    |       57 | 99.58%      | 99.97% | 99.33% | 98.12% | 21.59%          | 25.03%          | 10,625,319        | 2,294,062              | 349,342                       | 15,789               |
| CHS        | Pooled            |       51 | 99.90%      | 99.98% | 99.93% | 99.48% | 25.55%          | 28.18%          | 11,127,283        | 2,843,302              | 291,860                       | 0                    |
| CHS        | Single population |       51 | 99.89%      | 99.96% | 99.92% | 99.48% | 24.58%          | 27.21%          | 10,978,435        | 2,698,055              | 289,155                       | 0                    |
| CHS        | Ancestry tuned    |       51 | 99.89%      | 99.96% | 99.92% | 99.48% | 24.64%          | 27.27%          | 10,988,217        | 2,707,587              | 289,405                       | 0                    |

### All calls

| Population | Arm               | Children | Sensitivity | 5 Mb   | 2 Mb   | 1 Mb   | False discovery | Outside planted | Called bp / child | Unsupported bp / child | Low-heterozygosity bp / child | No-record bp / child |
|:-----------|:------------------|---------:|:------------|:-------|:-------|:-------|:----------------|:----------------|:------------------|:-----------------------|:------------------------------|:---------------------|
| ACB        | Pooled            |       20 | 99.96%      | 99.98% | 99.93% | 99.93% | 48.04%          | 51.17%          | 16,377,784        | 7,867,268              | 498,444                       | 15,000               |
| ACB        | Single population |       20 | 99.95%      | 99.98% | 99.89% | 99.93% | 45.10%          | 47.95%          | 15,363,339        | 6,928,690              | 438,315                       | 0                    |
| ACB        | Ancestry tuned    |       20 | 99.95%      | 99.98% | 99.89% | 99.92% | 44.88%          | 47.74%          | 15,301,708        | 6,866,986              | 438,444                       | 0                    |
| ASW        | Pooled            |       13 | 99.90%      | 99.99% | 99.70% | 99.86% | 48.81%          | 52.98%          | 16,995,632        | 8,295,806              | 661,483                       | 46,154               |
| ASW        | Single population |       13 | 99.84%      | 99.97% | 99.50% | 99.85% | 46.24%          | 49.86%          | 15,928,662        | 7,364,926              | 553,791                       | 23,077               |
| ASW        | Ancestry tuned    |       13 | 99.89%      | 99.98% | 99.70% | 99.85% | 46.60%          | 50.95%          | 16,291,966        | 7,592,746              | 661,483                       | 46,154               |
| CLM        | Pooled            |       35 | 99.90%      | 99.97% | 99.75% | 99.85% | 61.66%          | 64.81%          | 22,710,605        | 14,003,119             | 706,981                       | 8,571                |
| CLM        | Single population |       35 | 99.89%      | 99.97% | 99.72% | 99.86% | 60.19%          | 63.47%          | 21,877,138        | 13,167,298             | 709,838                       | 8,571                |
| CLM        | Ancestry tuned    |       35 | 99.90%      | 99.97% | 99.74% | 99.86% | 60.20%          | 63.47%          | 21,879,719        | 13,172,341             | 706,794                       | 8,571                |
| MXL        | Pooled            |       32 | 99.94%      | 99.98% | 99.93% | 99.77% | 62.52%          | 66.34%          | 23,751,840        | 14,848,827             | 898,464                       | 9,375                |
| MXL        | Single population |       32 | 99.93%      | 99.98% | 99.89% | 99.77% | 61.28%          | 65.21%          | 22,981,596        | 14,082,530             | 895,284                       | 9,375                |
| MXL        | Ancestry tuned    |       32 | 99.93%      | 99.98% | 99.89% | 99.77% | 61.26%          | 65.20%          | 22,974,246        | 14,075,119             | 895,284                       | 9,375                |
| PEL        | Pooled            |       35 | 99.67%      | 99.97% | 99.93% | 97.67% | 63.53%          | 71.70%          | 28,170,950        | 17,896,541             | 2,300,795                     | 0                    |
| PEL        | Single population |       35 | 99.67%      | 99.97% | 99.92% | 97.68% | 62.03%          | 70.51%          | 27,039,935        | 16,773,908             | 2,292,206                     | 0                    |
| PEL        | Ancestry tuned    |       35 | 99.67%      | 99.97% | 99.92% | 97.67% | 61.85%          | 70.38%          | 26,918,337        | 16,648,416             | 2,296,172                     | 0                    |
| PUR        | Pooled            |       35 | 99.90%      | 99.96% | 99.91% | 99.57% | 58.39%          | 61.62%          | 20,823,679        | 12,159,679             | 672,305                       | 0                    |
| PUR        | Single population |       35 | 99.89%      | 99.96% | 99.89% | 99.56% | 56.88%          | 60.21%          | 20,085,236        | 11,424,254             | 669,448                       | 0                    |
| PUR        | Ancestry tuned    |       35 | 99.87%      | 99.96% | 99.80% | 99.59% | 57.03%          | 60.35%          | 20,151,840        | 11,492,400             | 669,448                       | 0                    |
| YRI        | Pooled            |       56 | 99.73%      | 99.97% | 99.96% | 98.05% | 49.63%          | 54.29%          | 17,452,540        | 8,661,890              | 785,815                       | 26,786               |
| YRI        | Single population |       56 | 99.72%      | 99.97% | 99.93% | 98.03% | 45.39%          | 50.15%          | 16,003,501        | 7,263,512              | 746,518                       | 16,071               |
| YRI        | Ancestry tuned    |       56 | 99.72%      | 99.97% | 99.93% | 98.03% | 45.31%          | 50.08%          | 15,979,598        | 7,239,685              | 746,518                       | 16,071               |
| ESN        | Pooled            |       43 | 99.50%      | 99.87% | 99.55% | 97.56% | 50.13%          | 54.63%          | 17,543,508        | 8,793,816              | 754,717                       | 34,884               |
| ESN        | Single population |       43 | 99.24%      | 99.67% | 98.98% | 97.56% | 45.61%          | 50.32%          | 15,980,892        | 7,289,638              | 724,484                       | 27,907               |
| ESN        | Ancestry tuned    |       43 | 99.24%      | 99.67% | 98.98% | 97.56% | 45.51%          | 50.23%          | 15,949,960        | 7,258,725              | 724,484                       | 27,907               |
| CEU        | Pooled            |       57 | 99.58%      | 99.97% | 99.35% | 98.12% | 63.90%          | 66.73%          | 23,944,672        | 15,299,911             | 662,451                       | 15,789               |
| CEU        | Single population |       57 | 99.58%      | 99.97% | 99.33% | 98.12% | 61.58%          | 64.57%          | 22,483,499        | 13,844,376             | 657,131                       | 15,789               |
| CEU        | Ancestry tuned    |       57 | 99.58%      | 99.97% | 99.33% | 98.12% | 61.72%          | 64.71%          | 22,572,281        | 13,931,481             | 658,885                       | 15,789               |
| CHS        | Pooled            |       51 | 99.90%      | 99.98% | 99.93% | 99.48% | 67.50%          | 69.86%          | 26,514,576        | 17,896,282             | 626,172                       | 0                    |
| CHS        | Single population |       51 | 99.89%      | 99.96% | 99.92% | 99.48% | 64.78%          | 67.31%          | 24,443,545        | 15,833,586             | 618,734                       | 0                    |
| CHS        | Ancestry tuned    |       51 | 99.89%      | 99.96% | 99.92% | 99.48% | 64.93%          | 67.45%          | 24,553,919        | 15,943,687             | 619,007                       | 0                    |

### Calls of at least 1 Mb

| Population | Arm               | Children | Sensitivity | 5 Mb   | 2 Mb   | 1 Mb   | False discovery | Outside planted | Called bp / child | Unsupported bp / child | Low-heterozygosity bp / child | No-record bp / child |
|:-----------|:------------------|---------:|:------------|:-------|:-------|:-------|:----------------|:----------------|:------------------|:-----------------------|:------------------------------|:---------------------|
| ACB        | Pooled            |       20 | 98.09%      | 99.98% | 99.93% | 84.96% | 5.05%           | 8.52%           | 8,578,486         | 433,601                | 282,503                       | 15,000               |
| ACB        | Single population |       20 | 98.71%      | 99.98% | 99.89% | 89.94% | 3.72%           | 6.36%           | 8,432,444         | 313,469                | 222,503                       | 0                    |
| ACB        | Ancestry tuned    |       20 | 98.08%      | 99.98% | 99.89% | 84.95% | 3.67%           | 6.32%           | 8,376,209         | 307,117                | 222,503                       | 0                    |
| ASW        | Pooled            |       13 | 95.11%      | 99.99% | 99.70% | 61.49% | 6.95%           | 12.10%          | 8,655,937         | 601,332                | 400,000                       | 46,154               |
| ASW        | Single population |       13 | 95.04%      | 99.97% | 99.50% | 61.48% | 4.15%           | 7.97%           | 8,261,336         | 342,821                | 292,308                       | 23,077               |
| ASW        | Ancestry tuned    |       13 | 95.10%      | 99.98% | 99.70% | 61.48% | 6.85%           | 12.01%          | 8,646,558         | 592,559                | 400,000                       | 46,154               |
| CLM        | Pooled            |       35 | 96.69%      | 99.97% | 99.75% | 74.19% | 12.58%          | 16.70%          | 9,286,064         | 1,168,292              | 373,877                       | 8,571                |
| CLM        | Single population |       35 | 97.04%      | 99.97% | 99.72% | 77.04% | 12.73%          | 16.82%          | 9,333,372         | 1,187,697              | 373,877                       | 8,571                |
| CLM        | Ancestry tuned    |       35 | 97.05%      | 99.97% | 99.74% | 77.04% | 12.47%          | 16.55%          | 9,303,540         | 1,160,140              | 371,020                       | 8,571                |
| MXL        | Pooled            |       32 | 96.84%      | 99.98% | 99.93% | 74.94% | 13.43%          | 16.19%          | 9,243,378         | 1,241,230              | 245,969                       | 9,375                |
| MXL        | Single population |       32 | 96.44%      | 99.98% | 99.89% | 71.82% | 11.96%          | 14.78%          | 9,052,841         | 1,082,683              | 245,969                       | 9,375                |
| MXL        | Ancestry tuned    |       32 | 96.44%      | 99.98% | 99.89% | 71.80% | 12.06%          | 14.94%          | 9,070,077         | 1,093,656              | 252,219                       | 9,375                |
| PEL        | Pooled            |       35 | 97.46%      | 99.97% | 99.93% | 79.97% | 21.40%          | 26.46%          | 10,602,461        | 2,268,897              | 536,997                       | 0                    |
| PEL        | Single population |       35 | 97.10%      | 99.97% | 99.92% | 77.12% | 20.91%          | 26.00%          | 10,497,100        | 2,194,732              | 534,140                       | 0                    |
| PEL        | Ancestry tuned    |       35 | 97.10%      | 99.97% | 99.92% | 77.12% | 21.04%          | 26.29%          | 10,539,553        | 2,217,185              | 554,140                       | 0                    |
| PUR        | Pooled            |       35 | 96.01%      | 99.96% | 99.91% | 68.46% | 7.73%           | 13.07%          | 8,834,990         | 682,767                | 471,543                       | 0                    |
| PUR        | Single population |       35 | 96.01%      | 99.96% | 99.89% | 68.46% | 6.71%           | 12.10%          | 8,738,126         | 586,064                | 471,543                       | 0                    |
| PUR        | Ancestry tuned    |       35 | 95.99%      | 99.96% | 99.80% | 68.48% | 7.77%           | 13.11%          | 8,837,216         | 686,813                | 471,543                       | 0                    |
| YRI        | Pooled            |       56 | 93.94%      | 99.97% | 99.96% | 51.75% | 7.29%           | 12.91%          | 8,629,324         | 629,308                | 458,198                       | 26,786               |
| YRI        | Single population |       56 | 93.48%      | 99.97% | 99.93% | 48.16% | 6.20%           | 11.35%          | 8,436,691         | 522,953                | 418,913                       | 16,071               |
| YRI        | Ancestry tuned    |       56 | 93.48%      | 99.97% | 99.93% | 48.16% | 6.20%           | 11.36%          | 8,436,719         | 523,056                | 418,913                       | 16,071               |
| ESN        | Pooled            |       43 | 96.31%      | 99.87% | 99.55% | 72.05% | 10.41%          | 15.03%          | 9,067,674         | 944,037                | 383,790                       | 34,884               |
| ESN        | Single population |       43 | 96.05%      | 99.67% | 98.98% | 72.05% | 8.98%           | 13.16%          | 8,848,029         | 794,397                | 341,929                       | 27,907               |
| ESN        | Ancestry tuned    |       43 | 96.05%      | 99.67% | 98.98% | 72.05% | 8.84%           | 12.98%          | 8,829,812         | 780,851                | 337,278                       | 27,907               |
| CEU        | Pooled            |       57 | 96.52%      | 99.97% | 99.35% | 73.65% | 10.38%          | 12.04%          | 8,778,897         | 911,216                | 130,048                       | 15,789               |
| CEU        | Single population |       57 | 96.52%      | 99.97% | 99.33% | 73.65% | 9.80%           | 11.51%          | 8,725,908         | 854,992                | 133,557                       | 15,789               |
| CEU        | Ancestry tuned    |       57 | 96.52%      | 99.97% | 99.33% | 73.65% | 9.96%           | 11.66%          | 8,741,124         | 870,285                | 133,557                       | 15,789               |
| CHS        | Pooled            |       51 | 96.77%      | 99.98% | 99.93% | 74.47% | 13.19%          | 14.42%          | 9,046,628         | 1,193,121              | 111,530                       | 0                    |
| CHS        | Single population |       51 | 96.76%      | 99.96% | 99.92% | 74.47% | 12.69%          | 13.91%          | 8,991,617         | 1,140,967              | 109,569                       | 0                    |
| CHS        | Ancestry tuned    |       51 | 96.76%      | 99.96% | 99.92% | 74.47% | 12.71%          | 13.93%          | 8,993,514         | 1,142,864              | 109,569                       | 0                    |

## Criterion

The criterion **did not hold**: 6 of the 16 required checks pass, 10
fail and 0 are not evaluable. Differences are ancestry-tuned minus
comparator, as fractions, with their 95% paired-bootstrap intervals. The
values are in
[`criterion_by_population.csv`](results/roh-synthetic/criterion_by_population.csv).

| Population | Group   | Comparator        | Sensitivity difference       | False-discovery difference   | Better | Not worse | Required  | Pass  |
|:-----------|:--------|:------------------|:-----------------------------|:-----------------------------|:-------|:----------|:----------|:------|
| ACB        | admixed | Pooled            | -0.0001 \[-0.0003, +0.0001\] | -0.0124 \[-0.0378, +0.0010\] | FALSE  | TRUE      | better    | FALSE |
| ACB        | admixed | Single population | -0.0000 \[-0.0001, +0.0000\] | -0.0006 \[-0.0022, +0.0003\] | FALSE  | TRUE      | better    | FALSE |
| ASW        | admixed | Pooled            | -0.0001 \[-0.0002, +0.0000\] | -0.0008 \[-0.0018, -0.0001\] | TRUE   | TRUE      | better    | TRUE  |
| ASW        | admixed | Single population | +0.0006 \[+0.0000, +0.0018\] | +0.0266 \[-0.0000, +0.0722\] | FALSE  | FALSE     | better    | FALSE |
| CLM        | admixed | Pooled            | +0.0000 \[-0.0000, +0.0001\] | -0.0010 \[-0.0060, +0.0049\] | FALSE  | TRUE      | better    | FALSE |
| CLM        | admixed | Single population | +0.0001 \[-0.0000, +0.0002\] | -0.0020 \[-0.0095, +0.0045\] | FALSE  | TRUE      | better    | FALSE |
| MXL        | admixed | Pooled            | -0.0001 \[-0.0003, +0.0000\] | -0.0015 \[-0.0028, -0.0001\] | TRUE   | TRUE      | better    | TRUE  |
| MXL        | admixed | Single population | +0.0000 \[-0.0000, +0.0000\] | -0.0009 \[-0.0050, +0.0017\] | FALSE  | TRUE      | better    | FALSE |
| PEL        | admixed | Pooled            | +0.0000 \[-0.0001, +0.0002\] | -0.0032 \[-0.0069, +0.0007\] | FALSE  | TRUE      | better    | FALSE |
| PEL        | admixed | Single population | -0.0000 \[-0.0000, +0.0000\] | -0.0002 \[-0.0037, +0.0033\] | FALSE  | TRUE      | better    | FALSE |
| PUR        | admixed | Pooled            | -0.0002 \[-0.0007, +0.0001\] | -0.0019 \[-0.0058, +0.0008\] | FALSE  | TRUE      | better    | FALSE |
| PUR        | admixed | Single population | -0.0002 \[-0.0007, +0.0001\] | +0.0075 \[-0.0001, +0.0186\] | FALSE  | FALSE     | better    | FALSE |
| YRI        | control | Pooled            | -0.0001 \[-0.0002, -0.0000\] | -0.0118 \[-0.0245, -0.0024\] | TRUE   | TRUE      | none      | NA    |
| YRI        | control | Single population | -0.0000 \[-0.0000, +0.0000\] | -0.0009 \[-0.0029, +0.0000\] | FALSE  | TRUE      | not worse | TRUE  |
| ESN        | control | Pooled            | -0.0027 \[-0.0062, -0.0002\] | -0.0157 \[-0.0271, -0.0062\] | TRUE   | TRUE      | none      | NA    |
| ESN        | control | Single population | -0.0000 \[-0.0000, +0.0000\] | -0.0002 \[-0.0007, +0.0001\] | FALSE  | TRUE      | not worse | TRUE  |
| CEU        | control | Pooled            | -0.0000 \[-0.0001, +0.0000\] | -0.0006 \[-0.0039, +0.0026\] | FALSE  | TRUE      | none      | NA    |
| CEU        | control | Single population | -0.0000 \[-0.0000, +0.0000\] | -0.0000 \[-0.0024, +0.0022\] | FALSE  | TRUE      | not worse | TRUE  |
| CHS        | control | Pooled            | -0.0001 \[-0.0003, -0.0000\] | -0.0091 \[-0.0130, -0.0055\] | TRUE   | TRUE      | none      | NA    |
| CHS        | control | Single population | +0.0000 \[+0.0000, +0.0000\] | +0.0006 \[+0.0002, +0.0013\] | FALSE  | TRUE      | not worse | TRUE  |

### What the results show

These notes were written after the results. They describe the tables and
change nothing in the protocol.

- Sensitivity does not separate the arms. In the primary call set it
  lies between 99.24% and 99.96% for every population and arm, and the
  largest sensitivity difference between ancestry-tuned and a comparator
  is 0.0027. Every arm finds the planted segments.
- The false-discovery differences in the admixed populations lie between
  -0.0124 and +0.0266. In 2 of the 12 admixed comparisons the interval
  lies below zero. In 0 it lies above zero. In the other comparisons the
  interval contains zero, so this design shows no difference there.
- In ASW and PUR the point estimate of the ancestry-tuned arm is above
  that of the single-population arm, and the upper end of the interval
  is above the +0.01 margin.
- The control check passes in all four control populations: the
  ancestry-tuned arm is within the margin of the single-population arm.
- Called bases outside the planted segments are a large share of every
  arm’s calls: 10.27% to 40.52% of the called bases in the primary call
  set. The window rule excuses part of them. The rest counts as
  unsupported, but nothing here shows that it is false: the source
  children’s ROH between 500 kb and a few megabases is in this category
  when its windows hold two or more heterozygous calls. The planted
  truth is exact for sensitivity only. The false-discovery side still
  rests on the window heuristic.
- With the calls of at least 1 Mb, the sensitivity for the 1 Mb segments
  falls to between 48.16% and 89.94%. The calls that cover those
  segments in the primary call set are shorter than 1 Mb, so they leave
  that call set. This is the effect the 500 kb primary threshold was
  declared to avoid.

## Limitations

- A planted segment has no genotype error: every genotype in it is
  exactly homozygous. A real autozygous segment carries a few
  heterozygous calls from genotype error, so the sensitivities here are
  upper limits for real segments of the same lengths.
- The planted haplotype is a real haplotype of the child, so its allele
  frequencies are those of its population. Copying it makes a segment
  that is autozygous at zero generations; real autozygous segments of 1
  Mb to 5 Mb are older and have had time for mutation and gene
  conversion.
- The lengths are 1 Mb, 2 Mb and 5 Mb. The report says nothing about
  shorter segments.
- The category *source low heterozygosity* is a heuristic on a fixed 100
  kb grid. It uses no allele frequency, so it also excuses calls in
  windows that are homozygous for common haplotypes. A real ROH that
  covers only part of a window with two or more heterozygous calls
  counts as unsupported. The false-discovery ratio is therefore not an
  absolute error rate; the comparison between arms uses the same windows
  for every arm.
- The ancestry proportions are those of the source children (see the
  protocol).
- This is chr20 only, 377 children, one seed. A second seed or other
  chromosomes were not run.
- Every arm ran once. `benchmark_roh_ancestry.md` shows three identical
  repetitions of every arm on the source children; repeatability on the
  synthetic children was not measured again.
- The arms’ site set is the 108,757 sites of
  `roh_ancestry_chr20_af_sites`, while the planting changes all 849,143
  records. The truth is defined on bases, so it does not depend on the
  site set.
- No timing or memory claim is made. `benchmark_roh_scaling.md` holds
  the scaling evidence of the macros.

## Reproducing the evaluation

The drivers run from the repository root and resolve every input through
`r/duckhtsbench`. They load the extension from `build/release/`.

1.  `Rscript benchmarks/roh_synthetic/stage.R` stages the synthetic
    children and the truth and checks their registered identity.
2.  `Rscript benchmarks/roh_synthetic/verify_planting.R` checks the
    staged children against the source.
3.  `Rscript benchmarks/roh_synthetic/run_arm.R ARM` for `pooled`,
    `single_AFR`, `single_AMR`, `single_EUR`, `single_EAS` and
    `ancestry`.
4.  `Rscript benchmarks/roh_synthetic/summarise.R` writes the result
    tables.
5.  Render this report.
