Ancestry-tuned allele frequencies for ROH
================

## Evaluation protocol

This report evaluates whether ancestry-tuned allele frequencies reduce
unsupported long runs of homozygosity (ROH) in admixed 1000 Genomes trio
children without changing non-admixed controls. The acceptance criterion
is declared before calculating or comparing any arm:

- In admixed children, ancestry-tuned AF must yield fewer unsupported
  ROH of at least 1 Mb than both pooled and single-population AF.
- Its fraction of frequency-free truth-run length covered must be no
  more than 5% relatively below either comparator.
- In controls, both unsupported-ROH count and truth coverage must remain
  within 5% relative of the single-population arm. For a zero
  comparator, equality is required.

No amendment to the truth threshold has been made. A truth segment is a
merged run of windows (100 kb each) with at most one heterozygous call,
supported when at least 90% of its length meets that window criterion.
Truth calls use all biallelic SNVs in each full source VCF, not just
ancestry-reference sites. The window distribution and separation check
must be reviewed before any ROH-arm comparison; if inadequate, the
threshold may only be amended here before computing arm results. The
primary length threshold is 1 Mb, with 2 Mb and 5 Mb also reported. FROH
is total called ROH length divided by callable autosomal length.

## Data identity and processing status

The phased high-coverage VCF source is
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{N}.filtered.shapeit2-duohmm-phased.vcf.gz`.
The pedigree source is
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt`.
The ancestry reference is registry artifact
`ancestry_reference_grch38_parquet` (GRCh38, 5,012,826 sites).

| Input                     | Identity / SHA-256                                                                     | Status                                                                        |
|:--------------------------|:---------------------------------------------------------------------------------------|:------------------------------------------------------------------------------|
| chr20 phased VCF          | Source file SHA-256 `c2624c726c6cb27e288fe5908621125740379a3d523299dd07c13abf97645e4b` | Verified against the source file; not staged or processed for this evaluation |
| GRCh38 ancestry reference | SHA-256 `4982810aa69b9c71924e9a02dcd454ccedccc34bd8bcf71506ce9752589847fa`             | Available locally; not joined to the trio genotypes                           |
| Pedigree                  | URL above                                                                              | Not staged or checksummed                                                     |

The planned trio counts are ACB 20, ASW 13, CLM 35, MXL 32, PEL 35 and
PUR 35 (admixed), and YRI 56, ESN 43, CEU 57 and CHS 51 (controls).
Single-population AF uses AFR for ACB, ASW, YRI and ESN; AMR for CLM,
MXL, PEL and PUR; EUR for CEU; and EAS for CHS.

No chromosomes were processed. The chromosome-level children-only BCFs,
frequency-free truth inputs, data receipts, ancestry estimates, and ROH
arms have not been generated. The chr20 checksum identifies the supplied
file, but source identity, trio selection, and all requested output
denominators have not been verified by a benchmark run.

## Truth-window separation

Not measured. No heterozygote-window distribution is available, so the
truth definition has not been empirically assessed and no arm comparison
is reported.

## ROH results

No observations are available. All count and length summaries are
per-child quantities and must be aggregated by population only after
retaining each child’s observations.

| Population | Arm               | Unsupported ROH count / child (\>=1 Mb; \>=2 Mb; \>=5 Mb) | Unsupported length / child | Supported count / child | Supported length / child |         FROH | Truth-run length covered |
|:-----------|:------------------|----------------------------------------------------------:|---------------------------:|------------------------:|-------------------------:|-------------:|-------------------------:|
| ACB        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ACB        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ACB        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ASW        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ASW        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ASW        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CLM        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CLM        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CLM        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| MXL        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| MXL        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| MXL        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PEL        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PEL        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PEL        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PUR        | Pooled            |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PUR        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| PUR        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| YRI        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| YRI        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ESN        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| ESN        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CEU        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CEU        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CHS        | Single population |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |
| CHS        | Ancestry tuned    |                                              Not measured |               Not measured |            Not measured |             Not measured | Not measured |             Not measured |

Child-level q estimates and correlation status have not been calculated.
No arm result can establish the acceptance criterion. **Outcome: not
evaluated; neither improvement nor equivalence is supported by this
report.**

## Runtime and memory

No ROH arm was run. Runtime and peak RSS were not measured, no
fresh-process repetitions were performed, and thread count is not
applicable.

| Chromosome     | Pooled runtime / peak RSS | Single-population runtime / peak RSS | Ancestry-tuned runtime / peak RSS |        Threads |
|:---------------|:--------------------------|:-------------------------------------|:----------------------------------|---------------:|
| None processed | Not measured              | Not measured                         | Not measured                      | Not applicable |

The required per-chromosome fresh-process measurements remain
outstanding.
