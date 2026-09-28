Somalier find-sites: public 1000 Genomes slice
================

This report measures `rduckhts_somalier_find_sites()` on a retained,
unaltered three-sample 1000 Genomes chromosome-1 VCF slice distributed
by RBCFTools at [commit
9adeaf4](https://github.com/RGenomicsETL/RBCFTools/blob/9adeaf4cfcc3bff40efca6237749fefb53391678/inst/extdata/1000G_3samples.vcf.gz).
The 1,586-byte input has SHA-256
`d5b129eb4f7e431e0e8d9b73229fe3e304ed14863d2cb26c055f631eabb34110`. The
differential also uses the unaltered 9-record, 32,182-byte [gnomAD VCF
fixture](https://github.com/pysam-developers/pysam/blob/79293854c45889bedef09491ee7d20de3dbc76a6/tests/cbcf_data/gnomad.vcf)
(SHA-256
`afb6133e19dca997769ab5c163389a42a28c241243898e513b3c823b2371a73a`). The
binary oracle is Somalier v0.3.4 at
[`ff58fda`](https://github.com/brentp/somalier/tree/ff58fdade8f4f8293d904f10e0a4a13f1fac808d)
(SHA-256
`18717c205a9c4b65d479f1d2cf069a30b047a5378edf544d403d4081f06a3a78`).
Inputs are staged through the artifact registry without changing
physical records.

- Source revision: 57c0c042bd50aa027c71f8a54e2778dbe879cbf9.
- Workload: `--min-AN 6 --snp-dist 100`, default AF 0.15, target AF
  0.48, one DuckDB thread; five repeats of the caller-connection R
  wrapper after staging and connection setup.
- Denominators: 11 source rows, 3 gated candidates, and 3 selected rows.
  Somalier has 3 candidates and 3 final sites.
- Wall time: median 0.018 s, minimum 0.018 s. Process peak resident
  memory: 155.1 MiB (GNU time, includes R, DuckDB and staging; **not**
  isolated query memory).
- Differential on `(region, position, source_ref, source_alt)`: 0
  disagreements, all retained in
  `somalier_find_sites_disagreements.tsv`. The second public slice has 9
  source rows, 0 candidates, and 0 selected rows, equal to the pinned
  binary at `--min-AN 100 --snp-dist 100`.

| gate            | records |
|:----------------|--------:|
| src             |      11 |
| eligible        |      10 |
| pass_snv        |      10 |
| af_an           |       3 |
| annotation_gate |       3 |
| qc_gate         |       3 |
| interval_gate   |       3 |
| gated           |       3 |
| indels          |       0 |
| snps            |       8 |
| indel_clear     |       3 |
| neighbors       |       3 |
| kept            |       3 |
| selected        |       3 |

Source records and candidate gates. Indels and snps are
exclusion-neighbor inventories, not sequential pass counts.

The small public slice verifies the pipeline but does not predict memory
or time on a full population cohort. The nearest existing rendered
Somalier workload, `benchmark_somalier_site_extraction.md`, measures
different inputs and kernels; no comparable full-cohort selection
baseline is available.

Reproduce locally from the repository root with installed `Rduckhts` and
`duckhtsbench`, using the registry-staged pinned executable:

``` bash
/usr/bin/time -f '%M' -o benchmarks/somalier_find_sites_peak_kb.txt \
  Rscript test/scripts/somalier_find_sites_differential.R "$PWD"
cd benchmarks && Rscript -e "rmarkdown::render('benchmark_somalier_find_sites.Rmd')"
```
