# DuckVEP consequence tests

Ensembl VEP 116 is the behavioral oracle. The Rust prototype and fastVEP are not.

`data/so_consequences.tsv` is generated directly from VEP 116's
`%OVERLAP_CONSEQUENCES` table. `make duckvep-so-spec` verifies the pinned
`Constants.pm` checksum before replacing it; `make duckvep-so-spec-check` proves the
checked-in file is byte-identical. Behavioral rows are produced by executing VEP itself
through the differential runner below.

The runner uses the current `WangLabCSU/blit` command API (tested at
`940c2c1385ba6ad72f0c63b861e90abe8ae6e6f3`) to execute
`micromamba run -p "$VEP_PREFIX" vep ...` without a shell. The default prefix is
`/root/miniconda3/envs/vep`; override `VEP_PREFIX` for another VEP 116 environment.
Install that `blit` checkout with `R CMD INSTALL /path/to/blit`. A matching VEP
environment can be created with `micromamba create -p "$VEP_PREFIX"
bioconda::ensembl-vep=116.0`.

The tests have three independent jobs:

- `make test-duckvep-kernel` checks the pure C engine against brute-force interval,
  base-walk, genetic-code, edit-rebuild, and composition oracles. AddressSanitizer and
  UndefinedBehaviorSanitizer have separate targets.
- `make test-duckvep-kernel-statistical` repeats every randomized property 100,000 times.
  `DUCKVEP_PROP_TRIALS` and `DUCKVEP_PROP_SEED` make larger runs reproducible.
- `make test-duckvep-differential` generates boundary, splice, codon, and allele-shape
  witnesses, runs both engines on the same GFF and FASTA, and compares the exact SO term
  set for every `(variant, transcript)` pair.

For a large VCF, prepare an ordinary DuckDB database containing
`duckvep_sequence_regions`, `duckvep_transcripts`, `duckvep_exons`, and
`duckvep_transcript_names`, then run for example:

```sh
make duckvep-corpus-differential DUCKVEP_DIFFERENTIAL_ARGS="\
  --corpus giab \
  --vcf /data/HG002.vcf.gz \
  --gff /data/GRCh38.116.gff3.gz \
  --fasta /data/GRCh38.fa \
  --database /data/duckvep-model.duckdb \
  --model-sql '' \
  --sample-per-shape 50000"
```

Sampling is deterministic within allele type and length-change bin. Set
`--sample-per-shape 0` to retain every eligible biallelic record. Variants are sorted once
before the C engine runs. The output is a Parquet row set from both engines, a pair-level
Parquet difference, and CSV summaries with exact matches, unresolved engine rows, resolved
discordances, emission misses/extras, and exact binomial 95% upper bounds.

Ensembl also publishes release VCFs whose `VE` and `CSQ` fields contain every consequence
computed by its variation pipeline. DuckHTS reads their `Format=...` CSQ header directly;
the full typed record and the narrower consequence-oracle projection can be measured with:

```sh
make bench-duckvep-release-parquet DUCKVEP_RELEASE_PARQUET_ARGS="\
  --input /data/homo_sapiens_incl_consequences-chr22.vcf.gz \
  --source-url https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/homo_sapiens_incl_consequences-chr22.vcf.gz \
  --output-dir /data/duckvep-release-parquet \
  --release 116 --assembly GRCh38 --chromosome 22 --overwrite"
```

The VCF and generated Parquet files remain outside git. The benchmark ledger retains both
SHA-256 checksums, record/allele/CSQ cardinalities, exact byte sizes, compression settings,
thread count, and source revision.

`make duckvep-record-conformance` reruns the real VEP witnesses and records the current
source revision in `data/conformance_history.csv`. Rows include the complete consequence
set, individual SO terms, VEP impact, allele shape, unresolved reason, exact Ensembl
build, and annotation-artifact hash. `make bench-duckvep-throughput` records the sorted
stable-API path in `benchmarks/data/duckvep_throughput.csv`; its checked-in fixture has
one transcript and is not a whole-genome performance claim. Render both views with
`make duckvep-render-reports`. `make duckvep-record-properties` runs the pure-C
randomized suite and records every reported target, seed, trial count, and duplicate
count in `data/property_history.csv`; a failed suite writes no history row.

The corpus runner currently compares independent alleles. The pure C tests cover phased
edit grouping, same-codon interactions, open frameshifts, and restored frameshifts; a VEP
Haplosaurus differential belongs with the public phased input surface rather than being
faked in this runner.
