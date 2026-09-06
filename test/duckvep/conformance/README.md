# DuckVEP consequence tests

Ensembl VEP 116 is the behavioral oracle. The Rust prototype and fastVEP are not.

`data/so_consequences.tsv` is generated directly from VEP 116's
`%OVERLAP_CONSEQUENCES` table. `make duckvep-so-spec` verifies the pinned
`Constants.pm` checksum before replacing it; `make duckvep-so-spec-check` proves the
checked-in file is byte-identical. Behavioral rows are produced by executing VEP itself
through the differential runner below.

Variant-induced NMD uses the separate VEP Plugins release/116 `NMD.pm` policy pinned at
`0082591268417af618e03850c5ffdc7c09998a5d`. The pure-C fixed cases reproduce its four
escape rules and executable thresholds on both strands. They also preserve the plugin's
use of full VEP feature endpoints, which can differ from the minimized edit used for the
consequence itself. `NMD_transcript_variant` remains the independent VEP core biotype
consequence.

The runner uses the current `WangLabCSU/blit` command API (tested at
`940c2c1385ba6ad72f0c63b861e90abe8ae6e6f3`) to execute
`micromamba run -p "$VEP_PREFIX" vep ...` through a generated shell script whose dynamic
tokens are individually quoted. The default prefix is
`/root/miniconda3/envs/vep`; override `VEP_PREFIX` for another VEP 116 environment.
Install that `blit` checkout with `R CMD INSTALL /path/to/blit`. A matching VEP
environment can be created with `micromamba create -p "$VEP_PREFIX"
--file test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt` on
Linux-64. Before launching VEP, the runner requires the installed explicit package URL set
to equal that lock exactly; naming a lock file without matching the live environment is
not accepted as provenance.

The optional [`../../../pipelines/duckvep/`](../../../pipelines/duckvep/README.md)
workflow uses `{targets}` for campaign invalidation, branching, resume behavior, and saved
error workspaces. It builds the release extension once and passes a revision- and
byte-bound receipt to every campaign branch. `blit` remains the external VEP/micromamba
process layer; the targets graph does not reconstruct those command lines. Large corpus,
reference, model, cache-info, and cache-receipt paths are tracked inputs, while DuckVEP retains the
semantic receipts and full comparison denominators. The runner performs its explicit
receipt hashes whenever it executes, does not accept precomputed artifact digests, and
does not maintain an independent generic cache.

The validation gates have independent jobs:

- `make test-duckvep-kernel` checks the pure C engine against brute-force interval,
  base-walk, genetic-code, edit-rebuild, and composition oracles. AddressSanitizer and
  UndefinedBehaviorSanitizer have separate targets.
- `make test-duckvep-kernel-statistical` repeats every randomized property 100,000 times.
  `DUCKVEP_PROP_TRIALS` and `DUCKVEP_PROP_SEED` make larger runs reproducible.
- `make test-duckvep-differential` generates boundary, splice, codon, and allele-shape
  witnesses, runs both engines on the same GFF and FASTA, and compares the exact SO term
  set for every `(variant, transcript)` pair.
- `make test-duckvep-gvcf-differential` splits each ALT from a fixed mixed gVCF
  fixture into the same single-allele records given to DuckVEP, then compares
  those records with executable VEP 116. Both `T,<*>` and `<*>,T` source orders
  are present. Literal alleles retain their ordinary consequence; `<*>` receives
  VEP's generic coding or retained start/stop consequence; `<NON_REF>`, bare `*`,
  and `.` do not enter the alternate-overlap comparison. A 171-base REF whose
  `<*>` allele spans the complete transcript pins VEP 116's literal allele-length
  behavior: the executable emits `transcript_ablation`.
- `data/par_path_witnesses.vcf` is the release-116 GRCh38 exact-path PAR witness for a
  published-release/CLI UTR divergence and a sequence-dependent PLCXD1 start-loss/HGVS
  case on both X and Y. Run it through `corpus_differential.R` with the complete human
  cache/model, `--distance 0`, `--hgvs`, and `--sample-per-shape 0`. The retained
  acceptance denominator is 44 exact transcript pairs with no unresolved, missing, extra,
  consequence, HGVSc, or HGVSp difference.
- `--event-mode breakend --regulatory` adds source-derived RegulatoryFeature and
  MotifFeature start/mid/end points in both raw local (`point - 1`, before VEP anchor
  removal) and verbatim mate forms, plus same-object pairs whose local and mate points
  both hit one feature and pairs whose mate hits while the shifted local point is one
  base before it. It also generates shifted local points exactly 5,000 and 5,001 bases
  after a mate-discovered object, proving that VEP's fixed structural-allele admission
  cap is distinct from the caller-configurable transcript distance. This makes
  local-only, mate-only, both-exact, mate-exact/local-close, exact-cap, and beyond-cap
  states deterministic instead of relying on transcript-derived endpoints to hit a core
  interval by chance. Exact duplicate VEP object rows are collapsed and distinct
  allele-level consequence rows are unioned per `(event, object)` before pair
  denominators; conflicting status, reason, or NMD state fails the run.
- Run the same generated breakend corpus with `--distance 0`, `--distance 137`, the
  default `--distance 5000`, and `--distance 10000` before closing a BND semantic change.
  The non-default campaigns prevent the equal numeric defaults from becoming one hidden
  implementation authority.
- `make test-duckvep-state-exploration` runs every C property 100,000 times, then adds
  20,000 deterministic alleles concentrated around transcript boundaries and distributed
  across the transcript, and compares all of them with executable VEP 116. The generated
  VCF records the seed; pair-level Parquet keeps every disagreement and unresolved row.
  `DUCKVEP_STATE_MAX_LENGTH` controls the generated non-anchor allele length and the target
  passes that value plus the retained VCF anchor to the differential eligibility cap. A
  long-allele campaign therefore fails rather than quietly reverting to the runner's
  50-base default subset.
  `data/property_coverage_requirements.tsv` makes the rare-state denominator executable:
  a statistically required state must receive at least one observation, while an allowed
  zero count must name a fixed C witness. An undeclared zero counter, a missing required
  counter, or a statistically required zero fails and preserves the complete property log
  under `results/` with its seed and requested trial count.
  The acceptance path fails after writing those artifacts if it observes any consequence
  mismatch, missing/extra row, or DuckVEP unresolved state. It cannot turn a statistical
  failure into a successful campaign by excluding it from the resolved denominator.

The finite-state campaign manifest uses exact source revision, artifact digest, corpus,
model, oracle version, and oracle-build tokens. Its checker parses the selected CSV rows
and validates their all-row and mismatch counters; it never treats a corpus-name substring
as evidence. `make duckvep-state-current-check` is the stricter release audit: it fails
until every executable campaign named by the state relation has evidence at current HEAD.
Historical successful evidence remains useful for lineage, but cannot satisfy that
current-release gate after implementation or executable-conformance-input changes. A
campaign need not name the later evidence-only commit that records its result: the strict
check instead proves that its source commit is an ancestor of `HEAD` and that the DuckVEP
sources, build/catalog inputs, property harness, vendored property library, executable
conformance inputs, transition relation, and outcome relation are byte-unchanged since
measurement. Evidence histories, generated result artifacts, the campaign receipt
manifest, and documentation are excluded from that comparison so committing a receipt
does not immediately make the receipt stale. The strict check also rejects staged,
unstaged, or untracked semantic inputs; release evidence is evaluated only from a clean,
committed implementation and harness.

Property and distribution histories recover an interrupted journal before enforcing the
clean-worktree gate, then use one publication lock per destination directory,
which prevents concurrent publication of overlapping custom ledger pairs such as `(A, B)`
and `(A, C)`.
Publication is recoverable across process termination between the two POSIX renames;
removing the journal commits the pair before rollback copies are deleted. This is not a
claim of power-loss durability on filesystems that have not persisted file and directory
metadata. The campaign invokes the committed root `Makefile` explicitly, clears inherited
GNU Make control variables and optional property compiler flags, and rechecks all tracked
and untracked source inputs—including ignored compiler inputs in the complete C property
source/vendor closure—after the final evidence input has been consumed.

## Upstream Ensembl test lineage

The Ensembl suites are useful oracle checks, but they do not replace a DuckVEP
differential. Exact source files and per-file SHA-256 receipts are committed under
[`../upstream/`](../upstream/README.md); `perl ../upstream/check_sources.pl` rejects local
drift, omissions, and unreceipted additions without network access. Pin all repositories
to release 116 before using them:

- `Ensembl/ensembl-vep` release/116 commit
  `57ea5c52340acc1f156267f810ad162e26597082` supplies the VEP runner/parser,
  cache, regulatory/motif, structural-VCF, chromosome-alias, mitochondrial, and
  small Haplosaurus tests. Its bundled cache is primarily Ensembl 84 GRCh38 test
  data from small chromosome-21/22 regions; passing it proves that the Perl oracle
  environment is healthy, not that DuckVEP is conformant.
- `Ensembl/ensembl-variation` release/116 supplies the richer semantic authority.
  `modules/t/variation_effect.t` has 189 targeted consequence cases and
  `modules/t/hgvs_parser.t` has exact genomic, transcript, and protein HGVS cases;
  import provenance-labelled cases from those files and then compare both engines on
  the matching model. Do not describe this as DuckVEP directly passing Ensembl's
  internal suite.
- The legacy monolithic VEP suite is preserved in the Git history of
  `Ensembl/ensembl-tools/scripts/variant_effect_predictor`. Its last pre-removal
  snapshot is smaller than the current VEP suite. The older consequence and HGVS
  lineage is already present in `ensembl-variation` Git history, where
  `variation_effect.t` and `hgvs_parser.t` were added in 2011 and 2012. Ensembl
  migrated its CVS development history to Git for release 75; VEP's
  `VEP_SUB_VERSION` and “Subversion update” commit messages mean point releases such
  as 115.2, not the Apache Subversion version-control system. There is no separate
  public SVN conformance authority.

The upstream suites exercise core topology, start/stop, frameshift/in-frame,
incomplete-codon, mature-miRNA, frameshift-intron, selected HGVS shift,
regulatory/motif, structural-parser, mitochondrial-codon, and one small phased-SNV
case. They do not exhaust the named VEP-116 compatibility rules in
`design/duckvep_errata.md`, the 1,000-base shift cap, the full BND/SV predicate
matrix, long/pangenome alleles, or phased indel/MNV edit sets. Keep three separate
receipts: the pinned upstream self-test, extracted upstream semantic fixtures, and
the DuckVEP fixed/property/statistical/corpus conformance campaign.

The state exploration defaults are reproducible and can be widened without changing code:

```sh
DUCKVEP_STATE_CASES=100000 \
DUCKVEP_STATE_SEED=29 \
DUCKVEP_STATE_MAX_LENGTH=49 \
DUCKVEP_PROP_TRIALS=1000000 \
  make test-duckvep-state-exploration
```

A held-out pure-C run with seed `20260719` raised every randomized property to 100,000
cases. The ordinary, AddressSanitizer, and UndefinedBehaviorSanitizer targets each passed
175 tests and 204,759 assertions. That seed exposed an over-broad statistical oracle after
93,064 generated frame-changing cases: VEP's insertion-length-aware terminal-codon
reconstruction is a fallback only for an empty or `X`-containing local alternate peptide.
The minimized concrete-`*` scene now has a fixed regression, while the corrected random
oracle still observed 49 genuine fallback cases and 1,877 stop-gained frame changes.

Once a rare state is discovered, its generator receives a dedicated stratum instead of
depending on its accidental probability under a broad distribution. The terminal missing-
tail, terminal-duplication, and non-stop terminal-codon strata follow this rule; the
coverage manifest or direct property assertion makes their absence a failed run.

The independent executable differential with seed `20260716` wrote 100,268 variants (268
fixed witnesses plus 100,000 generated alleles) and matched all 100,268 VEP-116 transcript
pairs, with no unresolved, missing, or extra rows. The generated alleles included 384 SNVs,
33,032 MNVs, 30,328 insertions, 1,442 deletions, and 34,814 delins. These are measured counts
for that seed, not an assertion that the accepted cases are uniformly distributed:
duplicate rejection exhausts the small SNV/deletion state space sooner than the
longer-allele spaces.

## Registry-backed external chr22 corpora

`r/duckhtsbench/inst/benchmark_registry.tsv` is the sole authority for the HPRC v2
African-four carried-allele, Sniffles2 1KGP, and dbVar GRCh38 chr22 sources and
derivations. Stage all three with:

```sh
Rscript r/duckhtsbench/scripts/stage_duckvep_conformance_corpora.R --corpus all
```

Use `--corpus hprc-african4-chr22`, `--corpus sniffles2-chr22`, or
`--corpus dbvar-chr22` for one lane, and `--plan` to inspect its complete registry
closure. The complete cohort VCFs are read through registered HTTPS range sources;
the stage stores their pinned source indexes, source-identity receipts, deterministic
sites-only chr22 VCFs, tabix indexes, and adjacent provenance under
`DUCKHTS_CACHE_DIR`. The source receipts retain the HPRC S3 object versions, the
Sniffles2 ETag/release, and the dbVar dated publisher-manifest MD5. Derived VCF and
index checksums remain fail-closed.

`scripts/stage_duckvep_conformance_corpora.sh [OUTPUT_DIR] [CORPUS]` remains a
positional compatibility launcher, but it forwards to the R entry point and does not
own locators, identities, paths, or transformations. New callers should use the R CLI
and obtain an individual path from the registry when composing a campaign, for example:

```sh
DBVAR_CHR22=$(Rscript -e '
  Sys.setenv(DUCKHTSBENCH_REGISTRY="r/duckhtsbench/inst/benchmark_registry.tsv")
  source("r/duckhtsbench/R/registry.R")
  cat(duckhts_bench_artifact_path("duckvep_dbvar_grch38_20260127_chr22"))
')
```

Pass that VCF to the differential command below together with its receipt-matched
model, reference, and VEP oracle. Staging establishes source and derivation identity;
it does not by itself establish executable-VEP conformance.

For a large VCF, prepare an ordinary DuckDB database containing
`duckvep_sequence_regions`, `duckvep_transcripts`, `duckvep_exons`, and
`duckvep_transcript_names`. When the model carries Ensembl mature-miRNA
attributes, also provide `duckvep_mature_mirna` with transcript index and
inclusive genomic start/end columns. The runner loads that packed side relation
automatically. Then run, for example:

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

Use the matching indexed VEP cache for a release-level comparison with the
core-dump model:

```sh
make duckvep-corpus-differential DUCKVEP_DIFFERENTIAL_ARGS="\
  --corpus clinvar \
  --vcf /data/clinvar.vcf.gz \
  --cache-dir /data/vep-cache \
  --cache-info /data/vep-cache/homo_sapiens/116_GRCh38/info.txt \
  --cache-receipt /data/receipts/homo_sapiens-116-GRCh38.tsv \
  --fasta /data/GRCh38.fa \
  --database /data/duckvep-model.duckdb \
  --model-sql '' \
  --sample-per-shape 50000"
```

The corpus comparison is spillable rather than proportional to the complete annotation
relation. `--duckdb-memory-limit` (default `8GB`) and `--duckdb-threads` (default `4`)
control the DuckDB evidence-projection phase. The resident DuckVEP model is allocated
outside DuckDB's memory limit; the runner writes the two engines' unique
`(source, variant_id, transcript)` rows to Parquet, drops the resident model and in-memory
annotation tables, and then performs consequence and HGVS full outer joins from those
files. Duplicate pair keys fail before comparison instead of multiplying rows.

Add `--nmd-plugin-dir /path/to/pinned/VEP_plugins` to the same command to run
`NMD.pm`. The runner rejects any `NMD.pm` whose SHA-256 differs from the pinned
release/116 file. The companion `DuckVEPNMDState.pm` only serializes the parent
`TranscriptVariation` CDS coordinates that `NMD.pm` actually reads; this avoids
mistaking the different allele-level CDS range in VEP JSON for the plugin's
state. The runner stores both engines' `triggering`/`escaping`/`unresolved`
prediction on every transcript row and writes a separate
`*_nmd_conformance.csv` confusion table. Recording the run also adds
`nmd_prediction` rows to `data/conformance_history.csv`, which the rendered
report keeps separate from core SO-term agreement. A present consequence row
without an NMD observation remains `not_measured`, even if the other engine
omitted that consequence row. A missing consequence emission is
`not_comparable` only when the emission that does exist contains an eligible
NMD observation (`triggering`, `escaping`, or `unresolved`). Rows marked
`not_applicable` also stay out, so SO misses in older, plugin-free, or
NMD-ineligible rows cannot enter the NMD audit as observations.

`--gff` remains useful for small fixed fixtures and for auditing VEP's GFF
importer. It is not interchangeable with the indexed cache: VEP may skip GFF
feature types or parents that are present in the Ensembl core dump.

Sampling is deterministic within allele type and length-change bin. Set
`--sample-per-shape 0` to retain every eligible biallelic record. Variants are sorted once
before the C engine runs. The output is a Parquet row set from both engines, a pair-level
Parquet difference, and CSV summaries with exact matches, unresolved engine rows, resolved
discordances, emission misses/extras, and exact binomial 95% upper bounds.

Small-event sampling is not SNP-only. It retains four separate shapes: SNV, equal-length
MNV, insertion-like, and deletion-like. The last two include ordinary insertions/deletions
and non-empty delins; the pure-C properties independently generate all five semantic edit
classes, shared prefixes/suffixes, position-one right anchors, same- and cross-codon MNVs,
frame-changing edits, in-frame edits, and phased edit sets.

Use structural mode without `--vcf` to generate events from the loaded model itself. Each
seed samples real transcript-exact/containing/partial spans, coding segments of several
lengths, introns, exon-to-intron spans, both UTRs, start and stop codons, and insertion sites
in CDS, introns, and at exon/CDS edges. Span geometry is crossed with DEL, DUP, tandem DUP,
INV, and undirected CNV; insertion geometry uses INS. VEP 116 remains the only label oracle:
the generator supplies locations and operations, not expected consequences.

```sh
make duckvep-corpus-differential DUCKVEP_DIFFERENTIAL_ARGS="\
  --event-mode structural \
  --corpus sv_chr21_seed17 \
  --database /data/homo_sapiens_116_GRCh38.duckdb \
  --model-sql '' \
  --cache-dir /data/vep-cache \
  --cache-info /data/vep-cache/homo_sapiens/116_GRCh38/info.txt \
  --cache-receipt /data/receipts/homo_sapiens-116-GRCh38.tsv \
  --fasta /data/GRCh38.fa \
  --assembly GRCh38 --species homo_sapiens \
  --chrom 21 --seed 17 --sample-per-shape 100"
```

When `--model-sql` is empty, the runner loads DuckHTS in a private in-memory database and
attaches the prepared model read-only. Independent chromosome/seed jobs can therefore run
concurrently without copying the model or serializing on a DuckDB file lock. Every generated
VCF ID contains its state, source transcript ordinal, coordinates, and operation; retaining
`--sample-vcf` makes any discovered state directly reproducible.

Structural mode passes `--max_sv_size 10000000` to VEP by default and records that value
in the oracle receipt, so the executable does not silently skip an otherwise tested exact
span. This is deliberately different from VEP's own 5,000-base default. That same VEP
option also limits bounded `<CNV:TR>` literal expansion; it is unrelated to DuckVEP's
configurable transcript-distance default and VEP's separate fixed 5-kb BND admission rule.

The first eight-seed GRCh38 campaign covered chromosomes 1, 2, 6, 11, 17, 21, 22, and X:
40,375 generated events produced 2,140,911 transcript pairs, all exact against executable
indexed-cache VEP 116. Treat those denominators as the tested distribution. They do not
extend the claim to paired breakends, raw `<CNV:TR>` repeat reconstruction, producer-
specific symbolic encodings, or other species; add those as explicit event modes and
strata rather than silently widening this one. VEP 116's separate `CIPOS`/`CIEND` source
contract is that it preserves inner/outer coordinates while its registered consequence
predicates continue to use nominal `POS`/`END`.

The checked-in `data/structural_confidence_grch38.vcf` makes that confidence-coordinate
contract executable. It pairs nominal and `IMPRECISE;CIPOS;CIEND` records with identical
`POS`/`END` for CNV, DEL, DUP, tandem DUP, INV, and INS:

```sh
VEP_PREFIX=/opt/vep Rscript corpus_differential.R \
  --event-mode structural \
  --corpus sv_confidence_grch38 \
  --vcf data/structural_confidence_grch38.vcf \
  --database /data/homo_sapiens_116_GRCh38.duckdb \
  --model-sql '' \
  --cache-dir /data/vep-cache \
  --cache-info /data/vep-cache/homo_sapiens/116_GRCh38/info.txt \
  --cache-receipt /data/receipts/homo_sapiens-116-GRCh38.tsv \
  --fasta /data/GRCh38.fa \
  --assembly GRCh38 --species homo_sapiens \
  --chrom 21 --sample-per-shape 0 --fork 4
```

The pinned VEP 116 run emitted 466 transcript pairs. All 466 DuckVEP pairs were exact,
and the nominal/imprecise consequence multisets matched for all six event-kind pairs in
both engines. The sampled oracle VCF retains `IMPRECISE`, `CIPOS`, and `CIEND`; this is a
nominal-coordinate consequence test, not evidence that the uncertainty interval is exact.

Paired BNDs have their own generated mode because one event has two loci and cannot be
represented as one structural span. It crosses same- and cross-chromosome endpoint pairs
with all four bracket orientations and keeps raw ALT and orientation in the sampled VCF:

For source VCFs, the runner uses `duckvep_breakend_geometry()` for mate coordinates
and orientation instead of a separate regular-expression parser. The allocation-free
native parser follows [VCF 4.5 section 5.4](https://github.com/samtools/hts-specs/blob/e821e4f02ae25c2175f9a366edca1322d6a2de72/VCFv4.5.tex).
It preserves the replacement sequence (including retained local bases), exact contig
names and terminal coordinates. Annotation still requires a paired, positive,
model-addressable mate; parsing a single breakend or telomeric zero does not make that
event annotatable. `MATEID`, phase, reference validation and fusion reconstruction remain
separate from ALT syntax. `functions.yaml` defines the public result fields.

`data/breakend_default_witnesses.vcf` retains the seed-31 first-intron-base
counterexample and an ordinary intron control in all four orientations. With the
default minimal model, run `corpus_differential.R --event-mode breakend --vcf
test/duckvep/conformance/data/breakend_default_witnesses.vcf --corpus
bnd_default_witnesses --sample-per-shape 0` from the repo root. All eight pairs
must match, including the mate's default intergenic term in the four splice-site
cases; neither the pair union nor the report may suppress that term.

```sh
make duckvep-corpus-differential DUCKVEP_DIFFERENTIAL_ARGS="\
  --event-mode breakend \
  --corpus bnd_grch38_seed31 \
  --database /data/homo_sapiens_116_GRCh38.duckdb \
  --model-sql '' \
  --cache-dir /data/vep-cache \
  --cache-info /data/vep-cache/homo_sapiens/116_GRCh38/info.txt \
  --cache-receipt /data/receipts/homo_sapiens-116-GRCh38.tsv \
  --fasta /data/GRCh38.fa \
  --assembly GRCh38 --species homo_sapiens \
  --chrom 1,2,7,21,X --seed 31 --sample-per-shape 2"
```

The FASTA index supplies the deterministic VCF chromosome order. More importantly, the
runner forces VEP's BND `buffer_size` to one while retaining one Perl process. VEP 116
inserts mate positions from every buffered record into a coordinate-only interval tree;
larger BND buffers can therefore make one event gain or lose transcripts because of its
neighbors. The isolated five-chromosome run generated 1,004 events and matched all 91,428
transcript pairs, with no disagreement, extra row, or missing row.

Ensembl also publishes release VCFs whose `VE` and `CSQ` fields contain consequences
computed by its variation release pipeline. DuckHTS reads their `Format=...` CSQ header
directly; the full typed record and the narrower release-product projection can be
measured with:

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

Small pinned shards of these official release VCFs are suitable for ordinary CI because
they audit the published Ensembl Variation product without launching Perl VEP. They do not
certify VEP executable compatibility. In release 116, `X:276322 G>A` and `Y:276322 G>A`
have published `VE=intergenic_variant`, while cache-mode VEP with `--distance 0` emits
three path-specific `5_prime_UTR_variant` transcript rows on each chromosome. The
executable/cache combination remains the semantic authority; product differences remain
visible instead of being relabelled as DuckVEP conformance failures. Full scheduled
matrices should obtain the matching receipt-hashed DuckDB model from an external release
store rather than committing multi-gigabyte caches to git. The planned distribution
contract is one manifest entry per Ensembl/Ensembl Genomes release, species, assembly,
model ABI, transcript-filter policy, source-relation hash set, and artifact digest; a
Zenodo record can provide stable versioned storage while every downloaded model still
passes normal receipt/model-open validation.

`release_vcf_differential.R` makes that release-product audit executable. Its first
fail-closed stratum is literal SNVs. The release `VE` field retains
`Consequence|Index|Feature_type|Feature_id`; its zero-based Index maps each consequence to
the original GVF `Variant_seq` and corresponding VCF ALT, so multiallelic records remain
unambiguous. It aggregates the published VE consequence set per ALT/transcript,
runs the public rich `duckvep_annotate(...)` relation against the receipt-matched model, and
reports exact, missing, extra, and discordant transcript/object pairs without treating
the release relation as executable VEP. Both transcript distances are zero because the
variation database dump records overlapping feature consequences, not VEP CLI's optional
transcript flanks.
Do not instead map by CSQ allele text for indels. The Ensembl Variation
release-116 producer at `2fb834b987ede3824e200197a838ce11e91aeb4b` writes a GVF
`Variant_seq` and `Index` before `gvf2vcf.pl` asks `VariationFeature->to_VCF_record` for the
padded VCF alleles; the future non-SNV stratum must reproduce that indexed relation rather
than equate CSQ allele text with a transformed VCF ALT. `gvf2vcf.pl` also stores
`Consequence` in a hash keyed only by allele and feature while constructing CSQ, so repeated
VE terms overwrite one another there; CSQ is a useful typed presentation, not the complete
stored release-product consequence set.

Run it through the repository target so the exact input, model, release, assembly, and
output receipt remain visible in one command:

```sh
make test-duckvep-release-vcf DUCKVEP_RELEASE_DIFFERENTIAL_ARGS="\
  --input /data/homo_sapiens_incl_consequences-chr22.vcf.gz \
  --database /data/homo_sapiens_116_GRCh38.duckdb \
  --release 116 --assembly GRCh38 --chromosome 22 --threads 1 \
  --source-checksum sha256:HEX \
  --output test/duckvep/conformance/results/release_116_chr22_snv.csv"
```

`make duckvep-record-conformance` reruns the real VEP witnesses and records the current
source revision in `data/conformance_history.csv`. Rows include the complete consequence
set, individual SO terms, optional NMD-plugin predictions, VEP impact, allele shape,
unresolved reason, exact Ensembl build, and annotation-artifact hash. The same target runs
VEP with `--hgvs` and records exact HGVSc/HGVSp comparison counts in
`data/hgvs_history.csv`. The pair-level Parquet embeds the clean source revision,
vendored-htslib-distclean release-build binding, and SHA-256 receipts for the extension, model,
reference FASTA/index, source VCF, and exact sampled VCF passed to VEP. The history
writer rejects a diagnostic
artifact, a stale checkout, or non-constant receipts before it updates checked release
evidence.
`make bench-duckvep-throughput` records the sorted
stable-API path in `benchmarks/data/duckvep_throughput.csv`; its checked-in fixture has
one transcript and is not a whole-genome performance claim. Render both views with
`make duckvep-render-reports`. `make duckvep-record-properties` runs the pure-C
randomized suite and records every reported target, seed, trial count, and duplicate
count in `data/property_history.csv` plus named state-distribution counters in
`data/property_coverage_history.csv`. A failed suite writes no successful history row, but
retains its complete seed-specific log. The coverage-requirement manifest prevents a green
suite from silently omitting a declared rare state.

The fixed-event campaign is a closed release-regression gate for the declared independent
small-variant, typed DEL/DUP/tandem-DUP/INV/INS/CNV, paired-BND, transcript, mature-miRNA,
regulation/motif, supported codon-table/SeqEdit, and NMD-plugin surfaces. Those structural
event kinds have executable-VEP differentials. Structural `STR` additionally has pinned
VEP-source semantics, fixed SQL/R adapter tests, and randomized C coverage: an unexpanded
or oversized repeat uses the tandem-duplication gain/insertion fact algebra. Bounded
`<CNV:TR>` reconstruction from repeat metadata remains input preparation before the
literal small-variant path. `CIPOS`/`CIEND` remain relational evidence because strict VEP
116 consequence terms use nominal `POS`/`END`. VEP accepts a finite symbolic vocabulary
and rejects unrecognised kinds such as CPX; this gate does not promise arbitrary symbolic
parsing. Phased multi-record haplotypes and untested releases/species remain outside it. A
newly observed fixed-event mismatch reopens the gate.

The corpus runner currently compares independent alleles. The pure C tests cover phased
edit grouping, same-codon interactions, open frameshifts, and restored frameshifts; a VEP
Haplosaurus differential belongs with the public phased input surface rather than being
faked in this runner.

`make test-duckvep-haplotype-mechanics` runs the narrower executable prerequisite:
the existing pure-C edit application/translation helpers against the pinned Haplosaurus
parser, genomic-to-CDS mapper and transcript container. It generates two-exon transcripts
on both strands, ordinary genomic VCF records, cis/trans diploid carriers and shared
haplotypes. Same-codon substitutions, open/restored frameshifts, in-frame edits and seeded
random edit sets are compared by complete alternate CDS/protein sequence, frame flags,
every CDS contributor, per-sample counts and total carrier counts. For a larger campaign:

```sh
make test-duckvep-haplotype-mechanics DUCKVEP_HAPLOTYPE_ARGS="--cases 1000 --seed 173"
```

The direct-mutation oracle receives the generated original-CDS edit coordinates.
The separate carrier bridge receives genomic VCF alleles, transcript-ranked exons
and the borrowed reference CDS, but no projected edit coordinates. It uses the existing
event preparation and CDS projector once per event, feeds explicit carrier rows through
the native sparse prefix index, resumes between every carrier row, and rebuilds and
translates each occupied event path once. Their complete lane outputs must agree
before comparison with Haplosaurus. `carrier_metrics.csv` records input events/carriers,
peak active slots, completed event paths and translated bases; it is diagnostic work
accounting, not an execution-time or cohort-memory benchmark. Haplosaurus independently
parses and projects the same genomic alleles. The Perl observer
changes only serialization. Clean exact-commit VEP/Variation mirrors and the existing
exact-package environment lock are required. Generated inputs, both engines' complete
observations, mismatch rows, counts and byte receipts remain in a unique directory under
`results/`, including after a failed comparison. These diagnostic receipts do not update
release conformance histories or the `not_implemented` phased transition.
Seven deliberate corruptions exercise each comparison field; their rejection counts
are reported separately from the real engine comparisons.

This is **not** a public phased-executor certificate: it does not test DuckDB carrier
streaming, strict phase/PS interpretation, compound SO/HGVS, structural
composition or arbitrary ploidy. Haplosaurus exposes sequence differences and frame
flags, not a compound SO/HGVS oracle. Its offline container also defaults to two lanes
without inferring VCF ploidy; that behavior must not silently define DuckVEP's ploidy
contract. The existing independent-event and pure-C property campaigns remain separate
and unchanged.
