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
- `test/data/duckvep/hgvs_terminal_multiplication.vcf` retains the seed-27182818
  terminal `CGT>CCC` HGVS counterexample and copy-count controls. Run it with
  `corpus_differential.R --vcf test/data/duckvep/hgvs_terminal_multiplication.vcf
  --corpus hgvs_terminal_multiplication --sample-per-shape 0 --hgvs` on the
  default minimal model; do not add it to or change the frozen random generator.
- `make test-duckvep-projection ARGS='--vep-prefix /path/to/vep --consequences'`
  compares the typed presentation fields and, separately, native SO terms/status/reason
  on the derived transcript fixtures. Every physical VCF record is retained, including
  repeated alleles. Native pair artifacts and `consequence_summary.csv` remain available
  on failure; passing presentation fields alone do not certify native consequences.
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
`ERRATA.md`, the 1,000-base shift cap, the full BND/SV predicate
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

The native phased-protein campaign checks complete CDS replay, independent translation
and application of the typed protein operations. Its in-frame generator varies two to
five edits in a 36-base CDS. A separate reference-restoration stratum covers 1,344 cells:
seven coding locations, twelve distinct codon pairs, two transcript strands and all eight
orientations of three physical alleles. Every cell receives at least
`floor(DUCKVEP_PROP_TRIALS / 1344)` cases with seeded flanking sequence. These cases must
preserve changed local coding spans even when the complete CDS is unchanged. These two
generators use standard-table Ala/Gly/Asp/Pro body codons and intact start/stop codons.
A separate restored-frame generator samples twelve-base spans, conditioned on stop-free
reference/alternate translation with at least two differing runs separated by an unchanged
residue. It reports accepted cases, all proposals and rejected proposals separately;
every accepted case must split its protein operations and reproduce the complete protein.
The operations retain shared physical-edit provenance. A curated-reference stratum
requires draws in 80 cells: ten peptide-edit positions × two strands × four routes
(no DNA edits, ordinary in-frame edits, complete CDS restoration and frame restoration).
The reference has a single `U` peptide edit; the alternate is translated from literal
CDS replay. Every result must reconstruct that complete alternate protein from the
curated reference. Protein-only differences carry no invented physical source edit.
Every curated-reference draw also compares a borrowed prepared-protein operand
with the sparse model-edit view, checking complete replay, HGVS text, physical
spans and coding-context immutability. This paired check retains the original
case denominator; `prepared_views` counts the additional operand comparisons.
A terminal-reference property requires draws in 504 cells: three raw terminal
stop codons × 21 curated replacement residues × two strands × four edit routes.
The routes are no DNA edits, ordinary in-frame changes, complete CDS restoration
and frame restoration. Reference preparation adds a terminal marker after the
curated residue; the property checks complete edit replay and then compares both
displayed first-stop prefixes. Source spans and context immutability remain checked.
Generated/evaluated cases, passing cases and per-cell coverage are reported even
when a trial fails. A fixed cancellation witness retains a physical insertion
with equal prepared/alternate proteins; its repeat-deletion control must still
produce a deletion. Both strands and malformed physical-block rejection are checked.
This is sequence-mechanics evidence, not complete executable-VEP HGVS conformance.
The conditional frame-restoration route retains its proposal/rejection counts.
A separate terminal-repeat campaign targets `XY` motif insertion with a terminal
`X` reference peptide edit. Its 4,800 cells cross all twenty-by-twenty standard
amino-acid pairs, three stop codons, two transcript strands and two source-allele
strands. Seeded draws vary motif copy count, insertion placement and flanking
sequence. These stop-free bodies require exact complete-protein replay before any
first-stop display rule, so an inserted extra stop cannot be hidden by truncation.
All cells and failures are retained; this cohort has no conditional rejection.
The finite shifted-frame matrix crosses all 4,096 adjacent codon pairs, four
insertion positions, four inserted bases, a preceding in-frame insertion or deletion,
and both allele orientations: 262,144 cases. Each case also retains a left anchor
and a right anchor in both source replacements, yielding 524,288 equivalent-CDS
comparisons. Each anchored operation is also compared with its own isolated source
edit: complete HGVS text, shape, reference positions and termination distance must
agree. Sequence replay remains identical across equivalent source anchors. The matrix
counts 3,072 anchor-dependent presentation differences; none bypasses the per-record
operation checks. [The errata](../../../ERRATA.md#equivalent-dna-anchors-can-produce-contradictory-vep-116-protein-descriptions)
separates that pinned-VEP contract from HGVS nomenclature correctness.
These campaigns do not certify
complete HGVS nomenclature or executable-VEP compatibility. The all-property runner
retains its full denominator and the same failure controls.

The source-record HGVS differential invokes unmodified pinned VEP 116 at buffer sizes
1 and 5,000. Its complete matrix crosses three stop codons, all 64 following codons,
two transcript strands, four insertion positions, four inserted bases and both
retained-anchor forms: 384 models and 12,288 physical records. Literal right-retained
replacements remain unnormalized parser inputs. Both decoded strict calls and raw
`1|1` calls run at one and four DuckDB threads. Each query contains one source per
transcript, so raw-file duplicate retention cannot change the independent-record
comparison. Every input record participates in every route/thread configuration.
The same indexed FASTA and model also pass through `duckvep_annotate` as an independent-
event control, with a separate complete HGVSp comparison and retained output table.

```bash
Rscript test/duckvep/conformance/hgvs_anchor_differential.R --tail-codons 1
Rscript test/duckvep/conformance/hgvs_anchor_differential.R --tail-codons 64
```

`VEP_PREFIX` or `--vep-prefix` names the VEP environment. The smaller command keeps
all six stop/strand combinations and 192 records. Complete oracle VCFs, generated
models and records, nested native outputs, per-record comparisons and source hashes
are retained in a unique results directory. VCF percent escapes are decoded and the
prediction parentheses are removed only for string comparison; absent HGVSp remains
absent. The comparison routine must reject dropped records, duplicate IDs, changed
strings, missing HGVSp and invented HGVSp. Complete native outputs must agree between
thread counts. CDS replay and per-record HGVSp have independent verdicts. These
diagnostic runs fail on any HGVSp disagreement, including cases where VEP emits no
string. A passing terminal-anchor matrix does not certify compound or complete phased
HGVS, and these diagnostic runs do not enter release history.

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
automatically. Registry-built models instead store `model_regions` and nested
`model_transcripts`; shared read-only SQL projections expose the same flat
relations, including every exon, mature-miRNA segment and peptide edit. They do
not rewrite the artifact or recompute biological fields. Model-load queries
refer directly to the attached catalog so they also work from separate DuckDB
connections. The FastVEP benchmark worker uses the same projections.

Then run, for example:

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

The exact repeat preparer has a seeded SQL differential against base R's ordered
string repetition. Run from the repository root:

```sh
Rscript test/duckvep/conformance/repeat_sequence_differential.R --trials 100000 --seed 173
Rscript test/duckvep/conformance/repeat_sequence_differential.R --trials 100000 --seed 20260906
```

Seven required strata exercise exact multi-component sequences, summary-only evidence,
missing counts, missing units, fractional counts, empty sequences and unavailable component
lists. The run retains all input components and complete expected/observed rows in Parquet,
checks 32 generated capacity failures with recovery when enough exact scenes are available,
and rejects dropped rows, duplicate identities, changed sequence and changed status.
Receipts bind the seed, R RNG, extension bytes and source snapshots. These diagnostic runs
do not enter release history or claim VEP raw-STR parsing or biological population coverage.
SQL/R regressions separately compose exact pure and interrupted repeats through annotation
and haplotype inputs; summary metadata is not an exact-sequence oracle.

The native `breakend_parser_rejects_mutated_components` property generates a name,
replacement and 64-bit mate coordinate, then exercises 72 fault/form cells for each
draw: sixteen faults in each paired form and four in each single-breakend form.
Every unmodified form must parse correctly. Mutations target non-DNA replacement
bytes, embedded NUL, extra/missing/mismatched brackets, malformed mate names,
empty or nondecimal coordinates and UBIGINT overflow. Exact-length unterminated
buffers, unchanged input bytes, zeroed error output, output canaries and valid-record
recovery are checked. Coverage reports all attempted and passing mutations plus the
minimum draws per cell, including on failed runs. The 72 mutations share one generated
scene; they are not independent biological samples or evidence for fusion reconstruction.

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

`make test-duckvep-haplotype-mechanics` builds the release extension and runs the
narrower executable prerequisite:
the existing pure-C edit application/translation helpers against the pinned Haplosaurus
parser, genomic-to-CDS mapper and transcript container. It generates two-exon transcripts
on both strands, ordinary genomic VCF records, cis/trans diploid carriers and shared
haplotypes. Same-codon substitutions, open/restored frameshifts, in-frame edits and seeded
random edit sets are compared by complete alternate CDS/protein sequence, frame flags,
every CDS contributor, per-sample counts and total carrier counts. For a larger campaign:

```sh
make test-duckvep-haplotype-mechanics DUCKVEP_HAPLOTYPE_ARGS="--cases 1000 --seed 173"
```

The campaign also loads the built extension (`--extension`) and reads the generated
VCF through `read_geno()` and `read_bcf_samples()`. This third path derives allele
slots, sample identity and variant alleles from decoded calls, then runs the same
native projector/carrier/rebuild pipeline. Both generator-fed paths remain separate
checks, and the existing Haplosaurus verifier controls are unchanged. Six additional
routing controls require valid changed allele slots/ALT to change the haplotype and
missing/duplicate calls, a different model chromosome or reversed record order to fail.
The genotype path retains reader ordinals instead of sorting incorrect input. The complete typed
calls, metrics and controls are retained with the extension's SHA-256. Supplying
`--extension-receipt` validates the
existing clean-revision release-build receipt; runs without one are explicitly
`diagnostic_unbound`, not source-bound release evidence.

The direct-mutation oracle receives the generated original-CDS edit coordinates.
The separate carrier bridge receives genomic VCF alleles, transcript-ranked exons
and the borrowed reference CDS, but no projected edit coordinates. It uses the existing
event preparation and CDS projector once per event, passes the complete decoded diploid
call through the native phase reducer and sparse prefix index, resumes between calls, and rebuilds and
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

`haplotype_sql_differential.R --artifacts results/<haplotype-run>` is an additional
public-SQL lane. It verifies the original artifact hashes and successful controls,
loads the same model/genotypes, and runs `duckvep_haplotypes` under both named policies.
Required occupied carrier keys are checked against the complete independent paths,
not inferred from the rows that happened to survive. Reference lanes are explicitly
implicit in this SQL contract; their source-model paths complete the six-lane comparison
against the unchanged Haplosaurus observation. Full CDS/protein, frame flags, contributors,
sample/carrier counts and PS/ploidy are checked. Six additional output mutations must be
rejected. New outputs and receipts live in a separate result directory; original inputs,
generators, seeds, eligibility, denominators, failures and oracle remain untouched.
`--extension-receipt` applies the same clean-build binding as the original lane.
This tests the declared literal, diploid sequence-mechanics subset, not combined SO/HGVS,
structural composition, raw-parser emulation or broad missing/ploidy compatibility.

`haplotype_phase_differential.R` audits every GT over `0`, `1`, `2`, `.` at
ploidies 1–4, every intervening `/`/`|` combination, and absent, `/`, or `|`
leading prefixes: 7,020 profiles, without sampling. Each profile has a multiallelic
site and a homozygous second site in another PS. The registered 180-base replay
fixture supplies the reference; bcftools checks every REF before either engine runs.
The Haplosaurus observer and public SQL consume the same VCF/GFF/FASTA.
Complete CDS/protein multisets, source-record contributors and carrier counts are
compared; all raw observations, differences, decoded-GT collisions and receipts remain
in a separate result directory. Eighteen ordinary called diploid profiles and twelve
deliberate corruptions guard the verifier. This is a raw-input compatibility audit,
not a replacement for the original phased-sequence corpus or a population error rate.

An optional observer sidecar exposes the actual upstream retained genotype objects
and file-profile ploidy without changing its sequence output. Genotype and mapping
values are copied when the original upstream container constructor returns;
VariationFeature objects can be shared and remapped by another transcript before
output serialization. `replay_lanes` copies each sample's original mutator return
before equal-sequence grouping, including CDS, protein and applied-source keys.
It observes only samples entering that mutator; reference-only samples handled
by upstream's separate reference-haplotype path are not synthesized.
`Rscript test/duckvep/conformance/haplotype_observer_contract.R`
checks both strands against overlapping full-exon/two-exon models and verifies that
enabling the sidecar preserves complete oracle output. The same audit compares
the standalone native raw-GT parser on all 14,040 source calls: retained/omitted status,
parsed slot count, the two consumed allele ordinals, source ploidy and missingness.
Seven field-corruption controls guard that comparison. Parser results and failures
are retained separately; parser agreement does not turn full-replay failures into passes.
The native record-replay lane consumes the parser result and complete source alleles.
It compares CDS/protein multisets, counts and physical-edit source records with
Haplosaurus, and independently checks all 28,080 per-lane record observations against
the upstream genotype sidecar. Missing REF/omitted calls have conditional no-op
observations; undefined slots have conditional full-REF deletions. Four deliberate
observation corruptions guard source identity, evidence and sequence status. The
public decoded-call lane has its own verdict; its failures remain in the receipt.
The public `source_records` lane reads original GT text and complete ALT lists from
the same VCF fixture and checks full sequences, counts, physical-edit provenance
and per-lane record observations separately. Exact record IDs, regions, positions
and REF/ALT bytes are checked against the input, not just the repeated local site
labels. This two-site, single-exon grammar does not certify overlapping source
replacements, splicing, whole-haplotype SO/HGVS or broad phase compatibility.

```sh
Rscript test/duckvep/conformance/haplotype_phase_differential.R --max-ploidy 2
Rscript test/duckvep/conformance/haplotype_phase_differential.R
```

`--max-ploidy 2` is a 108-profile smoke test; the full command enumerates all 7,020
profiles in its finite grammar. Both preserve all comparisons and exit nonzero on
any disagreement. `--extension-receipt` requires clean-build evidence. Typed GT cannot
certify byte-level VEP parser behavior when distinct raw spellings decode identically.

`haplotype_record_differential.R` checks source-record geometry independently of
the raw-GT grammar audit. It uses the registered 180-base CDS, 144 fixed profiles
(12 record geometries × 6 GT patterns × 2 strands) and a declared seeded set of
overlapping pairs. Every profile has a separate homozygous anchor so both file
lanes are occupied. Full REF spans, equal-position file order, retained MNV bases,
insertions/deletions and known REF slots are part of the input, not normalized away.
Haplosaurus receives the exact VCF/GFF/FASTA; bcftools validates every REF first.
The public raw-record executor consumes the matching source relation. Full
CDS/protein multisets, carrier counts and applied-source identity sets are compared
with no excluded disagreements. The construction-time observer also captures both
sample/file lanes before upstream sequence grouping. Each lane's complete CDS,
protein and applied-source allele/record identity set is compared through the
shared lane comparator. These source sets do not certify physical-edit multiplicity.
Counts distinguish unavailable paths from sequence differences where every path
is available. Four corruptions guard the grouped comparator; eleven additional
controls guard lane identity and content, including a lane swap that leaves grouped
results unchanged. The multi-sample campaign retains its twelfth shared-source
control. Both campaigns check the complete native transcript domain, positive
carrier counts and exact list lengths before grouping output. Eight shared
controls reject extra/NULL transcript rows, zero-carrier rows, missing transcripts,
wrong/missing counts, duplicate transcript domains and total-count mismatches.
Eight group controls additionally check empty source sets, missing/null fields,
invalid identities, lost or invented sources, and dropped or duplicated rows.
Disjoint and adjacent records are positive controls.
`summary.csv` retains the grouped `equal` verdict and separately reports
`replay_lanes_equal` and `counts_equal`; `all_equal` requires all three. Receipts
retain every complete observation. A passing grouped comparison cannot hide
a failed lane or total-count comparison.

```sh
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 173 --random-cases 512
```

For mass rare-configuration trials, `--rare-per-stratum` requires that many draws
in every geometry × strand × GT-pattern/source-ploidy cell (1,128 cells). GT patterns
cover pipe/slash calls, mixed separators, leading separators, first/last missing
slots, all-missing calls and ALT only after the two consumed file lanes. Source
ploidies are 1, 2, 4, 8, 16 and 64; patterns require enough slots to express them.
Positions, lengths and inserted/replacement bases vary within each cell. A pair
shares its GT-pattern/ploidy class; callable allele slots are seeded within that
pattern. This does not enumerate every pairing of different classes.
The generator checks source ploidy and spelling, and `coverage.csv` must meet every
declared quota. Haplosaurus still consumes its two-lane file profile; source ploidy
is not a claim of arbitrary-ploidy output compatibility.

```sh
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 173 --rare-per-stratum 32
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 20260906 --rare-per-stratum 32
```

Each command includes the 144 fixed and 512 general-random cases: 36,752 profiles,
110,256 source records and 73,504 file lanes. Quotas establish cross-product
coverage, not exhaustive sequence coverage or independent biological observations.

`--pair-per-stratum` generates ordered pairs from all 47 eligible GT-pattern/ploidy
classes on both strands: 4,418 cells. Each record independently draws its callable
allele slots from its declared class. Geometry cycles through all twelve shapes,
so twelve draws per cell cover every geometry once. Source ploidy may change
between records; the homozygous anchor occupies both Haplosaurus file lanes.
Positions, lengths and replacement bases are seeded. The existing fixed, random,
same-class and source-context cohorts keep their seed stream and precede this cohort.

```sh
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 173 --pair-per-stratum 12
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 20260906 --pair-per-stratum 12
```

Each command has 53,016 paired profiles plus 656 fixed/general-random profiles:
161,016 source records and 107,344 file lanes. `pair_coverage.csv` and
`pair_geometry_coverage.csv` require the complete declared cross-product; input
checks verify both records' ploidy, separators and missing-slot placement. Five
GT corruptions supplement the four sequence/count/provenance controls. Complete
disagreements remain failures, including disjoint pairs with differing GT classes.
The 65,536-profile limit applies to the sum of all enabled cohorts. This finite
grammar does not establish arbitrary-ploidy output compatibility or population
error rates; every cell still has a much larger unsampled sequence space.

`--context-per-stratum` appends paired source-context trials. Identical edit
geometry, alleles and GTs are replayed with `0|0` records before, between or after
the edits. The 2,376 cells cross six geometries, six GT patterns, both strands,
three placements and eleven neutral-record counts through 36. Counts include
the neighbourhoods of the pinned interval tree's root changes. Each draw shares
its tested edits across all 33 contexts; these are paired observations, not
independent biological samples. The original fixed/general-random and optional
rare-GT cohorts retain their inputs and comparison rules.

```sh
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 173 --context-per-stratum 8
Rscript test/duckvep/conformance/haplotype_record_differential.R --seed 20260906 --context-per-stratum 8
```

Each command has 19,008 context profiles plus the 656 fixed/general-random profiles.
`context_coverage.csv` checks every quota; `context_summary.csv` retains complete
verdicts and paired changes relative to zero neutral records. The optional oracle
sidecar records the complete source buffer and each retained genotype's selected
CDS mapping. Buffered record IDs, coordinates and full alleles are checked against input,
and reference-only records must be absent from retained genotype objects. The
sidecar changes serialization only; sequence/count/provenance disagreements remain
failures. This single-exon context experiment does not certify multi-transcript
buffering, cross-exon record selection or complete phased annotation.

Use `--extension-receipt` for clean-build evidence. All inputs, outputs, comparisons
and failures are retained in the reported artifact directory; a mismatch exits
nonzero. This is a single-exon source-replacement audit, not whole-haplotype SO/HGVS,
splice prediction, population error rates or a claim that upstream behavior is wrong.

`haplotype_model_differential.R` gives each source region overlapping single-exon
and two-exon transcripts, with three diploid samples. Its 768-region baseline
preserves the seeded multi-exon diagnostic's input columns. Additional trials
require a quota in every one of 264 cells: four record geometries × three cohort
GT patterns × both strands × eleven neutral-record counts through 36. Positions
and spanning-record lengths vary within each cell; transcript geometry is fixed
at exons 11–70 and 101–190. Each region supplies two transcript comparisons.

```sh
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 173 --rare-per-stratum 32
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 20260906 --rare-per-stratum 32
```

Each command has 18,432 transcript cases, including 1,536 baseline cases, and
110,592 sample/file-lane observations. Counts are compared per sample and complete
CDS/protein group, so exchanging samples cannot pass through a pooled count.
Each sample/file lane separately compares complete CDS/protein and applied-source
identity sets against the upstream mutator. A homozygous alternate anchor puts
all six lanes per transcript through that observed path. Source keys retain the
full REF/ALT list and selected allele in transcript orientation; repeated physical
edit islands from one source remain a single identity in this comparison.
Native contributors separately retain exact
record identity, region, position, REF, interpreted ALT and occupied sample/lane.
The original source buffer, retained genotype multiplicity and transcript-owned
CDS coordinates are checked against the fixture, including exon-repeated and
unselected duplicate sources. Forty-two controls reject changed sequences,
reference models, samples, counts, source provenance, genotype observations, mapping coordinates,
lane identities and malformed observations. Lane swaps and missing lane-specific
sources must fail even when grouped sequences, sample counts and source sets agree.
The sidecar also records upstream's constructed reference CDS. Its comparison
to the loaded fixture has a separate failure count and a changed-sequence control.

`coverage.csv` requires every declared quota; `summary.csv` and `comparisons.rds`
retain all verdicts. The default quota is one; zero runs just the baseline.
`--extension-receipt` enforces clean-build binding. All sequence, provenance and
mapper or lane-association disagreements exit nonzero. This provides shared-transcript and
source-context coverage, not arbitrary exon/UTR geometry, phase/PS inference,
combined SO/HGVS, independent biological observations or population error rates.

`--geometry-per-stratum` adds generated exon/UTR models to the same executor and
comparison gates. It preserves the fixed-model input columns and seed stream,
then draws each extra model from 504 cells: 2/3/5/7 coding exons × first coding
split phase 0/1/2 × absent, intra-exon or separate-exon UTRs × both strands ×
seven source geometries. Coding exon lengths, intron lengths and intra-exon UTR
lengths vary; short internal exons can split a codon across three exons. Source
geometries cover coding substitutions, CDS start/end crossings, exon entry/exit
crossings, a whole exon with flanks and exon-end anchored insertions. Three samples
carry opposite lanes and a compacted missing call, with a homozygous ALT anchor.

```sh
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 173 --rare-per-stratum 32 --geometry-per-stratum 32
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 20260906 --rare-per-stratum 32 --geometry-per-stratum 32
```

Each command adds 16,128 models to the 18,432 fixed-model cases. All 34,560 cases
compare six sample/file lanes, source identities, input provenance and owned
mapper coordinates; all 42 corruption controls must pass. `geometry_coverage.csv`
checks the requested cell quotas. This grammar uses one registered 180-base CDS
and standard-code complete translations. It does not cover arbitrary biological
models, reference-only sample routes, general phase/PS inference or structural
composition. Failures remain in the combined summary and receipt.

`--interaction-per-stratum` appends a second edit to generated exon/UTR models.
Its 6,048 cells cross the 504 geometry cells with four partner allele shapes
(SNV, insertion, deletion and replacement) and three start locations: inside
the first record's REF span and transcript, within the same coding exon, or in
another coding exon. REF spans may cross exon–intron junctions or CDS endpoints;
every source span must intersect its explicitly selected transcript. Three
samples supply cis, trans and compacted-missing calls, with a homozygous ALT
anchor. The source records, including duplicate alleles, remain distinct.

```sh
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 173 --rare-per-stratum 0 --geometry-per-stratum 1 --interaction-per-stratum 4
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 20260906 --rare-per-stratum 0 --geometry-per-stratum 1 --interaction-per-stratum 4
```

Each command adds 24,192 interaction models after the original 1,536 fixed-model
and 504 geometry cases, for 26,232 models and 157,392 sample/file lanes.
`interaction_coverage.csv` enforces four draws in every cell. The original
cohorts retain their inputs and RNG stream; interactions consume later draws.
All sequence, count, reference-model, source, mapping and file-lane comparisons
use the existing gates. Empty applied-source arrays and native character vectors
represent the same empty identity set. Eight shared group controls reject absent
or null contributor fields, invalid identities, lost or invented sources, and
dropped or duplicated rows. This remains a standard-code 180-base CDS grammar,
not a certificate for arbitrary models, phase sets, consequences or HGVS.

`--length-per-stratum` crosses the 504 exon/UTR geometry cells with CDS lengths
36/37/38, 2,047/2,048/2,049 and 6,143/6,144/6,145 bases: 4,536 required cells.
These include short coding sequences, complete and partial terminal codons, and
sequence sizes around 2,048 bases and 2,048 codons. Each CDS retains the registered
reference's ATG start, repeats its internal codons to the requested length minus
six, and appends TAA. A raw TAA suffix is not necessarily an in-frame stop.
All derived sequences and source records are retained in the input receipt.

```sh
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 173 --rare-per-stratum 0 --length-per-stratum 1
Rscript test/duckvep/conformance/haplotype_model_differential.R --seed 20260909 --rare-per-stratum 0 --length-per-stratum 4
```

The commands cover 6,072 and 19,680 models respectively, including the 1,536
fixed-model cases, with 36,432 and 118,080 sample/file lanes. Length draws follow
any geometry/interaction cohorts without changing their inputs or seed stream.
`length_coverage.csv` enforces every requested quota. Complete sequence, count,
source, mapping, reference-model and lane checks retain all failures under the
same 42 corruption controls. The network-free generator test independently checks
genomic REF, spliced CDS, ranked exon phases and emitted GFF coordinates/phases.
This is a standard-code synthetic length grammar, not arbitrary genetic-code,
phase-set, compound consequence or HGVS certification.

`haplotype_reference_differential.R` exercises reference-only sample handling
through the original container JSON serializer. CDS and protein groups are compared
separately by sample and count; one CDS may link both a curated reference peptide
and a mutation peptide. The observer must preserve the complete unobserved JSON,
including array order. Mutation lanes compare complete sequences and applied sources;
all native contributors retain checked source geometry, allele and evidence.

```sh
Rscript test/duckvep/conformance/haplotype_reference_differential.R --seed 173 --rare-per-stratum 32
Rscript test/duckvep/conformance/haplotype_reference_differential.R --seed 20260906 --rare-per-stratum 32
```

Each command requires 32 draws in 360 cells: three start codons × internal-stop
presence × terminal-stop presence × strand × three all-missing GT spellings ×
five routes: missing-only, retained-REF-lane, later-retained-call, short-intronic
and long-intronic calls. Internal-stop codons, stop positions and two source
positions vary. Intronic cases have two exons and lengths sampled from 1–12 or
13–120 bases; the other routes use one exon. Models derive from the registered
180-base reference and use the standard genetic code. Twenty-nine
corruption controls guard groups, lanes, complete JSON, models, carrier keys and
source provenance. `--extension-receipt` requires a clean source-bound build.

Each of the 11,520 models has six upstream carrier memberships. DuckVEP emits four
memberships for the missing/retained samples; the two pure-reference memberships
remain implicit by API policy. Their upstream sequences/counts and absence from
native output are checked separately, without padding the native output with oracle
sequences. All expected groups, native rows, lane observations and disagreements are
retained. This is sequence/count/provenance conformance, not a comparison of every
native field to container JSON, arbitrary genetic-code/model coverage, general
ploidy inference or whole-haplotype consequences.

Add `--noncoding-contributors` to `haplotype_sql_differential.R` to run a separately
receipted augmented corpus after the original gate passes. It adds one homozygous
deep-intronic allele per transcript, reruns the pinned executable Haplosaurus, and
requires its complete observation to
remain unchanged. Public replay must retain the added source on all six diploid lanes,
including lanes with no coding edit, while preserving CDS/protein and coding flags.
Five additional corruption controls guard full carrier keys, sequence, provenance and
projection status. This checks literal replay, not splice prediction or combined SO.

`sequence_diff_differential.R --artifacts results/<haplotype-run>
--public-artifacts results/<haplotype-sql-run>` is a separate aligned-difference
lane. It reuses every CDS/protein sequence pair from the receipted six-lane
Haplosaurus cases, including unchanged reference lanes, and checks the native
alignment kernel against the pinned `TranscriptHaplotype::_get_raw_diffs` path.
The additional observer supplies the already validated reference/alternate
sequence operands; it does not replace the original parser/mapper/container
observer. The explicit environment must lack `Bio::Ext::Align`: this contract is
the release-116 pure-Perl NW score and tie order, not its optional alternative.
Public `cds_differences` and `protein_differences` from both policies must match the same complete expected
run lists. Missing/extra pairs and runs, either sequence coordinate, the alignment
coordinate, and either allele are corruption controls. Receipts retain every
pair and mismatch; no sampling or alteration of prior denominators occurs.
This proves alignment and public differences for the supplied standard,
uncurated generated transcripts, not general reference-peptide construction,
compound SO or HGVS. The pure-C suite separately compares the exact band against a full
matrix on exhaustive short pairs and exercises all declared allocation limits.

`reference_translation_differential.R` is a separate reference-peptide
investigation over all 125 ACGTN triplets in start/internal/terminal positions,
all 24 supported tables, all three trailing-partial lengths, and 13 explicit
start, stop, ambiguity, case and peptide-edit witnesses (27,013 cases). Real
Ensembl single-exon models feed `Transcript::translate`, `Translation::seq`
and `TranscriptHaplotypeContainer`; only the input slice's table attribute is
supplied without a database. No translation or container method is overridden.
The lane compares native reference-peptide preparation and complete/stop-truncated
alternate proteins against those methods. It separately retains the rejected
raw-translation-as-reference hypothesis and its mismatches; these are not
substituted for the actual native-reference comparison. The lane
retains all pairs, reference and alternate failure categories, exact module
hashes, package lock, source identity and missing/extra/sequence corruption
controls, and fails if any native reference or alternate comparison differs.
Original Haplosaurus inputs and expectations are not changed. This
lane does not certify public protein differences, compound SO or HGVS.

`compound_coding_audit.R --artifacts results/<haplotype-run>` observes the native
coding-context evaluator on every original edit set, including reference lanes.
It verifies complete CDS/displayed-protein replay and rejects false supported
consequences for compound indels, including net-zero sets. This is a support-limit
audit, not a combined-SO differential; unsupported rows remain explicit and are
never counted as biological agreement. The original corpus and assertions are unchanged.
It also opens every actual interaction block through the shared coding-window
API and checks its residues against the complete translated context. Separate
reference/alternate peptide offsets preserve earlier in-frame shifts; the local
length request uses the block's change, not the complete path's. Failed windows,
shifted-axis coverage and residue differences are retained as separate counters.
`blocks.csv` additionally retains every block's local predicate/support result,
physical span and earlier-stop observation. Substitution blocks reuse the
independent coding interpreter on their actual reference/alternate coordinates;
indel blocks reuse the independent length-change and start/terminal interpreters with
actual frame/stop geometry. Selected rebuilt block spans borrow unchanged reference
flanks; unrelated blocks do not become part of that local predicate's operands.
Single-record genomic insertion-length reach remains separate. The support audit
records supported endpoint blocks separately and rejects out-of-CDS
support, partial facts on failure and contradictory frame/in-frame facts.
The whole-context compound-indel substitution guard remains unchanged.
A later local missense predicate is not
evidence that translation reached that block. Neither local predicates nor
unsupported results are collapsed into a whole-haplotype consequence set.
Each block also records the full translation's `first_stop_position1` (zero when
absent), `frame_status`, and `stop_overlaps_displacement`. The last fact intersects
that stop's three alternate CDS bases with actual frame-changing/restoring edit
spans. A stop after restoration is distinct from one inside displaced bases;
these observations do not classify SO terms or establish protein rescue.
The public `haplotype_sql_differential.R` lane separately compares
`stop_in_displaced_frame` to rebuilt per-base coordinate markers for every occupied
carrier in both policies. Its original sequence/oracle assertions and denominators
are unchanged; frame comparisons, failures and a flipped-fact rejection control
are additional metrics, not VEP compound-consequence agreement.

The additional `bcftools csq -p a` observation builds bcftools/HTSlib 1.23 from
the exact `src/bcftools-1.23` tree in RBCFTools commit
`9adeaf4cfcc3bff40efca6237749fefb53391678`, exported from the local mirror into
the result directory. It does not reuse a prebuilt binary, modify the mirror,
disable assertions, or fetch dependencies. Receipts bind source archive, built
executable, command, complete native lane observations and upstream output/logs.
An upstream crash remains a failing audit with partial rows and missing carriers
counted; it does not interrupt observation of the complete native denominator.
bcftools supplies complementary compound-state evidence, not the VEP-116 oracle
or a replacement production classifier. In particular, seed 173's DHT000002
restores its DNA frame yet translates to `MGLS*`; downstream contributors remain
required. The native, SQL and R fixed witnesses preserve this case independently
of the optional external audit.
The assertion-enabled source build currently aborts at `csq.c:2433` on seed
173's DHT000895. Its four original records succeed with each sample separately
and with `cis,shared`, but fail with `cis,trans`. Retain that cohort-dependent
failure; an assertion-disabled local binary is not evidence that the source-built
lane passed.

`compound_sample_audit.R --audit-artifacts results/<compound-coding-run>` is an
additional sample-separability investigation. The original source-built binary
runs each sample separately on the unchanged complete VCF/GFF/FASTA. A separate
build applies `patches/bcftools_csq_prefix_local.patch`: this proposed correction
to Petr Danecek's bcftools `hap_add_csq` keeps leaf-specific rendering state local
instead of mutating a prefix node shared by another leaf. The original source,
binary, cohort failure and VEP oracle are untouched; the patched build is not
relabeled as upstream conformance. The audit checks the full occupied-carrier
domain and full row multisets, retains missing/extra rows and duplicate counts,
and fails on any discrepancy. On the existing two seeds the patch removes the
assertion and preserves every distinct consequence row, but cohort emission
still has extra duplicate rows. Those duplicates are unresolved evidence, not
discarded rows or a passing conformance result.

This is **not** a public phased-executor certificate: the R harness materializes
decoded calls, and it does not test native DuckDB carrier streaming, strict phase/PS
interpretation, compound SO/HGVS, structural
composition or arbitrary ploidy. Haplosaurus exposes sequence differences and frame
flags, not a compound SO/HGVS oracle. Its offline container also defaults to two lanes
without inferring VCF ploidy; that behavior must not silently define DuckVEP's ploidy
contract. The existing independent-event and pure-C property campaigns remain separate
and unchanged.
