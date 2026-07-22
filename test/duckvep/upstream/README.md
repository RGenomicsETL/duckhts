# Upstream Ensembl test authorities

Status: current source-lineage receipt. These files are exact, checksum-pinned mirrors of
public Ensembl tests used to identify semantic witnesses and to validate the executable
VEP oracle environment. Copying a Perl test here does **not** mean that DuckVEP directly
passes that Perl test.

## Why the sources are mirrored

DuckVEP is a compatibility rewrite. A fixed DuckVEP witness that originates in an Ensembl
test must remain reviewably connected to the exact upstream source, release, and expected
result. The repository therefore keeps the relevant public test authorities under paths
that preserve the upstream repository, release, and source path:

```text
ensembl-vep/release-116/t/...
ensembl-variation/release-116/modules/t/...
ensembl-tools/release-89/scripts/variant_effect_predictor/t/...
```

[`sources.tsv`](sources.tsv) is the machine-readable receipt. For every mirrored `.t`
file it records the upstream repository URL, ref, exact commit, upstream path, local path,
SHA-256, and one of these roles:

- `oracle_self_test`: exercises the pinned VEP executable, parsers, cache fixtures,
  regulatory/SV plumbing, output factories, or limited Haplosaurus behavior;
- `semantic_fixture_source`: contains targeted consequence or HGVS expectations from
  which a DuckVEP-native fixture may be extracted with source file, test name/line,
  assembly/model, input, and expected result retained;
- `historical_lineage`: preserves the smaller release-89 monolithic VEP suite so it is
  not mistaken for a separate, larger public conformance authority.

Run the offline integrity check with:

```bash
perl test/duckvep/upstream/check_sources.pl
```

That offline check rejects a changed byte relative to `sources.tsv`, an incomplete 49/2/5
authority inventory, an unreceipted `.t` file, or a role inconsistent with its authority.
For release-source verification, compare both bytes and complete public-suite inventories
against local Git objects at the declared commits:

```bash
DUCKVEP_ENSEMBL_VEP_GIT=/path/to/ensembl-vep \
DUCKVEP_ENSEMBL_VARIATION_GIT=/path/to/ensembl-variation \
DUCKVEP_ENSEMBL_TOOLS_GIT=/path/to/ensembl-tools \
  make -f Makefile duckvep-upstream-git-check
```

The Git-backed target requires local exact-commit mirrors and is therefore not part of the
network-free ordinary build. It is mandatory in the release evidence audit:

```bash
DUCKVEP_ENSEMBL_VEP_GIT=/path/to/ensembl-vep \
DUCKVEP_ENSEMBL_VARIATION_GIT=/path/to/ensembl-variation \
DUCKVEP_ENSEMBL_TOOLS_GIT=/path/to/ensembl-tools \
  make -f Makefile duckvep-release-conformance-audit
```

The Git-backed check prevents a coordinated local file/manifest rewrite from pretending
to be the pinned upstream commit. It also proves that all `.t` files in the VEP-116 and
legacy-tools suite directories are present; the Variation authority is deliberately the
two named semantic files rather than its entire module test tree. Updating a mirror is
an explicit compatibility-source change: copy from the new exact upstream commit,
regenerate `sources.tsv`, inspect the diff, and rerun the extracted DuckVEP witnesses and
executable differentials affected by the source change.

## Pinned authorities

| Authority | Ref | Commit | Use |
| --- | --- | --- | --- |
| [Ensembl VEP](https://github.com/Ensembl/ensembl-vep/tree/release/116/t) | `release/116` | `57ea5c52340acc1f156267f810ad162e26597082` | 49 executable upstream test files; oracle health and source cases for parser, cache, regulatory, SV, output, and limited Haplosaurus behavior. |
| [Ensembl Variation](https://github.com/Ensembl/ensembl-variation/tree/release/116/modules/t) | `release/116` | `2fb834b987ede3824e200197a838ce11e91aeb4b` | `variation_effect.t` (189 targeted consequence cases) and `hgvs_parser.t` (genomic, transcript, and protein HGVS cases); targeted semantic fixture sources. |
| [Ensembl Tools](https://github.com/Ensembl/ensembl-tools/tree/release/89/scripts/variant_effect_predictor/t) | `release/89` | `9d50991e514ea414981d3b5513256330daadce78` | Five historical monolithic VEP tests, retained only as lineage. |

The full upstream suites depend on Perl modules, databases, caches, and fixture files that
are intentionally not vendored wholesale here. The exact VEP release-116 source commit was
run with the Linux-64 Conda environment in
[`vep116_2026-07-22.conda-explicit.txt`](receipts/vep116_2026-07-22.conda-explicit.txt),
linked from [`self_test_receipts.tsv`](self_test_receipts.tsv). The lock records the complete
solved package environment, including VEP 116.0, Perl 5.32.1, and BioPerl 1.7.8; Git already
identifies this small manifest, so it is not wrapped in a redundant second digest. That run
passed all 49 files and 1,977
assertions; 293 assertions were explicitly skipped (281 without a local database, 11
without a local database or `Set::IntervalTree`, and one without `HTML::Lint`). Optional
Perl modules can change upstream's dynamic plan, so the receipt records the measured
denominator rather than treating an environment-dependent 1,979-assertion count as a
constant. The receipt names and hashes the committed verbose TAP output. The offline
checker requires the VEP authority and exact 49-file pinned inventory, rejects every
`not ok` or TAP bailout, validates per-file assertion numbering, plan and terminal state,
and recomputes the assertion, skip-reason, success-summary, and result counters.
This establishes oracle health only. DuckVEP conformance remains established by
DuckVEP-native fixed/property tests and by differentials that execute both engines on the
same declared model and allele corpus.

The recorded self-test was reproduced from a clean `t/` tree with:

```bash
VEP_GIT=/path/to/ensembl-vep-at-57ea5c52340acc1f156267f810ad162e26597082
VEP_ENV=/path/to/new/vep116-environment
micromamba create -y -p "$VEP_ENV" \
  --file test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt
(
  cd "$VEP_GIT"
  PERL5LIB="$PWD/modules:$VEP_ENV/share/ensembl-vep-116.0-0" \
    "$VEP_ENV/bin/prove" -lv t/*.t
)
```

The mirrored files are the upstream cases, byte for byte. A DuckVEP-native translation is
a separate test and must cite the repository, commit, source path, and upstream test
description or line in the SQL/C/R test itself. Until that translation exists, the source
is labelled `semantic_fixture_source`; its presence here is not counted as a passing
DuckVEP assertion.

## Licensing

The mirrored Ensembl test sources are distributed under the Apache License 2.0 and retain
their upstream copyright/license headers where present. The included
[`LICENSE-APACHE-2.0.txt`](LICENSE-APACHE-2.0.txt) is the upstream license text. DuckHTS's
own source remains under its repository license; the upstream mirrors are not relicensed.
