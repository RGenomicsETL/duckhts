# bcftools Conformance Harness

This harness standardizes DuckHTS vs upstream `bcftools` plugin comparisons for:

- `bcftools_score`
- `bcftools_munge_row` / `duckdb_munge`
- `bcftools_liftover` / `duckdb_liftover`

## Usage

```bash
bash scripts/run_bcftools_conformance.sh all
bash scripts/run_bcftools_conformance.sh score
bash scripts/run_bcftools_conformance.sh munge
bash scripts/run_bcftools_conformance.sh liftover
```

Outputs go under `.tmp/conformance/<scenario>/` by default and include:

- `raw/`: direct engine outputs and stderr captures
- `normalized/`: canonicalized tables used for comparison
- `report/summary.tsv`: case-level status and reason buckets
- `report/*.details.tsv`: row-level mismatches for each case

## Scenario Contract

Each scenario script in `scripts/conformance/scenarios/` must:

1. accept one argument: the scenario output directory
2. write `report/summary.tsv`
3. optionally write one or more `report/*.details.tsv`

The runner collates all scenario summaries into one top-level `summary.tsv`.

## Reason Buckets

The harness distinguishes raw mismatch status from the likely cause:

- `match`: canonical outputs agree
- `float_precision`: only small numeric drift remains
- `preset_alias_gap`: upstream found a valid preset alias path that DuckHTS did not
- `error_handling_divergence`: one engine emits structured rows where the other hard-errors or rejects the record differently
- `surface_difference`: both are behaving intentionally differently at the API/output-surface level
- `unclassified`: a real diff that still needs inspection

Reasons are scenario-specific on purpose. The framework standardizes artifacts and reporting, but it does not pretend every semantic difference can be inferred from a generic numeric diff alone.
