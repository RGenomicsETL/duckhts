# Reproducibility closure checklist

Status: open implementation checklist for PR #178. Remove completed items rather than
turning this file into a historical plan; code, tests, rendered reports, and git history
retain completed work.

## Merge blockers

- [ ] Let CI finish on the current head and fix every remaining platform failure.
- [ ] Request and receive a Codex review of the final head; address every actionable
  finding and repeat review when a fix changes the head.
- [ ] Implement `make stage-variantkey-providers` as a registry-backed executable
  stage command, replacing the temporary user-supplied-model workaround.
- [ ] Implement the public Ensembl 116 core/funcgen plus FASTA model producer, and
  derive the model regulatory/transcript Parquet outputs.
- [ ] Move DuckVEP HPRC, Sniffles, and dbVar corpus staging into `duckhtsbench`.

## Registry closure

- [ ] Add a completeness test: each active benchmark or conformance external input
  has a registry ID, one cache path, source/release/access/transform/consumer metadata,
  and an executable stage route.
- [ ] Audit browser, wasm, and webR assets against the same registry rather than only
  the shared cache helper convention.
- [ ] Verify each source with a publisher-provided digest or manifest when available;
  retain a known size guard only when it is the available verified identity.
- [ ] Add a network-free fake-source staging test for every migrated stage family.
- [ ] After each registry/stager migration, run the smallest representative real
  benchmark or conformance workload; retain its rendered report or say precisely why
  the declared public input cannot yet be staged.

## Documentation and cleanup

- [ ] Correct the stale root `NEWS.md` statement that VariantKey/FastVEP reports were
  retired; they were restored.
- [ ] Update the PR description to distinguish implemented registry stages from
  blocked model/corpus work.
- [ ] Keep the versioned DuckVEP model contract in `design/duckvep.md` aligned with
  the executable producer when it exists.
- [ ] Remove the temporary untracked VariantKey stager/test files after their
  registry-backed replacement is committed. Do not touch the unrelated HPRC
  investigation.

## Final verification

- [ ] `make test_release`
- [ ] `make test-benchmark-registry`
- [ ] Rduckhts bootstrap and tinytest
- [ ] `make docs`
- [ ] Relevant staged real workload and rendered benchmark checks
- [ ] `git diff --check`
- [ ] All CI green
- [ ] Current-head Codex approval
