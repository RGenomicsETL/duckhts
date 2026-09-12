## Submission

Rduckhts 1.5.2-0.1.5 packages DuckHTS 1.5.2 for the coordinated R and DuckDB
community-extension release. It adds selected typed FORMAT fields and original
genotype text to the genotype reader, improves reader allocation-failure
handling, and fixes independent-event consequence and protein-HGVS behavior
against pinned Ensembl VEP 116. The alpha haplotype interface supplies phased
sequence replay and provenance; complete compound consequence prediction is
outside this release's supported scope. Details are in `NEWS.md`.

This release also addresses the upcoming `duckdb` R release's reverse-dependency
failure. The connection test checks rejection and preservation of the existing
driver, connection and unsigned-extension policy without depending on DuckDB's
diagnostic wording. The original test reproduces the reported failure against
`duckdb` R commit `6a05aded03c52e1ef11cfa6804f146f86017502d`; the revised test
passes against that same build.

## Test environments

Local source-tarball installation and `R CMD check --as-cran`: Ubuntu 24.04.3
LTS, x86_64, R 4.6.0, GCC 13.3.0, and `duckdb` R 1.5.3.

Connection compatibility was also tested with `duckdb` R 1.5.5 and the pinned
release-candidate source above, which reports R package version 1.5.5.9013.36
and DuckDB engine v1.5.6-dev150 (`a3cd0deed1`).

## R CMD check results

0 errors, 0 warnings, 2 notes:

- CRAN incoming feasibility reports seven package updates in the preceding six
  months. This coordinated feature release also prepares for the upcoming
  `duckdb` R release.
- The local R installation supplies `-mno-omit-leaf-frame-pointer`, which the
  compilation-flags check reports as non-portable. Rduckhts does not add this
  flag; it preserves the installation's compiler settings.
