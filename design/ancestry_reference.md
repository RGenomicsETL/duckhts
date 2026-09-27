Status: current implementation guidance.

# Ancestry reference contract

The Figshare products from the bigsnpr ancestry vignette are the reference
frequencies at https://figshare.com/ndownloader/files/31620968 and the
16-component PC loadings at https://figshare.com/ndownloader/files/31620953.
They are paired products: chromosome, position, and reference/alternate allele
orientation must agree. Reference frequencies are one row per variant and group;
loadings are one row per variant and PC. Select biallelic, unambiguous SNVs with
complete group and PC coverage for a BAM/CRAM panel. Record the genome build,
reference release, and `duckhts_somalier_panel_sha256()` for every materialized
panel. `rduckhts_ancestry_panel()` orders selected loci and assigns zero-based
site indices; its default spacing is 5,000 bases with at most 17,000 sites.

`test/data/ancestry_correction_bigsnpr_1.12.21.tsv` is the 16-component
projection-shrinkage correction vector printed in the vignette at
https://privefl.github.io/bigsnpr/articles/ancestry.html. It is paired only
with those loadings. The duckhtsbench registry checks its committed bytes
before staging. The vector must not be applied to different PCA loadings.

The duckhtsbench registry pins the compressed Figshare files by exact byte
count and SHA-256. Automated fetches from some hosts receive HTTP 403 or a WAF
browser challenge. Place browser-downloaded files at their registry cache paths
and run `duckhts_bench_fetch(id)` to verify them without contacting Figshare.
A redirected or merged copy is not an acceptable identity. The epilepsy
input and public GRCh37 phase-3 1000 Genomes chr22 genotypes have independent
registry identities. On chr22, their allele agreement with
the reference products must be checked before site selection. GRCh38 30x CRAMs
require the registered GRCh37-to-GRCh38 chain and source/destination FASTAs;
retain the mapped, rejected, swapped, and duplicate-destination denominators
before comparing aligned reads to phased genotypes. The ancestry comparison
must retain its per-sample matching counts and correlation-gate failures.
