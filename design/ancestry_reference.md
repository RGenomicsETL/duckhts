Status: current implementation guidance; public reference staging remains open.

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

Figshare files require independently verified checksums before registration
in the benchmark registry. A redirected or merged copy is not an acceptable
identity for either original download. The published epilepsy and 1000 Genomes
comparisons require the paired Figshare files and matching sample inputs.
