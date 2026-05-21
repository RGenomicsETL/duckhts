# DuckVEP keying notes: VariantKey, RegionKey, and overlap kernels

## Purpose

This note records how DuckHTS's new VariantKey / RegionKey support should fit a future DuckVEP-style annotation layer.

The core split is:

- use `VariantKey` for exact small-variant identity joins
- use `RegionKey` and/or interval kernels for span/overlap annotation
- keep symbolic / SV / BND annotation as overlap-first plus exact metadata refinement

## 1. Exact variant annotation

For exact allele-based annotation stores such as dbSNP-, ClinVar-, or cohort-style joins:

- `variantkey(chrom, pos, ref, alt)` is the primary fast key
- reversible VariantKeys are exact for normalized biallelic variants
- nonreversible / hashed VariantKeys are still useful as a fast prefilter

Important consequence:

- hashed VariantKey rows must keep the original `REF` / `ALT`
- exact joins on hashed rows must verify `REF` and `ALT` after the key match
- DuckHTS does **not** need an NRVK reverse-lookup sidecar for SQL joins as long as the annotation relation still stores the original alleles

Recommended exact-match pattern:

1. join on `variantkey`
2. treat reversible rows as exact
3. for nonreversible rows, refine on `chrom`, `pos`, `ref`, and `alt`

This mirrors the intended DuckDB use case better than a decode-only design.

## 2. Region / overlap annotation

For interval-oriented annotation sources such as:

- genes
- transcripts
- exons / CDS / UTRs
- regulatory BED/GFF/GTF tracks
- OMIM-style locus intervals
- other span-based labels

use a span representation rather than VariantKey.

`regionkey(chrom, start, end[, strand])` is useful as a compact stored region identifier and for cheap same-span / coarse overlap-oriented work. However, for real overlap execution the primary engine should still be an interval kernel such as:

- `duckhts_cgranges_*`
- future DuckDB-native overlap primitives

So the practical layering is:

- `RegionKey` for compact persisted span keys and coarse filtering
- cgranges / interval kernels for actual overlap enumeration and candidate generation

## 3. Symbolic alleles, SVs, and breakends

VariantKey only encodes:

- `CHROM`
- `POS`
- `REF`
- `ALT`

It does **not** encode:

- `END`
- `SVLEN`
- mate breakend coordinates
- breakend orientation
- other SV metadata

RegionKey encodes a span, but a span alone is still not exact SV identity.

Therefore symbolic / SV / BND annotation should be treated as a separate path:

### 3a. Overlap-oriented SV annotation

Use span keys and interval kernels over:

- `chrom`
- `start`
- `end`
- optional strand/orientation summaries when meaningful

This is the correct shape for gene / regulatory / locus overlap annotation.

### 3b. Exact or near-exact SV matching

Do not rely on VariantKey or RegionKey alone.

Retain and refine on the relevant fields, for example:

- `svtype`
- `end`
- `svlen`
- `ref`
- `alt`
- mate chrom / mate pos
- breakend orientation fields

## 4. Suggested future DuckVEP storage model

A future annotation cache can legitimately store both:

- `vk` for exact allele joins
- `rk` for region/span joins

plus the original fields needed for exact refinement.

A practical row shape for variant-centric annotation tables is:

- `chrom`
- `pos`
- `ref`
- `alt`
- `vk`
- `vk_refalt_code`
- `vk_reversible`
- optional `start0`, `end0`, `rk` for span-aware annotations
- annotation payload columns

A practical row shape for region-centric annotation tables is:

- `chrom`
- `start0`
- `end0`
- optional `strand`
- `rk`
- feature / source metadata
- payload columns

## 5. Planned consequence-engine layering

A serious DuckVEP-style system should separate three steps:

1. **Exact variant annotation**
   - VariantKey-driven joins for small variants
   - hashed-key exact refinement where needed

2. **Region overlap annotation**
   - RegionKey / cgranges candidate generation for transcripts, exons, genes, regulatory regions, and locus tracks

3. **Exact consequence logic**
   - transcript, exon/CDS, strand, phase, codon, and allele-aware sequence consequence computation

RegionKey helps with candidate narrowing, but it is not itself a consequence engine.

## 6. Echtvar relation

`echtvar` is the relevant external performance reference for exact variant annotation archives.

For future benchmarking, compare:

- DuckHTS exact-key SQL joins on normalized biallelic variants
- DuckHTS mixed reversible/nonreversible fallback joins
- `echtvar` annotation throughput on equivalent exact-match workloads

But note that `echtvar` is an exact-variant archive design, not a transcript- or region-overlap consequence engine. For a full DuckVEP layer, exact-key and overlap/consequence workloads must be benchmarked separately.
