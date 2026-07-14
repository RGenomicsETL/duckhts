# DuckVEP compatibility errata

Status: current implementation guidance. This ledger records non-obvious VEP 116
behaviour that DuckVEP intentionally reproduces. Code and executable differential tests
remain authoritative; this is not an implementation diary.

VEP compatibility describes an observed result, not an endorsement of the underlying
biology or API design. When VEP's predicate ordering produces an unusual combination of
terms, DuckVEP must reproduce that state before offering a separately named alternative.

## An uploaded span can suppress a simpler coding edit

VEP 116 maps an equal-length uploaded `VariationFeature` to CDS as one span; it does not
default-minimize the substitution. If that complete feature starts in CDS and ends in
the transcript-oriented 3-prime UTR, one endpoint returned by
`BaseTranscriptVariation::cds_coords` is a `Bio::EnsEMBL::Mapper::Gap`. Consequently,
`TranscriptVariationAllele::codon` and `::peptide` are unavailable and
`VariationEffect::coding_unknown` emits `coding_sequence_variant` instead of a more
specific missense or stop consequence.

This is representation-dependent VEP behaviour. Semantic prefix/suffix trimming may
leave a single changed CDS base, but evaluating only that smaller edit produces a result
for a feature VEP was not asked to annotate. DuckVEP therefore retains both geometries:
the uploaded feature decides whether VEP can map peptide state, while the minimized edit
is used only after that mapping state permits sequence evaluation.

The paired VEP-116 witnesses isolate the distinction in the terminal `TAA` codon:

| Uploaded feature | Mapping state | Exact VEP terms |
| --- | --- | --- |
| `chrDuck:239 AA>CA` | both endpoints in CDS | `stop_lost` |
| `chrDuck:240 AA>TA` | final endpoint in 3-prime UTR | `3_prime_UTR_variant&coding_sequence_variant` |

Do not generalize this suppression to a feature crossing the 5-prime UTR boundary.
VEP has separate start-codon predicates that may still emit `start_lost` or
`start_retained_variant` without the ordinary peptide path. The executable witnesses and
held-out differentials are the authority for the precise boundary.

## Retained feature bases participate in peptide predicates

When an equal-length uploaded feature maps wholly into CDS, VEP 116 uses every codon
covered by that feature as the local peptide window. Common-prefix and suffix trimming
still identifies the bases to apply, but it does not shrink the window used by the start,
stop, synonymous, and missense predicates. Two records with the same changed CDS base can
therefore receive different consequences.

Paired VEP-116 witnesses isolate three forms of this representation dependence:

| Complete uploaded feature | One-base counterpart | Exact VEP distinction |
| --- | --- | --- |
| `chrDuck:122 GG>GA` | `chrDuck:123 G>A` | `start_lost&start_retained_variant` vs `missense_variant` |
| `chrDuck:237 GT>TT` | `chrDuck:237 G>T` | `stop_retained_variant` vs `missense_variant` |
| `chrDuck:237 GT>AT` | `chrDuck:237 G>A` | `missense_variant` vs `stop_gained` |

The first pair is especially non-obvious: the complete local peptides differ, so
`start_lost` is true, while the rebuilt first codon is still `ATG`, so the independently
evaluated `start_retained_variant` predicate is also true. In the third pair, the complete
window has a stop in both peptides at different local indexes, so none of VEP's
stop-lost, stop-retained, or stop-gained predicates succeeds and `missense_variant`
remains.

Keep one substitution evaluator over an explicitly selected peptide window. The uploaded
feature selects that window; the trimmed edit supplies the replacement bases. Do not
normalize the uploaded identity inside the consequence kernel or implement separate
start/stop rule copies for MNV-shaped input.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::start_lost`,
`::start_retained`, `::stop_lost`, `::stop_retained`, `::stop_gained`, and
`::missense`. Fixed witnesses must remain paired with large held-out differentials so a
special-case term substitution cannot satisfy conformance.

## ALT-only mismatch islands inherit an expanded intron cache

VEP 116 can emit `intron_variant` for a lengthening replacement whose REF-shaped
`VariationFeature` remains entirely exonic. The result comes from two stages that are
easy to miss when reading one predicate in isolation:

- the VCF parser removes one shared indel anchor but leaves the feature end based on REF;
- with `Set::IntervalTree`, `_overlapped_introns` preselects and caches introns over a
  three-base flank on both exon sides;
- `_intron_effects` later compares REF and ALT bytewise, and an ALT-only mismatch island
  may extend beyond the cached feature span into the intronic interior.

The intronic predicate still excludes the essential first and last two intron bases.
Consequently, an island can produce `intron_variant` together with coding, donor, and
donor-fifth-base terms. Moving the REF-shaped feature one base beyond the three-base
cache flank suppresses all intron-derived terms even when the ALT-only island reaches the
same bases. This is cache-dependent VEP behaviour, not a general interval-overlap rule.

The pinned VEP environment uses `Set::IntervalTree` 0.12. Its source expands the cached
intron interval in `BaseTranscriptVariation.pm::_create_intron_trees`, then consumes the
cached list in `BaseTranscriptVariationAllele.pm::_intron_effects`. Fixed witnesses at
`chrDuck:146 CGT>CACTGAGGGC` and `chrDuck:145 ACG>ACCTTCTGTGTA` pin the two sides of the
cache boundary. Large held-out differentials remain necessary because transcript 3-prime
shifting can move the predicate geometry before this cache is consulted.

## Terminal-stop insertions are predicate states, not modulo-three classes

VEP 116 can classify an insertion whose length is not divisible by three as
`inframe_insertion`, and can combine it with `coding_sequence_variant`, `stop_lost`, or
`stop_retained_variant`. The result follows the order and local state of
`stop_retained`, `stop_lost`, `frameshift`, `inframe_insertion`, and `coding_unknown`; it
cannot be reconstructed from insertion length alone.

The important upstream details are:

- `VariationEffect.pm::frameshift` returns false when the first affected reference
  peptide begins with `*`.
- `VariationEffect.pm::inframe_insertion` trims the alternate peptide after its first
  stop, then accepts a preserved reference prefix or suffix.
- `TranscriptVariationAllele.pm::peptide` appends `X` for a trailing partial codon,
  except when the complete translated peptide is exactly `*`.
- `VariationEffect.pm::coding_unknown` excludes frameshift and stop predicates, but not
  `inframe_insertion`, so both terms can be true.

Executable VEP-116 witnesses in the terminal `TAA` fixture include:

| Insertion site and local alternate codon | Exact VEP terms |
| --- | --- |
| inside `TAA`: `tAGCaa` (`*Q`) | `inframe_insertion&stop_retained_variant` |
| inside `TAA`: `tAGGTaa` (`*VX`) | `coding_sequence_variant&inframe_insertion` |
| inside `TAA`: `tCGATGTTATGAaa` (`SML*X`) | `inframe_insertion&stop_lost` |
| immediately before `TAA`: `TAAA` (`*`) | `inframe_insertion&stop_retained_variant` |

Source authority: Ensembl VEP 116, commit
`57ea5c52340acc1f156267f810ad162e26597082`. Preserve these cases as fixed witnesses and
also measure them through deterministic random differentials with held-out seeds. A
targeted witness alone is not evidence that the surrounding state space improved.

## Lengthening edits use two different peptide views

VEP 116 does not reduce every coding length increase to a single biological shape. Its
`inframe_insertion` and `protein_altering_variant` predicates inspect the same
codon-local edit through deliberately different peptide views:

- `inframe_insertion` truncates the alternate peptide at its first new stop, then asks
  whether the reference peptide is a prefix or suffix of that truncated value;
- `protein_altering_variant` tests the untruncated alternate peptide and suppresses the
  term when that complete value preserves the reference prefix or suffix.

Consequently, the local change `Y` -> `*QCY` is `stop_gained` alone. The reference `Y`
is absent from the stop-truncated alternate (`*`), so it is not an in-frame insertion;
it is present as the suffix of the untruncated alternate, so it is not a protein-altering
variant either. By contrast, `Y` -> `YRVALM*...` preserves `Y` before the stop and is
both `inframe_insertion` and `stop_gained`.

This is an upstream state-machine property, not a claim that the terms form a clean
biological partition. DuckVEP must reproduce the two predicate inputs and their unusual
combinations. Do not replace them with a single whole-protein diff, a net-length rule,
or one shared "in-frame versus protein-altering" classifier.

Source anchors: Ensembl Variation 116
`VariationEffect.pm::inframe_insertion` and
`VariationEffect.pm::protein_altering_variant`. The fixed VCF witness
`chrDuck:233:A:AACAATGTTA` has local codons `tac/taACAATGTTAc`, local peptides
`Y/*QCY`, and exact VEP consequence `stop_gained`. Keep a paired positive
`inframe_insertion&stop_gained` witness so a blanket suppression cannot satisfy the
test.

## In-frame deletion is not the inverse of in-frame insertion

VEP 116 tests `inframe_deletion` against the codon strings, not the peptide strings used
by `inframe_insertion`. It accepts a shorter alternate codon string when that string is
a reference prefix or suffix, or when common-edge trimming consumes the alternate and
leaves a reference remainder divisible by three. A net CDS change of minus three bases
is therefore insufficient.

This distinction matters after semantic allele trimming. A raw deletion-like delins
whose alternate is an exact nucleotide edge reduces to a pure deletion. A genuinely
non-empty local replacement may instead be `protein_altering_variant`, even though the
CDS becomes shorter by one codon. The fixed kernel scene replacing local peptide `PG`
with `A` is such a case.

Source anchors: Ensembl Variation 116
`VariationEffect.pm::inframe_deletion`,
`VariationEffect.pm::protein_altering_variant`, and
`Sequence.pm::trim_sequences`. Keep the deletion and insertion predicates separate;
making one the length-sign mirror of the other changes VEP results.
