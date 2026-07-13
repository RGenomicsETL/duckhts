# DuckVEP compatibility errata

Status: current implementation guidance. This ledger records non-obvious VEP 116
behaviour that DuckVEP intentionally reproduces. Code and executable differential tests
remain authoritative; this is not an implementation diary.

VEP compatibility describes an observed result, not an endorsement of the underlying
biology or API design. When VEP's predicate ordering produces an unusual combination of
terms, DuckVEP must reproduce that state before offering a separately named alternative.

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
