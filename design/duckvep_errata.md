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
