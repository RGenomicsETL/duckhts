# VEP 115 -> 116 behavioral change catalog: consequence predicates

Status: **reference catalog (2026-07-10).** Per-hunk behavioral catalog of the VEP 115→116
consequence-predicate delta (`VariationEffect.pm` / `Sequence.pm` / `TranscriptVariationAllele.pm`),
grounding the "pin VEP 116" decision in `duckvep_m6a_contract_gate.md` and telling M6b
implementers exactly which predicates moved. 17 meaningful hunks, 13 on the coding frontier;
the "Implications for pinning 116" section is the actionable output. Two hunks were escalated
to gpt-5.6-sol and resolved — see "Escalated hunks — RESOLVED" (dead `$consider_ins_len` param;
`c./n.` uses a separate unchanged 3′-shift routine from `g./m.`).

Sources diffed (byte diff, not `git diff`, since 115/116 live in unrelated trees):

- `VariationEffect.pm`
  - 115: `/root/ensembl-variation/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`
  - 116: `/root/miniconda3/envs/vep/share/ensembl-vep-116.0-0/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`
- `Sequence.pm`
  - 115: `/root/ensembl-variation/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm`
  - 116: `/root/miniconda3/envs/vep/share/ensembl-vep-116.0-0/Bio/EnsEMBL/Variation/Utils/Sequence.pm`
- `TranscriptVariationAllele.pm`
  - 115: `/root/ensembl-variation/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`
  - 116: `/root/miniconda3/envs/vep/share/ensembl-vep-116.0-0/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm` (actual 116
    path is directly under `Bio/EnsEMBL/Variation/...`, **not** under a `modules/` prefix as the task brief guessed —
    `modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm` does not exist in the 116 tree; verified with `find`.)

`diff -u` hunk counts: VariationEffect.pm = 14, Sequence.pm = 2, TranscriptVariationAllele.pm = 4. One hunk in each
of VariationEffect.pm and Sequence.pm and TranscriptVariationAllele.pm is only the `Copyright [2016-2025]` ->
`Copyright [2016-2026]` license-header bump; those are excluded from the behavioral tables below (noted once here
and not repeated as rows). That leaves **13 meaningful hunks in VariationEffect.pm, 1 in Sequence.pm, 3 in
TranscriptVariationAllele.pm** = **17 meaningful hunks total**.

---

## Table 1: VariationEffect.pm

| # | predicate/function | 115 behavior (quoted) | 116 behavior (quoted) | semantic change | SO term(s) affected | M6b frontier? |
|---|---|---|---|---|---|---|
| 1 | `inframe_insertion` | `VariationEffect.pm:1121-1126` (115): no guard between the stop-retained/start-retained checks and the stop-trim step; goes straight to `$alt_pep =~ s/\*.+/\*/;` then `return 1 if ($alt_pep =~ /^\Q$ref_pep\E/) \|\| ($alt_pep =~ /\Q$ref_pep\E$/);` | `VariationEffect.pm:1121` (116): adds `return 0 if $ref_pep eq "*" && $alt_pep eq "*"; # e.g. ref codon - TAG, alt codon - TAAG` immediately before the stop-trim/match step. | New early-exit: when the reference AA at the affected codon is already a stop (`*`) **and** the raw (untrimmed) alt peptide is also exactly `*` (i.e. the alt codon is itself a longer variant of a stop codon, e.g. `TAG`->`TAAG`), 116 refuses to call this `inframe_insertion`. In 115 this case would fall through to the trim-and-match step and, because `$alt_pep` (post `s/\*.+/\*/`) starts with `*` and `$ref_pep` is `*`, could still match `^\Q$ref_pep\E` and return 1 — i.e. 115 could label a pure "stop codon variant, both sides still `*`" event as `inframe_insertion`, which is really `stop_retained` territory. 116 explicitly excludes it here so it can be classified by `stop_retained`/`ref_eq_alt_sequence` instead. | `inframe_insertion` (and, by exclusion, `stop_retained`) | **Yes** — direct stop-codon / inframe-insertion boundary predicate. |
| 2 | `inframe_deletion` | `VariationEffect.pm:1172-1175` (115): after the length check, goes straight to `# simple string match` / `return 1 if ($ref_codon =~ /^\Q$alt_codon\E/) \|\| ...` | `VariationEffect.pm:1175` (116): adds `return 0 if $ref_pep eq "*"; # e.g. ref codon - TAG, alt codon - G` right after the length check, before the string-match return. | New guard: if the reference peptide at the affected codon is a stop (`*`), 116 never calls the deletion `inframe_deletion` — even if the raw codon strings would otherwise string-match (e.g. ref codon `TAG` / alt codon `G` would previously string-match `$ref_codon =~ /\Q$alt_codon\E$/` and return 1). This routes stop-codon-overlapping deletions to `stop_lost`/`stop_retained` logic instead of `inframe_deletion`. | `inframe_deletion` (and, by exclusion, `stop_lost`/`stop_retained`) | **Yes** — direct stop-codon / inframe-deletion boundary predicate. |
| 3 | `stop_gained` | `VariationEffect.pm:1217-1218` (115): `## check for inframe insertion before stop` / `return 0 if stop_retained(@_);` — no `stop_lost` guard. | `VariationEffect.pm:1221-1222` (116): `return 0 if stop_retained(@_);` / `return 0 if stop_lost(@_);` — comment dropped, new guard added. | 116 makes `stop_gained` and `stop_lost` **mutually exclusive by construction** (stop_lost now suppresses stop_gained). In 115 both predicates were computed independently from `($alt_pep =~ /\*/) and ($ref_pep !~ /\*/)` vs. `($alt_pep !~ /\*/) and ($ref_pep =~ /\*/)`, which are already largely disjoint on their own truth tables, but complex multi-codon-affecting alleles going through the `_ins_del_stop_altered` fallback path could in principle satisfy both independently-cached predicates. 116 forecloses that by explicit ordering: `stop_lost` is now evaluated (and takes priority) before `stop_gained` can fire. | `stop_gained`, `stop_lost` | **Yes** — call-ordering change at the core of the stop/inframe frontier. |
| 4 | `stop_lost` | `VariationEffect.pm:1236-1237` (115): `unless(exists($cache->{stop_lost})) { $cache->{stop_lost} = 0; ...` — no `partial_codon` check anywhere in this sub. | `VariationEffect.pm:1240` (116): `return 0 if partial_codon(@_);` inserted **before** the `unless(exists($cache->{stop_lost}))` cache block. | New guard, evaluated unconditionally on every call (not itself cached — it relies on `partial_codon`'s own cache): if the variant lands in an incomplete trailing codon (< 3 remaining CDS bases), `stop_lost` is now forced to 0 regardless of what the peptide-string comparison would have said. In 115, a partial/incomplete final codon could still be evaluated by the normal `($alt_pep !~ /\*/) and ($ref_pep =~ /\*/)` peptide comparison (where `_get_peptide_alleles` may return `X` for incomplete codons) and, before hunk #5 below, potentially return true. | `stop_lost` | **Yes** — directly named in the task's frontier list ("partial/incomplete codon"). |
| 5 | `stop_lost` | `VariationEffect.pm:1252-1254` (115, inside the TranscriptVariationAllele branch): `if(defined($ref_pep) && defined($alt_pep)) { $cache->{stop_lost} = ( ($alt_pep !~ /\*/) and ($ref_pep =~ /\*/) ); }` | `VariationEffect.pm:1259` (116): `if(defined($ref_pep) && defined($alt_pep) && $alt_pep !~ 'X') { ... same body ... }` | Adds `$alt_pep !~ 'X'` to the guard. `X` denotes an incomplete/ambiguous translated residue (see `ref_eq_alt_sequence`'s own comment "to account for incomplete coding terminal", unchanged in both versions at `VariationEffect.pm:1333`/`1337`). In 115, if `alt_pep` contained `X` the direct-comparison branch was still taken (and would report `stop_lost=0` since `X` doesn't match `/\*/` unless the alt also had a literal `*`, but this materially changes *which* codepath computes the answer). In 116, any `alt_pep` containing `X` is diverted to the `else` branch, i.e. `$cache->{stop_lost} = _ins_del_stop_altered(@_);` — a fundamentally different (indel/overlap-based) computation instead of the string-based one. | `stop_lost` | **Yes** — partial/incomplete-codon-adjacent stop_lost logic. |
| 6 | `stop_retained` | `VariationEffect.pm:1307-1313` (115): `if(defined($alt_pep) && $alt_pep ne '') { ## handle inframe insertion... $cache->{stop_retained} = ref_eq_alt_sequence(@_); } else { $cache->{stop_retained} = (...) && _overlaps_stop_codon(@_) && !_ins_del_stop_altered(@_); }` | `VariationEffect.pm:1313-1318` (116): `if(defined($alt_pep) && $alt_pep ne '' && $alt_pep !~ 'X') { $cache->{stop_retained} = ref_eq_alt_sequence(@_); } else { $cache->{stop_retained} = (...) && _overlaps_stop_codon_cil(@_) && !_ins_del_stop_altered_cil(@_); }` | Two combined changes: (a) same `X`-exclusion pattern as hunk #5 — alleles whose alt peptide contains `X` now fall into the `else` (overlap-based) branch instead of the string-comparison (`ref_eq_alt_sequence`) branch; (b) the `else` branch's overlap/alteration tests are swapped from the plain cDNA-coordinate-based `_overlaps_stop_codon`/`_ins_del_stop_altered` to the new genomic-coordinate, insertion-length-aware `_overlaps_stop_codon_cil`/`_ins_del_stop_altered_cil` (see hunks #10-#12). This changes which insertions/deletions near the stop codon are classified `stop_retained` for the "no clean peptide comparison available" fallback case. | `stop_retained` | **Yes** — core of the stop-codon predicate rewrite. |
| 7 | `ref_eq_alt_sequence` | `VariationEffect.pm:1327` (115): `my $mut_seq = $ref_seq;` placed immediately after `my $ref_seq = $bvfo->_peptide;`, i.e. before `$tl_start`/`$tl_end`/`_get_peptide_alleles` and before the two early `return 0`/`return 1` statements. | `VariationEffect.pm:1346` (116): the equivalent `my $mut_seq = $ref_seq;` is moved to immediately before the `substr($mut_seq, $tl_start-1, ...) = $alt_pep;` line, i.e. after all the early-return checks. | Pure reordering — nothing between the old and new positions mutates `$ref_seq`, and `$mut_seq` is not read before this point in either version. Net effect: no behavioral difference; only avoids one unnecessary scalar copy on the early-return paths (`$ref_pep eq "X" && $alt_pep eq "X"` short-circuit, and the "translation_start beyond ref_seq" short-circuit). | `stop_retained` (no-op contribution) | No — hygiene/perf only, not a behavior change. |
| 8 | `ref_eq_alt_sequence` | `VariationEffect.pm:1348-1354` (115): `my $final_stop_length = length($final_stop) if defined($final_stop) ne '';` then `return 1 if ( ($ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/) \|\| ($ref_seq eq $mut_substring && defined($final_stop_length) && $final_stop_length < 3) \|\| ( $ref_pep =~ /\*/ && (index($ref_pep, "*") + 1 == index($alt_pep, "*") + 1) ));` (three OR'd conditions) | `VariationEffect.pm:1351-1359` (116): `my $final_stop = substr($mut_seq, length($ref_seq)) if ...;` then `return 1 if ( ($ref_seq eq $mut_substring && defined($final_stop) && $final_stop =~ /^\Q*\E/) \|\| ($ref_pep =~ /\*/ && (index($ref_pep, "*") == index($alt_pep, "*"))) );` (two OR'd conditions) | This is the substantive rewrite of the `stop_retained` final-stop-length heuristic: (a) **condition 1 of 115 is deleted outright** — `($ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/)` ("ref AA equals first char of alt peptide, and alt peptide contains a stop anywhere") is gone in 116, so alleles that only satisfied this clause in 115 no longer return `stop_retained=1` from here (they may now be classified by a different predicate, or not at all). (b) **condition 2's test is replaced**: 115 required the "final stop" tail (the part of the mutated peptide beyond the original ref-peptide length, from an insertion) to have `length(...) < 3` — an arbitrary length threshold with no biological meaning tied to whether that tail actually starts with a stop; 116 instead requires the tail to literally start with a stop character, `$final_stop =~ /^\Q*\E/`, with **no length bound**. This is a materially different, more precise "does the extension terminate in a stop right away" test — it can both (i) accept longer 1‑extra-residue-then-stop tails that 115's `<3` bound happened to allow only by coincidence of length, and (ii) newly reject/accept cases where 115's `<3` characters happened not to start with `*` (115 would have returned 1 purely because the tail was short, even if it wasn't a stop) or where the tail is short with a stop at position>0. Net: 116's condition 2 is a correctness fix (checks for an actual stop) replacing a length-heuristic proxy. (c) **condition 3 is algebraically identical** — `index(...)+1 == index(...)+1` in 115 reduces to `index(...) == index(...)` in 116; this is a textual simplification with no behavioral effect. | `stop_retained`, indirectly `stop_gained`/`inframe_insertion` (whatever an allele falls back to when it no longer matches here) | **Yes** — this is the single largest semantic change in the file and squarely on the stop/inframe/partial-codon frontier named in the task brief. |
| 9 | `_overlaps_stop_codon` | `VariationEffect.pm:1372` (115): trailing whitespace on a blank line between the `cdna_start`/`cdna_end` guard and the `overlap(...)` call. | `VariationEffect.pm:1377` (116): same blank line, whitespace-only. | None — whitespace-only. | none | No. |
| 10 | `_overlaps_stop_codon_cil` (new) | Does not exist in 115. | `VariationEffect.pm:1389-1431` (116, new sub): computes overlap using **genomic** coordinates (`$bvf->seq_region_start`/`seq_region_end`) instead of pre-mapped cDNA coordinates, and when the variant looks like a zero-width insertion point on the given strand (`$v_end < $v_start` and the feature sequence is plain `ACTGN`), it extends `$v_start` by `length($vf_feature_seq)` in the direction implied by `$feat->strand` before testing overlap against the transcript's coding-region end (swapped for `strand == -1`), i.e. `overlap($v_start, $v_end, $t_end_s, $t_end_e)` where `($t_end_s,$t_end_e)` is the last/first 3 genomic bases of the CDS depending on strand. Also re-checks `cds_end_NF` the same way as the plain version. | New predicate, not a modification of an existing one, but it changes what `stop_retained`'s fallback branch (hunk #6) and `_ins_del_stop_altered_cil` (hunk #12) see as "overlaps the stop codon." The key behavioral delta vs. plain `_overlaps_stop_codon`: for an **insertion** represented as a zero-length genomic interval immediately before/after the stop codon, the plain cDNA-coordinate overlap test can fail to detect that the *inserted sequence*, once placed, extends into the stop codon (because the insertion point itself, as mapped to cDNA coords, may sit just outside the last-3-bases window) — `_overlaps_stop_codon_cil` compensates by "considering insertion length" and extending the tested window by the inserted sequence's length. | `stop_retained` (fallback path), transitively `stop_lost`/`stop_gained`/`inframe_insertion`/`inframe_deletion` for whichever predicate the allele falls through to if not `stop_retained`. | **Yes** — directly named in the task brief (`_overlaps_stop_codon_cil`). |
| 11 | `_ins_del_stop_altered` | `VariationEffect.pm:1382-1383` (115): four arguments; body calls `_overlaps_stop_codon(@_)`. | `VariationEffect.pm:1433-1434,1444` (116): adds `$consider_ins_len` and forwards it. | Full-tree resolution confirms a no-op: the dispatcher and only executable caller pass four arguments, while `_overlaps_stop_codon` destructures four and drops the forwarded extra. `stop_lost` behavior is unchanged by this hunk. The live insertion-length-aware behavior is the separate `_cil` pair in hunks 10/12. | `stop_lost` (no behavioral contribution from this hunk) | Do not port the dead parameter; implement and test the `_cil` pair. |
| 12 | `_ins_del_stop_altered_cil` (new) | Does not exist in 115. | `VariationEffect.pm:1488-1539` (116, new sub): identical body to `_ins_del_stop_altered` (same UTR/translateable-seq patching, same final-codon translate-and-compare-to-`*` logic) **except** its overlap gate is `return 0 unless _overlaps_stop_codon_cil(@_);` (line 1499) instead of `_overlaps_stop_codon(@_)`. | Same "consider insertion length" genomic-overlap gate as hunk #10, now used for the "was the stop actually altered by this indel" test. This is the function `stop_retained`'s new `else` branch calls (hunk #6): `!_ins_del_stop_altered_cil(@_)`. So for insertions/deletions overlapping the stop region (per the `_cil` genomic test), `stop_retained` in 116 asks "did the final translated codon change away from a stop" using the same downstream translate-and-compare logic as before, but gated by the more insertion-length-aware overlap test. | `stop_retained` | **Yes** — directly named in the task brief (`_ins_del_stop_altered_cil`). |
| 13 | `frameshift` | `VariationEffect.pm:1444-1447` (115): after `return 0 unless defined $bvfo->cds_start && defined $bvfo->cds_end;`, goes straight to `my $var_len = $bvfo->cds_end - $bvfo->cds_start + 1;` — no peptide check. | `VariationEffect.pm:1551-1556` (116): inserts `my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_); return 0 if defined $ref_pep && $ref_pep =~ /^\*/; # if the first base affected is the stop codon then it does no affect the reading frame` between the `cds_start`/`cds_end` guard and the `$var_len` computation. | New guard: if the reference peptide at the affected position **starts with** a stop codon (i.e. the variant begins at or after the annotated stop), 116 no longer calls it `frameshift` — regardless of what `abs($allele_len - $var_len) % 3` would have computed. In 115, an indel landing exactly at/after the stop (e.g. into the 3' extension/read-through region) could still satisfy the length-mod-3 arithmetic and be flagged `frameshift` purely by CDS-length bookkeeping even though it doesn't touch the coding reading frame in the intended biological sense. | `frameshift` (and by exclusion, whatever else the allele gets classified as, e.g. `stop_lost`/`stop_retained`/downstream/3'UTR terms) | **Yes** — directly named in the task brief ("frameshift" on the coding frontier). |

*(Trivial hunk 0, the `Copyright [2016-2025]` -> `[2016-2026]` line at `VariationEffect.pm:8` in both files, is license-only and omitted from numbering above.)*

---

## Table 2: Sequence.pm

| # | function | 115 behavior (quoted) | 116 behavior (quoted) | semantic change | HGVS impact | M6 frontier? |
|---|---|---|---|---|---|---|
| 1 | `get_3prime_seq_offset` | `Sequence.pm:801-802` (115): `my $check_length = length($post_seq) - length($seq_to_check);`. | `Sequence.pm:801-802` (116): `my $check_length = length($post_seq);`. | The unchanged loop consumes one downstream character per iteration while rotating the allele. Release 115 can therefore hit an artificial `flank_length - allele_length` cap before the first mismatch; 116 walks the full fetched flank. The executable call at `VariationFeature.pm:2003` feeds `VariationFeature::hgvs_genomic`. | Affects `g.` and, through that genomic path on mitochondrial sequence, `m.`. The import in `TranscriptVariationAllele.pm:76` is unused; `c./n.` uses its separate byte-identical `perform_shift` path. | M6c (`g.`) and M6d (`m.`); not M6b. |

*(Trivial hunk 0, the copyright-year bump at `Sequence.pm:8`, omitted.)*

---

## Table 3: TranscriptVariationAllele.pm

| # | function/area | 115 behavior (quoted) | 116 behavior (quoted) | semantic change | c./n. coordinate impact | M6 frontier? |
|---|---|---|---|---|---|---|
| 1 | `dbnsfp_alphamissense_prediction` / `dbnsfp_esm1b_prediction` / `dbnsfp_alphamissense_score` / `dbnsfp_esm1b_score` (new accessors) | Do not exist in 115. | `TranscriptVariationAllele.pm:1097-1116` (`dbnsfp_alphamissense_prediction`, `dbnsfp_esm1b_prediction`) and `:1217-1236` (`dbnsfp_alphamissense_score`, `dbnsfp_esm1b_score`) (116, two new blocks): four new thin wrapper subs, each delegating to the existing generic `_prediction('dbnsfp_alphamissense_prediction', ...)` / `_score('dbnsfp_alphamissense_score', ...)` machinery (same pattern as the pre-existing `dbnsfp_revel_prediction`/`dbnsfp_meta_lr_prediction` etc. immediately before/after them). | Purely additive: new named dbNSFP-backed missense-pathogenicity accessor methods (AlphaMissense, ESM1b) alongside the pre-existing REVEL/MetaLR/etc. accessors. No change to any existing method's logic; no change to SO-term computation or coordinate resolution. | none | Not a consequence predicate at all — dbNSFP score plumbing, not on the SO/consequence frontier. Listed for completeness since it's a real hunk in the diff. |
| 2 | `_var2transcript_slice_coords` | `TranscriptVariationAllele.pm:2655-2663` (115): `# Return undef if this VariationFeature does not fall within the supplied feature.` / `return undef if ( $vf_start < 1 \|\| $vf_end < 1 \|\| $vf_start > ($tr_end - $tr_start + 1) \|\| $vf_end > ($tr_end - $tr_start + 1) );` / `return( $vf_start , $vf_end, $self->_transcript_feature_Slice($tr));` | `TranscriptVariationAllele.pm:2719-2729` (116): `# Check for overlap before clamping` / `my $tr_length = $tr_end - $tr_start + 1;` / `# Variant is entirely before or after transcript` / `return undef if (($vf_start < 1 && $vf_end < 1) \|\| ($vf_start > $tr_length && $vf_end > $tr_length));` / `# Clamp coordinates to transcript boundaries` / `my $clamped_start = ($vf_start < 1 ? 1 : $vf_start > $tr_length ? $tr_length : $vf_start);` / `my $clamped_end = ($vf_end < 1 ? 1 : $vf_end > $tr_length ? $tr_length : $vf_end);` / `return ($clamped_start, $clamped_end, $self->_transcript_feature_Slice($tr));` | Boundary-overlap semantics change. 115 rejects (`return undef`) as soon as **either** endpoint of the transcript-mapped variant coordinates falls outside `[1, tr_length]` — this includes variants that only *partially* overlap the transcript (one end inside, one end past the transcript's 5' or 3' edge). 116 only rejects when the variant lies **entirely** on one side (`both` start and end < 1, or `both` > `tr_length`); for a partial-overlap variant it now **clamps** the out-of-range endpoint to the nearest boundary (`1` or `tr_length`) and returns those clamped coordinates plus the transcript feature slice, instead of `undef`. | This function has two call sites, both unchanged in 116: (a) `TranscriptVariationAllele.pm:167` inside `_return_3prime`'s RefSeq/`_rna_edit`-transcript shift-hash-reuse path (that whole code path is gated behind `return if ($tr->stable_id =~ /ENS/);` at `TranscriptVariationAllele.pm:150`, i.e. RefSeq-only) — in 115 an `undef` return falls into the `else { ... }` branch (`TranscriptVariationAllele.pm:201-206`) which just falls back to `$vf->{shift_hash}`/`shift_hash_reverse}` (may itself be undef); 116 instead gets a usable clamped slice and can attempt the flanking-sequence comparison instead of skipping straight to the fallback. (b) `TranscriptVariationAllele.pm:482` inside `look_for_slice_start`, which is called generically (not RefSeq-gated) from `hgvs_transcript` (`TranscriptVariationAllele.pm:1389`, `sub hgvs_transcript` starts at `:1311`) for **all** transcripts when building `c.`/`n.` HGVS strings. In 115, if `_var2transcript_slice_coords` returned `undef` for a boundary-straddling variant, `look_for_slice_start` would leave `$self->{_slice_start}` unset, and `hgvs_transcript` would hit `unless (defined $self->{_slice_start}) { ...; return undef; }` (`TranscriptVariationAllele.pm:1391-1394`) — i.e. **no `c.`/`n.` HGVS string at all** for that allele/transcript pair. In 116, the clamped coordinates let `look_for_slice_start` populate `_slice_start`/`_slice_end`/`_slice`, so `hgvs_transcript` proceeds to call `hgvs_variant_notation(...)` with boundary-clamped transcript-slice positions instead of bailing out. | **Yes** — directly changes whether/how `c.`/`n.` coordinates are resolved for transcript-boundary-crossing variants; this is exactly the cDNA-position resolution behavior named in the task brief. |

*(Trivial hunk 0, the copyright-year bump at `TranscriptVariationAllele.pm:8`, omitted.)*

Milestone assignment: preserve partial-overlap provenance in the M6b context builder; the
VEP-compatible `c./n.` clamp/render behavior is gated in M6c.

---

## Implications for pinning 116

A consequence engine written against 115 semantics and pointed at 116-shaped conformance data will get several
concrete cases wrong, concentrated exactly where the M6a note said to look (stop/inframe/partial-codon predicates),
plus two areas the headline didn't call out by name but that the diff shows are real: `frameshift`'s new
post-stop guard, and the `X` (incomplete/ambiguous peptide) routing change shared by `stop_lost` and
`stop_retained`.

1. **Stop-codon-overlapping indels need the `_cil` ("consider insertion length") genomic-coordinate overlap test**,
   not the plain cDNA-coordinate one, when computing `stop_retained`'s fallback path (VariationEffect.pm hunks
   #6, #10, #12). A 115-shaped engine that reuses one cDNA-span overlap test for both `stop_lost` and
   `stop_retained` will misclassify insertions whose insertion *point* maps just outside the last-3-bases window
   but whose *inserted sequence*, once placed, extends into the stop codon. Note `stop_lost` itself still uses the
   **plain** (non-`_cil`) `_ins_del_stop_altered` in 116 (hunk #11) — the `_cil` variant is deliberately
   `stop_retained`-only; an implementer must not blanket-apply `_cil` everywhere.
2. **The final-stop-length heuristic must be replaced, not just its threshold tuned.** 115's `ref_eq_alt_sequence`
   used `length($final_stop) < 3` as a proxy for "does the inserted tail terminate in a stop"; 116 checks the tail
   literally starts with `*` (hunk #8). An engine that hand-ports the `<3` cutoff, or drops the removed
   first-clause `($ref_pep eq substr($alt_pep,0,1) && $alt_pep =~ /\*/)` incorrectly, will diverge on
   `stop_retained` vs. `stop_gained`/`inframe_insertion` classification for stop-adjacent insertions.
3. **`stop_gained`, `stop_lost`, and `stop_retained` must be evaluated in the 116 call order**: `stop_gained` now
   explicitly defers to `stop_lost` (hunk #3), `stop_lost` now short-circuits on `partial_codon` unconditionally
   (hunk #4), and both `stop_lost` and `stop_retained` now divert alleles whose alt peptide contains `X`
   (incomplete/ambiguous translation) away from the direct string-comparison path and into the overlap/alteration
   fallback (hunks #5, #6). An engine that computes these three predicates independently/in parallel rather than
   with 116's dependency order will produce cases where more than one of `stop_gained`/`stop_lost`/`stop_retained`
   fires for the same allele, or where an `X`-containing allele takes the wrong branch.
4. **`inframe_insertion`/`inframe_deletion` must exclude alleles whose reference codon is already a stop** (hunks
   #1, #2) so that stop-codon-only-length-changing variants (`TAG`->`TAAG`, `TAG`->`G`, etc.) are left for the
   stop-specific predicates rather than being independently claimed by the inframe indel predicates.
5. **`frameshift` must not fire once the reference peptide at the affected position is itself a stop** (hunk #13).
   A 115-style engine that classifies frameshift purely from `abs(allele_len - cds_len) % 3` will wrongly tag
   post-stop indels as frameshift.
6. **Genomic HGVS 3' shifting (`get_3prime_seq_offset`) must not cap its walk at
   `flank_length - allele_length`.** A 115-ported implementation under-shifts indels before
   long repeats. Full-tree resolution proved this routine is executable for `g./m.` only;
   `c./n.` uses the separate, unchanged `perform_shift` path. Implement them separately in
   M6c/M6d, not in the M6b SO slice.
7. **`c.`/`n.` HGVS resolution for transcript-boundary-crossing variants must clamp, not reject.** 116's
   `_var2transcript_slice_coords` (Table 3 hunk #2) now returns clamped coordinates for variants that partially
   overlap a transcript's 5'/3' edge instead of `undef`. A 115-shaped engine will silently produce **no** `c.`/`n.`
   string at all for such alleles where 116 produces one (with boundary-clamped positions). This is a real
   coordinate-resolution behavior difference, not just an edge-case robustness fix, and it applies to **all**
   transcripts (not RefSeq-only) via the `hgvs_transcript` call path.
8. The M6a note's headline ("VariationEffect.pm differs by 155 lines, concentrated on stop_lost, partial_codon,
   stop_retained, `_overlaps_stop_codon`->`_overlaps_stop_codon_cil`, `_ins_del_stop_altered`->
   `_ins_del_stop_altered_cil`, and final-stop-length logic") is **directionally correct and not contradicted** by
   this catalog, but it undercounts scope in two ways worth flagging to implementers: it does not mention the
   `inframe_insertion`/`inframe_deletion` stop-codon guards (hunks #1-#2), the `stop_gained`/`stop_lost` mutual-
   exclusion ordering (hunk #3), or the `frameshift` post-stop guard (hunk #13) — all of which are separate,
   independently-landed predicate changes, not sub-parts of the four items it names. It also doesn't mention that
   `Sequence.pm` and `TranscriptVariationAllele.pm` both carry one real (non-cosmetic) hunk each, which the task
   brief's "reportedly changed" framing anticipated but the M6a note itself did not quantify.

---

## Escalated hunks — RESOLVED by gpt-5.6-sol (2026-07-10)

Both items were escalated to gpt-5.6-sol (max) on the persistent review session and resolved with full
call-chain evidence over the whole 115/116 tree:

**Q1 verdict — `$consider_ins_len` is DEAD scaffolding in VEP 116.** The only executable call to
`_ins_del_stop_altered` is `stop_lost`'s fallback, passing the bare 4-arg `@_`; the predicate dispatcher
(`BaseVariationFeatureOverlapAllele.pm:273`) invokes every predicate with exactly 4 args; the forwarded 5th
value reaches `_overlaps_stop_codon`, which destructures only 4 and silently drops extras. The live
insertion-length-aware logic is the *separate* `_cil` pair (`_overlaps_stop_codon_cil` /
`_ins_del_stop_altered_cil`) used by `stop_retained`. **Porting: implement the `_cil` path; ignore
`$consider_ins_len`.**

**Q2 verdict — the `get_3prime_seq_offset` cap change affects `g.`/`m.` HGVS ONLY, not `c./n.`** Transcript
HGVS uses a *separate, 115↔116-byte-identical* shift routine:
`hgvs_transcript → _return_3prime → _genomic_shift → perform_shift`, and `perform_shift` keeps its own
subtractive cap `length($post_seq) - $indel_length` (`TranscriptVariationAllele.pm:303`, unchanged). The only
executable call to `get_3prime_seq_offset` is `VariationFeature::hgvs_genomic` (`:2003`, the `g.` path); the
`TranscriptVariationAllele.pm:76` import is unused. **Porting rule: the two HGVS shift routines are DIFFERENT —
apply the changed `get_3prime_seq_offset` cap to `g./m.` only; `c./n.` retains the separate unchanged
`perform_shift` cap. An engine assuming one shared 3′-shift routine will get `c./n.` wrong.**

### Original escalation questions (historical context)

The two precisely located questions below explain why escalation was requested. The verdicts
above supersede their original uncertainty.

1. **VariationEffect.pm hunk #11 — `_ins_del_stop_altered`'s new `$consider_ins_len` parameter is unused dead code
   today, but I cannot confirm whether that is intentional.** `VariationEffect.pm:1433-1434` (116) adds a 5th
   positional parameter `$consider_ins_len` to `_ins_del_stop_altered`'s signature, and line `1444` forwards it as
   `_overlaps_stop_codon(@_, $consider_ins_len)`. `_overlaps_stop_codon` (`VariationEffect.pm:1363-1364`, both
   versions) only destructures `my ($bvfoa, $feat, $bvfo, $bvf) = @_;` — it never reads a 5th element, so the
   parameter is provably inert *given the source as written*. I traced every call site of `_ins_del_stop_altered`
   in both `VariationEffect.pm` files (`grep -n "_ins_del_stop_altered("`) and found none that pass a 5th
   argument — `stop_lost`'s call at `VariationEffect.pm:1263` (116) is still the bare `_ins_del_stop_altered(@_);`
   with the original 4-element `@_`. So as shipped in this exact 116.0-0 build, this is confirmed dead code, not a
   guess. What I could **not** resolve without wider access than this task scope (e.g. VEP's Runner/AF-plugin
   layer, or later 116.x patch releases, or upstream git history/commit messages I was not given access to) is
   *why* the parameter was added — whether it's leftover scaffolding from developing the `_cil` variants that
   should have been removed, or a forward-compat hook intentionally left for a caller outside the file I diffed
   (e.g. a plugin or CLI flag) that I have no visibility into. Recommend escalating only if the implementer's
   conformance harness sees `stop_lost` behavior on indels-near-stop that isn't explained by hunks #4/#5/#11 as
   documented here — that would indicate an external caller does pass a 5th arg somewhere I didn't find.

2. **Sequence.pm hunk — I could not find a direct call site proving `get_3prime_seq_offset`'s changed cap affects
   transcript-level (`c.`/`n.`) HGVS, only genomic (`g.`) HGVS.** `TranscriptVariationAllele.pm:76` (both versions)
   imports `get_3prime_seq_offset` via `use Bio::EnsEMBL::Variation::Utils::Sequence qw(hgvs_variant_notation
   format_hgvs_string get_3prime_seq_offset);`, but `grep -n "get_3prime_seq_offset("
   TranscriptVariationAllele.pm` returns **no matches** in either 115 or 116 — the import appears unused directly
   in this file (it may be called from a different transcript-HGVS module I was not asked to diff, or it may be
   genuinely vestigial). I confirmed the one concrete, provable call site is `VariationFeature.pm:2003` (both
   versions, `g.` HGVS path). I am not confident enough to state whether/how `c.`/`n.` shifting shares this same
   routine's changed loop bound, and I did not want to guess given Perl's `use ... qw(...)` doesn't itself prove
   usage. Recommend an implementer (or a follow-up grep across the full `ensembl-variation`/VEP module tree, which
   was out of this task's declared scope of exactly three files) confirm whether `c.`/`n.` 3'-shift math shares
   `get_3prime_seq_offset` or uses a separate transcript-level shift routine (`_return_3prime` /
   `perform_shift` / `create_shift_hash`, seen in `TranscriptVariationAllele.pm` around lines 118-465, were visible
   during this review but their internal shift-window-bound arithmetic was not diffed against 115 since they are
   textually identical between the two `TranscriptVariationAllele.pm` files — the diff tool confirmed no changes
   there).
