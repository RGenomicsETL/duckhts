# DuckVEP compatibility errata

Status: current implementation guidance. This ledger records non-obvious VEP 116
behaviour that DuckVEP intentionally reproduces. Code and executable differential tests
remain authoritative; this is not an implementation diary.

VEP compatibility describes an observed result, not an endorsement of the underlying
biology or API design. When VEP's predicate ordering produces an unusual combination of
terms, DuckVEP must reproduce that state before offering a separately named alternative.

## VEP removes EMAR rows before regulatory overlap evaluation

The Ensembl funcgen `regulatory_feature` table contains
`epigenetically_modified_region` rows whose feature-type name is EMAR. VEP 116 does not
turn those rows into `RegulatoryFeatureVariation` objects. Both its database annotation
source and the regulation-effect pipeline explicitly remove them before overlap
consequences are computed; the source comment says this avoids their very long names.
The indexed cache therefore cannot emit an EMAR stable ID even when an allele overlaps
the corresponding funcgen interval.

Loading every raw funcgen row produced a clean diagnostic signature: ordinary promoter,
enhancer, CTCF, open-chromatin, and motif overlaps agreed, while every extra DuckVEP row
was an EMAR `regulatory_region_variant`. The model compiler now removes EMAR rows before
dense feature ordinals, the resident SoA, and the cgranges index are built. This is source
selection, not a consequence suppression rule; keeping EMAR in the hot model and hiding
its output would waste memory and traversal work on features VEP never admits.

Source anchors: VEP 116 `AnnotationSource/Database/RegFeat.pm` and Ensembl Variation 116
`Pipeline/RegulationEffect.pm`.

## BND topology stays local while truncation observes both overlap alleles

VEP 116 does not evaluate a BND as two ordinary point variants. `BaseVCF4::get_start`
shifts the local VCF `POS` by one. `StructuralVariationFeature::_parse_breakends` keeps the
mate coordinate verbatim, and `StructuralVariationOverlap` can build an overlap allele for
each endpoint within its fixed 5 kb admission distance. Ordinary transcript predicates
still inspect the shifted local feature. Only `feature_truncation` inspects the current
overlap allele's endpoint.

The consequences are deliberately odd. A local intronic endpoint plus an intragenic mate
can produce `feature_truncation&intron_variant`. A transcript found only through an
interchromosomal mate produces `feature_truncation` with no local region. A mate near but
outside a transcript can create an internal overlap allele whose predicate list falls back
to `intergenic_variant`. The four VCF bracket orientations do not alter these transcript
consequence sets.

The fixed overlap-allele admission and the caller's directional window are independent.
For example, with `--distance 0`, a shifted local point that is outside a transcript but
within the fixed 5000-base admission range still creates a local overlap allele. Its
disabled upstream/downstream predicate leaves an empty predicate set, so that allele
defaults to `intergenic_variant`; an intragenic mate independently contributes
`feature_truncation`. VEP therefore returns
`feature_truncation&intergenic_variant`. At 5001 bases the local allele does not exist and
the same mate contributes only `feature_truncation`. Testing only the shared 5000-base
default hides this state.

The converse is equally non-obvious. If the mate creates the overlap object, predicates
on that mate allele still receive the local `StructuralVariationFeature`. A caller window
wider than 5000 bases can therefore emit a local upstream/downstream term even when the
local endpoint is too far away to create its own overlap allele. The fixed cap controls
allele construction; it does not clip the coordinates seen by ordinary predicates on an
allele constructed for the other endpoint.

Raw VEP output may contain two allele rows for one BND/transcript. DuckVEP supplies both
loci to one event, performs two cgranges candidate queries, evaluates ordinary topology
once from the local feature, applies mate-aware truncation, and emits the union once per
transcript. The rich SQL region is NULL when only the mate contributes; the compact region
mask is zero. Raw ALT, orientation, event identity, and provenance remain ordinary
relation columns for HGVS, fusion, and round-trip work.

VEP's executable oracle has a separate batching hazard. `InputBuffer::interval_tree`
inserts every mate position into the same coordinate tree as the local positions without
a chromosome key. When several cross-chromosome BNDs share a buffer, neighboring records
can therefore change another record's transcript set. This is not a per-event consequence
rule. The differential runs BNDs with `buffer_size=1` in one VEP process, preserving cache
reuse while isolating the event whose semantics are being compared. It also writes the VCF
with chromosomes contiguous and positions increasing, using the FASTA index only as the
deterministic chromosome order.

The broader seeded GRCh38 differential covers chromosomes 1, 2, 7, 21, and X, all four
bracket orientations, same- and cross-chromosome pairs, and transcript/exon/intron/CDS/flank
endpoint states. Its 1,004 BND events produced 91,428 transcript pairs, all exact against
isolated executable VEP 116 with no disagreement, extra row, or missing row. Source anchors
are VEP 116 `InputBuffer::interval_tree` and `BaseVCF4::get_start`, plus Ensembl Variation
116 `StructuralVariationFeature::_parse_breakends` and `StructuralVariationOverlap`.

## BND regulatory and motif overlap observes both endpoint points

The transcript asymmetry above does not mean regulation is local-only. VEP 116's
`AnnotationSource` fetches RegulatoryFeature and MotifFeature candidates around the local
and alternative breakend coordinates. `RegFeat::annotate_InputBuffer` then asks the input
buffer for overlapping variation features; that helper tests the shifted local point and
every mate point. Either point can therefore create a regulatory or motif overlap.

A BND has neither VEP's ordinary deletion predicate nor its copy-number-gain predicate,
but the result still depends on which endpoint found the feature. A shifted-local hit
retains the ordinary `regulatory_region_variant` or `TF_binding_site_variant` term.
`VariationEffect::feature_truncation` first checks `chromosome_breakpoint`; a mate-only
overlap reaches that branch and becomes generic HIGH-impact `feature_truncation` without
consulting deletion or copy loss.

There is a second, less obvious state. Once either endpoint discovers a feature exactly,
`StructuralVariationOverlap::_close_to_feature` admits the other endpoint when it is on
the same contig and no farther than 5000 bases from that feature. If the mate is inside
but the shifted local point is merely close, VEP emits two rows for the same stable
feature: local `intergenic_variant` and mate `feature_truncation`. DuckVEP's public unit
is their union, `feature_truncation&intergenic_variant`. This does not mean that a feature
5000 bases away is independently discoverable; one endpoint must still overlap it
exactly before the close endpoint is attached to that object.

This fixed value is not VEP's caller-configurable `--distance`, even though both default
to 5000. The former is compiled into
`StructuralVariationOverlap::_close_to_feature`; the latter controls directional
upstream/downstream transcript reach. DuckVEP therefore tests a 10,000-base caller window
with structural endpoints exactly 5,000 and 5,001 bases from the object. The first may
join the mate-discovered object; the second must not. Statistical sweep scenes also use
4,999, 5,000, 5,001, 10,000, and 65,535-base windows so the implementation is not trained
only on the shared default value.

If local and mate points both lie inside the same feature, the local base term wins. VEP
116 may emit several identical rows for that `(event, feature)` because
`InputBuffer::get_overlapping_vfs` can return the same StructuralVariationFeature through
both points. The executable witness at local raw `21:45553052`, shifted local
`21:45553053`, and mate `21:45553063` emits repeated
`regulatory_region_variant` for `ENSR21_5DP9TR`, but emits local
`intergenic_variant` plus mate `feature_truncation` for the nested motif
`ENSM00000587576`, whose exact interval starts at `21:45553055`. DuckVEP unions distinct
allele rows and returns each object once. It also rechecks contig plus point overlap, so a
matching numeric coordinate on the wrong contig is not an overlap.

Source anchors: VEP 116 `AnnotationSource.pm::get_all_features_by_InputBuffer`,
`AnnotationType/RegFeat.pm::annotate_InputBuffer`, and
`InputBuffer.pm::get_overlapping_vfs`; Ensembl Variation 116
`StructuralVariationOverlap::_close_to_feature` and
`Utils/VariationEffect.pm::{feature_truncation,within_regulatory_feature,within_motif_feature}`.

## Bounded tandem repeats change VEP object class; structural repeats mimic tandem duplication

For `<CNV:TR>`, VEP's VCF parser reads `RN` with `RUS`/`RUC` or `RB` and expands the repeat
to literal REF/ALT sequence when the resulting allele fits `max_sv_size`. The expanded
record becomes an ordinary `VariationFeature`; its consequence follows small-variant
normalization, sequence projection, and translation. An oversized or unexpanded record
remains a `StructuralVariationFeature` with structural class `tandem_repeat`.

That structural identity is observable metadata but not a distinct set of VEP consequence
predicates. `tandem_repeat` implies copy-number gain and insertion in
`VariationEffect.pm`, while the `duplication` predicate explicitly excludes tandem
repeats. The resulting transcript and interval-feature consequence terms are nevertheless
the same as a tandem duplication for the registered VEP-116 terms. DuckVEP therefore keeps
`STR` distinct in its event ABI but maps it to the same gain/insertion fact algebra; callers
retain repeat units and counts for provenance and HGVS.

Source anchors: VEP 116 `Parser/VCF.pm::_expand_tandem_repeat_allele_string` and Ensembl
Variation 116 `VariationEffect.pm::tandem_repeat`, `copy_number_gain`, `insertion`, and
`duplication`.

## Structural confidence and inserted-sequence payloads do not alter VEP 116 consequence terms

VEP's parser preserves `CIPOS`/`CIEND` as inner/outer structural coordinates, but its
release-116 candidate searches and registered `VariationEffect` predicates consume the
nominal start/end coordinates. Likewise, the structural branch of `inframe_insertion`
contains an explicit upstream TODO and returns false because it does not inspect inserted
sequence. These fields are still important evidence; they simply are not inputs to the
41-term consequence state machine being reproduced.

DuckVEP may therefore annotate nominal coordinates for strict VEP-116 consequence parity
while the surrounding relation retains confidence intervals and inserted sequence. It
must not promote the nominal span to exact experimental geometry, discard the payload, or
reuse it as HGVS, fusion, or round-trip geometry. Adding ignored payload parameters to the
hot consequence ABI would falsely imply semantics that VEP 116 does not have.

Source anchors: VEP 116 `Parser/VCF.pm` structural INFO handling and Ensembl Variation 116
`VariationEffect.pm::inframe_insertion` plus structural overlap predicates.

## Terminal-stop fallbacks bypass Translation SeqEdits

VEP applies supported Ensembl Translation SeqEdits to the reference peptide returned by
`TranscriptVariationAllele::peptide`, so its ordinary stop predicates inspect edited
local peptide text. The length-changing `X` or unavailable-peptide fallback is a separate
authority: `_overlaps_stop_codon_cil` tests the genomic coding endpoint, and
`_ins_del_stop_altered_cil` edits raw translateable DNA and translates the raw terminal
codon. It does not reapply Translation SeqEdits. A terminal deletion can consequently be
`stop_lost` through the CIL fallback even when a SeqEdit renders the corresponding local
reference residue `X`.

DuckVEP intentionally keeps the raw CDS terminal-stop gate for that fallback while using
the sparse edit overlay in direct peptide comparisons. Replacing both with one cleaner
peptide authority would change VEP 116 behavior. The focused C regression combines a
terminal `TAA`, an `X` reference-peptide edit, and a one-base terminal deletion to pin the
distinction.

## The insertion-length stop check is strand-asymmetric

VEP 116 has a second non-obvious state in the same fallback. For a minimized insertion,
the uploaded genomic interval is reversed: `(P+1,P)`. `_overlaps_stop_codon_cil` changes
the first coordinate by the inserted sequence length before testing the terminal codon.
On a reverse-strand transcript this turns an insertion just before the stop into an
ordinary interval that reaches the stop. On the forward strand the same arithmetic keeps
the interval reversed, so it does not create the corresponding upstream reach. This is an
observable strand asymmetry in the pinned implementation.

The subsequent `_ins_del_stop_altered_cil` path applies the insertion to translateable
CDS plus the 3-prime UTR, then retranslates the codon at the original CDS endpoint. The
extended overlap alone therefore does **not** imply stop retention: nearby insertions can
shift a non-stop codon into that position. The real GRCh38 witness
`7:148807690 C>CATCAGCCT` on reverse-strand transcript `ENST00001066230` happens to retain
the stop after that reconstruction. Its local peptides are `M/IG*X`, yet VEP emits
`protein_altering_variant&stop_retained_variant`, not
`frameshift_variant&stop_gained`. DuckVEP keeps transcript strand and uploaded-allele
orientation in both materialized and borrowed coding contexts, reproduces the exact
genomic coordinate calculation, and uses the shared CDS+tail endpoint reconstruction.
The reduced C witness uses the same inserted bases and peptide state; a plus-strand
“symmetry” rewrite or an overlap-implies-retention shortcut would be a compatibility
regression.

Source anchors: Ensembl Variation 116
`VariationEffect.pm::_overlaps_stop_codon_cil`, `::_ins_del_stop_altered_cil`,
`::stop_retained`, and `::frameshift`.

## `coding_unknown` and `missense` can coexist

VEP 116 does not treat an unknown peptide residue as a blanket failure. Its
`synonymous_variant` predicate rejects peptide strings containing `X`, but
`missense_variant` only asks whether the reference and alternate peptide strings have
equal length and differ, after the start/stop guards. Separately, `coding_unknown`
records that either local peptide contains `X`. Both predicates can therefore succeed
for the same allele.

The GRCh38 witness `17:75629084 GATGCCAGCAGA>TCTGCCTCTGGG` on
`ENST00000581825` has an incomplete first CDS codon. Ensembl represents its leading
bases as `NN`, giving local peptides `XLLAS/XPEAE`. VEP emits
`coding_sequence_variant&missense_variant`: the first residue remains unknown while the
later residues prove a missense change. DuckVEP permits `N` only in that declared
`cds_start_NF` first codon; ambiguity elsewhere still fails closed. The two resulting
facts remain independent rather than being cleaned into one biologically tidier label.

Source anchors: Ensembl Variation 116
`VariationEffect.pm::synonymous_variant`, `::missense_variant`, and `::coding_unknown`.

The converse matters at an incomplete terminal codon. When an uploaded feature starts
in a complete codon and reaches the trailing synthetic `X`, VEP may still prove a
missense change in the complete prefix, but it may **not** call an unchanged prefix
synonymous: `synonymous_variant` rejects the complete local peptide because it contains
`X`. The GRCh37 witness `2:228564238 CC>GG` on `ENST00000419059` therefore emits only
`coding_sequence_variant`, not `coding_sequence_variant&synonymous_variant`. Six
independently sampled GRCh37 transcript pairs exposed the same state. DuckVEP classifies
the complete-codon prefix, then applies the full-window `X` guard before retaining a
synonymous fact.

## A partial terminal codon is a sequence shape, not an attribute flag

VEP's `partial_codon` predicate is selected by the first affected peptide coordinate.
The durable model fact is the prepared CDS length modulo three, not merely Ensembl's
`cds_end_NF` attribute. A transcript whose prepared CDS ends in one or two bases can
therefore produce `incomplete_terminal_codon_variant` even when that attribute is absent;
this occurs in real mitochondrial and nuclear models.

Three nearby states must remain distinct:

- an edit beginning inside the one- or two-base terminal codon is
  `coding_sequence_variant&incomplete_terminal_codon_variant`;
- an edit beginning in a complete codon and continuing into the partial tail classifies
  the complete peptide prefix and independently records the trailing `X`;
- an equal-length feature that continues beyond the transcript's 3-prime edge still
  exposes its first mapped coding piece through `genomic2pep`, so the partial-codon term
  survives the outer mapper gap.

DuckVEP derives this from the prepared CDS byte length and the first mapped CDS position.
It does not promote every edit touching the final rounded codon to partial, and it does
not suppress the term just because one feature endpoint is a mapper gap. Fixed C cases
cover both strands, with and without `cds_end_NF`, and both internal and outer 3-prime
boundaries.

## A leading unknown codon does not suppress the frameshift predicate

For `cds_start_NF` transcripts, Ensembl may prefix the translateable sequence with one or
two synthetic `N` bases to preserve phase; the corresponding local peptide begins with
`X`. VEP still classifies a length-changing edit in that leading codon as a frameshift
from CDS edit geometry. Peptide-dependent start, missense, and stop facts remain absent,
but ambiguity is not a reason to discard the coordinate-level frameshift fact.

DuckVEP permits those synthetic `N` bases only in the declared incomplete first codon.
Its shared length-changing context resolves both a one-base insertion and one-base
deletion there to `frameshift_variant` on either strand, while ambiguity elsewhere still
fails closed. This is separate from the equal-length `coding_unknown&missense` state:
one concerns net frame displacement, the other compares equal-length peptide strings.

Source anchors: Ensembl Variation 116 `TranscriptVariationAllele.pm::peptide`,
`VariationEffect.pm::frameshift`, and the transcript `cds_start_NF` attribute path.

## An empty annotated UTR can overlap a spanning deletion

VEP 116 constructs its transcript-oriented before- and after-coding intervals as
`[transcript_start, cds_start - 1]` and `[cds_end + 1, transcript_end]`. When a CDS shares
the corresponding transcript endpoint, one of those intervals is inverted and contains
no annotated base. The generic `overlap()` predicate does not reject that empty interval;
a deletion spanning the endpoint can satisfy its four endpoint comparisons. The separate
cDNA-mapping requirement is also true, so VEP emits a UTR term together with
`coding_sequence_variant` even though the transcript has no UTR bases on that side.

The indexed VEP-116 cache exposed both orientations in the deterministic ClinVar chr21
sample. For example, `chr21:33602313 ATGTGCTAACTGAATAGCTATTG>A` on
`ENST00000417979` emits `3_prime_UTR_variant&coding_sequence_variant`. DuckVEP reproduces
the predicate geometry for a length-changing span; it does not manufacture a UTR sequence
or move the CDS endpoint. Complete feature-overlap states remain separate because VEP
routes them through feature ablation rather than these UTR predicates.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::_before_coding`,
`::_after_coding`, and `::overlap`. Keep the pure-C empty-interval witness and the real
cache differential together: either one alone is too easy to satisfy with an endpoint
special case that breaks the opposite strand or complete-overlap state.

Pure insertions have two additional explicit exceptions in those same predicates. VEP
stores an insertion after genomic base `P` as the reversed interval `(P+1,P)`.
`_before_coding` returns true when `P+1` is the CDS start, and `_after_coding` returns true
when `P` is the CDS end. The transcript mapper may accept the coding-side flank even when
the topology point used for the insertion is in the adjacent intron. The reverse-strand
GRCh37 witness `5:70356757 T>TG` on `ENST00000425596` consequently emits
`5_prime_UTR_variant&splice_region_variant`; treating the insertion as one ordinary point
loses the UTR term.

An equal-length feature can similarly cover the first CDS base and one base outside the
transcript. Its outer mapper endpoint is a `Gap`, so no start peptide exists, but the
inverted empty 5-prime UTR still overlaps and the mapped CDS base makes `within_cdna`
true. VEP's final result is supported
`5_prime_UTR_variant&coding_sequence_variant`, not a sequence failure. The full GRCh37
sample contained 36 such pairs on both strands. This outer 5-prime state must not be
generalized to the 3-prime edge: a partial terminal codon there remains classifiable as
described above.

## A mature miRNA term replaces the generic non-coding exon term

VEP 116 does not infer a mature miRNA from the transcript biotype alone. Ensembl stores
one or more `miRNA` transcript attributes whose values are inclusive cDNA ranges. VEP maps
each range through the transcript mapper, then tests the uploaded genomic span against the
resulting genomic segments. When the test succeeds, `mature_miRNA_variant` is selected at
an earlier consequence tier than the generic non-coding transcript and exon predicates;
those generic terms are therefore absent rather than added alongside it.

The deterministic ClinVar chromosome-21 cache sample exposed the final two differences:
`ENST00000290239` carries ranges `3-24` and `44-66`, while `ENST00000611994` carries
`41-61` and `6-28`. DuckVEP projects these cDNA ranges once during model preparation,
splits an interval at exon boundaries, and packs the resulting genomic segments as a
per-transcript side relation. The hot kernel checks that numeric slice; it does not parse
attributes or rebuild cDNA mappings per variant.

The cross-assembly GRCh37 cache run exposed an insertion-only boundary detail. VEP tests
the minimized `VariationFeature` coordinates, not the transcript mapper's chosen placement
point. A pure insertion after genomic base `P` has the reversed interval `(P+1,P)`;
therefore an insertion after the last mature-miRNA base does not overlap the mature range,
even though its retained VCF anchor does. Using the placement point produced 3,667 false
`mature_miRNA_variant` calls in the targeted GRCh37 corpus. The frozen boundary regression
keeps this distinction explicit for both transcript strands.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::within_mature_miRNA`, the
`miRNA` `attrib_type`, and the transcript mapper. Keep the rule-tier regression as well as
the importer-to-resident-model test: emitting both mature and generic exon terms would
find the right interval but still fail VEP compatibility.

## An empty predicate set becomes a transcript-associated intergenic result

VEP 116 does not drop a `TranscriptVariationAllele` merely because every consequence
predicate returns false. After evaluating the complete consequence list,
`BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences` replaces an empty list
with `$DEFAULT_OVERLAP_CONSEQUENCE`, whose term is `intergenic_variant`. The resulting
row still carries the transcript identifier and may carry coding coordinates and peptide
text. It is therefore not evidence that the uploaded variant lies between transcripts.

Seed 71 exposed 18 such terminal coding delins for `DUCK1-201`; for example,
`chrDuck:233 ACTGGTAA>AACCGGTTGACACTATTACTCATACCAATGGGTGC` has VEP codon and amino-acid
fields but the sole consequence `intergenic_variant`. DuckVEP assigns the same default
only after a real transcript candidate has exhausted its predicate set. This remains
separate from the SQL adapter's no-candidate behavior: only a receipt-backed complete
model may call a genuinely transcript-free input supported `intergenic_variant`; a
partial model returns `no_feature_in_loaded_model`.

Source anchors: Ensembl Variation 116
`BaseVariationFeatureOverlapAllele.pm::get_all_OverlapConsequences` and
`Utils/Constants.pm::$DEFAULT_OVERLAP_CONSEQUENCE`. Keep the transcript index in the
fixed regression. Replacing the missing row with a transcript-free synthetic row would
match the term while losing the upstream state.

## A mapper Gap can manufacture a retained stop

For a length-changing feature whose first or last transcript-mapper result is a `Gap`,
VEP cannot define both `cds_start` and `cds_end`, so ordinary codon and peptide alleles
do not exist. If any mapped part intersects CDS, the normal result is
`coding_sequence_variant` plus the independent topology and splice terms.

The terminal-stop path has a stranger outcome. `_overlaps_stop_codon_cil` ignores the
failed transcript endpoint and tests the uploaded genomic span directly. When that span
touches the terminal codon, `_ins_del_stop_altered_cil` then returns false because its
cDNA/CDS endpoint is missing. `stop_retained` negates that false result and emits
`stop_retained_variant` without ever editing or translating sequence. Held-out seed 71
exposed 2,237 instances of this state; `chrDuck:202` with a 50-base REF deletion is the
fixed witness. DuckVEP reproduces the predicate order rather than treating the missing
endpoint as proof that the biological stop was examined.

Source anchors: Ensembl Variation 116 `BaseTranscriptVariation::cds_start`,
`VariationEffect.pm::coding_unknown`, `::_overlaps_stop_codon_cil`,
`::_ins_del_stop_altered_cil`, and `::stop_retained`.

## A retained stop can still be protein-altering

VEP's consequence list is a set of independently evaluated predicates, not a clean
single-label hierarchy. For a length-changing insertion just before the terminal codon,
`_overlaps_stop_codon_cil` may reach the annotated endpoint and
`_ins_del_stop_altered_cil` may prove that the stop remains there. `stop_retained` then
suppresses `frameshift`, but it does not suppress `protein_altering_variant`. If the
local peptide strings have different lengths and the alternate preserves neither edge
of the reference peptide, VEP emits both `stop_retained_variant` and
`protein_altering_variant`.

The held-out randomized C run found reverse-strand examples with local peptide shapes
such as `M/IX`: the original endpoint still retranslates to stop, while the local
peptide shape satisfies the separate protein-altering predicate. DuckVEP keeps both
facts. The property oracle derives them separately from the local peptide strings; it
does not merely allow either result.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::stop_retained`,
`::frameshift`, and `::protein_altering_variant`.

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

## Ordinary literal deletions can ablate a complete transcript

VEP's `feature_ablation` predicate is not restricted to
`StructuralVariationFeatureOverlapAllele`. An ordinary VCF allele containing
only literal A/C/G/T bases remains a `VariationFeature` even when it is tens of
kilobases long. For length-changing alleles, the VCF parser removes the anchor
and `Parser::post_process_vfs` minimizes common allele edges; the resulting
ordinary deletion predicate is then combined with complete transcript overlap.
If both are true, the tier-1 result is `transcript_ablation`.

This state is common enough in graph/pangenome VCFs to matter. Treating all long
alleles as structural would change parser identity and other predicates, while
restricting ablation to the typed SV adapter loses the correct consequence.
DuckVEP therefore derives one normalized deletion fact for ordinary alleles and
sets the shared feature-ablation fact when that allele's VEP feature span
contains the complete transcript. Symbolic DEL/CNV loss continues through the
structural fact producer; the generated SO rule remains the single consumer.

Source anchors: VEP 116 `Parser/VCF.pm::create_VariationFeatures`,
`Parser.pm::post_process_vfs` / `::minimise`, and Ensembl Variation 116
`VariationEffect.pm::deletion` / `::feature_ablation`. Fixed tests must include
an ordinary shortened delins, not only a symbolic `<DEL>` record.

## Complete equal-length spans can expose both empty UTR predicates

VEP's `within_5_prime_utr` and `within_3_prime_utr` predicates call a generic
four-comparison overlap helper after requiring a cDNA mapping. They do not first
prove that the UTR interval contains a reference base. When a transcript and
its CDS share an endpoint, the corresponding UTR interval is inverted; an
equal-length uploaded feature containing the complete transcript can still
satisfy the overlap comparisons on both sides.

The HPRC long-literal witness at chromosome 22 position 22,919,681 has a
3,060-base REF and ALT and contains the 46-base `ENST00000390330` transcript,
whose transcript and CDS endpoints coincide. VEP 116 emits exactly
`3_prime_UTR_variant&5_prime_UTR_variant`. It emits neither
`coding_sequence_variant` nor a biological claim that UTR bases exist.
`VariationEffect::coding_unknown` is unconditionally false for complete feature
overlap, so an unknown-coding fact produced by the sequence view must not leak
back into the consequence set.

DuckVEP preserves this state in the shared span classifier on both transcript
strands, then makes complete-overlap suppression authoritative during effect
finalization. Do not replace the four comparisons with a non-empty-interval
test, and do not let an independently computed peptide/X state reintroduce
`coding_sequence_variant` afterward.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::within_5_prime_utr`,
`::within_3_prime_utr`, `::coding_unknown`, and the shared `overlap` helper.

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

The complete feature is also the reference-validation contract. Every retained REF base
must agree with the transcript CDS, and an ambiguous or incomplete codon anywhere in the selected
window suppresses the specific peptide predicate. Once that window is authoritative,
DuckVEP must preserve `reference_mismatch`, ambiguous-sequence, or `coding_sequence_variant`
state; retrying the smaller trimmed edit would annotate a different input representation.

Prepared CDS can prove that complete contract without genomic FASTA only when the complete
uploaded REF is also the minimized differing REF. If prefix/suffix trimming retained an
uploaded padding base or VCF anchor, validating only the differing CDS slice proves a
different, smaller assertion. A no-reference model must report `missing_reference` for
that HGVS row; it must not guess `reference_mismatch`, because the missing genomic base was
never observed. The regression `POS=124 REF=AA ALT=AC` over a transcript whose true bases
are `TA` is deliberately unresolved without FASTA and becomes an auditable mismatch only
when a reference provider is present.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::start_lost`,
`::start_retained`, `::stop_lost`, `::stop_retained`, `::stop_gained`, and
`::missense`. Fixed witnesses must remain paired with large held-out differentials so a
special-case term substitution cannot satisfy conformance.

## A five-prime boundary feature selects only the start codon

VEP 116 treats the opposite CDS boundary asymmetrically. When an equal-length uploaded
feature begins in the 5-prime UTR and continues into CDS, it does not suppress all peptide
predicates as the CDS-to-3-prime-UTR mapping gap does. Instead, the feature reaches the
start predicates, while the peptide view used there is limited to the first translated
codon. Coding bases later in the same uploaded feature do not produce an additional
missense or stop term through this path.

Four real VEP-116 witnesses separate the states:

| Uploaded feature | Coding change | Exact VEP terms |
| --- | --- | --- |
| `chrDuck:119 GA>AC` | start codon changes | `5_prime_UTR_variant&start_lost` |
| `chrDuck:118 CGA>GAA` | only UTR bases change; retained `A` overlaps CDS | `5_prime_UTR_variant&start_retained_variant` |
| `chrDuck:119 GATGGT>TATGAT` | `ATG` survives but the following coding bases change | `5_prime_UTR_variant&start_retained_variant` |
| `chrDuck:119 GATG>GTAA` | start codon becomes `TAA` | `5_prime_UTR_variant&start_lost`, without `stop_gained` |

The last two rows are the important guards: classifying the complete altered CDS would
add a coding or stop consequence that VEP does not emit. Conversely, treating the feature
like the 3-prime mapping gap would lose the start term. DuckVEP therefore retains the
uploaded topology, applies its contiguous CDS slice through the shared edit/context
engine, and selects codon one as the VEP peptide window. This is still one substitution
authority; there is no separate UTR consequence classifier.

Source anchors: Ensembl Variation 116
`BaseTranscriptVariation.pm::translation_coords`,
`TranscriptVariationAllele.pm::codon` / `::peptide`, and
`VariationEffect.pm::start_lost`, `::start_retained_variant`,
`::_snp_start_altered`. The fixed cases are generated and adjudicated by the real VEP
executable; randomized distributions remain the regression guard against overfitting.

## A length-changing start edit can be both lost and retained

For a feature crossing the 5-prime UTR/CDS boundary, VEP 116 evaluates two independent
string predicates. `_ins_del_start_altered` compares the edited UTR-plus-translateable
sequence with the original strings and feeds `start_retained_variant`. `start_lost` also
calls `_inv_start_altered`, even for an ordinary small allele; that helper asks whether
`ATG` remains at the original start offset. An edit can therefore preserve the translated
suffix while moving `ATG` away from that offset, making retained and lost true together.

The fixed VEP witness `chrDuck:117 ACGA>A` emits exactly
`start_lost&start_retained_variant&5_prime_UTR_variant`. This is why the resident model
keeps the complete spliced pre-CDS sequence and why the kernel evaluates the two predicates
separately over one allocation-free edited-sequence view. Do not collapse them into one
boolean or “correct” the combination. Models without complete transcript flanks return
`missing_transcript_flank`.

The same independence applies to `inframe_insertion`. VEP suppresses that term for a
start-retaining insertion only when the complete reference peptide is the *suffix* of the
alternate peptide, meaning that new residues were added before translation began. When the
reference peptide is the alternate *prefix*, the insertion follows the retained start and
`inframe_insertion&start_retained_variant` remains valid. Held-out seed 71 exposed the
minimal witness `chrDuck:121 T>TGTT`.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::_ins_del_start_altered`,
`::_inv_start_altered`, `::start_lost`, `::start_retained_variant`, and
`::inframe_insertion`. The reverse-strand kernel witness also pins VEP's genomic-left
anchor removal, which is not transcript-left on the negative strand.

## Length-changing predicates are independent

VEP 116 does not first assign a length-changing coding edit to one clean biological
class. `start_lost`, the three stop predicates, `frameshift`, `inframe_insertion`,
`inframe_deletion`, and `protein_altering_variant` independently inspect one local
codon/peptide state. DuckVEP therefore derives that state once and evaluates the same
predicate set; selecting one term from allele length or net frame loses valid results.

Four held-out VEP witnesses pin combinations that the older shape-specific paths left
unresolved:

| Uploaded feature | Exact VEP terms |
| --- | --- |
| `chrDuck:120 A>AAGTAAAATG` | `start_lost&stop_gained` |
| `chrDuck:129 ACGTACGTACGTACG>ACATAG` | `protein_altering_variant&stop_gained` |
| `chrDuck:202 CGTACGTACGTACGTACGTACGTACGTACGTACTGGTA>C` | `frameshift_variant&stop_lost` |
| `chrDuck:221 ACGT>A` | `inframe_deletion&stop_gained` |

The terminal shortening case follows another non-obvious guard in
`_ins_del_stop_altered`: VEP edits translateable CDS plus the complete 3-prime UTR and,
when that entire result is shorter than the original CDS, declares the original stop
altered without trying to read a replacement codon at the old endpoint. A complete
transcript flank proves this state; a model carrying only a short tail must fail closed.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::start_lost`, `::stop_gained`,
`::stop_lost`, `::frameshift`, `::inframe_deletion`,
`::protein_altering_variant`, and `::_ins_del_stop_altered`.

## Terminal overlap does not prove the stored codon is a stop

VEP's `_overlaps_stop_codon` is a coordinate predicate. For a complete annotated CDS it
asks whether the feature reaches the final three coding bases; it does not first require
the stored terminal codon or reference peptide residue to be `*`. The later
`_ins_del_stop_altered` comparison can consequently make `stop_lost` true for an edit at
an annotated coding endpoint whose prepared CDS actually ends in a non-stop codon.

The GRCh37 transcript `ENST00000599428` has a 60-base prepared CDS ending in `CCC`.
At `10:135341026`, the insertion `G>GT` is nevertheless
`frameshift_variant&stop_lost` in VEP 116. DuckVEP therefore separates “complete terminal
coordinate exists” from “reference terminal peptide is a stop.” Requiring `*` in the
coordinate gate loses this state; calling every endpoint edit stop-lost without running
the altered-endpoint comparison invents others.

Source anchors: Ensembl Variation 116 `VariationEffect.pm::_overlaps_stop_codon`,
`::_ins_del_stop_altered`, `::stop_lost`, and `::frameshift`.

## ALT-only mismatch islands inherit an expanded intron cache

VEP 116 can emit `intron_variant` for a lengthening replacement whose REF-shaped
`VariationFeature` remains entirely exonic. The result comes from two stages that are
easy to miss when reading one predicate in isolation:

- parser post-processing fully prefix/suffix-minimizes a length-changing biallelic pair,
  then defines the feature span from the remaining REF (which is empty for a pure
  insertion);
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

## The 12-base exon stretch is only a candidate set

`BaseTranscriptVariation::_overlapped_exons` stretches every exon by 12 bases when the
transcript contains any intron whose endpoint difference is at most 12. This is a coarse
lookup rule, not a replacement exon model. `non_coding_exon_variant` explicitly tests the
uploaded feature against each candidate's original exon coordinates before returning
true. A point inside the short intron is therefore neither a non-coding exon variant nor
an ordinary intron variant: `_intron_effects` recognizes the frameshift intron first and
skips its normal intron and splice predicates. The transcript-level non-coding fallback
remains.

The GRCh37 deletion `22:38616642 GG>G` minimizes to position `38616643`, inside a
four-base intron of `ENST00000541788`. VEP emits only
`non_coding_transcript_variant`. DuckVEP keeps the physical region intronic for reporting,
but keeps exact exon overlap and `within_intron` as separate predicate facts. The same
transcript-wide stretch still matters for consequence metadata whose `include` rule asks
for coarse exon exclusion, such as the polypyrimidine-tract gate; a remote short intron
can affect that candidate gate without turning every stretched base into an exon.

Source anchors: Ensembl Variation 116
`BaseTranscriptVariation.pm::_overlapped_exons`, `::_create_intron_trees`,
`BaseTranscriptVariationAllele.pm::_intron_effects`, and
`VariationEffect.pm::non_coding_exon_variant`.

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

The "consider insertion length" endpoint reconstruction is itself only a fallback.
`VariationEffect.pm::stop_retained` first asks whether the codon-local alternate peptide
is defined, non-empty, and contains no `X`. When it is, VEP calls
`ref_eq_alt_sequence` and does not consult `_ins_del_stop_altered_cil`, even if the
inserted length reaches the original CDS endpoint. A reverse-strand insertion can
therefore leave a stop at that original endpoint but expose the concrete local peptide
`*`; when the local reference peptide is `W`, VEP emits
`frameshift_variant&stop_gained` rather than `stop_retained_variant`. The held-out
statistical seed `20260719` found this distinction after 93,064 generated cases. Keep the
fixed `ATGGCATGGAGT`, CDS-position-9 `AT` insertion scene and the randomized oracle:
terminal-coordinate reconstruction applies only when the local alternate peptide is
empty or contains `X`.

Source anchors: Ensembl Variation 116
`VariationEffect.pm::stop_retained`, `::ref_eq_alt_sequence`,
`::_ins_del_stop_altered_cil`, `::stop_gained`, and `::frameshift`.

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

## The NMD plugin projects the full uploaded feature, not the minimized edit

Variant-induced NMD prediction in VEP Plugins release/116 does not consume the same
coordinates used to change CDS sequence. `NMD.pm` asks the parent
`TranscriptVariation` for CDS coordinates derived from the complete uploaded
`VariationFeature`, and reads the feature's genomic end for its exon-position rules.
Common-prefix/suffix minimization can therefore preserve the sequence change while
changing the NMD result.

Pure insertions expose the most surprising state. The parent transcript variation keeps an
empty feature as a reversed CDS interval such as `102,101`; it does not use the expanded
`101..102` allele range rendered in VEP JSON. At an exon edge, mapping one genomic flank
can be sufficient. Immediately before the first coding base, `1,0` is a valid result:
`NMD.pm` checks whether both coordinates are defined, not whether the numeric values are
nonzero, then treats the zero CDS end as an early-CDS escape.

DuckVEP consequently retains two projections:

- the minimized edit projection for REF validation, sequence application, translation,
  and core consequence predicates; and
- the full uploaded-feature projection for the NMD plugin's early-CDS, last-exon, and
  penultimate-exon-end rules.

A cached early-CDS fact from consequence prediction is reusable only when the full feature
and minimized edit cover identical genomic spans. Equal-length padded features and pure
insertions must use the full-feature projector. Broadening that cache condition changes
valid VEP results even when the underlying sequence edit is identical.

Source authority: VEP Plugins release/116 `NMD.pm`, commit
`0082591268417af618e03850c5ffdc7c09998a5d`, together with VEP 116
`BaseTranscriptVariation::cds_coords` and `TranscriptVariation` insertion mapping. Fixed
executable-plugin witnesses cover both strands, reversed insertion ranges, and the valid
`1,0` coordinate. The checked-in ClinVar chromosome-21 differential matches all 68,554
eligible transcript pairs, including 29,416 states both implementations leave unresolved.

## Structural coding predicates do not share one projection rule

VEP 116's structural consequence predicates reuse different pieces of mapper state. Five
cases found by seeded executable-VEP exploration must remain distinct:

- For an ordinary structural span, `_overlaps_start_codon` requires the first and last
  genomic positions to map to cDNA. A deletion or inversion can cover the genomic start
  codon but emit no start term when its other endpoint is intronic. When both endpoints are
  exonic, a structural allele has no inserted sequence, so VEP can emit the apparently
  contradictory pair `start_lost&start_retained_variant`.
- A structural insertion is different. `Mapper::map_insert` exposes its cDNA location as
  `start = end + 1` on both transcript strands. At an exon edge, its one mapped flank still
  determines which side of that base owns the insertion. Sorting the two coordinates or
  collapsing a one-flank result to a point invents start terms immediately before the CDS
  and after the third start-codon base. The reversed coordinate can nevertheless overlap a
  start codon split across exons, for example after its second base.
- structural `frameshift` requires a deletion wholly inside one exon and a length not
  divisible by three. It does not inspect `cds_coords`. Structural `inframe_deletion`
  additionally requires exactly one CDS `Coordinate` and no mapper `Gap`. A deletion that
  runs from CDS into 3-prime UTR can therefore be a frameshift when its length is not
  divisible by three, but only `coding_sequence_variant` when the same geometry has a
  length divisible by three.
- structural `stop_lost` calls the shared `partial_codon` predicate first. That predicate
  uses the event's transcript-oriented `translation_start`, not a transcript-wide
  "incomplete CDS" flag. On a 160-base prepared CDS, a deletion beginning in CDS base 158
  may still lose the stop because its first peptide coordinate names a complete codon;
  one beginning at CDS base 160 is suppressed as partial. If the event begins with a
  mapper `Gap`, `translation_start` is undefined and the partial-codon guard is false.
- `coding_transcript_variant` is restricted to the literal `protein_coding` biotype. Some
  Ensembl transcripts have a translation/CDS but another biotype; a neutral structural
  event containing such a transcript can exhaust the registry predicates and receive the
  transcript-associated default `intergenic_variant`.

These are compatibility facts, not a proposed biological cleanup. Fixed C witnesses pin
each state, while held-out generated structural events compare exact term sets against the
indexed VEP 116 cache. Source anchors: Ensembl Variation 116
`Mapper.pm::map_insert`, `BaseTranscriptVariation.pm::translation_start`,
`VariationEffect.pm::_overlaps_start_codon`, `partial_codon`, `stop_lost`, `frameshift`,
`inframe_deletion`, and the `coding_transcript_variant` registry entry in `Constants.pm`.

## Ensembl transcript HGVS first shifts against genomic sequence

VEP 116 does not normally discover the 3-prime `c.`/`n.` representation by walking the
spliced transcript. `TranscriptVariationAllele::_return_3prime` first calls
`_genomic_shift` in the transcript's strand direction. That routine fetches up to 1000
forward-reference bases on both sides of the event and calls `perform_shift`; the
resulting genomic offset and rotated allele are then used while the event is projected
back to transcript coordinates. A negative-strand transcript therefore walks toward
lower genomic coordinates, not toward the next byte in a stored CDS string.

The `perform_shift` limit is also not simply `min(1000, available bases)`. When the fetched
3-prime flank is at least as long as the inserted or deleted pattern, VEP compares at most
`flank_length - pattern_length + 1` positions. This allele-length-dependent remainder is
an executable compatibility detail. A byte-walk property covers both strands, rotated
multi-base alleles, duplications, and stops caused by that exact loop limit.

There is a separate sequence-region-end pathology. `_genomic_shift` expands the event to
1,000 bases on each genomic side, constrains that slice to the sequence region, and then
still defines `pre_seq` as the first 1,000 bytes and `post_seq` as the last 1,000 bytes of
the clipped slice. When the complete clipped slice is shorter than 2,000 bases, those two
substrings overlap. On a tiny contig, or close enough to either end of an ordinary contig,
`post_seq` can therefore begin before the event instead of at its 3-prime genomic flank.
The executable `chrDuck` fixture exposes shifted coordinates and rotated insertion payloads
that a clean non-overlapping flank implementation does not produce. This is a pinned VEP
implementation effect, not a general HGVS rule. DuckVEP must keep it as an explicit
compatibility case and must not train the ordinary interior-locus shift path on tiny-contig
output without modeling the clipped-slice construction itself.

Only RefSeq transcripts carrying relevant RNA-edit attributes may require a second
`perform_shift` against edited transcript sequence after VEP decides that a cached genomic
shift cannot be reused. Genomic `g.`/`m.` rendering is a third path:
`VariationFeature::get_all_hgvs_genomic` uses
`Utils::Sequence::get_3prime_seq_offset`. These paths share allele/event facts but must not
share one falsely generic sequence provider.

Source anchors: Ensembl Variation 116
`TranscriptVariationAllele::_return_3prime`, `_genomic_shift`, `perform_shift`, and
`VariationFeature::get_all_hgvs_genomic`, together with
`Utils::Sequence::get_3prime_seq_offset`.

The similarly named VEP options are different operations. `--shift_hgvs` controls this
HGVS-only 3-prime placement. `--shift_3prime` changes the internal transcript consequence
feature before consequence evaluation, and `--shift_genomic` rewrites the uploaded genomic
location toward the 3-prime representation. DuckVEP's independent-event HGVS surface
implements only the first operation; it must not move the consequence event or the caller's
variant identity as a side effect of printing HGVS.

## VEP HGVSc has an exonic-SNP-only phase fast path

VEP 116 does not use one coordinate routine for every transcript HGVS allele. Inside
`TranscriptVariationAllele::hgvs_transcript`, a literal `VariationFeature` whose
`var_class` is `SNP`, whose precomputed consequence state is exonic, and whose
`TranscriptVariation` has defined CDS endpoints takes a special path: both HGVS endpoints
are assigned directly from `TranscriptVariation::cds_start`. That value is already the
phase-aware CDS projection. Every other representation—including an intronic SNP, an
indel, an MNV, and an uploaded multi-base feature that trims to a one-base semantic
substitution—uses `_get_cDNA_position` instead.

The distinction is observable on `CDS_START_NF` transcripts whose first CDS exon has
positive phase. An exonic one-base `VariationFeature` includes that phase in its `c.`
coordinate, while an intronic SNP or a length-changing/multi-base feature at the same
transcript does not acquire a general phase offset. Applying phase in the generic
coordinate converter made a chromosome-21 ClinVar probe worse: the original 37 HGVSc
mismatches became 467, including shifted intronic and indel coordinates. The executable
state machine—not a uniform HGVS numbering rule—is therefore the compatibility authority.

DuckVEP keeps generic `_get_cDNA_position` coordinates phase-neutral and applies the
phase-aware coding projection only when the retained VEP feature is literally one
reference base and one alternate base, both endpoints are the same exonic coding base,
and coding projection succeeds. Fixed properties separately pin the literal-SNP and
multi-base-uploaded-feature cases so later cleanup cannot merge the two paths.

Source anchor: Ensembl Variation 116
`TranscriptVariationAllele::hgvs_transcript`, specifically its `var_class eq 'SNP'` plus
`pre_consequence_predicates->{exon}` branch and the following `_get_cDNA_position` branch.

## VEP caches a clamped pre-shift transcript slice for endpoint-crossing alleles

VEP 116 does not require every base of a deletion or delins to project inside a transcript
before it emits HGVSc/HGVSn. `TranscriptVariationAllele::_var2transcript_slice_coords`
first expresses the original `VariationFeature` relative to the transcript's genomic
slice. It returns no coordinates only when the complete feature is before or after that
slice. When a feature overlaps a transcript endpoint, it clamps the projected start and
end into `[1, transcript_length]`. Consequently, a four-base deletion with only the first
deleted base inside a transcript can be rendered as `n.1del`, and a longer deletion with
only two transcript bases retained can be rendered from those two bases rather than being
rejected as an invalid projection.

This state interacts non-obviously with HGVS 3-prime shifting. `hgvs_transcript` invokes
`_return_3prime` before lazily filling `_slice_start` and `_slice_end`, but
`_var2transcript_slice_coords` explicitly reads `unshifted_start` and `unshifted_end` when
the feature has been shifted. VEP then applies `_hgvs_offset` to the cached, clamped
coordinates when it calls `hgvs_variant_notation`; it does not replace the cache with a
fresh projection of the shifted feature. Executable VEP 116 therefore exposes the
original slice and the later shift at once:

- GRCh38 `21:45405086 CCCCGCCCCCTGCCCGGCCCCTGCCCGGCCCCTG>C` on
  `ENST00000859059` initially clamps the deleted feature to transcript slice bases 1 and
  2. VEP then reports an 18-base genomic HGVS shift, while the cached slice remains 1 and
  2, and renders `c.-95_-94del`.
- GRCh38 `21:31667246 TTTCA>T` on `ENST00000609934` leaves only one deleted base inside
  the transcript. The cached slice is `1..1`, the genomic shift is zero for this
  transcript, and VEP renders `n.1del`.
- GRCh38 `21:36930415 GAGGTAGTTCTAA>AGTTTGCACTATGTTGGAGTTTGCACTAT` on
  `ENST00000482273` clamps the feature to transcript slice `1..3` but retains the complete
  alternate allele state, producing
  `n.1_3delinsATAGTGCAAACTCCAACATAGTGCAAACT` on the negative strand.
- A pure insertion retains reversed transcript-slice coordinates. GRCh38
  `21:33262748 G>GT` is initially represented as `32781..32780` on
  `ENST00000980233`; after a one-base genomic shift the cached pair is unchanged and VEP
  renders `c.*860dup`.

The compatibility representation must therefore keep three concepts distinct: the
lossless uploaded feature, the semantic edit used by consequence and haplotype kernels,
and the clamped/cached transcript-slice facts used by VEP's HGVS formatter. Globally
clamping semantic edit endpoints, or teaching the shared transcript projector to accept
out-of-transcript coordinates, would change consequence, CDS, and phased-edit semantics.
DuckVEP therefore carries explicit transcript-slice start/end facts in its HGVS-facing
transcript edit, applies the separate HGVS shift offset during notation construction, and
leaves the semantic prepared event unchanged. Deterministic properties pin clipped
deletion, delins, and reversed insertion-slice states; the strict chromosome-21 ClinVar
differential exercises the same implementation without any HGVSc/HGVSp disagreement.

Source anchors: Ensembl Variation 116
`TranscriptVariationAllele::_var2transcript_slice_coords`, `look_for_slice_start`,
`_return_3prime`, and `TranscriptVariationAllele::hgvs_transcript`.

## Transcript-external consequence rows have no transcript HGVS

VEP can admit a transcript because an allele lies inside the configured upstream or
downstream distance while returning no `hgvsc` for that transcript allele. The absence is
not an HGVS projection error: neither edited genomic coordinate lies inside the transcript
span, so there is no transcript coordinate to print. A chromosome-19 GRCh38 GIAB probe
contained 767 such transcript pairs—404 `upstream_gene_variant` and 363
`downstream_gene_variant`—and executable VEP 116 returned no HGVSc for every pair.

DuckVEP consequently maps `DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT` to public HGVS
status `not_applicable` only when the complete semantic edit span—or both flanks of a pure
insertion—lies before or after the transcript. An event crossing a transcript endpoint is
different: its HGVS-facing transcript edit carries the clamped, cached slice facts from the
overlapping part of the original feature and may still render a supported string. A
printable VEP slice that is absent by design remains `not_applicable`; invalid exon/cDNA
projection, reference mismatch, missing reference sequence, and allocation-capacity
failures remain `unresolved`. Those states must not be collapsed merely because executable
VEP returned no string.

Source anchors: VEP 116
`TranscriptVariationAllele::hgvs_transcript` and the transcript-variation coordinate
projection used before it formats `c.` or `n.` notation.

## Protein HGVS preserves VEP's local-peptide state machine

VEP 116 protein HGVS is not a generic diff followed by a canonical formatter. Its helper
order makes several observable states that DuckVEP must retain even when a cleaner global
peptide comparison would choose another representation:

- HGVSp is attempted only when the pre-consequence `coding` predicate is true and both
  the first and last `TranscriptMapper::genomic2pep` items are coordinates. A feature may
  overlap CDS yet begin or end with a mapper Gap—for example at an intronic endpoint. VEP
  returns no HGVSp for that ordinary state; it is not a peptide-reconstruction failure.
- Equal local reference and alternate peptides return before shared-prefix/suffix clipping.
  The original translation start therefore remains the reported protein position.
- VEP converts a literal stop `*` to its internal `X`, BioPerl expands that symbol as
  `Xaa`, and the final formatter maps `Xaa` to `Ter`. A synthetic unknown residue and a
  translated stop can consequently share the rendered token even though their upstream
  facts differ.
- `_get_alternate_cds` calls `_trim_incomplete_codon` before appending the
  transcript-oriented 3-prime UTR used for continued translation, but that helper contains
  `if ($full_length = $keep_length)`. In executable Perl this assignment means edited CDS
  strings of one or two bases are trimmed to empty, while every length of at least three is
  returned wholly untrimmed—even when its length is not divisible by three. Replacing this
  with the biologically cleaner equality test changes frameshift/extension output.
- The first changed frameshift amino acid can occur after the first codon touched by the
  nucleotide edit. Its termination distance is measured from that first differing peptide
  residue to the new stop, not from the edit's first CDS coordinate.
- An insertion before the first translated residue, after the last translated residue, or
  immediately before the canonical terminal stop can have no protein HGVS. VEP's insertion
  representation requires two flanking residues, while its ordinary peptide string omits
  the terminal stop.
- `_shift_3prime` rotates a peptide insertion or deletion only when the complete changed
  peptide fits inside the remaining reference peptide. Its loop limit is
  `length(post_seq) - length(changed_peptide)`. A matching first residue is therefore not
  sufficient when the changed peptide is longer than the suffix; a cyclic-prefix search
  moves valid indels too far.
- Perl substring semantics leak into two start-of-peptide outputs. When clipping leaves an
  insertion at `start=1,end=0`, `substr(_peptide, -1, 2)` supplies the final reference
  residue and executable VEP can print strings such as
  `p.Trp0_Trp1insAsnThrAlaAlaThr`. A cached start-loss override over the same reversed
  interval can print `p.Trp1_?0`. These are VEP-116 compatibility strings, not valid HGVS
  positions, and the typed fact layer records them separately from semantic edit replay.
- The terminal stop becomes an insertion flank only for a pure insertion whose mapper
  coordinate lies inside the terminal codon. VEP may then print
  `p.Trp23_Ter24ins...`. An earlier coding insertion can clip to the same peptide-level
  insertion after the final non-stop residue, but VEP returns no HGVSp; globally appending
  a synthetic stop creates false protein strings.

The internal protein-HGVS facts keep equality, substitution, deletion, insertion, delins,
duplication, frameshift, start loss, and extension distinct before rendering. Randomized
properties replay the described peptide edit and publish counters for each observed shape,
strand, immediate-stop state, and unsupported terminal insertion. Executable-VEP strings
remain the final compatibility oracle.

VEP's ordinary `--hgvs` output renders the default protein suffix without parentheses,
for example `p.Ter394CysextTer9`. Parentheses are a presentation option enabled by
`--hgvsp_use_prediction`; they are not evidence of a different protein edit. The public
DuckVEP renderer follows the default unparenthesized form, while the internal renderer can
exercise both modes without changing the typed protein facts.

Source anchors: Ensembl Variation 116
`TranscriptVariationAllele::hgvs_protein`,
`BaseTranscriptVariation::translation_start`/`translation_coords`,
`TranscriptVariationAllele::_get_hgvs_protein_format`,
`_get_hgvs_protein_type`, `_get_hgvs_peptides`, `_clip_alleles`, `_get_fs_peptides`,
`_get_surrounding_peptides`, `_get_alternate_cds`, `_check_for_peptide_duplication`,
`_stop_loss_extra_AA`, and `_shift_3prime`.

## Consequence predicate sidecars preserve cached start- and stop-loss false states

The finalized SO set is not a closed-world serialization of every sequence predicate VEP
evaluated. Generalized edits that overlap splice logic can emit a splice/coding consequence
set without `frameshift_variant`, while the later shifted-CDS replay used by HGVSp still
proves a frameshift. Conversely, a consequence-side positive frameshift/start/stop
predicate remains useful when HGVS rebuilds the same edit.

DuckVEP therefore carries a compact predicate sidecar beside the consequence mask.
Frameshift and stop-retained bits are positive evidence: HGVS may OR them into its
independent CDS replay, but an absent bit must not clear a fact that the replay proved.
Deriving those HGVSp states only from the SO mask loses valid VEP combinations.

`start_lost` and `stop_lost` are the narrow exceptions. VEP evaluates and caches those
predicates while it constructs consequences, before `hgvs_protein()` performs its private
3-prime placement. The HGVSp call then reuses the cached values. The executable witness
`chrDuck:119 G>GAACT` is initially a `5_prime_UTR_variant`, shifts to
`c.1_2insACTA`, and renders `p.Met1AsnfsTer?`. Replaying the shifted CDS alone would infer
start loss and incorrectly render `p.Met1?`; VEP's cached `start_lost = false` prevents
that. The same cached complete-feature `coding` predicate also controls failure status. If
only that witness's uploaded anchor is replaced with a wrong base, genomic validation makes
HGVSc unresolved; HGVSp must be unresolved too because the feature was coding-admissible
even though its original region label was 5-prime UTR. Treating that region label as the
protein authority incorrectly returns not applicable. The same output keeps `fsTer?` even
when the shifted peptide replay reaches a stop:
`hgvs_protein()` deletes its private shift hash before `_stop_loss_extra_AA()` performs a
late alternate-CDS lookup, so the original UTR-side coordinates cannot prove a termination
distance. A deletion that removes the terminal peptide can likewise reconstruct stop loss
while the cached original predicate remains false, yielding an ordinary terminal
`p.Trp23_Ter24del` rather than an extension. Once the sidecar says sequence predicates are
valid, its start-loss and stop-loss bits are therefore closed-world even though its
frameshift and stop-retained bits are not.

The start-loss cache has a separate mapper trap. `Mapper::map_insert` maps the two genomic
flanks, drops an intronic `Gap`, and moves the surviving cDNA coordinate inward to retain
`start=end+1`. An insertion beside a start codon split across exons can consequently have
unshifted cDNA coordinates `3,2` and satisfy `_overlaps_start_codon` even though only one
genomic flank is exonic. Requiring two exonic endpoints clears a valid cached start-loss
fact. Fixed forward- and reverse-strand scenes cover both possible exon edges.

The pure-C regression deliberately restores a frameshift fact from a sidecar whose SO mask
contains only `coding_sequence_variant&splice_acceptor_variant`, then proves that an empty
sidecar cannot erase a separately reconstructed frameshift but does clear start and stop
loss. The public SQL regression pins the shifted 5-prime-UTR insertion and its exact
HGVSc/HGVSp.

This distinction is also a performance invariant. It is safe to synthesize the HGVS delta
from a valid sidecar when CDS length is unchanged, or when the sidecar positively records
frameshift. It is not safe to treat an absent frameshift bit as complete negative evidence
for a length-changing splice overlap. A strict chromosome-21 ClinVar differential exposed
three such rows: `21:44288343 AGG>A` on `ENST00000291582` and `ENST00000966178` must render
`p.Gly180AspfsTer36`, and `21:46117558 GCAGCCCAGCAGCCC>G` on `ENST00000984854` must render
`p.Ala359PhefsTer8`. Skipping the complete delta misrendered them as `p.Gly180Ter` and
`p.Ala359_His363delinsTer`. The centralized completeness predicate now admits the shortcut
only for an unchanged CDS length or positive frameshift evidence.

Source anchors: the separate VEP 116 consequence and `hgvs_protein` call paths,
`VariationEffect::start_lost`, `VariationEffect::stop_lost`, their `_predicate_cache`,
`Mapper::map_insert`, and DuckVEP
`duckvep_sequence_delta_consequence_flags` /
`duckvep_sequence_delta_apply_consequence_flags` /
`duckvep_sequence_delta_consequence_flags_complete_for_hgvs`.

## DNA duplication projects the copied source before requiring insertion flanks

VEP decides whether an inserted sequence is a tandem duplication before it formats an
ordinary insertion. For a duplication it prints the transcript coordinates of the copied
reference source. Only the fallback `ins` form requires both insertion flanks to project.
Reversing that order rejects a valid `dup` merely because the ordinary insertion notation
would lack one projected flank.

The low-level DuckVEP fact property pins this typed state with a copied final transcript
base rendering as `c.12dup`; it is a fact-layer reachability test, not a claim that the
current candidate sweep admits an insertion after a transcript endpoint. The executable
GRCh38 witness `21:33262748 G>GT` is end-adjacent rather than transcript-external: both
insertion flanks are in `ENST00000980233`, its copied source is projected after VEP's shift,
and the result is `c.*860dup`. Keep those two claims distinct when extending candidate
geometry.

Three genomic reference views must also remain distinct. `_genomic_shift` consumes the
exact sequence-region-constrained +/-1000 slice and takes its first and last 1000 bases as
the two shift strings. Complete uploaded VCF REF validation may require retained padding
outside that slice. Later, `hgvs_variant_notation` checks the transcript feature sequence
for an adjacent copied source whose length is the inserted allele length and is not capped
at 1000. DuckVEP performs one bounded faidx fetch but exposes the exact shift slice and the
wider assertion/duplication lookup as separate borrowed views. Widening the shift view for
retained REF changes VEP placement; limiting duplication lookup to it turns a valid
greater-than-1000-base `dup` into `ins`.

The randomized shift oracle exposed a second ownership lesson. Rare terminal duplications
can keep a printable copied-source DNA fact after the shifted insertion point has left a
synthetic transcript-sized reference view. Low-level CDS composition then depends on which
genomic VCF anchor that deliberately truncated view still contains, and forward/reverse
cases need not return the same low-level status. The public adapter first applies VEP's
shifted genomic-to-peptide endpoint precondition and does not compose that transcript-external
insertion into HGVSp. The statistical property therefore records the DNA-only state instead
of asserting an unreachable composition status; deterministic public tests own HGVSp
applicability.

Source anchors: Ensembl Variation 116 `TranscriptVariationAllele::hgvs_transcript`,
`TranscriptVariationAllele::_genomic_shift`, `hgvs_variant_notation`, and its duplication
check before insertion formatting.
