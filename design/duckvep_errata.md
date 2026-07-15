# DuckVEP compatibility errata

Status: current implementation guidance. This ledger records non-obvious VEP 116
behaviour that DuckVEP intentionally reproduces. Code and executable differential tests
remain authoritative; this is not an implementation diary.

VEP compatibility describes an observed result, not an endorsement of the underlying
biology or API design. When VEP's predicate ordering produces an unusual combination of
terms, DuckVEP must reproduce that state before offering a separately named alternative.

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

The complete feature is also the validation boundary. Every retained REF base must agree
with the transcript CDS, and an ambiguous or incomplete codon anywhere in the selected
window suppresses the specific peptide predicate. Once that window is authoritative,
DuckVEP must preserve `reference_mismatch`, ambiguous-sequence, or `coding_sequence_variant`
state; retrying the smaller trimmed edit would annotate a different input representation.

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
