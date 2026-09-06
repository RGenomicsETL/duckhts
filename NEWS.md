# DuckHTS Extension News

# duckhts 1.5.1.9000
- remove the redundant DuckVEP kernel model constructor in kernel 0.19.0: `duckvep_model_open`
  now accepts the optional interval-feature view directly. Native callers pass
  NULL for transcript-only models; model validation and SQL loading are unchanged

- require a pinned content checksum for VEP cache acquisition and verify the
  full compressed stream before publication; ETag/byte matches cannot substitute
  for content verification, and digest-process failures reject staging. Shared
  local-artifact SHA-256 and BSD-sum verification also requires a successful
  checksum process, even when failed commands print matching values

- share read-only projections of registry-built nested DuckVEP models between
  the conformance runner and FastVEP benchmark worker, removing their dependence
  on separately persisted flat views without rewriting model artifacts

- add registry-owned streaming acquisition of the complete VEP-116 GRCh38
  chromosome-21 cache for ClinVar acceptance campaigns, retaining all region
  shards and root metadata without storing the full upstream archive

- reuse validated DuckVEP indel CDS positions when the VEP feature selects the
  same genomic interval, avoiding repeated endpoint projection in consequence
  and HGVS context construction. One checked scalar phase converter serves
  base projections and indel positions; wider features retain full projection

- share DuckVEP's codon-rounded REF/ALT translation windows across substitutions
  and indels. Terminal SNVs/MNVs now borrow available 3-prime UTR bases for ALT,
  preserving independent partial-codon and stop-gained facts and recognizing
  stop retention when the full reference peptide is extended by a stop.
  Remove partial-substitution shortcuts and complete-codon-only clipping

- separate DuckVEP indel REF validation, VEP feature-CDS mutation, and unpadded
  cDNA start/stop predicates. Correct later-exon phase padding, known synthetic
  N codons, CDS-to-UTR edits and terminal-stop reconstruction on partial CDS.
  Standalone scratch calls now use the annotation dispatcher; physical edit
  arrays remain unchanged for sequence-edit and phased mechanics. CDS-to-UTR
  partial-codon shortcuts no longer bypass REF or required-flank validation

- separate complete MNV physical REF validation from VEP's feature-phase CDS
  edit, including retained bases and partial terminal codons. Share the
  independent start-offset test with SNVs and classify unchanged VEP peptide
  windows even when physical REF and ALT differ

- unify multi-base substitutions spanning 5-prime UTR and CDS with the shared
  transcript-string evaluator, correcting phase-padded start-loss/retention.
  Delete the clipped-CDS fallback and its start-only peptide mode: required
  missing UTR sequence is unresolved, and retained UTR REF bases are validated

- remove DuckVEP's insertion-only peptide special case inside a terminal
  partial codon. Keep codon-start insertions distinct from insertions between
  partial bases, correcting stop, insertion and protein-altering consequences
  while preserving the exact ClinVar HGVS witness. Codon windows include the
  borrowed 3-prime UTR and clamp REF and ALT requests independently

- separate DuckVEP's physical CDS reference checks from VEP feature/codon
  coordinates when a noncoding first exon precedes a phase-padded CDS. Correct
  native SNV body/stop consequences and terminal-codon admission without
  changing physical CDS edits or the SQL presentation reference. Do not reuse
  physical early-CDS NMD cache facts across a later coding-start phase

- retain DuckVEP's independent SNV start-codon predicate when phase-padded CDS
  sequence produces an unknown peptide, preserving physical REF validation.
  Add an opt-in native SO/status check to the R projection differential that
  retains every physical input record, including repeated alleles

- add SQL-first `duckvep_transcript_projection()` for typed VEP feature alleles,
  cDNA/CDS/protein ranges, exon/intron ordinals, transcript distance, quality
  flags and variable-length codon/amino-acid display from existing annotations
  and their prepared model; genetic-code tables remain kernel-owned and withheld
  coding sequence remains NULL rather than an empty-codon display

- connect the Haplosaurus differential to actual `read_geno()` calls while
  retaining both native comparison paths, the original verifier controls and
  source-bound build receipts; this is not a public phased-executor certificate

- add record-major `read_geno()` with typed arbitrary-ploidy GT/PS calls,
  per-allele phase flags and sparse non-reference calls. Share HTSlib sample
  selection with `read_bcf()` and expose stable original-header indices through
  `read_bcf_samples()`. Full genotype scans assign input-order, zero-based
  query-local ordinals; indexed region unions preserve physical duplicates

- check BCF output-list growth before accessing child storage and fail GT
  formatting on allocation errors instead of publishing a truncated genotype

- delegate BCF GT string formatting to HTSlib, retaining VCF 4.4 leading
  phase markers that the previous private formatter silently discarded

- stage an all-sample HPRC regional genotype cohort in VCF.gz and BCF through
  the benchmark registry, without changing consequence-only corpus derivations

- add `duckvep_breakend_geometry()` for lossless typed preparation of paired
  and single BND ALT alleles, preserving contig names, orientation, replacement
  sequence and telomeric mate coordinates with checked decimal parsing

- preserve VEP's per-endpoint default intergenic consequence when unioning
  paired BND transcript rows, including a close extragenic mate beside a local
  splice-site point with no ordinary structural consequence

- extend the R throughput driver with a raw BND ALT preparation workload,
  including the explicit mate-contig join and paired annotation in the timer

- bound native DuckVEP CDS projection by the borrowed transcript sequence
  slice before reading REF or insertion anchors, including hinted fallbacks.
  Add SQL/R checks for pre-publication rejection of short CDS models and feed
  genomic VCF alleles through the real projector in the carrier differential

- extend the R throughput benchmark with a checked-in-model indel workload
  that exercises CDS projection; record matched before/after bounds-check runs

- add a caller-owned native carrier index for phased-execution development:
  preserve explicit sample/phase-set/lane keys, share event prefixes, drain
  distinct occupied paths and recycle slots at transcript completion. Dense
  matrix properties and the Haplosaurus mechanics gate exercise the index;
  no public phased SQL function or genotype phase policy is introduced

- make the internal DuckVEP CDS rebuild reject coincident insertion sites,
  matching interaction partitioning instead of choosing ALT order from input.
  Exhaustive native tests compare both edit orders, both transcript strands,
  distinct scratch and exact-alias output, including unchanged buffers on conflict

- add a seeded executable Haplosaurus differential for native multi-edit CDS
  mutation and translation, comparing complete sequences, frame flags and
  contributor/sample counts. This developer gate does not certify the pending
  public phased executor, compound SO/HGVS or arbitrary-ploidy compatibility

- preserve literal colon-bearing contig names in indexed BCF/VCF full scans;
  internal contig tasks no longer parse dictionary names as region expressions.
  Share file/header/iterator ownership and decode diagnostics in a host-neutral
  scanner used by read_bcf and standalone tests, retaining worker-local mutable
  state and pending records across tidy output chunks

- stage the small BCF scan benchmark's committed inputs through duckhtsbench's
  registry, with pinned byte identities and a network-free cache-staging test

- retain parsed BCF/VCF indexes for the prepared plan, including region scans,
  instead of reloading them in workers. Index removal, corruption or valid
  replacement after bind now preserves the original rows and metadata counts.
  This supersedes the worker-reload error below; data and headers must remain
  unchanged, and the initial index must match the data. Reprepare to use a new
  index; sequential mode neither loads nor pins one

- plan indexed VCF full scans from the Tabix contig dictionary, preserving
  records on contigs absent from a partial header with a non-colocated index.
  Keep physical duplicate records and BCF header-ID semantics; VCF header-only
  empty contigs no longer create scan tasks

- fail index-dependent BCF/VCF parallel scans when a worker cannot reload its
  index, instead of silently returning zero rows. Preserve sequential fallback
  when no index was available at bind, and same-connection recovery after errors

- replace the BAM SIMD benchmark's Python driver with a shared R driver and
  fresh R/DuckDB processes per backend; retain workload controls and result
  checks, and invalidate cached measurements when binaries or inputs change

- report BAM CIGAR/AUX formatting and output-list growth failures as query
  errors instead of publishing missing, partial, or unwritable values. Share
  checked list growth and use reusable caller-owned text scratch with bounded
  native formatters with subtype-sized integer-array scratch; preserve genuine
  missing CIGAR/QUAL/tags and empty AUX strings

- restrict `read_bam` FILE_OFFSET to BGZF-compressed BAM: it is the virtual
  position immediately after the record, not its start. SAM (including compressed
  SAM), CRAM, and non-BGZF BAM return NULL instead of invalid transport-dependent
  values. Explicit ORDER BY FILE_OFFSET orders records within one unchanged BAM;
  scan arrival order is not guaranteed

- report allocation failures while constructing BAM/BCF/FASTA/FASTQ readers
  instead of crashing the host or publishing incomplete BCF schemas. Check
  array sizes and make partial cleanup safe, including FORMAT projection groups.
  Remove the redundant BAM sample-string copy and report failed SAMPLE_ID
  cache/HTSlib lookups instead of silently substituting NULL; genuinely absent
  samples remain NULL. Add isolated first-party allocation-failure/recovery tests

- reject empty region-list items and malformed known-contig intervals instead
  of silently widening or partially executing a request. Share one checked,
  host-neutral list parser across BAM, BCF, FASTA, and tabix readers; retain
  HTSlib coordinate semantics, unknown-contig handling, native record-union
  iteration, and FASTA's repeated-request rows. NULL/empty-string selection
  remains the no-filter API; list-allocation failure is a query error. FASTA
  region names come from the reader-owned index, preserving literal comma/colon
  header names instead of returning request-quoting braces or truncated names

- register the synthetic multi-region benchmark inputs and their network-free
  staging/test path; compare both reader revisions on the same staged artifacts

- remove experimental `read_bcf_appender` and its private connection,
  transaction, row formatter, and thread scheduler. Materialize `read_bcf`
  through DuckDB `CREATE TABLE AS`, `INSERT SELECT`, or `COPY`; retain canonical
  record/chunk ownership, sequential scanning, and native region semantics

- retain no-coordinate BAM/CRAM reads exactly once in automatic full-file
  contig scans, including across output chunks. Keep explicit region filtering
  unchanged; report read/iterator failures instead of silently ending a scan,
  and reject a lost worker index instead of duplicating sequential scans.
  Stream the tail when older BAI/CSI files omit its count and cannot locate it;
  do not use an ambiguous zero no-coordinate count for index-only counting

- remove the experimental `read_bcf_v2` SQL function and its generation-specific
  sample/schema selectors, text-only counting, and sample-run paths. Use
  `read_bcf`, SQL projection, and explicit `scan_mode := 'sequential'` for streaming;
  experimental sample-selection parameters are not accepted by `read_bcf`.
  Keep canonical header types, decode errors, worker-owned record caches, and
  HTSlib multi-region iteration unchanged

- cache Fedora Clang R dependencies under the exact build-profile inputs and
  use an explicit CRAN mirror for the Windows ARM64 dependency cache. Populate
  default-branch caches on main pushes; continue building and checking each
  current Rduckhts source tarball regardless of cache hits

- preserve CSQ/ANN/BCSQ values on every `read_bcf` tidy sample row, including
  records split across output chunks. Both BCF readers now retain one worker-owned
  decode cache; remove the duplicate per-chunk allocation/cleanup path and make
  partial cache initialization safe to destroy. Scan defaults remain unchanged

- share normalization and munging reference access through one host-neutral,
  thread-owned cache, retaining at most eight handles and eight 64-KiB sequence
  windows per thread. Oversized fetch results are not retained. Independent
  handles no longer serialize reference fetches; index creation remains locked.
  Keep caller-specific alias/end-of-contig, IUPAC, and remote-I/O tuning semantics;
  transport-tuning choices are part of cache identity. Test exact
  bytes, eviction, cleanup, and concurrent plain/BGZF access independently of SQL

- remove the unused DuckVEP variant-tile implementation, its isolated tests, and
  its extension/package build entries. The SQL annotation adapter already owns
  variant storage; public SQL behavior and the host-neutral kernel ABI are unchanged

- bound liftover clip-pad realignment to 4,194,304 Needleman–Wunsch cells,
  with checked matrix sizing and a recoverable DuckDB error instead of wrapped
  or unbounded input-driven allocation. Add replayable randomized cell-admission
  properties, malformed SQL/chain fuzzing that checks connection recovery, and
  complete-extension ASan/UBSan CI over the same focused host-safety corpus.
  Sanitizer CI installs a checksum-pinned DuckDB CLI, distribution tests
  consume the already built extension instead of rebuilding a relocated CMake
  cache, and MinGW jobs skip only the incompatible MSVC Python-wheel loader
  while retaining their native Windows contracts. The rendered liftover report
  keeps current synthetic evidence and the nearest comparable full-input
  real-callset baseline when current public inputs are not staged

- route malformed bcftools score-filter INFO/FORMAT storage and failed internal
  invariants through one recoverable filter authority instead of `exit()` or
  runtime assertions, unwind token stacks plus parser-owned index, file, and
  variable-argument temporaries after recovered failures, invoke `setjmp()`
  directly in a standards-conforming controlling expression, preserve the
  detailed filter cause in SQL errors, and remove the unused second recovery
  stack

- make Rduckhts SQL construction quote paths and option values through one
  DBI literal authority and table names through one identifier authority, while
  preserving explicitly documented SQL expression arguments. Quoted paths and
  identifiers can no longer terminate a wrapper-generated statement

- move HPRC v2, Sniffles2 1KGP, and dbVar GRCh38 chr22 DuckVEP corpus staging
  into the `duckhtsbench` artifact registry. The R-native CLI now owns pinned
  raw-source identities, remote indexes/manifests, deterministic derived VCFs
  and tabix indexes, cache paths, atomic publication, and adjacent receipts;
  retain the former positional shell command as a compatibility forwarder and
  cover version, ETag, manifest, path, checksum, and stale-cache behavior with
  a network-free fake-source integration test

- make the registered Ensembl 116 GRCh38 DuckVEP model reproducible from a clean
  cache: pin the public core and funcgen MySQL manifests, schemas, required table
  dumps, and primary-assembly FASTA in the benchmark registry; validate release,
  species, assembly, source manifests, reference sequence, and the deterministic
  DuckVEP receipt; preserve the imported source relation names; and build the
  model automatically before VariantKey provider exports. Add a network-free
  fixture-dump producer test and reject a mutated cached model

- make cgranges benchmark orchestration R-native with an `optparse` entry point,
  retain the former Python command as a compatibility launcher, and invalidate
  shared cgranges/GFFBase results whenever their source revision, extension,
  workload, or declared input identity changes. Resolve retained FastVEP input
  commands through the selected artifact registry; quote provider-stage paths;
  and reject an in-checkout wasm Docker work directory before recursive rsync.

- make the ordinary clone/build/test loop use the committed fixture closure: remove the
  obsolete instruction to regenerate test data from private mirrors, run SQLLogicTests
  through the cleanup-aware runner, and leave package builds network-free unless the
  developer explicitly chooses dependency bootstrap. Optional liftover, normalization, and
  VariantKey benchmark inputs are now explicit environment variables rather than one developer's paths.
  Explicit GIAB, liftover-reference, normalization, Riker, DuckVEP-corpus, local browser, and
  VariantKey-provider raw-source staging defaults now share `~/.cache/duckhts` (or
  `DUCKHTS_CACHE_DIR`), validate declared publisher identities before cache reuse, and record
  source/transformation provenance beside staged data. Retire only the unreferenced generic report;
  restore the whole-genome mosdepth report and a registry-backed `bench-mosdepth` entry point, and
  preserve the VariantKey and FastVEP reports as evidence while their declared public-provider
  staging work is made executable. Normalize embedded
  native source paths so equal clean builds can
  produce equal extension bytes.

- retire the duplicate R-side extension compiler, metadata serializer, and
  platform mapper; the package configure scripts and shared metadata writer are
  now the only build authority. Split the R package code into focused helper,
  connection, htslib, SIMD, catalog, and wrapper files

- make the Rduckhts binary architecture contract executable by compiling a
  minimal registered R native library. Binary repositories now derive a real
  target-specific `Built` field instead of classifying configure-built DuckHTS
  and htslib artifacts as architecture-independent

- make the nightly R-universe artifact sync validate the package `Built` field
  against the bundled DuckDB extension footer, publish only the receipted
  platform name, and fail on missing or contradictory architecture metadata

# duckhts 1.5.1
- mark Rduckhts as requiring compilation because its configure scripts build
  target-specific DuckHTS and htslib artifacts outside the conventional R
  `src/` directory. This prevents binary repositories from serving an amd64
  package artifact to Linux ARM64 consumers. Add a native Linux ARM64 clean-
  tarball install/load job alongside the Windows ARM64 package contract
- make Rduckhts package-owned DuckDB connections explicitly permit the locally
  compiled bundled extension while disabling implicit installation/loading of
  unrelated known extensions. This keeps package examples and tests executable
  with current `duckdb` builds that otherwise disable extension loading on
  Linux libc++ toolchains. File-backed connections reject a reused DuckDB
  driver whose creation policy cannot be changed, and failed new connections
  release their driver
- make native DuckHTS CMake builds and explicit CI profiles promote package-owned
  compiler diagnostics to errors while released R package builds preserve the
  warning policy supplied by R and the installation environment. Vendored htslib
  headers are system includes in Unix-like package builds, and the macOS package
  check now reproduces CRAN M1Mac's conversion-warning flags
- derive Windows package extension metadata from R's target and require the
  compiler target to agree, rather than trusting the MSYS shell architecture.
  Add native Windows ARM64 package coverage that loads the resulting
  `windows_arm64_mingw` extension and validates its bundled htslib receipt
- initialize the `bcftools_norm_row` shared out-of-memory cleanup count before
  any allocation can fail, avoiding an indeterminate cleanup-loop bound
- patch pinned libBigWig's `bwCleanup` definition to match its public
  `bwCleanup(void)` declaration, keeping Clang strict-prototype builds
  warning-clean
- reproduce VEP 116 on expanded gVCF input without per-row SQL regular
  expressions: omit `<NON_REF>`, bare `*`, and `.` ALT rows; retain real
  literal ALTs from mixed records; and annotate the distinct `<*>` catch-all
  allele as `coding_sequence_variant`, `start_retained_variant`, or
  `stop_retained_variant` without HGVS. Record-level `INFO/END` remains
  provenance rather than changing a literal ALT or `<*>` into a structural
  event. An executable differential splits and compares both `T,<*>` and
  `<*>,T` source orders against VEP 116 and pins its long-REF `<*>`
  allele-length quirk, which can emit an ablation term on complete feature
  overlap. The reusable native variant-tile builder accepts the same exact
  `<*>` event shape
- make pedantic extension builds warning-clean in DuckHTS-owned sources: initialize
  liftover alias storage, remove unused legacy helpers, retain portable allocation
  overflow checks, and omit ELF visibility attributes on MinGW
- recast the retained DuckVEP/FastVEP whole-genome comparison as an illustrated,
  evidence-first report. On the declared complete GIAB HG002 workload, DuckVEP
  records 2.55-fold and 2.12-fold same-core wall-time advantages at one and four
  pinned physical cores while matching all 56,998 held-out VEP 116 transcript
  pairs, consequence sets, and HGVSc/HGVSp suffixes. The report keeps its exact
  workload, output, run-count, and peak-RSS qualifications, directly benchmarks
  FastVEP fastSA and DuckDB on the same dated ClinVar payload and 4,095,611
  HG002 allele queries, and retains measured multi-provider AlphaMissense,
  assembly-correct REVEL, and phyloP BigWig workloads. It also explains the shared
  sweep/fact/rule design, late SQL presentation/provider joins, and layered test
  infrastructure, and plots the latest checked 5,100,500 randomized property
  trials plus 100,268-pair generated VEP differential. It does not extend the
  result to unmeasured predictors
- bounds-check htslib's in-place `BQ`/`ZQ` tag rename before writing the tag
  name, avoiding GCC 11's zero-size-region diagnostic while preserving BAQ
  behavior
- clarify the public DuckVEP model contract before release: the native kernel is
  transcript-source-neutral, while the supplied builder compiles the Ensembl core set;
  transcript-to-genome sequence corrections and circular-origin geometry remain
  unsupported, GENCODE Basic/Primary membership is retained for SQL filtering, and path
  projection or phased composition is never implicit
- make the Fedora CRAN-reproduction job install the suggested `Rtinycc` package before
  strict dependency checks instead of failing before Rduckhts is built
- document DuckVEP's exact sequence-path contract for X/Y pseudoautosomal regions,
  patches, and alternate haplotypes. Ensembl `assembly_exception` projections,
  `alt_allele` gene-equivalence metadata, and path-specific transcript rows are distinct;
  the current model loader performs no implicit path projection. Also reclassify official
  Ensembl Variation release-VCF `VE` comparisons as release-product audits rather than
  VEP executable oracles after a release-116 X/Y PAR witness showed the two products can
  differ while DuckVEP matches the pinned executable
- reproduce VEP 116's peptide-endpoint test before protein-HGVS 3-prime
  shifting. An in-frame insertion that copies the final translated residue can
  remain an HGVSp insertion rather than being promoted to a protein duplication,
  even when HGVSc is a transcript duplication. Add the minimized held-out
  differential witness and record the helper ordering in the compatibility errata
- make the executable rare-state campaign pass its generated maximum allele length
  through to the VEP differential, including the retained VCF anchor base. A 100-base
  campaign can no longer be reported from the runner's former 50-base default subset;
  the eligibility audit must retain the requested long-allele strata
- centralize the executable-language leaks required for Ensembl VEP 116 HGVS
  parity in one named compatibility policy instead of embedding release-specific
  codon-table and Perl-string behavior in individual algorithms. Consequence
  projection and HGVS now consume one borrowed per-event/transcript fact object,
  including one shared transcript edit and one shared terminal-partial-insertion
  state. Add `duckvep_allele_geometry(...)` so SQL provider joins can use the exact
  raw VCF span, VEP feature span, minimized edit span, and insertion point interpreted
  by the native engine. The held-out 100,000-trial UBSan campaign exposed an
  out-of-range read in the independent HGVS-shift test oracle; repair the oracle,
  preserve the failing seed, and require the same seed to pass before recording the
  campaign
- register the scalar and list forms of `duckhts_alt_to_list(...)` and
  `bcftools_norm_row(...)` as real DuckDB overload sets. Derived-query normalization
  can now consume the `VARCHAR[]` ALT column returned by `read_bcf()` instead of
  failing at bind time despite the documented list signature
- preflight the peak intermediate CDS length before applying an exact-alias
  multi-edit haplotype. A high-coordinate insertion followed by a lower-coordinate
  deletion can now return `BUFFER_TOO_SMALL` without mutating the caller buffer even
  when the final CDS length would fit. Also relabel the historical pair-level
  Clopper--Pearson calculation in the conformance report as descriptive under an
  independent-pair model; clustered transcript/object pairs and deliberately targeted
  generators are not presented as a deployment error-rate bound, while every observed
  mismatch, unresolved state, extra row, or missing row remains a release failure
- expose `duckhts_htslib_version()`, `duckhts_htslib_features()`, and
  `duckhts_htslib_feature_string()` as immutable runtime diagnostics for the htslib
  actually loaded with the extension. Rduckhts now publishes a versioned installed
  htslib linking receipt with exact headers/library paths, static dependencies,
  configured features, and source/build identity; its public accessor rejects
  source/header/runtime version disagreement and is exercised by an in-memory
  Rtinycc-built downstream C consumer that opens bundled BAM, CRAM, BCF, and VCF files,
  using the active SDK system headers on macOS. An omitted link mode selects the
  shared or static contract chosen when Rduckhts was configured; the compatibility
  `htslib_rpath()` helper follows that same configured contract and returns an empty
  loader path for static-only installations
- add `read_bigwig(path, region := NULL, blocks_per_iteration := 64)` for
  projection-aware scans of stored zero-based, half-open BigWig signal intervals.
  Comma-separated htslib-style regions are merged and deduplicated; full scans claim
  nonempty contigs and region scans claim merged ranges across DuckDB workers, with one
  mutable handle and iterator per worker. Vendor libBigWig 0.4.8 at commit
  `43c294ef1721a73b760803ca5e9410d581b98f17`, route its read-only I/O through DuckHTS's
  htslib `hFILE` transport for local, native remote, and browser-wasm access, and retain
  libBigWig's upstream correctness fixture and license. Oversized chromosome tables,
  malformed chromosome-tree IDs, truncated zoom headers, unreadable R-tree indexes,
  oversized block allocations, and truncated data blocks are rejected before indexed
  storage or interval decoding
- correct the real-WGS regulatory-provider plan after `EXPLAIN` showed that the
  documented chromosome equality forced a hash join with range residuals rather than
  IEJoin. Two packed RegionKey inequalities encode chromosome plus half-open coordinates
  exactly for the supported contigs and produce the identical 745,252 overlap pairs and
  414,813 matched alleles in 0.770 seconds instead of 83.960 seconds; the measured
  cgranges bulk alternative takes 1.144 seconds with lower peak RSS and remains the path
  for arbitrary contig names. Replace the 88-million-row retained consequence table and
  global ordered write with a direct four-writer Parquet stream under a 4 GB DuckDB memory
  limit. The full rich/HGVS plus ClinVar, ClinvArbitration, AlphaMissense, gene-constraint,
  and regulatory composition now writes 88,392,840 rows in 28.841 seconds with 5.31 GiB
  process high-water RSS; an order-independent full-row fingerprint matches the earlier
  46.87 GiB retained-intermediate run, which remains in the benchmark receipt as the
  diagnosed failure
- replace the root DuckVEP fixture-led narrative with a real human annotation flow over
  a public HG002 40x DeepVariant GRCh38 whole-genome callset, the complete Ensembl 116
  GRCh38 transcript/regulation model, collision-safe ClinVar and ClinvArbitration joins,
  AlphaMissense, gnomAD gene constraint, and an
  IEJoin-compatible regulatory interval provider. Show how each reusable provider Parquet
  relation is built from its declared release artifact, then present the case workflow as
  separately timed SQL stages and tabulate the resident-engine and FastVEP
  comparisons. Explain the continuing sorted transcript/feature sweep, immutable-model
  versus per-worker memory,
  relation-native supplementary-provider contract, and the fixed/property/statistical/
  executable/corpus conformance ladder. Record a timed whole-genome run over
  chromosomes 1--22, X, Y, and MT in the rendered throughput report, including its peak
  process memory, and distinguish extension-build, README-render, and optional VEP-oracle
  prerequisites
- document the end-to-end DuckVEP corpus and release-upgrade workflow: separate source,
  model, oracle, and comparison identities; preserve original-record-to-analysis-allele
  lineage and every pair denominator; inventory the generated, GIAB, ClinVar, GRCh37,
  non-human, regulatory, SV/BND, official Ensembl release-VCF, and HPRC lanes; and state
  exactly what is portable today versus what still needs a published external model/corpus
  pack. Limit large-artifact byte hashing to acquisition, publication, transfer, or an
  explicit release audit rather than presenting repeated full-file hashing as biological
  rigor. Add an optional `{targets}` campaign DAG that tracks
  corpus/model/reference/receipt paths, skips unchanged campaigns using its native file
  tracking, branches at coarse named campaigns, retains saved error workspaces, builds the
  clean release extension once, and reuses its revision-bound receipt across branches.
  Keep `blit` as the quoted shell-script execution layer for targets-to-runner and
  micromamba/VEP commands instead of growing a second command runner or generic digest
  cache; each executed evidence run hashes the artifacts it actually consumes, rejects
  caller-supplied precomputed digests, and publishes complete campaign outputs through an
  atomic directory replacement with interrupted-run recovery. Add focused contracts for
  traversal/symlink rejection, case-sensitive path identity, two-stage Windows batch
  quoting, atomic publication, and both recoverable interruption states. Require
  cache-backed campaigns to carry an acquisition/install receipt for the complete
  species/release/assembly cache leaf. The receipt binds the upstream checksum or immutable
  object identity plus `info.txt`, file count, byte total, and a deterministic
  relative-path/size/modification/change-time inventory. A cue-always cache-state target
  cheaply restates that inventory before reuse, detecting incomplete, replaced, or ordinarily edited shards
  without repeatedly hashing the 49 GB cache contents. Keep the staged
  annotation path separate from its stable published identity so methodology audits do
  not retain a temporary path after atomic publication
- close three rare consequence/HGVS states. VEP's start-loss predicate directly evaluates in-frame
  deletion before its generic start test; DuckVEP's staged classifier now establishes
  that dependent fact first, and the offset-based start test requires a 5-prime UTR. A pure
  insertion inside an incomplete terminal codon now exposes the empty-reference,
  insertion-only peptide to consequence predicates while preserving the codon-rounded
  edited-CDS view required by HGVSp. Fixed C witnesses plus exact public SQL and R
  regressions reproduce the remapped `ENST00000650713` state, including
  `c.280_281insAGT`, `p.Ter94=`, and its one-base 3-prime shift against a reference
  fixture that preserves the required downstream sequence. The large-corpus harness now
  spills pair evidence through Parquet,
  releases the resident model before comparison, asserts pair-key uniqueness, and keeps
  every HGVS discordance as an unconditional failure. The statistical executable-oracle
  audit now also fails on every consequence mismatch, missing/extra row, DuckVEP
  unresolved state, a status outside the explicit
  `supported`/`not_applicable`/`unresolved` vocabulary, or a status/value/reason
  inconsistency after preserving the complete pair and summary artifacts. Begin the
  machine-checkable whole-engine state model with a connected transition inventory that
  names observed dimensions, implementation and VEP authority, deterministic/property/
  executable evidence, and honest proof status. The generated check keeps haplotype
  mechanics and unimplemented combined classification visibly outside independent-event
  consequence claims. Make rare-state coverage requirements executable: selected
  statistical states must occur, zero counters require a named fixed C witness, and every
  failed campaign retains its complete seed-specific log. Treat held-out SV/BND and
  HPRC/pangenome campaigns as counterexample-guided state exploration rather than optional
  representativeness checks. Give newly discovered reachable rare classes dedicated
  generator strata instead of depending on accidental sampling. Document the pinned
  upstream VEP and Ensembl Variation test lineage as separate oracle-health and extracted-
  fixture gates. Preserve all 49 release-116 VEP tests, the two release-116 Variation
  consequence/HGVS authorities, and the five release-89 monolithic VEP tests under
  repository/ref-preserving paths with exact source commits and per-file SHA-256 receipts;
  the offline build check rejects changed, omitted, or unreceipted mirrors, while the
  release-source check compares every byte and complete suite inventory with the exact
  declared Git objects. Record the reproducible VEP-oracle self-test environment separately:
  the pinned 49-file suite passes 1,977 assertions with 293 explicit skips under a
  checked-in Linux-64 explicit Conda lock containing VEP 116.0, Perl 5.32.1, and BioPerl
  1.7.8; optional Perl modules can change the suite's dynamic assertion denominator, so the
  receipt preserves the complete solved environment without redundantly hashing its
  Git-tracked lock file.
  Add one explicit release-conformance audit target that requires the exact Git-backed
  mirror check as well as current executable-campaign evidence; the ordinary build remains
  network-free. Pin the two Variation semantic-source filenames rather than accepting any
  two tests from that release. Commit the verbose upstream TAP receipt and make its exact
  executed-file inventory, plans, assertions, per-file passes, skip reasons, and final
  result checksum- and content-verified; any `not ok` is fatal. Harden statistical
  history publication against excluded tests, untracked or mid-run source changes,
  concurrent campaigns with different ledger names, and process termination during
  journal cleanup; stale journals are recovered before the clean-worktree gate and one
  recoverable lock per destination directory prevents concurrent
  custom ledger pairs from replacing a shared output file. Official campaigns invoke the
  committed root Makefile with inherited GNU Make controls and optional property compiler
  flags cleared, then recheck source provenance after consuming the final evidence input. Dead
  publication owners are now detected from `pskill(..., signal = 0)`'s returned status,
  so their journals can actually be recovered.
  Define current campaign evidence by unchanged implementation and executable-conformance
  inputs since the measured ancestor revision, so the later commit that records a receipt
  does not make that receipt stale while any semantic source or harness change still does;
  the comparison now covers the complete extension/readers/vendor/build-helper closure,
  upstream semantic mirrors, and pinned build submodule. Staged, unstaged, compiler-like
  untracked inputs or a dirty build submodule make the strict audit fail.
  Checked executable campaigns likewise invoke explicit committed Makefiles with inherited
  GNU Make controls cleared, distclean vendored htslib, and reject ignored or untracked
  compiler inputs before rebuilding the receipt-bound extension.
  Clear inherited compiler, CMake, DuckDB, vcpkg, and toolchain controls from the checked
  build and regenerate its platform/version metadata through the committed configure helper
  before compilation.
  The smaller
  legacy monolithic VEP suite is not a hidden conformance authority, and VEP's
  “subversion” is its point-release number rather than SVN. Split
  alternative transcript, structural, interval, NMD, HGVS, and haplotype execution paths
  and require every implemented path's inputs and declared terminal state to be reachable.
  Resolve randomized properties, executable corpus campaigns, and explicit fail-closed
  error/status outcomes through checked manifests; report the unimplemented combined
  haplotype classifier as a structurally connected planned path, never as an executable
  one. Pin the rare complete-first-codon deletion and terminal partial-codon insertion
  through the public rich-plus-HGVS relation as well as the kernel and executable
  differential; the public reduction retains the no-UTR absence of `start_lost` and the
  apparently contradictory stop-gain consequence with `p.Ter5=` rendering.
- make HGVS differentials unconditionally fail closed: the corpus harness no longer
  accepts a flag that converts unresolved, mismatched, missing, or extra HGVS rows into a
  successful run. Reproduce two VEP 116 protein-formatting states exposed by the complete
  ClinVar mitochondrial shard: the late `fsTer`/`extTer` stop search uses BioPerl's
  implicit standard codon table even when the transcript peptide uses mitochondrial table
  2, and cached `stop_lost` plus a deletion-form peptide takes precedence over the later
  frameshift formatter. The fail-closed complete mitochondrial ClinVar rerun is exact for
  all 67,828 transcript pairs, including 3,294 present HGVSc and 2,354 present HGVSp
  values; the 25 initial HGVSp disagreements remain the discovery witness
- close rare HGVS states found by a 100,000-random-allele executable VEP 116
  differential, rather than by the C property suite alone. Shifted in-frame
  start-loss insertions perform VEP's sequence-dependent peptide-level 3-prime
  rotation before selecting flanking reference residues; a shifted frameshift
  at protein position 1 retains `Ter?`, while later positions can perform VEP's
  late stop search at the original CDS coordinate with the restored unrotated
  allele; and a
  minimized differing base just outside a transcript remains eligible for
  complete-feature clamping, allowing terminal `CG>CC` to render `c.*10dup`.
  Fixed SQL witnesses preserve the positive and negative rotation states and
  the exact HGVSc/HGVSp strings
- add a fail-closed official Ensembl release-VCF conformance runner for literal
  SNVs. It maps each lossless `VE` row to a VCF ALT through the producer's
  zero-based allele Index, aggregates
  consequence sets per transcript, runs the public rich DuckVEP relation against
  a receipt-matched model, and rejects mismatched, missing, extra, oracle-less,
  or unsupported rows. The documented non-SNV follow-up is pinned to Ensembl
  Variation release/116's GVF `Variant_seq`/`Index` plus `gvf2vcf.pl` lineage.
  Document that the release CSQ presentation overwrites repeated consequence
  terms for one allele/feature and therefore is not a full-set oracle; non-SNV
  ownership must not be guessed from padded VCF ALT strings
- let the deterministic DuckVEP corpus stager select every primary-contig tile
  for full GIAB/ClinVar performance gates, and record a caller-declared corpus
  name plus the dense-versus-full selection mode in its receipt. The historical
  density-ranked default remains unchanged
- extend the unified `duckvep_annotate(...)` relation with `rich := TRUE`.
  The fixed public schema retains the compact SO/region masks, IMPACT/status/
  reason codes, positions, amino-acid bytes, NMD codes, and object ordinals
  while adding VEP-facing consequence, IMPACT, region, amino-acid, NMD, and
  overlap-object text plus explicitly named `duckvep_status` and
  `duckvep_reason` audit fields. `rich := TRUE, hgvs := TRUE` fuses candidate
  discovery, consequence classification, and independent-event HGVSc/HGVSn/
  HGVSp in one native pass. Add mixed small/SV/BND and rich-plus-HGVS SQL
  regressions, extend the receipt-bound throughput driver to all four output
  choices, and replace the root README's synthetic-first DuckVEP overview with
  the Ensembl relation compiler, canonical event, conformance, performance,
  supplementary-join, and fair FastVEP comparison contracts
- make the unified `duckvep_annotate(...)` relation fail when a supported
  symbolic ALT disagrees with an explicit structural type, rather than
  silently reinterpreting the event, and make explicit `hgvs := NULL` follow
  the documented `FALSE` default instead of dropping small-variant rows
- replace the separate public small-variant, structural-event, breakend, compact,
  rich, and HGVS entry points with one pipe-friendly
  `duckvep_annotate(events_table, model_name, hgvs := false, ...)` relation.
  The canonical narrow event schema derives and validates each event family,
  dispatches to private native lanes, and returns one fixed compact schema with
  nullable HGVS fields. Add `duckvep_so_terms()` as the generated mask-decoding
  relation, public SQL/R regressions, and rendered root/package README examples
- measure all four public output contracts on all 4,095,611 model-addressable
  GIAB HG002 v4.2.1 literal alleles against the complete Ensembl 116 GRCh38
  model, including 1,383,580 regulatory/motif features. At checked revision
  `6c640890` the relation emits 47,835,851 rows at 964,354 compact, 346,030 rich,
  238,325 HGVS, and 180,495 fused rich-plus-HGVS alleles/s on one pinned core;
  four pinned cores reach 3,192,214, 1,513,530, 858,619, and 637,846 alleles/s.
  Five complete passes checksum the full output denominator and separate
  full-row fingerprints prove one-thread/four-thread equality for every output
  contract. See the rendered `benchmarks/duckvep_throughput.md` evidence
- extend the rendered VariantKey/RegionKey supplementary-annotation benchmark with
  isolated peak-RSS measurements and real dense exact providers. Official AlphaMissense
  v2 GRCh38 (71,697,556 rows), a receipted REVEL v1.3 GRCh37 projection (77,966,138
  rows), ClinvArbitration Zenodo record 16792026 GRCh38 (3,647,840 rows), ClinVar,
  GIAB HG002, Ensembl regulation, and gnomAD constraint are staged as typed Parquet
  projections. A deterministic genome-spanning REVEL key sample prevents favorable
  leading-row-group pruning. On one/four pinned i5-13500 performance cores, REVEL
  processes 1.27/4.44 million probe alleles/s, AlphaMissense 1.52/5.19 million/s, and
  ClinvArbitration 13.84/48.76 million/s. Fresh-process serving peaks remain 175--190
  MiB for those exact joins; one-time AlphaMissense and REVEL preparation peaks at
  1.10 and 1.25 GiB. The 643,528-interval string-labelled cgranges build plus full-GIAB
  probe peaks at 146 MiB versus 665/712 MiB for the per-contig, EXPLAIN-confirmed
  DuckDB IEJoin plan.
  Synthetic cgranges construction peaks at 129 MiB per million intervals without
  labels, 145 MiB with BIGINT labels, 197 MiB with VARCHAR labels, and 1.25 GiB for
  ten million BIGINT-labelled intervals. Narrow key-plus-payload storage projects to
  1.82--4.33 GB for 705,486,649 TOPMed Freeze 8 rows and 3.03--7.20 GB for
  1,172,689,405 live dbSNP Build 157 RefSNP rows. Twelve equally dense exact providers
  project to 27--32 minutes and 46--53 minutes respectively on four cores; these are
  explicit logical-time projections, not measured population runs. The accompanying
  design contract requires assembly/normalization receipts, collision-safe reversible
  and hashed lanes, sequential provider scans, and RSS-bounded chromosome/tile tasks
  rather than full-catalog cgranges indexes or twelve concurrently resident providers.
- add `duckhts_contig_key(...)` as the conservative chromosome/contig join authority for
  VCF, dbSNP, and DuckVEP preparation: remove one non-empty leading `chr` prefix,
  normalize mitochondrial `M`/`MT` spellings to `MT`, uppercase `X`/`Y`, and preserve all
  other suffixes. Numeric sex chromosomes, accessions, patches, and alternate loci are not
  guessed; DuckVEP callers must still reject model-side key collisions and the model
  compiler retains exact same-name, same-length reference matching
- add a rendered, receipt-pinned end-to-end DuckVEP/FastVEP comparison that
  includes compressed GIAB VCF decoding and real 17-column tab output. On one
  pinned i5-13500 performance core, DuckVEP also enforced the sorted-input
  contract inside the measured query and completed 4,048,342 records in a
  three-run median 64.54 seconds, versus 164.38 seconds for the current FastVEP
  native build. These are explicitly different native output contracts, not an
  engine-only speedup claim. The report separately controls FastVEP Rayon
  threads and DuckDB task threads. At four pinned physical performance cores,
  DuckVEP completes the same workload in a three-run median 32.06 seconds
  versus FastVEP's 68.09-second observation, a measured 2.12-fold
  same-core-count wall-time advantage. DuckVEP scales 2.01-fold from one to
  four threads at 251% median aggregate CPU use; measured process RSS rises
  only from 5.47 to 5.57 GiB, which rules out a fourfold total-RSS increase but
  is not per-allocation attribution. The report also records FastVEP's
  2.81-fold six-core speedup and 289% aggregate CPU use, retains output
  digests, and compares both engines with the same executable VEP 116 HGVS
  oracle instead of treating speed as evidence of semantic parity
- reset worker-local point and normalized-span exon ranks before reusing a workspace
  after a non-monotone prepared-event vector. A later variant can skip a transcript
  after an earlier forward jump, leaving that transcript's rank ahead even when the
  next DuckDB vector begins after the prior vector's final coordinate. A deterministic
  two-vector pure-C regression compares the reused workspace with a fresh workspace
  and preserves the exon consequence
- fix cumulative HGVS rendering for transcript or protein strings larger than
  the worker's initial scratch buffer. Renderers now report the exact required
  capacity and return `BUFFER_TOO_SMALL` before the adapter retries; the adapter
  also preserves DuckDB's checked UTF-8 rejection and validates every HGVS text-
  arena slice before materialization. Pure-C, SQL, and R regressions include a
  1,405-base inserted sequence matching the long-allele class found by the full
  ClinVar run and independent transcript/protein strings beyond the initial
  scratch capacity. At the exact fixed revision, pinned one-core complete-corpus runs
  process 4,438,467 literal ClinVar alleles at 391,848/s compact, 174,984/s rich,
  and 90,594/s with cumulative HGVS (126,733,707 rows), while 4,095,611 GIAB
  HG002 alleles run at 921,190/s, 447,266/s, and 244,222/s respectively
  (47,629,345 rows). The checked rendered report records five passes, the
  5,000-base transcript distance, CPU affinity, input and output denominators,
  and the absence of resident regulation/motif features. The throughput driver
  now also permits the cumulative HGVS surface to load the resident
  regulation/motif relation: interval-feature rows pass through the same
  consequence sweep with NULL HGVS fields, as pinned by SQL and R regressions,
  rather than requiring a second annotation scan. With all 1,383,580 resident
  RegulatoryFeature/MotifFeature intervals loaded, the full-corpus core-VEP
  rates are 387,944/s compact, 173,439/s rich, and 92,806/s HGVS for ClinVar,
  and 1,121,164/s, 468,230/s, and 249,004/s for GIAB. Each exact-head row now
  requires a vendored-htslib `distclean` followed by an in-tree release rebuild
  and retains its extension digest,
  physical/logical model, reference, staged-corpus, and original-VCF digests
  plus a full-public-row fingerprint. Executable HGVS pair artifacts carry the
  same source/binary/model/reference lineage, including separate source and exact
  sampled-VEP-input digests, and the checked history rejects
  unbound diagnostics. Capacity-one pure-C and greater-than-64-row SQL
  regressions prove that transcript, RegulatoryFeature, and MotifFeature
  observer callbacks resume without replay or loss when the result arena grows
- preserve declared source record IDs, `IMPRECISE`, `CIPOS`, and `CIEND` when the executable-
  VEP structural differential rebuilds its sampled VCF. A checked-in GRCh38 witness pairs
  nominal and imprecise CNV, DEL, DUP, tandem-DUP, INV, and INS records: VEP 116 and
  DuckVEP matched all 466 transcript pairs, and each engine produced identical nominal and
  imprecise consequence multisets for all six event kinds. This pins the intentional
  nominal-`POS`/`END` consequence contract while retaining confidence intervals as source
  evidence
- accelerate the shared DuckVEP consequence/HGVS execution path without adding a
  second biological authority: pack splice facts, use generated SO-impact masks,
  prepare normalized codons once, and let the cumulative HGVS adapter observe the
  consequence pass instead of replaying transcript classification. Stable active-set
  compaction now preserves `annotation_index` across DuckDB vector and disjoint ordered
  partition starts. The dense-region benchmark driver records arbitrary transcript
  distances, explicit ordered input partitions, full public-row multiset fingerprints, corpus
  receipts, and composition. On 517,097 annotation-dense GRCh38 alleles, one pinned
  i5-13500 P core emits 26,518,787 compact rows in 2.106 seconds at the 5,000-base
  distance (245,535 alleles/s); four pinned P cores and four disjoint ordered partitions
  take 0.632 seconds (818,191 alleles/s) with the same row-multiset fingerprint. Source
  ownership keeps the resident model immutable and shared; GNU `time -v` observed process
  peak RSS rising from 5,446,084 to 5,453,280 KiB, without a model-sized step, but that
  process measurement is not allocation attribution. Zero-, 10,000-, and 50,000-base
  runs retain exact one/four-worker row-multiset fingerprints; the
  50,000-base case emits 88,784,213 rows with an 8.2 MiB four-worker RSS premium. The
  fixed VEP breakend overlap-allele admission distance remains separate from the
  caller-configurable transcript search distance
- add `duckvep_annotate_hgvs(...)`, a cumulative independent-event surface that returns
  the compact consequence row together with transcript `c.`/`n.` and default-VEP protein
  `p.` HGVS suffixes, structured status/reason fields, and the applied transcript-direction
  shift without repeating candidate discovery. `duckvep_model_load(...)` can now bind an
  existing indexed reference FASTA through an exact ordinal/name/length relation; model
  validation never creates or modifies an index and pins open read descriptors for the
  validated FASTA, `.fai`, and optional `.gzi`. Linux workers reopen those descriptors,
  Windows retains a resolved source under deny-write sharing, and other POSIX workers open
  the resolved source independently with identity checks rather than sharing `/dev/fd` seek
  state. Each annotation worker owns
  its mutable faidx handle and a forward-read-ahead reference window reused by contained
  sorted alleles. External in-place mutation of a loaded source is outside the model
  contract and is rejected when descriptor metadata or the resolved source identity changes
  around a fetch. Explicit NULL
  for an optional model relation or reference is equivalent to omission. One prepared event/CDS edit authority
  feeds consequence evaluation and the HGVS-facing transcript edit; DNA placement, repeat
  rotation, and alternate-protein reconstruction then consume that same edit, while bounded
  allocation-free renderers keep VEP's default unparenthesized HGVSp distinct from its
  optional prediction display. Deterministic, SQL, R/DBI, and randomized pure-C coverage
  spans both strands, reference mismatch, missing-reference state, repeat shifting, mixed
  models in one DuckDB vector, ordinary peptide replay, and extended frameshift translation.
  Transcript rows admitted only through an upstream/downstream distance now report HGVS
  `not_applicable`, matching VEP's absent HGVSc, instead of presenting transcript-external
  edit geometry as an unresolved projection. HGVSc numbering also preserves VEP 116's
  representation-specific phase asymmetry: a literal exonic SNP uses the phase-aware CDS
  start, while intronic SNPs, indels, and multi-base features continue through VEP's
  phase-neutral generic transcript-coordinate path. Protein HGVS now returns
  `not_applicable` when VEP's first or last genomic-to-peptide mapper item is a Gap, instead
  of presenting an ordinary CDS-overlap/no-HGVSp state as a reconstruction failure. When
  transcript HGVS fails reference validation, protein status now uses VEP's cached
  complete-feature coding predicate rather than the original region label; a 5-prime-UTR
  insertion that could shift into coding coordinates remains unresolved instead of being
  mislabeled not applicable.
  A no-reference model now fails closed as `missing_reference` when retained uploaded REF
  padding or an anchor cannot be checked by the prepared CDS; it cannot validate only the
  minimized differing REF and emit supported HGVS for a different uploaded allele. Protein
  replay also reproduces VEP 116's `_trim_incomplete_codon` assignment bug: edited CDS
  strings of one or two bases are dropped before appending the 3-prime UTR, while every
  edited CDS of at least three bases remains wholly untrimmed. Reference-backed HGVS now
  keeps VEP's exact constrained +/-1000 genomic shift slice separate from the wider lookup
  view used to validate complete uploaded REF padding and identify adjacent copied sources;
  retained padding cannot change the shift, and insertions copied from more than 1000
  adjacent bases retain `dup` syntax instead of silently degrading to `ins`. Statistical runs publish
  named state-distribution counters so a passing total cannot hide unvisited edit shapes.
  A strict executable-VEP 116 chromosome-21 ClinVar differential matched all 56,998
  consequence rows and found zero HGVSc/HGVSp differences: 20,782 pairs matched both
  strings, 24,089 matched HGVSc with no protein string on either side, and 12,127 had
  neither string. This closes the exercised independent literal-small-variant matrix, not
  genomic HGVS, RefSeq RNA edits, exact structural/BND HGVS, or phased/compound HGVS;
  ferro-hgvs v0.9.0 is a pinned independent
  specification/corpus oracle, not a C/CRAN/Wasm dependency or replacement for VEP 116
  behavior. HGVS execution now consumes the kernel-prepared transcript, exon, and sequence
  views as its canonical model; derives and validates each reference CDS first-stop once at
  model open; projects the transcript edit before lazily attaching CDS edits; scans a
  virtual single-edit CDS without rebuilding the complete alternate sequence; and renders
  through one reusable worker buffer with retry only when the result exceeds its capacity.
  Consequence sidecars remain a cache rather than a second coding authority: an absent
  frameshift bit is conclusive only for length-preserving CDS edits, while length-changing
  splice overlaps run the full peptide delta unless the sidecar positively proves a
  frameshift. Exact VEP regressions include the ClinVar splice-overlapping deletions that
  render `p.Gly180AspfsTer36` and `p.Ala359PhefsTer8` despite lacking a consequence-side
  frameshift term. The production throughput driver accepts `chrom` as the sequence-region
  name column used by the final staged Ensembl model, so the HGVS workload is reproducible
  from that canonical staging relation rather than requiring a benchmark-only alias. On
  the pinned one-thread, 200,000-allele mixed-coding workload, cumulative HGVS execution
  falls from 32.106 to 4.579 seconds: 43,678 input alleles/s and 1,257,200 generated
  transcript-HGVS rows/s. The matched rich lane is 2.067 seconds, reducing the incremental
  HGVS tax from 15.79-fold to 2.22-fold while preserving the same 5,756,720 output rows and
  byte checksum. The sequence-pool first-stop cache extends the public standalone model
  view and advances the kernel ABI to 0.17.0
- match VEP 116 for complete transcript overlap by ordinary long literal
  alleles: a normalized deletion now enters the same tier-1
  `transcript_ablation` predicate as a symbolic deletion, while an equal-length
  uploaded span keeps VEP's four-comparison predicates on otherwise-empty UTR
  intervals and suppresses `coding_sequence_variant` even when sequence
  evaluation independently produced an unknown-coding fact. Add fixed C, SQL,
  and R witnesses plus checksum-pinned HPRC pangenome long-allele differential
  coverage; the configured 5,000-base upstream/downstream distance affects only
  candidate admission outside transcripts, not these overlapping states
- keep paired-BND result materialization linear in the dominant transcript-row
  count: transcript and interval-feature evaluators remain separate sorted
  streams, only the smaller feature stream is ordered by object kind, and a
  worker-owned merge buffer combines them in variant-major order. Models with
  no resident regulation features now retain the original transcript-only fast
  path instead of sorting an already ordered expansion
- tighten the DuckVEP statistical consequence oracle around VEP 116 terminal-codon
  insertions: the insertion-length-aware endpoint reconstruction is now modeled only for
  empty or `X`-containing local alternate peptides, while concrete local peptides retain
  VEP's `ref_eq_alt_sequence` path. A held-out reverse-strand regression discovered after
  93,064 generated cases now pins the resulting
  `frameshift_variant&stop_gained` state
- keep VEP's fixed 5000-base BND overlap-allele admission independent of the
  caller-selected transcript window in both directions: with a zero directional distance,
  an admitted local allele whose upstream/downstream predicate is disabled now contributes
  default `intergenic_variant` beside mate-derived `feature_truncation`. Add C, SQL, R, and
  executable-VEP regressions at zero, intermediate, default, and wider distances so the
  shared 5000-base default cannot hide this allele-row union
- close the remaining typed VEP-116 structural consequence inputs without inventing semantics for unused payloads: `duckvep_annotate_sv(...)` and its compact form now accept structural tandem repeats as `STR`/`TANDEM_REPEAT`/`CNV:TR`, preserve that event identity, and reproduce VEP's tandem-duplication gain/insertion predicates; bounded repeats expanded to literal alleles continue through the small-variant path. Paired-BND annotation now queries resident RegulatoryFeature and MotifFeature indexes at both the shifted local point and verbatim mate point and reproduces VEP's asymmetric feature semantics: local exact hits keep the base regulatory/motif term, mate-only exact hits use generic HIGH-impact `feature_truncation`, a same-contig shifted local point outside but within VEP's fixed 5000-base admission distance adds `intergenic_variant` to that mate-discovered object, and the local base term wins when both points hit one object exactly. Mixed transcript/regulation BND batches are restored to variant-major order before DuckDB list materialization, preserving one result list per input event. The executable oracle collapses byte-identical object rows, unions distinct allele-level consequences per object, and rejects conflicting non-consequence state. Add one shared explicit event/object-pair kernel path, SQL/R mixed-vector regressions, randomized all-feature BND/STR oracles, and executable-VEP breakend generation that targets feature starts, midpoints, and ends in raw-local, verbatim-mate, both-exact, mate-exact/local-close, exact-5000, and beyond-5000 forms. Fixed C/SQL/R regressions use a 10000-base caller window to prove it cannot widen VEP's separate fixed structural-allele cap; randomized sweep windows span 0 through 65535 bases instead of stopping at the shared default. Bump the standalone kernel ABI to 0.15.0. Exact VEP 116 source audit confirms that `CIPOS`/`CIEND` and structural inserted sequence do not alter its registered 41-term consequence predicates, so they remain relational provenance and future HGVS/round-trip inputs rather than ignored hot-kernel parameters
- keep unprojected FASTQ query names outside the BAM publishing limit: direct `read_fastq(...)` scans now accept long headers when `NAME` / `PAIR_ID` are neither returned nor needed for paired-file comparison, while name-producing scans retain the htslib-compatible limit. Also make an explicit `duckhts_fastq_qc(..., max_cycles)` value below 128 cap the aggregate's initial per-group cycle allocation instead of allocating the default 128-entry growth quantum
- add `duckhts_fastq_qc(sequence, quality [, max_cycles])`, a bounded aggregate over canonical FASTQ strings that emits exact global read/base/Q20/Q30/Q40/base-composition totals and compact per-cycle sufficient statistics without relationally expanding every base. The shared scalar kernel is the correctness oracle; AVX2, ARM NEON, and wasm SIMD128 slots use the existing per-kernel dispatch and preserve exact scalar results. Add SQL/R aggregate regressions, seeded standalone SIMD contracts, target-compilation checks, and a pinned one-core HG001 exome comparison against current fastp
- accelerate `read_fastq(...)` by parsing FASTQ records directly over htslib transport into only the projected DuckDB columns, removing the temporary `bam1_t` round trip while preserving multiline sequence/quality framing, htslib-style query names, paired/interleaved checks, string and packed nt16/Phred output, quality-encoding conversion, and count shortcuts. Truncated sequence or quality blocks now raise precise DuckDB errors instead of being mistaken for clean EOF. Add SQL/R multiline, packed-output, and truncation regressions plus a pinned one-core before/after benchmark and a separately labelled current-fastp QC comparison
- repair the installed R package's `duckhts_build(...)` fallback so it gets the DuckVEP kernel source list from the same helper as package bootstrap instead of referencing a bootstrap-local variable
- update the checksum-pinned vendored htslib dependency to 1.24 and use its native BCF/tabix multi-region iterators for indexed `read_bcf(...)`, `read_bcf_v2(...)`, `read_gff(...)`, `read_gtf(...)`, `read_tabix(...)`, and every `read_bcf_appender(...)` path. Overlapping or repeated region requests now emit each matching record once. The appender's `region_threads` partitions work only by primary contig and keeps every interval for that contig in one native iterator, so thread count no longer changes the target schema or row multiset; remove the thread-only `duckhts_region_idx` column and retain opt-in `FILE_OFFSET` as the file-order token. `read_bam(...)` already used the equivalent deduplicating `sam_itr_regarray(...)`; add a spanning-alignment regression across disjoint intervals. Keep BCF work claiming iterative and bounded to 16 workers, with a materializing regression over 256 leading header-only contigs plus appender checks across empty requested contigs. Build vendored htslib as GNU C17 so GCC 16 does not compile this pre-C23 dependency under its new GNU C23 default. Add SQL/R regressions and a Fedora GCC 16 R-devel check matching CRAN's `r-devel-linux-x86_64-fedora-gcc` environment
- emit Ensembl VEP 116's six regulatory-region and transcription-factor-binding-site consequences directly from the main small-variant and exact structural adapters. `duckvep_ensembl_regulation_features(...)` prepares a source-identified relation from release-matched funcgen tables and excludes EMAR rows exactly where VEP's source does; `duckvep_model_load(...)` keeps its five hot columns in a separate immutable SoA and cgranges seed index, while stable IDs and metadata remain relational. Transcript and feature sweeps advance together through the same generic sorted-interval helper and resumable output cursor; rich and compact rows identify the overlap object and its dense feature ordinal. The original eight-argument `duckvep_model_receipt(...)` call remains valid, with regulation_features_table as an optional named relation whose geometry is validated, counted, and hashed. The networked fixture recipe now pins and checksum-verifies exact Ensembl-116 core and funcgen manifests, folds them into a canonical source-manifest hash, and ships three real MotifFeature rows on the 1,143-base `KI270395.1` region for offline SQL/R tests. A brute-force all-pairs oracle checks one-row cursor resumption across 3,000 randomized small-allele scenes and 3,000 exact-SV scenes. Executable VEP 116 agrees on all 14,955 annotation-object pairs from 1,196 chromosome-21 GIAB sites and all 120,224 pairs from 2,700 generated exact structural events, covering every new SO term with no unresolved, extra, missing, or discordant row. The standalone kernel ABI is 0.14.0
- add dedicated paired-BND consequence functions `duckvep_annotate_breakend(...)` and `duckvep_annotate_breakend_compact(...)`. Both raw VCF loci enter one event; candidate discovery queries the resident cgranges index around both, while ordinary topology reproduces VEP 116's shifted local position and mate-aware `feature_truncation`. Results union VEP's endpoint-overlap consequence sets once per transcript; mate-only rows expose NULL rich region / zero compact region mask. The surrounding relation retains raw ALT, orientation, identity, and provenance. A seeded five-chromosome campaign covering all four bracket orientations generated 1,004 events and matched all 91,428 isolated executable-VEP transcript pairs with no disagreement, extra row, or missing row. The oracle uses one BND per VEP input buffer because VEP 116's chromosome-blind mate-coordinate interval tree otherwise lets neighboring breakends change the result. The new paired-event columns and explicit candidate-pair entry point make the standalone kernel ABI 0.13.0
- re-audit the two retained full-model GRCh38 ClinVar samples against the indexed VEP 116 cache after the cross-assembly consequence fixes: the coding-enriched sample now matches all 287,836 transcript pairs and the cross-chromosome sample all 316,397, with no unresolved, extra, missing, or discordant rows. Keep the previous nonzero runs in the append-only ledger so the rendered report shows that current evidence supersedes, rather than erases, the earlier frontier
- add typed exact structural-span annotation through `duckvep_annotate_sv(...)` and `duckvep_annotate_sv_compact(...)`: span operations use one-based inclusive coordinates, while INS uses `start = end = P` for the interbase site after reference base P; callers also supply DEL/DUP/TDUP/INV/INS/CNV/UNKNOWN operation and explicit copy direction, and contradictory metadata fails closed. Extend the pinned VEP-116 Sequence Ontology registry from the 34 transcript terms to all 41 registered terms, with regulatory-region and transcription-factor-binding-site facts evaluated by a separate interval-feature predicate producer rather than fake transcripts. Add seeded executable-VEP structural exploration across transcript/exon/intron/splice/UTR/CDS/start/stop states and allow independent campaigns to share one model database through read-only attachment. Eight GRCh38 chromosome/seed campaigns generated 40,375 events and matched all 2,140,911 executable-VEP transcript pairs. A separate held-out seed generated 100,000 random small events across SNV, MNV, insertion, deletion, and delins shapes; together with 268 fixed witnesses, all 100,268 transcript pairs matched executable VEP 116
- match the pinned VEP Plugins release/116 NMD implementation's separate feature geometry. Consequence and sequence changes still use the minimized edit, but NMD CDS coordinates and exon-end rules now use the complete VEP `VariationFeature`; retained bases in an equal-length allele can therefore change an early-CDS or penultimate-exon decision without changing the applied edit. Pure insertions reuse VEP's one-flank exon-edge projection and preserve the parent `TranscriptVariation`'s reversed CDS range, rather than the expanded allele range rendered in VEP JSON; the valid `1,0` range immediately before CDS base 1 remains defined and escapes by the plugin's early-CDS rule. Reuse the sequence classifier's cached `CDS end <= 101` fact only when the minimized edit and VEP feature span are identical. Pure-C tests pin both strands, exon edges, internal mapper gaps, and cached-versus-exhaustive behavior. The corpus differential executes checksum-pinned `NMD.pm` with a coordinate observer, retains the prediction confusion matrix, and renders it from the same append-only conformance ledger as SO terms; plugin-free and NMD-ineligible rows stay out of the NMD audit when one engine omitted the corresponding consequence row. The complete 48,967-allele ClinVar chromosome-21 run matched all 68,554 eligible NMD transcript predictions—including 29,416 mutually unresolved pairs—and all 1,331,664 core consequence sets
- avoid resolving the same exon-pair intron twice while evaluating VEP mismatch islands, and skip the shared intron predicate when neither the cached intron nor any splice/short-intron window can contribute. Keep that predicate inlined so the sorted point-classifier machine code remains unchanged; the checked-in benchmark reports the coding-SNV control separately as measurement noise rather than claiming an SNV speedup
- prepare each validated transcript's CDS-to-cDNA origin, end, starting exon, and phase once when opening an immutable DuckVEP model, so coding projection reuses those numeric relations instead of binary-searching both CDS endpoints for every variant/transcript pair. Topology-only transcripts whose CDS is not sequence-projectable retain the established exhaustive path. Keep one-time derived-model preparation out of the per-candidate instruction layout and align that hot loop to a 64-byte cache line on supporting compilers, preventing cold model changes from regressing the sorted-SNV lane. The extended standalone model ABI is version 0.12.0
- record the retained exact GRCh37 and *P. falciparum* VEP-cache artifacts in the checked-in DuckVEP conformance ledger, and render the newest tested source ancestor independently for each named corpus so a later run of one expensive corpus neither hides nor implicitly updates the evidence for another
- rebuild multi-edit DuckVEP CDS haplotypes in one linear reverse-coordinate pass over distinct reference and worker-scratch buffers, replacing one CDS-tail `memmove` per edit while preserving reverse-order coordinates, strand orientation, REF validation, exact input/output alias compatibility, and the independent randomized rebuild oracle. A second randomized gate proves that maximal differing islands, separate phased SNVs, and one equivalent uploaded MNV produce identical alternate CDS, peptides, and coding facts on both strands across start, body, and terminal codons
- match all 40,732 transcript consequence sets emitted by the indexed Ensembl Genomes VEP 63 cache for the retained 8,444-allele *P. falciparum* corpus, with every pair resolved and no disagreements, extra rows, or missing rows; the generic CDS-boundary fixes also closed the previous non-human frontier without organism-specific branches
- match all 486,464 transcript consequence sets emitted by the indexed VEP 116 cache for a deterministic 65,616-allele cross-chromosome GRCh37 sample, with every pair resolved and no disagreements, extra rows, or missing rows. This round reproduces VEP's independent start-lost/start-retained predicates; coordinate-only terminal-stop overlap; `cds_start_NF` frameshifts through a synthetic `X`; prepared-CDS partial terminal codons; the distinction between `coding_unknown&missense` and an `X`-suppressed synonymous call; exact-exon rechecking inside short frameshift introns; the inclusive short-intron splice/PPT edge; the unusual co-occurrence of retained-stop and protein-altering predicates; CDS-boundary insertion UTR exceptions; and supported generic coding results for equal-length features crossing the outer 5-prime transcript boundary
- expose VEP-style directional annotation windows on `duckvep_annotate(...)` and `duckvep_annotate_compact(...)`: five arguments retain the 5,000-base upstream/downstream default, a sixth value applies symmetrically, and a seventh value lets callers set upstream and downstream distances independently; zero disables the corresponding direction
- reproduce VEP 116's reverse-strand insertion-length terminal-stop fallback: an insertion whose added genomic span reaches the stop is reconstructed through the original CDS endpoint before deciding whether that stop is retained; materialized and borrowed coding contexts now preserve the same strand/orientation facts
- reproduce VEP 116's independent `coding_unknown` and `missense_variant` predicates for equal-length edits on `cds_start_NF` transcripts; a leading unknown codon may coexist with proven downstream amino-acid changes
- match VEP 116's mature-miRNA boundary predicate for pure insertions. VEP tests the minimized reversed `VariationFeature` interval `(P+1,P)`, rather than the transcript-relative placement point, so an insertion after the last mature-miRNA base now remains `non_coding_transcript_exon_variant`; the frozen GRCh37 indexed-cache differential removed all 3,667 false mature-miRNA calls exposed by this boundary case
- source DuckVEP codon tables from Ensembl `seq_region_attrib`, with VEP 116's table-1 default and the complete BioPerl-supported NCBI table set instead of chromosome-name inference. The Ensembl compiler now retains valid single-residue `initial_met`, `_selenocysteine`, `amino_acid_sub`, and `_stop_codon_rt` Translation SeqEdits as a sparse reference-peptide relation; the immutable model loader validates and packs that relation, and the consequence kernel applies it only to VEP's reference peptide while leaving alternate codon translation unchanged. Unsupported edit shapes and transcript RNA edits remain explicit sequence-withholding states. Attribute values are cast at the importer boundary, so official dumps whose numeric-only columns are inferred as integers remain valid inputs. Model receipts now hash and count mature-miRNA segments and peptide edits, the offline GRCh38/GRCh37 fixtures include their real sequence-region attributes and mitochondrial edit witnesses, and the extended standalone model ABI is version 0.11.0
- match every transcript consequence set emitted by the indexed VEP-116 cache for a deterministic 4,246-allele ClinVar chromosome-21 sample (seed 1, at most 1,000 alleles per shape/length bin): 126,320/126,320 exact pairs, with no unresolved, extra, or missing rows. The Ensembl builder now applies VEP's core-transcript selection before assigning model ordinals and projects each mature-miRNA cDNA attribute into genomic exon segments; `duckvep_model_load(...)` accepts those segments as an optional packed side relation, and the generated rule table gives `mature_miRNA_variant` VEP's precedence over generic non-coding exon terms. The same round fixes insertion placement at transcript/exon boundaries, VEP's empty-UTR endpoint predicate, and coding edits on `CDS_END_NF` transcripts with a one- or two-base partial terminal codon. Bump the standalone kernel ABI to 0.10.0. The cache differential remains a declared chromosome-21 workload rather than a claim of complete VEP behavior, while the corpus runner now records the cache/core/variation build provenance and the production throughput driver can load complete coverage, transcript flanks, and mature-miRNA segments
- add `duckvep_annotate_compact(...)`, the numeric materialization of the same sorted consequence kernel used by `duckvep_annotate(...)`, so transcript/gene ordinals, SO and region masks, status/reason codes, coordinates, amino-acid bytes, and NMD codes can be filtered or joined before text rendering. Cache DuckDB vector data/validity once per chunk, skip the coarse splice-region pass while applying the public per-call reaches directly to the exact generic VEP predicate, return default-window deep-intron states without rescanning irrelevant splice gaps, and generate a sparse fact-to-SO lookup from the existing rule ledger while retaining the rule interpreter as the test oracle. On the pinned one-thread Ensembl-116 GRCh38 topology/GIAB sample (644,292 transcripts, 5,078,384 exon memberships, 100,957 sorted alleles, 1,179,465 output rows), an already-ordered table scan measured 100 ms median across nine compact CLI passes (about 1.01 million input alleles/s and 11.79 million output rows/s). Restating `ORDER BY` at the timed SQL boundary measured 106 ms median in the CLI (about 0.95 million alleles/s) and 107 ms through the checked-in R/DBI driver (about 0.94 million alleles/s); the corresponding rich CLI projection measured 208 ms median. These are workload-specific materialization measurements, and the separate figures distinguish the sorted-kernel floor from SQL order enforcement rather than making a universal throughput claim
- make one edit/CDS/peptide context authoritative for length-changing small variants; production annotation no longer retries a failed context through a shape-specific classifier. Complete transcript-oriented sequence before and after the CDS is packed in one immutable cold byte pool, allowing the allocation-free evaluator to reproduce VEP 116's independent start/stop predicates, endpoint-mapper-Gap behavior, and transcript-associated default `intergenic_variant` without guessing from a short tail. The transcript loader accepts the original 11 columns, the legacy 12-column short-tail form, or a complete 13-column projection; older models fail closed as `missing_transcript_flank`. Bump the standalone kernel ABI to 0.9.0. The regenerated differential is exact for all 268 fixed pairs and all 401,072 pairs across four frozen 100,000-case distributions, with zero unresolved, missing, extra, or discordant rows
- match VEP 116 for equal-length uploaded features that cross from 5-prime UTR into the start codon. DuckVEP now applies the contiguous CDS slice through the shared coding context but evaluates only VEP's start predicates: changing the start emits `start_lost` without a co-occurring stop term, while preserving `ATG` emits `start_retained_variant` even when later coding bases in the same feature change. The 262-row fixed differential is exact for all 252 supported pairs. Pairwise comparison against the previous engine over three frozen 100,000-case distributions shows 3,016, 3,138, and 3,136 mismatch-to-match transitions, zero match-to-mismatch transitions, no remaining unresolved MNVs, and 14, 13, and 11 fewer supported disagreements
- match VEP 116 when retained bases in an equal-length uploaded feature widen the peptide-predicate window beyond its prefix/suffix-trimmed edit. Paired start- and terminal-stop witnesses now distinguish `start_lost&start_retained_variant`, `stop_retained_variant`, `missense_variant`, and `stop_gained` for records that apply the same one-base CDS change but carry different uploaded spans; retained REF mismatches and ambiguous widened windows fail closed instead of retrying the smaller edit. All 248 supported pairs in the expanded 258-row fixed VEP differential are exact, with 10 unsupported contexts remaining explicit. Three frozen 100,000-case distributions reduced resolved disagreements from 674 to 546, 669 to 550, and 682 to 578 through 351 mismatch-to-match transitions with no match-to-mismatch regression
- match VEP 116 when an equal-length uploaded feature crosses from CDS into the transcript-oriented 3-prime UTR. DuckVEP now preserves VEP's unavailable peptide-mapping state and emits `coding_sequence_variant` instead of translating a smaller prefix/suffix-trimmed edit; paired CDS-only and CDS-to-UTR witnesses on both strands prevent blanket stop-predicate suppression. The 252-row fixed differential is exact for all 242 supported pairs, and three independent 100,000-case distributions reduced supported disagreements from 742 to 661, 723 to 660, and 758 to 674
- match VEP 116 when a lengthening replacement's ALT-only mismatch island reaches an intron even though its REF-shaped feature remains exonic. DuckVEP now reproduces VEP's three-base interval-tree intron prefilter before applying the essential-site exclusion; paired fixed witnesses pin both sides of the cache boundary, and two independent 100,000-case held-out differentials reduced supported disagreements from 955 to 742 and from 933 to 723
- reproduce VEP 116's codon-local predicates for length-changing edits instead of classifying from a whole-protein diff or net length alone. `inframe_insertion` now sees the alternate peptide truncated at its first stop while `protein_altering_variant` sees the untruncated peptide, so a local `Y` to `*QCY` change is correctly `stop_gained` alone; shortening delins use VEP's separate codon-string deletion predicate and can be `protein_altering_variant` despite a minus-three-base change. Paired fixed witnesses and randomized delins properties prevent blanket term suppression or modulo-three classification
- reproduce VEP 116's terminal-stop insertion state machine instead of inferring the consequence from length modulo three. Insertions inside or immediately before the terminal codon can now emit VEP's unusual combinations of `inframe_insertion` with `coding_sequence_variant`, `stop_lost`, `stop_retained_variant`, or `stop_gained`; the local trailing-partial-codon `X` exception and predicate order are pinned in a compatibility errata ledger and executable VEP witnesses. The expanded fixed corpus is exact for all 239 supported pairs, with 11 further pairs remaining explicitly unresolved
- match VEP 116's three distinct allele geometries: preserve the uploaded VCF span, derive the parser's one-anchor-removed `VariationFeature` span for transcript candidate and splice predicates, and use the fully trimmed edit only for sequence application. Multi-base splice alleles now follow VEP's differing-region islands and predicate accumulation, start-overlapping length changes can emit `start_retained_variant` alongside their coding consequence, and terminal frame changes follow VEP's stop precedence. The standalone kernel ABI is 0.8.0
- add a deterministic large-state DuckVEP gate combining the independent pure-C property oracles with generated alleles concentrated around transcript boundaries and sampled throughout the transcript. Seeds, generator bounds, exact alleles, pair-level VEP 116 comparisons, unresolved reasons, and statistical summaries are retained per run; the R generator now uses preallocated typed columns so six-figure and larger campaigns do not spend their time growing one-row data frames. The rendered conformance report selects the nearest tested ancestor of the current commit, shows every corpus recorded there, and uses its largest corpus for detailed disagreement tables
- distinguish the existing VEP core `NMD_transcript_variant` biotype consequence from variant-induced nonsense-mediated-decay prediction. `duckvep_annotate(...)` now applies the pinned VEP Plugins release/116 NMD policy to stop-gained, frameshift, splice-donor, and splice-acceptor consequences, returning `triggering`, `escaping`, or explicit `unresolved` state plus independent intronless, early-CDS, last-exon, and penultimate-exon-end escape reasons. The standalone kernel ABI is 0.7.0; on the 64-bit release ABI the compact result remains 56 bytes by using former alignment padding
- `read_bcf()` and `read_bcf_v2()` now recognize the `Format=...` CSQ schema
  spelling used by Ensembl variation release VCFs, exposing their consequence
  fields through the same typed `VEP_*` columns as ordinary VEP output
- add a reproducible typed-Parquet measurement for official Ensembl consequence VCFs.
  Release-116 chr22 contains 14,920,904 records and 30,199,106 CSQ entries; Zstandard
  Parquet is 82.7% of the already BGZF-compressed source for all 51 reader columns and
  58.7% for the 14-column consequence projection
- retain the versioned RefSeq accession carried by each Ensembl MANE Select or MANE Plus Clinical attribute in the prepared transcript relation while keeping only its selection bit in the resident C model; reject empty or conflicting MANE mappings. Give `duckvep_model_receipt(...)` stable provenance column names, and add paired offline Ensembl-116 acceptance fixtures for GRCh38 and the archived GRCh37/GENCODE-19 model. The roughly 116 KiB fixtures exercise real MANE, ordinary coding, and mitochondrial fail-closed paths; explicit staging verifies both dump manifests, reference hashes, and deterministic model receipts while normal tests remain network-free
- document the current DuckVEP Ensembl-model build as a concept-to-code map: staged core relations plus FASTA, relation validation and sequence reconstruction, deterministic receipts, immutable model publication, failure behavior, ownership, and the exact boundary of work not yet implemented. Link the project architecture and DuckVEP contract from the generated README
- add `duckvep_ensembl_regions(...)`, `duckvep_ensembl_transcripts(...)`, and `duckvep_model_receipt(...)`: DuckDB now builds the resident consequence model directly from pinned Ensembl core relations plus tiled matching FASTA sequence, preserves stable/source identifiers and transcript attributes, verifies assembly/contig/exon/translation invariants, and emits a deterministic receipt over both region coverage and transcript content. Ensembl RNA/peptide edits not yet implemented by the kernel retain their flags and CDS span while withholding CDS bytes with an explicit reason, including `_rna_edit` records carried by either transcript or translation attributes. A complete local Ensembl-116/GRCh38-primary build validated 194 regions, 646,577 current transcripts, 5,087,789 exon memberships, and 370,580 sequence-backed coding transcripts through the resident model loader
- accept Ensembl exon phase `-1` for a sequence-backed CDS when translation starts after 5-prime UTR in the same exon. Ensembl and VEP treat that state as zero prepended phase bases; positive phases still prepend one or two `N` bases and the shared model validator still checks the resulting prepared-CDS length exactly
- match VEP 116's terminal-stop precedence for frame-changing edits: when the first affected reference codon is the terminal stop, DuckVEP now emits `stop_lost` or `stop_retained_variant` instead of `frameshift_variant`. The resident transcript query may add an optional twelfth `post_cds_bases` BLOB with up to three transcript-oriented reference bases; this costs three bytes per transcript and is sufficient to reconstruct the original translation endpoint for short terminal deletions. Models without the needed bases fail closed as `missing_transcript_tail`. Bump the standalone kernel ABI to 0.6.0; the real VEP differential is exact for all 210 resolved pairs, with 32 explicitly unresolved pairs and no supported disagreements
- match VEP 116 insertion placement at exon, CDS, and transcript boundaries. The consequence engine now projects a pure insertion through either genomic flank, so an intronic VCF padding base no longer hides a coding insertion at the next exon; splice predicates still use VEP's interbase interval. Real VEP witnesses now cover coding+splice, CDS-to-UTR, and transcript-to-downstream boundaries
- bump the standalone DuckVEP kernel ABI to 0.5.1 and preserve coding failure causes through the compact C result and SQL adapter: unresolved rows now distinguish missing or ambiguous sequence, reference mismatch, non-contiguous edits, unsupported compound effects, invalid projection, and internal capacity failures. Protein coordinates and scalar amino-acid values now have independent validity, so frameshifts and in-frame indels return NULL amino acids instead of empty strings
- make `duckvep_model_open(...)` the canonical fail-closed validator for both standalone and SQL-loaded resident models. It now rejects invalid transcript strand/span, exon order or envelope, discontinuous/non-length-preserving cDNA coordinates, invalid phase, non-exonic CDS boundaries, prepared-CDS length mismatches, unsupported codon tables, and non-ACGTN sequence before annotation can run
- make DuckVEP model coverage an explicit trust contract. Ordinary `duckvep_model_load(...)` calls are partial and now return an unresolved `sequence_variant` with reason `no_feature_in_loaded_model` when no loaded transcript is found. The named parameter `transcript_coverage_complete := TRUE` requires a `sequence_length` column in the sequence-region query, rejects out-of-contig spans, and only then permits supported `intergenic_variant` synthesis
- prepare each ordinary DuckVEP allele once per input row before transcript expansion. Event kind now follows the trimmed differing region, so padded one-base substitutions take the SNV path; an explicit interbase boundary and anchor side support VCF position-1 insertions/deletions whose mandatory padding base follows the event. The cursor retains its bounded-output contract, with C property, SQL, and R regressions covering padded SNVs and right-anchored indels
- consolidate the standalone DuckVEP kernel ABI as 0.5.0: remove ignored supplementary-track and HGVS arguments, placeholder consequence bytes, and impossible route states; keep the still-used direct coding fallback explicit and measured while the shared edit/CDS/peptide interpreter is completed
- accelerate sorted SNV annotation with one compact exon cursor per transcript and worker, retained across adjacent input chunks and reset on coordinate rewind. The point path shares the exact splice predicates with the span classifier and does not alter the existing MNV/indel/SV span classifiers or grouped-haplotype edit core. Sparse Sequence Ontology masks now avoid scanning all assigned terms during impact selection and single-term rendering
- make release and debug builds refresh extension metadata from the authoritative root `description.yml`, preventing a previously generated version file from leaking the prior release number into a new artifact
- emit VEP-116 `NMD_transcript_variant` alongside the ordinary regional consequence when a variant lies inside a transcript whose imported biotype is `nonsense_mediated_decay`; the upstream/downstream halo does not receive the term
- add `duckvep_model_load(...)`, `duckvep_model_drop(...)`, and `duckvep_annotate(...)`. Several named immutable transcript models can coexist; each is validated and narrowed once from three committed sorted DuckDB queries, keeps prepared coding sequences resident, and builds one cgranges transcript index. Annotation reads biallelic small-variant columns directly through the stable DuckDB extension API, reuses adapter and kernel storage, seeds each sorted run once, then advances through transcripts and exons. Missing sequence facts are returned as explicit unresolved rows instead of being counted as supported behavior
- port the pure-C consequence kernel and its testing framework from `duckvep-c` commit `9f922c8`. Fixed biological cases and randomized properties check interval sweeps, strand-aware projection, genetic codes, sequence edits, consequence composition, phased edit grouping, and bounded output against independent oracles; normal, AddressSanitizer, UndefinedBehaviorSanitizer, and 100,000-trial targets exercise rare states. Add deterministic boundary/codon witnesses and a large-VCF runner that compares the union of DuckVEP and VEP-116 `(variant, transcript)` rows, retains every mismatch/missing/extra/unresolved row, and reports exact binomial upper bounds by consequence and allele shape
- add append-only DuckVEP correctness, randomized-property, and throughput ledgers plus rendered R Markdown reports. Correctness rows come only from real VEP 116 runs, retain the exact Ensembl build and artifact hash, and split both by full consequence set and individual Sequence Ontology term. The property ledger records every theft target, seed, requested/passed/duplicate count, compiler, and suite result from a successful pure-C run. The stable-API throughput baseline consumes ten million sorted inputs and is explicitly labeled as a one-transcript hot-path measurement rather than a whole-genome claim
- make the vendored DuckVEP property runner build on MinGW by applying an explicit test-build patch that disables theft's optional POSIX fork/timeout path on Windows; the in-process property path used by DuckVEP remains identical across platforms, while requesting fork mode on Windows fails explicitly
- remove the discarded DuckVEP R/SQL model builder and replace its planning documents with the current source-data, resident-model, sorted-sweep, phased-edit, and statistical-test contract. Make extension metadata regeneration take its version from root `description.yml`, preventing release builds from retaining the stale generated version or falling back to a Git commit

# duckhts 1.4.0 (2026-07-10)
- implement `bam_nt16_counts` and `nt16_gc_counts` for AArch64 NEON and DuckDB-Wasm SIMD128, eliminating their scalar fallback on those targets while preserving the scalar contracts for odd packed tails, IUPAC buckets, invalid unpacked codes, and empty inputs. DuckDB-Wasm enables SIMD128 only for its backend translation unit; its real-browser test now requires concrete SIMD128 selection and compares both long unpacked nt16 input and a no-index BAM scan against scalar. Emscripten `bam_bin_counts(...)` disables unavailable htslib worker threads for that synchronous scan path. Replace the broken x86-only kernel microbenchmark with a normally linked backend contract runner covering deterministic boundary, tail, alignment, invalid-code, odd-padding, and randomized cases before benchmarking both nt16 kernels; remove signed-overflow-prone AVX2 loop bounds; and prevent retained Wasm builds from reusing a host curl-enabled htslib `config.status`
- add the private DuckVEP annotation-model foundation: a normalized 26-relation DuckDB/Parquet schema with an R validator that distinguishes structural, model-candidate, and fail-closed conformance checks; wide half-open edit primitives with semantic trimming and explicitly tagged VEP compatibility runs; deterministic 12,000-case property tests; and complete extension/R-package build wiring. No public SQL annotation surface is introduced yet
- overload the `cigar_*` family (`cigar_has_soft_clip`, `cigar_has_hard_clip`, `cigar_left_soft_clip`, `cigar_right_soft_clip`, `cigar_query_length`, `cigar_aligned_query_length`, `cigar_reference_length`, `cigar_has_op`) to accept a `UINTEGER[]` binary CIGAR — the packed `oplen<<4|op` list `read_bam(cigar_representation := 'binary')` emits — and overload `seq_hash_2bit(...)` to accept a `UTINYINT[]` of htslib nt16 codes, both via DuckDB scalar-function sets that resolve by argument type. The binary/nt16 paths decode the packed values directly (BAM op order `M=0,I=1,D=2,N=3,S=4,H=5,P=6,==7,X=8`) and are bit-identical to the text path (verified 0 mismatches across all cigar metrics and `seq_hash_2bit` over 3.2M reads, plus empty/`*`/non-ACGT edge cases). This completes the decode-free `seq_*`/`cigar_*` binary-overload family so BAM pipelines analyze projection-pushed binary columns without a text round-trip. `seq_kmers(...)` stays text-only (DuckDB table functions cannot overload by argument type); nt16 callers compose `seq_kmers(seq_decode_4bit(SEQ), k)`
- add a `nt16_gc_counts` logical SIMD kernel and route the nt16 `seq_gc_content(UTINYINT[])` overload through the capability-mask dispatch framework (snapshotting the immutable dispatch table once per chunk). Each row's contiguous nt16 codes are classified by the dispatched backend; a branchless scalar reference is the correctness oracle (results unchanged, byte-identical to the text path), and an AVX2 backend lands alongside it (`cmpeq` against A/C/G/T/N, bit-identical to scalar). Measured on an in-memory nt16 column (single-threaded, ~3.3 GB of codes): the classify runs ~7x faster on AVX2 than the scalar kernel (~1.5 -> ~10.6 GB/s); `duckhts_simd_kernel_info()` now reports the `nt16_gc_counts` kernel. The kernel runs for the common case (no NULL bases, i.e. a NULL child validity mask); rows whose list child carries a validity mask fall back to the scalar per-base path so a NULL base is treated as invalid. Follows the SIMD workflow (scalar oracle, backend-agnostic scalar-vs-auto SQL/R tests)
- optimize the nt16 `seq_gc_content(UTINYINT[])` overload by hoisting the per-base list-child validity probe out of the inner loop: the child validity mask is loop-invariant and BAM `SEQ` has no NULL bases, so the common case now folds the branch away instead of calling `duckdb_vector_get_validity()` per base. Results are byte-identical (0 mismatches over 3.2M reads); the gc_content-attributable time drops ~3.2x on an exome BAM (~540 -> ~1700 MB/s of nt16 codes, single-threaded). Also the prerequisite for a vectorized nt16 classify kernel (tracked as a follow-up); the same hoist applies to the `seq_revcomp`/`seq_canonical` nt16 paths once those land
- overload `seq_revcomp(...)` and `seq_canonical(...)` to accept BAM nt16 sequences: alongside the existing `VARCHAR`->`VARCHAR` signatures they now take a `UTINYINT[]` of htslib nt16 codes (from `read_bam(sequence_encoding := 'nt16')`) and return `UTINYINT[]`, via DuckDB scalar-function sets. The nt16 paths reverse-complement / canonicalize the codes directly (comparing by decoded base order for `seq_canonical`, since nt16 code order differs from base order at `N`) and reject any non-ACGTN code as NULL, exactly like the text path; both are bit-identical to the text result after decoding (verified 0 mismatches over 3.2M reads plus IUPAC/empty/rev-wins edge cases). BAM pipelines can now reverse-complement and canonicalize without leaving the nt16 encoding. Continues the `seq_*`/`cigar_*` binary-overload family
- overload `seq_gc_content(...)` to accept BAM nt16 sequences: alongside the existing `VARCHAR` signature it now takes a `UTINYINT[]` of htslib nt16 codes as produced by `read_bam(sequence_encoding := 'nt16')`, registered via a DuckDB scalar-function set so DuckDB resolves by argument type. The nt16 overload classifies codes directly against htslib's fixed `A=1,C=2,G=4,T=8` encoding (`called` ⇔ A/C/G/T, `gc` ⇔ C/G) and is bit-identical to the text path (verified 0 mismatches over 3.2M reads), so projection-pushed BAM pipelines compute GC without decoding sequences back to text. Scalar reference only for now; a nt16 SIMD backend and further `seq_*`/`cigar_*` binary overloads are tracked as follow-ups
- add `benchmarks/benchmark_riker_wgs.Rmd` (rendered `.md` + ggplot figures) plus the harness `scripts/stage_riker_wgs_bam.sh` (downloads + transcodes the input BAM) and `scripts/riker_wgs_benchmark.py` (cold-cache, core-pinned sweep), a harmonized three-way WGS coverage comparison of `riker 0.4.0` `wgs` vs native `duckhts_mosdepth(...)` vs upstream `mosdepth 0.3.13` on a 30x 1000 Genomes BAM (HG00188/ERR3240174), following the [fulcrumgenomics/riker benchmark pipeline](https://github.com/fulcrumgenomics/riker/tree/main/benchmark-pipeline) methodology (identical read selection, cold page cache per run, physical-core pinning, 3 reps). Records that `duckhts_mosdepth` is byte-identical to `mosdepth` under the harmonized flags, that riker wins single-threaded (one sequential SIMD pass, `u16` capped-250 depth) while `duckhts_mosdepth` scales to 1.70x faster at 4 threads via contig-parallel workers, and documents the decompressor (libdeflate vs zlib) and coverage-cap asymmetries plus the whole-contig `int32` memory that tiling (`design/coverage_memory_footprint.md`) would bound
- add `decode_error_policy := 'null'|'warn'|'error'` to `read_bcf(...)`, `read_bcf_v2(...)`, and `read_bcf_appender(...)` for corrupt BCF header-vs-payload FORMAT/INFO type clashes; the default `null` policy materializes NULLs, `warn` emits a DuckHTS warning and materializes NULLs, and `error` raises a DuckDB error. FORMAT and INFO decode now preflight BCF encoded types before htslib helpers under every policy (including the default `null`), so corrupt inputs neither terminate the process with `exit(1)` nor leak raw bytes; in particular an INFO field whose header claims `Type=String` over a numeric payload previously returned truncated garbage under the default policy because htslib's string path copies the payload without checking the encoded type, and now materializes NULL. Add corrupt-BCF regression fixtures (including a reverse String-header/numeric-payload clash) plus SQL/R coverage while keeping valid mixed-ploidy `Number=G` FORMAT records accepted
- add a headless Playwright smoke test for the duckdb-wasm extension build (`make wasm-playwright-test` locally and a `wasm-playwright` GitHub Actions job on pull requests), reusing the existing `scripts/start_duckdb_wasm_local_test.sh` Docker build path (now split into build-only via `SERVE=0`); the test loads the wasm extension in a real browser and asserts the SIMD dispatch kernels resolve on wasm (e.g. `bam_nt16_counts` falling back to scalar) plus an end-to-end `seq_gc_content(...)` check. The local Docker harness now excludes/removes stale host `configure/venv` copies before rebuilding inside the container and launches Chromium with `--no-sandbox` for root-run dev/CI containers. This supports a one-PR-per-optimization cadence where the full CI matrix (build pipeline + wasm builds + R CMD check + browser smoke test) validates each SIMD/reader optimization
- add a `bam_nt16_counts` logical SIMD kernel and route `bam_bin_counts(...)` GC base counting through the capability-mask dispatch framework (snapshotting the immutable dispatch table once per bind, shared read-only across worker threads); base classification follows htslib BAM nt16 semantics (`A=1,C=2,G=4,T=8,N=15`, reconciling to `seq_nt16_int` so `called` is exactly A+C+G+T as samtools/htslib would count). A branchless scalar backend is the correctness oracle (GC results unchanged) and an AVX2 backend lands alongside it (bit-identical to scalar; ~10x kernel throughput vs the scalar reference on an AVX2 host in the standalone microbenchmark). Add `make bench-simd-kernels` (`test/scripts/bench_simd_kernels.c`), a DuckDB/R-free C microbenchmark that gates every compiled backend against the scalar oracle and reports scalar-relative throughput per logical kernel. See `design/packed_state_kernels.md`
- surface malformed-record BCF/VCF scan failures from `read_bcf(...)`, projected `read_bcf_v2(...)` scans, and the non-parallel `read_bcf_appender(...)` path as DuckDB errors instead of treating htslib parse/read failures as EOF; run `read_bcf_appender(...)` writes in an internal transaction so malformed input rolls back target-table side effects; add a malformed-POS fixture plus SQL/R regression coverage
- restore the experimental `read_bcf_v2(...)` table function from recovered local work, preserving `read_bcf(...)` schema compatibility while adding sample pushdown, INFO/FORMAT/VEP field filters, projection-aware VCF unpacking, persistent decode caches, and count-only shortcuts for benchmarking; document the DuckVEP/`bcftools csq` implementation roadmap, `read_bcf` vs rendered/recovered `read_bcf_v2` optimization evidence, and TigerStyle/runtime-safety cleanup backlog under `design/`; add DuckVEP fixture plumbing for BCSQ parsing plus DuckDB exact/range annotation joins, add SQL/R coverage for malformed `additional_csq_column_types` errors, and keep generated R extension C sources synchronized from top-level `src/`
- add the initial experimental `read_bcf_appender(...)` benchmark helper, which read BCF/VCF and appended a narrow tidy table through DuckDB's chunk appender API for comparison against `read_bcf(...)` query pipelines; at introduction, the single-stream path ran inside an internal transaction while the `region_threads > 1` path was explicitly best-effort and lacked rollback (the current development entry above records its later transactional, thread-invariant contract)
- add the experimental `hts_region_union_query(reader, path, regions, params := '')` scalar macro, which generates a `UNION ALL BY NAME` query string reading one HTS file through separate per-region table-function scans (adding `filename`, `duckhts_region_shard_id`, and `duckhts_region_shard` columns); the generated query does not deduplicate records, so overlapping shards can duplicate boundary-spanning rows
- fix `read_bcf_appender(...)` so a region query that matches no records (absent contig or out-of-range coordinates) still creates and commits the empty target table instead of jumping to cleanup with the internal transaction open, which previously rolled back the `DROP`/`CREATE TABLE` and could leave a stale prior table under `overwrite := true` or never create the advertised table; add SQL coverage for the empty-region and overwrite-replacement cases
- harden the bcftools filter try/recovery shim so a `try_begin` env-stack overflow unwinds into the deepest active recovery frame instead of returning a sentinel that callers would mistake for a longjmp and then double-pop, which previously could unbalance the recovery stack on nested-at-max-depth use; extend the shim error-path probe with overflow-unwind coverage
- update the vendored DuckDB C API headers to DuckDB v1.5.3 while keeping the stable extension ABI metadata at v1.2.0; the refreshed headers include `DUCKDB_TYPE_GEOMETRY` and `DUCKDB_TYPE_VARIANT`
- add runtime DuckDB type-support probe macros `duckhts_duckdb_type_supported(...)`, `duckhts_duckdb_supports_variant()`, and `duckhts_duckdb_supports_geometry()` for version-safe SQL feature gating
- simplify the generated `duckhts_bcftools_norm(...)` site-preserving table-macro query shape by removing the extra correlated scalar `LATERAL` subquery around `bcftools_norm_row(...)`, eliminating the site-preserving `LEFT_DELIM_JOIN` plan overhead; preserve caller columns even when their names collide with DuckHTS helper-column names used internally by earlier macro forms, add regression coverage for DuckDB's suffixed behavior when callers already have normalized-output column names, and tune normalization benchmark/conformance drivers to use explicit full-file `read_bcf(..., scan_mode := 'sequential')` streaming by default
- add explicit `scan_mode := 'auto'|'sequential'` controls to the primary reader table functions (`read_bcf`, `read_bam`, `read_fasta`, `read_fastq`, `read_bed`, `read_gff`, `read_gtf`, and `read_tabix`), allowing callers to force full-file streaming/counting instead of index-backed count or parallel scan paths where applicable; sequential mode is rejected for region queries, which remain index-backed
- optimize `bcftools_norm_row(...)` for already-normalized plain ACGTN allele rows by skipping kstring left-realignment setup when right/left trim predicates prove the row is unchanged; avoid per-row FASTA path duplication after the vector-local cache is established, reuse larger bounded per-thread reference windows, and document/defensively serialize htslib FASTA fetches while keeping the normalization reference cache thread-local to avoid the faidx cache race class fixed in https://github.com/RGenomicsETL/duckhts/issues/17 / https://github.com/RGenomicsETL/duckhts/pull/18
- make `duckhts_bcftools_norm(...)` gVCF-aware for vt/vcfnorm-style row normalization: `<NON_REF>` and `<*>` reference-block alleles now pass through with `GVCFReferenceBlock`, and mixed real-plus-gVCF-symbolic alleles normalize the real alleles while preserving symbolic alleles and caller-supplied reference-block `END` in site-preserving output; mixed `*` plus real alleles now follows the same ignored-symbolic path, while `*`-only rows remain `SpanningDeletion`; add phased GT/PL/GP/DS/PS FORMAT fixtures, including haploid/triploid/tetraploid `Number=G` cardinality cases, and reader regression coverage that pins phase-separator preservation and documents that row normalization/liftover macros are not full-record genotype-remapping VCF rewrite pipelines; add a non-render-time-networked `stage-norm-1000g-dragen-gvcf` helper and optional benchmark section for a 10 Mb 1000 Genomes DRAGEN HG00096 chr22 gVCF slice
- add extension-owned Parquet conversion SQL builders `duckhts_bcf_convert_parquet_sql(...)`, `duckhts_bam_convert_parquet_sql(...)`, `duckhts_gff_convert_parquet_sql(...)`, and `duckhts_tabix_convert_parquet_sql(...)`, plus thin R/DBI conversion wrappers that execute the generated `COPY ... TO ... (FORMAT PARQUET)` statements with DuckHTS write-format metadata, preserved raw headers (`vcf_header`, `sam_header`, `gff_header`, `tabix_header`), optional corrected header text, SQL-filter provenance, selected-column/partition metadata, arbitrary user metadata via `metadata := map(...)`/R named lists, optional caller-managed JSON-file metadata when DuckDB's `json` extension is available, and partitioned-output support for DuckLake-style registration of premade Parquet files
- include the final VCF `#CHROM`/sample header line in `read_hts_header(..., mode := 'raw')`, so VCF/BCF header metadata is sufficient for future spec-compliant VCF/BCF regeneration

# duckhts 1.3.0 (2026-05-29)
- add a real-BAM SIMD GC-content benchmark (`scripts/benchmark_simd_bam_gc.py`, `benchmarks/benchmark_simd_bam_gc.Rmd`, and `make bench-simd-bam-gc`) that compares scalar/auto/concrete backend requests for both end-to-end `read_bam(...)` scans and materialized BAM `SEQ` strings
- replace the initial backend-wide SIMD ops table with a NumKong-style capability-mask dispatch design: logical kernels are generated from `duckhts_simd_kernels.def`, concrete backends register only the kernels they implement, initialization builds immutable dispatch tables for `auto`/`scalar`/concrete backend policies, runtime selection swaps an atomic table pointer, `auto` resolves the best implementation per kernel, and new `duckhts_simd_kernel_info()` diagnostics report the selected backend/capability for each logical kernel; add a design note plus SQL/R conformance coverage and benchmark reporting for the selected `seq_base_counts` kernel backend
- harden SIMD dispatch follow-up behavior: initialize the SIMD ops table with an atomic one-time state machine so an in-flight initialization cannot overwrite an explicit backend selection, keep per-kernel scalar fallback dispatch for future ops-table growth, clarify `selectable` versus `available` diagnostics, and render the generated SIMD backend catalog call with Markdown code ticks so ASCII SQL quotes survive README/function-catalog generation
- remove htslib autoconf `HAVE_*` macro guards from all SIMD backend translation units (`duckhts_simd_avx2.c`, `duckhts_simd_avx512.c`, `duckhts_simd_neon.c`, `duckhts_simd_wasm_simd128.c`, `duckhts_simd_dispatch.c`): the design intent is always-compile all ISA backends and dispatch at runtime via CPU feature detection, so static compile-time gating on autoconf-generated macros was both wrong by design and fragile across out-of-tree builds; the compile-time guard is now `defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(DUCKDB_WASM_EXTENSION)` for x86 backends and `defined(__aarch64__)` / `defined(__wasm_simd128__)` for NEON/wasm, which are universally available from the compiler without autoconf; `DUCKHTS_POPCOUNT32` now unconditionally expands to `__builtin_popcount` on x86; add scalar-vs-auto backend SQL correctness tests covering GC=0/0.5/1.0 sequences, embedded-N calling, and soft-masked lowercase bases

# duckhts 1.2.1 (2026-05-07)
- add an eager runtime SIMD dispatch scaffold for byte-oriented DuckHTS kernels, organized as a small dispatch core plus per-backend sources under `src/simd/`; expose backend-agnostic diagnostics and explicit selection through `duckhts_simd_backend()`, `duckhts_simd_requested_backend()`, `duckhts_simd_backend_available(...)`, and `duckhts_simd_set_backend(...)`; route `seq_gc_content(...)` through the new scalar/optional-AVX2 base-count vtable; add AVX-512, ARM NEON, and wasm SIMD128 backend translation units with runtime-gated availability; keep the R package bootstrap/manual rebuild paths wired to the SIMD sources; document the scalar/auto SIMD flow in the README; and add `scripts/benchmark_simd_seq_gc.py` plus `benchmarks/benchmark_simd_seq_gc.Rmd`/`make bench-simd-seq-gc` for scalar-vs-auto process-isolated benchmarking that skips unavailable platform-specific backends on ARM/wasm/scalar-only builds
- fix `bcftools_liftover(...)` / `duckdb_liftover(...)` source/destination FASTA contig alias resolution during reference validation and sequence fetches, and align spanning-deletion `*` allele handling with upstream `bcftools +liftover`: inputs such as `23`, `24`, `26`, `X`, `Y`, `MT`, and `chr*` aliases now follow the same canonicalization path instead of validating `23` only against `23`/`chr23`, removing spurious `SourceRefMismatch` rejects for X/Y/MT indels when the chain uses canonical contig names but the FASTA uses `chrX` / `chrY` / `chrM`; `*`-allele rows now follow upstream non-symbolic swap/ref-add semantics instead of being short-circuited as symbolic, full-file HG001/HG002/HG006 GIAB conformance is now exact against installed `bcftools +liftover`, and new SQL/R regression coverage pins the `23 -> chrX`, `SWAP=2`, and `SWAP=-1` cases
- vendor the official VariantKey / RegionKey C API (Nicola Asuni, 2018; <https://doi.org/10.1101/473744>) and add new scalar SQL helpers for both the high-level bcftools-style and raw upstream numeric surfaces: `variantkey(...)` now mirrors bcftools `%VKX` / `+add-variantkey` on 1-based VCF rows while preserving the official hashed nonreversible mode for large, ambiguous, and symbolic alleles, `regionkey(...)` adds official 0-based half-open span keys plus overlap helpers, raw encode/extract/decode functions are exposed for the underlying unsigned integer fields, SQL/R regression coverage now pins reversible, hashed, symbolic, and overlap cases, and new real-callset conformance tooling (`scripts/variantkey_conformance.sh`, `benchmarks/benchmark_variantkey_conformance.Rmd`) validates full downloaded GIAB whole-WGS VCFs against installed bcftools `%VKX`
- fix `duckhts_bcftools_norm(..., split_multiallelic := true)` so ref-only and empty-ALT rows are preserved instead of disappearing when the ALT side normalizes to an empty list: `ALT='.'`, NULL ALT values, and empty ALT lists now survive split-mode output with row-preserving `RefOnly` / `NullInput` status, SQL/R regression coverage pins the edge cases, `read_bcf(...)` now exposes `decompression_threads := 0` for explicit htslib worker-thread control on bgzipped VCF/BCF scans, the rendered normalization benchmark now uses fairer `bcftools norm -Ou` timing plus full-file local GIAB conformance cases, and the root/package READMEs now include concrete normalization examples
- return resolved default sidecar/output paths from helper writers when callers omit them: `fasta_index(...)` now reports the generated `.fai` path instead of an empty string, and regression coverage now also pins the default-path return behavior for FASTA indexing, BGZF compression/decompression, and tabix/BAM/BCF index builders
- add `bcftools_norm_row(...)` plus the table macro `duckhts_bcftools_norm(...)` for bcftools/vt-style FASTA-backed allele normalization over SQL tables or derived queries: the new surface preserves original columns, appends `pos_normed`, `end_pos_normed`, `ref_normed`, `alt_normed`, `normed`, and `norm_status`, accepts ALT as either comma-delimited `VARCHAR` or `VARCHAR[]`, optionally splits multiallelic rows before normalization, and handles symbolic `<DEL>`/`<DUP>` rows when `end_pos_col` / `svlen_col` are supplied; add SQL/R regression coverage and real-upstream conformance hooks for GIAB/gnomAD-style callsets
- fix `bcftools_liftover(...)` / `duckdb_liftover(...)` indel parity in two exact upstream rewrite points: (1) the source-side repeat-run extension path now keeps extending across the cached source-reference window boundary instead of stopping early, and (2) the clip-pad `Needleman-Wunsch` path now keeps the best shift even when all candidate alignment scores are negative, matching upstream `clip_pad()` / `nw()` instead of silently leaving the padded interval unshifted; together these restore upstream-compatible `REF/ALT` and `swap` results for the remaining repeat-heavy indel mismatches, add focused SQL/R regression fixtures, and bring real-callset conformance to exact parity against installed `bcftools +liftover` on GIAB HG001 chr20 plus the full HG006 GRCh37 benchmark VCF
- reject `bcftools_liftover(...)` / `duckdb_liftover(...)` rows whose source `REF` cannot be validated during the upstream indel/difficult-SNP source-fasta path: instead of fabricating padded lifted alleles or aborting the query, the scalar/table liftover surfaces now preserve the row schema with `mapped = false` and `reject_reason = 'SourceRefMismatch'`; update SQL/R regression coverage and README examples accordingly
- refresh liftover benchmarking/conformance tooling to use the installed `bcftools +liftover` plugin directly instead of an unpacked zip artifact, and add `scripts/liftover_conformance.sh` for side-by-side DuckHTS vs upstream comparisons on synthetic inputs or real callsets such as sliced GIAB HG001/HG006 VCFs; the conformance path now queries lifted BCF output directly instead of assuming the upstream output can always be indexed/sorted after liftover. Add a curated real-case table at `scripts/conformance_case_table.tsv` plus `scripts/liftover_conformance_batch.sh` so the same conformance workflow can be run reproducibly across multiple GIAB benchmark VCFs now and future known-sites/site-only VCF sets (for example gnomAD) later
- add `read_pileup(...)` for region-scoped BAM pileups backed by htslib's pileup engine, returning per-position `chrom`, `pos`, `depth`, `bases`, and `quals` after SAM-flag and MAPQ filtering; also add `read_bam(..., cigar_representation := 'binary')` to expose packed BAM CIGAR ops as `UINTEGER[]`, and add explicit `gzi_path := '...gzi'` support to `read_fasta(...)` and `fasta_nuc(...)` for bgzipped FASTA sidecars that are not colocated with the FASTA
- speed up `fasta_nuc(...)` nucleotide counting on capable x86_64 hosts with an AVX2+popcnt fast path selected via htslib-style runtime dispatch, while keeping scalar behavior as the fallback on non-AVX2 or non-x86 builds
- apply principled remote-I/O tuning across long-running HTS scan paths and `bam_index(...)`: native remote reads/index builds now raise htslib block/cache sizes for streaming vs indexed-region access patterns, while Emscripten/browser builds apply the same policy with smaller wasm-safe budgets; add a minimal vendored htslib patch exposing a pre-opened `sam_index_build4(...)` entry point so DuckHTS can tune `bam_index(...)` before the remote scan starts
- fix `read_bcf` scanner crash on certain VCF/BCF files where `FILTER` list elements were written without reserving list-vector capacity (`duckdb_list_vector_reserve` was missing before `duckdb_list_vector_set_size`), which could corrupt allocator state and abort with `double free`/`invalid pointer`; `FILTER` now reserves space before writing list entries and is stable for full-table scans including PASS and explicit filter tags
- compile DuckHTS extension sources with `-Wpedantic` in the CMake and R-package build paths while keeping vendored `htslib` on its own upstream warning flags
- keep the non-Emscripten `wasm_http_hfile.c` translation unit non-empty so native pedantic builds do not warn about an empty source file
- harden Windows R-package `libcurl` probing so htslib remote URL support is enabled only when `curl_easy_init` links successfully with the detected `pkg-config` dependency closure; otherwise the configure path disables libcurl/S3/GCS cleanly instead of passing a false-positive `-lcurl`

# duckhts 1.2.0 (2026-05-07)
- add richer parsed attribute outputs for `read_gff(...)` and `read_gtf(...)`: `attributes_list := true` returns `MAP(VARCHAR, VARCHAR[])` with grouped multi-values and GFF3 percent-decoding, while `attributes_pairs := true` returns `LIST<STRUCT(key VARCHAR, value VARCHAR, idx INTEGER)>` for exact key/value/index records; keep `attributes_map := true` as the backward-compatible raw scalar map and make it skip duplicate keys instead of constructing invalid DuckDB maps
- add `read_gff(..., strict := true)` for GFF3 structural validation during scans: strict mode rejects malformed rows with the wrong number of fields, invalid coordinates (including nonpositive concrete end coordinates), invalid strand/phase/CDS phase, invalid scores, empty seqids/features, whitespace-containing feature types, and malformed attribute segments; strict count-only scans validate records instead of using index metadata, while default `read_gff(...)` remains permissive for existing ingestion workflows
- add a reproducible GFF3 conformance/performance benchmark against GFFBase v0.1.0 (`scripts/gffbase_conformance_benchmark.py`, `benchmarks/benchmark_gffbase_conformance.Rmd`, `make bench-gffbase`), auditing GFFBase Rust-vs-Python parser fallback/parity and reporting DuckHTS default-vs-strict validation plus attribute-semantics gaps alongside direct scan/parser throughput, with optional `GFFBASE_HUMAN_GFF` real human-scale inputs
- extend `bcftools_score(...)` to accept multiple TSV/SSF summary files as a DuckDB LIST/array or via `summaries_list_file`, producing one PRS column per summary file in a single genotype scan; add `log_path` audit output with loaded, matched, allele-mismatch, and duplicate-marker counts per PRS; make `summaries_list_file` directory scans deterministic and sidecar-safe; validate generated score/count column names for uniqueness; align score accumulation with upstream `bcftools +score` float32 summation; and update score conformance/benchmark workflows to use the bundled `bcftools` and `score` plugin from `RGenomicsETL/RBCFTools`
- collapse the generated top-level README function catalog behind a disclosure widget so the README opens on the overview and examples instead of a long function table
- refresh release documentation for the community-extension path: generated README SQL examples now load `build/release/duckhts.duckdb_extension` by relative path instead of embedding a machine-local absolute build path
- add `duckhts_cgranges_overlaps_list(...)`, a vectorized scalar overlap expander that returns LIST-of-STRUCT hit records so provider rows from `read_bed(...)`, `read_bam(...)`, `read_bcf(...)`, or tables can be expanded with `UNNEST(...)` without generated bulk-probe SQL; SQL/R tests cover one-row-per-hit streaming expansion, and the existing `duckhts_cgranges_overlaps_bulk(...)` probe path now also handles DuckDB string vector lengths safely
- fix `duckhts_cgranges_from_query(...)` ingestion of DuckDB string vectors by respecting string lengths instead of assuming NUL-terminated buffers; this fixes bulk cgranges construction from providers such as `read_bed(...)` with longer chromosome names and adds SQL/R regression coverage
- rework the cgranges benchmark scripts and rendered benchmark report around the new streaming-provider scalar probe path: `scripts/cgranges_benchmark.py`, `scripts/cgranges_benchmark_real.py`, and `scripts/cgranges_benchmark_cli.sh` now avoid generated bulk-query SQL / one-query-per-probe loops, compare DuckHTS scalar filtering/counting against `bedtk flt`, include `bedtools intersect -u/-c` when available, and render refreshed benchmark results under `benchmarks/benchmark_cgranges.md`
- add vectorized scalar cgranges probe helpers `duckhts_cgranges_has_overlap(...)` and `duckhts_cgranges_count_overlaps(...)`, so rows from streaming providers such as `read_bed(...)`, `read_bam(...)`, or `read_bcf(...)` can be filtered/count-annotated directly against an already-finalized session cgranges index without routing through the materializing `overlaps_bulk` query-string path; add SQL/R coverage for overlap, contain, and NULL probe semantics
- add `duckhts_cgranges_overlaps_bulk(...)`: streams a SQL query of overlap probes through a finalized cgranges index in one table-function invocation, supports `mode := 'overlap'|'contain'`, accepts an optional `query_row_id_col`, and defaults `query_row_id` to the 1-based probe row ordinal when no id column is supplied; add SQL/R coverage for the new bulk path and its generated row ids
- expand cgranges conformance and documentation: add typed/null/contain/error SQL coverage for `duckhts_cgranges_*`, add a `bedtk`-based overlap-existence parity script under `test/scripts/cgranges_bedtk_conformance.py`, publish the cgranges entry points in `functions.yaml`, and add root/package README examples for row-wise and `from_query(...)` construction
- fix `fasta_nuc(...)` GC/AT percentage denominators for intervals containing `N`: `pct_gc` and `pct_at` now divide by informative `A/C/G/T` bases only instead of total sequence length, so ambiguous bases no longer dilute bin/interval composition percentages; add SQL and R regression coverage
- add C-built bulk cgranges population through `duckhts_cgranges_from_query(...)`, executing the source query on an extension-owned DuckDB connection and building the cgranges index in C before publishing it to the session registry; keep `duckhts_cgranges_from_table(...)` deferred for now and add SQL coverage for the new query-ingest path plus overlap queries
- vendor htslib 1.23.1 for its upstream CRAM decoder and GZI validation security fixes, which matter for DuckHTS's native and wasm/browser-exposed parsing paths
- add native `duckhts_bam_bed_coverage(...)`: computes samtools coverage-like per-BED-region summaries for indexed BAM/CRAM input, returning one row per target interval with pre/post-filter `numreads`, `covbases`, `coverage`, `meandepth`, `meanbaseq`, and `meanmapq`, plus read-mode strand-specific post-filter summaries (`numreads_fwd/rev_post`, `covbases_fwd/rev_post`, `coverage_fwd/rev_post`, `meandepth_fwd/rev_post`); add SQL/R coverage tests on bundled fixtures and reserve `fragment_mode` / `processing_threads` for later phases
- reduce `duckhts_bam_bed_coverage(...)` peak memory by allocating and freeing per-base working depth buffers during scan processing instead of holding them for the entire BED load, and tile large target intervals in fixed 1 Mb windows when computing covered-base breadth; region scans are now single-pass again after the tiled rewrite, `min_depth > 1` mean-depth behavior is pinned to match `samtools coverage`, strand-specific depth scratch is skipped when `strand_outputs := FALSE`, and a new `decompression_threads` parameter makes htslib decode worker counts explicit for BAM/CRAM input
- add native `duckhts_samtools_idxstats(...)`: writes samtools idxstats-compatible TAB-delimited output for BAM/CRAM/SAM inputs, uses BAM index metadata (`hts_idx_get_stat(...)`) on indexed BAM fast paths, falls back to a streaming scan for CRAM/SAM/unindexed BAM, supports explicit `output`, `index_path`, `threads`, and `overwrite`, and adds SQL/R conformance coverage for BAM fast-path, CRAM fallback, custom index paths, and overwrite errors
- clean repo hygiene for local development: fix stale ignore rules that still referenced the old `duckdb_extension` package path, ignore documented wasm/webR harness outputs plus benchmark/conformance byproducts, and add `make clean_local` to purge the reproducible local spill files
- add `processing_threads` parameter to `duckhts_mosdepth(...)` for parallel contig processing: each worker opens independent BAM file handles via index queries, claims contigs atomically, and writes output in header order using an ordered-completion pattern; on the NA12878 16 GB WGS benchmark (4 decompression threads, 2 processing threads), fast mode is 1.38x faster than mosdepth v0.3.13, default mode 1.40x faster, and fragment mode 1.61x faster, all with byte-identical output; default is now `processing_threads := 2`
- change `duckhts_mosdepth(...)` defaults to `threads := 2` (decompression) and `processing_threads := 2` (parallel contigs) for better out-of-the-box WGS performance; previous defaults were `threads := 0, processing_threads := 0` (single-threaded)
- ship htslib public headers (`htslib/*.h`) and `libhts.a` in the installed R package under `inst/duckhts_extension/htslib/{include,lib}/`; add `inst/htslib_config.R.in` template (substituted at configure time into `inst/htslib_config.R`) providing `htslib_cflags()`, `htslib_libs()`, `htslib_rpath()`, and `htslib_version()` helpers for downstream R packages that need to link against the bundled htslib
- fix `configure.win` to stage htslib headers and `libhts.a` into `include/htslib/` and `lib/` (matching the Unix configure), so Windows installs also ship the headers for downstream linking
- add peak RSS memory reporting to `mosdepth_benchmark.py` using GNU `/usr/bin/time -v` for accurate per-process measurement; the benchmark summary table now includes peak RSS alongside timing
- optimize `duckhts_mosdepth(...)` BED output formatting: replace `ksprintf` with direct `kstring.h` integer/string formatters (`kputs`, `kputc`, `kputw`, `kputll`) in all hot BED-writing paths, eliminating ~16% CPU overhead from printf format parsing; on the NA12878 16 GB WGS benchmark (4 threads), default mode is now at parity with upstream mosdepth v0.3.13 (1.02x faster), fragment mode is 1.07x faster, fast mode gap narrows from 1.35x to 1.15x, and fast+500bp windows from 1.26x to 1.07x
- fix `duckhts_mosdepth(...)` conformance and benchmark scripts to perform full BED interval comparison (including zero-depth regions) instead of filtering to non-zero bins only; verified byte-identical output against upstream mosdepth v0.3.13 on the 16 GB NA12878 low-coverage WGS BAM (172 M reads, 3.1 Gbp) across fast, default, fragment, and windowed modes
- change native `bam_bin_counts(...)` to emit a dense fixed-bin layout across each selected contig span, including zero-count bins up to the contig end instead of only observed bins; this makes the function a principled binning primitive for downstream CNV/sample serializers, keeps `start`/`end` coordinates aligned to the full header-defined span, and updates the SQL/R docs and conformance tests accordingly
- add native `bam_bin_counts(...)` for fixed-width BAM/CRAM read-start binning: returns per-bin `count_total`, `count_fwd`, and `count_rev`, supports `mapq`, `require_flags`, `exclude_flags`, and `rmdup := 'none'|'flag'|'streaming'` (including WisecondorX-style larp/larp2 streaming duplicate suppression), and can optionally emit one-pass per-bin GC and MAPQ sufficient statistics via `stats := 'gc'`, `'mq'`, or `'gc,mq'`; reference GC is included when `reference := '...'` is provided. Also bundle the tiny WisecondorX BAM/CRAM fixtures in `test/data/` and `r/Rduckhts/inst/extdata/`, add SQL/R conformance coverage for duplicate modes, flag filtering, CRAM, and GC/MQ outputs, and add top-level README examples for native fixed-bin counting
- add `duckhts_mosdepth(...)` examples to the top-level README, including windowed fragment coverage output and read-back of the generated BED.gz files, and refresh the generated function-catalog text so the mosdepth description matches the current v0.3.13 parity surface
- extend native `duckhts_mosdepth(...)` coverage parity to the pinned local `mosdepth 0.3.13` surface: `fragment_mode := TRUE` now matches upstream `--fragment-mode` full-fragment insert coverage for proper pairs, default mode (`fast_mode := FALSE`) performs CIGAR-aware coverage with mate-overlap correction, `read_groups := '...'` filters RG tags, `min_frag_len` / `max_frag_len` filter absolute template length, and `use_median := TRUE` switches `by := '<window|bed>'` outputs from mean to median; regenerate the public function catalog and add SQL/R/native conformance coverage for BAM and CRAM fixtures across fast, fragment, default, overlap, read-group, quantize, threshold, and median cases
- extend `duckhts_mosdepth(...)` fast-mode parity with `quantize := '...'`, writing mosdepth-style `.quantized.bed.gz` + CSI output, and add native/SQL/R conformance coverage for quantized output plus explicit `by := '<bed>'` validation against upstream mosdepth
- extend `duckhts_mosdepth(...)` fast-mode parity with `thresholds := '...'` for `by := '<window|bed>'`, writing mosdepth-style `.thresholds.bed.gz` + CSI output; also align window/BED mean accumulation and window-region distribution bucketing with upstream mosdepth's current implementation behavior, and add SQL/R/native-conformance coverage for the new outputs
- import upstream mosdepth edge-case fixtures (`big`, `empty-tids`, `overlapping-pairs`, `ovl`, `nanopore`, and related BED files) into `test/data/` and `r/Rduckhts/inst/extdata/`, add a native conformance-suite runner for them, and record Brent Pedersen as the original mosdepth author in the R package metadata/copyright bundle
- switch `scripts/mosdepth_conformance.py` and `scripts/mosdepth_benchmark.py` to the native `duckhts_mosdepth(...)` table function only; the old pure-SQL mosdepth reconstruction path is removed from the validation/benchmark workflow, and larger WGS validation runs can now be driven directly from real BAM/CRAM inputs through the native writer
- expand `duckhts_mosdepth(...)`: the native mosdepth-compatible fast-mode rewrite now accepts indexed CRAM input via `fasta := ...` when required by htslib, and exposes `precision_digits := 2` as an explicit function argument instead of relying on the `MOSDEPTH_PRECISION` environment variable; add bundled SQL/R coverage for BAM and CRAM outputs plus explicit precision validation
- expand the source READMEs with runnable compression/indexing examples covering `bgzip(...)`, `bgunzip(...)`, `bam_index(...)`, `bcf_index(...)`, and `tabix_index(...)`, and regenerate the rendered README outputs
- make `read_bam(...)` htslib decompression workers configurable via `decompression_threads := 2`; the previous hardcoded thread count is now the documented default, and `decompression_threads := 0` disables per-file htslib worker threads
- speed up zero-column `COUNT(*)` scans across the HTS readers: `read_bam(...)`, `read_bcf(...)`, `read_tabix(...)`, `read_gff(...)`, `read_gtf(...)`, and indexed `read_bed(...)` now use index metadata for full-file count-only scans when DuckDB projects no output columns; `read_fasta(...)` uses `faidx` sequence counts when an index is available and otherwise counts FASTA headers directly; `read_fastq(...)` continues to count raw FASTQ records directly when no projected columns are needed, while preserving paired/interleaved validation errors
- add `hts_union_query(reader, pattern, params)` SQL scalar macro: generates a `UNION ALL BY NAME` query string that reads every file matching a glob pattern through any named reader function (`read_bam`, `read_bcf`, etc.), with a `filename` column identifying each row's source file; pass the result to `query()` to execute
- add R multi-file reading wrappers: `rduckhts_bam_multi`, `rduckhts_bcf_multi`, `rduckhts_fastq_multi`, `rduckhts_fasta_multi`, `rduckhts_bed_multi`, `rduckhts_tabix_multi`, `rduckhts_gff_multi`, `rduckhts_gtf_multi`; each follows the standard `(con, table_name, files, ..., overwrite)` convention, creates a DuckDB table with a `filename` column, and accepts an optional `.params` data.frame for per-file parameter overrides (e.g. per-sample regions or index paths)
- clarify the package README's browser/webR documentation: `r/Rduckhts/README.Rmd` now documents the full `Module.duckhtsWasmHttpConfig` parameter set (`headers`, `allowHosts`, `enforceHostAllowlist`, `withCredentials`, `allowInsecureAuth`), explicitly notes that webR consumers can set that config from R via `webr::eval_js()` without editing the host page, and covers practical wasm/browser behaviors such as same-origin setup, CORS requirements, `.csi` to `.tbi` fallback, and non-fatal `Range` warnings under the local `http.server` harness
- fix DuckDB community-extension wasm side-module ABI parity with the working `webR` build: the top-level `Makefile` now preserves `-sWASM_BIGINT` alongside `-fwasm-exceptions`, and both wasm build paths now share one extension-owned Emscripten compat header from `src/include/`; the CMake duckdb-wasm path additionally remaps htslib and extension-level native i64 libc calls (`lseek`, `ftruncate`, `strtoll`, `strtoull`, `atoll`, `time`) onto duckdb-wasm's `orig$...` exports while the webR/package path keeps the default symbols. Together this removes the wasm-only socket drift and DuckDB host ABI mismatches that were causing browser `LOAD` failures (`indirect call signature mismatch`, then `env.ftruncate` / `env.strtoll` signature mismatch) without affecting native builds
- align the local duckdb-wasm browser harness with the current extension build line: `scripts/start_duckdb_wasm_local_test.sh` now defaults to `@duckdb/duckdb-wasm` `1.31.0` (the first DuckDB `v1.4.x` runtime line instead of the older `1.29.0` / DuckDB `v1.1.1` default), and `scripts/duckdb-wasm-local-test.html` now logs `SELECT version()` before `LOAD` so runtime/header mismatches show up immediately in the page log
- separate the two browser wasm harnesses so they no longer stage artifacts into overlapping code-owned areas: `scripts/start_duckdb_wasm_local_test.sh` now defaults to `.duckdb-wasm-local-artifacts/` for its built extension and served site, while `scripts/start_webr_local_test.sh` now defaults to `.webr-local-artifacts/` for its built package tarball, local repo, and served site
- make the bundled `Rduckhts` wasm extension self-contained with respect to `htslib`: the Emscripten/webR path in `r/Rduckhts/configure` now builds only `libhts.a`, links `duckhts.duckdb_extension` directly against that static archive, and skips bundling or requiring shared `libhts.so*` files on wasm
- add a browser-native wasm `http` / `https` backend in `src/wasm_http_hfile.c`: the DuckDB side-module path now bypasses libcurl sockets entirely and uses synchronous XHR from the webR worker, which makes same-origin and CORS-enabled remote HTS URLs readable in browser wasm builds while keeping vendored `htslib` `libcurl`, `S3`, and `GCS` disabled on wasm
- keep `libcurl` disabled on wasm in `r/Rduckhts/configure`: `r-wasm/webr` ships `/opt/webr/wasm/lib/libcurl.a`, but libcurl socket calls from a SIDE_MODULE still hit a webR Emscripten message-bus failure (`resolved is not a function`) on first network use, so the package-owned XHR backend is now the supported wasm HTTP path
- harden wasm browser HTTP range behavior at EOF in `src/wasm_http_hfile.c`: cache remote object size from `Content-Range`/`Content-Length`, clamp range ends when size is known, short-circuit reads at/after EOF, and add `SEEK_END` size discovery fallback via `GET Range: bytes=0-0` when `HEAD` metadata is unavailable; this avoids cross-origin 416 `Requested Range Not Satisfiable` failures on `.tbi` index probes (e.g. GTEx tabix in webR/browser)
- harden non-Range wasm/browser HTTP fallback in `src/wasm_http_hfile.c`: when a server answers ranged requests with `200 OK`, cache the full object per open handle and serve subsequent reads from that in-memory cache (instead of re-downloading on every seek/read), while also emitting one-time warnings that Range is ignored and that large fallback payloads (>=64 MiB) are being used
- add optional wasm/browser request-header configuration in `src/wasm_http_hfile.c` via `Module.duckhtsWasmHttpConfig`: supports custom headers (including bearer auth), host allowlisting, optional `withCredentials`, and a default HTTPS-only guard for `Authorization` headers unless explicitly overridden (`allowInsecureAuth`)
- extend `Module.duckhtsWasmHttpConfig` with `enforceHostAllowlist`: when enabled, wasm HTTP requests are denied for hosts not matched by `allowHosts` (in addition to scoped header injection), giving a hard outbound-host policy for browser runtimes
- add a dedicated duckdb-wasm local test harness (`scripts/start_duckdb_wasm_local_test.sh`, `scripts/duckdb-wasm-local-test.html`) to validate the generic community-extension wasm path separately from webR, including extension load, local HTTP smoke reads, and runtime `Module.duckhtsWasmHttpConfig` header/auth config wiring
- add a project-owned local duckdb-wasm build image (`scripts/docker/duckdb-wasm-local.Dockerfile`) and wire `scripts/start_duckdb_wasm_local_test.sh` to use it when `emcmake` is unavailable on the host, so local wasm harness testing follows extension-ci-tools style `emsdk` + `vcpkg` provisioning instead of relying on a generic Emscripten image
- improve local duckdb-wasm harness isolation and browser compatibility: `scripts/start_duckdb_wasm_local_test.sh` now builds in an isolated mirror worktree in Docker mode (leaving host native `build/` and `cmake_build/` trees untouched), stages duckdb-wasm runtime assets locally, and `scripts/duckdb-wasm-local-test.html` loads those local runtime files to avoid cross-origin worker security errors
- switch `scripts/start_duckdb_wasm_local_test.sh` to a Docker-only workflow and move more setup into `scripts/docker/duckdb-wasm-local.Dockerfile` for cache-friendly local iteration (pre-provisioned `emsdk`, pinned `vcpkg`, and duckdb-wasm runtime assets baked into the image)
- refine duckdb-wasm harness runtime loading to avoid both bare-specifier ESM failures and cross-origin worker/module issues: stage `duckdb-browser.mjs`/worker/wasm locally, add browser import-map resolution for `apache-arrow`, and load all runtime artifacts from the local harness origin
- clarify current wasm/browser support scope: local-file workflows now work end-to-end in a real browser webR worker, including `rduckhts_load()`, bundled FASTQ/FASTA smoke tests, the installed `tinytest` suite, and same-origin browser HTTP tabix/GFF reads served by the local harness; cross-origin remote URLs still depend on normal browser CORS policy
- fix the bundled `Rduckhts` wasm side-module final link in `r/Rduckhts/configure`: the extension link now preserves `${LDFLAGS}` so the webR/Emscripten `SIDE_MODULE` settings actually reach `duckhts.duckdb_extension`, and it exports `duckhts_init_c_api` explicitly for DuckDB's loader; this fixes browser/webR `LOAD` failures where the package shipped a wasm file containing the symbol name but no usable init export
- adjust the bundled `Rduckhts` extension metadata trailer for wasm/browser builds to encode `linux_i686_musl` instead of `wasm_mvp` when `r/Rduckhts/configure` detects the Emscripten/webR path
- harden the bundled `Rduckhts` WebAssembly build path: the package `configure` script now preserves `rwasm`/r-universe `NAME=VALUE` cache overrides, forwards explicit `--build` / `--host` triplets into the vendored `htslib` `./configure`, forwards webR's Emscripten port flags for `zlib`/`bzip2`, seeds wasm-safe Autoconf cache results for `zlib`/`bzip2`/socket probes, injects a tiny Emscripten-only socket compatibility shim for `recv`/`send`/`closesocket`, and disables the optional `htslib` features that are not available in the stock webR/r-universe wasm toolchain (`libcurl`, `S3`, `GCS`, `lzma`, `plugins`); this fixes the original `ac_cv_func_getrandom=no: command not found` failure and subsequent nested `htslib` cross-compile probe failures without changing non-Wasm targets
- fix top-level DuckDB wasm packaging without changing `extension-ci-tools`: the CMake wasm build now rebuilds `libduckhts.a` as a fat archive containing vendored `htslib` (and any static archive dependencies CMake can see), so the downstream DuckDB wasm packager no longer wraps a bare extension archive with unresolved `htslib` symbols such as `bcf_readrec`

## duckhts 1.1.6 (2026-04-09)

- implement real `read_hts_index_spans(...)` chunk expansion for binning indexes: replace the previous SQL macro placeholder rows with a native table function that enumerates CSI/TBI/BAI chunk records from the bundled `htslib` index structures, populating `bin`, `chunk_beg_vo`, `chunk_end_vo`, `chunk_bytes`, `seq_start`, and `seq_end` with real values for low-level index inspection; avoid spurious `tbx` probe errors on BCF indexes by loading BCF indexes directly
- add `FILE_OFFSET UBIGINT` column to `read_bam(...)`: exposes the BGZF virtual file offset after reading each record via `bgzf_tell()` (a zero-cost macro over already-open BGZF struct fields). The offset is monotonically increasing within a sequential BAM scan and correctly encodes block address and block offset as `(block_address << 16) | block_offset`. This enables `ORDER BY FILE_OFFSET` in SQL window functions to reproduce exact BAM file order, which is required to faithfully replicate position-based streaming deduplication algorithms (such as WisecondorX's larp/larp2 state machine) without writing C or R loops.

## duckhts 1.1.5 (2026-04-08)

- fix `bcftools_liftover(...)` / `duckdb_liftover(...)` cache and realignment hardening: bound the per-thread chain/FASTA context cache instead of retaining every distinct context for the lifetime of DuckDB worker threads, and restore failure-safe scalar realignment handling so left-alignment does not reuse stale traceback state after allocation/zero-length alignment failures
- fix indexed parallel full scans with leading empty contigs in `read_bam(...)` and `read_bcf(...)`: replace recursive contig claiming with iterative retry loops so worker threads do not recurse/starve across long runs of empty contigs, and keep `read_bam(...)` filling the current DuckDB chunk after a successful contig handoff instead of returning an empty chunk that terminates the scan early

- switch the top-level extension `README.Rmd` examples from `R`/`DBI` wrapper execution to a custom knitr `sql` engine that runs `duckdb -unsigned` directly, matching the built extension workflow used in `DuckTinyCC`; point the liftover example at bundled fixtures instead of generating temporary FASTA/chain files in `R`
- fix Windows GNU CI path selection in vendored `htslib` CMake builds: pass the explicit DuckDB platform (`windows_amd64_mingw` vs `windows_amd64_rtools`) into CMake and treat the two MinGW/Rtools paths differently during `htslib ./configure`; keep `windows_amd64_mingw` close to `configure.win` with a small configure-time `LIBS` set, but restore the fuller static-libcurl dependency closure for `windows_amd64_rtools`, which still needs it during `htslib` feature probes; keep `CURL_STATICLIB` on the actual MinGW build objects rather than on `./configure` test programs

- fix Windows `windows_amd64_rtools` CMake toolchain selection: the top-level `Makefile` now asks `R CMD config` for `CC`/`AR`/`RANLIB` and passes those into the initial CMake configure, preventing mixed `C:/mingw64` compiler vs `Rtools` library setups during vendored `htslib` builds; simplify the MinGW `htslib` `ExternalProject` path back to separate configure/build steps instead of inlining `./configure && make` into one shell command; define `CURL_STATICLIB` for MinGW static-libcurl builds so vendored `htslib` and the extension link against Rtools static `libcurl.a` without `__imp_curl_*` import-symbol failures

- fix `read_bcf(...)` fixed-count INFO/FORMAT arrays: `Number=2`, `Number=4`, and other exact-cardinality fields now map to DuckDB list types instead of being truncated to the first element as scalar values
- fix `read_bcf(...)` string FORMAT lists such as DRAGEN `FORMAT/LAA`: non-`GT` string FORMAT fields with `Number != 1` now populate `VARCHAR[]` columns instead of incorrectly writing scalar `VARCHAR` values into list vectors, which caused DuckDB internal assertion failures on queries such as `SELECT *` against DRAGEN multi-sample gVCFs

- fix `bcftools_munge_row(...)` / `duckdb_munge(...)` multithreaded FASTA fetch races: stop sharing one cached `faidx_t` across DuckDB worker threads, keep FASTA index handles thread-local, and synchronize FASTA fetches in `munge`, eliminating intermittent `fai_retrieve` errors and aborts seen with `fasta_ref` under `PRAGMA threads > 1`

- expand `bcftools_score(...)` test coverage: add FORMAT/AS integer dosage fixture, missing DS value fixture, and seven GWAS summary preset fixtures (REGENIE, SAIGE, BOLT, METAL, PGS, SSF/GWAS-SSF); add SQL test cases for TSV counts without threshold, GWAS-VCF with `q_score_thr`, mutual-exclusion error checks (`include`+`exclude`, `regions`+`regions_file`, `targets`+`targets_file`); expand conformance scenarios from 3 to 8 cases (DS, HDS, AP dosage modes, PLINK2 preset, `use_variant_id`)

- fix `scripts/conformance/scenarios/score.sh` unbound variable crash: replace bare `BCFTOOLS_PLUGINS` self-reference with `SCORE_PLUGIN_DIR="${SCORE_PLUGIN_DIR:-}"` so the script runs correctly under `set -u`

- add `no_left_align` BOOLEAN parameter to `bcftools_liftover(...)` and `duckdb_liftover(...)`: when `true`, skips post-liftover left-alignment (Step 8), mirroring `bcftools +liftover --no-left-align`; default is `false`, preserving existing behavior

- make `bcftools_score(...)` filter-expression error handling safer in embedded DuckDB/R contexts: `bcftools_filter.c` now traps parse/eval failures internally and returns structured status/error strings instead of unwinding across extension frames with an external shim `setjmp`/`longjmp` path; this targets CI-only Linux/Windows score test crashes
- expand `bcftools_score(...)` API with explicit bcftools-style filtering controls: `regions`, `regions_file`, `regions_overlap`, `targets`, `targets_file`, `targets_overlap`, `apply_filters`, plus `include`/`exclude` expression filtering over core VCF fields (`POS`, `QUAL`, `CHROM`, `ID`, `REF`, `ALT`, `FILTER`)
- improve `bcftools_liftover(...)` multiallelic semantic parity with upstream: run the full indel/reference/left-alignment pipeline over dynamic allele arrays (not only biallelic scalar slots), preserve all ALT alleles when a new reference is introduced (`swap=-1`), normalize duplicate/ref-equal ALT alleles from output, and treat `ALT='.'` as no alternate alleles; add SQL coverage for multiallelic forward/reverse cases
- add GWAS-VCF multi-PRS scoring to `bcftools_score(...)`: when the summary file is a VCF with `FORMAT/ES` and `FORMAT/LP` per sample, each GWAS-VCF sample is treated as an independent PRS, producing per-PRS output metrics alongside the existing TSV summary scoring path
- add INFO/END liftover support to `bcftools_liftover(...)` and `duckdb_liftover(...)`: when an `end_pos` column is provided, the END position is lifted through the chain alongside the primary position with contig/strand/ordering validation; the output gains a `dest_end` column
- add mitochondria passthrough to `bcftools_liftover(...)` and `duckdb_liftover(...)`: when `lift_mt := false` (default), MT/chrM contigs with matching source/destination lengths are renamed without chain lifting, matching upstream `bcftools +liftover` behavior; `lift_mt := true` forces chain-based lifting
- harden `rduckhts_liftover()` argument validation for CI consistency: reject negative `max_snp_gap`/`max_indel_inc` in R before SQL execution so tinytests see stable wrapper-level errors across platforms/runners
- add fai-only mode to `bcftools_munge_row(...)` and `duckdb_munge(...)`: when `fasta_ref` is NULL, alleles pass through as-is without reference matching or allele swapping, matching upstream `bcftools +munge --fai`-only behavior

- port full indel liftover pipeline from upstream `bcftools +liftover` (Tier 3 — Algorithmic Completeness): allele extension framework (`scalar_realign_t`, `extend_alleles`, `pad_left`/`pad_right`/`pad_from_left`/`pad_from_right`), Needleman-Wunsch re-alignment with BWA-MEM affine gap scoring (`scalar_nw`, `scalar_clip_pad`, `scalar_get_shift`), sophisticated 3-phase `find_reference` with interior matching, SNP rescue hack, and anchor update, left-alignment post-liftover via VT normalization (`trim_right`/`trim_left`/`shift_left`/`is_left_aligned`), and SNP retry-as-indel dispatch restructure — adds ~630 lines of new C code with edge-of-contig bounds checking not present in upstream
- restructure liftover dispatch to match upstream 8-step pipeline: try `liftover_bp()` for ALL variants first, retry failed SNPs ("difficult SNPs") via the indel path with allele extension — SNPs at chain boundaries can now be rescued by extending into chain range (**breaking**: some previously-rejected edge-of-chain SNPs now map successfully with `note = 'Padded'`)

- harden `bcftools_score(...)` preset column parity: add `SCORE_HDR_UNUSED` sentinel and 45 missing column entries across all 8 presets (PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF) so that non-scoring columns (A2, FRQ, SE, N, INFO, etc.) are recognized and silently consumed during header matching, exactly matching upstream `score.h` preset definitions
- add gzipped summary file support to `bcftools_score(...)`: replace `fopen`/`fgets` with `hts_open`/`hts_getline`/`kstring_t` for transparent `.gz`/`.bgz` decompression of PGS Catalog and public GWAS summary files, also removing the fixed 32KB line buffer limitation
- fix `bcftools_score(...)` PRS name derivation to strip multiple known extensions iteratively (`{gz,txt,tsv,vcf,bcf}`), so `PGS000001.txt.gz` → `PGS000001` instead of `PGS000001.txt`, matching upstream `score.c:729-737`
- add validation for explicitly requested `use` tag in `bcftools_score(...)`: when user specifies `use := 'DS'` (or GT/HDS/AP/GP/AS), verify the tag exists in the VCF header and error if absent instead of silently scoring zero, matching upstream `score.c:767-795`
- fix `bcftools_score(...)` CNT column naming to match upstream: when `q_score_thr` and `counts` are both active, count columns are now named `<prs>_CNT_p<thr>` (upstream pattern) instead of `<prs>_p<thr>_CNT` (**breaking** column name change for q_score_thr + counts queries)
- fix `bcftools_score(...)` threshold boundary precision to match upstream: parse `q_score_thr` values with `strtof`→double promotion before `-log10`, reproducing the exact float→double comparison asymmetry in upstream `bcftools +score`; markers at exact P-value boundaries may now be excluded where they were previously included
- add `scripts/score_conformance.sh` end-to-end conformance test that runs both DuckHTS and upstream `bcftools +score`, UNPIVOTs to long format, and FULL OUTER JOINs with numeric tolerance to catch any output divergence
- fix `bcftools_score(...)` wrong-result bugs: remove incorrect METAL `Zscore → P` and SSF `standard_error → P` column mappings that produced wrong PRS scores, add haploid GT support for chrX male samples, and fix memory leak when summary markers are skipped due to missing effect sizes
- fix `duckdb_liftover(...)` wrong-result bugs: skip reverse complement for symbolic alleles (`<DEL>`, `*`, multi-alt) that cannot be complemented, and detect insertions (1-bp ref, multi-bp alt) as indels so they route through the indel liftover path instead of the SNP path
- fix `duckdb_munge(...)` wrong-result bugs: emit rows with `filter := 'MissingContig'` instead of aborting the query when FASTA ref fetch fails for unknown chromosomes/scaffolds, and propagate NAN correctly when swapping AC with missing NS (matches upstream arithmetic)
- harden `bcftools_score(...)` upstream parity: port `bcf_hdr_name2id_flexible()` for flexible chromosome name resolution (chr prefix strip/prepend, 23→X, 24→Y, 26/MT→chrM aliases), port numerically stable `-log10(p)` parsing for very small p-values via mantissa/exponent splitting, and handle `NA`/`.` as missing values in summary stats fields (BETA, OR, P, LP)
- add comprehensive `bcftools_score(...)` test coverage: DS/HDS/AP/GP dosage modes with real FORMAT values, OR-to-beta conversion, PLINK2 preset with LOG10_P, custom columns_file mapping, allele mismatch (zero-match), missing genotype (`./. `) handling, NA in summary stats, chr-prefix flexible matching, auto-detection priority, and small p-value threshold precision
- tighten `bcftools_score(...)` TSV matching parity: enforce marker-ID column requirements for `use_variant_id := true`, keep SNP-ID fallback behavior explicit when CHR/BP are unavailable, and add SQL coverage for rsID-vs-CHR/BP matching paths
- add METAL meta-analysis output support: `duckdb_munge_metal(...)` macro emits `si` (imputation info), `i2` (Cochran's I²), `cq` (Cochran's Q -log10 p), and `ed` (effect direction) columns from METAL-specific input fields (`INFO`/`HET_I2`/`HET_P`/`HET_LP`/`DIRE`); `ed` direction strings are correctly flipped on allele swap; `cq` is computed from `HET_LP` (precedence) or `-log10(HET_P)`; the base `duckdb_munge(...)` macro omits these four columns for non-METAL presets
- clarify orientation semantics by renaming `swapped` to `alleles_swapped` in `bcftools_munge_row(...)` / `duckdb_munge(...)`
- add `bcftools_munge_row(...)` and `duckdb_munge(...)` entries to `functions.yaml` so munge APIs are documented in the generated function catalog and community extension descriptor
- add `benchmark_munge.Rmd` and `make bench-munge` to benchmark DuckHTS munge against `bcftools +munge` with normalized output-group parity checks
- harden liftover error handling and diagnostics: `duckdb_liftover(...)` and `bcftools_liftover(...)` now raise explicit SQL errors for invalid required inputs, liftover load failures surface the underlying chain/FASTA cause, macro registration no longer fails silently, and new SQL error-stress tests cover 1,000,000-row invalid-input paths
- harden liftover upstream parity: port full IUPAC complement table (`rev_nt()`) so reverse-complement handles all ambiguity codes (R/Y/S/W/K/M/B/V/D/H/N) instead of collapsing unknowns to N; add `replace_iupac_codes()` to sanitize fetched reference sequences (non-ACGTN→N) matching upstream `fetch_sequence()`; add contig aliases `25→X` and `26→MT` to match upstream `bcf_hdr_name2id_flexible()`; add chain ID dedup via khash so Ensembl-duplicated chains are skipped (chains with id=0 are never deduplicated); add chain coverage validation rejecting malformed chains where block sizes don't fully cover target/query intervals; change regidx payload from `int` to `uint64_t` matching upstream
- align liftover result semantics more closely with `bcftools +liftover` by splitting the old `warning` field into `reject_reason` for rejected rows and `note` for emitted rows with extra annotations, using upstream-style reject names such as `MissingContig`, `UnmappedAnchors`, `MismatchAnchors`, and `MissingFasta`
- add README liftover examples and expand impossible-liftover coverage so unmapped rows and invalid-input error paths are exercised in both SQL and package-level tests
- add a `liftover(...)` table macro for score-style variant rows (`chrom`, `pos`, optional `ref`/`alt`) backed by a new `liftover_variant(...)` scalar kernel that uses UCSC chain files plus destination/source FASTA references to return lifted coordinates, lifted alleles, reverse-complement/swap indicators, and liftover status annotations
- source the generated community extension descriptor `version` field from the repo-level `description.yml` instead of duplicating it in `functions.yaml`
- add `quality_representation := 'phred'` to `read_bam(...)` and `read_fastq(...)` so base qualities can be returned as `UTINYINT[]` raw Phred values instead of SAM/FASTQ text
- add `input_quality_encoding := 'phred33' | 'auto' | 'phred64' | 'solexa64'` to `read_fastq(...)`; default to modern `phred33`, with optional legacy decoding and canonical Phred output on read
- add `detect_quality_encoding(...)` to inspect FASTQ quality ASCII ranges and report compatible encodings plus a heuristic guessed encoding
- add `sequence_encoding := 'nt16'` parameter to `read_bam(...)`, `read_fasta(...)`, and `read_fastq(...)` to return sequence data as `UTINYINT[]` using htslib nt16 4-bit codes (`=ACMGRSVTWYHKDBN` → 0-15) instead of decoded `VARCHAR` strings; defaults to `'string'` for backward compatibility
- refactor `seq_encode_4bit(...)` / `seq_decode_4bit(...)` to use htslib's `seq_nt16_table[]` and `seq_nt16_str[]` directly instead of custom switch tables; `U` (RNA) now encodes as `T` (code 8) and code 0 (`=`) is accepted; **breaking**: unknown characters (e.g. `!`) now map to `N` (code 15) instead of returning NULL, matching htslib behavior — all encoding paths (UDF + reader `sequence_encoding := 'nt16'`) are now unified on the same shared code
- stop advertising unsupported `attributes_map := TRUE` on generic `read_tabix(...)`; parsed attribute maps remain available on `read_gff(...)` and `read_gtf(...)`
- centralize CSQ/ANN/BCSQ subfield typing for `read_bcf(...)` with builtin rules derived from `bcftools +split-vep`, conservative ANN defaults, and an `additional_csq_column_types := ...` override parameter using bcftools-style `PATTERN TYPE` entries
- add BGZF compression and decompression table functions: `bgzip(...)` and `bgunzip(...)`, both defaulting to preserving the source file unless `keep := FALSE` is requested
- add HTS index builders: `bam_index(...)`, `bcf_index(...)`, and `tabix_index(...)`
- add HTS metadata readers: `read_hts_header(...)`, `read_hts_index(...)`, `read_hts_index_spans(...)`, and `read_hts_index_raw(...)`
- add interval readers/helpers: `read_bed(...)` for BED3-BED12 input and `fasta_nuc(...)` for bedtools nuc-style FASTA interval composition over BED intervals or fixed-width bins
- add sequence helpers: `seq_encode_4bit(...)`, `seq_decode_4bit(...)`, `seq_gc_content(...)`, and `seq_kmers(...)`
- add SAM flag decoders: `sam_flag_bits(...)` for bulk struct output and `sam_flag_has(...)` for direct mask checks
- rename ambiguous SAM flag helpers to explicit/spec-aligned names: `is_paired(...)`, `is_proper_pair(...)`, `is_next_segment_unmapped(...)`, and `is_next_segment_reverse_complemented(...)`
- add `is_forward_aligned(...)` as a strand helper with SAM-spec semantics for mapped segments
- add `CIGAR utils` scalar helpers for soft clips, operator detection, and query/reference length calculations
- extend `read_bam(...)` with `standard_tags := TRUE` typed SAM tag columns and `auxiliary_tags := TRUE` for the remaining tags as a map
- improve `read_tabix(...)`, `read_gff(...)`, and `read_gtf(...)` with header-based column names, basic type inference, explicit column type overrides, and tabix metadata-aware parsing
- harden region-query behavior for files with incomplete contig metadata by returning empty results with a warning instead of failing
- add SQL coverage for BED reading, FASTA composition metrics, the new metadata readers, sequence/SAM-flag helpers, typed tabix parsing, and BGZF/index round-trips

## duckhts 0.1.1.9000 (2026-02-10)

- validate paired FASTQ mate files by exact QNAME match and equal record counts
- detect odd-length interleaved FASTQ input and raise an error
- return empty results (with warning) for region queries when contig headers are missing
- add non-conforming VCF fixture without ##contig for region query tests
- generic tabix reader now respects tabix header/meta configuration when inferring columns
- read_tabix supports header-based column names via header := true and header_names
- read_bam supports standard_tags (typed SAMtags columns) and auxiliary_tags (map of remaining tags)
- standard_tags + auxiliary_tags demo added to README R examples
- read_tabix supports auto_detect for basic numeric inference and column_types overrides
- added SQL/R demos and tests for tabix type inference and column overrides
