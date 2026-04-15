# Better Scan Planning Notes

This file captures the code changes needed for better indexed scan planning
without implementing them yet.

Current decision:

- Do implement later: reduce tail latency for indexed full-file scans by using
  weighted contig claiming instead of plain header-order claiming.
- Do not implement yet: expose a BCF/VCF record-offset column. That is a
  significant downstream-facing surface change and needs explicit sign-off on
  semantics first.

## 1. Tail latency: weighted contig claiming

### Goal

Keep the existing "one contig owned once" parallel scan semantics for indexed
full-file scans, but claim heavier contigs first so worker threads do not get
stuck on one late large contig while others go idle.

This is a scheduling change, not a query-semantics change.

### Why this is the right next step

- It preserves the current no-dup behavior of full indexed scans.
- It avoids reopening the known duplicate/overlap issues of multi-region BCF
  planning.
- It improves the long tail of parallel scans with minimal API impact.
- It keeps deterministic claim order if the weight sort uses a stable tie-break
  on original contig id.

### Scope

Apply this to the indexed full-file parallel paths only:

- `read_bam(...)` in `src/bam_reader.c`
- `read_bcf(...)` in `src/bcf_reader.c`

Do not change:

- explicit region scans
- multi-region chained scans
- count-only zero-column fast paths

### Required BAM changes

File: `src/bam_reader.c`

1. Extend `bam_bind_data_t` with planning metadata:
   - `uint64_t *contig_weights;`
   - `int *claim_order;`
   - `int n_claimable_contigs;`

2. Free the new arrays in `destroy_bam_bind()`.

3. After `sam_index_load3(...)` succeeds in `bam_read_bind()`:
   - allocate one weight per `tid`
   - call `hts_idx_get_stat(idx, tid, &mapped, &unmapped)` for each contig
   - set weight to `mapped + unmapped`
   - build `claim_order` as the list of contig ids sorted by:
     1. descending weight
     2. ascending original `tid` as the tie-break

4. Decide whether to keep zero-weight contigs:
   - preferred: omit zero-weight contigs from `claim_order` when stats are
     available and exact
   - fallback-safe option: keep them at the end with weight `0`

5. In `bam_read_global_init()`:
   - use `bind->n_claimable_contigs` instead of raw `bind->n_contigs` for
     `global->n_contigs` and `duckdb_init_set_max_threads(...)`

6. Extend `bam_global_init_data_t` with:
   - `int *claim_order;`

7. In `bam_read_global_init()`, point `global->claim_order` at
   `bind->claim_order`.

8. In `claim_next_contig()`:
   - replace `next` as direct `tid`
   - use `tid = global->claim_order[next]`
   - keep `assigned_contig = tid`
   - continue using `sam_itr_queryi(local->idx, tid, 0, HTS_POS_MAX)`

9. Keep the existing zero-indexed-read retry loop, but it should almost never
   trigger once zero-weight contigs are filtered or placed last.

### Required BCF changes

File: `src/bcf_reader.c`

1. Extend `bcf_bind_data_t` with:
   - `uint64_t *contig_weights;`
   - `int *claim_order;`
   - `int n_claimable_contigs;`

2. Free the new arrays in `destroy_bind_data()`.

3. In the index-probing section of `bcf_read_bind()`:
   - after index load succeeds and contig names are known
   - compute one weight per `tid` with `hts_idx_get_stat(idx_or_tbx->idx, tid, &mapped, &unmapped)`
   - use `mapped + unmapped` as the scheduling weight
   - build a stable `claim_order` sorted by descending weight, then ascending
     original `tid`

4. In `bcf_global_init_data_t`, add:
   - `int *claim_order;`

5. In `bcf_read_global_init()`:
   - use `bind->n_claimable_contigs` instead of `bind->n_contigs`
   - point `global->claim_order` to `bind->claim_order`
   - cap threads against `n_claimable_contigs`

6. In `claim_next_contig()`:
   - fetch the next slot from the atomic counter
   - resolve `tid = global->claim_order[next]`
   - use `tid` to select `contig_name` and build the iterator

7. Do not change the multi-region BCF path. The duplicate-risk boundary stays
   exactly where it is today: chained region iterators.

### Shared helper recommendation

To avoid duplicating sort logic, add a small helper in each file or one shared
internal helper module that:

- accepts `n_contigs`
- accepts a `uint64_t weights[]`
- fills `claim_order[]`
- sorts by `(weight desc, tid asc)`

If this is shared across readers, wire the new source file through:

- `CMakeLists.txt`
- `r/Rduckhts/bootstrap.R`
- `r/Rduckhts/configure`
- `r/Rduckhts/configure.win`

### Testing needed

SQL:

- keep `test/sql/parallel_empty_contigs.test`
- add one regression that exercises a weighted-claim path on a file with
  uneven contig densities and verifies full-row count correctness
- avoid tests that assert physical output order; this is a scheduling change

R:

- extend existing BAM/BCF end-to-end reader tests only if a deterministic
  correctness assertion is useful
- no wrapper signature changes are needed for weighted claiming alone

Benchmarking:

- add or update a benchmark script that logs per-scan wall clock under
  `PRAGMA threads > 1`
- success criterion is reduced tail time, not different row contents

### Non-goals

- no new user parameter
- no change to full-scan result set
- no region splitting below contig level
- no use of overlapping region strings for planning

## 2. Deferred: BCF/VCF record offset

### Why this is deferred

Adding a BCF/VCF record-offset column is a significant downstream-facing change.
Before implementing it, we need agreement on:

- column name
- exact semantics
- stability guarantees across sequential vs indexed scans
- stability guarantees across BCF vs text VCF paths
- whether the value is a BGZF virtual offset, iterator next-offset, or some
  weaker best-effort location marker

This should not be merged as an incidental part of scan-planning work.

### Intended outcome

Provide a stable ordering key analogous to BAM `FILE_OFFSET`, so downstream SQL
can request deterministic file-order semantics explicitly when needed.

### Why it is harder than BAM

BCF/VCF has multiple scan paths:

- sequential BCF reads
- indexed BCF iterator reads
- indexed tabix VCF text reads followed by `vcf_parse1(...)`

Those paths do not obviously expose the same "record offset after this record"
semantics that BAM currently gets from `bgzf_tell(...)`.

### Questions to settle before coding

1. Should the column be named `FILE_OFFSET`, `VIRTUAL_OFFSET`, or something
   BCF-specific?
2. Must the value mean the same thing for BCF and text VCF?
3. Is it acceptable for the value to be present only for BGZF-backed paths and
   NULL otherwise?
4. For indexed iterator scans, should the value reflect:
   - the start of the current record
   - the position after the current record
   - the iterator's next seek point
5. Can remote scans expose the same semantics reliably?

### Likely files to change later

- `src/bcf_reader.c`
- `test/sql/duckhts.test`
- `r/Rduckhts/inst/tinytest/test_seq_ops.R`
- `functions.yaml`
- `NEWS.md`
- `r/Rduckhts/NEWS.md`
- possibly `r/Rduckhts/R/duckhts.R` only if wrapper docs need explicit mention

### Recommended rollout

1. write down the chosen semantics first
2. implement only after explicit approval
3. add schema tests plus monotonicity/order tests
4. document clearly that downstream users should `ORDER BY` the offset column
   when exact file order matters

## 3. Things not to do

- Do not replace contig claiming with arbitrary overlapping region planning for
  full-file scans.
- Do not mix tail-latency scheduling work with a public BCF offset column in
  one commit.
- Do not make parallel output order promises unless a stable ordering column is
  exposed and the query explicitly orders by it.
