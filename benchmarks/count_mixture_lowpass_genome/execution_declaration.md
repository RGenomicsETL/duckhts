Status: pre-run execution contract for the genome-wide low-pass campaign.

# Memory-bounded execution declaration

The campaign retains the declaration in `declaration.R`: 90 strata, 18,000
replicates, 75 primary strata, both sequencing-error values, and the unchanged
model, calibration grid, independent oracle, statistical criteria, memory
ceiling formula and no-spill requirement. Every replicate uses its declared
seed and all 4,883,749 panel sites to simulate depth.

Fourteen long-lived, one-thread workers each load all chromosome panels once.
Input identities are verified once at startup. The measured microbatch sizes
are 4, 2 and 1 replicates at depths 0.1, 0.2 and 0.4 respectively. Each call uses
one bulk-loaded counts table with `sample_id = '<stratum>-<replicate>'`.
Fresh-process timing and RSS remain confined to the separate scaling series.

Only covered sites enter the campaign counts table. The simulator retains the
same random draws and counts as the all-site implementation. Full-panel sites,
zero-depth sites, excluded sites and unavailable sites are recorded separately
for every replicate; native covered-table QC is not presented as full-panel QC.

`results/covered_equivalence.csv` retains three replicates per depth. All nine
all-site/covered-site comparisons are bit-identical for `share`, `share_lower`,
`share_upper`, `share_se`, `singles_share`, `pair_residual`, `status`,
`relation_diagnostic`, `blocks` and `pair_blocks`. Simulated covered rows and QC
also match exactly. All-site controls use seven threads and registered frames;
covered controls use one thread and bulk-loaded tables, additionally checking
those execution differences. Both retain the configured 6 GB limit and zero
spill. The full-site controls' reported buffer peaks reach 7.60 GB, within the
original overhead-plus-three-times-live-state ceiling; the configured limit is
not a measured peak-memory guarantee for those controls.

`results/microbatch_probe.csv` measures the chosen covered batches in one
long-lived connection: 1,857,208 / 1,772,959 / 1,611,064 counts rows, 1.66 GB
maximum DuckDB buffer and zero spill. Simulation, bulk loading and fitting
cost 1.50 / 2.60 / 4.16 seconds per replicate in this probe. The depth-weighted
18,000-replicate projection is about 59 minutes on fourteen workers; allowing
startup, independent checks and concurrency overhead, expected campaign wall
time is about 1.3 hours, subject to shared memory-bandwidth contention. A full
200-replicate stratum is timed before concurrent execution.

Before collection, one replicate per stratum is checked against the independent
R moment oracle, refitted alone with bit-identical rows, and compared at one
and four DuckDB threads. No gate is calibrated from campaign outcomes.

Probe implementation SHA-256:

- `batch.R`: `7880abe9ab728637d9146822f1157142f8af0a16faaad82bffd41125e1c5b386`
- `probe.R`: `1605807acd713aeca0e693b87770c12fb508f1e78260c0fd88c4be497efe2cbe`
