Revision: f1498fbe6a2b8d67fdbcdd9f08a48f90ca50d4f5. Workload:
deterministic biallelic synthetic reference frequencies and PC loadings;
three independently perturbed samples; 4 DuckDB threads; two reference
groups; two PCs; `min_cor = 0.4`. Each variant is present in both groups
and PCs, and every input row matches. This report measures query and
materialization together. The public-reference workload is measured
separately below.

The correction vector is staged from the committed, checksummed bigsnpr
vignette product. `VmHWM` is cumulative process high-water RSS (MiB),
including package loading and all preceding cases; the baseline before
the first run was 178 MiB. It is not per-query allocated memory.

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 15%" />
<col style="width: 8%" />
<col style="width: 13%" />
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 6%" />
<col style="width: 7%" />
<col style="width: 12%" />
<col style="width: 15%" />
<col style="width: 6%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">samples</th>
<th style="text-align: right;">variants_per_sample</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">output_group_rows</th>
<th style="text-align: right;">groups</th>
<th style="text-align: right;">pcs</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">repeat_id</th>
<th style="text-align: right;">elapsed_seconds</th>
<th style="text-align: right;">process_peak_rss_mb</th>
<th style="text-align: left;">failures</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">3000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.089</td>
<td style="text-align: right;">219.887</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: right;">3</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">3000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">0.086</td>
<td style="text-align: right;">221.449</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">3000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">0.083</td>
<td style="text-align: right;">223.793</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: right;">3</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">15000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.101</td>
<td style="text-align: right;">226.137</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">15000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">0.099</td>
<td style="text-align: right;">227.074</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: right;">3</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">15000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">0.098</td>
<td style="text-align: right;">229.262</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">51000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">0.143</td>
<td style="text-align: right;">234.574</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: right;">3</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">51000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">0.143</td>
<td style="text-align: right;">238.168</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">51000</td>
<td style="text-align: right;">6</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">3</td>
<td style="text-align: right;">0.146</td>
<td style="text-align: right;">242.699</td>
<td style="text-align: left;"></td>
</tr>
</tbody>
</table>

Failed gates are retained in the `failures` column for each repeat. The
nearest identical synthetic workload is the earlier rendered
`benchmarks/benchmark_ancestry.md` at repository revision `8ab38b99`
(source revision `56aad9cf`): its three 17,000-site repeats had median
0.078 seconds, versus 0.143 seconds here. These are separate R processes
with uncontrolled background load; this difference cannot establish a
regression or equivalence. There is no ancestry-projection baseline on
`develop`; the synthetic rows are not evidence about real-reference
loading or BAM/CRAM depth sensitivity.

## Public reference and GRCh37 chr22 genotypes

Two checksum-verified Figshare CSV products each contain 5,816,590 rows:
21 reference-frequency groups and 16 PC loadings. A staged, key-sorted
Parquet join holds the four site-key columns, 16 DOUBLE loadings and 21
DOUBLE frequencies. Its SHA-256 is
`a0bf01ae16db605ff0add59e5ebb747c35454406900c11fd8596f23c167fe02b`; the
staging receipt binds this digest to both registered CSV digests and
DuckDB 1.5.5. The one-time conversion is excluded from query times.
Converting each compressed CSV to narrow typed Parquet before the keyed
sort keeps the CSV scanner and sort from holding wide decompressed rows
simultaneously; reading just the compressed header also avoids a full
R-side CSV scan.

<table>
<colgroup>
<col style="width: 13%" />
<col style="width: 26%" />
<col style="width: 5%" />
<col style="width: 9%" />
<col style="width: 7%" />
<col style="width: 8%" />
<col style="width: 5%" />
<col style="width: 8%" />
<col style="width: 7%" />
<col style="width: 8%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">implementation</th>
<th style="text-align: left;">source_revision</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">reference_rows</th>
<th style="text-align: right;">pc_columns</th>
<th style="text-align: right;">group_columns</th>
<th style="text-align: right;">seconds</th>
<th style="text-align: right;">peak_rss_mib</th>
<th style="text-align: right;">budget_mib</th>
<th style="text-align: left;">within_budget</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">wide_csv_join</td>
<td
style="text-align: left;">29f3e695489e172dbec04a7d5323159cf78a38c7</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">322.26</td>
<td style="text-align: right;">15269.00</td>
<td style="text-align: right;">2560</td>
<td style="text-align: left;">FALSE</td>
</tr>
<tr class="even">
<td style="text-align: left;">astra_narrow_parquet</td>
<td style="text-align: left;">independent_prototype</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">24.92</td>
<td style="text-align: right;">2491.38</td>
<td style="text-align: right;">2560</td>
<td style="text-align: left;">TRUE</td>
</tr>
<tr class="odd">
<td style="text-align: left;">keyed_narrow_parquet</td>
<td
style="text-align: left;">5494ebdf5d80d324b6dc8efee749b87361b7be55</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">28.93</td>
<td style="text-align: right;">2480.53</td>
<td style="text-align: right;">2560</td>
<td style="text-align: left;">TRUE</td>
</tr>
</tbody>
</table>

All three staging rows use the same two public CSV products and four
threads. The independently run Astra prototype is an external
comparator, not a DuckHTS revision. Source conversion and digest
validation are included in staging RSS; these figures are not query
memory. The epilepsy summary input has 4,880,492 records. A separate
public phase-3 chr22 VCF supplies 11 individuals from AFR (3), EUR (2),
EAS (2), SAS (2) and AMR (2). Each workload uses four DuckDB threads and
the bigsnpr vignette’s 16 correction coefficients. These are
complete-process peak RSS measurements, including package loading,
CSV-based bigsnpr oracles and preceding operations in that process. The
SQL time includes scanning the staged Parquet view, matching sites and
aggregation, but excludes loading the oracle CSVs.

<table style="width:100%;">
<colgroup>
<col style="width: 17%" />
<col style="width: 4%" />
<col style="width: 11%" />
<col style="width: 6%" />
<col style="width: 10%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 8%" />
<col style="width: 13%" />
<col style="width: 6%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">workload</th>
<th style="text-align: right;">samples</th>
<th style="text-align: right;">variants_per_sample</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">output_group_rows</th>
<th style="text-align: right;">groups</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">reference_rows</th>
<th style="text-align: right;">reference_load_seconds</th>
<th style="text-align: right;">sql_seconds</th>
<th style="text-align: right;">process_peak_rss_mib</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">epilepsy (chr22)</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">25.51</td>
<td style="text-align: right;">2.24</td>
<td style="text-align: right;">7595.44</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy (genome-wide)</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">25.55</td>
<td style="text-align: right;">4.91</td>
<td style="text-align: right;">7740.84</td>
</tr>
<tr class="odd">
<td style="text-align: left;">1000G phase-3 chr22 genotypes</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">26.84</td>
<td style="text-align: right;">2.54</td>
<td style="text-align: right;">4871.47</td>
</tr>
</tbody>
</table>

Epilepsy matching audit: bigsnpr matched 3343235 / 4880492 records
genome-wide (6 ambiguous SNPs excluded, 1,576,097 reversed, 0
strand-flipped). DuckHTS matched 3343235 / 3343235 pre-matched inputs,
after matching against the real reference product. The independently
selected, position-unique, unambiguous chr22 subset is 17000 / 38061
eligible chr22 matches, with 17000 / 17000 matched by DuckHTS. Both
genome-wide solvers consume the same allele-aligned sites.

<table>
<thead>
<tr class="header">
<th style="text-align: left;">group_id</th>
<th style="text-align: right;">proportion</th>
<th style="text-align: right;">oracle</th>
<th style="text-align: right;">absolute_difference</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Africa (East)</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Africa (North)</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Africa (South)</td>
<td style="text-align: right;">0.0028022</td>
<td style="text-align: right;">0.0028022</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Africa (West)</td>
<td style="text-align: right;">0.0065084</td>
<td style="text-align: right;">0.0065084</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Ashkenazi</td>
<td style="text-align: right;">0.0175426</td>
<td style="text-align: right;">0.0175426</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Asia (East)</td>
<td style="text-align: right;">0.0308440</td>
<td style="text-align: right;">0.0308440</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Bangladesh</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Europe (North East)</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Europe (South East)</td>
<td style="text-align: right;">0.1394387</td>
<td style="text-align: right;">0.1394387</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Europe (South West)</td>
<td style="text-align: right;">0.0213679</td>
<td style="text-align: right;">0.0213679</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Finland</td>
<td style="text-align: right;">0.0650006</td>
<td style="text-align: right;">0.0650006</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Ireland</td>
<td style="text-align: right;">0.0873397</td>
<td style="text-align: right;">0.0873397</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Italy</td>
<td style="text-align: right;">0.0339252</td>
<td style="text-align: right;">0.0339252</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Japan</td>
<td style="text-align: right;">0.0025044</td>
<td style="text-align: right;">0.0025044</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle East</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Pakistan</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Philippines</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Scandinavia</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">South America</td>
<td style="text-align: right;">0.0004969</td>
<td style="text-align: right;">0.0004969</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">Sri Lanka</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0.0000000</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">United Kingdom</td>
<td style="text-align: right;">0.5922292</td>
<td style="text-align: right;">0.5922292</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

Genome-wide `cor_pred`: DuckHTS 0.999629260942, bigsnpr 0.999629260942.
The conservative seven-decimal QP rounding correlation bound is
0.000008450; the 0.4 gate is farther away. On chr22, both return
0.9995518 at seven decimals with group difference 0; its rounding bound
is 0.000008978. All 21 genome-wide coefficients agree at seven decimals;
the script checks rounded equality.

The individual input is restricted from 1055452 unique biallelic SNVs to
70371 same-assembly, allele-agreeing, unambiguous reference overlaps
before spaced site selection. Both solvers see the same per-sample
genotypes. bigsnpr’s imperfect-match warning (`cor_pred < 0.99`) is
retained; passing the 0.4 gate does not imply that the frequencies match
perfectly. The 11 individuals’ 231 coefficients agree with bigsnpr at
seven decimals.

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 6%" />
<col style="width: 8%" />
<col style="width: 8%" />
<col style="width: 13%" />
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 4%" />
<col style="width: 8%" />
<col style="width: 9%" />
<col style="width: 8%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">sample_id</th>
<th style="text-align: left;">superpopulation</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">duckhts_sites</th>
<th style="text-align: right;">bigsnpr_sites</th>
<th style="text-align: right;">max_group_difference</th>
<th style="text-align: right;">cor_pred</th>
<th style="text-align: right;">oracle_cor_pred</th>
<th style="text-align: left;">status</th>
<th style="text-align: left;">oracle_error</th>
<th style="text-align: right;">rounding_bound</th>
<th style="text-align: right;">gate_distance</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: left;">EUR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6617542</td>
<td style="text-align: right;">0.6617542</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.5e-06</td>
<td style="text-align: right;">0.2617542</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: left;">EAS</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.7073448</td>
<td style="text-align: right;">0.7073448</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">7.6e-06</td>
<td style="text-align: right;">0.3073448</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00731</td>
<td style="text-align: left;">AMR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6567316</td>
<td style="text-align: right;">0.6567316</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.6e-06</td>
<td style="text-align: right;">0.2567316</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG01112</td>
<td style="text-align: left;">AMR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6480976</td>
<td style="text-align: right;">0.6480976</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.7e-06</td>
<td style="text-align: right;">0.2480976</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG01879</td>
<td style="text-align: left;">AFR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6832799</td>
<td style="text-align: right;">0.6832799</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.0e-06</td>
<td style="text-align: right;">0.2832799</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG02057</td>
<td style="text-align: left;">EAS</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6770518</td>
<td style="text-align: right;">0.6770518</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">7.6e-06</td>
<td style="text-align: right;">0.2770518</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG02561</td>
<td style="text-align: left;">AFR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6978570</td>
<td style="text-align: right;">0.6978570</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">7.8e-06</td>
<td style="text-align: right;">0.2978570</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG03742</td>
<td style="text-align: left;">SAS</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6385715</td>
<td style="text-align: right;">0.6385715</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.5e-06</td>
<td style="text-align: right;">0.2385715</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG03866</td>
<td style="text-align: left;">SAS</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6878663</td>
<td style="text-align: right;">0.6878663</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.5e-06</td>
<td style="text-align: right;">0.2878663</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA12878</td>
<td style="text-align: left;">EUR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.6653880</td>
<td style="text-align: right;">0.6653880</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">8.6e-06</td>
<td style="text-align: right;">0.2653880</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: left;">AFR</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">20000</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0.7353511</td>
<td style="text-align: right;">0.7353511</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">7.5e-06</td>
<td style="text-align: right;">0.3353511</td>
</tr>
</tbody>
</table>

## Keyed Parquet profile and input scaling

The pre-conversion genome-wide DuckDB JSON profile, on revision
`0a5fd898c2bb8449259e01dd7aff161bb4dc3bc8`, measured 30.27 s and 40,232
MiB peak buffer memory. Its complete R process peaked at 38,783 MiB RSS
while retaining the CSV-based bigsnpr oracle. The 53,491,760-row
materialized `basis` and 1,123,326,960-row `x` join dominate
cardinality; `pred_sites` separately groups 3,343,235 variants after a
long-reference join. The old `list()` calls assemble only the small
post-aggregation solver vectors; there is no per-variant list
construction. The profile also shows hash joins and grouping over the
long reference, but DuckDB does not attribute its single reported peak
buffer value to individual operators. CSV parsing occurs in R before
this SQL profile, not in a DuckDB scan operator. Sorting is limited to
the input window partition and final group ordering.

The keyed query scans projected Parquet columns and aggregates 336
`sum(PC * group_frequency)` values, 16 corrected
`sum(PC * aligned_frequency)` values and group correlations per sample.
The 336-element solver vector and 21 group correlations are assembled
**after** aggregation. The earlier keyed plan grouped all 5,816,590
reference loci into a site hash and retained 3,343,235 classified rows
while the correlation branch rebuilt aligned input. Its four-thread
fresh-process run took 3.43 s, 1,872 MiB RSS and 3,175 MiB DuckDB peak
buffer. The old JSON profile attributes 0.80 s to the reference grouping
and has two large join branches; profile peak memory is query-wide, not
operator-attributed.

Astra’s external plan converts CSVs to narrow Parquet, then builds its
solver hash from the 3,343,235 matched inputs, not the 5,816,590
reference sites. Its solver-only process ran 2.39–2.43 s with 587 MiB
RSS and zero spill under `memory_limit='1GB'`;
classification/preparation ran separately. DuckHTS keeps classification
and its complete audit in the measured wrapper call: it writes a
temporary 3,343,235-row, approximately 22.6 MiB aligned Parquet
relation, then materializes only per-sample moments, solver output and
correlation. The single-sample reductions use ungrouped aggregates. For
multiple samples, the 1..N sample index from the small audit table is
the group key: a string `GROUP BY sample_id` led DuckDB to estimate
440,000 groups for 22 samples and spill 951 MiB with one thread, while
the indexed 44-sample query has zero spill. The temporary Parquet file
is removed before returning; its bytes are explicit bounded intermediate
I/O, not a DuckDB spill.

The hash-build input is 95.65 MiB of decoded key/frequency payload for
the full epilepsy run (3,343,235 rows, measured by
`benchmark_ancestry_query_stages.R`). With a fixed 512 MiB R/DuckDB
overhead and 3× this payload, the process-RSS budget is 799 MiB, further
capped at 1 GiB under `memory_limit='1GB'`. Its measured 4-thread RSS
below is within both limits. Every stage’s JSON profile reported zero
temp spill; the largest query-session DuckDB buffer estimate was 1,218
MiB. DuckDB buffer accounting is not resident process RSS and must not
be added to it. Genotype 44-sample build payload is 27.69 MiB (880,000
rows, including the four-byte sample index used for aggregation), for a
595 MiB budget; the one-thread process measured 476 MiB RSS with zero
spill.

`benchmark_ancestry_reference_load.R` scans every numeric column in each
source CSV and the staged Parquet product with DuckDB; neither scan
materializes 5.8 million rows in R. These timings differ from the CSV
`fread()` oracle-loading times above and have their own denominator.

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 16%" />
<col style="width: 17%" />
<col style="width: 26%" />
<col style="width: 30%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">threads</th>
<th style="text-align: right;">reference_rows</th>
<th style="text-align: right;">numeric_columns</th>
<th style="text-align: right;">csv_full_column_seconds</th>
<th style="text-align: right;">parquet_full_column_seconds</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">37</td>
<td style="text-align: right;">19.62</td>
<td style="text-align: right;">1.57</td>
</tr>
<tr class="even">
<td style="text-align: right;">4</td>
<td style="text-align: right;">5816590</td>
<td style="text-align: right;">37</td>
<td style="text-align: right;">14.28</td>
<td style="text-align: right;">0.44</td>
</tr>
</tbody>
</table>

`benchmark_ancestry_scaling.R` runs three independent R processes per
point with one or four DuckDB threads and the same registry-staged
`ancestry_reference_parquet`. Files are reused from the local cache
across processes; the OS page cache is not cleared between repetitions.
Input Parquet is emitted by `test/scripts/ancestry_bigsnpr_real.R` and
`test/scripts/ancestry_1000g_real.R` from registry-staged
`ancestry_epilepsy`, `ancestry_1000g_chr22` and the 1000G sample panel;
the scaling driver accepts those output paths as arguments. Epilepsy
1×/2×/4× contains nested prefixes of 835,809/1,671,618/3,343,235
physical input variants (one sample). The sample series keeps 20,000
chr22 sites per person while selecting 11/22/44 distinct phase-3 VCF
individuals. The joint series selects 5,000/10,000/20,000 sites for
11/22/44 individuals: both axes grow, with 55,000/220,000/880,000 input
rows. No sample IDs or genotypes are cloned. Every input row
participates; each group produces one output row. Query time covers the
public wrapper, including reference validation, aligned file I/O,
solver, quality gates and result materialisation, but excludes one-time
input staging, source checksums and bigsnpr oracle runs. Peak RSS
includes R and DuckDB; DuckDB buffer peaks are not process RSS. The
aligned Parquet size is reported separately from DuckDB temporary spill.
The single-thread 1× floor determines whether a workload has a timing
verdict; sub-five-second series are memory profiles only. Each
repetition has a 1 GB DuckDB memory limit and zero DuckDB spill
allowance. The declared process-RSS ceilings are 512 MiB fixed overhead
plus three times the full-scale decoded required-column payload: 798.95
MiB for epilepsy (95.65 MiB decoded), and 595.07 MiB for both sample and
joint cases (27.69 MiB decoded), capped at 1 GiB. These ceilings apply
to every point, including the joint series. The keyed reference scans
5,816,590 loci with 16 DOUBLE PC and 21 DOUBLE group columns (296
decoded numeric bytes per row plus four keys). The admitted aligned
relation holds at most the physical input row count: 3,343,235 dense
summary-statistic sites or 880,000 joint genotype rows, not sites × PCs
× groups. The largest per-sample aggregate consumes 3,343,235 sites and
emits 16 × (21 + 1) numeric moments; each fitted sample returns 21 rows.
The aligned file is written once and scanned for moments and prediction;
its reported size is logical file I/O, not DuckDB spill or measured
physical disk traffic. An optional single-run TEMP-table probe retained
identical proportions and zero spill, but its full-epilepsy peak was
approximately 865 MiB, above the 798.95 MiB ceiling. It is not a
repeated comparable workload; the measured path here keeps the aligned
Parquet file. The 21 full-epilepsy coefficients and 231 individual
coefficients agree with bigsnpr at seven decimals.

<table>
<caption>Keyed-query baseline at 29f3e695; genotype 22/44 repeated
sample IDs</caption>
<colgroup>
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 11%" />
<col style="width: 10%" />
<col style="width: 6%" />
<col style="width: 10%" />
<col style="width: 13%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">workload</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">scale</th>
<th style="text-align: right;">samples</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">used_variants</th>
<th style="text-align: right;">output_rows</th>
<th style="text-align: right;">seconds</th>
<th style="text-align: right;">peak_rss_mib</th>
<th style="text-align: right;">peak_buffer_mib</th>
<th style="text-align: right;">peak_temp_mib</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">3.89</td>
<td style="text-align: right;">877.57</td>
<td style="text-align: right;">1210.64</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">6.35</td>
<td style="text-align: right;">1005.39</td>
<td style="text-align: right;">1559.88</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">11.60</td>
<td style="text-align: right;">1769.77</td>
<td style="text-align: right;">2996.53</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">1.35</td>
<td style="text-align: right;">1066.92</td>
<td style="text-align: right;">1404.67</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">1.88</td>
<td style="text-align: right;">1106.41</td>
<td style="text-align: right;">1740.59</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">3.43</td>
<td style="text-align: right;">1872.41</td>
<td style="text-align: right;">3174.94</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">1.37</td>
<td style="text-align: right;">695.71</td>
<td style="text-align: right;">849.84</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">1.72</td>
<td style="text-align: right;">695.96</td>
<td style="text-align: right;">890.09</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">2.40</td>
<td style="text-align: right;">750.57</td>
<td style="text-align: right;">1081.69</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">0.53</td>
<td style="text-align: right;">655.62</td>
<td style="text-align: right;">962.30</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">0.70</td>
<td style="text-align: right;">698.44</td>
<td style="text-align: right;">1029.94</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">1.05</td>
<td style="text-align: right;">827.65</td>
<td style="text-align: right;">1222.04</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

Keyed-query baseline at 29f3e695; genotype 22/44 repeated sample IDs

<table>
<caption>Single-run bounded-engine baseline at 58fb3264; reference
admission contract differs</caption>
<colgroup>
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 11%" />
<col style="width: 10%" />
<col style="width: 6%" />
<col style="width: 10%" />
<col style="width: 13%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">workload</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">scale</th>
<th style="text-align: right;">samples</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">used_variants</th>
<th style="text-align: right;">output_rows</th>
<th style="text-align: right;">seconds</th>
<th style="text-align: right;">peak_rss_mib</th>
<th style="text-align: right;">peak_buffer_mib</th>
<th style="text-align: right;">peak_temp_mib</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">2.97</td>
<td style="text-align: right;">301.81</td>
<td style="text-align: right;">247.05</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">5.52</td>
<td style="text-align: right;">394.24</td>
<td style="text-align: right;">419.52</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">10.92</td>
<td style="text-align: right;">586.05</td>
<td style="text-align: right;">799.24</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">1.01</td>
<td style="text-align: right;">407.61</td>
<td style="text-align: right;">473.30</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">1.75</td>
<td style="text-align: right;">504.02</td>
<td style="text-align: right;">639.62</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">3.26</td>
<td style="text-align: right;">690.68</td>
<td style="text-align: right;">979.23</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">0.69</td>
<td style="text-align: right;">258.81</td>
<td style="text-align: right;">174.30</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">1.17</td>
<td style="text-align: right;">332.45</td>
<td style="text-align: right;">299.12</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">2.10</td>
<td style="text-align: right;">471.25</td>
<td style="text-align: right;">546.02</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">0.39</td>
<td style="text-align: right;">286.07</td>
<td style="text-align: right;">255.58</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">0.58</td>
<td style="text-align: right;">356.35</td>
<td style="text-align: right;">373.84</td>
<td style="text-align: right;">0</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">0.96</td>
<td style="text-align: right;">494.19</td>
<td style="text-align: right;">614.06</td>
<td style="text-align: right;">0</td>
</tr>
</tbody>
</table>

Single-run bounded-engine baseline at 58fb3264; reference admission
contract differs

<table>
<caption>Three fresh processes per point; spread is IQR and range;
exponent is log median-time ratio / log input-row ratio to preceding
scale</caption>
<colgroup>
<col style="width: 7%" />
<col style="width: 5%" />
<col style="width: 4%" />
<col style="width: 3%" />
<col style="width: 4%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 4%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 6%" />
<col style="width: 7%" />
<col style="width: 8%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 4%" />
<col style="width: 7%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: left;">workload</th>
<th style="text-align: right;">threads</th>
<th style="text-align: right;">scale</th>
<th style="text-align: right;">samples</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">output_rows</th>
<th style="text-align: right;">decoded_mib</th>
<th style="text-align: right;">median_s</th>
<th style="text-align: right;">iqr_s</th>
<th style="text-align: right;">min_s</th>
<th style="text-align: right;">max_s</th>
<th style="text-align: right;">max_rss_mib</th>
<th style="text-align: right;">max_buffer_mib</th>
<th style="text-align: right;">aligned_file_mib</th>
<th style="text-align: right;">spill_mib</th>
<th style="text-align: right;">budget_mib</th>
<th style="text-align: right;">exponent</th>
<th style="text-align: left;">timing_verdict</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">epilepsy.1.1</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">23.91</td>
<td style="text-align: right;">10.01</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">10.01</td>
<td style="text-align: right;">10.05</td>
<td style="text-align: right;">558.45</td>
<td style="text-align: right;">746.04</td>
<td style="text-align: right;">5.60</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy.1.2</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">47.83</td>
<td style="text-align: right;">12.00</td>
<td style="text-align: right;">0.03</td>
<td style="text-align: right;">11.98</td>
<td style="text-align: right;">12.04</td>
<td style="text-align: right;">558.40</td>
<td style="text-align: right;">772.39</td>
<td style="text-align: right;">11.18</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">0.26</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy.1.4</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">95.65</td>
<td style="text-align: right;">16.12</td>
<td style="text-align: right;">0.03</td>
<td style="text-align: right;">16.09</td>
<td style="text-align: right;">16.16</td>
<td style="text-align: right;">654.05</td>
<td style="text-align: right;">931.67</td>
<td style="text-align: right;">22.41</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">0.43</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy.4.1</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">835809</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">23.91</td>
<td style="text-align: right;">2.94</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">2.91</td>
<td style="text-align: right;">2.96</td>
<td style="text-align: right;">475.47</td>
<td style="text-align: right;">675.84</td>
<td style="text-align: right;">5.65</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epilepsy.4.2</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1671618</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">47.83</td>
<td style="text-align: right;">3.48</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">3.47</td>
<td style="text-align: right;">3.51</td>
<td style="text-align: right;">547.62</td>
<td style="text-align: right;">740.28</td>
<td style="text-align: right;">11.30</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">0.24</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">epilepsy.4.4</td>
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">95.65</td>
<td style="text-align: right;">4.71</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">4.70</td>
<td style="text-align: right;">4.72</td>
<td style="text-align: right;">742.08</td>
<td style="text-align: right;">1124.17</td>
<td style="text-align: right;">22.59</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">798.95</td>
<td style="text-align: right;">0.43</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes.1.1</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">6.92</td>
<td style="text-align: right;">8.41</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">8.39</td>
<td style="text-align: right;">8.43</td>
<td style="text-align: right;">558.56</td>
<td style="text-align: right;">773.39</td>
<td style="text-align: right;">0.39</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes.1.2</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">13.85</td>
<td style="text-align: right;">8.96</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">8.89</td>
<td style="text-align: right;">8.99</td>
<td style="text-align: right;">558.69</td>
<td style="text-align: right;">773.39</td>
<td style="text-align: right;">0.76</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.09</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes.1.4</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">27.69</td>
<td style="text-align: right;">9.95</td>
<td style="text-align: right;">0.04</td>
<td style="text-align: right;">9.93</td>
<td style="text-align: right;">10.01</td>
<td style="text-align: right;">558.33</td>
<td style="text-align: right;">779.01</td>
<td style="text-align: right;">1.32</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.15</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes.4.1</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">6.92</td>
<td style="text-align: right;">2.51</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">2.50</td>
<td style="text-align: right;">2.52</td>
<td style="text-align: right;">474.15</td>
<td style="text-align: right;">649.24</td>
<td style="text-align: right;">0.40</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">genotypes.4.2</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">440000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">13.85</td>
<td style="text-align: right;">2.71</td>
<td style="text-align: right;">0.00</td>
<td style="text-align: right;">2.70</td>
<td style="text-align: right;">2.71</td>
<td style="text-align: right;">469.12</td>
<td style="text-align: right;">676.83</td>
<td style="text-align: right;">0.76</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.11</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes.4.4</td>
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">27.69</td>
<td style="text-align: right;">3.15</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">3.14</td>
<td style="text-align: right;">3.17</td>
<td style="text-align: right;">567.76</td>
<td style="text-align: right;">733.14</td>
<td style="text-align: right;">1.30</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.22</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">joint.1.1</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">55000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">1.73</td>
<td style="text-align: right;">8.00</td>
<td style="text-align: right;">0.05</td>
<td style="text-align: right;">7.99</td>
<td style="text-align: right;">8.09</td>
<td style="text-align: right;">559.06</td>
<td style="text-align: right;">781.27</td>
<td style="text-align: right;">0.10</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">joint.1.2</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">6.92</td>
<td style="text-align: right;">8.43</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">8.40</td>
<td style="text-align: right;">8.44</td>
<td style="text-align: right;">559.26</td>
<td style="text-align: right;">781.27</td>
<td style="text-align: right;">0.39</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.04</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">joint.1.4</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">27.69</td>
<td style="text-align: right;">10.02</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">9.99</td>
<td style="text-align: right;">10.03</td>
<td style="text-align: right;">559.07</td>
<td style="text-align: right;">781.27</td>
<td style="text-align: right;">1.32</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.13</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">joint.4.1</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">55000</td>
<td style="text-align: right;">231</td>
<td style="text-align: right;">1.73</td>
<td style="text-align: right;">2.39</td>
<td style="text-align: right;">0.02</td>
<td style="text-align: right;">2.36</td>
<td style="text-align: right;">2.40</td>
<td style="text-align: right;">468.70</td>
<td style="text-align: right;">676.74</td>
<td style="text-align: right;">0.10</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">NA</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="odd">
<td style="text-align: left;">joint.4.2</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">22</td>
<td style="text-align: right;">220000</td>
<td style="text-align: right;">462</td>
<td style="text-align: right;">6.92</td>
<td style="text-align: right;">2.63</td>
<td style="text-align: right;">0.01</td>
<td style="text-align: right;">2.61</td>
<td style="text-align: right;">2.63</td>
<td style="text-align: right;">466.92</td>
<td style="text-align: right;">676.48</td>
<td style="text-align: right;">0.38</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.07</td>
<td style="text-align: left;">timed</td>
</tr>
<tr class="even">
<td style="text-align: left;">joint.4.4</td>
<td style="text-align: left;">joint</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">44</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">924</td>
<td style="text-align: right;">27.69</td>
<td style="text-align: right;">3.16</td>
<td style="text-align: right;">0.03</td>
<td style="text-align: right;">3.16</td>
<td style="text-align: right;">3.22</td>
<td style="text-align: right;">564.49</td>
<td style="text-align: right;">731.61</td>
<td style="text-align: right;">1.30</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">595.07</td>
<td style="text-align: right;">0.13</td>
<td style="text-align: left;">timed</td>
</tr>
</tbody>
</table>

Three fresh processes per point; spread is IQR and range; exponent is
log median-time ratio / log input-row ratio to preceding scale

<table>
<caption>Nearest identical input baseline (single run, 58fb3264) versus
repeated checked-reference wrapper; absolute differences are not a
no-regression claim</caption>
<colgroup>
<col style="width: 10%" />
<col style="width: 12%" />
<col style="width: 9%" />
<col style="width: 13%" />
<col style="width: 8%" />
<col style="width: 14%" />
<col style="width: 15%" />
<col style="width: 15%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">workload</th>
<th style="text-align: right;">input_rows</th>
<th style="text-align: right;">median_s</th>
<th style="text-align: right;">max_rss_mib</th>
<th style="text-align: right;">seconds</th>
<th style="text-align: right;">peak_rss_mib</th>
<th style="text-align: right;">extra_seconds</th>
<th style="text-align: right;">extra_rss_mib</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">epilepsy</td>
<td style="text-align: right;">3343235</td>
<td style="text-align: right;">4.71</td>
<td style="text-align: right;">742.08</td>
<td style="text-align: right;">3.26</td>
<td style="text-align: right;">690.68</td>
<td style="text-align: right;">1.46</td>
<td style="text-align: right;">51.39</td>
</tr>
<tr class="even">
<td style="text-align: left;">genotypes</td>
<td style="text-align: right;">880000</td>
<td style="text-align: right;">3.15</td>
<td style="text-align: right;">567.76</td>
<td style="text-align: right;">0.96</td>
<td style="text-align: right;">494.19</td>
<td style="text-align: right;">2.18</td>
<td style="text-align: right;">73.57</td>
</tr>
</tbody>
</table>

Nearest identical input baseline (single run, 58fb3264) versus repeated
checked-reference wrapper; absolute differences are not a no-regression
claim

Full-reference duplicate and finite-column admission scans account for
additional fixed work across all input scales, especially the smaller
sample runs. The comparison uses the same staged input denominators but
different validation contracts; the latency/RSS trade-off requires
explicit review before merging.

## Indexed 30x CRAM sensitivity

These retained measurements were produced at
`0a5fd898c2bb8449259e01dd7aff161bb4dc3bc8` with the existing long-form
BAM ancestry path, not remeasured for the keyed Parquet projection.
Three public 30x CRAMs (NA18507/AFR, HG00188/EUR, HG00403/EAS) were read
via their indexes at chr22:15,335,303–50,800,284 (GRCh38). The
registered GRCh37-to-GRCh38 chain maps 19,981/20,000 selected GRCh37
loci onto chr22; 430 mapped loci require reference swaps. The 19,551
unswapped, unique SNVs provide a spaced 17,000-site panel. No whole CRAM
was downloaded. Each 25% read subsample uses samtools seed 42; median
depth and allele-read totals are measured after downsampling. Each mode
sees the same `min_depth = 7` criterion and the same VCF site subset.
Each query emits 21 group rows with four DuckDB threads and two htslib
workers. Peak RSS is the R process high-water mark (MiB) after reference
loading and regional CRAM staging; it excludes samtools child-process
RSS and is not incremental SQL memory. The `seconds` column is the
ancestry BAM query alone, excluding CRAM staging and depth audit. The
native and quarter-depth columns retain every status, including failed
gates and their `cor_pred` values.

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">sites</th>
<th style="text-align: left;">panel_sha256</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1000</td>
<td
style="text-align: left;">0c0f0ec14cc1f94e0394787839d5361142545da58557b4aa2dd18279c801328f</td>
</tr>
<tr class="even">
<td style="text-align: right;">5000</td>
<td
style="text-align: left;">49f46df0de88a5df86edf45d6be8f661c77b404e8963eb9b2d959d301ffc3cb0</td>
</tr>
<tr class="odd">
<td style="text-align: right;">17000</td>
<td
style="text-align: left;">9490246d8372fe5bbc0e3b20ee12dc648f3a57b37e6b3ba6ef7ee48e040798b6</td>
</tr>
</tbody>
</table>

<table style="width:100%;">
<colgroup>
<col style="width: 4%" />
<col style="width: 2%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 7%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 4%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 7%" />
<col style="width: 6%" />
<col style="width: 3%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">sample_id</th>
<th style="text-align: right;">sites</th>
<th style="text-align: left;">depth</th>
<th style="text-align: left;">mode</th>
<th style="text-align: right;">vcf_used</th>
<th style="text-align: right;">bam_used</th>
<th style="text-align: right;">depth_eligible</th>
<th style="text-align: right;">median_depth</th>
<th style="text-align: right;">allele_reads</th>
<th style="text-align: right;">vcf_cor_pred</th>
<th style="text-align: right;">cor_pred</th>
<th style="text-align: left;">status</th>
<th style="text-align: right;">max_group_difference</th>
<th style="text-align: right;">rounding_bound</th>
<th style="text-align: right;">gate_distance</th>
<th style="text-align: right;">seconds</th>
<th style="text-align: right;">process_peak_rss_mib</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">994</td>
<td style="text-align: right;">994</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">35276</td>
<td style="text-align: right;">0.7368</td>
<td style="text-align: right;">0.7236</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3236</td>
<td style="text-align: right;">1.925</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">966</td>
<td style="text-align: right;">994</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">35276</td>
<td style="text-align: right;">0.7368</td>
<td style="text-align: right;">0.7315</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3315</td>
<td style="text-align: right;">1.950</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">784</td>
<td style="text-align: right;">784</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">8933</td>
<td style="text-align: right;">0.7368</td>
<td style="text-align: right;">0.7272</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3272</td>
<td style="text-align: right;">0.699</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">726</td>
<td style="text-align: right;">784</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">8933</td>
<td style="text-align: right;">0.7368</td>
<td style="text-align: right;">0.7530</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3530</td>
<td style="text-align: right;">0.697</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">992</td>
<td style="text-align: right;">992</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">32929</td>
<td style="text-align: right;">0.6637</td>
<td style="text-align: right;">0.6520</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.1791</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2520</td>
<td style="text-align: right;">2.064</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">967</td>
<td style="text-align: right;">992</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">32929</td>
<td style="text-align: right;">0.6637</td>
<td style="text-align: right;">0.6621</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.1098</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2621</td>
<td style="text-align: right;">1.811</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">686</td>
<td style="text-align: right;">686</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8230</td>
<td style="text-align: right;">0.6637</td>
<td style="text-align: right;">0.6361</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.3586</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2361</td>
<td style="text-align: right;">0.620</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">620</td>
<td style="text-align: right;">686</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8230</td>
<td style="text-align: right;">0.6637</td>
<td style="text-align: right;">0.6676</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.4107</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2676</td>
<td style="text-align: right;">0.617</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">991</td>
<td style="text-align: right;">991</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">33339</td>
<td style="text-align: right;">0.7023</td>
<td style="text-align: right;">0.6933</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0335</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2933</td>
<td style="text-align: right;">1.906</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">960</td>
<td style="text-align: right;">991</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">33339</td>
<td style="text-align: right;">0.7023</td>
<td style="text-align: right;">0.6966</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0439</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2966</td>
<td style="text-align: right;">1.886</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">716</td>
<td style="text-align: right;">716</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8410</td>
<td style="text-align: right;">0.7023</td>
<td style="text-align: right;">0.6888</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0155</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2888</td>
<td style="text-align: right;">0.639</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">1000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">1000</td>
<td style="text-align: right;">645</td>
<td style="text-align: right;">716</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">8410</td>
<td style="text-align: right;">0.7023</td>
<td style="text-align: right;">0.7200</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0278</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3200</td>
<td style="text-align: right;">0.633</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4961</td>
<td style="text-align: right;">4961</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">174362</td>
<td style="text-align: right;">0.7279</td>
<td style="text-align: right;">0.7229</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3229</td>
<td style="text-align: right;">2.749</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4811</td>
<td style="text-align: right;">4961</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">174362</td>
<td style="text-align: right;">0.7279</td>
<td style="text-align: right;">0.7300</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3300</td>
<td style="text-align: right;">2.556</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3762</td>
<td style="text-align: right;">3762</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">43452</td>
<td style="text-align: right;">0.7279</td>
<td style="text-align: right;">0.7175</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3175</td>
<td style="text-align: right;">0.799</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3445</td>
<td style="text-align: right;">3762</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">43452</td>
<td style="text-align: right;">0.7279</td>
<td style="text-align: right;">0.7552</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3552</td>
<td style="text-align: right;">0.797</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4966</td>
<td style="text-align: right;">4966</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">164947</td>
<td style="text-align: right;">0.6557</td>
<td style="text-align: right;">0.6490</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0912</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2490</td>
<td style="text-align: right;">2.317</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4798</td>
<td style="text-align: right;">4966</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">164947</td>
<td style="text-align: right;">0.6557</td>
<td style="text-align: right;">0.6553</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0684</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2553</td>
<td style="text-align: right;">2.336</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3470</td>
<td style="text-align: right;">3470</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">41323</td>
<td style="text-align: right;">0.6557</td>
<td style="text-align: right;">0.6327</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0892</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2327</td>
<td style="text-align: right;">0.748</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3189</td>
<td style="text-align: right;">3470</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">41323</td>
<td style="text-align: right;">0.6557</td>
<td style="text-align: right;">0.6613</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0892</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2613</td>
<td style="text-align: right;">0.753</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4963</td>
<td style="text-align: right;">4963</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">165959</td>
<td style="text-align: right;">0.7108</td>
<td style="text-align: right;">0.7075</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0250</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3075</td>
<td style="text-align: right;">2.475</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">4822</td>
<td style="text-align: right;">4963</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">165959</td>
<td style="text-align: right;">0.7108</td>
<td style="text-align: right;">0.7154</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0348</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3154</td>
<td style="text-align: right;">2.478</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3562</td>
<td style="text-align: right;">3562</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">41694</td>
<td style="text-align: right;">0.7108</td>
<td style="text-align: right;">0.6983</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0109</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2983</td>
<td style="text-align: right;">0.775</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">5000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">5000</td>
<td style="text-align: right;">3274</td>
<td style="text-align: right;">3562</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">41694</td>
<td style="text-align: right;">0.7108</td>
<td style="text-align: right;">0.7316</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0205</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3316</td>
<td style="text-align: right;">0.779</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16888</td>
<td style="text-align: right;">16888</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">593948</td>
<td style="text-align: right;">0.7338</td>
<td style="text-align: right;">0.7260</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3260</td>
<td style="text-align: right;">3.183</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16335</td>
<td style="text-align: right;">16888</td>
<td style="text-align: right;">35</td>
<td style="text-align: right;">593948</td>
<td style="text-align: right;">0.7338</td>
<td style="text-align: right;">0.7340</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3340</td>
<td style="text-align: right;">2.924</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">12832</td>
<td style="text-align: right;">12832</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">148596</td>
<td style="text-align: right;">0.7338</td>
<td style="text-align: right;">0.7175</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3175</td>
<td style="text-align: right;">0.947</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">NA18507</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">11795</td>
<td style="text-align: right;">12832</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">148596</td>
<td style="text-align: right;">0.7338</td>
<td style="text-align: right;">0.7520</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0000</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3520</td>
<td style="text-align: right;">1.016</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16879</td>
<td style="text-align: right;">16879</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">559198</td>
<td style="text-align: right;">0.6575</td>
<td style="text-align: right;">0.6490</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0021</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2490</td>
<td style="text-align: right;">2.587</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16373</td>
<td style="text-align: right;">16879</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">559198</td>
<td style="text-align: right;">0.6575</td>
<td style="text-align: right;">0.6567</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0021</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2567</td>
<td style="text-align: right;">2.591</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">11794</td>
<td style="text-align: right;">11794</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">139660</td>
<td style="text-align: right;">0.6575</td>
<td style="text-align: right;">0.6327</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0129</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2327</td>
<td style="text-align: right;">0.884</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00188</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">10813</td>
<td style="text-align: right;">11794</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">139660</td>
<td style="text-align: right;">0.6575</td>
<td style="text-align: right;">0.6661</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0007</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2661</td>
<td style="text-align: right;">0.874</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16881</td>
<td style="text-align: right;">16881</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">564975</td>
<td style="text-align: right;">0.7041</td>
<td style="text-align: right;">0.6975</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0243</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2975</td>
<td style="text-align: right;">2.709</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">native</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">16343</td>
<td style="text-align: right;">16881</td>
<td style="text-align: right;">33</td>
<td style="text-align: right;">564975</td>
<td style="text-align: right;">0.7041</td>
<td style="text-align: right;">0.7062</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0173</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3062</td>
<td style="text-align: right;">2.721</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="odd">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">allele_fraction</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">12132</td>
<td style="text-align: right;">12132</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">141852</td>
<td style="text-align: right;">0.7041</td>
<td style="text-align: right;">0.6844</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0699</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.2844</td>
<td style="text-align: right;">0.905</td>
<td style="text-align: right;">4871.648</td>
</tr>
<tr class="even">
<td style="text-align: left;">HG00403</td>
<td style="text-align: right;">17000</td>
<td style="text-align: left;">quarter</td>
<td style="text-align: left;">called_genotype</td>
<td style="text-align: right;">17000</td>
<td style="text-align: right;">11156</td>
<td style="text-align: right;">12132</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">141852</td>
<td style="text-align: right;">0.7041</td>
<td style="text-align: right;">0.7189</td>
<td style="text-align: left;">ok</td>
<td style="text-align: right;">0.0716</td>
<td style="text-align: right;">0.0011</td>
<td style="text-align: right;">0.3189</td>
<td style="text-align: right;">0.895</td>
<td style="text-align: right;">4871.648</td>
</tr>
</tbody>
</table>

0 / 36 correlation gates failed. The minimum distance from the 0.4 gate
among the measurable cases is 0.2326775; the largest
coefficient-rounding correlation bound across the three site counts is
0.0010957. Called-genotype frequency excludes sites whose balance rule
makes no call; its `bam_used` denominator is therefore smaller than the
depth-eligible count. Maximum group differences compare the same-site
VCF result with each BAM/CRAM mode; they are sensitivity, not an
equality assertion.
