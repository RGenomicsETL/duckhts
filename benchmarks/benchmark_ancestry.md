Revision: c235ef57fefa94b8ea90a2d154b76193079cc861. Workload:
deterministic biallelic synthetic reference frequencies and PC loadings;
three independently perturbed samples; 4 DuckDB threads; two reference
groups; two PCs; `min_cor = 0.4`. Each variant is present in both groups
and PCs, and every input row matches. This report measures query and
materialization together. The public-reference workload is measured
separately below.

The correction vector is staged from the committed, checksummed bigsnpr
vignette product. `VmHWM` is cumulative process high-water RSS (MiB),
including package loading and all preceding cases; the baseline before
the first run was 183.9 MiB. It is not per-query allocated memory.

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
<td style="text-align: right;">0.039</td>
<td style="text-align: right;">222.660</td>
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
<td style="text-align: right;">0.040</td>
<td style="text-align: right;">232.504</td>
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
<td style="text-align: right;">0.041</td>
<td style="text-align: right;">241.879</td>
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
<td style="text-align: right;">0.049</td>
<td style="text-align: right;">245.160</td>
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
<td style="text-align: right;">0.049</td>
<td style="text-align: right;">247.660</td>
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
<td style="text-align: right;">0.049</td>
<td style="text-align: right;">250.941</td>
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
<td style="text-align: right;">0.085</td>
<td style="text-align: right;">267.660</td>
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
<td style="text-align: right;">0.081</td>
<td style="text-align: right;">271.566</td>
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
<td style="text-align: right;">0.082</td>
<td style="text-align: right;">275.473</td>
<td style="text-align: left;"></td>
</tr>
</tbody>
</table>

Failed gates are retained in the `failures` column for each repeat. The
nearest identical synthetic workload is the earlier rendered
`benchmarks/benchmark_ancestry.md` at repository revision `8ab38b99`
(source revision `56aad9cf`): its three 17,000-site repeats had median
0.078 seconds, versus 0.082 seconds here. These are separate R processes
with uncontrolled background load; this difference cannot establish a
regression or equivalence. There is no ancestry-projection baseline on
`develop`; the synthetic rows are not evidence about real-reference
loading or BAM/CRAM depth sensitivity.

## Public reference and GRCh37 chr22 genotypes

Two checksum-verified Figshare products each contain 5,816,590 rows: 21
reference-frequency groups and 16 PC loadings. The epilepsy summary
input has 4,880,492 records. A separate public phase-3 chr22 VCF
supplies 11 individuals from AFR (3), EUR (2), EAS (2), SAS (2) and AMR
(2). Each workload uses four DuckDB threads and the bigsnpr vignette’s
16 correction coefficients. These are complete-process peak RSS
measurements, including package loading and preceding operations in that
process; the SQL time excludes reading CSVs, matching sites, and
materializing reference relations.

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
<td style="text-align: right;">25.45</td>
<td style="text-align: right;">0.22</td>
<td style="text-align: right;">7596.47</td>
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
<td style="text-align: right;">25.48</td>
<td style="text-align: right;">31.09</td>
<td style="text-align: right;">37337.41</td>
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
<td style="text-align: right;">26.64</td>
<td style="text-align: right;">1.10</td>
<td style="text-align: right;">4871.58</td>
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

Genome-wide `cor_pred`: DuckHTS 0.999629260844, bigsnpr 0.999629260942.
The conservative seven-decimal QP rounding correlation bound is
0.000008450; the 0.4 gate is farther away. On chr22, both return
0.9995518 at seven decimals with group difference 0; its rounding bound
is 0.000008978. Tolerance: 2e-4 per group.

The individual input is restricted from 1055452 unique biallelic SNVs to
70371 same-assembly, allele-agreeing, unambiguous reference overlaps
before spaced site selection. Both solvers see the same per-sample
genotypes. bigsnpr’s imperfect-match warning (`cor_pred < 0.99`) is
retained; passing the 0.4 gate does not imply that the frequencies match
perfectly.

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
<td style="text-align: right;">0.6770519</td>
<td style="text-align: right;">0.6770518</td>
<td style="text-align: left;">ok</td>
<td style="text-align: left;"></td>
<td style="text-align: right;">7.6e-06</td>
<td style="text-align: right;">0.2770519</td>
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

## Indexed 30x CRAM sensitivity

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
