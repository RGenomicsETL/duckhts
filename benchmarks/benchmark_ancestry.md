Revision: 03536163eec7de307cc74375bf1009fa7d885ecb. Workload:
deterministic biallelic synthetic reference frequencies and PC loadings;
three independently perturbed samples; 4 DuckDB threads; two reference
groups; two PCs; `min_cor = 0.4`. Each variant is present in both groups
and PCs, and every input row matches. This report measures query and
materialization together, not the 850 MB Figshare products.

The correction vector is staged from the committed, checksummed bigsnpr
vignette product. `VmHWM` is cumulative process high-water RSS (MiB),
including package loading and all preceding cases; the baseline before
the first run was 179.7 MiB. It is not per-query allocated memory.

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
<td style="text-align: right;">222.242</td>
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
<td style="text-align: right;">0.038</td>
<td style="text-align: right;">226.461</td>
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
<td style="text-align: right;">0.038</td>
<td style="text-align: right;">228.492</td>
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
<td style="text-align: right;">0.051</td>
<td style="text-align: right;">239.117</td>
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
<td style="text-align: right;">0.048</td>
<td style="text-align: right;">241.930</td>
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
<td style="text-align: right;">0.048</td>
<td style="text-align: right;">245.367</td>
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
<td style="text-align: right;">0.083</td>
<td style="text-align: right;">259.742</td>
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
<td style="text-align: right;">0.082</td>
<td style="text-align: right;">265.055</td>
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
<td style="text-align: right;">277.086</td>
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
`develop` and neither report measures public reference products or
BAM/CRAM depth sensitivity.
