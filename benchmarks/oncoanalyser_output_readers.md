
<!-- Render from the repository root with: Rscript -e 'rmarkdown::render("benchmarks/oncoanalyser_output_readers.Rmd")' -->

# Reading Oncoanalyser outputs from S3 with DuckHTS

This report is a runnable DuckDB script wrapped in R Markdown. It
demonstrates a minimal loader plan for the public nf-core Oncoanalyser
2.3.0 example output.

The report runs in this order:

1.  start DuckDB through `duckknit` and load `httpfs` plus DuckHTS;
2.  define one S3 root macro and thin reader macros;
3.  list the S3 prefix with `glob()` and classify files by suffix;
4.  run representative reads for tables, JSON, HTS files, report blobs,
    and provenance text;
5.  summarize which reader covers each output family.

Source data:

``` text
s3://nf-core-awsmegatests/oncoanalyser/results-234fd82acc16a3beb01bf301900d83346b6ec812/
```

All file reads below use `s3://` paths. Large BAM examples use
`read_bam(...) LIMIT ...` for row previews and `read_hts_header()` for
metadata; use unbounded scans such as `count(*)` only when full
streaming is intended.

The setup chunk loads:

``` sql
LOAD httpfs;
SET s3_region='eu-west-1';
SET allow_extensions_metadata_mismatch=true;
LOAD 'build/release/duckhts.duckdb_extension';
```

## Reader macros

``` sql
CREATE OR REPLACE MACRO onco_root() AS
  's3://nf-core-awsmegatests/oncoanalyser/results-234fd82acc16a3beb01bf301900d83346b6ec812/';
```

``` sql
CREATE OR REPLACE MACRO onco_rel(path) AS onco_root() || path;
```

``` sql
CREATE OR REPLACE MACRO onco_file_suffix(path) AS
CASE
  WHEN lower(path) LIKE '%.vcf.gz.tbi' THEN '.vcf.gz.tbi'
  WHEN lower(path) LIKE '%.vcf.gz' THEN '.vcf.gz'
  WHEN lower(path) LIKE '%.bam.bai' THEN '.bam.bai'
  WHEN lower(path) LIKE '%.tsv.gz' THEN '.tsv.gz'
  WHEN lower(path) LIKE '%.layout.gz' THEN '.layout.gz'
  WHEN lower(path) LIKE '%.tar.gz' THEN '.tar.gz'
  WHEN strpos(regexp_extract(path, '[^/]+$', 0), '.') > 0
    THEN '.' || regexp_extract(lower(path), '[.]([^.]+)$', 1)
  ELSE '(none)'
END;
```

``` sql
CREATE OR REPLACE MACRO onco_file_kind(path) AS
CASE
  WHEN onco_file_suffix(path) IN ('.vcf', '.vcf.gz') THEN 'variant_vcf'
  WHEN onco_file_suffix(path) = '.vcf.gz.tbi' THEN 'variant_index'
  WHEN onco_file_suffix(path) = '.bam' THEN 'alignment_bam'
  WHEN onco_file_suffix(path) = '.bam.bai' THEN 'alignment_index'
  WHEN onco_file_suffix(path) IN ('.tsv', '.tsv.gz', '.csv', '.pcf', '.circos') THEN 'tabular_text'
  WHEN onco_file_suffix(path) = '.json' THEN 'json_document'
  WHEN onco_file_suffix(path) IN ('.html', '.png', '.pdf') THEN 'report_blob'
  WHEN onco_file_suffix(path) IN ('.log', '.out', '.err', '.run', '.sh', '.conf', '.txt', '.version', '.qc', '.metrics', '.yml') THEN 'provenance_text'
  WHEN onco_file_suffix(path) = '.layout.gz' THEN 'gzip_text'
  ELSE 'other'
END;
```

``` sql
CREATE OR REPLACE MACRO onco_read_tsv(path) AS TABLE
SELECT *
FROM read_csv(
  path,
  delim = '\t',
  header = true,
  union_by_name = true,
  filename = true,
  all_varchar = true,
  nullstr = ['', 'NA', 'null'],
  strict_mode = false,
  null_padding = true,
  ignore_errors = true
);
```

``` sql
CREATE OR REPLACE MACRO onco_read_csv(path) AS TABLE
SELECT *
FROM read_csv(
  path,
  header = true,
  union_by_name = true,
  filename = true,
  all_varchar = true,
  nullstr = ['', 'NA', 'null'],
  strict_mode = false,
  null_padding = true,
  ignore_errors = true
);
```

``` sql
CREATE OR REPLACE MACRO onco_read_json(path) AS TABLE
SELECT *
FROM read_json_auto(
  path,
  filename = true,
  union_by_name = true,
  maximum_object_size = 134217728
);
```

``` sql
CREATE OR REPLACE MACRO onco_read_vcf(path) AS TABLE
SELECT * FROM read_bcf(path);
```

``` sql
CREATE OR REPLACE MACRO onco_read_bam(path) AS TABLE
SELECT * FROM read_bam(path);
```

``` sql
CREATE OR REPLACE MACRO onco_read_blob(path) AS TABLE
SELECT
  filename,
  size,
  md5(content) AS md5,
  left(hex(content), 32) AS magic_hex,
  content
FROM read_blob(path);
```

``` sql
CREATE OR REPLACE MACRO onco_read_text_lines(path) AS TABLE
WITH docs AS (
  SELECT filename, string_split(content, chr(10)) AS lines
  FROM read_text(path)
)
SELECT
  filename,
  i AS line_no,
  regexp_replace(lines[i], chr(13) || '$', '') AS line
FROM docs, range(1, len(lines) + 1) AS r(i)
WHERE regexp_replace(lines[i], chr(13) || '$', '') <> '';
```

``` sql
CREATE OR REPLACE MACRO onco_read_gzip_lines(path) AS TABLE
SELECT filename, row_number() OVER (PARTITION BY filename) AS line_no, column0 AS line
FROM read_csv(
  path,
  header = false,
  columns = {'column0': 'VARCHAR'},
  compression = 'gzip',
  delim = '\t',
  strict_mode = false,
  ignore_errors = true,
  max_line_size = 10000000,
  filename = true
);
```

``` sql
CREATE OR REPLACE MACRO onco_read_key_value(path) AS TABLE
SELECT
  filename,
  line_no,
  trim(CASE WHEN strpos(line, chr(9)) > 0 THEN split_part(line, chr(9), 1) ELSE split_part(line, '=', 1) END) AS key,
  trim(CASE WHEN strpos(line, chr(9)) > 0 THEN split_part(line, chr(9), 2) ELSE regexp_extract(line, '^[^=]+=(.*)$', 1) END) AS value
FROM onco_read_text_lines(path)
WHERE line NOT LIKE '#%' AND line NOT LIKE '<%' AND (strpos(line, chr(9)) > 0 OR strpos(line, '=') > 0);
```

``` sql
CREATE OR REPLACE MACRO onco_read_nextflow_trace(path) AS TABLE
SELECT *
FROM read_csv(
  path,
  delim = '\t',
  header = true,
  filename = true,
  all_varchar = true,
  union_by_name = true,
  strict_mode = false,
  null_padding = true,
  ignore_errors = true
);
```

``` sql
CREATE OR REPLACE MACRO onco_read_conf_assignments(path) AS TABLE
SELECT
  filename,
  line_no,
  trim(split_part(line, '=', 1)) AS key,
  trim(regexp_replace(regexp_extract(line, '^[^=]+=(.*)$', 1), '#.*$', '')) AS value
FROM onco_read_text_lines(path)
WHERE regexp_matches(line, '^[A-Za-z0-9_.-]+[[:space:]]*=');
```

## S3 inventory without downloading objects

`glob()` lists the public S3 prefix through DuckDB/httpfs. This is key
inventory only; it does not materialize the BAM/VCF/table/report
payloads.

``` sql
CREATE OR REPLACE TEMP TABLE onco_s3_files AS
SELECT
  file AS s3_uri,
  regexp_replace(file, '^s3://nf-core-awsmegatests/oncoanalyser/results-[^/]+/', '') AS rel_path,
  onco_file_suffix(file) AS suffix,
  onco_file_kind(file) AS kind
FROM glob(onco_root() || '**/*');
```

``` sql
SELECT suffix, kind, count(*) AS files
FROM onco_s3_files
GROUP BY suffix, kind
ORDER BY files DESC, suffix;
```

    ┌─────────────┬─────────────────┬───────┐
    │   suffix    │      kind       │ files │
    │   varchar   │     varchar     │ int64 │
    ├─────────────┼─────────────────┼───────┤
    │ .png        │ report_blob     │   145 │
    │ .tsv        │ tabular_text    │    91 │
    │ .log        │ provenance_text │    35 │
    │ .sh         │ provenance_text │    35 │
    │ .err        │ provenance_text │    34 │
    │ .out        │ provenance_text │    34 │
    │ .run        │ provenance_text │    34 │
    │ .vcf.gz     │ variant_vcf     │    17 │
    │ .vcf.gz.tbi │ variant_index   │    17 │
    │ .tsv.gz     │ tabular_text    │    12 │
    │ .html       │ report_blob     │    11 │
    │ .bam        │ alignment_bam   │     9 │
    │ .csv        │ tabular_text    │     9 │
    │ .circos     │ tabular_text    │     8 │
    │ .bam.bai    │ alignment_index │     7 │
    │ .json       │ json_document   │     6 │
    │ .txt        │ provenance_text │     6 │
    │ .version    │ provenance_text │     5 │
    │ .pcf        │ tabular_text    │     3 │
    │ .conf       │ provenance_text │     2 │
    │ .layout.gz  │ gzip_text       │     2 │
    │ .qc         │ provenance_text │     2 │
    │ .metrics    │ provenance_text │     1 │
    │ .pdf        │ report_blob     │     1 │
    │ .vcf        │ variant_vcf     │     1 │
    │ .yml        │ provenance_text │     1 │
    └─────────────┴─────────────────┴───────┘
      26 rows                     3 columns

Large BAMs are visible in the inventory but are not read wholesale here:

``` sql
SELECT rel_path, suffix
FROM onco_s3_files
WHERE kind = 'alignment_bam'
ORDER BY rel_path;
```

    ┌─────────────────────────────────────────────────────┬─────────┐
    │                      rel_path                       │ suffix  │
    │                       varchar                       │ varchar │
    ├─────────────────────────────────────────────────────┼─────────┤
    │ HCC1395/alignments/dna/HCC1395.normal_dna.redux.bam │ .bam    │
    │ HCC1395/alignments/dna/HCC1395.tumor_dna.redux.bam  │ .bam    │
    │ HCC1395/alignments/rna/HCC1395.tumor_rna.md.bam     │ .bam    │
    │ HCC1395/cider/HCC1395.tumor_dna.cider.bam           │ .bam    │
    │ HCC1395/cider/HCC1395.tumor_rna.cider.bam           │ .bam    │
    │ HCC1395/esvee/HCC1395.normal_dna.esvee.prep.bam     │ .bam    │
    │ HCC1395/esvee/HCC1395.tumor_dna.esvee.prep.bam      │ .bam    │
    │ HCC1395/teal/HCC1395.normal_dna.teal.telbam.bam     │ .bam    │
    │ HCC1395/teal/HCC1395.tumor_dna.teal.telbam.bam      │ .bam    │
    └─────────────────────────────────────────────────────┴─────────┘

## Table-like outputs

`cider` in names such as `*.cider.vdj.tsv.gz` is the Oncoanalyser/CIDER
tool prefix, not a separate file format. The files themselves are
ordinary BAM, TSV/TSV.GZ, and gzipped text outputs for immune-receptor /
V(D)J calls and supporting alignments/layouts.

``` sql
SELECT 'tsv' AS reader, count(*) AS rows
FROM onco_read_tsv(onco_rel('HCC1395/bamtools/HCC1395.normal_dna/HCC1395.normal_dna.bam_metric.summary.tsv'))
UNION ALL
SELECT 'tsv.gz' AS reader, count(*) AS rows
FROM onco_read_tsv(onco_rel('HCC1395/cider/HCC1395.tumor_rna.cider.alignment_match.tsv.gz'))
UNION ALL
SELECT 'csv' AS reader, count(*) AS rows
FROM onco_read_csv(onco_rel('HCC1395/isofox/HCC1395.tumor_rna.isf.summary.csv'))
UNION ALL
SELECT 'pcf' AS reader, count(*) AS rows
FROM onco_read_tsv(onco_rel('HCC1395/amber/HCC1395.tumor_dna.amber.baf.pcf'))
UNION ALL
SELECT 'circos' AS reader, count(*) AS rows
FROM onco_read_tsv(onco_rel('HCC1395/purple/circos/HCC1395.tumor_dna.link.circos'));
```

    ┌─────────┬───────┐
    │ reader  │ rows  │
    │ varchar │ int64 │
    ├─────────┼───────┤
    │ tsv     │     1 │
    │ tsv.gz  │    12 │
    │ csv     │     1 │
    │ pcf     │  2651 │
    │ circos  │  1274 │
    └─────────┴───────┘

Key/value QC files are line-oriented text:

``` sql
SELECT regexp_extract(filename, '[^/]+$', 0) AS file, key, value
FROM onco_read_key_value(onco_rel('HCC1395/purple/HCC1395.tumor_dna.purple.qc'))
ORDER BY line_no
LIMIT 10;
```

    ┌─────────────────────────────┬───────────────────────────────┬─────────┐
    │            file             │              key              │  value  │
    │           varchar           │            varchar            │ varchar │
    ├─────────────────────────────┼───────────────────────────────┼─────────┤
    │ HCC1395.tumor_dna.purple.qc │ QCStatus                      │ PASS    │
    │ HCC1395.tumor_dna.purple.qc │ Method                        │ NORMAL  │
    │ HCC1395.tumor_dna.purple.qc │ CopyNumberSegments            │ 2743    │
    │ HCC1395.tumor_dna.purple.qc │ UnsupportedCopyNumberSegments │ 5       │
    │ HCC1395.tumor_dna.purple.qc │ Purity                        │ 1.0000  │
    │ HCC1395.tumor_dna.purple.qc │ AmberGender                   │ MALE    │
    │ HCC1395.tumor_dna.purple.qc │ CobaltGender                  │ MALE    │
    │ HCC1395.tumor_dna.purple.qc │ DeletedGenes                  │ 147     │
    │ HCC1395.tumor_dna.purple.qc │ Contamination                 │ 0.0     │
    │ HCC1395.tumor_dna.purple.qc │ GermlineAberrations           │ NONE    │
    └─────────────────────────────┴───────────────────────────────┴─────────┘
      10 rows                                                     3 columns

Picard-style metrics can be landed as text first; table splitting is a
second-stage, tool-specific view.

``` sql
SELECT regexp_extract(filename, '[^/]+$', 0) AS file, line_no, left(line, 100) AS line
FROM onco_read_text_lines(onco_rel('HCC1395/alignments/rna/HCC1395.tumor_rna.md.metrics'))
ORDER BY line_no
LIMIT 8;
```

    ┌──────────────────────────────┬─────────┬───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │             file             │ line_no │                                                   line                                                    │
    │           varchar            │  int64  │                                                  varchar                                                  │
    ├──────────────────────────────┼─────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────┤
    │ HCC1395.tumor_rna.md.metrics │       1 │ ## htsjdk.samtools.metrics.StringHeader                                                                   │
    │ HCC1395.tumor_rna.md.metrics │       2 │ # MarkDuplicates --INPUT HCC1395.tumor_rna.HCC1395__tumour_wts__2891075264__SRR892423__subsampled__l      │
    │ HCC1395.tumor_rna.md.metrics │       3 │ ## htsjdk.samtools.metrics.StringHeader                                                                   │
    │ HCC1395.tumor_rna.md.metrics │       4 │ # Started on: Mon May 11 23:28:43 GMT 2026                                                                │
    │ HCC1395.tumor_rna.md.metrics │       6 │ ## METRICS CLASS\tpicard.sam.DuplicationMetrics                                                           │
    │ HCC1395.tumor_rna.md.metrics │       7 │ LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUN │
    │ HCC1395.tumor_rna.md.metrics │       8 │ Unknown Library\t2649505\t57711259\t13597250\t3548113\t2293242\t8452952\t0\t0.162605\t177260073           │
    │ HCC1395.tumor_rna.md.metrics │      10 │ ## HISTOGRAM\tjava.lang.Double                                                                            │
    └──────────────────────────────┴─────────┴───────────────────────────────────────────────────────────────────────────────────────────────────────────┘

`*.cider.layout.gz` files are CIDER diagnostic layouts: compressed
plain-text read pileups around CIDER immune-receptor calls. They
interleave call headers such as `type: IGHV`, ANSI-coloured
nucleotide/quality strings, and `read:` alignment lines. Treat them as
provenance/debug text; use `*.cider.vdj.tsv.gz` for structured V(D)J
calls.

They are directly readable from S3 too:

``` sql
SELECT regexp_extract(filename, '[^/]+$', 0) AS file, line_no, left(line, 120) AS line
FROM onco_read_gzip_lines(onco_rel('HCC1395/cider/HCC1395.tumor_dna.cider.layout.gz'))
ORDER BY line_no
LIMIT 5;
```

    ┌───────────────────────────────────┬─────────┬────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │               file                │ line_no │                                                                line                                                                │
    │              varchar              │  int64  │                                                              varchar                                                               │
    ├───────────────────────────────────┼─────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
    │ HCC1395.tumor_dna.cider.layout.gz │       1 │                                                                 \e[31mT\e[35mC\e[35mC\e[32mA\e[33mG\e[33mG\e[33mG\e[32mA\e[32mA\e[ │
    │ HCC1395.tumor_dna.cider.layout.gz │       2 │                                                                 AAFFFKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK           │
    │ HCC1395.tumor_dna.cider.layout.gz │       3 │     read: WGS_IL_T_1__SRR7890856.37946297 1/2 151b aligned to chr16:33803043-33803193.                                             │
    │ HCC1395.tumor_dna.cider.layout.gz │       4 │                                                                   \e[35mC\e[32mA\e[33mG\e[33mG\e[33mG\e[32mA\e[32mA\e[33mG\e[33mG  │
    │ HCC1395.tumor_dna.cider.layout.gz │       5 │                                                                   AAFFFKKKKKKKKKKAKKKKKKKKKKKKKKKKKKKKKKAFKKKKKKKKKKKKKF           │
    └───────────────────────────────────┴─────────┴────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

## JSON outputs

``` sql
SELECT regexp_extract(filename, '[^/]+$', 0) AS file, count(*) AS rows
FROM onco_read_json(onco_rel('pipeline_info/params_2026-05-11_23-24-56.json'))
GROUP BY filename;
```

    ┌─────────────────────────────────┬───────┐
    │              file               │ rows  │
    │             varchar             │ int64 │
    ├─────────────────────────────────┼───────┤
    │ params_2026-05-11_23-24-56.json │     1 │
    └─────────────────────────────────┴───────┘

## HTS outputs through DuckHTS

VCF and BAM examples use the `s3://` URLs directly. `LIMIT` is enough
for a BAM row preview; do not turn this into an unbounded scan on the
80+ GB BAMs unless you actually want to stream them.

``` sql
SELECT CHROM, POS, REF, ALT, FILTER
FROM onco_read_vcf(onco_rel('HCC1395/sage/somatic/HCC1395.tumor_dna.sage.somatic.vcf.gz'))
LIMIT 5;
```

    ┌─────────┬───────┬─────────┬───────────┬──────────────────────────────────────────────┐
    │  CHROM  │  POS  │   REF   │    ALT    │                    FILTER                    │
    │ varchar │ int64 │ varchar │ varchar[] │                  varchar[]                   │
    ├─────────┼───────┼─────────┼───────────┼──────────────────────────────────────────────┤
    │ chr1    │ 13613 │ T       │ [A]       │ [maxGermlineVAF, minMQF]                     │
    │ chr1    │ 13684 │ C       │ [T]       │ [maxGermlineVAF, minMQF]                     │
    │ chr1    │ 13813 │ T       │ [G]       │ [maxGermlineRelQual, maxGermlineVAF, minMQF] │
    │ chr1    │ 13838 │ C       │ [T]       │ [maxGermlineRelQual, maxGermlineVAF, minMQF] │
    │ chr1    │ 13912 │ G       │ [A]       │ [maxGermlineVAF, minGermlineDepth, minMQF]   │
    └─────────┴───────┴─────────┴───────────┴──────────────────────────────────────────────┘

``` sql
SELECT qname, flag, rname, pos, mapq, cigar
FROM onco_read_bam(onco_rel('HCC1395/alignments/dna/HCC1395.normal_dna.redux.bam'))
LIMIT 5;
```

    ┌──────────────────────────────────┬────────┬─────────┬───────┬───────┬─────────────┐
    │              QNAME               │  FLAG  │  RNAME  │  POS  │ MAPQ  │    CIGAR    │
    │             varchar              │ uint16 │ varchar │ int64 │ int32 │   varchar   │
    ├──────────────────────────────────┼────────┼─────────┼───────┼───────┼─────────────┤
    │ WGS_IL_N_1__SRR7945633.241510432 │    163 │ chr1    │ 10048 │    60 │ 88M         │
    │ WGS_IL_N_1__SRR7945634.179122439 │    163 │ chr1    │ 10048 │    25 │ 130M1I11M8S │
    │ WGS_IL_N_1__SRR7945633.241510432 │     83 │ chr1    │ 10048 │    60 │ 88M         │
    │ WGS_IL_N_1__SRR7945634.179122439 │     83 │ chr1    │ 10048 │    25 │ 130M1I14M   │
    │ WGS_IL_N_1__SRR7945634.14921884  │     99 │ chr1    │ 10056 │     8 │ 50M2D35M4S  │
    └──────────────────────────────────┴────────┴─────────┴───────┴───────┴─────────────┘

For metadata, use header or region/index-aware peeks instead of
whole-file reads:

``` sql
SELECT record_type, id, length, key_values
FROM read_hts_header(onco_rel('HCC1395/alignments/dna/HCC1395.normal_dna.redux.bam'))
WHERE record_type IN ('HD', 'SQ')
LIMIT 6;
```

    ┌─────────────┬─────────┬───────────┬─────────────────────────┐
    │ record_type │   id    │  length   │       key_values        │
    │   varchar   │ varchar │   int64   │  map(varchar, varchar)  │
    ├─────────────┼─────────┼───────────┼─────────────────────────┤
    │ HD          │ NULL    │      NULL │ {VN=1.6, SO=coordinate} │
    │ SQ          │ chr1    │ 248956422 │ {SN=chr1, LN=248956422} │
    │ SQ          │ chr2    │ 242193529 │ {SN=chr2, LN=242193529} │
    │ SQ          │ chr3    │ 198295559 │ {SN=chr3, LN=198295559} │
    │ SQ          │ chr4    │ 190214555 │ {SN=chr4, LN=190214555} │
    │ SQ          │ chr5    │ 181538259 │ {SN=chr5, LN=181538259} │
    └─────────────┴─────────┴───────────┴─────────────────────────┘

## Reports, images, and other blobs

``` sql
SELECT
  onco_file_suffix(filename) AS suffix,
  size,
  md5,
  magic_hex
FROM onco_read_blob(onco_rel('HCC1395/linx/somatic_plots/all/HCC1395.tumor_dna.chrX.071.png'));
```

    ┌─────────┬─────────┬──────────────────────────────────┬──────────────────────────────────┐
    │ suffix  │  size   │               md5                │            magic_hex             │
    │ varchar │  int64  │             varchar              │             varchar              │
    ├─────────┼─────────┼──────────────────────────────────┼──────────────────────────────────┤
    │ .png    │ 1021812 │ 31d14fa20688686f2125be996f139dbb │ 89504E470D0A1A0A0000000D49484452 │
    └─────────┴─────────┴──────────────────────────────────┴──────────────────────────────────┘

HTML reports can be read as text when metadata is enough:

``` sql
SELECT
  regexp_extract(filename, '[^/]+$', 0) AS file,
  regexp_extract(content, '<title>([^<]+)</title>', 1) AS title
FROM read_text(onco_rel('pipeline_info/execution_report_2026-05-11_23-23-53.html'));
```

    ┌───────────────────────────────────────────┬─────────────────────────────────────────────┐
    │                   file                    │                    title                    │
    │                  varchar                  │                   varchar                   │
    ├───────────────────────────────────────────┼─────────────────────────────────────────────┤
    │ execution_report_2026-05-11_23-23-53.html │ [cheesy_babbage_4] Nextflow Workflow Report │
    └───────────────────────────────────────────┴─────────────────────────────────────────────┘

## Nextflow provenance and config-like files

``` sql
SELECT
  regexp_extract(filename, '[^/]+$', 0) AS trace_file,
  count(*) AS tasks,
  count(*) FILTER (WHERE status = 'COMPLETED') AS completed,
  count(*) FILTER (WHERE status = 'CACHED') AS cached
FROM onco_read_nextflow_trace(onco_rel('pipeline_info/execution_trace_2026-05-11_23-23-53.txt'))
GROUP BY filename;
```

    ┌─────────────────────────────────────────┬───────┬───────────┬────────┐
    │               trace_file                │ tasks │ completed │ cached │
    │                 varchar                 │ int64 │   int64   │ int64  │
    ├─────────────────────────────────────────┼───────┼───────────┼────────┤
    │ execution_trace_2026-05-11_23-23-53.txt │    53 │        44 │      9 │
    └─────────────────────────────────────────┴───────┴───────────┴────────┘

``` sql
SELECT onco_file_suffix(filename) AS suffix, count(*) AS lines
FROM onco_read_text_lines([
  onco_rel('HCC1395/logs/HCC1395.teal_pipeline.command.log'),
  onco_rel('HCC1395/logs/HCC1395.teal_pipeline.command.out'),
  onco_rel('HCC1395/logs/HCC1395.teal_pipeline.command.err'),
  onco_rel('HCC1395/logs/HCC1395.teal_pipeline.command.run'),
  onco_rel('HCC1395/logs/HCC1395.teal_pipeline.command.sh')
])
GROUP BY suffix
ORDER BY lines DESC;
```

    ┌─────────┬───────┐
    │ suffix  │ lines │
    │ varchar │ int64 │
    ├─────────┼───────┤
    │ .log    │ 10208 │
    │ .out    │ 10206 │
    │ .run    │   337 │
    │ .sh     │    20 │
    │ .err    │     2 │
    └─────────┴───────┘

``` sql
SELECT regexp_extract(filename, '[^/]+$', 0) AS file, key, value
FROM onco_read_conf_assignments(onco_rel('HCC1395/purple/circos/HCC1395.tumor_dna.circos.conf'))
ORDER BY line_no
LIMIT 10;
```

    ┌───────────────────────────────┬─────────────────────────────┬─────────────────────────────────────────┐
    │             file              │             key             │                  value                  │
    │            varchar            │           varchar           │                 varchar                 │
    ├───────────────────────────────┼─────────────────────────────┼─────────────────────────────────────────┤
    │ HCC1395.tumor_dna.circos.conf │ show_ticks                  │ no                                      │
    │ HCC1395.tumor_dna.circos.conf │ show_tick_labels            │ no                                      │
    │ HCC1395.tumor_dna.circos.conf │ karyotype                   │ data/karyotype/karyotype.human.hg38.txt │
    │ HCC1395.tumor_dna.circos.conf │ chromosomes_units           │ 1000000                                 │
    │ HCC1395.tumor_dna.circos.conf │ chromosomes_display_default │ yes                                     │
    │ HCC1395.tumor_dna.circos.conf │ chromosomes                 │ -hsZ                                    │
    └───────────────────────────────┴─────────────────────────────┴─────────────────────────────────────────┘

## Coverage summary

``` sql
SELECT *
FROM (VALUES
  ('BAM / BAI', 'read_bam(), read_hts_header()', 'remote S3 works; avoid whole-file scans on huge BAMs unless intended'),
  ('VCF / VCF.GZ / TBI', 'read_bcf()', 'remote S3 works through DuckHTS/htslib'),
  ('TSV / TSV.GZ / CSV', 'read_csv()', 'remote S3 works through DuckDB/httpfs'),
  ('JSON', 'read_json_auto()', 'remote S3 works; ORANGE may need larger maximum_object_size'),
  ('PCF / QC / VERSION', 'read_csv() or key/value macro', 'covered by simple SQL macros'),
  ('CIRCOS / CONF', 'read_csv() or text assignment macro', 'covered by simple SQL macros'),
  ('LOG / OUT / ERR / RUN / SH', 'read_text() line macro', 'covered as provenance text'),
  ('LAYOUT.GZ', 'read_csv(... compression=gzip) line macro', 'covered as gzip text'),
  ('PNG / PDF / HTML', 'read_blob(); read_text() for HTML', 'covered as artefact blobs plus metadata')
) AS t(output_family, reader, status);
```

    ┌────────────────────────────┬───────────────────────────────────────────┬──────────────────────────────────────────────────────────────────────┐
    │       output_family        │                  reader                   │                                status                                │
    │          varchar           │                  varchar                  │                               varchar                                │
    ├────────────────────────────┼───────────────────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
    │ BAM / BAI                  │ read_bam(), read_hts_header()             │ remote S3 works; avoid whole-file scans on huge BAMs unless intended │
    │ VCF / VCF.GZ / TBI         │ read_bcf()                                │ remote S3 works through DuckHTS/htslib                               │
    │ TSV / TSV.GZ / CSV         │ read_csv()                                │ remote S3 works through DuckDB/httpfs                                │
    │ JSON                       │ read_json_auto()                          │ remote S3 works; ORANGE may need larger maximum_object_size          │
    │ PCF / QC / VERSION         │ read_csv() or key/value macro             │ covered by simple SQL macros                                         │
    │ CIRCOS / CONF              │ read_csv() or text assignment macro       │ covered by simple SQL macros                                         │
    │ LOG / OUT / ERR / RUN / SH │ read_text() line macro                    │ covered as provenance text                                           │
    │ LAYOUT.GZ                  │ read_csv(... compression=gzip) line macro │ covered as gzip text                                                 │
    │ PNG / PDF / HTML           │ read_blob(); read_text() for HTML         │ covered as artefact blobs plus metadata                              │
    └────────────────────────────┴───────────────────────────────────────────┴──────────────────────────────────────────────────────────────────────┘

Conclusion: no new C reader is needed for the initial Oncoanalyser
loader. A small, organized SQL macro pack on top of DuckHTS plus
DuckDB/httpfs covers the observed output families directly on S3.
