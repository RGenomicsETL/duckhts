# DuckHTS

A [DuckDB](https://duckdb.org/) extension for reading high-throughput sequencing (HTS) file formats using [htslib](https://github.com/samtools/htslib).

Query VCF, BCF, BAM, CRAM, FASTA, FASTQ, GTF, GFF, and tabix-indexed files directly from SQL.

## Functions

| Function | Description | Schema |
|---|---|---|
| `read_bcf(path, [region, tidy_format])` | Read VCF/BCF files | CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\_\*, FORMAT\_\* |
| `read_bam(path, [region])` | Read SAM/BAM/CRAM files | QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL |
| `read_fasta(path)` | Read FASTA files | NAME, DESCRIPTION, SEQUENCE |
| `read_fastq(path)` | Read FASTQ files | NAME, DESCRIPTION, SEQUENCE, QUALITY |
| `read_gff(path, [region])` | Read GFF3 files | seqname, source, feature, start, end, score, strand, frame, attributes |
| `read_gtf(path, [region])` | Read GTF files | seqname, source, feature, start, end, score, strand, frame, attributes |
| `read_tabix(path, [region])` | Read any tabix-indexed file | column0, column1, … (auto-detected) |

## Examples

```sql
LOAD 'duckhts';

-- Read a VCF file
SELECT CHROM, POS, REF, ALT, QUAL
FROM read_bcf('variants.vcf.gz')
WHERE CHROM = '1' AND POS > 1000000;

-- Region query on an indexed VCF (requires .tbi or .csi index)
SELECT * FROM read_bcf('variants.bcf', region := 'chr1:1000000-2000000');

-- Read a BAM file
SELECT QNAME, FLAG, RNAME, POS, MAPQ, CIGAR
FROM read_bam('alignments.bam')
WHERE FLAG & 4 = 0;  -- mapped reads only

-- Region query on an indexed BAM
SELECT count(*) FROM read_bam('alignments.bam', region := 'chr1:1-1000000');

-- Read FASTA sequences
SELECT NAME, length(SEQUENCE) as seq_length
FROM read_fasta('reference.fa');

-- Read FASTQ and compute average quality
SELECT NAME, length(SEQUENCE) as read_length
FROM read_fastq('reads.fq.gz');

-- Read GFF3 annotations
SELECT seqname, feature, start, "end", attributes
FROM read_gff('annotations.gff3')
WHERE feature = 'gene';

-- Read a tabix-indexed BED file
SELECT * FROM read_tabix('regions.bed.gz', region := 'chr1:1-1000000');
```

## Building

### Prerequisites

- C compiler (GCC or Clang)
- CMake ≥ 3.5
- Make
- Python 3 + venv
- Git
- htslib build dependencies: zlib, libbz2, liblzma, libdeflate, libcurl, libcrypto (OpenSSL)

On Debian/Ubuntu:
```bash
sudo apt install build-essential cmake python3 python3-venv git \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-openssl-dev libssl-dev
```

On macOS:
```bash
brew install cmake htslib xz libdeflate
```

### Vendor htslib

```bash
./scripts/vendor_htslib.sh
```

This downloads and verifies htslib 1.23 into `third_party/htslib/`.

### Build

```bash
make configure    # one-time setup (Python venv, platform detection)
make release      # build optimised extension
```

The build uses CMake's `ExternalProject` to run htslib's own `./configure && make` automatically — htslib's configure detects all available libraries (zlib, bz2, lzma, libdeflate, libcurl, OpenSSL) and enables them.

The extension binary is written to `build/release/duckhts.duckdb_extension`.

### Debug build

```bash
make debug
```

## Loading

```sql
-- Unsigned extensions must be loaded with -unsigned flag:
-- duckdb -unsigned

LOAD '/path/to/duckhts.duckdb_extension';
```

## Testing

SQL tests live in `test/sql/` using DuckDB's SQLLogicTest format.

Before running tests for the first time, prepare the indexed test data:

```bash
./test/scripts/prepare_test_data.sh   # requires samtools, bcftools, bgzip, tabix
```

This copies files from the vendored htslib test suite into `test/data/` and
builds the required indexes (BAI, CSI, TBI, FAI) so region queries work.

Then run:

```bash
make test_release
```

## Project Structure

```
src/
  duckhts.c          # Extension entry point
  bcf_reader.c       # VCF/BCF reader (read_bcf)
  bam_reader.c       # SAM/BAM/CRAM reader (read_bam)
  seq_reader.c       # FASTA/FASTQ reader (read_fasta, read_fastq)
  tabix_reader.c     # Tabix/GTF/GFF reader (read_tabix, read_gtf, read_gff)
  vep_parser.c       # VEP/CSQ annotation parser
  include/
    vcf_types.h
    vep_parser.h
third_party/
  htslib/            # Vendored htslib 1.23 (built automatically)
test/
  sql/               # SQL logic tests
duckdb_capi/
  duckdb.h           # DuckDB C API headers
  duckdb_extension.h
```

## License

MIT