#!/usr/bin/env python3
"""Exercise the release extension under a DuckDB 2.0 Python wheel."""

from pathlib import Path
import re
import shutil
import tempfile

import duckdb


root = Path(__file__).resolve().parents[2]
extension = root / "build/release/duckhts.duckdb_extension"
requirements = (root / "test/scripts/duckdb-v2-requirements.txt").read_text(encoding="utf-8")
pinned = re.search(r"^duckdb==(\S+)", requirements, re.M).group(1)
if duckdb.__version__ != pinned:
    raise SystemExit(f"Expected the pinned DuckDB {pinned}, found {duckdb.__version__}")
if not extension.is_file():
    raise SystemExit(f"Build the release extension first: {extension}")

con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
con.execute(f"LOAD '{extension}'")


def check(sql, expected):
    actual = con.execute(sql).fetchall()
    if actual != expected:
        raise AssertionError(f"{sql}\nexpected {expected!r}; got {actual!r}")


vcf = root / "test/data/mapping_number_families.vcf"
check(f"SELECT count(*) FROM read_bcf('{vcf}')", [(2,)])
con.execute("""
CREATE TEMP TABLE somalier_v2_panel AS
SELECT * FROM (VALUES
  ('GRCh38', 0::UBIGINT, 'chr1', 100::UBIGINT, 'A', 'C'),
  ('GRCh38', 1::UBIGINT, 'chr1', 200::UBIGINT, 'A', 'G'),
  ('GRCh38', 2::UBIGINT, 'chr1', 300::UBIGINT, 'C', 'T')
) AS p(assembly, site_index, region, position, allele_a, allele_b)
""")
check(f"""
SELECT sample_id, site_index, a, b, other, status
FROM duckhts_somalier_vcf_counts('{vcf}', 'somalier_v2_panel')
ORDER BY source_sample_index, site_index
""", [
    ("S1", 0, 9, 3, 0, "measured"),
    ("S1", 1, 5, 9, 4, "measured"),
    ("S1", 2, None, None, None, "unavailable_no_record"),
    ("S2", 0, 0, 20, 0, "measured"),
    ("S2", 1, 0, 10, 5, "measured"),
    ("S2", 2, None, None, None, "unavailable_no_record"),
])


# A DuckHTS scalar that reports an error must raise its own error. DuckDB 2.0
# builds before duckdb/duckdb#25664 turned it into an INTERNAL error, because
# C API v1 cannot mark a function fallible. Each case gets its own connection:
# on such a build the error also invalidated the database.
def check_error(sql, message):
    case = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    case.execute(f"LOAD '{extension}'")
    try:
        case.execute(sql).fetchall()
    except duckdb.InvalidInputException as error:
        if message not in str(error):
            raise AssertionError(f"{sql}\nexpected {message!r}; got {error}") from error
    else:
        raise AssertionError(f"{sql}\nexpected an InvalidInputException")
    finally:
        case.close()


check_error("SELECT duckhts_somalier_spacing([100, 110]::UBIGINT[], 0)", "duckhts_somalier_spacing")
check_error("SELECT duckhts_cgranges_index('no_such_index')", "unknown index name")


# DuckDB 2.0 binds a prepared statement again at every EXECUTE. A writer that
# worked at bind wrote at PREPARE, and the first EXECUTE then found its own
# output (https://github.com/RGenomicsETL/duckhts/issues/358). R prepares and
# then executes every query, so this is the path DBI users take.
work = Path(tempfile.mkdtemp(prefix="duckhts_v2_writers_"))
try:
    bam = work / "range.bam"
    shutil.copy(root / "test/data/range.bam", bam)
    writers = {
        "bam_index": (f"bam_index('{bam}', index_path := '{work}/range.bai', threads := 1)",
                      work / "range.bai"),
        "tabix_index": (f"tabix_index('{root}/test/data/bcf_scalar_counts.vcf.gz', index_path := '{work}/scalar.tbi', threads := 1)",
                        work / "scalar.tbi"),
        "fasta_index": (f"fasta_index('{root}/test/data/ce.fa', index_path := '{work}/ce.fai')",
                        work / "ce.fai"),
        "duckhts_samtools_idxstats": (f"duckhts_samtools_idxstats('{bam}', output := '{work}/idxstats.txt')",
                                      work / "idxstats.txt"),
        "duckhts_mosdepth": (f"duckhts_mosdepth('{work}/md', '{bam}', fast_mode := TRUE, no_per_base := TRUE)",
                             work / "md.mosdepth.summary.txt"),
    }
    for number, (name, (call, output)) in enumerate(writers.items()):
        con.execute(f"PREPARE writer_{number} AS SELECT success FROM {call}")
        if output.exists():
            raise AssertionError(f"{name}: PREPARE wrote {output}")
        check(f"EXECUTE writer_{number}", [(True,)])
        if not output.exists():
            raise AssertionError(f"{name}: EXECUTE did not write {output}")
finally:
    shutil.rmtree(work, ignore_errors=True)
print(f"DuckDB {duckdb.__version__}: LOAD, read_bcf, Somalier, scalar errors and prepared writers passed")
