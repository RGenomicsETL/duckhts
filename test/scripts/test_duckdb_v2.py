#!/usr/bin/env python3
"""Exercise the release extension under a DuckDB 2.0 Python wheel."""

from pathlib import Path
import re

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
print(f"DuckDB {duckdb.__version__}: LOAD, read_bcf, Somalier and scalar errors passed")
