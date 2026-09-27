#!/usr/bin/env python3
"""Build equivalent feature stores and reject mismatched keyed query results."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import re
import subprocess
import sys
import tempfile
import time

import duckdb
import gffbase
import pyarrow as pa

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "r/duckhtsbench/inst/benchmark_registry.tsv"
SQL = ROOT / "scripts/gffbase_featuredb.sql"
FIELDS = ("anchor", "descendant_id", "seqid", "source", "featuretype", "start",
          "end", "score", "strand", "frame", "file_order", "depth")
PARITY_KEYS = ("features", "anchors", "windows", "anchors_sha256", "windows_sha256",
               "descendants", "descendants_sha256", "region_hits", "region_sha256")


def digest_rows(table):
    rows = sorted(zip(*(table.column(field).to_pylist() for field in FIELDS)))
    digest = hashlib.sha256()
    for row in rows:
        digest.update(json.dumps(row, separators=(",", ":"), ensure_ascii=False).encode())
        digest.update(b"\n")
    return len(rows), digest.hexdigest()


def input_identity(artifact_id):
    with REGISTRY.open(newline="") as source:
        row = next((r for r in csv.DictReader(source, delimiter="\t") if r["id"] == artifact_id), None)
    if row is None:
        raise ValueError(f"unregistered input: {artifact_id}")
    cache = os.environ.get("DUCKHTS_CACHE_DIR")
    if not cache:
        cache = str(Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "duckhts")
    path = Path(cache).expanduser() / row["cache_relpath"]
    if not path.is_file():
        raise FileNotFoundError(f"stage registered input first: {artifact_id}: {path}")
    identity = dict(part.split("=", 1) for part in row["supplier_identity"].split(";") if part)
    if path.stat().st_size != int(identity["bytes"]):
        raise ValueError(f"incorrect size: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != identity["sha256"]:
        raise ValueError(f"incorrect SHA-256: {path}")
    return path, digest


def resolve_duplicate_ids(db):
    """Assign GFFBase create_unique IDs to rows whose names can collide."""
    rows = db.execute(r"""
        WITH f AS (
          SELECT rid, id, regexp_replace(id, '(_[0-9]+)+$', '') AS stem
          FROM features WHERE NOT starts_with(id, 'row:')
        ), duplicated AS (
          SELECT DISTINCT stem FROM f WHERE id IN (SELECT id FROM f GROUP BY id HAVING count(*) > 1)
        )
        SELECT rid, id FROM f WHERE stem IN (SELECT stem FROM duplicated) ORDER BY rid
    """).fetchall()
    taken, counters, renames = set(), {}, []
    for rid, raw in rows:
        name = raw
        while name in taken:
            counters[raw] = counters.get(raw, 0) + 1
            name = f"{raw}_{counters[raw]}"
        taken.add(name)
        if name != raw:
            renames.append((rid, name))
    if renames:
        # MANE renames ~185k repeated CDS IDs: hand them over as Arrow (binding
        # Python lists as parameters costs seconds) and rebuild the sorted table.
        rids, ids = zip(*renames)
        db.register("renames", pa.table({"rid": pa.array(rids, pa.int64()), "id": ids}))
        db.execute("""
            CREATE OR REPLACE TABLE features AS
            SELECT f.* REPLACE (coalesce(r.id, f.id) AS id)
            FROM features f LEFT JOIN renames r USING (rid)
            ORDER BY seqid, start, "end"
        """)
        db.unregister("renames")


def worker(engine, input_path, database, extension, threads, limit):
    if database.exists():
        database.unlink()
    started = time.perf_counter()
    if engine == "DuckHTS":
        db = duckdb.connect(str(database), config={
            "allow_unsigned_extensions": "true", "memory_limit": limit, "threads": str(threads)
        })
        db.execute("LOAD " + "'" + str(extension).replace("'", "''") + "'")
        sql = SQL.read_text().replace("{input}", str(input_path).replace("'", "''"))
        before, after = sql.split("-- @resolve-duplicate-ids", 1)
        for statement in before.split(";\n"):
            if statement.strip():
                db.execute(statement)
        resolve_duplicate_ids(db)
        for statement in after.split(";\n"):
            if statement.strip() and not all(line.startswith("--") for line in statement.strip().splitlines()):
                db.execute(statement)
    else:
        db = gffbase.create_db(str(input_path), str(database), force=True,
                               merge_strategy="create_unique", force_gff=True,
                               pragmas={"threads": threads}, validation=None)
    build_s = time.perf_counter() - started

    genes = [r[0] for r in db.execute(
        "SELECT id FROM features WHERE featuretype = 'gene' ORDER BY id").fetchall()]
    random.Random(1).shuffle(genes)
    anchors = genes[:5000]
    regions = db.execute(
        'SELECT seqid, start, "end" FROM features WHERE featuretype = \'gene\' ORDER BY id').fetchall()
    random.Random(2).shuffle(regions)
    regions = regions[:1000]
    started = time.perf_counter()
    if engine == "DuckHTS":
        db.execute("CREATE TEMP TABLE anchors AS SELECT unnest($1::VARCHAR[]) AS id", [anchors])
        descendants = db.execute("""
            SELECT a.id AS anchor, f.id AS descendant_id, f.seqid, f.source,
                   f.featuretype, f.start, f."end",
                   coalesce(cast(f.score AS VARCHAR), '.') AS score,
                   f.strand, f.frame, f.rid AS file_order, c.level AS depth
            FROM anchors a JOIN closure c ON c.ancestor = a.id
            JOIN features f ON f.id = c.descendant
        """).fetch_arrow_table()
    else:
        descendants = db.children_batched(anchors, format="arrow")
    descendant_s = time.perf_counter() - started
    count, keyed_digest = digest_rows(descendants)

    started = time.perf_counter()
    if engine == "DuckHTS":
        region_ids = [[row[0] for row in db.execute(
            'SELECT id FROM features WHERE seqid = ? AND start <= ? AND "end" >= ? ORDER BY id',
            [seqid, end, start]).fetchall()] for seqid, start, end in regions]
    else:
        region_ids = [sorted(feature.id for feature in db.region(
            seqid=seqid, start=start, end=end)) for seqid, start, end in regions]
    region_s = time.perf_counter() - started
    output = {
        "features": db.execute("SELECT count(*) FROM features").fetchone()[0],
        "anchors": len(anchors), "windows": len(regions),
        "anchors_sha256": hashlib.sha256(json.dumps(anchors).encode()).hexdigest(),
        "windows_sha256": hashlib.sha256(json.dumps(regions).encode()).hexdigest(),
        "descendants": count, "descendants_sha256": keyed_digest,
        "region_hits": sum(map(len, region_ids)),
        "region_sha256": hashlib.sha256(json.dumps(region_ids).encode()).hexdigest(),
        "build_s": build_s, "descendant_s": descendant_s, "region_s": region_s,
    }
    db.close()
    output["db_bytes"] = database.stat().st_size
    print(json.dumps(output))


def measured(engine, input_path, extension, threads, limit, database, scope):
    command = ["/usr/bin/time", "-v", "-o", str(database) + ".time",
               sys.executable, str(Path(__file__).resolve()), "--worker", engine,
               "--input", str(input_path), "--database", str(database),
               "--extension", str(extension), "--threads", str(threads), "--limit", limit]
    if scope and engine == "gffbase":
        command = ["systemd-run", "--quiet", "--scope", "-p", "MemoryMax=45G"] + command
    output = subprocess.run(command, capture_output=True, text=True)
    if output.returncode:
        raise RuntimeError(f"{engine} failed: {output.stdout}\n{output.stderr}")
    result = json.loads(output.stdout)
    timing = Path(str(database) + ".time").read_text()
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", timing)
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", timing)
    if not rss or not wall:
        raise RuntimeError(f"missing /usr/bin/time measurements: {timing}")
    units = [float(part) for part in wall.group(1).split(":")]
    result["wall_s"] = sum(unit * 60 ** index for index, unit in enumerate(reversed(units)))
    result["peak_rss_kb"] = int(rss.group(1))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("DuckHTS", "gffbase"))
    parser.add_argument("--input", type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--extension", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--limit", default="8GB")
    parser.add_argument("--mane-passes", type=int, default=3)
    parser.add_argument("--gencode-passes", type=int, default=1)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "benchmarks/data")
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.input, args.database, args.extension, args.threads, args.limit)
        return
    if args.threads < 1 or args.mane_passes < 3 or args.gencode_passes < 1:
        parser.error("threads must be positive; MANE needs at least 3 passes; GENCODE needs at least 1")
    if gffbase.__version__ != "0.2.1" or duckdb.__version__ != "1.5.2":
        raise RuntimeError("this comparison pins gffbase 0.2.1 and DuckDB 1.5.2")
    if not args.extension.is_file():
        raise FileNotFoundError(args.extension)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    extension_sha = hashlib.sha256(args.extension.read_bytes()).hexdigest()
    cpu = next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")), platform.processor())
    host = platform.node()
    checks = []
    timings = []
    with tempfile.TemporaryDirectory(prefix="gffbase-featuredb-") as temp:
        for artifact, passes in (("mane_v15_ensembl_gff3", args.mane_passes),
                                 ("gencode_v49_basic_gff3", args.gencode_passes)):
            input_path, input_sha = input_identity(artifact)
            for pass_number in range(1, passes + 1):
                results = {}
                for engine in ("DuckHTS", "gffbase"):
                    database = Path(temp) / f"{artifact}-{pass_number}-{engine}.duckdb"
                    result = measured(engine, input_path, args.extension, args.threads,
                                      args.limit, database, artifact == "gencode_v49_basic_gff3")
                    results[engine] = result
                    database.unlink()
                duck, base = results["DuckHTS"], results["gffbase"]
                disagreements = {key: [duck[key], base[key]] for key in PARITY_KEYS
                                 if duck[key] != base[key]}
                if disagreements:
                    raise ValueError(f"parity failure for {artifact} pass {pass_number}: "
                                     + json.dumps(disagreements))
                checks.append(dict(dataset=artifact, pass_number=pass_number,
                                   features=duck["features"], anchors=duck["anchors"],
                                   windows=duck["windows"], descendants=duck["descendants"],
                                   region_hits=duck["region_hits"],
                                   descendants_sha256=duck["descendants_sha256"],
                                   region_sha256=duck["region_sha256"], parity="PASS"))
                for engine, result in results.items():
                    timings.append(dict(dataset=artifact, pass_number=pass_number, engine=engine,
                                        input_sha256=input_sha, revision=revision,
                                        extension_sha256=extension_sha, host=host,
                                        os=platform.platform(), cpu=cpu, threads=args.threads,
                                        duckdb_version=duckdb.__version__,
                                        gffbase_version=gffbase.__version__,
                                        duckhts_memory_limit=args.limit,
                                        gffbase_memory_cap="45G" if artifact == "gencode_v49_basic_gff3" else "none",
                                        validation="not run",
                                        **{key: result[key] for key in PARITY_KEYS},
                                        **{key: result[key] for key in (
                                            "build_s", "wall_s", "peak_rss_kb", "db_bytes",
                                            "descendant_s", "region_s")}))
                print(f"{artifact} pass {pass_number}/{passes}: parity PASS", flush=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("gffbase_featuredb_parity.csv", checks),
                           ("gffbase_featuredb_timings.csv", timings)):
        with (args.out_dir / filename).open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
