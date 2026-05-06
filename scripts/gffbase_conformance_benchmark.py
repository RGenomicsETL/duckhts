#!/usr/bin/env python3
"""Benchmark DuckHTS read_gff against GFFBase for GFF3 conformance.

The conformance cases are adapted from GFFBase's Apache-2.0 NCBI/GFF3
compliance tests (Copyright 2026 Kuan-Hao Chao) at:

    https://github.com/Kuanhao-Chao/gffbase.git
    commit 78714cf30a9d799eab544e00a79a4da9754987ca

This script is intentionally skeptical about GFFBase as well as DuckHTS: it
records whether GFFBase is using its Rust or Python parser, checks Rust/Python
parity when both are available, and then compares DuckHTS default and
``strict := true`` reads against the observed GFFBase strict-parser behavior.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator


GFFBASE_UPSTREAM_URL = "https://github.com/Kuanhao-Chao/gffbase.git"
GFFBASE_UPSTREAM_COMMIT = "78714cf30a9d799eab544e00a79a4da9754987ca"


@dataclass(frozen=True)
class ConformanceCase:
    name: str
    text: bytes
    focus: str
    expected_status: str
    expected_kind: str = ""


GOOD = b"chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=a\n"
CONFORMANCE_CASES: tuple[ConformanceCase, ...] = (
    ConformanceCase(
        "valid_minimal",
        GOOD,
        "well-formed 9-column GFF3 feature",
        "accepted",
    ),
    ConformanceCase(
        "too_few_fields_7",
        b"chr1\tsrc\texon\t1\t10\t.\t+\n",
        "wrong column count: 7 fields",
        "rejected",
        "TooFewFields",
    ),
    ConformanceCase(
        "too_few_fields_8",
        b"chr1\tsrc\texon\t1\t10\t.\t+\t.\n",
        "wrong column count: 8 fields without attributes",
        "rejected",
        "TooFewFields",
    ),
    ConformanceCase(
        "non_integer_start",
        b"chr1\tsrc\texon\tabc\t10\t.\t+\t.\tID=x\n",
        "start coordinate must be an integer",
        "rejected",
        "InvalidCoordinate",
    ),
    ConformanceCase(
        "non_integer_end",
        b"chr1\tsrc\texon\t1\txyz\t.\t+\t.\tID=x\n",
        "end coordinate must be an integer",
        "rejected",
        "InvalidCoordinate",
    ),
    ConformanceCase(
        "negative_start",
        b"chr1\tsrc\texon\t-5\t10\t.\t+\t.\tID=x\n",
        "start coordinate must be >= 1",
        "rejected",
        "InvalidCoordinate",
    ),
    ConformanceCase(
        "zero_start",
        b"chr1\tsrc\texon\t0\t10\t.\t+\t.\tID=x\n",
        "start coordinate must be 1-based",
        "rejected",
        "InvalidCoordinate",
    ),
    ConformanceCase(
        "end_less_than_start",
        b"chr1\tsrc\texon\t100\t10\t.\t+\t.\tID=x\n",
        "end coordinate must be >= start",
        "rejected",
        "InvalidCoordinate",
    ),
    ConformanceCase(
        "invalid_strand_char",
        b"chr1\tsrc\texon\t1\t10\t.\t@\t.\tID=x\n",
        "strand must be one of +, -, ., ?",
        "rejected",
        "InvalidStrand",
    ),
    ConformanceCase(
        "multichar_strand",
        b"chr1\tsrc\texon\t1\t10\t.\t+-\t.\tID=x\n",
        "strand must be a single character",
        "rejected",
        "InvalidStrand",
    ),
    ConformanceCase(
        "strand_question_mark",
        b"chr1\tsrc\texon\t1\t10\t.\t?\t.\tID=x\n",
        "question-mark strand is accepted by NCBI GFF3 guidance",
        "accepted",
    ),
    ConformanceCase(
        "invalid_phase_value",
        b"chr1\tsrc\texon\t1\t10\t.\t+\t5\tID=x\n",
        "phase must be 0, 1, 2, or .",
        "rejected",
        "InvalidPhase",
    ),
    ConformanceCase(
        "cds_phase_dot",
        b"chr1\tsrc\tCDS\t1\t10\t.\t+\t.\tID=c1\n",
        "CDS phase must be concrete",
        "rejected",
        "InvalidPhase",
    ),
    ConformanceCase(
        "cds_phase_concrete",
        b"chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1\n",
        "concrete CDS phase is accepted",
        "accepted",
    ),
    ConformanceCase(
        "empty_seqid",
        b"\tsrc\texon\t1\t10\t.\t+\t.\tID=x\n",
        "seqid must not be empty",
        "rejected",
        "EmptySeqid",
    ),
    ConformanceCase(
        "empty_featuretype",
        b"chr1\tsrc\t\t1\t10\t.\t+\t.\tID=x\n",
        "feature type must not be empty",
        "rejected",
        "EmptyFeaturetype",
    ),
    ConformanceCase(
        "featuretype_whitespace",
        b"chr1\tsrc\texon foo\t1\t10\t.\t+\t.\tID=x\n",
        "feature type must not contain whitespace",
        "rejected",
        "InvalidFeaturetype",
    ),
    ConformanceCase(
        "invalid_score",
        b"chr1\tsrc\texon\t1\t10\tabc\t+\t.\tID=x\n",
        "score must be numeric or .",
        "rejected",
        "InvalidScore",
    ),
    ConformanceCase(
        "truncated_final_line",
        b"chr1\tsrc\texon\t1",
        "truncated final line must still have 9 fields",
        "rejected",
        "TooFewFields",
    ),
    ConformanceCase(
        "directives_only",
        b"##gff-version 3\n# comment\n##sequence-region chr1 1 10\n",
        "directives/comments without features are accepted as zero records",
        "accepted",
    ),
    ConformanceCase(
        "all_dots_row",
        b"chr1\tsrc\tregion\t.\t.\t.\t.\t.\tID=r\n",
        "real-world region rows may use . for unknown start/end",
        "accepted",
    ),
    ConformanceCase(
        "dbxref_multivalue",
        b"chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g;Dbxref=GeneID:1,HGNC:HGNC:1\n",
        "GFF3 comma-separated Dbxref values are multiple attribute values",
        "accepted",
    ),
    ConformanceCase(
        "garbage_attribute",
        b"chr1\tsrc\texon\t1\t10\t.\t+\t.\tjustsomegarbage\n",
        "attribute field must be key=value pairs or .",
        "rejected",
        "InvalidAttribute",
    ),
    ConformanceCase(
        "dot_attribute",
        b"chr1\tsrc\texon\t1\t10\t.\t+\t.\t.\n",
        "dot attribute field is accepted",
        "accepted",
    ),
    ConformanceCase(
        "percent_decode_attribute",
        b"chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g;Note=hello%20world\n",
        "GFF3 attribute values are URL-decoded by GFFBase",
        "accepted",
    ),
)


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def sql_string(value: str | Path) -> str:
    return str(value).replace("'", "''")


def git_rev(path: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return ""


def collect_system_specs() -> dict[str, Any]:
    specs: dict[str, Any] = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "cpu_logical": os.cpu_count(),
    }
    try:
        specs["cpu_affinity"] = len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except Exception:
        specs["cpu_affinity"] = None

    cpu_model = ""
    cpu_mhz = ""
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        try:
            for line in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
                if not cpu_model and (line.startswith("model name") or line.startswith("Hardware") or line.startswith("Processor")):
                    cpu_model = line.split(":", 1)[1].strip() if ":" in line else line.strip()
                if not cpu_mhz and line.startswith("cpu MHz") and ":" in line:
                    cpu_mhz = line.split(":", 1)[1].strip()
                if cpu_model and cpu_mhz:
                    break
        except Exception:
            pass
    specs["cpu_model"] = cpu_model or platform.processor()
    specs["cpu_mhz_first_core"] = cpu_mhz

    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        try:
            for line in meminfo.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    specs["mem_total_gib"] = round(kb / (1024 * 1024), 2)
                    specs["mem_total_kib"] = kb
                    break
        except Exception:
            pass

    os_release = Path("/etc/os-release")
    if os_release.exists():
        try:
            for line in os_release.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("PRETTY_NAME="):
                    specs["os_pretty_name"] = line.split("=", 1)[1].strip().strip('"')
                    break
        except Exception:
            pass
    return specs


@contextlib.contextmanager
def suppress_c_stderr() -> Iterator[None]:
    """Temporarily silence C-library stderr noise from optional tabix probes."""
    saved = os.dup(2)
    try:
        with open(os.devnull, "w", encoding="utf-8") as devnull:
            os.dup2(devnull.fileno(), 2)
            yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)


def import_gffbase(repo_root: Path):
    try:
        import gffbase  # type: ignore

        return gffbase
    except ImportError as first_error:
        source_path = repo_root / ".sync" / "gffbase" / "python"
        if source_path.exists():
            sys.path.insert(0, str(source_path))
            try:
                import gffbase  # type: ignore

                return gffbase
            except ImportError:
                pass
        raise RuntimeError(
            "Could not import gffbase. Install it, e.g.\n"
            "  python3 -m pip install --target .tmp/gffbase_site gffbase==0.1.0\n"
            "  PYTHONPATH=.tmp/gffbase_site python3 scripts/gffbase_conformance_benchmark.py\n"
            "or install from the pinned upstream mirror under .sync/gffbase."
        ) from first_error


def run_gffbase_case(gffbase: Any, case: ConformanceCase, engine: str) -> dict[str, Any]:
    try:
        iterator = gffbase.parse_bytes(case.text, strict=True, engine=engine)
        records = list(iterator)
        first = records[0] if records else None
        attrs: dict[str, list[str]] = {}
        if first is not None:
            attrs = first.attributes_dict()
        return {
            "status": "accepted",
            "rows": len(records),
            "error_kind": "",
            "message": "",
            "attributes": attrs,
            "start": getattr(first, "start", None),
            "end": getattr(first, "end", None),
        }
    except Exception as exc:  # noqa: BLE001 - report exact parser error kind
        return {
            "status": "rejected",
            "rows": 0,
            "error_kind": getattr(exc, "kind", exc.__class__.__name__),
            "message": str(exc),
            "attributes": {},
            "start": None,
            "end": None,
        }


def run_duckhts_case(con: Any, path: Path, *, strict: bool = False) -> dict[str, Any]:
    strict_arg = ", strict := true" if strict else ""
    query = (
        "SELECT seqname, source, feature, start, \"end\", score, strand, frame, "
        "attributes, attributes_list, attributes_pairs "
        f"FROM read_gff('{sql_string(path)}', attributes_list := true, attributes_pairs := true{strict_arg})"
    )
    try:
        with suppress_c_stderr():
            rows = con.execute(query).fetchall()
        first = rows[0] if rows else None
        attrs = first[9] if first is not None else None
        if attrs is None:
            attrs = {}
        return {
            "status": "accepted",
            "rows": len(rows),
            "error_kind": "",
            "message": "",
            "first_row": first,
            "attributes": attrs,
            "start": first[3] if first is not None else None,
            "end": first[4] if first is not None else None,
            "score": first[5] if first is not None else None,
            "strand": first[6] if first is not None else None,
            "frame": first[7] if first is not None else None,
        }
    except Exception as exc:  # noqa: BLE001 - report exact DuckDB/extension error
        return {
            "status": "rejected",
            "rows": 0,
            "error_kind": exc.__class__.__name__,
            "message": str(exc),
            "first_row": None,
            "attributes": {},
            "start": None,
            "end": None,
            "score": None,
            "strand": None,
            "frame": None,
        }


def gffbase_expectation_ok(case: ConformanceCase, observed: dict[str, Any]) -> bool:
    if observed["status"] != case.expected_status:
        return False
    if case.expected_kind and observed["error_kind"] != case.expected_kind:
        return False
    return True


def classify_duckhts(case: ConformanceCase, gffbase_result: dict[str, Any], duckhts_result: dict[str, Any]) -> tuple[str, str]:
    if gffbase_result["status"] == "rejected":
        if duckhts_result["status"] == "rejected":
            return "match", f"rejected ({duckhts_result['error_kind']})"
        detail = (
            f"accepted {duckhts_result['rows']} row(s); "
            f"start={duckhts_result.get('start')!r}, end={duckhts_result.get('end')!r}, "
            f"score={duckhts_result.get('score')!r}, strand={duckhts_result.get('strand')!r}, "
            f"frame={duckhts_result.get('frame')!r}"
        )
        return "gap", detail

    if duckhts_result["status"] != "accepted":
        return "gap", f"rejected ({duckhts_result['error_kind']}): {duckhts_result['message']}"

    if duckhts_result["rows"] != gffbase_result["rows"]:
        return "gap", f"row count {duckhts_result['rows']} != GFFBase {gffbase_result['rows']}"

    attrs = duckhts_result.get("attributes") or {}
    if case.name == "all_dots_row":
        if duckhts_result.get("start") is not None or duckhts_result.get("end") is not None:
            return "gap", (
                "accepted, but unknown . coordinates became "
                f"start={duckhts_result.get('start')!r}, end={duckhts_result.get('end')!r}"
            )
    if case.name == "dbxref_multivalue":
        dbxref = attrs.get("Dbxref") if isinstance(attrs, dict) else None
        if dbxref != ["GeneID:1", "HGNC:HGNC:1"]:
            return "gap", f"accepted, but Dbxref map value is {dbxref!r} rather than two values"
    if case.name == "percent_decode_attribute":
        note = attrs.get("Note") if isinstance(attrs, dict) else None
        if note not in ("hello world", ["hello world"]):
            return "gap", f"accepted, but Note remains {note!r} rather than URL-decoded"

    return "match", f"accepted {duckhts_result['rows']} row(s)"


def run_gffbase_audit(gffbase: Any, out_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    native = bool(gffbase.native_available())
    rows.append({
        "check": "import",
        "status": "ok",
        "detail": getattr(gffbase, "__file__", ""),
        "value": getattr(gffbase, "__version__", ""),
    })
    rows.append({
        "check": "native_available",
        "status": "ok" if native else "fallback",
        "detail": "auto engine resolves to rust" if native else "auto engine resolves to python fallback",
        "value": str(native),
    })

    if native:
        mismatches = 0
        for case in CONFORMANCE_CASES:
            rust = run_gffbase_case(gffbase, case, "rust")
            py = run_gffbase_case(gffbase, case, "python")
            same = (rust["status"], rust["rows"], rust["error_kind"]) == (py["status"], py["rows"], py["error_kind"])
            if not same:
                mismatches += 1
        rows.append({
            "check": "rust_python_parser_parity",
            "status": "ok" if mismatches == 0 else "mismatch",
            "detail": f"{len(CONFORMANCE_CASES) - mismatches}/{len(CONFORMANCE_CASES)} strict cases matched by status,row_count,error_kind",
            "value": str(mismatches),
        })

    audit_dir = out_dir / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    gff3_path = audit_dir / "dialect.gff3"
    gtf_path = audit_dir / "dialect.gtf"
    gff3_path.write_bytes(b"##gff-version 3\nchr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g;Name=alpha\n")
    gtf_path.write_bytes(b"chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n")
    for label, path, expected_fmt in (("gff3", gff3_path, "gff3"), ("gtf", gtf_path, "gtf")):
        try:
            dialect = gffbase.detect_dialect(str(path), engine="rust" if native else "python")
            fmt = dialect.get("fmt", "")
            rows.append({
                "check": f"detect_dialect_{label}",
                "status": "ok" if fmt == expected_fmt else "mismatch",
                "detail": json.dumps(dialect, sort_keys=True),
                "value": fmt,
            })
        except Exception as exc:  # noqa: BLE001 - audit output
            rows.append({
                "check": f"detect_dialect_{label}",
                "status": "error",
                "detail": str(exc),
                "value": "",
            })
    return rows


def write_case_files(out_dir: Path) -> dict[str, Path]:
    cases_dir = out_dir / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for case in CONFORMANCE_CASES:
        path = cases_dir / f"{case.name}.gff3"
        path.write_bytes(case.text)
        paths[case.name] = path
    return paths


def run_conformance(gffbase: Any, con: Any, out_dir: Path, engine: str) -> list[dict[str, Any]]:
    paths = write_case_files(out_dir)
    rows: list[dict[str, Any]] = []
    for case in CONFORMANCE_CASES:
        gffbase_result = run_gffbase_case(gffbase, case, engine)
        duckhts_default = run_duckhts_case(con, paths[case.name], strict=False)
        default_verdict, default_detail = classify_duckhts(case, gffbase_result, duckhts_default)
        duckhts_strict = run_duckhts_case(con, paths[case.name], strict=True)
        strict_verdict, strict_detail = classify_duckhts(case, gffbase_result, duckhts_strict)
        rows.append(
            {
                "case": case.name,
                "focus": case.focus,
                "expected_gffbase_status": case.expected_status,
                "expected_gffbase_kind": case.expected_kind,
                "gffbase_status": gffbase_result["status"],
                "gffbase_rows": gffbase_result["rows"],
                "gffbase_error_kind": gffbase_result["error_kind"],
                "gffbase_expectation_ok": gffbase_expectation_ok(case, gffbase_result),
                "duckhts_default_status": duckhts_default["status"],
                "duckhts_default_rows": duckhts_default["rows"],
                "duckhts_default_error_kind": duckhts_default["error_kind"],
                "duckhts_default_vs_gffbase": default_verdict,
                "duckhts_default_detail": default_detail,
                "duckhts_strict_status": duckhts_strict["status"],
                "duckhts_strict_rows": duckhts_strict["rows"],
                "duckhts_strict_error_kind": duckhts_strict["error_kind"],
                "duckhts_strict_vs_gffbase": strict_verdict,
                "duckhts_strict_detail": strict_detail,
            }
        )
    return rows


def write_synthetic_gff3(path: Path, rows: int) -> int:
    chroms = ("chr1", "chr2", "chr3", "chr4")
    feature_cycle = ("gene", "mRNA", "exon", "CDS", "exon", "CDS")
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")
        handle.write("##sequence-region chr1 1 100000000\n")
        for i in range(rows):
            chrom = chroms[i % len(chroms)]
            feature = feature_cycle[i % len(feature_cycle)]
            start = i * 7 + 1
            end = start + 50 + (i % 200)
            strand = "+" if i % 2 == 0 else "-"
            phase = str(i % 3) if feature == "CDS" else "."
            ident = f"{feature}_{i}"
            parent = f"Parent=gene_{max(0, i - (i % 17))};" if feature not in {"gene"} else ""
            attrs = (
                f"ID={ident};{parent}Name={feature}-{i};"
                f"Dbxref=GeneID:{i},HGNC:HGNC:{i};Note=hello%20world"
            )
            handle.write(
                f"{chrom}\tduckhts_bench\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"
            )
    return rows


def timed(label: str, fn: Callable[[], Any], passes: int, warmup: int) -> tuple[Any, list[float]]:
    result: Any = None
    for _ in range(warmup):
        result = fn()
    times: list[float] = []
    for _ in range(passes):
        start = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - start)
    return result, times


def summarize_timing(tool: str, variant: str, rows: int, bytes_read: int, passes: int, times: list[float], extra: dict[str, Any] | None = None) -> dict[str, Any]:
    median_sec = statistics.median(times)
    row: dict[str, Any] = {
        "tool": tool,
        "variant": variant,
        "rows": rows,
        "bytes_read": bytes_read,
        "passes": passes,
        "min_sec": min(times),
        "median_sec": median_sec,
        "mean_sec": statistics.mean(times),
        "max_sec": max(times),
        "rows_per_sec": rows / median_sec if median_sec > 0 else None,
        "mb_per_sec": (bytes_read / (1024 * 1024)) / median_sec if median_sec > 0 else None,
    }
    if extra:
        row.update(extra)
    return row


def run_performance(gffbase: Any, con: Any, synthetic_path: Path, rows: int, passes: int, warmup: int, engine: str, include_create_db: bool) -> list[dict[str, Any]]:
    bytes_read = synthetic_path.stat().st_size
    synthetic_sql = sql_string(synthetic_path)
    results: list[dict[str, Any]] = []
    dataset_extra = {"dataset": synthetic_path.name, "path": str(synthetic_path)}

    def duckhts_count() -> int:
        with suppress_c_stderr():
            return con.execute(f"SELECT count(*) FROM read_gff('{synthetic_sql}')").fetchone()[0]

    value, times = timed("duckhts_count", duckhts_count, passes, warmup)
    observed_rows = int(value)
    results.append(summarize_timing("DuckHTS", "read_gff COUNT(*)", observed_rows, bytes_read, passes, times, dataset_extra))

    def duckhts_strict_count() -> int:
        with suppress_c_stderr():
            return con.execute(f"SELECT count(*) FROM read_gff('{synthetic_sql}', strict := true)").fetchone()[0]

    value, times = timed("duckhts_strict_count", duckhts_strict_count, passes, warmup)
    results.append(summarize_timing("DuckHTS", "read_gff strict COUNT(*)", int(value), bytes_read, passes, times, dataset_extra))

    def duckhts_projection() -> tuple[int, int, int]:
        with suppress_c_stderr():
            return con.execute(
                "SELECT count(*), sum(\"end\" - start + 1), "
                "sum(CASE WHEN strand = '+' THEN 1 ELSE 0 END) "
                f"FROM read_gff('{synthetic_sql}') WHERE feature IN ('exon', 'CDS')"
            ).fetchone()

    value, times = timed("duckhts_projection", duckhts_projection, passes, warmup)
    results.append(
        summarize_timing(
            "DuckHTS",
            "read_gff projected filter/sum",
            observed_rows,
            bytes_read,
            passes,
            times,
            {**dataset_extra, "selected_rows": int(value[0])},
        )
    )

    def duckhts_attributes() -> tuple[int, int]:
        with suppress_c_stderr():
            return con.execute(
                "SELECT count(*), count(map_extract_value(attributes_map, 'ID')) "
                f"FROM read_gff('{synthetic_sql}', attributes_map := true)"
            ).fetchone()

    value, times = timed("duckhts_attributes", duckhts_attributes, passes, warmup)
    results.append(
        summarize_timing(
            "DuckHTS",
            "read_gff attributes_map",
            int(value[0]),
            bytes_read,
            passes,
            times,
            {**dataset_extra, "selected_rows": int(value[1])},
        )
    )

    def duckhts_attributes_list() -> tuple[int, int]:
        with suppress_c_stderr():
            return con.execute(
                "SELECT count(*), count(map_extract_value(attributes_list, 'ID')) "
                f"FROM read_gff('{synthetic_sql}', attributes_list := true)"
            ).fetchone()

    value, times = timed("duckhts_attributes_list", duckhts_attributes_list, passes, warmup)
    results.append(
        summarize_timing(
            "DuckHTS",
            "read_gff attributes_list",
            int(value[0]),
            bytes_read,
            passes,
            times,
            {**dataset_extra, "selected_rows": int(value[1])},
        )
    )

    def duckhts_attributes_pairs() -> tuple[int, int]:
        with suppress_c_stderr():
            return con.execute(
                "SELECT count(*), sum(list_count(attributes_pairs)) "
                f"FROM read_gff('{synthetic_sql}', attributes_pairs := true)"
            ).fetchone()

    value, times = timed("duckhts_attributes_pairs", duckhts_attributes_pairs, passes, warmup)
    results.append(
        summarize_timing(
            "DuckHTS",
            "read_gff attributes_pairs",
            int(value[0]),
            bytes_read,
            passes,
            times,
            {**dataset_extra, "selected_rows": int(value[1] or 0)},
        )
    )

    def duckhts_attributes_all() -> tuple[int, int, int, int]:
        with suppress_c_stderr():
            return con.execute(
                "SELECT count(*), count(map_extract_value(attributes_map, 'ID')), "
                "count(map_extract_value(attributes_list, 'ID')), sum(list_count(attributes_pairs)) "
                f"FROM read_gff('{synthetic_sql}', attributes_map := true, attributes_list := true, attributes_pairs := true)"
            ).fetchone()

    value, times = timed("duckhts_attributes_all", duckhts_attributes_all, passes, warmup)
    results.append(
        summarize_timing(
            "DuckHTS",
            "read_gff attributes_map+list+pairs",
            int(value[0]),
            bytes_read,
            passes,
            times,
            {**dataset_extra, "selected_rows": int(value[3] or 0)},
        )
    )

    def gffbase_parse() -> int:
        n = 0
        for _feature in gffbase.parse_gff(str(synthetic_path), strict=True, engine=engine):
            n += 1
        return n

    value, times = timed("gffbase_parse", gffbase_parse, passes, warmup)
    results.append(
        summarize_timing(
            "GFFBase",
            f"parse_gff strict ({engine})",
            int(value),
            bytes_read,
            passes,
            times,
            dataset_extra,
        )
    )

    if include_create_db:
        db_dir = synthetic_path.parent / "gffbase_dbs"
        db_dir.mkdir(parents=True, exist_ok=True)
        counter = {"i": 0}

        def gffbase_create_db() -> int:
            counter["i"] += 1
            db_path = db_dir / f"synthetic_{counter['i']}.duckdb"
            if db_path.exists():
                db_path.unlink()
            db = gffbase.create_db(str(synthetic_path), str(db_path), force=True)
            try:
                count = db.conn.execute("SELECT count(*) FROM features").fetchone()[0]
            finally:
                try:
                    db.conn.close()
                except Exception:
                    pass
            return int(count)

        value, times = timed("gffbase_create_db", gffbase_create_db, passes, 0)
        results.append(
            summarize_timing(
                "GFFBase",
                "create_db ingest",
                int(value),
                bytes_read,
                passes,
                times,
                dataset_extra,
            )
        )

    return results


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str]) -> argparse.Namespace:
    root = repo_root_from_script()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extension",
        type=Path,
        default=root / "build" / "release" / "duckhts.duckdb_extension",
        help="Path to duckhts.duckdb_extension",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / ".tmp" / "gffbase_conformance",
        help="Directory for CSV/JSON outputs and generated fixtures",
    )
    parser.add_argument("--rows", type=int, default=200_000, help="Synthetic GFF3 feature rows for timing")
    parser.add_argument("--passes", type=int, default=3, help="Timed passes per benchmark")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup passes before timing")
    parser.add_argument("--duckdb-threads", type=int, default=4, help="PRAGMA threads for DuckDB")
    parser.add_argument(
        "--gffbase-engine",
        default="auto",
        choices=("auto", "rust", "python"),
        help="GFFBase parser engine",
    )
    parser.add_argument(
        "--human-gff",
        action="append",
        default=[],
        type=Path,
        help="Additional real human-scale GFF3/GFF file to benchmark. May be supplied multiple times.",
    )
    parser.add_argument(
        "--include-create-db",
        action="store_true",
        help="Also time gffbase.create_db ingestion. This is slower and compares a different workload shape.",
    )
    parser.add_argument(
        "--skip-performance",
        action="store_true",
        help="Only run the conformance matrix",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    repo_root = repo_root_from_script()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not args.extension.exists():
        raise SystemExit(f"DuckHTS extension not found: {args.extension}")

    gffbase = import_gffbase(repo_root)
    effective_gffbase_engine = args.gffbase_engine
    if effective_gffbase_engine == "auto":
        effective_gffbase_engine = "rust" if gffbase.native_available() else "python"

    import duckdb  # type: ignore

    con = duckdb.connect(":memory:", config={"allow_unsigned_extensions": "true"})
    con.execute(f"LOAD '{sql_string(args.extension.resolve())}'")
    con.execute(f"PRAGMA threads={int(args.duckdb_threads)}")

    audit_rows = run_gffbase_audit(gffbase, args.out_dir)
    write_csv(args.out_dir / "gffbase_audit.csv", audit_rows)

    conformance_rows = run_conformance(gffbase, con, args.out_dir, effective_gffbase_engine)
    write_csv(args.out_dir / "conformance.csv", conformance_rows)

    performance_rows: list[dict[str, Any]] = []
    synthetic_path = args.out_dir / f"synthetic_{args.rows}.gff3"
    if not args.skip_performance:
        write_synthetic_gff3(synthetic_path, args.rows)
        performance_rows = run_performance(
            gffbase,
            con,
            synthetic_path,
            args.rows,
            args.passes,
            args.warmup,
            effective_gffbase_engine,
            args.include_create_db,
        )
        for human_path in args.human_gff:
            if not human_path.exists():
                raise SystemExit(f"human-scale GFF path not found: {human_path}")
            performance_rows.extend(run_performance(
                gffbase,
                con,
                human_path,
                0,
                args.passes,
                args.warmup,
                effective_gffbase_engine,
                args.include_create_db,
            ))
        write_csv(args.out_dir / "performance.csv", performance_rows)

    gffbase_sync = repo_root / ".sync" / "gffbase"
    metadata = {
        "duckhts_git_rev": git_rev(repo_root),
        "duckhts_extension": str(args.extension.resolve()),
        "duckhts_extension_size": args.extension.stat().st_size,
        "gffbase_version": getattr(gffbase, "__version__", ""),
        "gffbase_native_available": bool(gffbase.native_available()),
        "gffbase_engine_requested": args.gffbase_engine,
        "gffbase_engine": effective_gffbase_engine,
        "gffbase_upstream_url": GFFBASE_UPSTREAM_URL,
        "gffbase_upstream_commit": GFFBASE_UPSTREAM_COMMIT,
        "gffbase_sync_rev": git_rev(gffbase_sync) if gffbase_sync.exists() else "",
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "server_specs": collect_system_specs(),
        "duckdb_python_version": getattr(duckdb, "__version__", ""),
        "rows": args.rows,
        "passes": args.passes,
        "warmup": args.warmup,
        "duckdb_threads": args.duckdb_threads,
        "synthetic_path": str(synthetic_path.resolve()) if synthetic_path.exists() else "",
        "synthetic_bytes": synthetic_path.stat().st_size if synthetic_path.exists() else 0,
        "human_gff": [str(p.resolve()) for p in args.human_gff],
        "case_source": ".sync/gffbase/tests/test_ncbi_compliance.py",
    }
    (args.out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    gaps = sum(1 for row in conformance_rows if row["duckhts_strict_vs_gffbase"] == "gap")
    print(f"wrote {args.out_dir / 'gffbase_audit.csv'} ({len(audit_rows)} audit rows)")
    print(f"wrote {args.out_dir / 'conformance.csv'} ({len(conformance_rows)} cases; {gaps} DuckHTS strict gaps)")
    if performance_rows:
        print(f"wrote {args.out_dir / 'performance.csv'} ({len(performance_rows)} timing rows)")
    print(f"wrote {args.out_dir / 'metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
