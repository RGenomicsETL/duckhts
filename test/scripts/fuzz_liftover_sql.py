#!/usr/bin/env python3
"""Replayable malformed-input fuzzing for bcftools_liftover()."""

import argparse
from pathlib import Path
import random
import tempfile

import duckdb

SQL = """SELECT bcftools_liftover(
    ?::VARCHAR, ?::BIGINT, ?::VARCHAR, ?::VARCHAR,
    ?::VARCHAR, ?::VARCHAR, ?::VARCHAR,
    ?::INTEGER, ?::INTEGER, false, NULL::BIGINT, true)"""
LIMIT_ERROR = "clip-pad alignment exceeds the 4194304-cell limit"


def execute_case(connection, values, seed, trial, expected_error=None):
    try:
        connection.execute(SQL, values).fetchone()
        if expected_error:
            raise AssertionError(f"seed={seed} trial={trial}: expected {expected_error!r}")
    except duckdb.Error as error:
        if expected_error and expected_error not in str(error):
            raise AssertionError(f"seed={seed} trial={trial}: {error}") from error
    if connection.execute("SELECT 42").fetchone() != (42,):
        raise AssertionError(f"seed={seed} trial={trial}: connection did not recover")


def random_text(rng, maximum):
    alphabet = "ACGTNacgtn,*<>[]:/_-' \t"
    return "".join(rng.choice(alphabet) for _ in range(rng.randrange(maximum + 1)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extension", required=True)
    parser.add_argument("--seed", type=int, default=169)
    parser.add_argument("--trials", type=int, default=250)
    parser.add_argument("--data-dir", default="test/data")
    args = parser.parse_args()
    data = Path(args.data_dir).resolve()
    extension = Path(args.extension).resolve()
    source = data / "liftover_nw_limit_src.fa"
    destination = data / "liftover_nw_limit_dst.fa"
    chain = data / "liftover_nw_limit.chain"
    reference = source.read_text().splitlines()[1][1:2201]
    rng = random.Random(args.seed)

    connection = duckdb.connect(config={"allow_unsigned_extensions": "true"})
    connection.execute(f"LOAD '{str(extension).replace(chr(39), chr(39) * 2)}'")
    base = ["chrS", 2, reference, "C", str(chain), str(destination), str(source), 1, 250]
    execute_case(connection, base, args.seed, "oversized", LIMIT_ERROR)

    fixed = [
        [None, 2, "C", "T", *base[4:]],
        ["", 2, "C", "T", *base[4:]],
        ["chrS", 0, "C", "T", *base[4:]],
        ["chrS", 2, "", "", *base[4:]],
    ]
    for trial, values in enumerate(fixed):
        execute_case(connection, values, args.seed, f"fixed-{trial}")

    malformed_templates = [
        "", "not a chain\n", "chain 1 chrS x + 0 1 chrD 1 + 0 1 1\n1\n",
        "chain 1 chrS 10 - 0 10 chrD 10 + 0 10 1\n10\n",
        "chain 1 chrS 10 + 0 10 chrD 10 + 0 10 1\n9\n",
    ]
    with tempfile.TemporaryDirectory(prefix="duckhts_liftover_fuzz_") as directory:
        malformed = Path(directory) / "input.chain"
        for trial in range(args.trials):
            if trial < len(malformed_templates):
                malformed.write_text(malformed_templates[trial])
            else:
                malformed.write_text(random_text(rng, 1024) + "\n")
            chain_path = str(malformed) if rng.randrange(3) else str(chain)
            chrom = rng.choice([None, "", "chrS", "chrD", random_text(rng, 32)])
            pos = rng.choice([None, -1, 0, 1, 2, 4005, 2**63 - 1])
            ref = rng.choice([None, "", "C", "*", "<DEL>", random_text(rng, 5000)])
            alt = rng.choice([None, "", "T", "*", "<INS>", random_text(rng, 5000)])
            values = [chrom, pos, ref, alt, chain_path, str(destination), str(source),
                      rng.choice([-1, 0, 1, 127]), rng.choice([-1, 0, 1, 250, 4096])]
            execute_case(connection, values, args.seed, trial)
    print(f"liftover SQL fuzz: OK ({args.trials} trials, seed={args.seed})")


if __name__ == "__main__":
    main()
