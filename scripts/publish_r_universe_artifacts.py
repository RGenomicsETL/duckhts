#!/usr/bin/env python3
"""Publish validated R-universe bundles without deleting existing binaries.

The caller validates package receipts before invoking this script. Existing ZIPs
are immutable: a repeated upload must contain the same package bytes. GitHub has
no atomic asset replacement, so only the small version pointer may be clobbered,
and only after every ZIP has been downloaded from the release and verified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


RELEASE_TAG = "r-universe-binaries"
POINTER_NAME = "latest-version.txt"


def run_gh(*args: str) -> str:
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"GitHub command failed: {result.stderr.strip()}")
    return result.stdout


def release_assets(repository: str, release_id: int) -> dict[str, dict]:
    pages = json.loads(run_gh(
        "api", f"repos/{repository}/releases/{release_id}/assets?per_page=100",
        "--paginate", "--slurp",
    ))
    assets = {}
    for page in pages:
        for asset in page:
            name = asset["name"]
            if name in assets:
                raise ValueError(f"Duplicate release asset: {name}")
            assets[name] = asset
    return assets


def package_digest(path: Path, expected: tuple | None = None) -> tuple[str, int, str]:
    """Compare package bytes, not the outer ZIP's download-time metadata."""
    with zipfile.ZipFile(path) as archive:
        members = archive.infolist()
        if len(members) != 1 or members[0].is_dir():
            raise ValueError(f"Expected one package in {path.name}")
        member = members[0]
        if Path(member.filename).name != member.filename or member.file_size == 0:
            raise ValueError(f"Invalid package member in {path.name}")
        if expected is not None and (member.filename, member.file_size) != expected[:2]:
            raise ValueError(f"Existing binary differs: {path.name}; use a new package version")
        digest = hashlib.sha256()
        with archive.open(member) as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        receipt = (member.filename, member.file_size, digest.hexdigest())
        if expected is not None and receipt != expected:
            raise ValueError(f"Existing binary differs: {path.name}; use a new package version")
        return receipt


def verify_binary(repository: str, asset: dict, expected: tuple) -> None:
    if asset.get("state") != "uploaded":
        raise ValueError(f"Release asset is not fully uploaded: {asset['name']}")
    with tempfile.TemporaryDirectory() as tmp:
        run_gh(
            "release", "download", RELEASE_TAG, "--repo", repository,
            "--pattern", asset["name"], "--dir", tmp,
        )
        package_digest(Path(tmp) / asset["name"], expected)


def check_asset_names(assets: dict, expected: dict, prefix: str) -> None:
    unexpected = sorted(name for name in assets if name.startswith(prefix) and name not in expected)
    if unexpected:
        raise ValueError(
            "Unvalidated assets already exist for this package version: "
            + ", ".join(unexpected)
            + "; inspect them manually or use a new package version"
        )


def publish(directory: Path, repository: str) -> None:
    pointer = directory / POINTER_NAME
    version = pointer.read_text(encoding="ascii").strip()
    if not re.fullmatch(r"[0-9]+(?:[.-][0-9]+)+", version):
        raise ValueError("Invalid package version")
    prefix = f"Rduckhts_{version}_"
    binaries = sorted((directory / "zips").glob("*.zip"))
    if not binaries:
        raise ValueError("No binary ZIPs to publish")
    expected = {}
    for path in binaries:
        if not path.name.startswith(prefix) or not re.fullmatch(r"[A-Za-z0-9_.-]+", path.name):
            raise ValueError(f"Unexpected binary name: {path.name}")
        expected[path.name] = package_digest(path)

    # A missing release or any preflight error is read-only. Do not silently
    # create a release or remove historical assets to make publication succeed.
    release = json.loads(run_gh("api", f"repos/{repository}/releases/tags/{RELEASE_TAG}"))
    assets = release_assets(repository, release["id"])
    check_asset_names(assets, expected, prefix)
    for name, receipt in expected.items():
        if name in assets:
            verify_binary(repository, assets[name], receipt)

    for path in binaries:
        if path.name not in assets:
            # No --clobber: even a concurrent upload must not delete an existing ZIP.
            run_gh("release", "upload", RELEASE_TAG, str(path), "--repo", repository)

    assets = release_assets(repository, release["id"])
    check_asset_names(assets, expected, prefix)
    for name, receipt in expected.items():
        if name not in assets:
            raise ValueError(f"Published binary is missing: {name}")
        verify_binary(repository, assets[name], receipt)

    # A failed pointer replacement can leave the pointer absent, but cannot
    # advertise a partially uploaded set. Retrying reuses verified binaries.
    run_gh("release", "upload", RELEASE_TAG, str(pointer), "--repo", repository, "--clobber")
    with tempfile.TemporaryDirectory() as tmp:
        run_gh(
            "release", "download", RELEASE_TAG, "--repo", repository,
            "--pattern", POINTER_NAME, "--dir", tmp,
        )
        if (Path(tmp) / POINTER_NAME).read_bytes() != pointer.read_bytes():
            raise ValueError("Published version pointer does not match the verified binaries")
    print(f"Verified {len(binaries)} binaries and published version {version}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("dist"))
    parser.add_argument("--repository", required=True)
    args = parser.parse_args()
    try:
        publish(args.directory, args.repository)
    except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as exc:
        print(f"R-universe publication failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
