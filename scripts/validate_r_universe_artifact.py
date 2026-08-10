#!/usr/bin/env python3
"""Validate an R-universe Rduckhts binary package receipt.

The R package DESCRIPTION and the bundled DuckDB extension footer are two
independent records of the build target.  Artifact publishing must use their
intersection, not the coarse operating-system label returned by R-universe.
"""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FOOTER_SIZE = 534
FOOTER_MAGIC = b"\x00\x93\x04\x10"
FOOTER_SIGNATURE = b"duckdb_signature"
FOOTER_HEADER = b"\x80\x04"

# These names are the DuckDB extension platform values emitted by the package
# configure scripts.  linux_i686_musl is the existing receipt value for the
# Emscripten/webR package path.
SUPPORTED_PLATFORMS = frozenset(
    {
        "linux_amd64",
        "linux_arm64",
        "linux_i686_musl",
        "osx_amd64",
        "osx_arm64",
        "windows_amd64_mingw",
        "windows_arm64_mingw",
    }
)
API_OS_PLATFORMS = {
    "linux": frozenset({"linux_amd64", "linux_arm64"}),
    "mac": frozenset({"osx_amd64", "osx_arm64"}),
    "win": frozenset({"windows_amd64_mingw", "windows_arm64_mingw"}),
    "wasm": frozenset({"linux_i686_musl"}),
}

DESCRIPTION_PATH = "Rduckhts/DESCRIPTION"
EXTENSION_SUFFIX = "Rduckhts/duckhts_extension/build/duckhts.duckdb_extension"


class ArtifactValidationError(ValueError):
    """A package does not contain a self-consistent build receipt."""


@dataclass(frozen=True)
class ArtifactReceipt:
    package_version: str
    built_r_version: str
    built_target: str
    platform: str
    duckdb_version: str


def _find_unique_member(names: Iterable[str], expected: str) -> str:
    matches = [name for name in names if name == expected or name.endswith("/" + expected)]
    if len(matches) != 1:
        if not matches:
            raise ArtifactValidationError(f"package is missing {expected}")
        raise ArtifactValidationError(f"package contains multiple copies of {expected}")
    return matches[0]


def _read_tar_member(package: Path, expected: str) -> bytes:
    try:
        with tarfile.open(package, mode="r:*") as archive:
            member_name = _find_unique_member((m.name for m in archive.getmembers()), expected)
            member = archive.getmember(member_name)
            if not member.isfile():
                raise ArtifactValidationError(
                    f"package member is not a regular file: {member_name}"
                )
            stream = archive.extractfile(member)
            if stream is None:
                raise ArtifactValidationError(f"could not read package member: {member_name}")
            return stream.read()
    except ArtifactValidationError:
        raise
    except (OSError, tarfile.TarError) as exc:
        raise ArtifactValidationError(f"not a readable tar package: {package}") from exc


def _read_zip_member(package: Path, expected: str) -> bytes:
    try:
        with zipfile.ZipFile(package) as archive:
            member_name = _find_unique_member(archive.namelist(), expected)
            info = archive.getinfo(member_name)
            if info.is_dir():
                raise ArtifactValidationError(
                    f"package member is not a regular file: {member_name}"
                )
            return archive.read(info)
    except ArtifactValidationError:
        raise
    except (OSError, zipfile.BadZipFile, KeyError) as exc:
        raise ArtifactValidationError(f"not a readable zip package: {package}") from exc


def _read_member(package: Path, expected: str) -> bytes:
    try:
        return _read_tar_member(package, expected)
    except ArtifactValidationError as tar_error:
        try:
            return _read_zip_member(package, expected)
        except ArtifactValidationError as zip_error:
            raise ArtifactValidationError(
                f"could not read {expected} from {package}: {tar_error}; {zip_error}"
            ) from zip_error


def _dcf_fields(description: bytes) -> dict[str, str]:
    try:
        text = description.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ArtifactValidationError("DESCRIPTION is not UTF-8") from exc

    fields: dict[str, str] = {}
    current: str | None = None
    for line in text.splitlines():
        if line[:1].isspace():
            if current is not None:
                fields[current] += " " + line.strip()
            continue
        if ":" not in line:
            current = None
            continue
        name, value = line.split(":", 1)
        current = name.strip()
        fields[current] = value.strip()
    return fields


def _built_info(fields: dict[str, str]) -> tuple[str, str]:
    built = fields.get("Built", "")
    parts = [part.strip() for part in built.split(";")]
    if len(parts) < 2 or not parts[0].startswith("R ") or not parts[1]:
        raise ArtifactValidationError("DESCRIPTION has no valid Built target")
    built_r_version = parts[0][2:].strip()
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)+", built_r_version):
        raise ArtifactValidationError(
            f"DESCRIPTION has an invalid Built R version: {built_r_version}"
        )
    return built_r_version, parts[1]


def target_to_platform(target: str) -> str:
    """Map an R Built target to the package's exact DuckDB platform value."""
    normalized = target.strip().lower()
    if normalized == "wasm32-unknown-emscripten":
        return "linux_i686_musl"
    if normalized.startswith(("aarch64-", "arm64-")) and "-linux-" in normalized:
        return "linux_arm64"
    if normalized.startswith("x86_64-") and "-linux-" in normalized:
        return "linux_amd64"
    if normalized.startswith(("aarch64-", "arm64-")) and "-apple-darwin" in normalized:
        return "osx_arm64"
    if normalized.startswith("x86_64-") and "-apple-darwin" in normalized:
        return "osx_amd64"
    if normalized.startswith(("aarch64-", "arm64-")) and normalized.endswith("-w64-mingw32"):
        return "windows_arm64_mingw"
    if normalized.startswith("x86_64-") and normalized.endswith("-w64-mingw32"):
        return "windows_amd64_mingw"
    raise ArtifactValidationError(f"unsupported R Built target: {target}")


def _field(footer: bytes, offset: int, name: str) -> str:
    raw = footer[offset : offset + 32]
    if len(raw) != 32:
        raise ArtifactValidationError(f"DuckDB footer is missing {name}")
    value, separator, padding = raw.partition(b"\0")
    if not separator:
        padding = b""
    if any(byte != 0 for byte in padding):
        raise ArtifactValidationError(f"DuckDB footer has non-zero {name} padding")
    try:
        decoded = value.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ArtifactValidationError(f"DuckDB footer {name} is not ASCII") from exc
    if not decoded:
        raise ArtifactValidationError(f"DuckDB footer has an empty {name}")
    return decoded


def _footer_receipt(extension: bytes) -> tuple[str, str, str]:
    if len(extension) < FOOTER_SIZE:
        raise ArtifactValidationError("DuckDB extension is shorter than its metadata footer")
    footer = extension[-FOOTER_SIZE:]
    if footer[:4] != FOOTER_MAGIC or footer[4:20] != FOOTER_SIGNATURE:
        raise ArtifactValidationError("DuckDB extension has no recognized metadata footer")
    if footer[20:22] != FOOTER_HEADER:
        raise ArtifactValidationError("DuckDB extension metadata footer has an invalid header")
    extension_version = _field(footer, 150, "extension version")
    duckdb_version = _field(footer, 182, "DuckDB API version")
    platform = _field(footer, 214, "platform")
    if platform not in SUPPORTED_PLATFORMS:
        raise ArtifactValidationError(f"unsupported DuckDB footer platform: {platform}")
    return extension_version, duckdb_version, platform


def validate_artifact(
    package: str | Path,
    api_os: str,
    expected_version: str,
    expected_r_version: str | None = None,
) -> ArtifactReceipt:
    """Validate one R-universe package and return its exact receipt."""
    package_path = Path(package)
    if not package_path.is_file():
        raise ArtifactValidationError(f"package does not exist: {package_path}")
    if api_os not in API_OS_PLATFORMS:
        raise ArtifactValidationError(f"unsupported R-universe operating system: {api_os}")
    if not expected_version:
        raise ArtifactValidationError("expected package version is empty")

    fields = _dcf_fields(_read_member(package_path, DESCRIPTION_PATH))
    package_version = fields.get("Version", "")
    if package_version != expected_version:
        raise ArtifactValidationError(
            f"package Version {package_version!r} disagrees with API Version {expected_version!r}"
        )
    if fields.get("NeedsCompilation") != "yes":
        raise ArtifactValidationError("binary package does not declare NeedsCompilation: yes")

    built_r_version, built_target = _built_info(fields)
    if expected_r_version is not None and built_r_version != expected_r_version:
        raise ArtifactValidationError(
            f"package Built R version {built_r_version!r} disagrees with API R version "
            f"{expected_r_version!r}"
        )
    built_platform = target_to_platform(built_target)
    extension_version, duckdb_version, footer_platform = _footer_receipt(
        _read_member(package_path, EXTENSION_SUFFIX)
    )
    if extension_version != expected_version:
        raise ArtifactValidationError(
            f"DuckDB footer extension version {extension_version!r} disagrees with API Version "
            f"{expected_version!r}"
        )
    if built_platform != footer_platform:
        raise ArtifactValidationError(
            f"DESCRIPTION Built target {built_target!r} resolves to {built_platform!r}, "
            f"but DuckDB footer records {footer_platform!r}"
        )
    if footer_platform not in API_OS_PLATFORMS[api_os]:
        raise ArtifactValidationError(
            f"R-universe os {api_os!r} disagrees with DuckDB footer platform {footer_platform!r}"
        )
    return ArtifactReceipt(
        expected_version,
        built_r_version,
        built_target,
        footer_platform,
        duckdb_version,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--api-os", required=True, choices=sorted(API_OS_PLATFORMS))
    parser.add_argument("--expected-version", required=True)
    parser.add_argument("--expected-r-version", required=True)
    args = parser.parse_args(argv)
    try:
        receipt = validate_artifact(
            args.package,
            args.api_os,
            args.expected_version,
            args.expected_r_version,
        )
    except ArtifactValidationError as exc:
        print(f"artifact receipt validation failed: {exc}", file=sys.stderr)
        return 1
    print(receipt.platform)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
