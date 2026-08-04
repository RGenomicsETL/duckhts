#!/usr/bin/env python3
"""Mocked package receipt tests for the R-universe synchronizer."""

from __future__ import annotations

import io
import sys
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from validate_r_universe_artifact import (  # noqa: E402
    ArtifactValidationError,
    FOOTER_HEADER,
    FOOTER_MAGIC,
    FOOTER_SIGNATURE,
    FOOTER_SIZE,
    validate_artifact,
)


VERSION = "1.5.1.9000-0.1.5"
DESCRIPTION_PATH = "Rduckhts/DESCRIPTION"
EXTENSION_PATH = "Rduckhts/duckhts_extension/build/duckhts.duckdb_extension"


def _field(value: str) -> bytes:
    encoded = value.encode("ascii")
    if len(encoded) > 32:
        raise ValueError(value)
    return encoded + b"\0" * (32 - len(encoded))


def _extension(platform: str, version: str = VERSION) -> bytes:
    footer = b"".join(
        (
            FOOTER_MAGIC,
            FOOTER_SIGNATURE,
            FOOTER_HEADER,
            _field(""),
            _field(""),
            _field(""),
            _field("C_STRUCT"),
            _field(version),
            _field("v1.2.0"),
            _field(platform),
            _field("4"),
            b"\0" * 256,
        )
    )
    if len(footer) != FOOTER_SIZE:
        raise AssertionError(len(footer))
    return b"mock extension payload\0" + footer


def _description(
    target: str,
    version: str = VERSION,
    needs_compilation: str = "yes",
    r_version: str = "4.6.1",
) -> bytes:
    return (
        f"Package: Rduckhts\nVersion: {version}\nNeedsCompilation: {needs_compilation}\n"
        f"Built: R {r_version}; {target}; 2026-08-04 22:00:00 UTC; unix\n"
    ).encode("utf-8")


def _write_package(
    directory: Path,
    *,
    archive_kind: str,
    target: str,
    platform: str,
    version: str = VERSION,
    r_version: str = "4.6.1",
    include_extension: bool = True,
) -> Path:
    path = directory / ("package.tar.gz" if archive_kind == "tar" else "package.zip")
    members = {DESCRIPTION_PATH: _description(target, version, r_version=r_version)}
    if include_extension:
        members[EXTENSION_PATH] = _extension(platform, version)
    if archive_kind == "tar":
        with tarfile.open(path, "w:gz") as archive:
            for name, content in members.items():
                info = tarfile.TarInfo(name)
                info.size = len(content)
                archive.addfile(info, io.BytesIO(content))
    else:
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, content in members.items():
                archive.writestr(name, content)
    return path


class RUniverseArtifactReceiptTests(unittest.TestCase):
    def test_accepts_linux_arm64_tar_and_returns_exact_footer_platform(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="aarch64-unknown-linux-gnu",
                platform="linux_arm64",
            )
            receipt = validate_artifact(package, "linux", VERSION, "4.6.1")
        self.assertEqual(receipt.platform, "linux_arm64")

    def test_accepts_windows_amd64_zip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="zip",
                target="x86_64-w64-mingw32",
                platform="windows_amd64_mingw",
            )
            receipt = validate_artifact(package, "win", VERSION, "4.6.1")
        self.assertEqual(receipt.platform, "windows_amd64_mingw")

    def test_accepts_wasm_receipt_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="wasm32-unknown-emscripten",
                platform="linux_i686_musl",
                r_version="4.6.0",
            )
            receipt = validate_artifact(package, "wasm", VERSION, "4.6.0")
        self.assertEqual(receipt.platform, "linux_i686_musl")

    def test_rejects_built_and_footer_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="x86_64-pc-linux-gnu",
                platform="linux_arm64",
            )
            with self.assertRaisesRegex(ArtifactValidationError, "resolves to"):
                validate_artifact(package, "linux", VERSION, "4.6.1")

    def test_rejects_coarse_api_os_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="zip",
                target="aarch64-unknown-linux-gnu",
                platform="linux_arm64",
            )
            with self.assertRaisesRegex(ArtifactValidationError, "R-universe os"):
                validate_artifact(package, "win", VERSION, "4.6.1")

    def test_rejects_missing_extension_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="x86_64-pc-linux-gnu",
                platform="linux_amd64",
                include_extension=False,
            )
            with self.assertRaisesRegex(ArtifactValidationError, "missing"):
                validate_artifact(package, "linux", VERSION, "4.6.1")

    def test_rejects_package_version_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="x86_64-pc-linux-gnu",
                platform="linux_amd64",
                version="1.5.1.9000-0.1.4",
            )
            with self.assertRaisesRegex(ArtifactValidationError, "package Version"):
                validate_artifact(package, "linux", VERSION, "4.6.1")

    def test_rejects_r_version_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = _write_package(
                Path(tmp),
                archive_kind="tar",
                target="x86_64-pc-linux-gnu",
                platform="linux_amd64",
                r_version="4.5.3",
            )
            with self.assertRaisesRegex(ArtifactValidationError, "Built R version"):
                validate_artifact(package, "linux", VERSION, "4.6.1")


if __name__ == "__main__":
    unittest.main()
