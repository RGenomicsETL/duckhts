#!/usr/bin/env python3
"""Network-free regression tests for failure-safe R-universe publication."""

from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from publish_r_universe_artifacts import POINTER_NAME, publish, release_assets  # noqa: E402


VERSION = "2.0.2.9000.0.1.1"
PREFIX = f"Rduckhts_{VERSION}_"
FIRST = PREFIX + "linux_amd64_r4.6.1.zip"
SECOND = PREFIX + "linux_arm64_r4.6.1.zip"
OLDER = "Rduckhts_2.0.2.0.1.1_linux_amd64_r4.6.1.zip"
OLD_POINTER = b"2.0.2.0.1.1\n"


def bundle(payload: bytes, year: int = 2026) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        member = zipfile.ZipInfo("Rduckhts.tar.gz", date_time=(year, 1, 1, 0, 0, 0))
        archive.writestr(member, payload)
    return output.getvalue()


class FakeGitHub:
    def __init__(self) -> None:
        self.files = {OLDER: bundle(b"older package"), POINTER_NAME: OLD_POINTER}
        self.states = {}
        self.calls = []
        self.fail_upload = None
        self.corrupt_download = None
        self.hide_asset = None

    def run(self, *args: str) -> str:
        self.calls.append(args)
        if args[0] == "api":
            if "/tags/" in args[1]:
                return json.dumps({"id": 17})
            if args[-2:] != ("--paginate", "--slurp"):
                raise AssertionError("Asset listing must include all pages")
            assets = [
                {"name": name, "state": self.states.get(name, "uploaded")}
                for name in self.files if name != self.hide_asset
            ]
            # Multiple pages are intentional, including a final empty page.
            return json.dumps([assets[:1], assets[1:], []])
        if args[:2] == ("release", "upload"):
            path = Path(args[3])
            if path.name in self.files:
                if "--clobber" not in args:
                    raise RuntimeError("Asset already exists")
                del self.files[path.name]
            if self.fail_upload == path.name:
                raise RuntimeError("Injected upload failure")
            self.files[path.name] = path.read_bytes()
            return ""
        if args[:2] == ("release", "download"):
            name = args[args.index("--pattern") + 1]
            directory = Path(args[args.index("--dir") + 1])
            content = self.files[name]
            if self.corrupt_download == name:
                content = bundle(b"corrupt package")
            (directory / name).write_bytes(content)
            return ""
        raise AssertionError(f"Unexpected GitHub command: {args}")

    def uploads(self) -> list[str]:
        return [Path(args[3]).name for args in self.calls if args[:2] == ("release", "upload")]


class PublicationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        (self.directory / "zips").mkdir()
        (self.directory / POINTER_NAME).write_text(VERSION + "\n", encoding="ascii")
        (self.directory / "zips" / FIRST).write_bytes(bundle(b"amd64 package"))
        (self.directory / "zips" / SECOND).write_bytes(bundle(b"arm64 package"))
        self.github = FakeGitHub()
        self.original_files = dict(self.github.files)
        self.mock = patch("publish_r_universe_artifacts.run_gh", self.github.run)
        self.mock.start()
        self.addCleanup(self.mock.stop)

    def publish(self) -> None:
        publish(self.directory, "owner/repo")

    def test_uploads_and_verifies_all_binaries_before_pointer(self) -> None:
        self.publish()
        self.assertEqual(self.github.uploads(), [FIRST, SECOND, POINTER_NAME])
        self.assertEqual(self.github.files[OLDER], self.original_files[OLDER])
        self.assertEqual(self.github.files[POINTER_NAME], (VERSION + "\n").encode())
        pointer_call = next(i for i, args in enumerate(self.github.calls)
                            if args[:2] == ("release", "upload") and Path(args[3]).name == POINTER_NAME)
        prior_downloads = [args[args.index("--pattern") + 1] for args in self.github.calls[:pointer_call]
                           if args[:2] == ("release", "download")]
        self.assertCountEqual(prior_downloads, [FIRST, SECOND])
        for args in self.github.calls:
            if args[:2] == ("release", "upload") and Path(args[3]).name != POINTER_NAME:
                self.assertNotIn("--clobber", args)

    def test_initial_upload_failure_preserves_every_existing_asset(self) -> None:
        self.github.fail_upload = FIRST
        with self.assertRaisesRegex(RuntimeError, "upload failure"):
            self.publish()
        self.assertEqual(self.github.files, self.original_files)
        self.assertNotIn(POINTER_NAME, self.github.uploads())

    def test_partial_upload_is_retryable_without_deleting_successful_binary(self) -> None:
        self.github.fail_upload = SECOND
        with self.assertRaises(RuntimeError):
            self.publish()
        first_bytes = self.github.files[FIRST]
        self.assertEqual(self.github.files[POINTER_NAME], OLD_POINTER)
        self.github.fail_upload = None
        self.github.calls.clear()
        self.publish()
        self.assertEqual(self.github.uploads(), [SECOND, POINTER_NAME])
        self.assertEqual(self.github.files[FIRST], first_bytes)

    def test_same_package_with_different_zip_timestamp_is_reused(self) -> None:
        self.github.files[FIRST] = bundle(b"amd64 package", year=2025)
        original = self.github.files[FIRST]
        self.publish()
        self.assertEqual(self.github.uploads(), [SECOND, POINTER_NAME])
        self.assertEqual(self.github.files[FIRST], original)

    def test_conflicting_existing_binary_fails_before_any_upload(self) -> None:
        self.github.files[SECOND] = bundle(b"wrong package")
        before = dict(self.github.files)
        with self.assertRaisesRegex(ValueError, "Existing binary differs"):
            self.publish()
        self.assertEqual(self.github.files, before)
        self.assertEqual(self.github.uploads(), [])

    def test_unvalidated_same_version_asset_is_not_silently_retained_or_deleted(self) -> None:
        self.github.files[PREFIX + "obsolete_platform.zip"] = bundle(b"unsafe mapping")
        before = dict(self.github.files)
        with self.assertRaisesRegex(ValueError, "Unvalidated assets"):
            self.publish()
        self.assertEqual(self.github.files, before)
        self.assertEqual(self.github.uploads(), [])

    def test_incomplete_existing_upload_is_rejected_without_mutation(self) -> None:
        self.github.files[FIRST] = b""
        self.github.states[FIRST] = "starter"
        with self.assertRaisesRegex(ValueError, "not fully uploaded"):
            self.publish()
        self.assertEqual(self.github.uploads(), [])

    def test_corrupted_download_does_not_advance_pointer(self) -> None:
        self.github.corrupt_download = SECOND
        with self.assertRaises(ValueError):
            self.publish()
        self.assertEqual(self.github.files[POINTER_NAME], OLD_POINTER)
        self.assertNotIn(POINTER_NAME, self.github.uploads())

    def test_missing_published_asset_does_not_advance_pointer(self) -> None:
        self.github.hide_asset = SECOND
        with self.assertRaisesRegex(ValueError, "Published binary is missing"):
            self.publish()
        self.assertEqual(self.github.files[POINTER_NAME], OLD_POINTER)

    def test_pointer_failure_leaves_binaries_intact_and_retry_repairs_pointer(self) -> None:
        self.github.fail_upload = POINTER_NAME
        with self.assertRaises(RuntimeError):
            self.publish()
        self.assertNotIn(POINTER_NAME, self.github.files)
        binary_bytes = dict(self.github.files)
        self.github.fail_upload = None
        self.github.calls.clear()
        self.publish()
        self.assertEqual(self.github.uploads(), [POINTER_NAME])
        for name, content in binary_bytes.items():
            self.assertEqual(self.github.files[name], content)

    def test_empty_set_is_read_only(self) -> None:
        for path in (self.directory / "zips").iterdir():
            path.unlink()
        with self.assertRaisesRegex(ValueError, "No binary ZIPs"):
            self.publish()
        self.assertEqual(self.github.calls, [])

    def test_wrong_version_filename_is_read_only(self) -> None:
        (self.directory / "zips" / OLDER).write_bytes(bundle(b"wrong version"))
        with self.assertRaisesRegex(ValueError, "Unexpected binary name"):
            self.publish()
        self.assertEqual(self.github.calls, [])

    def test_asset_listing_uses_all_pages(self) -> None:
        self.assertEqual(set(release_assets("owner/repo", 17)), set(self.github.files))


if __name__ == "__main__":
    unittest.main()
