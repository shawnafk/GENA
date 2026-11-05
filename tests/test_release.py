"""Release archives preserve provenance, omit engineering files and detect damage."""

import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from scripts import build_release, verify_release


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / "libs").mkdir()
        resource = b"small synthetic kernel"
        (self.root / "libs/de432s.bsp").write_bytes(resource)
        resources = {"libs/de432s.bsp": {
            "size": len(resource), "sha256": hashlib.sha256(resource).hexdigest(),
        }}
        files = {name: b"" for name in build_release.RELEASE_FILES}
        files.update({
            "release-resources.json": json.dumps(resources).encode(),
            "gena/_version.py": b'__version__ = "0.2.0rc1"\n',
            "load/loaders/gena.py": b'"""Namespace package loader."""\n',
            "gena/raw/doc/internal.pdf": b"excluded engineering document",
            "analysis/private.py": b"excluded analysis",
        })
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w") as archive:
            for name, data in files.items():
                entry = tarfile.TarInfo(name)
                entry.size = len(data)
                archive.addfile(entry, io.BytesIO(data))
        self.archive = buffer.getvalue()
        self.commit = "a" * 40

    def git(self, *args):
        return {
            "status": b"", "rev-parse": self.commit.encode(), "archive": self.archive,
        }[args[0]]

    def build(self, destination):
        with patch.object(build_release, "ROOT", self.root), patch.object(
            build_release, "git", side_effect=self.git,
        ):
            build_release.build("HEAD", destination)
        return destination / "gena-code-0.2.0rc1-full.zip"

    def test_archive_is_repeatable_and_contains_only_selected_files(self):
        first = self.build(self.root / "first")
        second = self.build(self.root / "second")
        self.assertEqual(first.read_bytes(), second.read_bytes())
        with ZipFile(first) as archive:
            names = archive.namelist()
            self.assertIn("gena-code-0.2.0rc1/load/loaders/gena.py", names)
            self.assertFalse(any("internal.pdf" in name or "private.py" in name for name in names))
            archive.extractall(self.root / "extracted")
        extracted = self.root / "extracted/gena-code-0.2.0rc1"
        verify_release.verify(extracted)
        manifest = json.loads((extracted / "release-manifest.json").read_text())
        self.assertEqual(manifest["commit"], self.commit)
        (extracted / "libs/de432s.bsp").write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            verify_release.verify(extracted)

    def test_changed_resource_fails_before_creating_an_archive(self):
        (self.root / "libs/de432s.bsp").write_bytes(b"unexpected replacement")
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            self.build(self.root / "output")
        self.assertFalse((self.root / "output").exists())

    def test_existing_release_is_preserved(self):
        destination = self.root / "output"
        archive = self.build(destination)
        original = archive.read_bytes()
        with self.assertRaises(FileExistsError):
            self.build(destination)
        self.assertEqual(archive.read_bytes(), original)
