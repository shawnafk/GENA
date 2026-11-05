"""Check the extracted release files before installation, using only the stdlib."""

import hashlib
import json
from pathlib import Path, PurePosixPath


def verify(root):
    root = root.resolve()
    checksums = {}
    for line in (root / "SHA256SUMS.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        path = PurePosixPath(name)
        target = (root / name).resolve()
        if path.is_absolute() or ".." in path.parts or not target.is_relative_to(root):
            raise ValueError(f"Invalid checksum path: {name}")
        with target.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != expected:
            raise ValueError(f"SHA256 mismatch: {name}")
        checksums[name] = expected
    manifest = json.loads((root / "release-manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["files"].items():
        if checksums.get(name) != expected:
            raise ValueError(f"Incomplete checksum list: {name}")
    if (root / "VERSION").read_text().strip() != manifest["version"]:
        raise ValueError("VERSION differs from release-manifest.json")
    print(f"Verified GENA-CODE {manifest['version']} ({manifest['commit']})")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parents[1])
