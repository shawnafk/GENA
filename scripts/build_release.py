"""Build a full release ZIP from a clean Git commit and pinned local resources."""

import argparse
import ast
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PACKAGES = {"common", "gena", "libs", "load", "plots"}
RELEASE_FILES = {
    "pyproject.toml", "readme.md", "readme.zh.md", "RELEASE.md", "NOTICE",
    "install.ps1", "release-resources.json", "scripts/verify_release.py",
}


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def build(ref, output):
    if git("status", "--porcelain").strip():
        raise ValueError("Commit or stash working-tree changes before building a release")
    commit = git("rev-parse", "--verify", f"{ref}^{{commit}}").decode().strip()
    archive = git("archive", "--format=tar", commit)
    files = {}
    with tarfile.open(fileobj=io.BytesIO(archive)) as source:
        for member in source.getmembers():
            path = PurePosixPath(member.name)
            selected = member.name in RELEASE_FILES or (
                path.parts[0] in PACKAGES and path.suffix == ".py"
            )
            if member.isfile() and selected:
                files[member.name] = source.extractfile(member).read()
    missing = RELEASE_FILES - files.keys()
    if missing:
        raise ValueError(f"Commit is missing release files: {sorted(missing)}")
    tree = ast.parse(files["gena/_version.py"])
    version = next(
        ast.literal_eval(node.value) for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets)
    )
    if ref.startswith("v") and ref != f"v{version}":
        raise ValueError(f"Tag {ref} does not match source version {version}")
    resources = json.loads(files["release-resources.json"])
    for name, expected in resources.items():
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or path.parts[0] != "libs":
            raise ValueError(f"Invalid resource path: {name}")
        # Metadata is versioned with the code. Large binaries live outside Git.
        data = git("show", f"{commit}:{name}") if name.endswith(".json") else (ROOT / name).read_bytes()
        if len(data) != expected["size"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
            raise ValueError(f"Resource differs from the committed fingerprint: {name}")
        files[name] = data
    manifest = {
        "version": version,
        "commit": commit,
        "variant": "full",
        "resources": resources,
        "files": {name: hashlib.sha256(data).hexdigest() for name, data in sorted(files.items())},
    }
    files["release-manifest.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    files["VERSION"] = (version + "\n").encode()
    # Include generated version metadata in the file checksum list as well.
    sums = "".join(f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in sorted(files.items()))
    files["SHA256SUMS.txt"] = sums.encode()
    output.mkdir(parents=True, exist_ok=True)
    destination = output / f"gena-code-{version}-full.zip"
    checksum_path = destination.with_suffix(".zip.sha256")
    if destination.exists() or checksum_path.exists():
        raise FileExistsError(f"Release output already exists: {destination}")
    temporary = destination.with_suffix(".zip.tmp")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as release:
            for name, data in sorted(files.items()):
                info = zipfile.ZipInfo(f"gena-code-{version}/{name}", date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                release.writestr(info, data, compresslevel=9)
        with temporary.open("rb") as stream:
            checksum = hashlib.file_digest(stream, "sha256").hexdigest()
        temporary.replace(destination)
        checksum_path.write_text(f"{checksum}  {destination.name}\n", encoding="ascii")
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Built {destination}\nCommit {commit}\nSHA256 {checksum}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", default="HEAD", help="Committed revision or matching version tag")
    parser.add_argument("--output", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    build(args.ref, args.output.resolve())
