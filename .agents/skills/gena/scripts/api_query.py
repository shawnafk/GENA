#!/usr/bin/env python3
"""Query Python symbols in GENA-CODE-Dev without importing project modules."""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from pathlib import Path
import re
import sys


PACKAGE_ROOTS = ("gena", "common", "load", "plots", "libs")
SKIP_PARTS = {"__pycache__", "deprecated"}
LEGACY_FILES = {"filter_od.py", "jinfa_ref.py"}


@dataclass
class Symbol:
    module: str
    qualname: str
    signature: str
    kind: str
    path: Path
    line: int
    doc: str


def find_repo_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "gena").is_dir() and (candidate / "common").is_dir():
            return candidate
    raise SystemExit("Could not locate a GENA-CODE-Dev root containing gena/ and common/.")


def module_name(root: Path, path: Path) -> str:
    relative = path.relative_to(root).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def format_args(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    try:
        return ast.unparse(node.args)
    except Exception:
        return "..."


def first_doc_line(node: ast.AST) -> str:
    doc = ast.get_docstring(node, clean=True) or ""
    return doc.splitlines()[0] if doc else ""


def collect_symbols(
    root: Path,
    include_private: bool = False,
    include_legacy: bool = False,
) -> list[Symbol]:
    symbols: list[Symbol] = []

    def visit(path: Path, module: str, body: list[ast.stmt], parent: str = "") -> None:
        for node in body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            if not include_private and node.name.startswith("_"):
                continue
            qualname = f"{parent}.{node.name}" if parent else node.name
            if isinstance(node, ast.ClassDef):
                signature = ""
                kind = "class"
            else:
                signature = format_args(node)
                kind = "async function" if isinstance(node, ast.AsyncFunctionDef) else "function"
            symbols.append(
                Symbol(module, qualname, signature, kind, path, node.lineno, first_doc_line(node))
            )
            if isinstance(node, ast.ClassDef):
                visit(path, module, node.body, qualname)

    for package in PACKAGE_ROOTS:
        package_dir = root / package
        if not package_dir.is_dir():
            continue
        for path in sorted(package_dir.rglob("*.py")):
            if any(part in SKIP_PARTS for part in path.parts):
                continue
            if not include_legacy and path.name in LEGACY_FILES:
                continue
            try:
                tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            except (OSError, SyntaxError, UnicodeDecodeError) as exc:
                print(f"warning: skipped {path}: {exc}", file=sys.stderr)
                continue
            visit(path, module_name(root, path), tree.body)
    return symbols


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find GENA-CODE-Dev functions/classes by name, module, signature, or docstring."
    )
    parser.add_argument("query", nargs="?", default="", help="case-insensitive text or regex")
    parser.add_argument("--root", type=Path, help="repository root; auto-detected by default")
    parser.add_argument("--module", help="restrict to a module prefix, e.g. gena.post")
    parser.add_argument("--regex", action="store_true", help="interpret query as a regular expression")
    parser.add_argument("--private", action="store_true", help="include underscore-prefixed symbols")
    parser.add_argument("--legacy", action="store_true", help="include historical raw modules")
    parser.add_argument("--list-modules", action="store_true", help="list scanned modules and symbol counts")
    args = parser.parse_args()

    root = (args.root.expanduser().resolve() if args.root else find_repo_root(Path.cwd().resolve()))
    symbols = collect_symbols(
        root,
        include_private=args.private,
        include_legacy=args.legacy,
    )
    if args.module:
        symbols = [item for item in symbols if item.module.startswith(args.module)]

    if args.list_modules:
        counts: dict[str, int] = {}
        for item in symbols:
            counts[item.module] = counts.get(item.module, 0) + 1
        for module, count in sorted(counts.items()):
            print(f"{module}\t{count}")
        return 0

    if args.query:
        if args.regex:
            try:
                pattern = re.compile(args.query, re.IGNORECASE)
            except re.error as exc:
                parser.error(f"invalid regex: {exc}")
            matches = lambda value: bool(pattern.search(value))
        else:
            needle = args.query.casefold()
            matches = lambda value: needle in value.casefold()
        symbols = [
            item
            for item in symbols
            if matches(" ".join((item.module, item.qualname, item.signature, item.doc)))
        ]

    for item in symbols:
        suffix = f"({item.signature})" if item.kind != "class" else ""
        location = item.path.relative_to(root)
        print(f"{item.module}:{item.qualname}{suffix}  [{location}:{item.line}]")
        if item.doc:
            print(f"  {item.doc}")
    return 0 if symbols else 1


if __name__ == "__main__":
    raise SystemExit(main())
