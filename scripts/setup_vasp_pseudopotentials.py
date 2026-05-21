#!/usr/bin/env python3
"""Prepare a licensed VASP PAW library for MCP_IC_Tool.

The tool expects:

    PSEUDO_PATH/
      PAWPickList
      paw_pbe/<label>/POTCAR

This script ingests an already licensed VASP potential directory or tarball,
then creates that layout using symlinks by default.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import tarfile
from pathlib import Path


ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)($|[_-])")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare licensed VASP PAW potentials")
    parser.add_argument("--source", required=True, help="Licensed VASP potential dir or .tar/.tar.gz/.tgz")
    parser.add_argument("--dest", default="~/pseudopotential", help="Destination PSEUDO_PATH")
    parser.add_argument("--family", default="paw_pbe", help="Destination family directory name")
    parser.add_argument("--copy", action="store_true", help="Copy POTCAR files instead of symlinking")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing family directory and PAWPickList",
    )
    return parser.parse_args()


def expand(path: str) -> Path:
    return Path(path).expanduser().resolve()


def extract_if_needed(source: Path, dest: Path) -> Path:
    if source.is_dir():
        return source
    if not tarfile.is_tarfile(source):
        raise SystemExit(f"source is not a directory or tar archive: {source}")

    extract_root = dest / "_official_sources" / source.stem.replace(".tar", "")
    if extract_root.exists():
        return extract_root
    extract_root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(source) as archive:
        archive.extractall(extract_root)
    return extract_root


def find_potential_root(root: Path) -> Path:
    candidates: list[tuple[int, Path]] = []
    for path in (root, *root.rglob("*")):
        if not path.is_dir():
            continue
        count = sum(1 for child in path.iterdir() if child.is_dir() and (child / "POTCAR").is_file())
        if count:
            candidates.append((count, path))
    if not candidates:
        raise SystemExit(
            "No VASP potential root found. Expected element directories containing POTCAR files. "
            "If your files are POTCAR.Z, uncompress them first."
        )
    candidates.sort(reverse=True, key=lambda item: item[0])
    return candidates[0][1]


def element_from_label(label: str) -> str | None:
    match = ELEMENT_RE.match(label)
    if not match:
        return None
    return match.group(1)


def choose_mapping(labels: list[str]) -> dict[str, str]:
    by_element: dict[str, list[str]] = {}
    for label in labels:
        element = element_from_label(label)
        if element:
            by_element.setdefault(element, []).append(label)

    mapping: dict[str, str] = {}
    for element, element_labels in sorted(by_element.items()):
        element_labels = sorted(element_labels)
        if element in element_labels:
            mapping[element] = element
        else:
            mapping[element] = element_labels[0]
    return mapping


def link_or_copy(source: Path, dest: Path, *, copy: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if copy:
        shutil.copy2(source, dest)
    else:
        os.symlink(source, dest)


def main() -> None:
    args = parse_args()
    source = expand(args.source)
    dest = expand(args.dest)
    family_dir = dest / args.family

    if not source.exists():
        raise SystemExit(f"source does not exist: {source}")
    if family_dir.exists() and not args.force:
        raise SystemExit(f"destination already exists: {family_dir}; use --force to replace")

    dest.mkdir(parents=True, exist_ok=True)
    root = extract_if_needed(source, dest)
    potential_root = find_potential_root(root)

    if family_dir.exists():
        shutil.rmtree(family_dir)
    family_dir.mkdir(parents=True)

    labels = sorted(
        child.name
        for child in potential_root.iterdir()
        if child.is_dir() and (child / "POTCAR").is_file()
    )
    for label in labels:
        link_or_copy(potential_root / label / "POTCAR", family_dir / label / "POTCAR", copy=args.copy)

    mapping = choose_mapping(labels)
    pick_list = dest / "PAWPickList"
    if pick_list.exists() and not args.force:
        raise SystemExit(f"PAWPickList already exists: {pick_list}; use --force to replace")
    with pick_list.open("w", encoding="utf-8") as handle:
        for element, label in mapping.items():
            handle.write(f"{element} {label}\n")

    print(f"source root: {potential_root}")
    print(f"destination: {dest}")
    print(f"family: {family_dir}")
    print(f"potentials linked: {len(labels)}")
    print(f"PAWPickList entries: {len(mapping)}")
    print("Set PSEUDO_PATH to:", dest)


if __name__ == "__main__":
    main()
