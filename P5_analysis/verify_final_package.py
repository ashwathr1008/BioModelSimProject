"""Verify hashes of the frozen final P5 report/PPT artifacts (standard library)."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "FINAL_FREEZE_MANIFEST.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def verify() -> list[str]:
    if not MANIFEST.is_file():
        return ["Missing FINAL_FREEZE_MANIFEST.json"]
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    issues = []
    for row in manifest["files"]:
        rel = Path(row["path"])
        if rel.is_absolute() or ".." in rel.parts:
            issues.append(f"Unsafe manifest path: {rel}")
            continue
        file = ROOT / rel
        if not file.is_file():
            issues.append(f"MISSING: {rel}")
        elif file.stat().st_size != row["bytes"]:
            issues.append(f"SIZE CHANGED: {rel}")
        elif sha256(file) != row["sha256"]:
            issues.append(f"HASH CHANGED: {rel}")
    return issues


if __name__ == "__main__":
    errors = verify()
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        sys.exit(1)
    n = len(json.loads(MANIFEST.read_text(encoding="utf-8"))["files"])
    print(f"PASS: {n} frozen P5 files match checksums")
