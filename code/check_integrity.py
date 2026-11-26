#!/usr/bin/env python3
"""Check the listed files, sizes and hashes in a release ZIP or Git clone."""
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
manifest = "FILE_MANIFEST.json"
records = json.loads((ROOT / manifest).read_text())
expected = {row["file"] for row in records}
if len(expected) != len(records):
    raise SystemExit("Duplicate paths in the release manifest")

actual = set()
for directory, subdirs, filenames in os.walk(ROOT):
    subdirs[:] = [name for name in subdirs if name not in {".git", "__pycache__"}]
    for name in filenames:
        relative = (Path(directory) / name).relative_to(ROOT).as_posix()
        if relative not in {manifest, ".git"} and name != ".DS_Store":
            actual.add(relative)
if expected != actual:
    raise SystemExit(f"Missing: {sorted(expected - actual)}; unexpected: {sorted(actual - expected)}")
for row in records:
    data = (ROOT / row["file"]).read_bytes()
    if len(data) != row["bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
        raise SystemExit(f"Content mismatch: {row['file']}")
print(f"Integrity verified for {len(expected)} files")
