"""Parsing the input list: one ligand per line, '#' starts a comment."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .errors import LigprepError

KINDS = ("name", "cid", "smiles")


@dataclass(frozen=True)
class Entry:
    kind: str  # name | cid | smiles
    value: str
    raw: str

    @property
    def label(self) -> str:
        base = self.value if self.kind != "smiles" else "smiles"
        return re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("_")[:40] or "ligand"


def parse_line(line: str) -> Entry | None:
    text = line.split("#", 1)[0].strip()
    if not text:
        return None
    prefix, sep, rest = text.partition(":")
    if sep and prefix.lower() in KINDS and rest.strip():
        kind, value = prefix.lower(), rest.strip()
    elif text.isdigit():
        kind, value = "cid", text
    else:
        kind, value = "name", text
    if kind == "cid" and not value.isdigit():
        raise LigprepError(f"CID must be an integer: {value!r}")
    return Entry(kind, value, text)


def read_list(path: str | Path) -> list[Entry]:
    entries = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        try:
            entry = parse_line(line)
        except LigprepError as e:
            raise LigprepError(f"{path}:{n}: {e}") from None
        if entry:
            entries.append(entry)
    return entries
