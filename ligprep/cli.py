"""Command line entry point: ligprep --input list.txt --out out/"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__
from .convert import find_obabel, sdf_to_pdbqt
from .errors import LigprepError
from .inputs import Entry, read_list
from .pubchem import PubChem
from .validate import read_sdf, validate

MANIFEST_FIELDS = ["input", "kind", "cid", "smiles", "formula", "heavy_atoms", "torsdof",
                   "charge_sum", "pdbqt", "status", "messages"]


def prepare(entry: Entry, client: PubChem, out: Path, obabel: str, allow_flat: bool,
            keep_sdf: bool) -> dict:
    row = {k: "" for k in MANIFEST_FIELDS}
    row.update(input=entry.raw, kind=entry.kind, status="failed")
    try:
        cid = client.cid_for(entry)
        row["cid"] = cid
        props = client.properties(cid)
        row["smiles"] = props.get("SMILES", "")
        sdf = client.sdf_3d(cid)
        mol = read_sdf(sdf)
        pdbqt = sdf_to_pdbqt(sdf, obabel)
        rep = validate(mol, props["MolecularFormula"], pdbqt, allow_flat)
        row.update(formula=rep.formula, heavy_atoms=rep.heavy_atoms, torsdof=rep.torsdof,
                   charge_sum=rep.charge_sum)
        row["messages"] = "; ".join([f"ERROR: {m}" for m in rep.errors] + rep.warnings)
        if rep.ok:
            name = f"{entry.label}_CID{cid}"
            (out / f"{name}.pdbqt").write_text(pdbqt)
            if keep_sdf:
                (out / f"{name}.sdf").write_text(sdf)
            row.update(pdbqt=f"{name}.pdbqt", status="ok")
    except LigprepError as e:
        row["messages"] = f"ERROR: {e}"
    return row


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ligprep", description=__doc__)
    ap.add_argument("--input", "-i", required=True, help="text file, one name, CID or SMILES per line")
    ap.add_argument("--out", "-o", required=True, help="output directory")
    ap.add_argument("--allow-flat", action="store_true", help="do not fail on planar (z = 0) structures")
    ap.add_argument("--keep-sdf", action="store_true", help="also keep the downloaded 3D SDF files")
    ap.add_argument("--delay", type=float, default=0.25, help="seconds between PubChem requests")
    ap.add_argument("--version", action="version", version=f"ligprep {__version__}")
    a = ap.parse_args(argv)

    try:
        entries = read_list(a.input)
        obabel = find_obabel()
    except (LigprepError, OSError) as e:
        print(f"ligprep: {e}", file=sys.stderr)
        return 2
    if not entries:
        print("ligprep: input list is empty", file=sys.stderr)
        return 2

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    client = PubChem(delay=a.delay)
    rows = []
    for e in entries:
        row = prepare(e, client, out, obabel, a.allow_flat, a.keep_sdf)
        rows.append(row)
        tag = "ok  " if row["status"] == "ok" else "FAIL"
        extra = f"  {row['messages']}" if row["messages"] else ""
        print(f"[{tag}] {e.raw}  CID {row['cid'] or '-'}  TORSDOF {row['torsdof'] if row['torsdof'] != '' else '-'}{extra}")

    with open(out / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_ok = sum(r["status"] == "ok" for r in rows)
    print(f"{n_ok}/{len(rows)} ligands prepared -> {out}/ (manifest.csv)")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
