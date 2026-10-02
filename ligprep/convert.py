"""SDF -> PDBQT with Open Babel (Gasteiger charges)."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from .errors import LigprepError


def find_obabel() -> str:
    exe = shutil.which("obabel")
    if not exe:
        raise LigprepError("Open Babel ('obabel') not found on PATH; install it first")
    return exe


def sdf_to_pdbqt(sdf_text: str, obabel: str | None = None) -> str:
    """Return PDBQT text. Hydrogens in the 3D conformer are kept; non-polar ones are merged by obabel."""
    exe = obabel or find_obabel()
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "in.sdf", Path(d) / "out.pdbqt"
        src.write_text(sdf_text)
        cmd = [exe, str(src), "-O", str(dst), "--partialcharge", "gasteiger"]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)
        except subprocess.TimeoutExpired:
            raise LigprepError("obabel timed out after 120 s") from None
        if proc.returncode != 0 or not dst.exists() or dst.stat().st_size == 0:
            raise LigprepError(f"obabel failed: {(proc.stderr or proc.stdout).strip()[:200]}")
        return dst.read_text()
