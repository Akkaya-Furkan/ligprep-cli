import subprocess

import pytest

from ligprep import convert
from ligprep.errors import LigprepError
from ligprep.validate import pdbqt_stats

from .conftest import needs_obabel


@needs_obabel
def test_aspirin_converts_to_a_valid_pdbqt(aspirin_sdf):
    stats = pdbqt_stats(convert.sdf_to_pdbqt(aspirin_sdf))
    assert stats["heavy"] == 13  # 21 atoms in the SDF; non-polar hydrogens are merged
    assert stats["torsdof"] is not None
    assert abs(stats["charge"]) < 0.05


@needs_obabel
def test_garbage_input_is_a_ligprep_error():
    with pytest.raises(LigprepError, match="obabel failed"):
        convert.sdf_to_pdbqt("not a molecule")


def test_missing_obabel_is_reported(monkeypatch):
    monkeypatch.setattr(convert.shutil, "which", lambda name: None)
    with pytest.raises(LigprepError, match="Open Babel"):
        convert.find_obabel()


def test_obabel_timeout_is_a_ligprep_error(monkeypatch, aspirin_sdf):
    def boom(*a, **kw):
        raise subprocess.TimeoutExpired(cmd="obabel", timeout=120)
    monkeypatch.setattr(convert.subprocess, "run", boom)
    with pytest.raises(LigprepError, match="timed out"):
        convert.sdf_to_pdbqt(aspirin_sdf, obabel="/bin/true")
