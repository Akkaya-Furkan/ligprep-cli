# Ligand preparation CLI

**Question.** A command-line tool that turns names, CIDs or SMILES into docking-ready PDBQT files — with built-in checks.

## Install
Requires Python >= 3.10 and [Open Babel](https://openbabel.org/) (`obabel` on `PATH`).
```bash
git clone https://github.com/Akkaya-Furkan/ligprep-cli && cd ligprep-cli
python -m venv .venv && . .venv/bin/activate
pip install -e .
```

## Usage
```bash
ligprep --input examples/drugs.txt --out out/
```
One ligand per line; `#` starts a comment. Bare text is a name, a bare integer is a PubChem CID; prefixes
`name:`, `cid:` and `smiles:` make it explicit.
```
aspirin
cid:3672
smiles:CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```
Output in `out/`: one `<label>_CID<cid>.pdbqt` per valid ligand and `manifest.csv` (input, CID, SMILES, formula,
heavy atoms, TORSDOF, charge sum, status, messages). Options: `--keep-sdf`, `--allow-flat`, `--delay`.
Exit code: 0 all prepared, 1 at least one failed, 2 bad input or Open Babel missing.

## Pipeline and checks
name / CID / SMILES -> PubChem CID -> 3D SDF -> Open Babel (Gasteiger charges) -> PDBQT. A PDBQT is written only if:

| Check | Fails when |
| --- | --- |
| Formula | SDF formula differs from PubChem's molecular formula |
| 3D | every z coordinate is 0 (a 2D structure; skip with `--allow-flat`) |
| Heavy atoms | PDBQT heavy-atom count differs from the formula |
| TORSDOF | line missing (warning above 20 torsions) |
| Charges | warning if the Gasteiger sum is more than 0.1 from the formal charge |

Failed ligands are reported in `manifest.csv` and do not stop the run.

## Data
PubChem PUG REST 3D conformers (no key needed, requests are rate-limited to 4 per second). Examples use public
drugs: aspirin, ibuprofen, caffeine.

## Tests
```bash
pip install -e ".[dev]"
pytest -q        # 50 tests, no network: PubChem is mocked with `responses`
ruff check .
```
Open Babel must be installed for the conversion and end-to-end tests (they skip without it). CI runs both on
Python 3.10 and 3.12. The suite covers input parsing, PubChem lookups (retries, 404, network errors, the
User-Agent), each validation check including a formula mismatch and a flat structure, the Open Babel step and
the CLI exit codes. Each check was also switched off by hand to confirm a test fails without it.

## Checked by hand
`ligprep --input examples/drugs.txt --out out/` prepares 3/3 molecules against live PubChem. The three PDBQT
files were also read and docked by AutoDock Vina 1.2.5 (a parse check only: exhaustiveness 1, so the scores mean
nothing).

## Limitations
- Quality is that of PubChem's computed 3D conformer; no protonation at a given pH and no tautomer enumeration.
- Open Babel merges non-polar hydrogens, as AutoDock expects; the TORSDOF counts rotatable bonds the way Open Babel does.
- Ligands are prepared one at a time; PubChem rate limits make very large lists slow.

## Status
v0.1.0: packaging, tests and CI done.
Part of a 30-project computational biology and scientific software portfolio.
