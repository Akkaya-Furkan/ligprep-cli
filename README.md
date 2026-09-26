# Ligand preparation CLI

**Question.** A command-line tool that turns names, CIDs or SMILES into docking-ready PDBQT files — with built-in checks.

## Data
PubChem PUG REST 3D conformers; examples use public drugs (aspirin, ibuprofen, caffeine).

## Approach
1. name/CID/SMILES → PubChem 3D → Open Babel Gasteiger → PDBQT
2. Built-in validation: formula match, non-flat 3D (z ≠ 0), TORSDOF report, CID + SMILES manifest CSV
3. Offline tests with mocked HTTP, GitHub Actions (pytest + ruff), v0.1.0 release

## Deliverables
`ligprep --input list.txt --out out/` produces validated files; CI green; release published.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
