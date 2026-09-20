# Reproducibility guide

## Repository policy

This repository contains custom code and small text configuration files only.
The following files are excluded because they are large, generated, or contain
simulation-specific topology information:

- GROMACS trajectories and run files (`*.xtc`, `*.trr`, `*.tpr`, `*.edr`,
  `*.cpt` and `*.gro`).
- PyMOL sessions, rendered figures and spreadsheet exports.
- gmx_MMPBSA result archives and temporary Amber/GROMACS products.

For publication, archive the source-data tables, representative PDB files and
any minimally necessary processed trajectories in a research-data repository
and cite its DOI separately from the code DOI.

## Standard directory layout

The plotting scripts accept an external data root. A compatible local layout is:

```text
data/
  CC_1/
    disulfide_geometry/
    mmpbsa_ABC_vs_DEF/
  mC/
    disulfide_geometry/
    mmpbsa_ABC_vs_DEF/
  NC/
    disulfide_geometry/
    mmpbsa_ABC_vs_DEF/
  P10/
    mmpbsa_ABC_vs_DEF/
  cys_screening/
  rmsd_rg_fel/
  ensemble_dccm/
    source_data/
    structures/
```

Set the root once per shell:

```bash
export COLLAGEN_DATA_ROOT=/absolute/path/to/data
export COLLAGEN_RESULTS_ROOT=/absolute/path/to/results
```

On Windows PowerShell:

```powershell
$env:COLLAGEN_DATA_ROOT = (Resolve-Path ".\data").Path
$env:COLLAGEN_RESULTS_ROOT = (Resolve-Path ".\results").Path
```

## Module-specific variables

| Variable | Purpose | Default |
|---|---|---|
| `COLLAGEN_DATA_ROOT` | Primary MD/MM-PBSA input root | `data/` |
| `COLLAGEN_RESULTS_ROOT` | Generated figures and tables | `results/` |
| `COLLAGEN_CACHE_DIR` | ASCII-only temporary trajectory cache | OS temporary directory |
| `COLLAGEN_CYS_SCREENING_DATA` | PyRosetta candidate tables | `data/cys_screening/` |
| `COLLAGEN_FEL_SOURCE` | Shared-grid RMSD-Rg CSV files | `data/rmsd_rg_fel/` |
| `COLLAGEN_FEL_TIMESERIES` | Per-system RMSD/Rg time series | `data/rmsd_rg_fel/` |
| `COLLAGEN_DCCM_SOURCE` | DCCM matrices and manifest | `data/ensemble_dccm/` |
| `COLLAGEN_ENSEMBLE_STRUCTURES` | Ensemble and putty PDBs | `data/ensemble_dccm/structures/` |
| `COLLAGEN_MMPBSA_SOURCE` | Consolidated MM/PBSA source tables | `data/mmpbsa_analysis/` |
| `COLLAGEN_DOCKING_SOURCE` | Starting docking PDBs | `data/docking_structures/` |
| `COLLAGEN_PYMOL_SESSIONS` | Optional local PyMOL sessions | `data/pymol_sessions/` |
| `PYMOL_BIN` | PyMOL executable | `pymol` on `PATH` |

## Trajectory conventions

`repair_docking_md.sh` treats all six protein chains as one complex during PBC
clustering and global fitting. It produces:

- `full_clustered.xtc`: complete solvated system for interface hydration,
  water bridges and sulfur-coordinate extraction.
- `protein_complex_fit.xtc`: fitted protein-only trajectory for RMSD, PCA,
  RMSF, DCCM and conformational ensembles.
- `reference_protein.gro`: the repaired protein reference at 0 ps.

The general docking analysis uses 20--200 ns by default for equilibrium
structural analyses. The sulfur-proximity and MM/PBSA runners use 0--200 ns by
default, matching the displayed time-course panels. Override their start/end
variables when applying a different equilibration policy, and report the
chosen window explicitly.

## Six-chain index split

Chains A--C are treated as trimer ABC and chains D--F as trimer DEF. The
validated defaults in `analyze_docking_md.sh` are system specific. For a new
topology, inspect its zero-based GROMACS `resindex` range and set:

```bash
ABC_LAST=<last_resindex_of_ABC> \
DEF_FIRST=<first_resindex_of_DEF> \
DEF_LAST=<last_resindex_of_DEF> \
bash scripts/02_trajectory_preprocessing/analyze_docking_md.sh \
  /path/to/simulation /path/to/repaired
```

Always confirm group atom counts before calculating RMSD, PCA or MM/PBSA.

## PyMOL recipes

Run repository-relative PML files from the repository root so their relative
`data/` and `results/` paths resolve correctly. The Python PyMOL renderers can
instead be directed to external structures through the variables above.

## Numerical interpretation

- Free-energy maps are projected relative free energies derived from sampled
  RMSD/Rg occupancies; they are not absolute folding free energies.
- DCCM entries range from -1 (anti-correlated) through 0 (weakly correlated) to
  +1 (correlated).
- An SG-SG distance below 0.40 nm is labelled a proximity event only.
- MM/PBSA is calculated as complex minus receptor minus ligand. The displayed
  total does not include configurational entropy.
