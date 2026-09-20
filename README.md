# Collagen Cys MD Analysis

Analysis and figure-generation code for position-resolved cysteine screening and
six-chain collagen docking molecular dynamics. The repository is organised by
figure panel and contains the custom code, PyMOL rendering recipes and
MM/PBSA configuration files used in the study.

## Scope

This is a code-only release. Large trajectories, GROMACS run files, intermediate
MM/PBSA files and rendered figures are intentionally excluded from GitHub. The
expected input layout is documented in
[`docs/reproducibility.md`](docs/reproducibility.md). Source data should be
archived separately in a suitable research-data repository when the manuscript
is published.

## Figure modules

| Panel | Analysis | Code |
|---|---|---|
| a | Position-resolved Cys permissibility landscape | `scripts/01_cys_screening/` |
| b | Conformational ensembles and Cys mapping | `scripts/04_ensemble_dccm/`, `pymol/` |
| c | RMSD-Rg free-energy landscapes | `scripts/03_rmsd_rg_fel/` |
| d | Dynamic cross-correlation matrices | `scripts/04_ensemble_dccm/` |
| e | Nearest inter-trimer S-gamma--S-gamma distance | `scripts/05_sg_distance/` |
| f | Frame-wise MM/PBSA association-energy estimates | `scripts/06_mmpbsa/` |
| g | Representative docking conformations and local Cys interfaces | `scripts/07_interface_rendering/`, `pymol/` |
| h | MM/PBSA association-energy component profile | `scripts/06_mmpbsa/` |

The detailed script-to-panel map is in
[`docs/figure_code_map.md`](docs/figure_code_map.md).

## Requirements

The analysis was developed around:

- GROMACS 2024.5 for trajectory preprocessing and coordinate extraction.
- Python 3.11 with NumPy, pandas, Matplotlib, SciPy, Pillow and MDAnalysis.
- gmx_MMPBSA and AmberTools for endpoint MM/PBSA calculations.
- PyMOL for structural rendering.
- PyRosetta for the cysteine-placement screen.

Install the open Python dependencies with either:

```bash
conda env create -f environment.yml
conda activate collagen-cys-md
```

or:

```bash
python -m pip install -r requirements.txt
```

GROMACS, gmx_MMPBSA, PyMOL and PyRosetta should be installed separately under
their respective licences. See
[`docs/software_environment.md`](docs/software_environment.md).

## Configuration

No script contains a user-specific filesystem path. The principal locations can
be set with environment variables:

```bash
export COLLAGEN_DATA_ROOT=/path/to/input-data
export COLLAGEN_RESULTS_ROOT=/path/to/results
export COLLAGEN_CACHE_DIR=/path/to/cache
export PYMOL_BIN=/path/to/pymol
```

More specialised variables are documented beside the relevant commands in
[`docs/reproducibility.md`](docs/reproducibility.md). By default, scripts look
for inputs under `data/` and write generated files under `results/`; both
directories are excluded from version control except for their README files.

## Core workflows

Repair a six-chain trajectory and calculate docking-MD metrics:

```bash
bash scripts/02_trajectory_preprocessing/repair_docking_md.sh \
  /path/to/simulation /path/to/repaired

bash scripts/02_trajectory_preprocessing/analyze_docking_md.sh \
  /path/to/simulation /path/to/repaired
```

Calculate all cross-trimer cysteine sulfur distances over 0--200 ns:

```bash
START_NS=0 END_NS=200 STRIDE_NS=0.1 \
bash scripts/05_sg_distance/analyze_sg_sg_geometry.sh \
  /path/to/simulation /path/to/repaired
```

Run ABC-versus-DEF single-trajectory MM/PBSA over 0--200 ns:

```bash
START_NS=0 END_NS=200 FRAME_INTERVAL_NS=1 \
bash scripts/06_mmpbsa/run_mmpbsa_abc_def.sh \
  /path/to/simulation /path/to/repaired
```

The MM/PBSA runner resolves `COMPLEX_Protein`, `TRIMER_ABC` and `TRIMER_DEF`
from the GROMACS index by name. Its default input file is
[`config/mmpbsa_abc_def.in`](config/mmpbsa_abc_def.in).

## Interpretation boundaries

- S-gamma--S-gamma distances below 0.40 nm are geometric proximity events and
  do not establish formation of a covalent disulfide bond.
- MM/PBSA values are endpoint association-energy estimates. Configurational
  entropy is not included, and trajectory frames are not independent biological
  replicates.
- DCCM and free-energy landscapes are descriptive outputs from one trajectory
  per construct unless independent repeats are supplied separately.

## Citation

Use the repository citation provided in [`CITATION.cff`](CITATION.cff). Before
article publication, create a versioned GitHub release, archive it with Zenodo
and add the resulting DOI to both the manuscript and `CITATION.cff`.

## Licence

The custom code in this repository is released under the
[MIT License](LICENSE). Third-party software and input data retain their own
licences and are not redistributed here.
