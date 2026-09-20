# Figure-to-code map

This document maps each main-figure panel to the custom code retained in the
repository. Generated files and large numerical inputs are intentionally not
tracked by Git.

## Panel a: position-resolved Cys permissibility

- `scripts/01_cys_screening/pyrosetta_cys_centered_multisite_screen.py`
  performs the centered local PyRosetta screen and writes candidate-level
  metrics.
- `scripts/01_cys_screening/make_six_metric_cys_landscape.py` converts the
  tabulated screen into the position-resolved heat map.

## Panel b: conformational ensembles and Cys mapping

- `scripts/04_ensemble_dccm/make_md_panels_ab_ensemble_dccm.py` aligns sampled
  trajectory frames, calculates C-alpha RMSF and writes ensemble/putty PDBs.
- `scripts/04_ensemble_dccm/render_pymol_transparent_ensemble_plate.py` and
  `scripts/04_ensemble_dccm/render_pymol_medoid_cartoon_cloud.py` render the
  ensemble representations.
- `pymol/*_transparent_ensemble.pml` stores reusable PyMOL styling recipes.

## Panel c: RMSD-Rg free-energy landscapes

- `scripts/03_rmsd_rg_fel/make_md_fel_rmsd_rg.py` calculates the projected
  RMSD/Rg distributions and relative free energies.
- `scripts/03_rmsd_rg_fel/make_additional_md_fel.py` handles additional
  trajectories, including the P10 control.
- `scripts/03_rmsd_rg_fel/make_combined_p10_primary_fel.py` places systems on
  a shared grid.
- `scripts/03_rmsd_rg_fel/plot_rmsd_rg_energy_basins_screening_palette.py`
  generates the final harmonised panel.

## Panel d: dynamic cross-correlation matrices

- `scripts/04_ensemble_dccm/make_md_panels_ab_ensemble_dccm.py` calculates the
  C-alpha DCCM after common-core alignment.
- `scripts/04_ensemble_dccm/plot_dccm_screening_palette.py` renders the four
  matrices on a common -1 to +1 colour scale.

## Panel e: inter-trimer sulfur proximity

- `scripts/05_sg_distance/analyze_sg_sg_geometry.sh` selects all Cys SG atoms
  from a repaired full-atom trajectory.
- `scripts/05_sg_distance/sg_sg_geometry_from_xvg.py` computes all SG-SG pair
  distances and cross-trimer contact occupancies.
- `scripts/05_sg_distance/plot_full_trajectory_sg_sg_proximity.py` plots the
  per-frame nearest cross-trimer distance and its rolling-median version.

## Panels f and h: MM/PBSA

- `scripts/06_mmpbsa/run_mmpbsa_abc_def.sh` runs the ABC-versus-DEF endpoint
  calculation.
- `scripts/06_mmpbsa/run_gmx_mmpbsa_general.sh` is the group-name-driven
  general runner.
- `scripts/06_mmpbsa/analyze_plot_mmpbsa.py` extracts per-frame totals,
  gas-solvent compensation, block stability and component means.
- `scripts/06_mmpbsa/plot_mmpbsa_components_with_total.py` renders the final
  component heat map with the reported total.
- `scripts/06_mmpbsa/plot_sg_sg_mmpbsa_timecourse_pair.py` provides matched
  distance and MM/PBSA time-course panels.
- `config/` contains the PB/MM-PBSA input files.

## Panel g: representative docking conformations and local interfaces

- `scripts/07_interface_rendering/analyze_collagen_docking.py` audits chain
  lengths and inter-trimer contacts in six-chain PDB files.
- `scripts/07_interface_rendering/crop_collagen_complexes.py` creates defined
  interface-centred structural variants without coordinate transformation.
- `scripts/07_interface_rendering/render_local_panels.py` and the local
  `pymol/*.pml` recipes render the interface close-ups.
- `scripts/07_interface_rendering/make_p10_p10_docking_control.py` constructs
  the P10:P10 control while preserving the source rigid-body arrangement.

## Shared preprocessing and utilities

- `scripts/02_trajectory_preprocessing/repair_docking_md.sh` repairs PBC and
  creates full-system and fitted protein trajectories.
- `scripts/02_trajectory_preprocessing/analyze_docking_md.sh` generates RMSD,
  RMSF, Rg, PCA, interface, hydration and snapshot outputs.
- `scripts/08_figure_utilities/make_shared_energy_colorbar.py` creates the
  shared energy scale used during figure assembly.
