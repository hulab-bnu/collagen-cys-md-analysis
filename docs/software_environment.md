# Software environment

## Recorded primary software

- GROMACS 2024.5.
- Python 3.11 for the gmx_MMPBSA analysis environment.

## Python dependencies

The open Python dependencies are listed in `environment.yml` and
`requirements.txt`. Exact historical package build identifiers were not
preserved in the analysis records; the repository therefore does not claim an
unverified lockfile. For a final archival release, run the commands below in
the validated environment and retain the outputs with the release assets:

```bash
python --version
python -m pip freeze > software_versions_pip.txt
gmx --version > software_versions_gromacs.txt
gmx_MMPBSA --version > software_versions_gmx_mmpbsa.txt
pymol -cq -d "print(cmd.get_version()); quit" > software_versions_pymol.txt
```

PyRosetta is distributed separately and must be installed under its own
academic or commercial licence. Record the PyRosetta build string printed by
the screening run in the archived analysis log.

## Platform notes

GROMACS shell workflows target a POSIX shell. Python plotting and table scripts
are platform independent after configuring data paths. PyMOL render appearance
can vary slightly by version and graphics backend; the PML files preserve the
camera-independent molecular styling and output dimensions.
