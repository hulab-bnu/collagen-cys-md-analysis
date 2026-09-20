# MM/PBSA configuration

- `mmpbsa_abc_def.in` is the input used by the six-chain ABC-versus-DEF runner.
- `mmpbsa_charmm_pb_general.in` is the corresponding general CHARMM/PB setup.

The calculation uses the receptor and ligand group names resolved from the
GROMACS index. Confirm ionic strength, PB radii and decomposition settings
against the manuscript before reusing these files for another system.
