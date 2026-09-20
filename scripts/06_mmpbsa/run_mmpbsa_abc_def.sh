#!/usr/bin/env bash
# Run single-trajectory MM/PBSA for a six-chain docking system.
# Receptor: trimer ABC. Ligand: trimer DEF.
# The analysis assumes that the repaired full_clustered.xtc trajectory is whole,
# PBC-corrected and contains the full simulated system.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    cat <<'USAGE'
Usage:
  bash run_mmpbsa_abc_def.sh <simulation_dir> <repaired_dir> [output_dir]

Example:
  bash run_mmpbsa_abc_def.sh \
    /path/to/simulation/CC_1 \
    /path/to/repaired/CC_1

Environment overrides:
  GMX_BIN=gmx                    GROMACS executable (default: gmx)
  MMPBSA_BIN=gmx_MMPBSA          gmx_MMPBSA executable
  MPI_RANKS=1                    Set >1 to use mpirun
  START_NS=0 END_NS=200            Production analysis window
  FRAME_INTERVAL_NS=1             Sampling interval after repair
  NDX_FILE=/path/docking_groups.ndx
  MMPBSA_INPUT=/path/mmpbsa_abc_def.in
  USE_REFERENCE_PDB=0             Set to 1 only with a residue-matched reference PDB
  REFERENCE_PDB=/path/complex_reference.pdb
USAGE
    exit 2
fi

SIM_DIR=$(cd "$1" && pwd)
REPAIRED_DIR=$(cd "$2" && pwd)
OUT_DIR=${3:-"$REPAIRED_DIR/mmpbsa_ABC_vs_DEF"}
OUT_DIR=$(mkdir -p "$OUT_DIR" && cd "$OUT_DIR" && pwd)

GMX_BIN=${GMX_BIN:-gmx}
MMPBSA_BIN=${MMPBSA_BIN:-gmx_MMPBSA}
MPI_RANKS=${MPI_RANKS:-1}
START_NS=${START_NS:-0}
END_NS=${END_NS:-200}
FRAME_INTERVAL_NS=${FRAME_INTERVAL_NS:-1}
USE_REFERENCE_PDB=${USE_REFERENCE_PDB:-0}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MMPBSA_INPUT=${MMPBSA_INPUT:-"$SCRIPT_DIR/../../config/mmpbsa_abc_def.in"}

TPR="$SIM_DIR/md_prod.tpr"
TOP="$SIM_DIR/topol.top"
TRAJ="$REPAIRED_DIR/full_clustered.xtc"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -f "$TPR" ]] || die "missing run-input file: $TPR"
[[ -f "$TOP" ]] || die "missing topology file: $TOP"
[[ -f "$TRAJ" ]] || die "missing repaired full-system trajectory: $TRAJ"
[[ -f "$MMPBSA_INPUT" ]] || die "missing gmx_MMPBSA input file: $MMPBSA_INPUT"
command -v "$GMX_BIN" >/dev/null 2>&1 || die "cannot find GROMACS executable: $GMX_BIN"
command -v "$MMPBSA_BIN" >/dev/null 2>&1 || die "cannot find gmx_MMPBSA executable: $MMPBSA_BIN"

if [[ -z "${AMBERHOME:-}" ]]; then
    die "AMBERHOME is unset. Activate the gmx_MMPBSA conda environment and run: export AMBERHOME=\$CONDA_PREFIX"
fi
[[ -d "$AMBERHOME" ]] || die "AMBERHOME does not exist: $AMBERHOME"

if [[ -n "${NDX_FILE:-}" ]]; then
    NDX_FILE=$(cd "$(dirname "$NDX_FILE")" && pwd)/$(basename "$NDX_FILE")
else
    NDX_FILE=$(find "$REPAIRED_DIR" -type f -name 'docking_groups.ndx' -print | sort -V | head -n 1 || true)
fi
[[ -n "$NDX_FILE" && -f "$NDX_FILE" ]] || die "docking_groups.ndx was not found. Set NDX_FILE=/absolute/path/to/docking_groups.ndx"

index_group_id() {
    local index_file=$1
    local requested_name=$2
    awk -v requested_name="$requested_name" '
        BEGIN { group_id = 0 }
        /^\[/ {
            group_name = $0
            sub(/^\[[[:space:]]*/, "", group_name)
            sub(/[[:space:]]*\][[:space:]]*$/, "", group_name)
            if (group_name == requested_name) {
                print group_id
                exit
            }
            group_id++
        }
    ' "$index_file"
}

list_groups() {
    awk '/^\[/ { line=$0; sub(/^\[[[:space:]]*/, "", line); sub(/[[:space:]]*\][[:space:]]*$/, "", line); print line }' "$1"
}

COMPLEX_GROUP=$(index_group_id "$NDX_FILE" "COMPLEX_Protein")
ABC_GROUP=$(index_group_id "$NDX_FILE" "TRIMER_ABC")
DEF_GROUP=$(index_group_id "$NDX_FILE" "TRIMER_DEF")

if [[ -z "$COMPLEX_GROUP" || -z "$ABC_GROUP" || -z "$DEF_GROUP" ]]; then
    echo "Available groups in $NDX_FILE:" >&2
    list_groups "$NDX_FILE" >&2
    die "the index must contain COMPLEX_Protein, TRIMER_ABC and TRIMER_DEF"
fi

if grep -qi 'charmm' "$TOP"; then
    if grep -Eq '^[[:space:]]*PBRadii[[:space:]]*=[[:space:]]*7' "$MMPBSA_INPUT" \
        && grep -Eq '^[[:space:]]*radiopt[[:space:]]*=[[:space:]]*0' "$MMPBSA_INPUT"; then
        echo "[INFO] CHARMM topology detected; using CHARMM PB radii (PBRadii=7, radiopt=0)."
    else
        echo "WARNING: a CHARMM topology was detected, but the input does not specify PBRadii=7 and radiopt=0." >&2
        echo "         Update $MMPBSA_INPUT before treating PB values as final." >&2
    fi
fi

SAMPLED_TRAJ="$OUT_DIR/complex_${START_NS}-${END_NS}ns_every_${FRAME_INTERVAL_NS}ns.xtc"
if [[ ! -s "$SAMPLED_TRAJ" ]]; then
    echo "[1/3] Sampling repaired full-system trajectory: ${START_NS}-${END_NS} ns, every ${FRAME_INTERVAL_NS} ns"
    # No custom index is supplied: default group 0 is System, preserving the atom count of md_prod.tpr.
    printf '0\n' | "$GMX_BIN" trjconv \
        -s "$TPR" \
        -f "$TRAJ" \
        -o "$SAMPLED_TRAJ" \
        -b "$START_NS" \
        -e "$END_NS" \
        -dt "$FRAME_INTERVAL_NS" \
        -tu ns
else
    echo "[1/3] Reusing existing sampled trajectory: $SAMPLED_TRAJ"
fi

if (( USE_REFERENCE_PDB == 1 )); then
    if [[ -n "${REFERENCE_PDB:-}" ]]; then
        REFERENCE_PDB=$(cd "$(dirname "$REFERENCE_PDB")" && pwd)/$(basename "$REFERENCE_PDB")
    else
        REFERENCE_PDB="$OUT_DIR/complex_reference_${START_NS}ns.pdb"
    fi

    if [[ ! -s "$REFERENCE_PDB" ]]; then
        echo "[2/3] Creating an ABC+DEF protein reference structure for residue and chain mapping"
        printf '%s\n' "$COMPLEX_GROUP" | "$GMX_BIN" trjconv \
            -s "$TPR" \
            -f "$TRAJ" \
            -n "$NDX_FILE" \
            -o "$REFERENCE_PDB" \
            -b "$START_NS" \
            -e "$START_NS" \
            -tu ns
    else
        echo "[2/3] Reusing existing reference PDB: $REFERENCE_PDB"
    fi
else
    echo "[2/3] No external reference PDB: using the complex-derived single-trajectory mapping"
fi

RESULT_DAT="$OUT_DIR/FINAL_RESULTS_MMPBSA.dat"
RESULT_CSV="$OUT_DIR/FINAL_RESULTS_MMPBSA.csv"
DECOMP_DAT="$OUT_DIR/FINAL_DECOMP_MMPBSA.dat"
DECOMP_CSV="$OUT_DIR/FINAL_DECOMP_MMPBSA.csv"
LOG_FILE="$OUT_DIR/gmx_MMPBSA.log"

CMD=(
    "$MMPBSA_BIN" -O -nogui
    -i "$MMPBSA_INPUT"
    -cs "$TPR"
    -ct "$SAMPLED_TRAJ"
    -ci "$NDX_FILE"
    -cg "$ABC_GROUP" "$DEF_GROUP"
    -cp "$TOP"
    -o "$RESULT_DAT"
    -eo "$RESULT_CSV"
    -do "$DECOMP_DAT"
    -deo "$DECOMP_CSV"
)

if (( USE_REFERENCE_PDB == 1 )); then
    CMD+=( -cr "$REFERENCE_PDB" )
fi

echo "[3/3] Running single-trajectory MM/PBSA"
echo "        receptor = TRIMER_ABC (index $ABC_GROUP); ligand = TRIMER_DEF (index $DEF_GROUP)"
echo "        selected complex = COMPLEX_Protein (index $COMPLEX_GROUP)"
echo "        input = $MMPBSA_INPUT"

if (( MPI_RANKS > 1 )); then
    command -v mpirun >/dev/null 2>&1 || die "MPI_RANKS=$MPI_RANKS but mpirun is not available"
    python -c 'from mpi4py import MPI; print(MPI.Get_library_version())' >/dev/null 2>&1 || die "MPI_RANKS=$MPI_RANKS but mpi4py is unavailable"
    mpirun -np "$MPI_RANKS" "${CMD[@]}" 2>&1 | tee "$LOG_FILE"
else
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE"
fi

cat <<EOF

Completed: $OUT_DIR
Primary outputs:
  FINAL_RESULTS_MMPBSA.dat/.csv       Complex, receptor, ligand and binding terms
  FINAL_DECOMP_MMPBSA.dat/.csv        Per-residue decomposition within 6 A
  gmx_MMPBSA.log                      Full run log
  complex_*ns_every_*ns.xtc           Sampled PBC-corrected trajectory used for analysis

Interpretation:
  DeltaG_bind = G(complex) - G(trimer ABC) - G(trimer DEF).
  Compare systems using the same window, salt concentration, radii settings and frame interval.
  These are endpoint MM/PBSA estimates, not a direct covalent S-S bond-formation energy.
EOF
