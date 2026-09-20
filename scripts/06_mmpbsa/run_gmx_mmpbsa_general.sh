#!/usr/bin/env bash
# Generic single-trajectory gmx_MMPBSA runner.
# Receptor and ligand are resolved by group name, not hard-coded group number.

set -euo pipefail

if [[ $# -ne 7 ]]; then
    cat <<'USAGE'
Usage:
  bash run_gmx_mmpbsa_general.sh \
    <simulation_dir> <repaired_full_trajectory> <index_file> \
    <receptor_group_name> <ligand_group_name> <mmpbsa_input> <output_dir>

Example:
  START_NS=0 END_NS=200 FRAME_INTERVAL_NS=1 \
  bash run_gmx_mmpbsa_general.sh \
    /path/system /path/full_clustered.xtc /path/mmpbsa_groups.ndx \
    RECEPTOR LIGAND /path/mmpbsa_charmm_pb_general.in /path/results

Optional environment variables:
  GMX_BIN=gmx
  MMPBSA_BIN=gmx_MMPBSA
  TPR_NAME=md_prod.tpr
  TOP_NAME=topol.top
  START_NS=0
  END_NS=200
  FRAME_INTERVAL_NS=1
  MPI_RANKS=1
  REFERENCE_PDB=/path/exactly_matched_complex_reference.pdb
USAGE
    exit 2
fi

SIM_DIR=$(cd "$1" && pwd)
TRAJ=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
NDX_FILE=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
RECEPTOR_NAME=$4
LIGAND_NAME=$5
MMPBSA_INPUT=$(cd "$(dirname "$6")" && pwd)/$(basename "$6")
OUT_DIR=$(mkdir -p "$7" && cd "$7" && pwd)

GMX_BIN=${GMX_BIN:-gmx}
MMPBSA_BIN=${MMPBSA_BIN:-gmx_MMPBSA}
TPR_NAME=${TPR_NAME:-md_prod.tpr}
TOP_NAME=${TOP_NAME:-topol.top}
START_NS=${START_NS:-0}
END_NS=${END_NS:-200}
FRAME_INTERVAL_NS=${FRAME_INTERVAL_NS:-1}
MPI_RANKS=${MPI_RANKS:-1}

TPR="$SIM_DIR/$TPR_NAME"
TOP="$SIM_DIR/$TOP_NAME"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

for required_file in "$TPR" "$TOP" "$TRAJ" "$NDX_FILE" "$MMPBSA_INPUT"; do
    [[ -f "$required_file" ]] || die "missing required file: $required_file"
done

command -v "$GMX_BIN" >/dev/null 2>&1 || die "cannot find GROMACS executable: $GMX_BIN"
command -v "$MMPBSA_BIN" >/dev/null 2>&1 || die "cannot find gmx_MMPBSA executable: $MMPBSA_BIN"

if [[ -z "${AMBERHOME:-}" && -n "${CONDA_PREFIX:-}" ]]; then
    export AMBERHOME="$CONDA_PREFIX"
fi
[[ -n "${AMBERHOME:-}" ]] || die "AMBERHOME is unset; activate the gmxMMPBSA environment"

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

RECEPTOR_GROUP=$(index_group_id "$NDX_FILE" "$RECEPTOR_NAME")
LIGAND_GROUP=$(index_group_id "$NDX_FILE" "$LIGAND_NAME")

if [[ -z "$RECEPTOR_GROUP" || -z "$LIGAND_GROUP" ]]; then
    echo "Available groups in $NDX_FILE:" >&2
    list_groups "$NDX_FILE" >&2
    die "could not resolve receptor='$RECEPTOR_NAME' and ligand='$LIGAND_NAME'"
fi

if [[ "$RECEPTOR_GROUP" == "$LIGAND_GROUP" ]]; then
    die "receptor and ligand resolve to the same index group"
fi

SAMPLED_TRAJ="$OUT_DIR/complex_${START_NS}-${END_NS}ns_every_${FRAME_INTERVAL_NS}ns.xtc"
if [[ ! -s "$SAMPLED_TRAJ" ]]; then
    echo "[1/2] Sampling ${START_NS}-${END_NS} ns every ${FRAME_INTERVAL_NS} ns"
    printf '0\n' | "$GMX_BIN" trjconv \
        -s "$TPR" \
        -f "$TRAJ" \
        -o "$SAMPLED_TRAJ" \
        -b "$START_NS" \
        -e "$END_NS" \
        -dt "$FRAME_INTERVAL_NS" \
        -tu ns
else
    echo "[1/2] Reusing sampled trajectory: $SAMPLED_TRAJ"
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
    -cg "$RECEPTOR_GROUP" "$LIGAND_GROUP"
    -cp "$TOP"
    -o "$RESULT_DAT"
    -eo "$RESULT_CSV"
)

if grep -Eq '^[[:space:]]*&decomp' "$MMPBSA_INPUT"; then
    CMD+=( -do "$DECOMP_DAT" -deo "$DECOMP_CSV" )
fi

if [[ -n "${REFERENCE_PDB:-}" ]]; then
    [[ -f "$REFERENCE_PDB" ]] || die "missing REFERENCE_PDB: $REFERENCE_PDB"
    CMD+=( -cr "$REFERENCE_PDB" )
fi

echo "[2/2] Running gmx_MMPBSA"
echo "      receptor: $RECEPTOR_NAME (group $RECEPTOR_GROUP)"
echo "      ligand:   $LIGAND_NAME (group $LIGAND_GROUP)"

if (( MPI_RANKS > 1 )); then
    command -v mpirun >/dev/null 2>&1 || die "MPI_RANKS=$MPI_RANKS but mpirun is unavailable"
    python -c 'from mpi4py import MPI; print(MPI.Get_library_version())' >/dev/null 2>&1 \
        || die "MPI_RANKS=$MPI_RANKS but mpi4py is unavailable"
    mpirun -np "$MPI_RANKS" "${CMD[@]}" 2>&1 | tee "$LOG_FILE"
else
    "${CMD[@]}" 2>&1 | tee "$LOG_FILE"
fi

cat <<EOF

Completed: $OUT_DIR
Check these files:
  FINAL_RESULTS_MMPBSA.dat
  FINAL_RESULTS_MMPBSA.csv
  gmx_MMPBSA.log

If &decomp was enabled, also check:
  FINAL_DECOMP_MMPBSA.dat
  FINAL_DECOMP_MMPBSA.csv

Before interpretation, inspect the generated fixed complex, receptor and ligand
structures and verify that the selected groups are correct.
EOF
