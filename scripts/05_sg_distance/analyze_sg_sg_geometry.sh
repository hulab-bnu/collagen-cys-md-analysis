#!/usr/bin/env bash
# Extract all Cys SG coordinates from a repaired docking trajectory and run an
# unbiased SG-SG proximity analysis. Run this script on the GROMACS server.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    cat <<'USAGE'
Usage:
  bash analyze_sg_sg_geometry.sh <simulation_dir> <repaired_dir> [output_dir]

Example:
  bash analyze_sg_sg_geometry.sh \
    /path/to/simulation/CC_1 \
    /path/to/repaired/CC_1

Environment overrides:
  GMX_BIN=/path/to/gmx  PYTHON_BIN=python3  START_NS=0  END_NS=200  STRIDE_NS=0.1
USAGE
    exit 2
fi

SIM_DIR=$(cd "$1" && pwd)
REPAIRED_DIR=$(cd "$2" && pwd)
OUT_DIR=${3:-"$REPAIRED_DIR/disulfide_geometry"}
OUT_DIR=$(mkdir -p "$OUT_DIR" && cd "$OUT_DIR" && pwd)

GMX_BIN=${GMX_BIN:-gmx}
PYTHON_BIN=${PYTHON_BIN:-python3}
START_NS=${START_NS:-0}
END_NS=${END_NS:-200}
STRIDE_NS=${STRIDE_NS:-0.1}

TPR="$SIM_DIR/md_prod.tpr"
TRAJ="$REPAIRED_DIR/full_clustered.xtc"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_PY="$SCRIPT_DIR/sg_sg_geometry_from_xvg.py"

if [[ ! -f "$TPR" ]]; then
    echo "ERROR: missing run-input file: $TPR" >&2
    exit 1
fi
if [[ ! -f "$TRAJ" ]]; then
    echo "ERROR: missing repaired full-atom trajectory: $TRAJ" >&2
    echo "Run the repair script first; SG coordinates cannot be extracted from a CA-only trajectory." >&2
    exit 1
fi
if [[ ! -f "$ANALYSIS_PY" ]]; then
    echo "ERROR: missing companion Python script: $ANALYSIS_PY" >&2
    exit 1
fi

# Older analysis runs may keep snapshot PDBs inside docking_metrics/, while
# newer runs may not retain them at all. A single protein PDB is needed only
# to map the GROMACS SG atom indices to chain/residue labels.
REFERENCE_PDB=$(find "$REPAIRED_DIR" -maxdepth 2 -type f -name 'snapshot_*ps.pdb' -print | sort -V | head -n 1 || true)
if [[ -z "$REFERENCE_PDB" ]]; then
    REFERENCE_PDB="$OUT_DIR/sg_reference_protein.pdb"
    echo "No saved snapshot PDB found; generating a protein reference at ${START_NS} ns"
    printf 'Protein\n' | "$GMX_BIN" trjconv \
        -s "$TPR" \
        -f "$TRAJ" \
        -o "$REFERENCE_PDB" \
        -b "$START_NS" \
        -e "$START_NS" \
        -tu ns
fi
if [[ ! -s "$REFERENCE_PDB" ]]; then
    echo "ERROR: could not create a protein reference PDB for SG atom mapping." >&2
    exit 1
fi

echo "[1/3] Selecting all cysteine sulfur atoms (SG)"
"$GMX_BIN" select \
    -s "$TPR" \
    -select 'name SG and (resname CYS or resname CYX or resname CYM)' \
    -seltype atom \
    -on "$OUT_DIR/cys_sg.ndx"

if ! grep -q '[0-9]' "$OUT_DIR/cys_sg.ndx"; then
    echo "ERROR: no SG atoms were selected. Check Cys residue/atom names in the topology." >&2
    exit 1
fi

if "$GMX_BIN" trajectory -h >/dev/null 2>&1; then
    TRAJECTORY_COMMAND=trajectory
else
    TRAJECTORY_COMMAND=traj
fi

echo "[2/3] Extracting SG coordinates every ${STRIDE_NS} ns from ${START_NS} to ${END_NS} ns"
# cys_sg.ndx contains exactly one group; therefore its index is 0.
printf '0\n' | "$GMX_BIN" "$TRAJECTORY_COMMAND" \
    -s "$TPR" \
    -f "$TRAJ" \
    -n "$OUT_DIR/cys_sg.ndx" \
    -ox "$OUT_DIR/cys_sg_coordinates.xvg" \
    -b "$START_NS" \
    -e "$END_NS" \
    -dt "$STRIDE_NS" \
    -tu ns \
    -xvg none

echo "[3/3] Calculating every SG-SG distance and potential disulfide-contact occupancy"
"$PYTHON_BIN" "$ANALYSIS_PY" \
    --coordinates "$OUT_DIR/cys_sg_coordinates.xvg" \
    --index "$OUT_DIR/cys_sg.ndx" \
    --reference-pdb "$REFERENCE_PDB" \
    --output-dir "$OUT_DIR" \
    --cutoff-nm 0.40

cat <<EOF

Completed: $OUT_DIR
Key files:
  sg_atom_mapping.csv                         SG identity and chain mapping
  all_sg_sg_pair_summary.csv                  median/IQR/minimum/occupancy for every pair
  cross_trimer_sg_sg_pair_summary.csv         ABC-versus-DEF pairs, ranked by median distance
  sg_sg_pair_distance_timeseries.csv          SG-SG distance at every sampled time
  potential_disulfide_contact_count.csv       geometric contact-pair counts over time

Interpretation: <0.40 nm denotes an SG-SG proximity event, not a formed S-S bond.
EOF
