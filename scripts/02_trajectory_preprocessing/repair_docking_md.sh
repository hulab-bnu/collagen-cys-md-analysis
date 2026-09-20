#!/usr/bin/env bash
# Repair a GROMACS docking MD trajectory without modifying the raw simulation.
# The system contains two collagen-like trimers (six chains). PBC operations
# always treat all protein chains as one complex; no individual chain is fitted.
# Usage:
#   bash repair_docking_md.sh /path/to/simulation_dir /path/to/analysis_repaired
# Optional: GMX_BIN=gmx_mpi bash repair_docking_md.sh ...

set -euo pipefail

GMX_BIN="${GMX_BIN:-gmx}"
SRC_DIR="${1:-$PWD}"
OUT_DIR="${2:-${SRC_DIR}/analysis_repaired}"

TPR="${SRC_DIR}/md_prod.tpr"
XTC="${SRC_DIR}/md_prod.xtc"

for file in "$TPR" "$XTC"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file" >&2
        exit 1
    fi
done

if ! command -v "$GMX_BIN" >/dev/null 2>&1; then
    echo "ERROR: Cannot find GROMACS executable: $GMX_BIN" >&2
    echo "Set GMX_BIN, for example: GMX_BIN=gmx_mpi bash repair_docking_md.sh ..." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

if [[ -s "$OUT_DIR/whole_system.xtc" ]]; then
    echo "[1/5] Reusing existing whole_system.xtc..."
else
    echo "[1/5] Making molecules whole while retaining water and ions..."
    printf "System\n" | "$GMX_BIN" trjconv \
        -s "$TPR" \
        -f "$XTC" \
        -o "$OUT_DIR/whole_system.xtc" \
        -pbc mol
fi

# trjconv asks, in order: clustering group, centering group, output group.
# Protein -> Protein -> System keeps the complete solvated complex together.
echo "[2/5] Clustering and centering the complete system for solvent/interface analysis..."
printf "Protein\nProtein\nSystem\n" | "$GMX_BIN" trjconv \
    -s "$TPR" \
    -f "$OUT_DIR/whole_system.xtc" \
    -o "$OUT_DIR/full_clustered.xtc" \
    -pbc cluster \
    -center

echo "[3/5] Creating an all-protein fitted trajectory for RMSD, PCA, and visualization..."
printf "Backbone\nProtein\n" | "$GMX_BIN" trjconv \
    -s "$TPR" \
    -f "$OUT_DIR/full_clustered.xtc" \
    -o "$OUT_DIR/protein_complex_fit.xtc" \
    -fit rot+trans

echo "[4/5] Writing the time-zero protein reference structure..."
printf "Protein\n" | "$GMX_BIN" trjconv \
    -s "$TPR" \
    -f "$OUT_DIR/full_clustered.xtc" \
    -o "$OUT_DIR/reference_protein.gro" \
    -dump 0

echo "[5/5] Checking the fitted trajectory..."
"$GMX_BIN" check -f "$OUT_DIR/protein_complex_fit.xtc" \
    > "$OUT_DIR/gmx_check_protein_fit.txt" 2>&1

cat > "$OUT_DIR/README.txt" <<EOF
Generated from: $SRC_DIR

full_clustered.xtc
  Full system (protein, water, ions), PBC repaired and centered.
  Use for interface hydration, water bridges, solvent-mediated hydrogen bonds,
  and protein-water contacts.

protein_complex_fit.xtc
  Protein-only trajectory, fitted as one complete six-chain complex.
  Use for complex RMSD, PCA, RMSF, DCCM, conformational ensembles, and figures.

reference_protein.gro
  Protein reference at t = 0 ps after PBC repair.

Raw md_prod.xtc is not changed by this script.
EOF

echo
echo "Finished. Repaired files are in: $OUT_DIR"
