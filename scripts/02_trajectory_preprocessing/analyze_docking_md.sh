#!/usr/bin/env bash
# Quantitative docking-MD analysis for six-chain collagen-like complexes.
# Input trajectories must be produced by repair_docking_md.sh:
#   full_clustered.xtc: protein + water + ions, for interface/hydration
#   protein_complex_fit.xtc: protein-only and globally fitted, for PCA/RMSD
#
# Usage:
#   bash analyze_docking_md.sh /path/to/simulation_dir /path/to/analysis_repaired_v2/CC_1
#
# Optional environment variables:
#   GMX_BIN=gmx_mpi BEGIN_PS=20000 STRIDE_PS=100 bash analyze_docking_md.sh ...
#   ABC_LAST=104 DEF_FIRST=105 DEF_LAST=206 bash analyze_docking_md.sh ...

set -euo pipefail

GMX_BIN="${GMX_BIN:-gmx}"
SIM_DIR="${1:?Usage: bash analyze_docking_md.sh SIMULATION_DIR REPAIRED_DIR}"
REPAIRED_DIR="${2:?Usage: bash analyze_docking_md.sh SIMULATION_DIR REPAIRED_DIR}"
SYSTEM_NAME="$(basename "${SIM_DIR}")"
BEGIN_PS="${BEGIN_PS:-20000}"
STRIDE_PS="${STRIDE_PS:-100}"
SNAPSHOT_STEP_PS="${SNAPSHOT_STEP_PS:-10000}"
BEGIN_NS="$(awk -v value="$BEGIN_PS" 'BEGIN {printf "%.6f", value / 1000}')"
STRIDE_NS="$(awk -v value="$STRIDE_PS" 'BEGIN {printf "%.6f", value / 1000}')"

TPR="${SIM_DIR}/md_prod.tpr"
FULL_TRAJ="${REPAIRED_DIR}/full_clustered.xtc"
PROTEIN_TRAJ="${REPAIRED_DIR}/protein_complex_fit.xtc"
REF="${REPAIRED_DIR}/reference_protein.gro"
OUT="${REPAIRED_DIR}/docking_metrics"
NDX="${OUT}/docking_groups.ndx"

for file in "$TPR" "$FULL_TRAJ" "$PROTEIN_TRAJ" "$REF"; do
    [[ -f "$file" ]] || { echo "ERROR: Missing input: $file" >&2; exit 1; }
done
command -v "$GMX_BIN" >/dev/null 2>&1 || {
    echo "ERROR: GROMACS executable not found: $GMX_BIN" >&2
    exit 1
}

# Residue indices are zero-based and unique across the six-chain topology.
# Chains A-C form trimer ABC; chains D-F form trimer DEF.
case "$SYSTEM_NAME" in
    CC_1) : "${ABC_LAST:=104}"; : "${DEF_FIRST:=105}"; : "${DEF_LAST:=206}" ;;
    NC)   : "${ABC_LAST:=98}";  : "${DEF_FIRST:=99}";  : "${DEF_LAST:=194}" ;;
    mC|mC_ok)
          : "${ABC_LAST:=89}";  : "${DEF_FIRST:=90}";  : "${DEF_LAST:=179}" ;;
    P10|P10_0)
          # Validated split for the 180-residue P10:P10 topology.
          : "${ABC_LAST:=89}";  : "${DEF_FIRST:=90}";  : "${DEF_LAST:=179}" ;;
    *)
          : "${ABC_LAST:?Set ABC_LAST for this construct}";
          : "${DEF_FIRST:?Set DEF_FIRST for this construct}";
          : "${DEF_LAST:?Set DEF_LAST for this construct}" ;;
esac

mkdir -p "$OUT"

HAS_CYS=0
if grep -Eiq '(^|[[:space:]])(CYS|CYX|CYM)([[:space:]]|$)' "$SIM_DIR"/topol_Protein_chain_*.itp; then
    HAS_CYS=1
fi

echo "[1/10] Building static protein, trimer, backbone, C-alpha, water, and optional Cys groups..."
SELECTIONS=(
    '"COMPLEX_Protein" group "Protein"'
    '"COMPLEX_Backbone" group "Backbone"'
    '"COMPLEX_CA" group "C-alpha"'
    "\"TRIMER_ABC\" group \"Protein\" and resindex 0 to ${ABC_LAST}"
    "\"TRIMER_DEF\" group \"Protein\" and resindex ${DEF_FIRST} to ${DEF_LAST}"
    "\"ABC_Backbone\" group \"Backbone\" and resindex 0 to ${ABC_LAST}"
    "\"DEF_Backbone\" group \"Backbone\" and resindex ${DEF_FIRST} to ${DEF_LAST}"
    "\"ABC_CA\" group \"C-alpha\" and resindex 0 to ${ABC_LAST}"
    "\"DEF_CA\" group \"C-alpha\" and resindex ${DEF_FIRST} to ${DEF_LAST}"
    '"WATER" resname SOL'
    '"WATER_O" resname SOL and name OW OH2 O'
)
if (( HAS_CYS )); then
    SELECTIONS+=(
        "\"CYS_ABC\" group \"Protein\" and resname CYS CYX CYM and resindex 0 to ${ABC_LAST}"
        "\"CYS_DEF\" group \"Protein\" and resname CYS CYX CYM and resindex ${DEF_FIRST} to ${DEF_LAST}"
    )
else
    echo "No cysteine residues detected: Cys-specific interface counts will be recorded as zero."
fi
"$GMX_BIN" select -s "$TPR" -on "$NDX" -select "${SELECTIONS[@]}"

echo "[2/10] Global and trimer-internal backbone RMSD..."
printf "COMPLEX_Backbone\nCOMPLEX_Backbone\n" | "$GMX_BIN" rms \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -o "$OUT/rmsd_complex_backbone.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -mw no -pbc no -xvg none
printf "ABC_Backbone\nABC_Backbone\n" | "$GMX_BIN" rms \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -o "$OUT/rmsd_abc_internal.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -mw no -pbc no -xvg none
printf "DEF_Backbone\nDEF_Backbone\n" | "$GMX_BIN" rms \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -o "$OUT/rmsd_def_internal.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -mw no -pbc no -xvg none

echo "[3/10] Per-residue C-alpha RMSF..."
printf "COMPLEX_CA\n" | "$GMX_BIN" rmsf \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -o "$OUT/rmsf_ca_per_residue.xvg" -res -fit -b "$BEGIN_PS" -dt "$STRIDE_PS" -xvg none

echo "[4/10] Complex and individual-trimer radius of gyration..."
"$GMX_BIN" gyrate -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -sel 'group "COMPLEX_Protein"' \
    -o "$OUT/rg_complex.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" gyrate -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -sel 'group "TRIMER_ABC"' \
    -o "$OUT/rg_abc.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" gyrate -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -sel 'group "TRIMER_DEF"' \
    -o "$OUT/rg_def.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none

echo "[5/10] Inter-trimer COM separation and interface contact-residue counts..."
"$GMX_BIN" distance -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" \
    -select 'com of group "TRIMER_ABC" plus com of group "TRIMER_DEF"' \
    -oall "$OUT/abc_def_com_distance.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" select -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" -seltype res_com \
    -select \
    '"ABC_contact_residues" group "TRIMER_ABC" and same residue as within 0.45 of group "TRIMER_DEF"' \
    '"DEF_contact_residues" group "TRIMER_DEF" and same residue as within 0.45 of group "TRIMER_ABC"' \
    -os "$OUT/interface_contact_residue_counts.xvg" \
    -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none

echo "[6/10] Interface and protein-water hydrogen-bond time series..."
"$GMX_BIN" hbond -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" \
    -r 'group "TRIMER_ABC"' -t 'group "TRIMER_DEF"' \
    -num "$OUT/hbonds_abc_def.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" hbond -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" \
    -r 'group "TRIMER_ABC"' -t 'group "WATER"' \
    -num "$OUT/hbonds_abc_water.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" hbond -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" \
    -r 'group "TRIMER_DEF"' -t 'group "WATER"' \
    -num "$OUT/hbonds_def_water.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none

echo "[7/10] Interfacial bridging-water proxy and optional Cys-interface occupancy..."
"$GMX_BIN" select -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" -seltype res_com \
    -select '"bridging_water_proxy" same residue as (group "WATER_O" and within 0.35 of group "TRIMER_ABC") and same residue as (group "WATER_O" and within 0.35 of group "TRIMER_DEF")' \
    -os "$OUT/bridging_water_proxy_count.xvg" \
    -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
if (( HAS_CYS )); then
    "$GMX_BIN" select -s "$TPR" -f "$FULL_TRAJ" -n "$NDX" -seltype res_com \
        -select \
        '"ABC_Cys_interface" group "CYS_ABC" and same residue as within 0.45 of group "TRIMER_DEF"' \
        '"DEF_Cys_interface" group "CYS_DEF" and same residue as within 0.45 of group "TRIMER_ABC"' \
        -os "$OUT/cys_interface_residue_counts.xvg" \
        -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
else
    awk '!/^[[:space:]]*[@#]/ && NF {printf "%s 0 0\n", $1}' \
        "$OUT/interface_contact_residue_counts.xvg" \
        > "$OUT/cys_interface_residue_counts.xvg"
fi

echo "[8/10] SASA and buried interface area..."
"$GMX_BIN" sasa -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -surface 'group "TRIMER_ABC"' -output 'group "TRIMER_ABC"' \
    -o "$OUT/sasa_abc_isolated.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" sasa -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -surface 'group "TRIMER_DEF"' -output 'group "TRIMER_DEF"' \
    -o "$OUT/sasa_def_isolated.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
"$GMX_BIN" sasa -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -fgroup 'group "COMPLEX_Protein"' -surface 'group "COMPLEX_Protein"' -output 'group "COMPLEX_Protein"' \
    -o "$OUT/sasa_complex.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none

strip_xvg() { awk '!/^[[:space:]]*[@#]/ && NF {print}' "$1"; }
strip_xvg "$OUT/sasa_abc_isolated.xvg" > "$OUT/.sasa_abc.dat"
strip_xvg "$OUT/sasa_def_isolated.xvg" > "$OUT/.sasa_def.dat"
strip_xvg "$OUT/sasa_complex.xvg" > "$OUT/.sasa_complex.dat"
paste "$OUT/.sasa_abc.dat" "$OUT/.sasa_def.dat" "$OUT/.sasa_complex.dat" | \
    awk 'BEGIN {OFS=","; print "time_ns,buried_sasa_nm2"} {print $1,$2+$4-$6}' \
    > "$OUT/interface_buried_sasa_nm2.csv"
rm -f "$OUT/.sasa_abc.dat" "$OUT/.sasa_def.dat" "$OUT/.sasa_complex.dat"

echo "[9/10] PCA input/projection and C-alpha coordinates for DCCM..."
printf "COMPLEX_CA\nCOMPLEX_CA\n" | "$GMX_BIN" covar \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -o "$OUT/pca_eigenvalues.xvg" -v "$OUT/pca_eigenvectors.trr" \
    -av "$OUT/pca_average_ca.gro" -xpm "$OUT/pca_covariance.xpm" \
    -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
# gmx anaeig requests the covar fitting group and then the eigenvector group.
printf "COMPLEX_CA\nCOMPLEX_CA\n" | "$GMX_BIN" anaeig \
    -s "$REF" -f "$PROTEIN_TRAJ" -n "$NDX" -v "$OUT/pca_eigenvectors.trr" \
    -2d "$OUT/pca_pc1_pc2.xvg" -first 1 -last 2 \
    -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none
printf "COMPLEX_CA\n" | "$GMX_BIN" traj \
    -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
    -ox "$OUT/ca_coordinates.xvg" -b "$BEGIN_NS" -dt "$STRIDE_NS" -tu ns -xvg none

echo "[10/10] Representative protein snapshots and download archive..."
for time_ps in $(seq "$BEGIN_PS" "$SNAPSHOT_STEP_PS" 200000); do
    printf "COMPLEX_Protein\n" | "$GMX_BIN" trjconv \
        -s "$TPR" -f "$PROTEIN_TRAJ" -n "$NDX" \
        -o "$OUT/snapshot_${time_ps}ps.pdb" -dump "$time_ps"
done

cat > "$OUT/README.txt" <<EOF
System: ${SYSTEM_NAME}
Sampled interval: ${BEGIN_PS} to 200000 ps, every ${STRIDE_PS} ps
ABC (chains A-C): residue indices 0-${ABC_LAST}
DEF (chains D-F): residue indices ${DEF_FIRST}-${DEF_LAST}

Core readouts
- rmsd_*: global-complex and individual-trimer backbone RMSD (nm)
- rmsf_ca_per_residue: C-alpha RMSF (nm)
- rg_complex, rg_abc, rg_def: complex and individual-trimer radii of gyration (nm)
- abc_def_com_distance: inter-trimer COM separation (nm)
- interface_contact_residue_counts: residues within 0.45 nm across the interface
- hbonds_abc_def: direct inter-trimer hydrogen bonds
- hbonds_*_water: hydration hydrogen bonds for each trimer
- bridging_water_proxy_count: water oxygen within 0.35 nm of both trimers
- cys_interface_residue_counts: Cys residues within 0.45 nm of the opposite trimer (all zero for Cys-free controls)
- interface_buried_sasa_nm2.csv: SASA_ABC + SASA_DEF - SASA_complex
- pca_pc1_pc2 and pca_eigenvalues: PCA landscape inputs
- ca_coordinates: sampled C-alpha coordinates for DCCM calculation
- snapshot_*ps.pdb: protein-only conformational ensemble, every ${SNAPSHOT_STEP_PS} ps
EOF

tar -czf "${REPAIRED_DIR}/${SYSTEM_NAME}_docking_metrics.tar.gz" -C "$OUT" .
echo "Finished: ${REPAIRED_DIR}/${SYSTEM_NAME}_docking_metrics.tar.gz"
