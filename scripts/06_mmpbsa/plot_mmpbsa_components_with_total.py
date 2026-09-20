"""Add the reported mean TOTAL to the existing MM/PBSA component heat map."""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = Path(os.environ.get("COLLAGEN_MMPBSA_SOURCE", REPO_ROOT / "data" / "mmpbsa_analysis"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "mmpbsa_components_with_total"
ORDER = ["P10", "CC", "mC", "NC"]
FIELDS = ["VDWAALS", "EEL", "EPB", "NONPOLAR", "TOTAL"]
LABELS = [
    r"$\Delta E_{\mathrm{vdW}}$",
    r"$\Delta E_{\mathrm{ele}}$",
    r"$\Delta G_{\mathrm{PB}}$",
    r"$\Delta G_{\mathrm{nonpolar}}$",
    r"$\Delta G_{\mathrm{TOTAL}}$",
]


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(SOURCE_DIR / "mmpbsa_delta_per_frame.csv")
    archived = pd.read_csv(SOURCE_DIR / "mmpbsa_system_summary.csv").set_index("system")
    frame = frame[frame["system"].isin(ORDER)].copy()
    frame["NONPOLAR"] = frame["ENPOLAR"] + frame["EDISPER"]
    for system in ORDER:
        subset = frame.loc[frame["system"] == system].sort_values("time_ns")
        assert len(subset) == 201, (system, len(subset))
        np.testing.assert_allclose(subset["time_ns"].to_numpy(), np.arange(201))
        assert np.isfinite(subset[FIELDS].to_numpy()).all()
        np.testing.assert_allclose(
            subset["TOTAL"].mean(),
            archived.loc[system, "delta_total_mean_kcal_mol"],
            atol=1e-8,
        )
    # Preserve the original non-polar row; all dispersion entries are zero here.
    assert np.allclose(frame["EDISPER"], 0)
    means = frame.groupby("system")[FIELDS].mean().reindex(ORDER)
    matrix = means.to_numpy().T
    residual = frame["TOTAL"] - frame["GGAS"] - frame["GSOLV"]
    assert residual.abs().max() < 0.031

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "axes.unicode_minus": True,
    })
    cmap = LinearSegmentedColormap.from_list(
        "screening_matched_energy",
        [(0.00, "#2F7FA8"), (0.36, "#67B6B0"), (0.50, "#F2DE7D"),
         (0.74, "#EC9447"), (1.00, "#C84E3B")],
    )
    fig = plt.figure(figsize=(3.50, 3.35), facecolor="white")
    ax = fig.add_axes([0.245, 0.105, 0.535, 0.84])
    cax = fig.add_axes([0.835, 0.16, 0.025, 0.73])
    norm = TwoSlopeNorm(vmin=-450, vcenter=0, vmax=450)
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(4), ORDER)
    ax.set_yticks(range(5), LABELS)
    ax.tick_params(axis="x", length=3, width=0.6, pad=3)
    ax.tick_params(axis="y", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axhline(3.5, color="white", linewidth=1.6)
    ax.get_yticklabels()[-1].set_fontweight("bold")
    for row in range(5):
        for col in range(4):
            value = matrix[row, col]
            label = f"{value:.1f}" if row == 4 else f"{value:.0f}"
            color = "white" if abs(value) > 0.56 * 450 else "#20262D"
            ax.text(
                col, row, label, ha="center", va="center", color=color,
                fontsize=7.6, fontweight="bold" if row == 4 else "normal",
            )
    cb = fig.colorbar(im, cax=cax, ticks=[-400, -200, 0, 200, 400])
    cb.ax.tick_params(labelsize=7, length=2.5, pad=2)
    cb.outline.set_linewidth(0.5)
    cb.set_label(r"Energy (kcal mol$^{-1}$)", labelpad=5, fontsize=7.5)
    cb.ax.yaxis.set_label_position("left")

    stem = OUTPUT_DIR / "Fig_MMPBSA_components_with_total"
    fig.savefig(stem.with_suffix(".png"), dpi=600, transparent=True)
    fig.savefig(stem.with_suffix(".svg"), transparent=True)
    fig.savefig(stem.with_suffix(".pdf"), transparent=True)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, transparent=True,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUTPUT_DIR / "Fig_MMPBSA_components_with_total_preview.png", dpi=300,
                facecolor="white", transparent=False)

    rows = []
    for component in FIELDS:
        for system in ORDER:
            value = means.loc[system, component]
            rows.append({
                "system": system, "component": component,
                "mean_kcal_mol": value,
                "display_value": f"{value:.1f}" if component == "TOTAL" else f"{value:.0f}",
                "n_frames": 201, "start_ns": 0, "end_ns": 200,
            })
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "heatmap_source_data.csv", index=False)
    qa = {
        "conclusion": "Reported endpoint TOTAL is displayed alongside component means; no new mechanistic claim.",
        "archetype": "quantitative grid",
        "backend": "Python matplotlib",
        "units": "kcal/mol",
        "statistic": "Arithmetic mean of 201 exported frames, 0-200 ns inclusive.",
        "independent_trajectories_per_system": 1,
        "total_source": "Exported TOTAL field, not the sum of the four displayed components.",
        "max_total_minus_gas_minus_solv_kcal_mol": float(residual.abs().max()),
        "total_matches_archived_summary": True,
        "original_four_component_rows_unchanged": True,
        "color_limits": [-450, 450],
        "raw_data_modified": False,
        "caution": "Existing nonzero internal delta terms require topology-consistency review; results are not experimental binding free energies.",
    }
    (OUTPUT_DIR / "figure_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    (OUTPUT_DIR / "figure_legend.txt").write_text(
        "MM/PBSA association-energy components and total endpoint estimate. Values are "
        "arithmetic means over 201 frames at 1-ns intervals from 0 to 200 ns. All values "
        "are complex minus receptor minus ligand, in kcal/mol. The TOTAL row uses the "
        "reported TOTAL field and is not reconstructed from the four displayed components. "
        "Configurational entropy is not included. One trajectory was analyzed per system; "
        "frames are not independent simulation replicates. Existing nonzero internal delta "
        "terms require topology-consistency review before definitive interpretation.\n",
        encoding="utf-8",
    )
    plt.close(fig)
    print(means.to_string())
    print(stem)


if __name__ == "__main__":
    main()
