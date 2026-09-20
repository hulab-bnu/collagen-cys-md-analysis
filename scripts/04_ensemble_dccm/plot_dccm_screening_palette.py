#!/usr/bin/env python3
"""Replot the triple-helix DCCM matrices with the screening heatmap palette.

The numerical DCCM remains unchanged: -1 is anti-correlated, 0 is weakly
correlated, and +1 is correlated. Only the colour language is harmonised with
the position-screening heatmap.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SOURCE_DIR = Path(os.environ.get("COLLAGEN_DCCM_SOURCE", REPO_ROOT / "data" / "ensemble_dccm"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "ensemble_dccm"

SYSTEMS = (
    ("CC", "CC", "#236FAE"),
    ("mC ok", "mC ok", "#267BB7"),
    ("NC", "NC", "#25875C"),
    ("P10_0", "P10_0", "#B25945"),
)

# Same low-to-high language used by the site-screening heatmap.
SCREENING_CMAP = LinearSegmentedColormap.from_list(
    "screening_palette",
    ["#2C79AD", "#54A4AE", "#91C6AE", "#D4DA8B", "#F3C35D", "#E98449", "#C73931"],
    N=256,
)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.0,
            "axes.labelsize": 7.0,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 6.0,
            "ytick.labelsize": 6.0,
            "axes.linewidth": 0.72,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def load_matrices() -> tuple[dict[str, np.ndarray], dict[str, dict[str, object]]]:
    matrix_table = pd.read_csv(SOURCE_DIR / "source_data" / "dccm_matrices.csv")
    manifest = json.loads((SOURCE_DIR / "analysis_manifest.json").read_text(encoding="utf-8"))
    metadata = {row["system"]: row for row in manifest["summary_rows"]}
    matrices: dict[str, np.ndarray] = {}
    for system, _, _ in SYSTEMS:
        subset = matrix_table.loc[matrix_table["system"] == system]
        size = int(np.sqrt(len(subset)))
        if size * size != len(subset):
            raise ValueError(f"DCCM source data for {system} is not square")
        matrices[system] = subset["DCCM"].to_numpy(dtype=float).reshape(size, size)
    return matrices, metadata


def cys_positions(metadata: dict[str, object]) -> list[int]:
    value = str(metadata["Cys_positions_per_chain"])
    if value.lower() == "none":
        return []
    return [int(item) for item in value.split(",")]


def draw_matrix(ax: plt.Axes, matrix: np.ndarray, metadata: dict[str, object], label: str, color: str, panel: str) -> mpl.image.AxesImage:
    residues_per_chain = int(metadata["segment_residues_per_chain"])
    image = ax.imshow(matrix, cmap=SCREENING_CMAP, norm=Normalize(-1, 1), interpolation="nearest", origin="upper", rasterized=True)
    for boundary in (residues_per_chain - 0.5, 2 * residues_per_chain - 0.5):
        ax.axvline(boundary, color="#465A64", linewidth=0.68, linestyle="--", zorder=3)
        ax.axhline(boundary, color="#465A64", linewidth=0.68, linestyle="--", zorder=3)
    for position in cys_positions(metadata):
        for chain_index in range(3):
            coordinate = chain_index * residues_per_chain + position - 0.5
            ax.axvline(coordinate, color="#F3AE22", linewidth=0.65, alpha=0.98, zorder=4)
            ax.axhline(coordinate, color="#F3AE22", linewidth=0.65, alpha=0.98, zorder=4)
    centres = [residues_per_chain / 2 - 0.5, 1.5 * residues_per_chain - 0.5, 2.5 * residues_per_chain - 0.5]
    ax.set_xticks(centres, ["A", "B", "C"])
    ax.set_yticks(centres, ["A", "B", "C"])
    ax.set_xlabel("Chain-residue block")
    ax.set_ylabel("Chain-residue block")
    cys_text = "no Cys" if not cys_positions(metadata) else "Cys " + ", ".join(map(str, cys_positions(metadata)))
    ax.set_title(f"{label} | {residues_per_chain} aa/chain; {cys_text}", loc="left", color=color, fontweight="bold", pad=3)
    ax.text(-0.16, 1.04, panel, transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="bottom")
    ax.tick_params(length=2.5, pad=1.5)
    return image


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(OUTPUT_DIR / f"{stem}.{extension}", dpi=500 if extension == "png" else None, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(
        OUTPUT_DIR / f"{stem}.tiff",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.04,
        pil_kwargs={"compression": "tiff_lzw"},
    )


def main() -> None:
    configure_style()
    matrices, metadata = load_matrices()
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 6.0))
    fig.subplots_adjust(left=0.095, right=0.845, top=0.865, bottom=0.105, hspace=0.34, wspace=0.30)
    for panel, ax, (key, label, color) in zip("abcd", axes.flat, SYSTEMS):
        draw_matrix(ax, matrices[key], metadata[key], label, color, panel)

    colorbar_axis = fig.add_axes([0.875, 0.245, 0.022, 0.52])
    colorbar = mpl.colorbar.ColorbarBase(colorbar_axis, cmap=SCREENING_CMAP, norm=Normalize(-1, 1), orientation="vertical")
    colorbar.set_ticks([-1, 0, 1])
    colorbar.set_ticklabels(["Anti-correlated\n-1", "Weakly correlated\n0", "Correlated\n+1"])
    colorbar.ax.tick_params(labelsize=6.0, length=2.6, pad=3)
    colorbar.set_label("Dynamic cross-correlation", fontsize=7.0, labelpad=8)

    fig.suptitle("Panel B | Dynamic cross-correlation matrices", x=0.45, y=0.972, fontsize=11.0, fontweight="bold")
    fig.text(0.45, 0.94, "Screening-palette colours are used only to harmonise visual language; DCCM values remain signed from -1 to +1.", ha="center", fontsize=6.15, color="#465B67")
    fig.text(0.45, 0.020, "Dashed lines delimit chains A, B and C. Gold guides mark every Cys residue in each chain.", ha="center", fontsize=6.0, color="#465B67")
    save_figure(fig, "Fig_DCCM_screening_palette")
    plt.close(fig)
    print(f"Saved DCCM figure to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
