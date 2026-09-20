#!/usr/bin/env python3
"""Recolour the four RMSD-Rg free-energy basins with the screening palette."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SOURCE_DIR = Path(os.environ.get("COLLAGEN_FEL_SOURCE", REPO_ROOT / "data" / "rmsd_rg_fel"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "rmsd_rg_fel"

SYSTEMS = (
    ("CC_shared_grid.csv", "CC", "#236FAE"),
    ("mC_shared_grid.csv", "mC", "#267BB7"),
    ("NC_shared_grid.csv", "NC", "#25875C"),
    # GPP10 is the experimental control name for the P10_0 trajectory.
    ("P10_shared_grid.csv", "GPP10", "#B25945"),
)
FREE_ENERGY_CAP = 12.5

# Match the site-screening and recoloured DCCM figures: low -> high.
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
            "axes.titlesize": 8.2,
            "xtick.labelsize": 6.1,
            "ytick.labelsize": 6.1,
            "axes.linewidth": 0.72,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def read_landscape(filename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[float, float]]:
    table = pd.read_csv(SOURCE_DIR / filename)
    x = np.sort(table["core_backbone_rmsd_nm"].unique())
    y = np.sort(table["core_radius_of_gyration_nm"].unique())
    energy = table.pivot(index="core_backbone_rmsd_nm", columns="core_radius_of_gyration_nm", values="delta_G_kJ_mol").reindex(index=x, columns=y).to_numpy(dtype=float)
    supported = table.pivot(index="core_backbone_rmsd_nm", columns="core_radius_of_gyration_nm", values="display_supported").reindex(index=x, columns=y).to_numpy(dtype=bool)
    x_grid, y_grid = np.meshgrid(x, y, indexing="ij")
    valid = np.isfinite(energy) & supported
    peak_index = np.unravel_index(np.nanargmin(np.where(valid, energy, np.nan)), energy.shape)
    peak = (float(x[peak_index[0]]), float(y[peak_index[1]]))
    return x_grid, y_grid, energy, valid, peak


def draw_landscape(ax: plt.Axes, filename: str, label: str, label_color: str, panel: str, show_x: bool, show_y: bool) -> None:
    x_grid, y_grid, energy, supported, peak = read_landscape(filename)
    displayed = np.ma.masked_where(~supported, np.minimum(energy, FREE_ENERGY_CAP))
    levels = np.linspace(0, FREE_ENERGY_CAP, 11)
    ax.contourf(x_grid, y_grid, displayed, levels=levels, cmap=SCREENING_CMAP, norm=Normalize(0, FREE_ENERGY_CAP), extend="max")
    ax.contour(x_grid, y_grid, displayed, levels=[2.5, 5.0, 7.5, 10.0], colors="#435861", linewidths=0.34, alpha=0.62)
    ax.scatter(*peak, marker="*", s=32, color="white", edgecolor="#20313A", linewidth=0.5, zorder=5)
    ax.set_facecolor("#F2F7F8")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.16, 1.04, panel, transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="bottom")
    ax.set_title(label, loc="left", color=label_color, fontweight="bold", pad=3)
    if show_x:
        ax.set_xlabel("RMSD (nm)")
    else:
        ax.tick_params(labelbottom=False)
    if show_y:
        ax.set_ylabel(r"$R_g$ (nm)")
    else:
        ax.tick_params(labelleft=False)


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(OUTPUT_DIR / f"{stem}.{extension}", dpi=500 if extension == "png" else None, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTPUT_DIR / f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.04, pil_kwargs={"compression": "tiff_lzw"})


def main() -> None:
    configure_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.20, 5.85))
    fig.subplots_adjust(left=0.10, right=0.82, top=0.86, bottom=0.12, hspace=0.32, wspace=0.27)
    for panel, ax, (filename, label, color) in zip("abcd", axes.flat, SYSTEMS):
        draw_landscape(ax, filename, label, color, panel, show_x=ax in axes[1], show_y=ax in axes[:, 0])

    colorbar_axis = fig.add_axes([0.865, 0.235, 0.023, 0.54])
    colorbar = mpl.colorbar.ColorbarBase(colorbar_axis, cmap=SCREENING_CMAP, norm=Normalize(0, FREE_ENERGY_CAP), orientation="vertical")
    colorbar.set_ticks([0, FREE_ENERGY_CAP / 2, FREE_ENERGY_CAP])
    colorbar.set_ticklabels(["Low\n0", "Moderate\n6.25", "High\n12.5"])
    colorbar.ax.tick_params(labelsize=6.0, length=2.6, pad=3)
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", fontsize=7.0, labelpad=7)

    fig.suptitle(r"Panel C | RMSD-$R_g$ energy basins", x=0.42, y=0.985, fontsize=11.0, fontweight="bold")
    fig.text(0.42, 0.936, "All coordinates use the common 30-residue segment per chain; GPP10 denotes the P10_0 control trajectory.", ha="center", fontsize=6.15, color="#465B67")
    fig.text(0.42, 0.025, "White stars mark the dominant sampled basin; uncoloured regions were not sufficiently sampled for free-energy estimation.", ha="center", fontsize=6.0, color="#465B67")
    save_figure(fig, "Fig_RMSD_Rg_energy_basins_screening_palette")
    plt.close(fig)
    print(f"Saved RMSD-Rg basins to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
