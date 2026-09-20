#!/usr/bin/env python3
"""Export the shared free-energy colourbar used by the RMSD-Rg landscapes."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "docking_MD_results"
STEM = "Fig_shared_free_energy_colorbar_transparent_300dpi"
DCCM_STEM = "Fig_shared_DCCM_colorbar_transparent_300dpi"
FREE_ENERGY_CAP = 12.5

SCREENING_CMAP = LinearSegmentedColormap.from_list(
    "screening_palette",
    ["#2C79AD", "#54A4AE", "#91C6AE", "#D4DA8B", "#F3C35D", "#E98449", "#C73931"],
    N=256,
)


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.2,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "none",
            "figure.facecolor": "none",
        }
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(1.20, 3.00))
    figure.patch.set_alpha(0)
    axis = figure.add_axes([0.08, 0.08, 0.23, 0.84])
    colorbar = mpl.colorbar.ColorbarBase(
        axis,
        cmap=SCREENING_CMAP,
        norm=Normalize(0.0, FREE_ENERGY_CAP),
        orientation="vertical",
    )
    colorbar.set_ticks([0.0, FREE_ENERGY_CAP / 2, FREE_ENERGY_CAP])
    colorbar.set_ticklabels(["Low\n0", "Moderate\n6.25", "High\n12.5"])
    colorbar.ax.tick_params(labelsize=6.2, length=2.7, pad=3.0)
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", fontsize=7.2, labelpad=8.0)

    for extension in ("svg", "pdf", "png"):
        figure.savefig(
            OUTPUT_DIR / f"{STEM}.{extension}",
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            pad_inches=0.015,
            transparent=True,
            facecolor="none",
        )
    figure.savefig(
        OUTPUT_DIR / f"{STEM}.tiff",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.015,
        transparent=True,
        facecolor="none",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)

    dccm_figure = plt.figure(figsize=(1.20, 3.00))
    dccm_figure.patch.set_alpha(0)
    dccm_axis = dccm_figure.add_axes([0.08, 0.08, 0.23, 0.84])
    dccm_colorbar = mpl.colorbar.ColorbarBase(
        dccm_axis,
        cmap=SCREENING_CMAP,
        norm=Normalize(-1.0, 1.0),
        orientation="vertical",
    )
    dccm_colorbar.set_ticks([-1.0, 0.0, 1.0])
    dccm_colorbar.set_ticklabels(["Anti-correlated\n-1", "Weakly correlated\n0", "Correlated\n+1"])
    dccm_colorbar.ax.tick_params(labelsize=6.2, length=2.7, pad=3.0)
    dccm_colorbar.set_label("Dynamic cross-correlation", fontsize=7.2, labelpad=8.0)
    for extension in ("svg", "pdf", "png"):
        dccm_figure.savefig(
            OUTPUT_DIR / f"{DCCM_STEM}.{extension}",
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            pad_inches=0.015,
            transparent=True,
            facecolor="none",
        )
    dccm_figure.savefig(
        OUTPUT_DIR / f"{DCCM_STEM}.tiff",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.015,
        transparent=True,
        facecolor="none",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(dccm_figure)
    print(f"Saved {STEM} and {DCCM_STEM} to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
