#!/usr/bin/env python3
"""Create matched SG-SG proximity and MM/PBSA time-course panels.

The left panel reports the instantaneous nearest *cross-trimer* SG-SG pair.
It is a geometric proximity readout only and does not establish formation of a
covalent disulfide bond. Distances are converted from nm to Angstrom.
"""

from __future__ import annotations

import csv
import os
from io import StringIO
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_ROOT = Path(os.environ.get("COLLAGEN_DATA_ROOT", REPO_ROOT / "data"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "geometry_energy_timecourse"

SG_SYSTEMS = ("mC", "CC_1", "NC")
ENERGY_SYSTEMS = ("P10", "mC", "CC_1", "NC")
LABELS = {"P10": "P10", "CC_1": "CC", "mC": "mC", "NC": "NC"}
COLORS = {"P10": "#C84E3B", "CC_1": "#7652A8", "mC": "#2879B8", "NC": "#238A56"}
# A shorter window preserves transient nearest-pair exchange events while
# retaining the thin-raw/thick-smoothed visual language of the energy panel.
SG_SMOOTHING_FRAMES = 21  # 2.1 ns for the retained 0.1-ns SG-SG sampling.

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.linewidth": 0.85,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def read_nearest_cross_trimer_distance(system: str) -> pd.DataFrame:
    """Return the nearest SG-SG distance at each frame, in Angstrom."""
    path = INPUT_ROOT / system / "disulfide_geometry" / "sg_sg_pair_distance_timeseries.csv"
    table = pd.read_csv(path)
    table = table[
        table["pair_class"].eq("between_trimers")
        & table["time_ns"].between(0.0, 200.0, inclusive="both")
    ].copy()
    nearest = table.loc[table.groupby("time_ns")["sg_sg_distance_nm"].idxmin()].sort_values("time_ns")
    nearest = nearest[["time_ns", "pair", "sg_sg_distance_nm"]].copy()
    nearest["sg_sg_distance_A"] = 10.0 * nearest.pop("sg_sg_distance_nm")
    nearest["rolling_median_sg_sg_distance_A"] = nearest["sg_sg_distance_A"].rolling(
        window=SG_SMOOTHING_FRAMES,
        center=True,
        min_periods=1,
    ).median()
    nearest.insert(0, "system", LABELS[system])
    return nearest


def read_delta_total(system: str) -> pd.DataFrame:
    """Read the Delta Energy Terms block from gmx_MMPBSA framewise output."""
    path = INPUT_ROOT / system / "mmpbsa_ABC_vs_DEF" / "FINAL_RESULTS_MMPBSA.csv"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next(index for index, line in enumerate(lines) if line.strip() == "Delta Energy Terms")
    records: list[str] = []
    for line in lines[start + 1 :]:
        if not line.strip():
            break
        records.append(line)
    table = pd.read_csv(StringIO("\n".join(records)))[["Frame #", "TOTAL"]].copy()
    table["time_ns"] = table["Frame #"].astype(float) - 1.0
    table = table[table["time_ns"].between(0.0, 200.0, inclusive="both")]
    table.rename(columns={"TOTAL": "delta_g_mmpbsa_kcal_mol"}, inplace=True)
    table.insert(0, "system", LABELS[system])
    return table[["system", "time_ns", "delta_g_mmpbsa_kcal_mol"]]


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.06, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    distance_data = [read_nearest_cross_trimer_distance(system) for system in SG_SYSTEMS]
    energy_data = [read_delta_total(system) for system in ENERGY_SYSTEMS]
    pd.concat(distance_data, ignore_index=True).to_csv(
        OUTPUT_DIR / "source_nearest_cross_trimer_sg_sg_distance_A.csv",
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        float_format="%.5f",
    )
    pd.concat(energy_data, ignore_index=True).to_csv(
        OUTPUT_DIR / "source_mmpbsa_total_timecourse.csv",
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        float_format="%.5f",
    )

    fig, (ax_distance, ax_energy) = plt.subplots(1, 2, figsize=(10.4, 3.45))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.80, bottom=0.25, wspace=0.34)

    for frame in distance_data:
        system_label = frame["system"].iat[0]
        system = next(key for key, value in LABELS.items() if value == system_label)
        ax_distance.plot(
            frame["time_ns"],
            frame["sg_sg_distance_A"],
            color=COLORS[system],
            lw=0.55,
            alpha=0.15,
        )
        ax_distance.plot(
            frame["time_ns"],
            frame["rolling_median_sg_sg_distance_A"],
            color=COLORS[system],
            lw=1.55,
            label=system_label,
        )
    ax_distance.axhline(4.0, color="#66717E", lw=0.85, ls="--", zorder=0)
    ax_distance.text(0.99, 0.95, "P10: no Cys", transform=ax_distance.transAxes, ha="right", va="top", fontsize=7, color=COLORS["P10"])
    ax_distance.text(0.99, 0.86, "Dashed: 4 Å threshold", transform=ax_distance.transAxes, ha="right", va="top", fontsize=6.8, color="#596574")
    ax_distance.set_title("Nearest inter-trimer sulfur proximity", loc="left", fontsize=10, fontweight="bold")
    ax_distance.set_ylabel("Nearest cross-trimer SG-SG distance (Å)")
    ax_distance.legend(loc="upper left", ncol=3, fontsize=7.2, handlelength=2.1, columnspacing=1.0)
    ax_distance.text(
        0.985,
        0.07,
        "Thin: 0.1-ns frames\nThick: 2.1-ns rolling median",
        transform=ax_distance.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        color="#596574",
    )
    add_panel_label(ax_distance, "a")

    for frame in energy_data:
        system_label = frame["system"].iat[0]
        system = next(key for key, value in LABELS.items() if value == system_label)
        rolling = frame["delta_g_mmpbsa_kcal_mol"].rolling(window=11, center=True, min_periods=1).mean()
        ax_energy.plot(frame["time_ns"], frame["delta_g_mmpbsa_kcal_mol"], color=COLORS[system], alpha=0.16, lw=0.65)
        ax_energy.plot(frame["time_ns"], rolling, color=COLORS[system], lw=1.55, label=system_label)
    ax_energy.axhline(0, color="#9AA1A8", lw=0.75, ls="--", zorder=0)
    ax_energy.set_title("Time-resolved association estimate", loc="left", fontsize=10, fontweight="bold")
    ax_energy.set_ylabel(r"$\Delta G_{\mathrm{MM/PBSA}}$ (kcal mol$^{-1}$)")
    ax_energy.legend(loc="upper left", ncol=4, fontsize=7.2, handlelength=2.0, columnspacing=0.8)
    ax_energy.text(0.985, 0.07, "Thin: 1-ns frames\nThick: 11-ns rolling mean", transform=ax_energy.transAxes, ha="right", va="bottom", fontsize=6.8, color="#596574")
    add_panel_label(ax_energy, "b")

    for axis in (ax_distance, ax_energy):
        axis.set_xlim(0, 200)
        axis.set_xticks([0, 50, 100, 150, 200])
        axis.set_xlabel("Simulation time (ns)")
        axis.tick_params(labelsize=8)
    ax_distance.set_ylim(2.5, 12.0)
    ax_distance.set_yticks([4, 6, 8, 10, 12])

    fig.suptitle("Inter-trimer sulfur proximity and association energy", x=0.085, ha="left", fontsize=13, fontweight="bold")
    fig.text(
        0.085,
        0.07,
        "Panel a reports the nearest eligible SG-SG pair in each frame; pair identity may change. A distance below 4 Å denotes geometric proximity, not a formed covalent S-S bond. "
        "Panel b shows endpoint MM/PBSA estimates from one trajectory per system.",
        ha="left",
        fontsize=7,
        color="#596574",
    )
    save_figure(fig, OUTPUT_DIR / "Fig_intertrimer_SG_SG_and_MMPBSA_timecourses")
    plt.close(fig)


if __name__ == "__main__":
    main()
