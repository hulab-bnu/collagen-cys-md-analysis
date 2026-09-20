#!/usr/bin/env python3
"""Create a publication-style figure for unbiased SG-SG proximity screening."""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


warnings.filterwarnings("ignore", message="Pandas requires version")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.3,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

SYSTEMS = ("CC_1", "mC", "NC")
DISPLAY = {"CC_1": "CC", "mC": "mC", "NC": "NC"}
# Match the construct palette used across the docking-MD figure set.
COLORS = {"CC_1": "#7250A5", "mC": "#2878B8", "NC": "#23864F"}
BLOCK_STARTS = np.arange(0.0, 200.0, 20.0)
WINDOW_START_NS = 0.0
WINDOW_END_NS = 200.0
SMOOTH_WINDOW_FRAMES = 21  # 2.1 ns at the retained 0.1-ns sampling interval


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.035, 1.07, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def block_medians(time: pd.Series, values: pd.Series) -> np.ndarray:
    block = pd.cut(
        time,
        bins=np.append(BLOCK_STARTS, 200.0001),
        labels=BLOCK_STARTS.astype(int),
        right=False,
        include_lowest=True,
    )
    return pd.DataFrame({"block": block, "value": values}).groupby("block", observed=True)["value"].median().to_numpy()


def load_system(root: Path, name: str) -> dict[str, object]:
    folder = root / name / "disulfide_geometry"
    time_series = pd.read_csv(folder / "sg_sg_pair_distance_timeseries.csv")
    cross_trimer = time_series[
        time_series["pair_class"].eq("between_trimers")
        & time_series["time_ns"].between(WINDOW_START_NS, WINDOW_END_NS, inclusive="both")
    ].copy()
    # Select the fixed candidate over the same full trajectory window shown in the figure.
    candidate_stats = (
        cross_trimer.groupby("pair", as_index=False)["sg_sg_distance_nm"]
        .agg(
            median_sg_sg_nm="median",
            q1_sg_sg_nm=lambda values: values.quantile(0.25),
            q3_sg_sg_nm=lambda values: values.quantile(0.75),
            minimum_sg_sg_nm="min",
        )
        .sort_values(["median_sg_sg_nm", "minimum_sg_sg_nm"])
    )
    top = candidate_stats.iloc[0]
    trajectory = cross_trimer[cross_trimer["pair"].eq(top["pair"])].sort_values("time_ns").copy()
    nearest_index = cross_trimer.groupby("time_ns")["sg_sg_distance_nm"].idxmin()
    instantaneous_minimum = cross_trimer.loc[nearest_index, ["time_ns", "pair", "sg_sg_distance_nm"]].sort_values("time_ns")
    count = pd.read_csv(folder / "potential_disulfide_contact_count.csv")
    count = count[count["time_ns"].between(WINDOW_START_NS, WINDOW_END_NS, inclusive="both")].copy()
    if not np.isclose(count["cutoff_nm"].iloc[0], 0.40):
        raise ValueError(f"{name}: expected a 0.40 nm contact cutoff")
    return {
        "name": name,
        "display": DISPLAY[name],
        "color": COLORS[name],
        "top": top,
        "trajectory": trajectory,
        "instantaneous_minimum": instantaneous_minimum,
        "top_occupancy": float((trajectory["sg_sg_distance_nm"] < 0.40).mean()),
        "any_occupancy": float((count["max_nonoverlapping_cross_trimer_pairs_lt_cutoff"] > 0).mean()),
        "mean_distinct_contacts": float(count["max_nonoverlapping_cross_trimer_pairs_lt_cutoff"].mean()),
    }


def write_source_data(systems: list[dict[str, object]], output_dir: Path) -> None:
    rows: list[pd.DataFrame] = []
    quant_rows: list[dict[str, object]] = []
    for system in systems:
        top = system["top"]
        trajectory = system["trajectory"][["time_ns", "sg_sg_distance_nm"]].copy()
        trajectory.insert(0, "pair", top["pair"])
        trajectory.insert(0, "system", system["display"])
        rows.append(trajectory)
        instantaneous = system["instantaneous_minimum"].copy()
        instantaneous.insert(0, "system", system["display"])
        instantaneous.to_csv(
            output_dir / f"SG_SG_instantaneous_minimum_{system['display']}.csv", index=False, float_format="%.5f"
        )
        quant_rows.append(
            {
                "system": system["display"],
                "closest_cross_trimer_pair": top["pair"],
                "median_sg_sg_nm": top["median_sg_sg_nm"],
                "q1_sg_sg_nm": top["q1_sg_sg_nm"],
                "q3_sg_sg_nm": top["q3_sg_sg_nm"],
                "minimum_sg_sg_nm": top["minimum_sg_sg_nm"],
                "top_pair_fraction_lt_0p40_nm": system["top_occupancy"],
                "any_cross_trimer_proximity_fraction": system["any_occupancy"],
                "mean_distinct_cross_trimer_contacts": system["mean_distinct_contacts"],
            }
        )
    pd.concat(rows, ignore_index=True).to_csv(output_dir / "SG_SG_top_candidate_source_data.csv", index=False, float_format="%.5f")
    pd.DataFrame(quant_rows).to_csv(output_dir / "SG_SG_top_candidate_summary.csv", index=False, float_format="%.5f")


def draw(systems: list[dict[str, object]], output_dir: Path, *, smoothed: bool, stem: str) -> None:
    fig, (ax_time, ax_occupancy) = plt.subplots(
        1,
        2,
        figsize=(6.85, 2.72),
        gridspec_kw={"width_ratios": [1.36, 0.72]},
    )
    fig.subplots_adjust(left=0.105, right=0.975, top=0.735, bottom=0.275, wspace=0.50)

    displayed_values: list[pd.Series] = []
    for system in systems:
        trajectory = system["instantaneous_minimum"]
        raw = trajectory["sg_sg_distance_nm"]
        values = raw.rolling(window=SMOOTH_WINDOW_FRAMES, center=True, min_periods=1).median() if smoothed else raw
        displayed_values.append(values)
        ax_time.plot(
            trajectory["time_ns"],
            values,
            color=system["color"],
            lw=1.35 if smoothed else 0.82,
            alpha=0.98,
            label=system["display"],
            zorder=2,
        )
    ax_time.axhline(0.40, color="#566474", lw=0.75, ls="--", zorder=0)
    ax_time.set_xlim(WINDOW_START_NS, WINDOW_END_NS)
    ymax = max(1.18, max(values.max() for values in displayed_values) + 0.05)
    ax_time.set_ylim(0.25, ymax)
    ax_time.set_xticks([0, 50, 100, 150, 200])
    ax_time.set_xlabel("Simulation time (ns)")
    ax_time.set_ylabel("Nearest cross-trimer SG-SG distance (nm)")
    ax_time.set_title("Rolling-median nearest Cys pair" if smoothed else "Per-frame nearest Cys pair", loc="left", fontsize=8.2, fontweight="bold", pad=8)
    ax_time.legend(loc="upper right", ncol=3, fontsize=5.9, handlelength=1.55, columnspacing=0.85, frameon=False, borderaxespad=0.20)
    panel_label(ax_time, "a")

    labels = [system["display"] for system in systems]
    fixed = np.array([system["top_occupancy"] for system in systems]) * 100
    any_pair = np.array([system["any_occupancy"] for system in systems]) * 100
    y = np.arange(len(systems))[::-1]
    for y_value, system, fixed_value, any_value in zip(y, systems, fixed, any_pair):
        color = system["color"]
        ax_occupancy.hlines(y_value, fixed_value, any_value, color=color, lw=1.55, zorder=1)
        ax_occupancy.scatter(fixed_value, y_value, s=35, color=color, edgecolor="white", linewidth=0.65, zorder=3)
        ax_occupancy.scatter(
            any_value,
            y_value,
            s=44,
            facecolor="none" if np.isclose(fixed_value, any_value) else "white",
            edgecolor=color,
            linewidth=1.2,
            zorder=4,
        )
        if np.isclose(fixed_value, any_value):
            ax_occupancy.text(any_value + 2.2, y_value, f"{any_value:.0f}%", ha="left", va="center", fontsize=6.3, color=color, fontweight="bold")
        else:
            ax_occupancy.text(fixed_value, y_value + 0.19, f"{fixed_value:.0f}%", ha="center", va="bottom", fontsize=5.8, color=color)
            ax_occupancy.text(any_value, y_value - 0.22, f"{any_value:.0f}%", ha="center", va="top", fontsize=5.8, color=color)
    ax_occupancy.set_xlim(0, 80)
    ax_occupancy.set_ylim(-0.15, 2.45)
    ax_occupancy.set_xticks([0, 20, 40, 60, 80])
    ax_occupancy.set_yticks(y, labels)
    ax_occupancy.set_xlabel("Frames with SG-SG <0.40 nm (%)")
    ax_occupancy.grid(axis="x", color="#DFE6EC", lw=0.55, zorder=0)
    ax_occupancy.set_title("Proximity occupancy", loc="left", fontsize=8.2, fontweight="bold", pad=8)
    ax_occupancy.tick_params(axis="y", length=0)
    panel_label(ax_occupancy, "b")

    fig.suptitle("Inter-trimer Cys proximity in docking MD", fontsize=10.5, fontweight="bold", y=0.985)
    time_subtitle = (
        f"Displayed traces: centered 2.1-ns rolling medians (21 frames) of the per-frame nearest Cys pair ({WINDOW_START_NS:.0f}-{WINDOW_END_NS:.0f} ns)."
        if smoothed
        else f"Each trace reports the nearest cross-trimer Cys pair at every retained 0.1-ns frame ({WINDOW_START_NS:.0f}-{WINDOW_END_NS:.0f} ns)."
    )
    fig.text(
        0.105,
        0.805,
        time_subtitle,
        ha="left",
        fontsize=5.7,
        color="#536578",
    )
    fig.text(
        ax_occupancy.get_position().x0,
        0.805,
        "Filled: fixed candidate pair. Open: any cross-trimer pair.",
        ha="left",
        fontsize=5.7,
        color="#536578",
    )
    bottom_note = (
        f"At each 0.1-ns frame, the nearest SG-SG pair is re-identified across trimers ABC and DEF. All occupancy values use {WINDOW_START_NS:.0f}-{WINDOW_END_NS:.0f} ns. The 0.40-nm cutoff denotes geometric proximity only; it does not establish a covalent S-S bond."
        if not smoothed
        else f"Distances were calculated at each 0.1-ns frame; panel a shows a centered 2.1-ns rolling median. All occupancy values use {WINDOW_START_NS:.0f}-{WINDOW_END_NS:.0f} ns. The 0.40-nm cutoff denotes geometric proximity only; it does not establish a covalent S-S bond."
    )
    fig.text(
        0.5,
        0.055,
        bottom_note,
        ha="center",
        fontsize=5.65,
        color="#46566a",
    )
    base = output_dir / stem
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base.with_suffix(".png"), dpi=500, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(os.environ.get("COLLAGEN_DATA_ROOT", repo_root / "data")))
    parser.add_argument("--output", type=Path, default=Path(os.environ.get("COLLAGEN_RESULTS_ROOT", repo_root / "results")) / "sg_distance")
    args = parser.parse_args()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    systems = [load_system(args.input.resolve(), system) for system in SYSTEMS]
    write_source_data(systems, output_dir)
    draw(systems, output_dir, smoothed=False, stem="Fig_full_trajectory_SG_SG_proximity")
    draw(systems, output_dir, smoothed=True, stem="Fig_full_trajectory_SG_SG_proximity_smoothed")
    print(f"Wrote SG-SG figures and source data to: {output_dir}")


if __name__ == "__main__":
    main()
