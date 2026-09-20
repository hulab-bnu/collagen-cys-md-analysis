#!/usr/bin/env python3
"""Summarize and plot single-trajectory CHARMM-MM/PBSA results for ABC vs DEF."""

from __future__ import annotations

import argparse
import csv
import os
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


# Keep SVG text editable for later layout adjustment.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)

SYSTEMS = ("P10", "CC_1", "mC", "NC")
SYSTEM_LABELS = {"P10": "P10", "CC_1": "CC", "mC": "mC", "NC": "NC"}
# Global construct palette shared with the screening, PCA, DCCM, and structure panels.
SYSTEM_COLORS = {"P10": "#C84E3B", "CC_1": "#7652A8", "mC": "#2879B8", "NC": "#238A56"}
GAS_COLOR = "#3977B8"
SOLV_COLOR = "#E49B3C"
NEUTRAL = "#545E6A"
FAVOURABLE = "#2F7FA8"
UNFAVOURABLE = "#C84E3B"

ENERGY_COLUMNS = [
    "BOND",
    "ANGLE",
    "DIHED",
    "UB",
    "IMP",
    "CMAP",
    "VDWAALS",
    "EEL",
    "1-4 VDW",
    "1-4 EEL",
    "EPB",
    "ENPOLAR",
    "EDISPER",
    "GGAS",
    "GSOLV",
    "TOTAL",
]
INTERNAL_COLUMNS = ["BOND", "ANGLE", "DIHED", "UB", "IMP", "CMAP"]
MECHANISM_COLUMNS = ["VDWAALS", "EEL", "EPB", "ENPOLAR"]


def parse_delta_section(path: Path) -> pd.DataFrame:
    """Read only the Delta (complex - receptor - ligand) block from gmx_MMPBSA CSV."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == "Delta Energy Terms")
    except StopIteration as exc:
        raise ValueError(f"Delta Energy Terms section not found in {path}") from exc

    header_index = start + 1
    data_lines: list[str] = []
    for line in lines[header_index:]:
        if not line.strip():
            break
        data_lines.append(line)

    delta = pd.read_csv(StringIO("\n".join(data_lines)))
    expected = ["Frame #", *ENERGY_COLUMNS]
    missing = sorted(set(expected).difference(delta.columns))
    if missing:
        raise ValueError(f"Missing fields in {path}: {', '.join(missing)}")
    return delta[expected].copy()


def add_time_and_blocks(delta: pd.DataFrame) -> pd.DataFrame:
    """Map frame 1..201 to 0..200 ns and assign contiguous 40-ns stability blocks."""
    result = delta.copy()
    result["time_ns"] = result["Frame #"].astype(float) - 1.0
    bins = [0, 40, 80, 120, 160, 201]
    labels = ["0-40", "40-80", "80-120", "120-160", "160-200"]
    result["block"] = pd.cut(
        result["time_ns"],
        bins=bins,
        labels=labels,
        right=False,
        include_lowest=True,
    )
    return result


def bootstrap_mean_ci(values: np.ndarray, n_boot: int = 5000, seed: int = 20260807) -> tuple[float, float]:
    """Descriptive bootstrap interval for five contiguous block means."""
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(n_boot, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def build_summaries(data_by_system: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries = []
    block_rows = []
    per_frame_rows = []

    for system, data in data_by_system.items():
        for block, subset in data.groupby("block", observed=False):
            if subset.empty:
                continue
            block_rows.append(
                {
                    "system": SYSTEM_LABELS[system],
                    "block": str(block),
                    "n_frames": len(subset),
                    "mean_delta_total_kcal_mol": subset["TOTAL"].mean(),
                    "sd_delta_total_kcal_mol": subset["TOTAL"].std(ddof=1),
                }
            )

        block_means = (
            data.groupby("block", observed=False)["TOTAL"].mean().dropna().to_numpy(dtype=float)
        )
        ci_low, ci_high = bootstrap_mean_ci(block_means)
        summary = {
            "system": SYSTEM_LABELS[system],
            "n_frames": len(data),
            "n_time_blocks": len(block_means),
            "delta_total_mean_kcal_mol": data["TOTAL"].mean(),
            "delta_total_frame_sd_kcal_mol": data["TOTAL"].std(ddof=1),
            "delta_total_block_sd_kcal_mol": block_means.std(ddof=1),
            "delta_total_block_bootstrap_ci_low": ci_low,
            "delta_total_block_bootstrap_ci_high": ci_high,
            "delta_ggas_mean_kcal_mol": data["GGAS"].mean(),
            "delta_gsolv_mean_kcal_mol": data["GSOLV"].mean(),
            "internal_residual_kcal_mol": data[INTERNAL_COLUMNS].sum(axis=1).mean(),
        }
        summary["internal_residual_fraction_of_abs_total"] = (
            abs(summary["internal_residual_kcal_mol"])
            / abs(summary["delta_total_mean_kcal_mol"])
        )
        for column in MECHANISM_COLUMNS:
            summary[f"delta_{column.lower().replace('-', '_').replace(' ', '_')}_mean_kcal_mol"] = data[column].mean()
        summaries.append(summary)

        frame_out = data[["Frame #", "time_ns", "block", *ENERGY_COLUMNS]].copy()
        frame_out.insert(0, "system", SYSTEM_LABELS[system])
        per_frame_rows.append(frame_out)

    return pd.DataFrame(summaries), pd.DataFrame(block_rows), pd.concat(per_frame_rows, ignore_index=True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.04, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def plot_mmpbsa_figure(
    data_by_system: dict[str, pd.DataFrame],
    summary: pd.DataFrame,
    blocks: pd.DataFrame,
    out_base: Path,
) -> None:
    """Create an evidence-ranked publication figure and a separate QC panel."""
    fig = plt.figure(figsize=(10.6, 6.7))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.75, 1.0],
        height_ratios=[1.12, 1.0],
        left=0.075,
        right=0.985,
        top=0.885,
        bottom=0.135,
        hspace=0.68,
        wspace=0.52,
    )
    ax_time = fig.add_subplot(grid[0, 0])
    ax_balance = fig.add_subplot(grid[0, 1])
    ax_blocks = fig.add_subplot(grid[1, 0])
    ax_heat = fig.add_subplot(grid[1, 1])

    fig.text(
        0.075,
        0.945,
        "Inter-trimer association energy: P10 baseline and Cys-topology series",
        fontsize=14,
        fontweight="bold",
        ha="left",
    )
    fig.text(
        0.075,
        0.913,
        "P10-P10 is the Cys-free baseline; ABC receptor versus DEF ligand; 201 frames across 0-200 ns; single-trajectory CHARMM-PB protocol",
        fontsize=8.5,
        color=NEUTRAL,
        ha="left",
    )

    # a. Time-resolved endpoint association estimates.
    for system in SYSTEMS:
        data = data_by_system[system]
        color = SYSTEM_COLORS[system]
        rolling = data["TOTAL"].rolling(window=11, center=True, min_periods=1).mean()
        ax_time.plot(data["time_ns"], data["TOTAL"], color=color, alpha=0.15, lw=0.7)
        ax_time.plot(data["time_ns"], rolling, color=color, lw=1.8, label=SYSTEM_LABELS[system])
    ax_time.axhline(0, color="#9AA1A8", lw=0.7, ls="--", zorder=0)
    ax_time.set_xlim(0, 200)
    ax_time.set_xlabel("Simulation time (ns)")
    ax_time.set_ylabel(r"$\Delta G_{\mathrm{MM/PBSA}}$ (kcal mol$^{-1}$)")
    ax_time.set_title("Time-resolved association estimate", loc="left", fontsize=9.5, fontweight="bold")
    ax_time.legend(ncol=4, loc="lower left", fontsize=7.3, handlelength=2.2, columnspacing=0.9)
    ax_time.text(
        0.99,
        0.06,
        "Thin: 1-ns frames\nThick: 11-ns rolling mean",
        transform=ax_time.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color=NEUTRAL,
    )
    add_panel_label(ax_time, "a")

    # b. Gas versus solvent compensation. Separate offsets avoid overlapping bars.
    system_labels = [SYSTEM_LABELS[system] for system in SYSTEMS]
    label_to_system = {label: system for system, label in SYSTEM_LABELS.items()}
    summary_plot = summary.set_index("system").loc[system_labels].reset_index()
    y = np.arange(len(summary_plot))[::-1]
    for yi, row in zip(y, summary_plot.itertuples(index=False)):
        ax_balance.barh(
            yi + 0.19,
            row.delta_ggas_mean_kcal_mol,
            color=GAS_COLOR,
            height=0.26,
            label=r"$\Delta G_{\mathrm{gas}}$" if yi == y[0] else None,
        )
        ax_balance.barh(
            yi - 0.19,
            row.delta_gsolv_mean_kcal_mol,
            color=SOLV_COLOR,
            height=0.26,
            label=r"$\Delta G_{\mathrm{solv}}$" if yi == y[0] else None,
        )
        ax_balance.scatter(
            row.delta_total_mean_kcal_mol,
            yi,
            s=30,
            color=SYSTEM_COLORS[label_to_system[row.system]],
            edgecolor="white",
            linewidth=0.7,
            zorder=4,
        )
        ax_balance.text(
            row.delta_total_mean_kcal_mol + 8,
            yi,
            f"{row.delta_total_mean_kcal_mol:.1f}",
            va="center",
            fontsize=7,
        )
    ax_balance.axvline(0, color="#9AA1A8", lw=0.7)
    ax_balance.set_yticks(y, summary_plot["system"])
    ax_balance.set_xlabel(r"Energy contribution (kcal mol$^{-1}$)")
    ax_balance.set_title("Gas-solvent compensation", loc="left", fontsize=9.5, fontweight="bold", pad=23)
    ax_balance.text(
        0.02,
        1.07,
        r"$\blacksquare$ $\Delta G_{\mathrm{gas}}$",
        transform=ax_balance.transAxes,
        fontsize=7.1,
        color=GAS_COLOR,
        ha="left",
        va="bottom",
    )
    ax_balance.text(
        0.52,
        1.07,
        r"$\blacksquare$ $\Delta G_{\mathrm{solv}}$",
        transform=ax_balance.transAxes,
        fontsize=7.1,
        color=SOLV_COLOR,
        ha="left",
        va="bottom",
    )
    add_panel_label(ax_balance, "b")

    # c. Five contiguous time-block means. They show stability, not independent replicates.
    block_order = ["0-40", "40-80", "80-120", "120-160", "160-200"]
    x = np.arange(len(system_labels))
    for xi, system_label in enumerate(system_labels):
        subset = blocks[blocks["system"] == system_label].set_index("block").loc[block_order]
        values = subset["mean_delta_total_kcal_mol"].to_numpy()
        color = SYSTEM_COLORS[label_to_system[system_label]]
        jitter = np.linspace(-0.13, 0.13, len(values))
        ax_blocks.plot(xi + jitter, values, color=color, alpha=0.45, lw=0.9, zorder=1)
        ax_blocks.scatter(xi + jitter, values, s=22, color=color, edgecolor="white", linewidth=0.5, zorder=2)
        mean = values.mean()
        lo, hi = values.min(), values.max()
        ax_blocks.plot([xi, xi], [lo, hi], color=color, lw=2.2, zorder=0)
        ax_blocks.scatter(xi, mean, s=46, color="white", edgecolor=color, linewidth=1.4, zorder=3)
        ax_blocks.text(xi, hi + 2.2, f"{mean:.1f}", ha="center", va="bottom", fontsize=7, color=color, fontweight="bold")
    ax_blocks.axhline(0, color="#9AA1A8", lw=0.7, ls="--")
    ax_blocks.set_xticks(x, system_labels)
    ax_blocks.set_ylabel(r"Block mean $\Delta G$ (kcal mol$^{-1}$)")
    ax_blocks.set_title("Stability across five 40-ns blocks", loc="left", fontsize=9.5, fontweight="bold")
    add_panel_label(ax_blocks, "c")

    # d. Mechanistic components use a signed, labelled heat map.
    heat_labels = [r"$\Delta E_{\mathrm{vdW}}$", r"$\Delta E_{\mathrm{ele}}$", r"$\Delta G_{\mathrm{PB}}$", r"$\Delta G_{\mathrm{nonpolar}}$"]
    matrix = np.array(
        [
            [
                summary_plot.loc[summary_plot["system"] == system_label, f"delta_{column.lower().replace('-', '_').replace(' ', '_')}_mean_kcal_mol"].iloc[0]
                for system_label in system_labels
            ]
            for column in MECHANISM_COLUMNS
        ]
    )
    limit = float(np.ceil(np.abs(matrix).max() / 25.0) * 25.0)
    # Match the screening heat-map language: blue/teal is lower and more
    # favourable, while yellow/orange/red is higher and more unfavourable.
    cmap = LinearSegmentedColormap.from_list(
        "screening_matched_energy",
        [
            (0.00, "#2F7FA8"),
            (0.36, "#67B6B0"),
            (0.50, "#F2DE7D"),
            (0.74, "#EC9447"),
            (1.00, "#C84E3B"),
        ],
    )
    im = ax_heat.imshow(matrix, cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect="auto")
    ax_heat.set_xticks(range(len(system_labels)), system_labels)
    ax_heat.set_yticks(range(4), heat_labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            text_color = "white" if abs(matrix[i, j]) > 0.56 * limit else "#20262D"
            ax_heat.text(j, i, f"{matrix[i, j]:.0f}", ha="center", va="center", fontsize=7.5, color=text_color)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.055, pad=0.045)
    cbar.set_label(r"Favourable $\leftarrow$ contribution $\rightarrow$ unfavourable (kcal mol$^{-1}$)", fontsize=6.6)
    cbar.ax.tick_params(labelsize=6.5)
    ax_heat.set_title("Mechanistic energy signature", loc="left", fontsize=9.5, fontweight="bold")
    add_panel_label(ax_heat, "d")

    fig.text(
        0.075,
        0.040,
        "Dots: contiguous 40-ns blocks; open dot: overall mean; range bar: min-max. Endpoint estimates derive from one trajectory per system; "
        "blocks are a time-stability check, not independent biological replicates. "
        "Use cross-system ranking and energy balance rather than the absolute magnitude as the principal inference.",
        fontsize=6.5,
        color=NEUTRAL,
        ha="left",
    )

    for suffix, kwargs in {
        "svg": {},
        "pdf": {},
        "png": {"dpi": 350},
        "tiff": {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}},
    }.items():
        fig.savefig(out_base.with_suffix(f".{suffix}"), bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path(os.environ.get("COLLAGEN_DATA_ROOT", repo_root / "data")))
    parser.add_argument("--output-dir", type=Path, default=Path(os.environ.get("COLLAGEN_RESULTS_ROOT", repo_root / "results")) / "mmpbsa_analysis")
    args = parser.parse_args()

    data_by_system: dict[str, pd.DataFrame] = {}
    for system in SYSTEMS:
        source = args.input_root / system / "mmpbsa_ABC_vs_DEF" / "FINAL_RESULTS_MMPBSA.csv"
        if not source.exists():
            raise FileNotFoundError(f"Missing MM/PBSA result: {source}")
        data_by_system[system] = add_time_and_blocks(parse_delta_section(source))

    summary, blocks, per_frame = build_summaries(data_by_system)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "mmpbsa_system_summary.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)
    blocks.to_csv(args.output_dir / "mmpbsa_block_summary.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)
    per_frame.to_csv(args.output_dir / "mmpbsa_delta_per_frame.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)

    plot_mmpbsa_figure(data_by_system, summary, blocks, args.output_dir / "Fig_MMPBSA_intertrimer_comparison")

    print("MM/PBSA summary")
    print(summary.round(3).to_string(index=False))
    print(f"\nOutputs: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
