#!/usr/bin/env python3
"""Integrate MM/PBSA endpoint estimates with cross-trimer Cys SG geometry.

Inputs are the three completed docking-MD result folders (CC_1, mC, NC).
The script creates publication-ready SVG/PDF/TIFF/PNG figures and the exact
source-data tables used for plotting. All MM/PBSA values are reported in
kcal/mol and all SG-SG distances in nm.
"""

from __future__ import annotations

import csv
import io
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("COLLAGEN_DATA_ROOT", REPO_ROOT / "data"))
OUT = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "mmpbsa_cys_geometry"

SYSTEMS = {
    "CC_1": {"label": "CC", "color": "#7450A8"},
    "mC": {"label": "mC", "color": "#2C78B7"},
    "NC": {"label": "NC", "color": "#278D5B"},
}
SYSTEM_ORDER = list(SYSTEMS)
BLOCK_EDGES = np.array([0, 40, 80, 120, 160, 201])
ENERGY_COLUMNS = ["VDWAALS", "EEL", "EPB", "ENPOLAR", "GGAS", "GSOLV", "TOTAL"]
HEATMAP_COLUMNS = ["VDWAALS", "EEL", "EPB", "ENPOLAR", "TOTAL"]
INTERNAL_COMPONENTS = ["BOND", "ANGLE", "DIHED", "UB", "IMP", "CMAP"]

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def read_delta_energy(path: Path) -> pd.DataFrame:
    """Read only the Delta Energy Terms block from gmx_MMPBSA output."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "Delta Energy Terms")
    header = start + 1
    end = next((i for i in range(header + 1, len(lines)) if not lines[i].strip()), len(lines))
    frame = pd.read_csv(io.StringIO("\n".join(lines[header:end])))
    frame = frame.apply(pd.to_numeric, errors="coerce").dropna(subset=["Frame #"])
    frame["time_ns"] = frame["Frame #"].astype(float) - 1.0
    return frame


def read_delta_decomposition(path: Path) -> pd.DataFrame:
    """Read the DELTAS Total Decomposition Contribution block only."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    delta_start = next(i for i, line in enumerate(lines) if line.strip() == "DELTAS:")
    header = next(i for i in range(delta_start, len(lines)) if lines[i].startswith("Frame #,"))
    end = next(
        i
        for i in range(header + 1, len(lines))
        if lines[i].strip() == "Sidechain Decomposition Contribution (SDC)"
    )
    frame = pd.read_csv(io.StringIO("\n".join(lines[header:end])))
    for col in frame.columns:
        if col not in {"Residue"}:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame.dropna(subset=["Frame #", "TOTAL"])


def read_geometry(system_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    geometry_dir = system_dir / "disulfide_geometry"
    pairs = pd.read_csv(geometry_dir / "cross_trimer_sg_sg_pair_summary.csv")
    contacts = pd.read_csv(geometry_dir / "potential_disulfide_contact_count.csv")
    distances = pd.read_csv(geometry_dir / "sg_sg_pair_distance_timeseries.csv")
    return pairs, contacts, distances


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.10, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames_all: list[pd.DataFrame] = []
    blocks_all: list[pd.DataFrame] = []
    energy_summary: list[dict[str, object]] = []
    geometry_summary: list[dict[str, object]] = []
    residue_summary: list[pd.DataFrame] = []
    closest_timeseries: dict[str, pd.DataFrame] = {}

    for system in SYSTEM_ORDER:
        system_dir = ROOT / system
        metadata = SYSTEMS[system]

        delta = read_delta_energy(system_dir / "mmpbsa_ABC_vs_DEF" / "FINAL_RESULTS_MMPBSA.csv")
        delta["system"] = metadata["label"]
        frames_all.append(delta)

        for block_id, (start, stop) in enumerate(zip(BLOCK_EDGES[:-1], BLOCK_EDGES[1:]), start=1):
            subset = delta[(delta["time_ns"] >= start) & (delta["time_ns"] < stop)]
            for component in ENERGY_COLUMNS:
                blocks_all.append(
                    {
                        "system": metadata["label"],
                        "block": block_id,
                        "start_ns": start,
                        "end_ns": stop - 1,
                        "component": component,
                        "mean_kcal_mol": subset[component].mean(),
                        "sd_kcal_mol": subset[component].std(ddof=1),
                        "n_frames": len(subset),
                    }
                )

        for component in [*ENERGY_COLUMNS, *INTERNAL_COMPONENTS]:
            energy_summary.append(
                {
                    "system": metadata["label"],
                    "component": component,
                    "mean_kcal_mol": delta[component].mean(),
                    "sd_frames_kcal_mol": delta[component].std(ddof=1),
                    "sem_frames_kcal_mol": delta[component].sem(ddof=1),
                }
            )

        decomp = read_delta_decomposition(system_dir / "mmpbsa_ABC_vs_DEF" / "FINAL_DECOMP_MMPBSA.csv")
        residue = (
            decomp.groupby("Residue", as_index=False)
            .agg(
                mean_total_kcal_mol=("TOTAL", "mean"),
                sd_total_kcal_mol=("TOTAL", "std"),
                n_frames=("TOTAL", "count"),
            )
            .sort_values("mean_total_kcal_mol")
        )
        residue["system"] = metadata["label"]
        residue_summary.append(residue)

        pairs, contacts, distances = read_geometry(system_dir)
        best_pair = pairs.sort_values("median_sg_sg_nm", kind="stable").iloc[0]
        best_name = str(best_pair["pair"])
        best_ts = distances[(distances["pair_class"] == "between_trimers") & (distances["pair"] == best_name)].copy()
        best_ts["rolling_median_nm"] = best_ts["sg_sg_distance_nm"].rolling(21, center=True, min_periods=1).median()
        closest_timeseries[system] = best_ts

        geometry_summary.append(
            {
                "system": metadata["label"],
                "fixed_pair": best_name,
                "median_fixed_pair_nm": best_pair["median_sg_sg_nm"],
                "q1_fixed_pair_nm": best_pair["q1_sg_sg_nm"],
                "q3_fixed_pair_nm": best_pair["q3_sg_sg_nm"],
                "fixed_pair_lt_0p40_pct": 100.0 * best_pair["fraction_lt_0p40_nm"],
                "any_cross_trimer_lt_0p40_pct": 100.0
                * (contacts["raw_cross_trimer_pairs_lt_cutoff"] > 0).mean(),
                "mean_cross_trimer_contact_pairs": contacts["raw_cross_trimer_pairs_lt_cutoff"].mean(),
            }
        )

    frames = pd.concat(frames_all, ignore_index=True)
    blocks = pd.DataFrame(blocks_all)
    energy = pd.DataFrame(energy_summary)
    geometry = pd.DataFrame(geometry_summary)
    residues = pd.concat(residue_summary, ignore_index=True)

    energy_wide = energy.pivot(index="system", columns="component", values="mean_kcal_mol")
    block_total = blocks[blocks["component"] == "TOTAL"].groupby("system")["mean_kcal_mol"].agg(["mean", "min", "max"])
    qc = geometry.set_index("system").join(energy_wide[["TOTAL", *INTERNAL_COMPONENTS]]).join(block_total, how="left", rsuffix="_block")
    qc["internal_residual_kcal_mol"] = qc[INTERNAL_COMPONENTS].sum(axis=1)
    qc["abs_internal_over_abs_total_pct"] = 100.0 * qc["internal_residual_kcal_mol"].abs() / qc["TOTAL"].abs()
    qc = qc.reset_index().rename(
        columns={
            "TOTAL": "delta_g_bind_kcal_mol",
            "mean": "block_mean_kcal_mol",
            "min": "block_min_kcal_mol",
            "max": "block_max_kcal_mol",
        }
    )

    frames.to_csv(OUT / "source_framewise_mmpbsa_delta.csv", index=False)
    blocks.to_csv(OUT / "source_mmpbsa_40ns_blocks.csv", index=False)
    energy.to_csv(OUT / "source_mmpbsa_energy_summary.csv", index=False)
    geometry.to_csv(OUT / "source_cross_trimer_cys_geometry.csv", index=False)
    residues.to_csv(OUT / "source_residue_decomposition_summary.csv", index=False)
    qc.to_csv(OUT / "source_integrated_mmpbsa_geometry_qc_summary.csv", index=False)

    # Main figure: energy quality and SG geometry are deliberately shown together.
    fig = plt.figure(figsize=(10.4, 7.65))
    grid = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.05], wspace=0.42, hspace=0.48)
    fig.subplots_adjust(left=0.09, right=0.965, top=0.88, bottom=0.13)
    ax_time = fig.add_subplot(grid[0, 0])
    ax_blocks = fig.add_subplot(grid[0, 1])
    ax_heat = fig.add_subplot(grid[1, 0])
    ax_sg = fig.add_subplot(grid[1, 1])

    # Panel A: all 201 framewise endpoint estimates plus a 10-ns rolling average.
    for system in SYSTEM_ORDER:
        meta = SYSTEMS[system]
        subset = frames[frames["system"] == meta["label"]].sort_values("time_ns")
        smooth = subset["TOTAL"].rolling(11, center=True, min_periods=1).mean()
        ax_time.plot(subset["time_ns"], subset["TOTAL"], color=meta["color"], lw=0.5, alpha=0.20)
        ax_time.plot(subset["time_ns"], smooth, color=meta["color"], lw=1.8, label=meta["label"])
    ax_time.axhline(0, color="#6F7680", lw=0.75, zorder=0)
    ax_time.set(xlim=(0, 200), xlabel="Simulation time (ns)", ylabel=r"$\Delta G_{\mathrm{PB}}$ (kcal mol$^{-1}$)")
    ax_time.set_title("Frame-resolved MM/PBSA estimate", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_time.legend(loc="lower left", ncol=3, handlelength=2.4, columnspacing=1.2)
    add_panel_label(ax_time, "a")

    # Panel B: five contiguous 40-ns blocks, the appropriate temporal replicate display.
    system_labels = [SYSTEMS[s]["label"] for s in SYSTEM_ORDER]
    for xpos, system in enumerate(SYSTEM_ORDER):
        meta = SYSTEMS[system]
        values = blocks[(blocks["system"] == meta["label"]) & (blocks["component"] == "TOTAL")]["mean_kcal_mol"].to_numpy()
        offsets = np.linspace(-0.13, 0.13, len(values))
        ax_blocks.scatter(np.full(len(values), xpos) + offsets, values, s=25, color=meta["color"], edgecolor="white", linewidth=0.5, zorder=3)
        mean = values.mean()
        lo, hi = values.min(), values.max()
        ax_blocks.errorbar(xpos, mean, yerr=[[mean - lo], [hi - mean]], fmt="D", ms=5.2, color="#222222", mfc="white", capsize=3, lw=1.0, zorder=4)
        ax_blocks.text(xpos, np.max(values) + 3.8, f"{mean:.1f}", ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax_blocks.axhline(0, color="#6F7680", lw=0.75, zorder=0)
    ax_blocks.set(xticks=np.arange(len(system_labels)), xticklabels=system_labels, ylabel=r"Block mean $\Delta G_{\mathrm{PB}}$ (kcal mol$^{-1}$)")
    ax_blocks.set_title("Temporal block robustness", loc="left", fontsize=10, fontweight="bold", pad=8)
    add_panel_label(ax_blocks, "b")

    # Panel C: physical component fingerprint. Blue is favorable (negative), red unfavorable (positive).
    heat = energy.pivot(index="system", columns="component", values="mean_kcal_mol").reindex(index=system_labels, columns=HEATMAP_COLUMNS)
    vmax = np.ceil(np.abs(heat.to_numpy()).max() / 50) * 50
    image = ax_heat.imshow(heat.to_numpy(), cmap="RdYlBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax_heat.set(xticks=np.arange(len(HEATMAP_COLUMNS)), xticklabels=[r"$\Delta E_{vdW}$", r"$\Delta E_{ele}$", r"$\Delta G_{PB}$", r"$\Delta G_{nonpol}$", r"$\Delta G_{bind}$"], yticks=np.arange(len(system_labels)), yticklabels=system_labels)
    ax_heat.tick_params(axis="x", rotation=28)
    for row in range(heat.shape[0]):
        for col in range(heat.shape[1]):
            value = heat.iat[row, col]
            color = "white" if abs(value) > 0.56 * vmax else "#1D232D"
            ax_heat.text(col, row, f"{value:.1f}", ha="center", va="center", fontsize=7.2, color=color, fontweight="bold")
    cbar = fig.colorbar(image, ax=ax_heat, fraction=0.046, pad=0.04)
    cbar.set_label("Mean energy (kcal mol$^{-1}$)")
    cbar.set_ticks([-vmax, 0, vmax])
    cbar.set_ticklabels(["Favorable", "0", "Unfavorable"])
    ax_heat.set_title("Energetic fingerprint", loc="left", fontsize=10, fontweight="bold", pad=8)
    add_panel_label(ax_heat, "c")

    # Panel D: closest persistent inter-trimer Cys pair, preserving all 0.1-ns sampled frames in the raw data.
    for system in SYSTEM_ORDER:
        meta = SYSTEMS[system]
        row = geometry[geometry["system"] == meta["label"]].iloc[0]
        ts = closest_timeseries[system]
        ax_sg.plot(ts["time_ns"], ts["sg_sg_distance_nm"], color=meta["color"], lw=0.35, alpha=0.18)
        ax_sg.plot(
            ts["time_ns"],
            ts["rolling_median_nm"],
            color=meta["color"],
            lw=1.55,
            label=f"{meta['label']}: {row['fixed_pair_lt_0p40_pct']:.0f}% <0.40 nm",
        )
    ax_sg.axhline(0.40, color="#C7821D", lw=1.0, ls="--")
    ax_sg.text(2, 0.415, "SG proximity threshold (0.40 nm)", fontsize=7, color="#9A6113", va="bottom")
    ax_sg.set(xlim=(0, 200), ylim=(0.25, 1.18), xlabel="Simulation time (ns)", ylabel="Fixed-pair SG-SG distance (nm)")
    ax_sg.set_title("Persistent cross-trimer Cys geometry", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_sg.legend(loc="upper left", fontsize=7, handlelength=2.0, frameon=True, facecolor="white", edgecolor="#D9DEE5")
    add_panel_label(ax_sg, "d")

    fig.suptitle("Docking MD couples topology-specific energetics to Cys proximity", fontsize=13, fontweight="bold", y=0.975)
    fig.text(
        0.09,
        0.035,
        "a, 1-ns endpoint estimates (thin) and 11-ns rolling means (thick). b, Five contiguous 40-ns means; diamond, mean; whiskers, block range. c, Mean energy terms; negative is favorable. d, Fixed cross-trimer SG pair with lowest median distance. SG-SG <0.40 nm denotes proximity, not a formed S-S bond.",
        fontsize=7,
        color="#56616D",
    )
    save_figure(fig, OUT / "Fig_MMPBSA_Cys_geometry")
    plt.close(fig)

    # Supplementary figure: top interface residue contributions, retaining the native gmx_MMPBSA identifiers.
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 5.3), sharex=False)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.79, bottom=0.14, wspace=0.30)
    for ax, system in zip(axes, SYSTEM_ORDER):
        meta = SYSTEMS[system]
        table = residues[residues["system"] == meta["label"]].nsmallest(12, "mean_total_kcal_mol").sort_values("mean_total_kcal_mol")
        labels = table["Residue"].str.replace(r"^[RL]:", "", regex=True)
        is_cys = labels.str.contains("CYS")
        colors = np.where(is_cys, "#C7821D", meta["color"])
        ax.barh(np.arange(len(table)), table["mean_total_kcal_mol"], color=colors, height=0.72)
        ax.axvline(0, color="#6F7680", lw=0.7)
        left_limit = min(-0.5, table["mean_total_kcal_mol"].min() * 1.12)
        ax.set(xlim=(left_limit, 0), yticks=np.arange(len(table)), yticklabels=labels, xlabel="Mean contribution (kcal mol$^{-1}$)")
        ax.tick_params(axis="y", labelsize=7)
        ax.set_title(meta["label"], color=meta["color"], fontweight="bold", fontsize=11)
    fig.text(0.16, 0.95, "Interface-residue MM/PBSA decomposition", fontsize=13, fontweight="bold", va="top")
    fig.text(0.16, 0.90, "Twelve most favorable residues per system; cysteines are highlighted in amber. Identifiers follow the native gmx_MMPBSA mapping; each panel has an independent x-axis scale.", fontsize=8, color="#56616D")
    save_figure(fig, OUT / "Fig_MMPBSA_interface_residue_decomposition")
    plt.close(fig)

    # Human-readable method note.
    note = """MM/PBSA and Cys-geometry analysis summary

Systems: CC, mC, NC.
Trajectory window: 0-200 ns, 201 frames for MM/PBSA (1-ns stride).
Temporal uncertainty display: five contiguous 40-ns block means.
MM/PBSA method: single-trajectory CHARMM-PB endpoint calculation (PBRadii=7, radiopt=0).
Geometry definition: SG-SG <0.40 nm is a proximity event, not a formed covalent disulfide bond.
Interpretation boundary: use cross-system ranking and physical component trends; do not interpret endpoint values as experimental absolute binding free energies.
"""
    (OUT / "README.txt").write_text(note, encoding="utf-8")


if __name__ == "__main__":
    main()
