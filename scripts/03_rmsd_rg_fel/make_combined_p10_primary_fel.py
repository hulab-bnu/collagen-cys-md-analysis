#!/usr/bin/env python3
"""Combine P10 with CC, mC, NC, and NC3 on one RMSD-Rg FEL grid."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import make_md_fel_rmsd_rg as core

np = core.np
mpl = core.mpl
plt = core.plt
Normalize = core.Normalize


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
INPUT_ROOT = Path(os.environ.get("COLLAGEN_FEL_TIMESERIES", REPO_ROOT / "data" / "rmsd_rg_fel"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "md_fel_combined_p10_primary"
SOURCE_DATA_DIR = OUTPUT_DIR / "source_data"

DATASETS = [
    {
        "key": "P10",
        "label": "P10",
        "path": INPUT_ROOT / "P10_0_rmsd_rg_timeseries.csv",
        "color": "#2b7fad",
    },
    {
        "key": "CC",
        "label": "CC",
        "path": INPUT_ROOT / "CC_rmsd_rg_timeseries.csv",
        "color": "#7251a3",
    },
    {
        "key": "mC",
        "label": "mC",
        "path": INPUT_ROOT / "mC_rmsd_rg_timeseries.csv",
        "color": "#2d6ea3",
    },
    {
        "key": "NC",
        "label": "NC",
        "path": INPUT_ROOT / "NC_rmsd_rg_timeseries.csv",
        "color": "#20854e",
    },
    {
        "key": "NC3",
        "label": "NC3",
        "path": INPUT_ROOT / "NC3_rmsd_rg_timeseries.csv",
        "color": "#b66a18",
    },
]


def read_timeseries(spec: dict) -> dict:
    time = []
    rmsd = []
    rg = []
    include = []
    with spec["path"].open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            time.append(float(row["time_ns"]))
            rmsd.append(float(row["core_backbone_rmsd_nm"]))
            rg.append(float(row["core_radius_of_gyration_nm"]))
            include.append(bool(int(row["included_after_burn_in"])))
    result = {
        "key": spec["key"],
        "label": spec["label"],
        "color": spec["color"],
        "source_file": spec["path"].name,
        "time_ns": np.asarray(time, dtype=float),
        "rmsd_nm": np.asarray(rmsd, dtype=float),
        "rg_nm": np.asarray(rg, dtype=float),
        "include": np.asarray(include, dtype=bool),
    }
    if len(result["time_ns"]) != 20001 or result["time_ns"][-1] < 199.9:
        raise ValueError(f"{spec['label']} is not a complete 200-ns comparison trajectory")
    return result


def occupancy_overlap(probability_a: np.ndarray, probability_b: np.ndarray) -> float:
    return float(np.minimum(probability_a, probability_b).sum())


def summarize(result: dict, p10_probability: np.ndarray) -> dict:
    include = result["include"]
    rmsd = result["rmsd_nm"][include]
    rg = result["rg_nm"][include]
    landscape = result["landscape"]
    early = (result["time_ns"] >= 20.0) & (result["time_ns"] < 110.0)
    late = result["time_ns"] >= 110.0
    return {
        "system": result["label"],
        "analyzed_frames": int(include.sum()),
        "median_rmsd_nm": float(np.median(rmsd)),
        "rmsd_q25_nm": float(np.percentile(rmsd, 25)),
        "rmsd_q75_nm": float(np.percentile(rmsd, 75)),
        "median_rg_nm": float(np.median(rg)),
        "rg_q25_nm": float(np.percentile(rg, 25)),
        "rg_q75_nm": float(np.percentile(rg, 75)),
        "dominant_basin_rmsd_nm": landscape["peak_rmsd_nm"],
        "dominant_basin_rg_nm": landscape["peak_rg_nm"],
        "low_energy_area_nm2_at_2.5_kJ_mol": landscape["low_energy_area_nm2"],
        "occupancy_overlap_with_P10": occupancy_overlap(
            p10_probability, landscape["probability"]
        ),
        "late_minus_early_median_rmsd_nm": float(
            np.median(result["rmsd_nm"][late]) - np.median(result["rmsd_nm"][early])
        ),
        "late_minus_early_median_rg_nm": float(
            np.median(result["rg_nm"][late]) - np.median(result["rg_nm"][early])
        ),
    }


def draw_map(ax, result: dict, panel_label: str, show_xlabel: bool, show_ylabel: bool) -> None:
    landscape = result["landscape"]
    z = np.ma.masked_where(~landscape["supported"], landscape["display_energy"])
    levels = np.linspace(0, core.FREE_ENERGY_CAP_KJ_MOL, 11)
    norm = Normalize(vmin=0, vmax=core.FREE_ENERGY_CAP_KJ_MOL)
    ax.contourf(
        landscape["rmsd_grid"],
        landscape["rg_grid"],
        z,
        levels=levels,
        cmap=core.ENERGY_CMAP,
        norm=norm,
        extend="max",
    )
    ax.contour(
        landscape["rmsd_grid"],
        landscape["rg_grid"],
        z,
        levels=[2.5, 5.0, 7.5, 10.0],
        colors="#29313a",
        linewidths=0.32,
        alpha=0.55,
    )
    ax.scatter(
        landscape["peak_rmsd_nm"],
        landscape["peak_rg_nm"],
        marker="*",
        s=30,
        color="white",
        edgecolor="#1f252b",
        linewidth=0.45,
        zorder=5,
    )
    ax.set_facecolor("#f4f6f7")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.15, 1.04, panel_label, transform=ax.transAxes, fontsize=8.7, fontweight="bold")
    ax.set_title(
        result["label"],
        loc="left",
        fontsize=8.7,
        fontweight="bold",
        color=result["color"],
        pad=3,
    )
    ax.set_xlabel("Core RMSD (nm)" if show_xlabel else "")
    ax.set_ylabel(r"Core $R_g$ (nm)" if show_ylabel else "")
    if not show_xlabel:
        ax.tick_params(labelbottom=False)
    if not show_ylabel:
        ax.tick_params(labelleft=False)


def draw_interval_panel(ax, results: list[dict], summaries: dict[str, dict], metric: str) -> None:
    y = np.arange(len(results))[::-1]
    if metric == "rmsd":
        median_key, q25_key, q75_key = "median_rmsd_nm", "rmsd_q25_nm", "rmsd_q75_nm"
        xlabel = "Core RMSD (nm)"
        title = "Conformational displacement"
    else:
        median_key, q25_key, q75_key = "median_rg_nm", "rg_q25_nm", "rg_q75_nm"
        xlabel = r"Core $R_g$ (nm)"
        title = "Core compactness"

    for position, result in zip(y, results):
        summary = summaries[result["key"]]
        median = summary[median_key]
        low = summary[q25_key]
        high = summary[q75_key]
        ax.plot([low, high], [position, position], color="#929aa3", linewidth=2.2, solid_capstyle="round")
        ax.scatter(
            median,
            position,
            s=28,
            color=result["color"],
            edgecolor="#333940",
            linewidth=0.45,
            zorder=3,
        )
    ax.set_yticks(y, [result["label"] for result in results])
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold", pad=4)
    ax.grid(axis="x", color="#e1e4e7", linewidth=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def draw_overlap_panel(ax, results: list[dict], summaries: dict[str, dict]) -> None:
    comparison = [result for result in results if result["key"] != "P10"]
    comparison.sort(key=lambda item: summaries[item["key"]]["occupancy_overlap_with_P10"])
    y = np.arange(len(comparison))
    values = [summaries[item["key"]]["occupancy_overlap_with_P10"] for item in comparison]
    colors = [item["color"] for item in comparison]
    ax.barh(y, values, color=colors, height=0.56, alpha=0.90)
    ax.set_yticks(y, [item["label"] for item in comparison])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Occupancy overlap with P10")
    ax.set_title("Whole-landscape similarity", loc="left", fontsize=7.5, fontweight="bold", pad=4)
    for position, value in zip(y, values):
        ax.text(value + 0.02, position, f"{value:.2f}", va="center", fontsize=5.8)
    ax.grid(axis="x", color="#e1e4e7", linewidth=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def draw_combined_2d(
    results: list[dict],
    summaries: dict[str, dict],
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 7.0))
    grid = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.0, 1.0, 0.78],
        left=0.09,
        right=0.96,
        bottom=0.08,
        top=0.90,
        wspace=0.24,
        hspace=0.34,
    )
    map_positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    panel_labels = ["a", "b", "c", "d", "e"]
    for index, (result, position, panel_label) in enumerate(
        zip(results, map_positions, panel_labels)
    ):
        ax = fig.add_subplot(grid[position])
        draw_map(
            ax,
            result,
            panel_label,
            show_xlabel=position[0] == 1,
            show_ylabel=position[1] == 0,
        )
        ax.set_xlim(*rmsd_limits)
        ax.set_ylim(*rg_limits)

    info_ax = fig.add_subplot(grid[1, 2])
    info_ax.axis("off")
    norm = Normalize(vmin=0, vmax=core.FREE_ENERGY_CAP_KJ_MOL)
    colorbar_ax = info_ax.inset_axes([0.08, 0.70, 0.84, 0.09])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=core.ENERGY_CMAP),
        cax=colorbar_ax,
        orientation="horizontal",
    )
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", labelpad=2)
    colorbar.outline.set_linewidth(0.6)
    info_ax.text(0.08, 0.93, "Shared analysis", fontsize=7.7, fontweight="bold")
    info_ax.text(
        0.08,
        0.53,
        "20-200 ns  |  18,001 frames/system\n72 x 72 common grid  |  300 K\nidentical Gaussian smoothing",
        fontsize=6.0,
        color="#46515d",
        va="top",
        linespacing=1.45,
    )
    info_ax.text(
        0.08,
        0.25,
        r"$\Delta G=-RT\ln(P/P_{max})$" "\nwhite star: dominant basin",
        fontsize=5.8,
        color="#46515d",
        va="top",
        linespacing=1.5,
    )

    rmsd_ax = fig.add_subplot(grid[2, 0])
    rg_ax = fig.add_subplot(grid[2, 1])
    overlap_ax = fig.add_subplot(grid[2, 2])
    draw_interval_panel(rmsd_ax, results, summaries, "rmsd")
    draw_interval_panel(rg_ax, results, summaries, "rg")
    draw_overlap_panel(overlap_ax, results, summaries)
    rmsd_ax.text(-0.18, 1.04, "f", transform=rmsd_ax.transAxes, fontsize=8.7, fontweight="bold")
    rg_ax.text(-0.18, 1.04, "g", transform=rg_ax.transAxes, fontsize=8.7, fontweight="bold")
    overlap_ax.text(-0.18, 1.04, "h", transform=overlap_ax.transAxes, fontsize=8.7, fontweight="bold")

    fig.suptitle(
        r"P10 benchmarked against CC, mC, NC, and NC3 conformational landscapes",
        x=0.51,
        y=0.975,
        fontsize=10.8,
        fontweight="bold",
    )
    fig.text(
        0.51,
        0.938,
        "All landscapes use the same common-core coordinates, analysis window, grid and free-energy scale",
        ha="center",
        fontsize=6.4,
        color="#46515d",
    )
    fig.text(
        0.51,
        0.018,
        "Points and horizontal lines show medians and interquartile ranges; overlap = shared probability mass with P10.",
        ha="center",
        fontsize=5.9,
        color="#46515d",
    )
    return fig


def draw_combined_3d(
    results: list[dict],
    summaries: dict[str, dict],
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 6.1))
    fig.subplots_adjust(left=0.035, right=0.92, bottom=0.07, top=0.89, wspace=0.02, hspace=0.08)
    norm = Normalize(vmin=0, vmax=core.FREE_ENERGY_CAP_KJ_MOL)

    for index, result in enumerate(results, start=1):
        ax = fig.add_subplot(2, 3, index, projection="3d")
        landscape = result["landscape"]
        z = np.ma.masked_where(~landscape["supported"], landscape["display_energy"])
        ax.plot_surface(
            landscape["rmsd_grid"],
            landscape["rg_grid"],
            z,
            cmap=core.ENERGY_CMAP,
            norm=norm,
            linewidth=0,
            antialiased=True,
            alpha=0.96,
        )
        ax.contourf(
            landscape["rmsd_grid"],
            landscape["rg_grid"],
            z,
            zdir="z",
            offset=-0.75,
            levels=np.linspace(0, core.FREE_ENERGY_CAP_KJ_MOL, 11),
            cmap=core.ENERGY_CMAP,
            norm=norm,
        )
        ax.scatter(
            [landscape["peak_rmsd_nm"]],
            [landscape["peak_rg_nm"]],
            [-0.68],
            s=17,
            marker="*",
            color="white",
            edgecolor="#202020",
            linewidth=0.4,
            depthshade=False,
        )
        ax.set_xlim(*rmsd_limits)
        ax.set_ylim(*rg_limits)
        ax.set_zlim(-0.75, core.FREE_ENERGY_CAP_KJ_MOL)
        ax.set_zticks([0, 5, 10])
        ax.set_xlabel("RMSD (nm)", labelpad=1)
        ax.set_ylabel(r"$R_g$ (nm)", labelpad=1)
        ax.set_zlabel(r"$\Delta G$", labelpad=1)
        ax.tick_params(pad=0, labelsize=5.2)
        ax.view_init(elev=27, azim=-128)
        ax.set_box_aspect((1.08, 1.0, 0.80))
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((1, 1, 1, 0))
            axis.pane.set_edgecolor("#d5d9df")
        ax.text2D(
            0.02, 0.96, chr(96 + index), transform=ax.transAxes,
            fontsize=8.7, fontweight="bold", va="top"
        )
        ax.text2D(
            0.13, 0.96, result["label"], transform=ax.transAxes,
            fontsize=8.8, fontweight="bold", color=result["color"], va="top"
        )
        summary = summaries[result["key"]]
        ax.text2D(
            0.13,
            0.89,
            f"RMSD {summary['median_rmsd_nm']:.3f} nm  |  Rg {summary['median_rg_nm']:.3f} nm",
            transform=ax.transAxes,
            fontsize=5.1,
            color="#46515d",
            va="top",
        )

    info_ax = fig.add_subplot(2, 3, 6)
    info_ax.axis("off")
    colorbar_ax = info_ax.inset_axes([0.10, 0.83, 0.78, 0.06])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=core.ENERGY_CMAP),
        cax=colorbar_ax,
        orientation="horizontal",
    )
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", labelpad=2)
    colorbar.outline.set_linewidth(0.6)
    ranking = sorted(
        [result for result in results if result["key"] != "P10"],
        key=lambda result: summaries[result["key"]]["occupancy_overlap_with_P10"],
        reverse=True,
    )
    info_ax.text(0.10, 0.58, "Similarity to P10", fontsize=7.7, fontweight="bold")
    for row, result in enumerate(ranking):
        overlap = summaries[result["key"]]["occupancy_overlap_with_P10"]
        info_ax.text(0.10, 0.47 - row * 0.10, result["label"], fontsize=6.5, color=result["color"], fontweight="bold")
        info_ax.text(0.78, 0.47 - row * 0.10, f"{overlap:.2f}", fontsize=6.5, ha="right")
    info_ax.text(
        0.10,
        0.05,
        "Overlap compares the full\nprobability landscape, not\nonly the median RMSD.",
        fontsize=5.5,
        color="#46515d",
        linespacing=1.4,
    )

    fig.suptitle(
        r"P10, CC, mC, NC, and NC3 RMSD-$R_g$ free-energy landscapes",
        x=0.49,
        y=0.97,
        fontsize=10.8,
        fontweight="bold",
    )
    fig.text(
        0.49,
        0.93,
        "Shared common-core axes and relative free-energy scale; 20-200 ns",
        ha="center",
        fontsize=6.3,
        color="#46515d",
    )
    return fig


def write_outputs(results: list[dict], summary_rows: list[dict], limits) -> None:
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    core.write_csv(
        SOURCE_DATA_DIR / "combined_summary_metrics.csv",
        list(summary_rows[0].keys()),
        summary_rows,
    )
    for result in results:
        landscape = result["landscape"]
        rows = []
        for i, rmsd in enumerate(landscape["rmsd_centers"]):
            for j, rg in enumerate(landscape["rg_centers"]):
                energy = landscape["free_energy"][i, j]
                rows.append(
                    {
                        "core_backbone_rmsd_nm": f"{rmsd:.7f}",
                        "core_radius_of_gyration_nm": f"{rg:.7f}",
                        "smoothed_probability": f"{landscape['probability'][i, j]:.10g}",
                        "delta_G_kJ_mol": "" if not np.isfinite(energy) else f"{energy:.7f}",
                        "display_supported": int(landscape["supported"][i, j]),
                    }
                )
        core.write_csv(
            SOURCE_DATA_DIR / f"{result['key']}_shared_grid.csv",
            [
                "core_backbone_rmsd_nm",
                "core_radius_of_gyration_nm",
                "smoothed_probability",
                "delta_G_kJ_mol",
                "display_supported",
            ],
            rows,
        )

    manifest = {
        "analysis": "Combined P10/CC/mC/NC/NC3 common-core RMSD-Rg landscapes",
        "systems": [result["label"] for result in results],
        "settings": {
            "analysis_window_ns": [20.0, 200.0],
            "frames_per_system": 18001,
            "common_core_residues_per_chain": 30,
            "bins_per_axis": core.N_BINS,
            "smoothing_sigma_bins": core.SMOOTH_SIGMA_BINS,
            "temperature_K": core.TEMPERATURE_K,
            "free_energy_cap_kJ_mol": core.FREE_ENERGY_CAP_KJ_MOL,
            "shared_rmsd_limits_nm": limits[0],
            "shared_rg_limits_nm": limits[1],
            "occupancy_overlap": "sum over common grid of min(P10 probability, comparison probability)",
        },
    }
    (OUTPUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    legend = """Figure X | P10 benchmarked against CC, mC, NC, and NC3 conformational landscapes.

a-e, Two-dimensional common-core RMSD-Rg projected free-energy landscapes for P10, CC, mC, NC, and NC3. All systems were analyzed over 20-200 ns using 18,001 frames, an identical 72 x 72 grid, and the same Gaussian smoothing and color scale. Relative free energy was calculated as DeltaG(x,y) = -RT ln[P(x,y)/Pmax] at 300 K. White stars mark the most populated basin. f,g, Median and interquartile range of common-core backbone RMSD (f) and mass-weighted radius of gyration (g). h, Probability overlap between each construct and P10, calculated as the summed minimum probability over the common RMSD-Rg grid. Larger values indicate more similar conformational sampling. One trajectory was analyzed per construct.
"""
    (OUTPUT_DIR / "figure_legend.txt").write_text(legend, encoding="utf-8")


def main() -> None:
    core.configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    results = [read_timeseries(spec) for spec in DATASETS]
    limits = core.common_limits(results)
    for result in results:
        core.build_landscape(result, *limits)

    p10_probability = results[0]["landscape"]["probability"]
    summary_rows = [summarize(result, p10_probability) for result in results]
    summaries = {
        result["key"]: summary for result, summary in zip(results, summary_rows)
    }
    write_outputs(results, summary_rows, limits)
    core.save_figure(
        draw_combined_2d(results, summaries, *limits),
        OUTPUT_DIR / "combined_p10_cc_mc_nc_nc3_fel_2d_quantitative",
    )
    core.save_figure(
        draw_combined_3d(results, summaries, *limits),
        OUTPUT_DIR / "combined_p10_cc_mc_nc_nc3_fel_3d",
    )
    print(f"Done: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
