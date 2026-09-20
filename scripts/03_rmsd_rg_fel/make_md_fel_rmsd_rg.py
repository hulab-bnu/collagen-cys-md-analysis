#!/usr/bin/env python3
"""Build comparable RMSD-Rg projected free-energy landscapes from four MD runs.

The four source models have different residue counts. To avoid confounding Rg
with construct length, the analysis uses a centered 30-residue core from each
of the three collagen chains in every system.
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "md_fel_rmsd_rg"
SOURCE_DATA_DIR = OUTPUT_DIR / "source_data"
# MDAnalysis' XDR reader cannot open the private-use character in the CC
# source directory. A temporary ASCII-only cache is used for that system.
CACHE_DIR = Path(os.environ.get("COLLAGEN_CACHE_DIR", Path(tempfile.gettempdir()) / "collagen_cys_md_cache"))
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".mplconfig"))
(OUTPUT_DIR / ".mplconfig").mkdir(parents=True, exist_ok=True)

LOCAL_DEPS = HERE / ".mdanalysis-extracted"
if LOCAL_DEPS.exists():
    # Keep the locally bundled binary stack internally consistent. In
    # particular, do not mix its NumPy with unrelated global binary packages.
    sys.path = [
        str(LOCAL_DEPS),
        *[path for path in sys.path if "site-packages" not in path.lower()],
    ]

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.ndimage import gaussian_filter, maximum_filter

import MDAnalysis as mda


SOURCE_ROOT = Path(os.environ.get("COLLAGEN_DATA_ROOT", REPO_ROOT / "data"))

SYSTEMS = [
    {"key": "CC", "label": "CC", "folder_prefix": "CC"},
    {"key": "mC", "label": "mC", "folder": "mC ok"},
    {"key": "NC", "label": "NC", "folder": "NC"},
    {"key": "NC3", "label": "NC3", "folder": "NC3_C2"},
]

CORE_RESIDUES_PER_CHAIN = 30
BURN_IN_NS = 20.0
TEMPERATURE_K = 300.0
GAS_CONSTANT_KJ_MOL_K = 0.008314462618
RT = GAS_CONSTANT_KJ_MOL_K * TEMPERATURE_K
N_BINS = 72
SMOOTH_SIGMA_BINS = 1.15
FREE_ENERGY_CAP_KJ_MOL = 12.5
LOW_ENERGY_CUTOFF_KJ_MOL = 2.5

ENERGY_CMAP = LinearSegmentedColormap.from_list(
    "energy_sequential",
    ["#382a84", "#2878b8", "#45b6b0", "#f0d264", "#ef8738", "#c93435"],
    N=256,
)

SYSTEM_COLORS = {
    "mC": "#2d6ea3",
    "CC": "#6f4aa1",
    "NC": "#20854e",
    "NC3": "#b66a18",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.0,
            "axes.labelsize": 7.0,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def resolve_source_folder(spec: dict[str, str]) -> Path:
    if "folder" in spec:
        folder = SOURCE_ROOT / spec["folder"]
        if not folder.is_dir():
            raise FileNotFoundError(folder)
        return folder

    matches = sorted(
        path for path in SOURCE_ROOT.iterdir()
        if path.is_dir() and path.name.startswith(spec["folder_prefix"])
    )
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one folder beginning with {spec['folder_prefix']!r}; found {matches}"
        )
    return matches[0]


def contains_private_use_character(path: Path) -> bool:
    return any(0xE000 <= ord(char) <= 0xF8FF for char in str(path))


def prepare_input_files(spec: dict[str, str]) -> tuple[Path, Path, Path]:
    source_folder = resolve_source_folder(spec)
    topology = source_folder / "1_processed.gro"
    trajectory = source_folder / "fit.xtc"
    if not topology.is_file() or not trajectory.is_file():
        raise FileNotFoundError(f"Missing 1_processed.gro or fit.xtc in {source_folder}")

    # The CC directory contains a private-use character that the XDR C library
    # cannot open. Cache only that pair under a conventional path.
    if contains_private_use_character(source_folder):
        local = CACHE_DIR / spec["key"]
        local.mkdir(parents=True, exist_ok=True)
        local_topology = local / "topology.gro"
        local_trajectory = local / "trajectory.xtc"
        for source, destination in (
            (topology, local_topology),
            (trajectory, local_trajectory),
        ):
            if not destination.exists() or destination.stat().st_size != source.stat().st_size:
                shutil.copy2(source, destination)
        topology, trajectory = local_topology, local_trajectory

    return source_folder, topology, trajectory


def ordered_chain_residue_groups(universe: mda.Universe) -> list:
    residues = universe.residues
    if len(residues) % 3 != 0:
        raise ValueError(f"Residue count {len(residues)} is not divisible into three chains")
    chain_length = len(residues) // 3
    return [residues[i * chain_length : (i + 1) * chain_length] for i in range(3)]


def build_common_core(universe: mda.Universe) -> tuple:
    chain_groups = ordered_chain_residue_groups(universe)
    chain_length = len(chain_groups[0])
    if chain_length < CORE_RESIDUES_PER_CHAIN:
        raise ValueError(
            f"Chain length {chain_length} is shorter than requested core "
            f"({CORE_RESIDUES_PER_CHAIN})"
        )
    trim_left = (chain_length - CORE_RESIDUES_PER_CHAIN) // 2
    selected_residues = []
    for residues in chain_groups:
        selected_residues.extend(
            residues[trim_left : trim_left + CORE_RESIDUES_PER_CHAIN]
        )
    atom_indices = np.concatenate([residue.atoms.indices for residue in selected_residues])
    core = universe.atoms[atom_indices]
    backbone = core.select_atoms("backbone")
    if len(backbone) == 0:
        raise ValueError("Backbone atom selection is empty")
    return core, backbone, chain_length, trim_left


def weighted_kabsch_rmsd(
    mobile: np.ndarray,
    reference_centered: np.ndarray,
    weights: np.ndarray,
) -> float:
    weight_sum = weights.sum()
    mobile_center = np.sum(mobile * weights[:, None], axis=0) / weight_sum
    mobile_centered = mobile - mobile_center
    covariance = (mobile_centered * weights[:, None]).T @ reference_centered
    left, _, right_t = np.linalg.svd(covariance)
    if np.linalg.det(left @ right_t) < 0:
        left[:, -1] *= -1
    rotation = left @ right_t
    delta = mobile_centered @ rotation - reference_centered
    return float(np.sqrt(np.sum(weights * np.sum(delta * delta, axis=1)) / weight_sum))


def radius_of_gyration(coordinates: np.ndarray, masses: np.ndarray) -> float:
    mass_sum = masses.sum()
    center = np.sum(coordinates * masses[:, None], axis=0) / mass_sum
    squared_distance = np.sum((coordinates - center) ** 2, axis=1)
    return float(np.sqrt(np.sum(masses * squared_distance) / mass_sum))


def extract_timeseries(spec: dict[str, str]) -> dict:
    source_folder, topology, trajectory = prepare_input_files(spec)
    universe = mda.Universe(str(topology), str(trajectory))
    core, backbone, chain_length, trim_left = build_common_core(universe)

    if not np.all(np.isfinite(core.masses)) or np.any(core.masses <= 0):
        raise ValueError(f"Missing or invalid atomic masses in {spec['key']}")

    universe.trajectory[0]
    backbone_weights = backbone.masses.astype(np.float64)
    reference = backbone.positions.astype(np.float64).copy()
    reference_center = np.sum(reference * backbone_weights[:, None], axis=0) / backbone_weights.sum()
    reference_centered = reference - reference_center
    core_masses = core.masses.astype(np.float64)

    frame_count = len(universe.trajectory)
    time_ns = np.empty(frame_count, dtype=np.float64)
    rmsd_nm = np.empty(frame_count, dtype=np.float64)
    rg_nm = np.empty(frame_count, dtype=np.float64)

    for index, timestep in enumerate(universe.trajectory):
        time_ns[index] = timestep.time / 1000.0
        rmsd_nm[index] = weighted_kabsch_rmsd(
            backbone.positions.astype(np.float64),
            reference_centered,
            backbone_weights,
        ) / 10.0
        rg_nm[index] = radius_of_gyration(
            core.positions.astype(np.float64), core_masses
        ) / 10.0

    include = time_ns >= BURN_IN_NS
    if include.sum() < 1000:
        raise ValueError(f"Too few post-equilibration frames for {spec['key']}")

    return {
        "key": spec["key"],
        "label": spec["label"],
        "source_folder_name": spec.get("folder", spec["key"]),
        "topology_name": "1_processed.gro",
        "trajectory_name": "fit.xtc",
        "chain_length": chain_length,
        "trim_left": trim_left,
        "core_residues": CORE_RESIDUES_PER_CHAIN * 3,
        "core_atoms": len(core),
        "backbone_atoms": len(backbone),
        "frames": frame_count,
        "dt_ps": float(universe.trajectory.dt),
        "time_ns": time_ns,
        "rmsd_nm": rmsd_nm,
        "rg_nm": rg_nm,
        "include": include,
    }


def common_limits(results: list[dict]) -> tuple[tuple[float, float], tuple[float, float]]:
    all_rmsd = np.concatenate([item["rmsd_nm"][item["include"]] for item in results])
    all_rg = np.concatenate([item["rg_nm"][item["include"]] for item in results])

    def padded(values: np.ndarray) -> tuple[float, float]:
        low = float(np.min(values))
        high = float(np.max(values))
        pad = max((high - low) * 0.035, 0.005)
        return low - pad, high + pad

    return padded(all_rmsd), padded(all_rg)


def build_landscape(
    result: dict,
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> dict:
    rmsd = result["rmsd_nm"][result["include"]]
    rg = result["rg_nm"][result["include"]]
    histogram, rmsd_edges, rg_edges = np.histogram2d(
        rmsd,
        rg,
        bins=N_BINS,
        range=[rmsd_limits, rg_limits],
    )
    probability = gaussian_filter(histogram.astype(float), sigma=SMOOTH_SIGMA_BINS)
    probability /= probability.sum()
    probability_ratio = probability / probability.max()
    with np.errstate(divide="ignore", invalid="ignore"):
        free_energy = -RT * np.log(probability_ratio)
    supported = np.isfinite(free_energy) & (free_energy <= FREE_ENERGY_CAP_KJ_MOL)
    display_energy = np.where(
        np.isfinite(free_energy),
        np.minimum(free_energy, FREE_ENERGY_CAP_KJ_MOL),
        FREE_ENERGY_CAP_KJ_MOL,
    )

    rmsd_centers = 0.5 * (rmsd_edges[:-1] + rmsd_edges[1:])
    rg_centers = 0.5 * (rg_edges[:-1] + rg_edges[1:])
    rmsd_grid, rg_grid = np.meshgrid(rmsd_centers, rg_centers, indexing="ij")
    peak_index = np.unravel_index(np.argmax(probability), probability.shape)

    peak_mask = (
        (probability == maximum_filter(probability, size=5, mode="nearest"))
        & supported
        & (free_energy <= 5.0)
    )
    peak_candidates = np.argwhere(peak_mask)
    peak_candidates = sorted(
        peak_candidates,
        key=lambda index: probability[tuple(index)],
        reverse=True,
    )
    separated_peaks = []
    for index in peak_candidates:
        if all(np.linalg.norm(index - chosen) >= 5 for chosen in separated_peaks):
            separated_peaks.append(index)
        if len(separated_peaks) == 5:
            break

    cell_area = (rmsd_edges[1] - rmsd_edges[0]) * (rg_edges[1] - rg_edges[0])
    low_energy_area = float(np.sum(free_energy <= LOW_ENERGY_CUTOFF_KJ_MOL) * cell_area)

    landscape = {
        "histogram": histogram,
        "probability": probability,
        "free_energy": free_energy,
        "display_energy": display_energy,
        "supported": supported,
        "rmsd_centers": rmsd_centers,
        "rg_centers": rg_centers,
        "rmsd_grid": rmsd_grid,
        "rg_grid": rg_grid,
        "peak_rmsd_nm": float(rmsd_centers[peak_index[0]]),
        "peak_rg_nm": float(rg_centers[peak_index[1]]),
        "basin_count": len(separated_peaks),
        "low_energy_area_nm2": low_energy_area,
    }
    result["landscape"] = landscape
    return landscape


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q))


def summarize(result: dict) -> dict:
    rmsd = result["rmsd_nm"][result["include"]]
    rg = result["rg_nm"][result["include"]]
    midpoint_ns = 0.5 * (BURN_IN_NS + float(result["time_ns"][-1]))
    early = (result["time_ns"] >= BURN_IN_NS) & (result["time_ns"] < midpoint_ns)
    late = result["time_ns"] >= midpoint_ns
    early_rmsd = float(np.median(result["rmsd_nm"][early]))
    late_rmsd = float(np.median(result["rmsd_nm"][late]))
    early_rg = float(np.median(result["rg_nm"][early]))
    late_rg = float(np.median(result["rg_nm"][late]))
    landscape = result["landscape"]
    return {
        "system": result["label"],
        "source_folder": result["source_folder_name"],
        "trajectory": result["trajectory_name"],
        "total_frames": result["frames"],
        "analyzed_frames": int(result["include"].sum()),
        "trajectory_dt_ps": result["dt_ps"],
        "burn_in_ns": BURN_IN_NS,
        "analysis_end_ns": float(result["time_ns"][-1]),
        "original_residues_per_chain": result["chain_length"],
        "trimmed_residues_each_N_side": result["trim_left"],
        "common_core_residues_total": result["core_residues"],
        "common_core_atoms": result["core_atoms"],
        "common_core_backbone_atoms": result["backbone_atoms"],
        "median_rmsd_nm": float(np.median(rmsd)),
        "rmsd_q25_nm": percentile(rmsd, 25),
        "rmsd_q75_nm": percentile(rmsd, 75),
        "median_rg_nm": float(np.median(rg)),
        "rg_q25_nm": percentile(rg, 25),
        "rg_q75_nm": percentile(rg, 75),
        "early_20_110_ns_median_rmsd_nm": early_rmsd,
        "late_110_200_ns_median_rmsd_nm": late_rmsd,
        "late_minus_early_median_rmsd_nm": late_rmsd - early_rmsd,
        "early_20_110_ns_median_rg_nm": early_rg,
        "late_110_200_ns_median_rg_nm": late_rg,
        "late_minus_early_median_rg_nm": late_rg - early_rg,
        "dominant_basin_rmsd_nm": landscape["peak_rmsd_nm"],
        "dominant_basin_rg_nm": landscape["peak_rg_nm"],
        "low_energy_area_nm2_at_2.5_kJ_mol": landscape["low_energy_area_nm2"],
        "estimated_basin_count_below_5_kJ_mol": landscape["basin_count"],
    }


def write_csv(path: Path, fieldnames: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_source_data(results: list[dict], summaries: list[dict]) -> None:
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    for result in results:
        rows = (
            {
                "time_ns": f"{time:.5f}",
                "core_backbone_rmsd_nm": f"{rmsd:.7f}",
                "core_radius_of_gyration_nm": f"{rg:.7f}",
                "included_after_burn_in": int(include),
            }
            for time, rmsd, rg, include in zip(
                result["time_ns"],
                result["rmsd_nm"],
                result["rg_nm"],
                result["include"],
            )
        )
        write_csv(
            SOURCE_DATA_DIR / f"{result['key']}_rmsd_rg_timeseries.csv",
            [
                "time_ns",
                "core_backbone_rmsd_nm",
                "core_radius_of_gyration_nm",
                "included_after_burn_in",
            ],
            rows,
        )

        landscape = result["landscape"]
        grid_rows = []
        for i, rmsd in enumerate(landscape["rmsd_centers"]):
            for j, rg in enumerate(landscape["rg_centers"]):
                value = landscape["free_energy"][i, j]
                grid_rows.append(
                    {
                        "core_backbone_rmsd_nm": f"{rmsd:.7f}",
                        "core_radius_of_gyration_nm": f"{rg:.7f}",
                        "smoothed_probability": f"{landscape['probability'][i, j]:.10g}",
                        "delta_G_kJ_mol": "" if not np.isfinite(value) else f"{value:.7f}",
                        "display_supported": int(landscape["supported"][i, j]),
                    }
                )
        write_csv(
            SOURCE_DATA_DIR / f"{result['key']}_free_energy_grid.csv",
            [
                "core_backbone_rmsd_nm",
                "core_radius_of_gyration_nm",
                "smoothed_probability",
                "delta_G_kJ_mol",
                "display_supported",
            ],
            grid_rows,
        )

    write_csv(
        SOURCE_DATA_DIR / "summary_metrics.csv",
        list(summaries[0].keys()),
        summaries,
    )


def format_panel_subtitle(summary: dict) -> str:
    return (
        f"median RMSD {summary['median_rmsd_nm']:.3f} nm  |  "
        f"median Rg {summary['median_rg_nm']:.3f} nm"
    )


def draw_3d_figure(
    results: list[dict],
    summaries_by_key: dict[str, dict],
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 6.35), constrained_layout=False)
    fig.subplots_adjust(left=0.045, right=0.89, bottom=0.075, top=0.91, wspace=0.02, hspace=0.10)
    norm = Normalize(vmin=0, vmax=FREE_ENERGY_CAP_KJ_MOL)
    panel_labels = ["a", "b", "c", "d"]

    for panel_index, result in enumerate(results, start=1):
        ax = fig.add_subplot(2, 2, panel_index, projection="3d")
        landscape = result["landscape"]
        z = np.ma.masked_where(~landscape["supported"], landscape["display_energy"])
        ax.plot_surface(
            landscape["rmsd_grid"],
            landscape["rg_grid"],
            z,
            cmap=ENERGY_CMAP,
            norm=norm,
            rstride=1,
            cstride=1,
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
            levels=np.linspace(0, FREE_ENERGY_CAP_KJ_MOL, 11),
            cmap=ENERGY_CMAP,
            norm=norm,
            antialiased=True,
        )
        ax.scatter(
            [landscape["peak_rmsd_nm"]],
            [landscape["peak_rg_nm"]],
            [-0.68],
            s=18,
            marker="*",
            color="white",
            edgecolor="#202020",
            linewidth=0.45,
            depthshade=False,
            zorder=10,
        )
        ax.set_xlim(*rmsd_limits)
        ax.set_ylim(*rg_limits)
        ax.set_zlim(-0.75, FREE_ENERGY_CAP_KJ_MOL)
        ax.set_zticks([0, 5, 10])
        ax.set_xlabel("Core RMSD (nm)", labelpad=3)
        ax.set_ylabel(r"Core $R_g$ (nm)", labelpad=3)
        ax.set_zlabel(r"$\Delta G$ (kJ mol$^{-1}$)", labelpad=3)
        ax.tick_params(axis="both", which="major", pad=0, length=2)
        ax.view_init(elev=27, azim=-128)
        ax.set_box_aspect((1.10, 1.00, 0.82))
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((1, 1, 1, 0))
            axis.pane.set_edgecolor("#d5d9df")
        ax.text2D(
            0.03,
            0.94,
            panel_labels[panel_index - 1],
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
        )
        ax.text2D(
            0.12,
            0.94,
            result["label"],
            transform=ax.transAxes,
            fontsize=9.5,
            fontweight="bold",
            color=SYSTEM_COLORS[result["key"]],
            va="top",
        )
        ax.text2D(
            0.12,
            0.885,
            format_panel_subtitle(summaries_by_key[result["key"]]),
            transform=ax.transAxes,
            fontsize=5.8,
            color="#4a5562",
            va="top",
        )

    colorbar_ax = fig.add_axes([0.915, 0.25, 0.018, 0.47])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=ENERGY_CMAP), cax=colorbar_ax)
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", labelpad=5)
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.outline.set_linewidth(0.6)

    fig.suptitle(
        r"RMSD-$R_g$ projected free-energy landscapes of collagen-like constructs",
        x=0.47,
        y=0.975,
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.47,
        0.942,
        "Common 30-residue core per chain; 20-200 ns; identical bins, smoothing and color scale",
        ha="center",
        va="center",
        fontsize=6.8,
        color="#46515d",
    )
    fig.text(
        0.47,
        0.018,
        r"$\Delta G(x,y)=-RT\ln[P(x,y)/P_{max}]$ at 300 K; white star marks the most populated basin.",
        ha="center",
        va="bottom",
        fontsize=6.2,
        color="#46515d",
    )
    return fig


def draw_2d_figure(
    results: list[dict],
    summaries_by_key: dict[str, dict],
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.75), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.095, right=0.88, bottom=0.11, top=0.89, wspace=0.16, hspace=0.23)
    levels = np.linspace(0, FREE_ENERGY_CAP_KJ_MOL, 11)
    panel_labels = ["a", "b", "c", "d"]
    norm = Normalize(vmin=0, vmax=FREE_ENERGY_CAP_KJ_MOL)

    for ax, result, panel_label in zip(axes.flat, results, panel_labels):
        landscape = result["landscape"]
        z = np.ma.masked_where(~landscape["supported"], landscape["display_energy"])
        ax.contourf(
            landscape["rmsd_grid"],
            landscape["rg_grid"],
            z,
            levels=levels,
            cmap=ENERGY_CMAP,
            norm=norm,
            extend="max",
        )
        ax.contour(
            landscape["rmsd_grid"],
            landscape["rg_grid"],
            z,
            levels=[2.5, 5.0, 7.5, 10.0],
            colors="#27313a",
            linewidths=0.38,
            alpha=0.58,
        )
        ax.scatter(
            landscape["peak_rmsd_nm"],
            landscape["peak_rg_nm"],
            marker="*",
            s=34,
            color="white",
            edgecolor="#202020",
            linewidth=0.5,
            zorder=5,
        )
        ax.set_xlim(*rmsd_limits)
        ax.set_ylim(*rg_limits)
        ax.set_facecolor("#f5f6f7")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.13, 1.04, panel_label, transform=ax.transAxes, fontsize=9, fontweight="bold")
        ax.set_title(
            result["label"],
            loc="left",
            color=SYSTEM_COLORS[result["key"]],
            fontweight="bold",
            fontsize=9.5,
            pad=4,
        )
        summary = summaries_by_key[result["key"]]
        ax.text(
            0.99,
            1.02,
            f"basin: {summary['dominant_basin_rmsd_nm']:.3f}, {summary['dominant_basin_rg_nm']:.3f} nm",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.8,
            color="#4a5562",
        )

    for ax in axes[-1, :]:
        ax.set_xlabel("Common-core backbone RMSD (nm)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Common-core $R_g$ (nm)")

    colorbar_ax = fig.add_axes([0.905, 0.20, 0.018, 0.57])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=ENERGY_CMAP), cax=colorbar_ax)
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", labelpad=5)
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.outline.set_linewidth(0.6)

    fig.suptitle(
        r"Comparable RMSD-$R_g$ projected free-energy landscapes",
        x=0.49,
        y=0.965,
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.49,
        0.925,
        "Common 30-residue core per chain; equilibrated 20-200 ns segment",
        ha="center",
        fontsize=6.8,
        color="#46515d",
    )
    fig.text(
        0.49,
        0.025,
        "Low-energy basins indicate frequently sampled conformations; blank regions are not sufficiently sampled.",
        ha="center",
        fontsize=6.2,
        color="#46515d",
    )
    return fig


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_methods_and_manifest(
    results: list[dict],
    summaries: list[dict],
    rmsd_limits: tuple[float, float],
    rg_limits: tuple[float, float],
) -> None:
    manifest = {
        "analysis": "RMSD-Rg projected free-energy landscape",
        "software": {
            "python": sys.version.split()[0],
            "MDAnalysis": mda.__version__,
            "numpy": np.__version__,
            "matplotlib": mpl.__version__,
        },
        "systems": [
            {
                "display_name": result["label"],
                "source_folder_name": result["source_folder_name"],
                "topology": result["topology_name"],
                "trajectory": result["trajectory_name"],
            }
            for result in results
        ],
        "settings": {
            "common_core_residues_per_chain": CORE_RESIDUES_PER_CHAIN,
            "chains": 3,
            "rmsd_atoms": "backbone atoms of centered common core",
            "rmsd_reference": "first trajectory frame (0 ns)",
            "rmsd_alignment": "mass-weighted Kabsch superposition on the same backbone atoms",
            "rg_atoms": "all atoms of centered common core",
            "rg_weighting": "atomic mass",
            "burn_in_ns": BURN_IN_NS,
            "temperature_K": TEMPERATURE_K,
            "bins_per_axis": N_BINS,
            "gaussian_smoothing_sigma_bins": SMOOTH_SIGMA_BINS,
            "free_energy_cap_kJ_mol": FREE_ENERGY_CAP_KJ_MOL,
            "low_energy_basin_cutoff_kJ_mol": LOW_ENERGY_CUTOFF_KJ_MOL,
            "shared_rmsd_limits_nm": rmsd_limits,
            "shared_rg_limits_nm": rg_limits,
            "free_energy_equation": "DeltaG = -RT ln(P/Pmax)",
        },
        "summary_metrics_file": "source_data/summary_metrics.csv",
    }
    with (OUTPUT_DIR / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    methods = f"""# RMSD-Rg projected free-energy landscape

## Analysis definition

All four 200-ns trajectories were sampled every 10 ps. The first {BURN_IN_NS:.0f} ns were excluded as equilibration, leaving {int((200 - BURN_IN_NS) * 100 + 1):,} frames per system. Because the source models contain different numbers of residues per chain, each chain was center-trimmed to a common {CORE_RESIDUES_PER_CHAIN}-residue core ({CORE_RESIDUES_PER_CHAIN * 3} residues across the collagen triple helix).

Core-backbone RMSD was calculated after mass-weighted Kabsch superposition to the 0-ns common-core backbone:

`RMSD = sqrt[sum_i m_i ||R(r_i-r_cm) - (r_i^ref-r_cm^ref)||^2 / sum_i m_i]`

The core radius of gyration used all atoms in the same common core:

`Rg = sqrt[sum_i m_i ||r_i-r_cm||^2 / sum_i m_i]`

The joint RMSD-Rg distribution was estimated on a shared {N_BINS} x {N_BINS} grid and smoothed with a Gaussian kernel (sigma = {SMOOTH_SIGMA_BINS:.2f} bins). Relative free energy was calculated at {TEMPERATURE_K:.0f} K:

`DeltaG(x,y) = -RT ln[P(x,y)/Pmax]`

All systems use identical RMSD/Rg limits, grid dimensions, smoothing, and a 0-{FREE_ENERGY_CAP_KJ_MOL:.1f} kJ mol^-1 color scale. Regions above the display cap are treated as insufficiently sampled and left blank in the 2D map.

Early-versus-late median RMSD and Rg values (20-110 ns versus 110-200 ns) are reported in `source_data/summary_metrics.csv` as a simple stationarity diagnostic.

## Interpretation limits

This is a two-dimensional projection of finite MD sampling, not an absolute binding free energy. Basin positions and occupancy breadth can be compared under the shared protocol; apparent barriers and basin counts should be interpreted cautiously because only one trajectory is available for each construct. Rg is made comparable by using the same centered core, but sequence/topology differences can still affect the coordinate projection.

## Recommended manuscript use

Use the 2D map as the main comparative panel because basin position and occupancy are easier to read without perspective distortion. The 3D surface matches the requested Figure-f visual style and is well suited to a slide or supplementary visualization.
"""
    (OUTPUT_DIR / "methods_and_interpretation.md").write_text(methods, encoding="utf-8")


def main() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    for spec in SYSTEMS:
        print(f"Reading and analyzing {spec['label']}...", flush=True)
        results.append(extract_timeseries(spec))

    rmsd_limits, rg_limits = common_limits(results)
    for result in results:
        build_landscape(result, rmsd_limits, rg_limits)

    summaries = [summarize(result) for result in results]
    summaries_by_key = {result["key"]: summary for result, summary in zip(results, summaries)}
    write_source_data(results, summaries)
    write_methods_and_manifest(results, summaries, rmsd_limits, rg_limits)

    print("Rendering 3D surface figure...", flush=True)
    save_figure(
        draw_3d_figure(results, summaries_by_key, rmsd_limits, rg_limits),
        OUTPUT_DIR / "figure_f_style_rmsd_rg_fel_3d",
    )
    print("Rendering 2D comparison figure...", flush=True)
    save_figure(
        draw_2d_figure(results, summaries_by_key, rmsd_limits, rg_limits),
        OUTPUT_DIR / "rmsd_rg_fel_2d_comparison",
    )
    print(f"Done: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
