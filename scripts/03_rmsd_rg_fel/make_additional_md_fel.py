#!/usr/bin/env python3
"""Analyze additional MD trajectories discovered beside the primary four runs."""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import make_md_fel_rmsd_rg as core
from MDAnalysis.lib.distances import minimize_vectors

np = core.np
mpl = core.mpl
plt = core.plt
Normalize = core.Normalize


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "md_fel_additional"
SOURCE_DATA_DIR = OUTPUT_DIR / "source_data"

ADDITIONAL_SYSTEMS = [
    {"key": "mC_original", "label": "mC (original)", "folder": "mC"},
    {"key": "P10_0", "label": "P10_0", "folder": "P10_0"},
    {"key": "P10_1", "label": "P10_1", "folder": "P10_1"},
]

COLORS = {
    "mC_original": "#a63d40",
    "P10_0": "#2878a8",
    "P10_1": "#8a64a8",
}


def configure_output() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    core.configure_style()


def extract_clustered_timeseries(spec: dict[str, str]) -> dict:
    """Recluster the three chains by minimum image before RMSD/Rg analysis."""
    source_folder, topology, trajectory = core.prepare_input_files(spec)
    universe = core.mda.Universe(str(topology), str(trajectory))
    residue_groups = core.ordered_chain_residue_groups(universe)
    chain_length = len(residue_groups[0])
    trim_left = (chain_length - core.CORE_RESIDUES_PER_CHAIN) // 2
    chain_cores = [
        residues[trim_left : trim_left + core.CORE_RESIDUES_PER_CHAIN].atoms
        for residues in residue_groups
    ]
    chain_backbones = [atoms.select_atoms("backbone") for atoms in chain_cores]
    core_masses = np.concatenate([atoms.masses for atoms in chain_cores]).astype(float)
    backbone_masses = np.concatenate([atoms.masses for atoms in chain_backbones]).astype(float)
    if np.any(core_masses <= 0) or np.any(backbone_masses <= 0):
        raise ValueError(f"Invalid masses in {spec['key']}")

    def clustered_coordinates() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        box = universe.trajectory.ts.dimensions
        if box is None or not np.all(np.isfinite(box)):
            raise ValueError(f"Missing unit-cell dimensions in {spec['key']}")
        coms = np.asarray([atoms.center_of_mass() for atoms in chain_cores])
        shifts = [np.zeros(3, dtype=float)]
        for chain_index in (1, 2):
            vector = np.asarray(coms[chain_index] - coms[0], dtype=np.float32).reshape(1, 3)
            minimum_vector = minimize_vectors(vector, box)[0].astype(float)
            shifts.append(coms[0] + minimum_vector - coms[chain_index])
        corrected_core = np.concatenate(
            [atoms.positions.astype(float) + shift for atoms, shift in zip(chain_cores, shifts)]
        )
        corrected_backbone = np.concatenate(
            [atoms.positions.astype(float) + shift for atoms, shift in zip(chain_backbones, shifts)]
        )
        corrected_coms = coms + np.asarray(shifts)

        return corrected_core, corrected_backbone, coms, corrected_coms

    universe.trajectory[0]
    _, reference_backbone, _, _ = clustered_coordinates()
    reference_center = np.sum(
        reference_backbone * backbone_masses[:, None], axis=0
    ) / backbone_masses.sum()
    reference_centered = reference_backbone - reference_center

    frame_count = len(universe.trajectory)
    time_ns = np.empty(frame_count, dtype=float)
    rmsd_nm = np.empty(frame_count, dtype=float)
    rg_nm = np.empty(frame_count, dtype=float)
    original_max_separation_nm = np.empty(frame_count, dtype=float)
    clustered_max_separation_nm = np.empty(frame_count, dtype=float)

    for index, timestep in enumerate(universe.trajectory):
        corrected_core, corrected_backbone, original_coms, corrected_coms = clustered_coordinates()
        time_ns[index] = timestep.time / 1000.0
        rmsd_nm[index] = core.weighted_kabsch_rmsd(
            corrected_backbone, reference_centered, backbone_masses
        ) / 10.0
        rg_nm[index] = core.radius_of_gyration(corrected_core, core_masses) / 10.0
        original_distances = [
            np.linalg.norm(original_coms[0] - original_coms[1]),
            np.linalg.norm(original_coms[0] - original_coms[2]),
            np.linalg.norm(original_coms[1] - original_coms[2]),
        ]
        corrected_distances = [
            np.linalg.norm(corrected_coms[0] - corrected_coms[1]),
            np.linalg.norm(corrected_coms[0] - corrected_coms[2]),
            np.linalg.norm(corrected_coms[1] - corrected_coms[2]),
        ]
        original_max_separation_nm[index] = max(original_distances) / 10.0
        clustered_max_separation_nm[index] = max(corrected_distances) / 10.0

    include = time_ns >= core.BURN_IN_NS
    return {
        "key": spec["key"],
        "label": spec["label"],
        "source_folder_name": spec["folder"],
        "topology_name": "1_processed.gro",
        "trajectory_name": "fit.xtc",
        "chain_length": chain_length,
        "trim_left": trim_left,
        "core_residues": core.CORE_RESIDUES_PER_CHAIN * 3,
        "core_atoms": int(sum(len(atoms) for atoms in chain_cores)),
        "backbone_atoms": int(sum(len(atoms) for atoms in chain_backbones)),
        "frames": frame_count,
        "dt_ps": float(universe.trajectory.dt),
        "time_ns": time_ns,
        "rmsd_nm": rmsd_nm,
        "rg_nm": rg_nm,
        "include": include,
        "original_max_interchain_com_nm": original_max_separation_nm,
        "clustered_max_interchain_com_nm": clustered_max_separation_nm,
    }


def summarize_additional(result: dict) -> dict:
    include = result["include"]
    time = result["time_ns"]
    rmsd = result["rmsd_nm"][include]
    rg = result["rg_nm"][include]
    midpoint = 0.5 * (core.BURN_IN_NS + float(time[-1]))
    early = (time >= core.BURN_IN_NS) & (time < midpoint)
    late = time >= midpoint
    landscape = result["individual_landscape"]
    original_separation = result["original_max_interchain_com_nm"]
    clustered_separation = result["clustered_max_interchain_com_nm"]
    status = "PBC clustering corrected"
    if time[-1] < 199.9:
        status += "; processed trajectory truncated"
    return {
        "system": result["label"],
        "source_folder": result["source_folder_name"],
        "total_frames": result["frames"],
        "analyzed_frames": int(include.sum()),
        "analysis_start_ns": core.BURN_IN_NS,
        "analysis_end_ns": float(time[-1]),
        "half_split_ns": midpoint,
        "median_rmsd_nm": float(np.median(rmsd)),
        "rmsd_q25_nm": float(np.percentile(rmsd, 25)),
        "rmsd_q75_nm": float(np.percentile(rmsd, 75)),
        "median_rg_nm": float(np.median(rg)),
        "rg_q25_nm": float(np.percentile(rg, 25)),
        "rg_q75_nm": float(np.percentile(rg, 75)),
        "late_minus_early_median_rmsd_nm": float(
            np.median(result["rmsd_nm"][late]) - np.median(result["rmsd_nm"][early])
        ),
        "late_minus_early_median_rg_nm": float(
            np.median(result["rg_nm"][late]) - np.median(result["rg_nm"][early])
        ),
        "dominant_basin_rmsd_nm": landscape["peak_rmsd_nm"],
        "dominant_basin_rg_nm": landscape["peak_rg_nm"],
        "low_energy_area_nm2_at_2.5_kJ_mol": landscape["low_energy_area_nm2"],
        "initial_original_max_interchain_com_nm": float(original_separation[0]),
        "median_original_max_interchain_com_nm": float(np.median(original_separation[include])),
        "median_clustered_max_interchain_com_nm": float(np.median(clustered_separation[include])),
        "qc_status": status,
    }


def write_additional_source_data(results: list[dict], summaries: list[dict]) -> None:
    for result in results:
        rows = (
            {
                "time_ns": f"{time:.5f}",
                "core_backbone_rmsd_nm": f"{rmsd:.7f}",
                "core_radius_of_gyration_nm": f"{rg:.7f}",
                "original_max_interchain_com_distance_nm": f"{original_sep:.7f}",
                "clustered_max_interchain_com_distance_nm": f"{clustered_sep:.7f}",
                "included_after_burn_in": int(include),
            }
            for time, rmsd, rg, original_sep, clustered_sep, include in zip(
                result["time_ns"],
                result["rmsd_nm"],
                result["rg_nm"],
                result["original_max_interchain_com_nm"],
                result["clustered_max_interchain_com_nm"],
                result["include"],
            )
        )
        core.write_csv(
            SOURCE_DATA_DIR / f"{result['key']}_rmsd_rg_timeseries.csv",
            [
                "time_ns",
                "core_backbone_rmsd_nm",
                "core_radius_of_gyration_nm",
                "original_max_interchain_com_distance_nm",
                "clustered_max_interchain_com_distance_nm",
                "included_after_burn_in",
            ],
            rows,
        )

        landscape = result["individual_landscape"]
        grid_rows = []
        for i, rmsd in enumerate(landscape["rmsd_centers"]):
            for j, rg in enumerate(landscape["rg_centers"]):
                energy = landscape["free_energy"][i, j]
                grid_rows.append(
                    {
                        "core_backbone_rmsd_nm": f"{rmsd:.7f}",
                        "core_radius_of_gyration_nm": f"{rg:.7f}",
                        "smoothed_probability": f"{landscape['probability'][i, j]:.10g}",
                        "delta_G_kJ_mol": "" if not np.isfinite(energy) else f"{energy:.7f}",
                        "display_supported": int(landscape["supported"][i, j]),
                    }
                )
        core.write_csv(
            SOURCE_DATA_DIR / f"{result['key']}_individual_scale_free_energy_grid.csv",
            [
                "core_backbone_rmsd_nm",
                "core_radius_of_gyration_nm",
                "smoothed_probability",
                "delta_G_kJ_mol",
                "display_supported",
            ],
            grid_rows,
        )

    core.write_csv(
        SOURCE_DATA_DIR / "additional_summary_metrics.csv",
        list(summaries[0].keys()),
        summaries,
    )


def draw_additional_3d(results: list[dict], summaries_by_key: dict[str, dict]) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 3.55))
    fig.subplots_adjust(left=0.035, right=0.90, bottom=0.13, top=0.82, wspace=0.02)
    norm = Normalize(vmin=0, vmax=core.FREE_ENERGY_CAP_KJ_MOL)

    for index, result in enumerate(results, start=1):
        ax = fig.add_subplot(1, 3, index, projection="3d")
        landscape = result["individual_landscape"]
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
            s=18,
            marker="*",
            color="white",
            edgecolor="#202020",
            linewidth=0.45,
            depthshade=False,
            zorder=10,
        )
        rmsd_limits, rg_limits = result["individual_limits"]
        ax.set_xlim(*rmsd_limits)
        ax.set_ylim(*rg_limits)
        ax.set_zlim(-0.75, core.FREE_ENERGY_CAP_KJ_MOL)
        ax.set_zticks([0, 5, 10])
        ax.set_xlabel("Core RMSD (nm)", labelpad=2)
        ax.set_ylabel(r"Core $R_g$ (nm)", labelpad=2)
        ax.set_zlabel(r"$\Delta G$ (kJ mol$^{-1}$)", labelpad=2)
        ax.tick_params(pad=0, labelsize=5.5)
        ax.view_init(elev=27, azim=-128)
        ax.set_box_aspect((1.08, 1.00, 0.82))
        ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((1, 1, 1, 0))
            axis.pane.set_edgecolor("#d5d9df")
        ax.text2D(
            0.03, 0.96, chr(96 + index), transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top"
        )
        ax.text2D(
            0.13, 0.96, result["label"], transform=ax.transAxes,
            fontsize=8.8, fontweight="bold", color=COLORS[result["key"]], va="top"
        )
        summary = summaries_by_key[result["key"]]
        ax.text2D(
            0.13,
            0.89,
            f"PBC-clustered  |  {summary['analysis_end_ns']:.2f} ns\n"
            f"median RMSD {summary['median_rmsd_nm']:.3f} nm",
            transform=ax.transAxes,
            fontsize=5.3,
            color="#28734f" if summary["analysis_end_ns"] >= 199.9 else "#8a5a20",
            va="top",
        )

    colorbar_ax = fig.add_axes([0.925, 0.22, 0.016, 0.50])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=core.ENERGY_CMAP), cax=colorbar_ax
    )
    colorbar.set_label(r"Relative free energy, $\Delta G$ (kJ mol$^{-1}$)", labelpad=4)
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.outline.set_linewidth(0.6)

    fig.suptitle(
        r"PBC-corrected RMSD-$R_g$ landscapes of additional simulations",
        x=0.48,
        y=0.975,
        fontsize=10.5,
        fontweight="bold",
    )
    fig.text(
        0.48,
        0.905,
        "Three chains clustered by minimum image; common 30-residue core; 20 ns discarded",
        ha="center",
        fontsize=6.2,
        color="#46515d",
    )
    fig.text(
        0.48,
        0.025,
        "Independent panel scales reveal basin shape but must not be used to compare absolute spread.",
        ha="center",
        fontsize=5.8,
        color="#8a3030",
    )
    return fig


def draw_additional_2d(results: list[dict], shared: bool) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.82), sharex=shared, sharey=shared)
    fig.subplots_adjust(left=0.075, right=0.89, bottom=0.19, top=0.79, wspace=0.16)
    norm = Normalize(vmin=0, vmax=core.FREE_ENERGY_CAP_KJ_MOL)
    levels = np.linspace(0, core.FREE_ENERGY_CAP_KJ_MOL, 11)

    for index, (ax, result) in enumerate(zip(axes, results), start=1):
        landscape = result["shared_landscape"] if shared else result["individual_landscape"]
        z = np.ma.masked_where(~landscape["supported"], landscape["display_energy"])
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
            levels=[2.5, 5, 7.5, 10],
            colors="#27313a",
            linewidths=0.35,
            alpha=0.58,
        )
        ax.scatter(
            landscape["peak_rmsd_nm"],
            landscape["peak_rg_nm"],
            marker="*",
            s=28,
            color="white",
            edgecolor="#202020",
            linewidth=0.45,
            zorder=5,
        )
        if not shared:
            rmsd_limits, rg_limits = result["individual_limits"]
            ax.set_xlim(*rmsd_limits)
            ax.set_ylim(*rg_limits)
        ax.set_facecolor("#f5f6f7")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.14, 1.05, chr(96 + index), transform=ax.transAxes, fontsize=9, fontweight="bold")
        ax.set_title(
            result["label"], loc="left", color=COLORS[result["key"]],
            fontsize=8.8, fontweight="bold", pad=3
        )
        ax.text(
            0.99,
            1.03,
            "PBC-clustered",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.2,
            color="#28734f",
        )
        ax.set_xlabel("Core RMSD (nm)")
        if index == 1:
            ax.set_ylabel(r"Core $R_g$ (nm)")

    colorbar_ax = fig.add_axes([0.915, 0.22, 0.015, 0.50])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=core.ENERGY_CMAP), cax=colorbar_ax
    )
    colorbar.set_label(r"$\Delta G$ (kJ mol$^{-1}$)", labelpad=3)
    colorbar.set_ticks([0, 2.5, 5, 7.5, 10, 12.5])
    colorbar.outline.set_linewidth(0.6)

    scale_text = "shared RMSD/Rg axes" if shared else "independently scaled RMSD/Rg axes"
    fig.suptitle(
        r"PBC-corrected RMSD-$R_g$ landscapes",
        x=0.48,
        y=0.97,
        fontsize=10.2,
        fontweight="bold",
    )
    fig.text(
        0.48,
        0.875,
        f"Minimum-image chain clustering; 20 ns discarded; {scale_text}",
        ha="center",
        fontsize=6.2,
        color="#46515d",
    )
    if not shared:
        fig.text(
            0.48,
            0.035,
            "Use this view to inspect basin shape; use the shared-axis view for cross-system displacement.",
            ha="center",
            fontsize=5.7,
            color="#8a3030",
        )
    return fig


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if len(values) < window:
        return values.copy()
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(values, kernel, mode="same")


def draw_qc_timeseries(results: list[dict]) -> plt.Figure:
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 5.65), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.80, bottom=0.10, top=0.89, hspace=0.18)
    window = 501  # approximately 5 ns for 10-ps sampling

    for result in results:
        color = COLORS[result["key"]]
        time = result["time_ns"]
        for ax, values in zip(axes, (result["rmsd_nm"], result["rg_nm"])):
            ax.plot(time[::10], values[::10], color=color, alpha=0.13, linewidth=0.45)
            smoothed = moving_average(values, window)
            valid = np.arange(len(values)) >= window // 2
            valid &= np.arange(len(values)) < len(values) - window // 2
            ax.plot(time[valid], smoothed[valid], color=color, linewidth=1.35, label=result["label"])
        original_sep = result["original_max_interchain_com_nm"]
        clustered_sep = result["clustered_max_interchain_com_nm"]
        axes[2].plot(
            time[::10],
            original_sep[::10],
            color=color,
            linewidth=0.7,
            alpha=0.45,
            linestyle="--",
        )
        clustered_smooth = moving_average(clustered_sep, window)
        valid = np.arange(len(clustered_sep)) >= window // 2
        valid &= np.arange(len(clustered_sep)) < len(clustered_sep) - window // 2
        axes[2].plot(
            time[valid],
            clustered_smooth[valid],
            color=color,
            linewidth=1.35,
        )

    for ax in axes:
        ax.axvspan(0, core.BURN_IN_NS, color="#d9dde2", alpha=0.55, linewidth=0)
        ax.axvline(core.BURN_IN_NS, color="#69727c", linewidth=0.7, linestyle="--")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#dfe3e7", linewidth=0.45)

    axes[0].set_ylabel("Core backbone RMSD (nm)")
    axes[1].set_ylabel(r"Core $R_g$ (nm)")
    axes[2].set_ylabel("Maximum interchain\nCOM distance (nm)")
    axes[2].set_xlabel("Simulation time (ns)")
    axes[2].set_xlim(0, 200)
    axes[2].set_ylim(0, 12.2)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    axes[0].text(
        0.015, 0.92, "a", transform=axes[0].transAxes, fontsize=9, fontweight="bold"
    )
    axes[1].text(
        0.015, 0.92, "b", transform=axes[1].transAxes, fontsize=9, fontweight="bold"
    )
    axes[2].text(
        0.015, 0.90, "c", transform=axes[2].transAxes, fontsize=9, fontweight="bold"
    )
    axes[0].text(
        0.02,
        0.07,
        "excluded",
        transform=axes[0].transAxes,
        fontsize=5.6,
        color="#626a73",
    )
    p10_1_end = next(item["time_ns"][-1] for item in results if item["key"] == "P10_1")
    for ax in axes:
        ax.axvline(p10_1_end, color=COLORS["P10_1"], linewidth=0.75, linestyle=":")
    axes[2].text(
        p10_1_end - 1.5,
        0.55,
        f"P10_1 ends\n{p10_1_end:.2f} ns",
        ha="right",
        va="bottom",
        fontsize=5.6,
        color=COLORS["P10_1"],
    )
    axes[2].text(
        0.99,
        0.94,
        "dashed: stored fit.xtc\nsolid: minimum-image clustered",
        transform=axes[2].transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
        color="#46515d",
    )
    fig.suptitle(
        "Quality-control trajectories for additional simulations",
        x=0.44,
        y=0.97,
        fontsize=10.5,
        fontweight="bold",
    )
    fig.text(
        0.44,
        0.92,
        "RMSD and Rg use PBC-clustered coordinates; panel c documents the minimum-image correction",
        ha="center",
        fontsize=6.2,
        color="#46515d",
    )
    return fig


def write_documentation(results: list[dict], summaries: list[dict]) -> None:
    manifest = {
        "analysis": "Additional RMSD-Rg projected free-energy landscapes",
        "relationship_to_primary_figure": "Quality-control/additional simulations; not pooled into primary ranking",
        "settings": {
            "common_core_residues_per_chain": core.CORE_RESIDUES_PER_CHAIN,
            "periodic_boundary_correction": "minimum-image clustering of chain 2 and chain 3 around chain 1 for every frame",
            "burn_in_ns": core.BURN_IN_NS,
            "temperature_K": core.TEMPERATURE_K,
            "bins_per_axis": core.N_BINS,
            "gaussian_smoothing_sigma_bins": core.SMOOTH_SIGMA_BINS,
            "free_energy_cap_kJ_mol": core.FREE_ENERGY_CAP_KJ_MOL,
        },
        "systems": [
            {
                "name": item["label"],
                "source_folder": item["source_folder_name"],
                "frames": item["frames"],
                "last_time_ns": float(item["time_ns"][-1]),
            }
            for item in results
        ],
    }
    (OUTPUT_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    notes = f"""# Additional MD dataset audit and interpretation

Three additional processed trajectories were discovered: `mC`, `P10_0`, and `P10_1`.

- `mC` contains a complete 200-ns processed trajectory, but its original production `md.xtc` is no longer present in the folder. The retained `fit.xtc` is readable and was analyzed.
- `P10_0` contains a complete 200-ns processed trajectory.
- `P10_1` contains {results[2]['frames']:,} processed frames and ends at {results[2]['time_ns'][-1]:.2f} ns rather than 200 ns.

All three systems contain 30 residues per chain, so no length trimming was required beyond selecting the same 30-residue core definition used in the primary analysis. The first {core.BURN_IN_NS:.0f} ns were excluded.

The stored `fit.xtc` trajectories for `mC` and `P10_1` placed one chain approximately 11 nm from the other two, although the input topology placed all three chains within approximately 0.6 nm. This is a periodic-boundary reconstruction artifact. Before calculating RMSD and Rg, every frame was therefore corrected by translating chains 2 and 3 to the nearest periodic image around chain 1 using the triclinic unit-cell vectors. The correction reduced the maximum interchain center-of-mass distance to the expected clustered range without changing intrachain coordinates.

The independently scaled 3D and 2D landscapes are intended for inspecting each basin. The shared-axis 2D figure is the valid view for comparing overall conformational displacement. `P10_1` remains incomplete at {results[2]['time_ns'][-1]:.2f} ns and should be labelled as a truncated trajectory. The corrected additional simulations are not pooled into the primary CC/mC-ok/NC/NC3 figure because their construct roles have not been assigned in the primary experimental comparison.
"""
    (OUTPUT_DIR / "dataset_audit_and_interpretation.md").write_text(notes, encoding="utf-8")

    legend = f"""Supplementary Figure X | Projected free-energy landscapes and quality-control traces for additional simulations.

a-c, RMSD-Rg projected free-energy landscapes for the original mC run, P10_0, and P10_1. Before analysis, the three protein chains were reconstructed into the same periodic image by minimum-image clustering in every frame. The first {core.BURN_IN_NS:.0f} ns were excluded. RMSD and mass-weighted Rg were calculated for the same 30-residue-per-chain core used in the primary analysis. Relative free energy was calculated as DeltaG(x,y) = -RT ln[P(x,y)/Pmax] at {core.TEMPERATURE_K:.0f} K. The 3D and independently scaled 2D panels show basin shape, whereas the shared-axis panel enables direct comparison using identical axes. P10_1 ends at {results[2]['time_ns'][-1]:.2f} ns. One trajectory was available for each system.
"""
    (OUTPUT_DIR / "figure_legend.txt").write_text(legend, encoding="utf-8")


def main() -> None:
    configure_output()
    results = []
    for spec in ADDITIONAL_SYSTEMS:
        print(f"Reading and analyzing {spec['label']}...", flush=True)
        result = extract_clustered_timeseries(spec)
        individual_limits = core.common_limits([result])
        individual_landscape = core.build_landscape(result, *individual_limits)
        result["individual_limits"] = individual_limits
        result["individual_landscape"] = individual_landscape
        results.append(result)

    shared_limits = core.common_limits(results)
    for result in results:
        result["shared_landscape"] = core.build_landscape(result, *shared_limits)
        result["landscape"] = result["individual_landscape"]

    summaries = [summarize_additional(result) for result in results]
    summaries_by_key = {
        result["key"]: summary for result, summary in zip(results, summaries)
    }
    write_additional_source_data(results, summaries)
    write_documentation(results, summaries)

    print("Rendering additional 3D landscapes...", flush=True)
    core.save_figure(
        draw_additional_3d(results, summaries_by_key),
        OUTPUT_DIR / "additional_simulations_fel_3d_individual_axes",
    )
    print("Rendering additional 2D landscapes...", flush=True)
    core.save_figure(
        draw_additional_2d(results, shared=False),
        OUTPUT_DIR / "additional_simulations_fel_2d_individual_axes",
    )
    core.save_figure(
        draw_additional_2d(results, shared=True),
        OUTPUT_DIR / "additional_simulations_fel_2d_shared_axes",
    )
    print("Rendering trajectory quality-control figure...", flush=True)
    core.save_figure(
        draw_qc_timeseries(results),
        OUTPUT_DIR / "additional_simulations_rmsd_rg_qc",
    )
    print(f"Done: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
