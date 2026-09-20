#!/usr/bin/env python3
"""Plot a six-metric Cys-placement perturbation landscape."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("COLLAGEN_CYS_SCREENING_DATA", REPO_ROOT / "data" / "cys_screening"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "cys_screening"

SEQUENCE_LENGTH = 183
POSITION_XLIM = (0.5, SEQUENCE_LENGTH + 0.5)
SELECTED_SITES = {
    32: "P32C",
    54: "A54C",
    60: "Q60C",
    72: "A72C",
    90: "P90C",
    105: "V105C",
    117: "A117C",
    129: "A129C",
    150: "Q150C",
    167: "P167C",
    183: "P183C",
}

COLORS = {
    "ink": "#20242A",
    "muted": "#68727D",
    "grid": "#D8DDE2",
    "blue": "#276D9C",
    "blue_soft": "#DDECF5",
    "gold": "#C38122",
    "gold_soft": "#F5E6C5",
    "green": "#238451",
    "green_soft": "#DDEEE3",
    "purple": "#76559A",
    "purple_soft": "#E8E0EF",
    "red": "#B93632",
    "candidate": "#72A9C7",
}

REGIONS = [
    (1, 30, "N-terminal flank", COLORS["blue_soft"]),
    (31, 150, "Functional region", COLORS["gold_soft"]),
    (151, 183, "C-terminal flank", COLORS["green_soft"]),
]

RAW_METRICS = [
    "mutation_context_ca_rmsd_A",
    "gly_triangle_max_abs_delta_A2",
    "delta_local_interchain_hbond_bb_bb_energy",
    "delta_total_constraint_energy",
]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def burden_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        "perturbation_burden",
        ["#2C7FB8", "#7FCDBB", "#F4D66D", "#EF8A4C", "#C83E3A"],
    )


def positive(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").clip(lower=0)


def bounded(series: pd.Series, scale: float) -> pd.Series:
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"Invalid normalization scale: {scale}")
    return (positive(series) / scale).clip(0, 1)


def load_and_calibrate() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    replicates = pd.read_csv(DATA_DIR / "screening_replicates.csv")
    candidates = pd.read_csv(DATA_DIR / "candidate_library.csv")

    required = {
        "Job_Type",
        "Center_Candidate",
        "Center_Full_Position",
        *RAW_METRICS,
    }
    missing = required - set(replicates.columns)
    if missing:
        raise ValueError(f"Missing replicate columns: {sorted(missing)}")

    single_replicates = replicates.loc[replicates["Job_Type"].eq("single")].copy()
    grouped = single_replicates.groupby(
        ["Center_Candidate", "Center_Full_Position"],
        as_index=False,
    )
    medians = grouped[RAW_METRICS].median()
    replicate_counts = grouped.size().rename(columns={"size": "Replicate_N"})
    data = medians.merge(
        replicate_counts,
        on=["Center_Candidate", "Center_Full_Position"],
        validate="one_to_one",
    )
    data = data.sort_values("Center_Full_Position").reset_index(drop=True)

    if len(data) != 118:
        raise ValueError(f"Expected 118 single-site candidates, found {len(data)}")
    if not data["Replicate_N"].eq(3).all():
        bad = data.loc[~data["Replicate_N"].eq(3), "Center_Candidate"].tolist()
        raise ValueError(f"Expected three paired relax runs for every candidate: {bad}")
    if data[RAW_METRICS].isna().any().any():
        raise ValueError("One or more raw metric medians are missing")

    scales = {
        "Local_backbone_RMSD_fail_A": 1.5,
        "Triple_chain_geometry_fail_A2": 2.5,
        "HBond_energy_loss_positive_Q95_REU": float(
            positive(data["delta_local_interchain_hbond_bb_bb_energy"]).quantile(0.95)
        ),
        "Relax_constraint_positive_Q95_REU": float(
            positive(data["delta_total_constraint_energy"]).quantile(0.95)
        ),
    }

    data["Burden_Local_backbone_RMSD"] = bounded(
        data["mutation_context_ca_rmsd_A"],
        scales["Local_backbone_RMSD_fail_A"],
    )
    data["Burden_Triple_chain_geometry"] = bounded(
        data["gly_triangle_max_abs_delta_A2"],
        scales["Triple_chain_geometry_fail_A2"],
    )
    data["Burden_Backbone_HBond_energy_loss"] = bounded(
        data["delta_local_interchain_hbond_bb_bb_energy"],
        scales["HBond_energy_loss_positive_Q95_REU"],
    )
    data["Burden_Relax_constraint"] = bounded(
        data["delta_total_constraint_energy"],
        scales["Relax_constraint_positive_Q95_REU"],
    )

    # Give every displayed component the same contribution to the composite.
    data["Composite_perturbation"] = (
        0.25 * data["Burden_Local_backbone_RMSD"]
        + 0.25 * data["Burden_Triple_chain_geometry"]
        + 0.25 * data["Burden_Backbone_HBond_energy_loss"]
        + 0.25 * data["Burden_Relax_constraint"]
    )
    data["Experimental_site"] = data["Center_Full_Position"].isin(SELECTED_SITES)
    return data, candidates, scales


def position_edges(positions: np.ndarray) -> np.ndarray:
    midpoints = (positions[:-1] + positions[1:]) / 2
    return np.concatenate(([0.5], midpoints, [SEQUENCE_LENGTH + 0.5]))


def draw_flow(ax: plt.Axes, data: pd.DataFrame, candidates: pd.DataFrame) -> None:
    ax.set_axis_off()
    allowed = int(candidates["Allowed"].astype(bool).sum())
    boxes = [
        (0.02, 0.16, 0.18, "183 residues", "Complete collagen domain", COLORS["blue_soft"], COLORS["blue"]),
        (0.275, 0.16, 0.19, f"{allowed} candidates", "Allowed X/Y-register sites", COLORS["gold_soft"], COLORS["gold"]),
        (0.535, 0.16, 0.19, "Four-component landscape", "Four components + composite", COLORS["green_soft"], COLORS["green"]),
        (0.80, 0.08, 0.18, "11 experimental sites", "Topology and density series", COLORS["purple_soft"], COLORS["purple"]),
    ]
    for x, y, width, title, subtitle, face, edge in boxes:
        height = 0.70 if x < 0.80 else 0.86
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            transform=ax.transAxes,
            boxstyle="round,pad=0.010,rounding_size=0.015",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.2,
        )
        ax.add_patch(patch)
        ax.text(
            x + width / 2,
            y + height * 0.61,
            title,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10.5,
            fontweight="bold",
            color=edge,
        )
        ax.text(
            x + width / 2,
            y + height * 0.28,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7,
            color=COLORS["ink"],
        )
    for left, right in ((0.205, 0.27), (0.47, 0.53), (0.73, 0.795)):
        ax.add_patch(
            FancyArrowPatch(
                (left, 0.51),
                (right, 0.51),
                transform=ax.transAxes,
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.3,
                color=COLORS["muted"],
            )
        )
    ax.text(
        0.5,
        1.18,
        "Four-component structural perturbation landscape for rational Cys-site selection",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        0.5,
        1.06,
        "All metrics are derived from the median of three paired WT/mutant relax runs",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.4,
        color=COLORS["blue"],
        fontweight="bold",
    )


def add_region_background(ax: plt.Axes, alpha: float) -> None:
    for start, end, _, color in REGIONS:
        ax.axvspan(start - 0.5, end + 0.5, color=color, alpha=alpha, lw=0)


def draw_candidate_track(ax: plt.Axes, data: pd.DataFrame, candidates: pd.DataFrame) -> None:
    # This limit is shared with the heatmap so a star, a dashed guide and the
    # center of the corresponding heatmap column always use the same x value.
    ax.set_xlim(*POSITION_XLIM)
    ax.set_ylim(0, 1)
    add_region_background(ax, alpha=0.72)
    for start, end, label, _ in REGIONS:
        ax.text(
            (start + end) / 2,
            0.84,
            label,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )

    allowed = candidates.loc[candidates["Allowed"].astype(bool)]
    ax.scatter(
        allowed["Full_Position"],
        np.full(len(allowed), 0.29),
        s=11,
        color=COLORS["candidate"],
        alpha=0.65,
        edgecolor="none",
        zorder=2,
    )
    selected = data.loc[data["Experimental_site"]]
    ax.scatter(
        selected["Center_Full_Position"],
        np.full(len(selected), 0.29),
        marker="*",
        s=92,
        facecolor=COLORS["red"],
        edgecolor="white",
        linewidth=0.6,
        zorder=4,
    )
    for row in selected.itertuples():
        position = int(row.Center_Full_Position)
        offset = (-4, 11) if position >= 178 else (0, 11)
        alignment = "right" if position >= 178 else "center"
        ax.annotate(
            row.Center_Candidate,
            xy=(position, 0.29),
            xytext=offset,
            textcoords="offset points",
            rotation=61,
            ha=alignment,
            va="bottom",
            fontsize=6.1,
            color=COLORS["red"],
            annotation_clip=False,
        )

    ax.text(1, 0.06, "Allowed candidates", ha="left", va="bottom", fontsize=7, color=COLORS["muted"])
    ax.text(
        183,
        0.06,
        "Red stars: experimental sites",
        ha="right",
        va="bottom",
        fontsize=7,
        color=COLORS["red"],
    )
    ax.set_xticks([1, 30, 60, 90, 120, 150, 183])
    ax.set_yticks([])
    ax.tick_params(axis="x", length=3)
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(COLORS["muted"])


def draw_heatmap(ax: plt.Axes, data: pd.DataFrame, show_group_strip: bool = False) -> pd.DataFrame:
    definitions = [
        ("Composite perturbation", "Composite_perturbation", "summary"),
        ("Local backbone RMSD", "Burden_Local_backbone_RMSD", "structure"),
        ("Triple-chain geometry", "Burden_Triple_chain_geometry", "structure"),
        (
            "Backbone H-bond energy-loss burden",
            "Burden_Backbone_HBond_energy_loss",
            "qc",
        ),
        ("Relax constraint compatibility", "Burden_Relax_constraint", "qc"),
    ]
    matrix = np.vstack([data[column].to_numpy(float) for _, column, _ in definitions])
    positions = data["Center_Full_Position"].to_numpy(float)
    mesh = ax.pcolormesh(
        position_edges(positions),
        np.arange(len(definitions) + 1),
        matrix,
        cmap=burden_cmap(),
        vmin=0,
        vmax=1,
        shading="flat",
    )
    ax.invert_yaxis()
    ax.set_xlim(*POSITION_XLIM)
    ax.set_xticks([1, 30, 60, 90, 120, 150, 183])
    ax.set_xlabel("Residue position in the 183-aa collagen-like domain")
    ax.set_yticks(np.arange(len(definitions)) + 0.5)
    ax.set_yticklabels([label for label, _, _ in definitions])
    ax.tick_params(axis="y", length=0, pad=6)

    for position in SELECTED_SITES:
        ax.axvline(position, color=COLORS["red"], lw=0.65, alpha=0.62, linestyle="--")
    for boundary in (30.5, 150.5):
        ax.axvline(boundary, color=COLORS["ink"], lw=0.9, alpha=0.5)
    for y in range(1, len(definitions)):
        ax.axhline(y, color="white", lw=0.8)

    if show_group_strip:
        groups = [
            (0, 1, "Summary", COLORS["purple"]),
            (1, 3, "Structure", COLORS["blue"]),
            (3, 5, "QC", COLORS["green"]),
        ]
        for start, end, label, color in groups:
            y0 = 1 - end / len(definitions)
            height = (end - start) / len(definitions)
            ax.add_patch(
                Rectangle(
                    (-0.175, y0),
                    0.020,
                    height,
                    transform=ax.transAxes,
                    facecolor=color,
                    edgecolor="none",
                    clip_on=False,
                )
            )
            ax.text(
                -0.190,
                y0 + height / 2,
                label,
                transform=ax.transAxes,
                rotation=90,
                ha="center",
                va="center",
                fontsize=6.2,
                fontweight="bold",
                color=color,
            )

    colorbar = ax.figure.colorbar(mesh, ax=ax, fraction=0.019, pad=0.016)
    colorbar.set_ticks([0, 0.5, 1])
    colorbar.set_ticklabels(["Low", "Moderate", "High"])
    colorbar.set_label("Calibrated perturbation burden", rotation=270, labelpad=12)
    colorbar.outline.set_linewidth(0.7)

    note = (
        "Composite = 25% RMSD + 25% triple-chain geometry "
        "+ 25% H-bond energy loss + 25% relax-constraint penalty. "
        "Blue denotes preserved structure; red denotes greater perturbation."
    )
    ax.text(
        0,
        -0.29,
        note,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        color=COLORS["muted"],
    )
    for spine in ax.spines.values():
        spine.set_visible(False)

    source = data[
        [
            "Center_Candidate",
            "Center_Full_Position",
            "Replicate_N",
            *RAW_METRICS,
            *[column for _, column, _ in definitions],
            "Experimental_site",
        ]
    ].copy()
    return source


def metric_definitions(scales: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Displayed_metric": "Composite perturbation",
                "Raw_field": "derived",
                "Aggregation": "weighted sum of four burdens",
                "Normalization": "already bounded 0-1",
                "Composite_weight": 1.0,
            },
            {
                "Displayed_metric": "Local backbone RMSD",
                "Raw_field": "mutation_context_ca_rmsd_A",
                "Aggregation": "median of 3 paired relax runs",
                "Normalization": f"clip(value/{scales['Local_backbone_RMSD_fail_A']:.4f} A, 0, 1)",
                "Composite_weight": 0.25,
            },
            {
                "Displayed_metric": "Triple-chain geometry",
                "Raw_field": "gly_triangle_max_abs_delta_A2",
                "Aggregation": "median of 3 paired relax runs",
                "Normalization": f"clip(value/{scales['Triple_chain_geometry_fail_A2']:.4f} A2, 0, 1)",
                "Composite_weight": 0.25,
            },
            {
                "Displayed_metric": "Backbone H-bond energy-loss burden",
                "Raw_field": "delta_local_interchain_hbond_bb_bb_energy",
                "Aggregation": "median of 3 paired relax runs",
                "Normalization": f"clip(max(value,0)/{scales['HBond_energy_loss_positive_Q95_REU']:.6f} REU [positive Q95], 0, 1)",
                "Composite_weight": 0.25,
            },
            {
                "Displayed_metric": "Relax constraint compatibility",
                "Raw_field": "delta_total_constraint_energy",
                "Aggregation": "median of 3 paired relax runs",
                "Normalization": f"clip(max(value,0)/{scales['Relax_constraint_positive_Q95_REU']:.6f} REU [positive Q95], 0, 1)",
                "Composite_weight": 0.25,
            },
        ]
    )


def write_legend(scales: dict[str, float]) -> None:
    legend = f"""Figure | Four-component structural perturbation landscape for rational Cys-site selection.

Of 183 residues, 118 X/Y-register positions were evaluated. Red stars mark the 11 sites used in the experimental topology-control and internal-Cys density series. Every heatmap cell is calculated from the median of three independent paired WT/mutant relax runs.

The four component rows are local backbone RMSD, triple-chain geometry, local interchain backbone-backbone H-bond energy-loss burden and relax-constraint compatibility. Triple-chain geometry is represented by the maximum absolute change in the cross-sectional Gly triangle area, providing a direct local readout of three-chain packing geometry.

Composite perturbation is the arithmetic mean of the four component burdens: 0.25 x RMSD burden + 0.25 x triple-chain geometry burden + 0.25 x backbone H-bond energy-loss burden + 0.25 x relax-constraint burden. Mutation total energy, steric repulsion and the overlapping Gly-register readout are not included. RMSD and triple-chain geometry use fixed scales of 1.5 A and 2.5 A2. Positive H-bond energy loss and positive relax-constraint change use candidate-wide Q95 values of {scales['HBond_energy_loss_positive_Q95_REU']:.4f} and {scales['Relax_constraint_positive_Q95_REU']:.4f} REU, respectively; non-positive changes receive zero burden.

Lower values indicate greater structural preservation. Composite perturbation is a transparent screening summary, not a PyRosetta energy term, free energy, disulfide-formation probability or MD persistence measurement.
"""
    (OUTPUT_DIR / "figure_legend.txt").write_text(legend, encoding="utf-8")


def build_figure() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    source_dir = OUTPUT_DIR / "source_data"
    source_dir.mkdir(exist_ok=True)

    data, candidates, scales = load_and_calibrate()

    # Minimal two-part panel: the sequence-position track and its aligned
    # five-row screening landscape.  No workflow or side panels are retained.
    fig = plt.figure(figsize=(12.4, 4.65))
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.88, 2.58],
        hspace=0.075,
        left=0.198,
        right=0.935,
        top=0.945,
        bottom=0.145,
    )
    track_ax = fig.add_subplot(grid[0])
    heatmap_ax = fig.add_subplot(grid[1], sharex=track_ax)

    draw_candidate_track(track_ax, data, candidates)
    track_ax.tick_params(axis="x", labelbottom=False)
    source = draw_heatmap(heatmap_ax, data, show_group_strip=False)

    # Assert the visual contract in code: both panels have identical residue axes.
    if tuple(track_ax.get_xlim()) != tuple(heatmap_ax.get_xlim()):
        raise RuntimeError("Candidate track and heatmap residue axes are not aligned")

    base = OUTPUT_DIR / "four_component_screening_landscape_compact_aligned"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)

    source.to_csv(source_dir / "four_component_candidate_landscape.csv", index=False)
    metric_definitions(scales).to_csv(source_dir / "metric_definitions.csv", index=False)
    pd.DataFrame(
        [{"Calibration_parameter": key, "Value": value} for key, value in scales.items()]
    ).to_csv(source_dir / "normalization_scales.csv", index=False)
    write_legend(scales)


if __name__ == "__main__":
    build_figure()
