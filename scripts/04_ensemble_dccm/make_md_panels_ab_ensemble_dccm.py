#!/usr/bin/env python3
"""Generate Panel A (ensemble/RMSF putty) and Panel B (DCCM) for four MD runs.

The centered common 30-residue segment in each chain is used solely for
PBC-aware alignment. RMSF, rendered structures and DCCM retain every residue
present in each simulation segment, including terminal Cys sites. The output
is descriptive for one trajectory per construct; it does not treat time frames
as independent repeats.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "ensemble_dccm"
SOURCE_DIR = OUTPUT_DIR / "source_data"
STRUCTURE_DIR = OUTPUT_DIR / "structures"
RENDER_DIR = OUTPUT_DIR / "pymol_renders"
SOURCE_ROOT = Path(os.environ.get("COLLAGEN_DATA_ROOT", REPO_ROOT / "data"))
CACHE_DIR = Path(os.environ.get("COLLAGEN_CACHE_DIR", Path(tempfile.gettempdir()) / "collagen_cys_md_cache"))
PYMOL = Path(os.environ.get("PYMOL_BIN", shutil.which("pymol") or "pymol"))

LOCAL_DEPS = HERE / ".mdanalysis-extracted"
if LOCAL_DEPS.exists():
    sys.path = [str(LOCAL_DEPS), *[path for path in sys.path if "site-packages" not in path.lower()]]

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors


SYSTEMS = [
    {"key": "CC", "label": "CC", "folder_prefix": "CC", "color": "#7250a5"},
    {"key": "mC", "label": "mC ok", "folder": "mC ok", "color": "#2878b8"},
    {"key": "NC", "label": "NC", "folder": "NC", "color": "#23864f"},
    {"key": "P10_0", "label": "P10_0", "folder": "P10_0", "color": "#b74744"},
]
REFERENCE_KEY = "P10_0"
CORE_RESIDUES_PER_CHAIN = 30
BURN_IN_NS = 20.0
ENSEMBLE_FRAMES = 16
RMSF_CMAP = LinearSegmentedColormap.from_list(
    "rmsf_putty", ["#1f4e99", "#2f8fc0", "#66c2a5", "#f6d55c", "#d73027"], N=256
)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.0,
            "axes.labelsize": 7.0,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "axes.linewidth": 0.65,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def resolve_folder(spec: dict) -> Path:
    if "folder" in spec:
        path = SOURCE_ROOT / spec["folder"]
        if not path.is_dir():
            raise FileNotFoundError(path)
        return path
    matches = sorted(path for path in SOURCE_ROOT.iterdir() if path.is_dir() and path.name.startswith(spec["folder_prefix"]))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one folder beginning with {spec['folder_prefix']!r}; found {matches}")
    return matches[0]


def has_private_character(path: Path) -> bool:
    return any(0xE000 <= ord(char) <= 0xF8FF for char in str(path))


def prepare_input(spec: dict) -> tuple[Path, Path, Path]:
    source_folder = resolve_folder(spec)
    topology = source_folder / "1_processed.gro"
    trajectory = source_folder / "fit.xtc"
    if not topology.is_file() or not trajectory.is_file():
        raise FileNotFoundError(f"Expected 1_processed.gro and fit.xtc in {source_folder}")
    if has_private_character(source_folder):
        cached = CACHE_DIR / spec["key"]
        cached.mkdir(parents=True, exist_ok=True)
        cached_topology = cached / "topology.gro"
        cached_trajectory = cached / "trajectory.xtc"
        for source, destination in ((topology, cached_topology), (trajectory, cached_trajectory)):
            if not destination.exists() or destination.stat().st_size != source.stat().st_size:
                shutil.copy2(source, destination)
        topology, trajectory = cached_topology, cached_trajectory
    return source_folder, topology, trajectory


def select_cys_aware_segments(universe: mda.Universe) -> dict:
    """Keep every simulated residue while using the common 30-aa segment only as an alignment anchor."""
    residues = universe.residues
    if len(residues) % 3:
        raise ValueError("Residues cannot be divided into three chains")
    chain_length = len(residues) // 3
    if chain_length < CORE_RESIDUES_PER_CHAIN:
        raise ValueError("Chain is shorter than the requested common core")
    trim_left = (chain_length - CORE_RESIDUES_PER_CHAIN) // 2
    all_residues = [
        residues[chain * chain_length : (chain + 1) * chain_length]
        for chain in range(3)
    ]
    anchor_residues = [
        residues[chain * chain_length + trim_left : chain * chain_length + trim_left + CORE_RESIDUES_PER_CHAIN]
        for chain in range(3)
    ]
    anchor_atoms = [group.atoms for group in anchor_residues]
    anchor_ca_atoms = [group.atoms.select_atoms("name CA") for group in anchor_residues]
    all_atoms = [group.atoms for group in all_residues]
    all_ca_atoms = [group.atoms.select_atoms("name CA") for group in all_residues]
    if any(len(group) != CORE_RESIDUES_PER_CHAIN for group in anchor_ca_atoms):
        raise ValueError("Every core residue must contain exactly one CA atom")
    if any(len(group) != chain_length for group in all_ca_atoms):
        raise ValueError("Every simulated residue must contain exactly one CA atom")
    cys_positions = [
        local_index
        for local_index, residue in enumerate(all_residues[0], start=1)
        if str(residue.resname).upper() in {"CYS", "CYX"}
    ]
    return {
        "chain_length": chain_length,
        "all_residues": all_residues,
        "anchor_atoms": anchor_atoms,
        "anchor_ca_atoms": anchor_ca_atoms,
        "all_atoms": all_atoms,
        "all_ca_atoms": all_ca_atoms,
        "cys_positions": cys_positions,
    }


def pbc_clustered_coordinates(universe: mda.Universe, anchor_atoms: list, anchor_ca_atoms: list, all_atoms: list, all_ca_atoms: list) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    box = universe.trajectory.ts.dimensions
    if box is None or not np.all(np.isfinite(box)):
        raise ValueError("Trajectory frame lacks finite periodic-box dimensions")
    centers = np.asarray([atoms.center_of_mass() for atoms in anchor_atoms])
    shifts = [np.zeros(3, dtype=float)]
    for chain in (1, 2):
        vector = np.asarray(centers[chain] - centers[0], dtype=np.float32).reshape(1, 3)
        nearest = minimize_vectors(vector, box)[0].astype(float)
        shifts.append(centers[0] + nearest - centers[chain])
    anchor_ca = np.concatenate([atoms.positions.astype(float) + shift for atoms, shift in zip(anchor_ca_atoms, shifts)])
    full_atoms = np.concatenate([atoms.positions.astype(float) + shift for atoms, shift in zip(all_atoms, shifts)])
    full_ca = np.concatenate([atoms.positions.astype(float) + shift for atoms, shift in zip(all_ca_atoms, shifts)])
    return anchor_ca, full_atoms, full_ca


def kabsch_transform(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    mobile_centered = mobile - mobile_center
    reference_centered = reference - reference_center
    left, _, right_t = np.linalg.svd(mobile_centered.T @ reference_centered)
    if np.linalg.det(left @ right_t) < 0:
        left[:, -1] *= -1
    return mobile_center, reference_center, left @ right_t


def apply_transform(coordinates: np.ndarray, mobile_center: np.ndarray, reference_center: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    return (coordinates - mobile_center) @ rotation + reference_center


def build_atom_records(all_residues: list, chain_length: int) -> list[dict]:
    records = []
    serial = 1
    for chain_index, residues in enumerate(all_residues):
        for local_index, residue in enumerate(residues, start=1):
            for atom in residue.atoms:
                name = str(atom.name).strip()
                first = name[0].upper() if name else "C"
                element = first if first in {"C", "N", "O", "S", "H", "P"} else "C"
                records.append(
                    {
                        "serial": serial,
                        "name": name,
                        "resname": str(residue.resname)[:3],
                        "chain": "ABC"[chain_index],
                        "local_resid": local_index,
                        "source_resid": int(residue.resid),
                        "rmsf_index": chain_index * chain_length + local_index - 1,
                        "element": element,
                    }
                )
                serial += 1
    return records


def write_pdb(path: Path, models: list[np.ndarray], records: list[dict], rmsf: np.ndarray) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for model_index, coordinates in enumerate(models, start=1):
            if len(models) > 1:
                handle.write(f"MODEL     {model_index:4d}\n")
            for record, coordinate in zip(records, coordinates):
                bfactor = float(rmsf[record["rmsf_index"]])
                handle.write(
                    f"ATOM  {record['serial']:5d} {record['name']:<4.4s} {record['resname']:>3.3s} {record['chain']}{record['local_resid']:4d}    "
                    f"{coordinate[0]:8.3f}{coordinate[1]:8.3f}{coordinate[2]:8.3f}{1.00:6.2f}{bfactor:6.2f}          {record['element']:>2s}\n"
                )
            handle.write("TER\n")
            if len(models) > 1:
                handle.write("ENDMDL\n")
        handle.write("END\n")


def analyze_system(spec: dict) -> dict:
    source_folder, topology, trajectory = prepare_input(spec)
    universe = mda.Universe(str(topology), str(trajectory))
    segments = select_cys_aware_segments(universe)
    records = build_atom_records(segments["all_residues"], segments["chain_length"])
    frame_count = len(universe.trajectory)
    universe.trajectory[0]
    reference_anchor_ca, _, _ = pbc_clustered_coordinates(
        universe,
        segments["anchor_atoms"],
        segments["anchor_ca_atoms"],
        segments["all_atoms"],
        segments["all_ca_atoms"],
    )
    expected_times = float(universe.trajectory.ts.time / 1000.0) + np.arange(frame_count) * float(universe.trajectory.dt / 1000.0)
    post_indices = np.flatnonzero(expected_times >= BURN_IN_NS)
    if len(post_indices) < 1000:
        raise ValueError(f"Too few post-equilibration frames for {spec['label']}")
    sample_indices = np.unique(np.rint(np.linspace(post_indices[0], post_indices[-1], ENSEMBLE_FRAMES)).astype(int))
    sample_lookup = {int(index) for index in sample_indices}
    representative_index = int(post_indices[len(post_indices) // 2])
    n_ca = 3 * segments["chain_length"]
    dccm_positions = np.empty((len(post_indices), n_ca, 3), dtype=np.float32)
    mean = np.zeros((n_ca, 3), dtype=float)
    m2 = np.zeros((n_ca, 3), dtype=float)
    sample_models = []
    representative_model = None
    time_ns = np.empty(len(post_indices), dtype=float)
    post_cursor = 0

    for frame_index, timestep in enumerate(universe.trajectory):
        if frame_index < post_indices[0]:
            continue
        anchor_ca, full_atoms, full_ca = pbc_clustered_coordinates(
            universe,
            segments["anchor_atoms"],
            segments["anchor_ca_atoms"],
            segments["all_atoms"],
            segments["all_ca_atoms"],
        )
        mobile_center, reference_center, rotation = kabsch_transform(anchor_ca, reference_anchor_ca)
        aligned_ca = apply_transform(full_ca, mobile_center, reference_center, rotation)
        aligned_atoms = apply_transform(full_atoms, mobile_center, reference_center, rotation)
        if frame_index in sample_lookup:
            sample_models.append(aligned_atoms.astype(np.float32))
        if frame_index == representative_index:
            representative_model = aligned_atoms.astype(np.float32)
        dccm_positions[post_cursor] = aligned_ca.astype(np.float32)
        time_ns[post_cursor] = timestep.time / 1000.0
        count = post_cursor + 1
        delta = aligned_ca - mean
        mean += delta / count
        m2 += delta * (aligned_ca - mean)
        post_cursor += 1

    if representative_model is None:
        representative_model = sample_models[len(sample_models) // 2]
    if post_cursor != len(post_indices):
        raise RuntimeError("Unexpected number of post-equilibration frames")
    rmsf = np.sqrt(np.sum(m2 / max(post_cursor - 1, 1), axis=1))
    fluctuations = dccm_positions.astype(float) - dccm_positions.mean(axis=0)
    covariance = np.einsum("tip,tjp->ij", fluctuations, fluctuations) / len(fluctuations)
    denominator = np.sqrt(np.clip(np.diag(covariance), 1e-12, None))
    dccm = covariance / np.outer(denominator, denominator)
    dccm = np.clip(dccm, -1, 1)
    return {
        "spec": spec,
        "source_folder": source_folder.name,
        "records": records,
        "rmsf": rmsf,
        "dccm": dccm,
        "ensemble_models": sample_models,
        "representative_model": representative_model,
        "frames": int(post_cursor),
        "analysis_start_ns": float(time_ns[0]),
        "analysis_end_ns": float(time_ns[-1]),
        "chain_length": segments["chain_length"],
        "cys_positions": segments["cys_positions"],
    }


def write_structures(result: dict, rmsf_min: float, rmsf_max: float) -> tuple[Path, Path, Path]:
    key = result["spec"]["key"]
    putty_pdb = STRUCTURE_DIR / f"{key}_rmsf_putty.pdb"
    ensemble_pdb = STRUCTURE_DIR / f"{key}_ensemble_{len(result['ensemble_models'])}_frames.pdb"
    pml = RENDER_DIR / f"{key}_ensemble_putty.pml"
    png = RENDER_DIR / f"{key}_ensemble_putty.png"
    write_pdb(putty_pdb, [result["representative_model"]], result["records"], result["rmsf"])
    write_pdb(ensemble_pdb, result["ensemble_models"], result["records"], result["rmsf"])
    pml.write_text(
        "\n".join(
            [
                "reinitialize",
                f'load "structures/{putty_pdb.name}", putty',
                f'load "structures/{ensemble_pdb.name}", ensemble',
                "set all_states, on",
                "hide everything, all",
                "show cartoon, ensemble",
                "color gray85, ensemble",
                "set cartoon_transparency, 0.84, ensemble",
                "show cartoon, putty",
                "cartoon putty, putty",
                "set cartoon_putty_transform, 0, putty",
                "set cartoon_putty_scale_min, 0.40, putty",
                "set cartoon_putty_scale_max, 1.65, putty",
                f"spectrum b, blue_white_red, putty, minimum={rmsf_min:.4f}, maximum={rmsf_max:.4f}",
                "show spheres, putty and resn CYS",
                "color yellow, putty and resn CYS",
                "set sphere_scale, 0.34, putty and resn CYS",
                "set orthoscopic, on",
                "set ray_opaque_background, off",
                "set antialias, 2",
                "bg_color white",
                "orient putty",
                "zoom putty, 1.32",
                f"png pymol_renders/{png.name}, 1280, 980, 300, 1",
                "quit",
            ]
        ),
        encoding="ascii",
    )
    return putty_pdb, ensemble_pdb, pml


def run_pymol(pml: Path) -> None:
    if not PYMOL.is_file():
        raise FileNotFoundError(f"PyMOL executable was not found: {PYMOL}")
    completed = subprocess.run(
        [str(PYMOL), "-cq", str(pml)],
        capture_output=True,
        text=True,
        timeout=180,
        cwd=str(OUTPUT_DIR),
    )
    if completed.returncode:
        raise RuntimeError(f"PyMOL rendering failed for {pml.name}:\n{completed.stdout}\n{completed.stderr}")
    expected_png = RENDER_DIR / f"{pml.stem}.png"
    if not expected_png.is_file():
        raise RuntimeError(f"PyMOL did not create the expected image: {expected_png}")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=450, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTPUT_DIR / f"{stem}.svg", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTPUT_DIR / f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.04, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def trim_render_margin(image_path: Path) -> np.ndarray:
    """Remove transparent PyMOL canvas margin without altering molecular geometry."""
    image = plt.imread(image_path)
    if image.ndim == 3 and image.shape[-1] == 4:
        mask = image[..., 3] > 0.03
    else:
        mask = np.any(image[..., :3] < 0.985, axis=-1)
    rows, columns = np.where(mask)
    if not len(rows):
        return image
    pad_y, pad_x = 24, 36
    y0, y1 = max(rows.min() - pad_y, 0), min(rows.max() + pad_y + 1, image.shape[0])
    x0, x1 = max(columns.min() - pad_x, 0), min(columns.max() + pad_x + 1, image.shape[1])
    return image[y0:y1, x0:x1]


def draw_panel_a(results: list[dict], rmsf_min: float, rmsf_max: float) -> plt.Figure:
    fig, axes = plt.subplots(4, 1, figsize=(8.65, 7.20))
    fig.subplots_adjust(left=0.13, right=0.84, bottom=0.075, top=0.865, hspace=0.16)
    for letter, ax, result in zip("abcd", axes, results):
        image_path = RENDER_DIR / f"{result['spec']['key']}_ensemble_putty.png"
        ax.imshow(trim_render_margin(image_path))
        ax.axis("off")
        ax.text(-0.105, 0.93, letter, transform=ax.transAxes, fontsize=9, fontweight="bold")
        ax.text(-0.070, 0.93, result["spec"]["label"], transform=ax.transAxes, fontsize=9.5, color=result["spec"]["color"], fontweight="bold")
    colorbar_ax = fig.add_axes([0.875, 0.19, 0.018, 0.59])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(rmsf_min, rmsf_max), cmap=RMSF_CMAP), cax=colorbar_ax)
    colorbar.set_label(r"Core C$\alpha$ RMSF ($\AA$)", fontsize=6.7, labelpad=4)
    colorbar.ax.tick_params(labelsize=5.8, length=2)
    fig.suptitle("Panel A | Conformational ensembles and residue-level flexibility", x=0.47, y=0.972, fontsize=10.8, fontweight="bold")
    fig.text(0.47, 0.942, f"{ENSEMBLE_FRAMES} PBC-clustered, aligned frames from 20-200 ns are overlaid in gray; a representative structure is rendered as RMSF putty", ha="center", fontsize=6.2, color="#46515d")
    fig.text(0.47, 0.018, "RMSF color and tube-radius ranges are fixed across all constructs; yellow spheres denote every Cys residue in the simulated segment.", ha="center", fontsize=5.9, color="#46515d")
    return fig


def draw_dccm(ax: plt.Axes, result: dict, norm, cmap, letter: str) -> mpl.image.AxesImage:
    matrix = result["dccm"]
    chain_length = result["chain_length"]
    image = ax.imshow(matrix, cmap=cmap, norm=norm, interpolation="nearest", origin="upper", rasterized=True)
    for boundary in (chain_length - 0.5, 2 * chain_length - 0.5):
        ax.axvline(boundary, color="#202830", linewidth=0.65, linestyle="--")
        ax.axhline(boundary, color="#202830", linewidth=0.65, linestyle="--")
    for cys_position in result["cys_positions"]:
        for chain_index in range(3):
            coordinate = chain_index * chain_length + cys_position - 0.5
            ax.axvline(coordinate, color="#d4a20a", linewidth=0.55, alpha=0.95)
            ax.axhline(coordinate, color="#d4a20a", linewidth=0.55, alpha=0.95)
    midpoints = [chain_length / 2 - 0.5, 1.5 * chain_length - 0.5, 2.5 * chain_length - 0.5]
    ax.set_xticks(midpoints, ["A", "B", "C"])
    ax.set_yticks(midpoints, ["A", "B", "C"])
    ax.set_xlabel("Chain-residue block")
    ax.set_ylabel("Chain-residue block")
    cys_label = "no Cys" if not result["cys_positions"] else "Cys " + ", ".join(map(str, result["cys_positions"]))
    ax.set_title(f"{result['spec']['label']} | {chain_length} aa/chain; {cys_label}", loc="left", color=result["spec"]["color"], fontweight="bold", pad=3)
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes, fontsize=9, fontweight="bold")
    return image


def draw_panel_b_dccm(results: list[dict]) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(7.45, 6.1))
    fig.subplots_adjust(left=0.10, right=0.86, bottom=0.10, top=0.86, hspace=0.33, wspace=0.28)
    norm = Normalize(vmin=-1, vmax=1)
    for letter, ax, result in zip("abcd", axes.flat, results):
        draw_dccm(ax, result, norm, "RdBu_r", letter)
    colorbar_ax = fig.add_axes([0.89, 0.22, 0.018, 0.53])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=colorbar_ax)
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    colorbar.set_label("Dynamic cross-correlation", fontsize=6.7, labelpad=4)
    colorbar.ax.tick_params(labelsize=5.8, length=2)
    fig.suptitle("Panel B | Dynamic cross-correlation matrices", x=0.46, y=0.972, fontsize=11.0, fontweight="bold")
    fig.text(0.46, 0.94, "All simulated C-alpha residues are analyzed after PBC clustering, common-core alignment and exclusion of the first 20 ns", ha="center", fontsize=6.1, color="#46515d")
    fig.text(0.46, 0.015, "Red: correlated motion; blue: anti-correlated motion. Dashed lines separate A/B/C; gold guides mark every Cys residue in each chain.", ha="center", fontsize=5.9, color="#46515d")
    return fig


def write_source_data(results: list[dict]) -> None:
    rmsf_rows = []
    dccm_rows = []
    for result in results:
        label = result["spec"]["label"]
        chain_length = result["chain_length"]
        for residue_index, rmsf in enumerate(result["rmsf"], start=1):
            rmsf_rows.append(
                {
                    "system": label,
                    "chain": "ABC"[(residue_index - 1) // chain_length],
                    "segment_residue_position": (residue_index - 1) % chain_length + 1,
                    "is_Cys": (residue_index - 1) % chain_length + 1 in result["cys_positions"],
                    "RMSF_A": f"{rmsf:.8f}",
                }
            )
        for i in range(result["dccm"].shape[0]):
            for j in range(result["dccm"].shape[1]):
                dccm_rows.append(
                    {
                        "system": label,
                        "residue_i": i + 1,
                        "chain_i": "ABC"[i // chain_length],
                        "residue_j": j + 1,
                        "chain_j": "ABC"[j // chain_length],
                        "DCCM": f"{result['dccm'][i, j]:.8f}",
                    }
                )
    with (SOURCE_DIR / "rmsf_per_residue.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rmsf_rows[0]))
        writer.writeheader()
        writer.writerows(rmsf_rows)
    with (SOURCE_DIR / "dccm_matrices.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dccm_rows[0]))
        writer.writeheader()
        writer.writerows(dccm_rows)

    summary_rows = []
    for result in results:
        matrix = result["dccm"]
        chain_length = result["chain_length"]
        block_means = {}
        for first, second, name in ((0, 1, "AB"), (0, 2, "AC"), (1, 2, "BC")):
            first_slice = slice(first * chain_length, (first + 1) * chain_length)
            second_slice = slice(second * chain_length, (second + 1) * chain_length)
            block_means[f"mean_crosschain_DCCM_{name}"] = float(matrix[first_slice, second_slice].mean())
        cys_indices = [chain * chain_length + position - 1 for chain in range(3) for position in result["cys_positions"]]
        cys_rmsf = float(np.mean(result["rmsf"][cys_indices])) if cys_indices else float("nan")
        summary_rows.append(
            {
                "system": result["spec"]["label"],
                "segment_residues_per_chain": chain_length,
                "Cys_positions_per_chain": ",".join(map(str, result["cys_positions"])) or "none",
                "analysis_frames": result["frames"],
                "analysis_start_ns": result["analysis_start_ns"],
                "analysis_end_ns": result["analysis_end_ns"],
                "mean_RMSF_A": float(result["rmsf"].mean()),
                "max_RMSF_A": float(result["rmsf"].max()),
                "mean_Cys_RMSF_A": cys_rmsf,
                **block_means,
            }
        )
    with (OUTPUT_DIR / "dynamic_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    return summary_rows


def write_documentation(results: list[dict], rmsf_min: float, rmsf_max: float, summary_rows: list[dict]) -> None:
    manifest = {
        "analysis": "Panel A conformational ensembles/RMSF putty and Panel B DCCM",
        "systems": [item["spec"]["label"] for item in results],
        "common_settings": {
            "chain_count": 3,
            "alignment_anchor_residues_per_chain": CORE_RESIDUES_PER_CHAIN,
            "reported_segment_residues_per_chain": {item["spec"]["label"]: item["chain_length"] for item in results},
            "Cys_positions_per_chain": {item["spec"]["label"]: item["cys_positions"] for item in results},
            "burn_in_ns": BURN_IN_NS,
            "ensemble_frames_per_system": ENSEMBLE_FRAMES,
            "pbc_handling": "minimum-image clustering of chains 2 and 3 around chain 1 before alignment",
            "alignment": "unweighted Kabsch alignment of the centered common 30-residue CA anchor to the initial clustered frame of the same construct; all simulated residues are then transformed and analyzed",
            "rmsf_color_range_A": [float(rmsf_min), float(rmsf_max)],
            "dccm": "Cij=<dri dot drj>/sqrt(<|dri|^2><|drj|^2>) using aligned CA fluctuations of every simulated residue, including terminal Cys residues",
            "statistical_scope": "one trajectory per construct; trajectory-level descriptive dynamics only",
        },
        "summary_rows": summary_rows,
    }
    (OUTPUT_DIR / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    legend = f"""Figure X | Ensemble convergence, residue-level flexibility and dynamic coupling of four collagen-like constructs.

Panel A, {ENSEMBLE_FRAMES} evenly spaced frames from {BURN_IN_NS:.0f}-200 ns were PBC-clustered and aligned through a common centered 30-residue anchor before transparent overlay. Every residue present in each simulated segment, including terminal and internal Cys residues, is rendered; the representative backbone is RMSF-encoded by color and putty radius. The same RMSF color range ({rmsf_min:.2f}-{rmsf_max:.2f} A) is used for every construct. Panel B, DCCM matrices were calculated from aligned CA fluctuations of all simulated residues (CC, 34 aa/chain; mC ok, 30 aa/chain; NC, 32 aa/chain; P10_0, 30 aa/chain); dashed boundaries delimit chains A, B and C, and gold guide lines identify Cys positions. One trajectory was available per construct, so the displays are descriptive trajectory-level analyses rather than replicate-based inference.
"""
    (OUTPUT_DIR / "figure_legend.txt").write_text(legend, encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    configure_style()
    results = []
    for spec in SYSTEMS:
        print(f"Analyzing {spec['label']}...", flush=True)
        results.append(analyze_system(spec))
    all_rmsf = np.concatenate([result["rmsf"] for result in results])
    rmsf_min, rmsf_max = np.quantile(all_rmsf, [0.05, 0.95])
    rmsf_max = max(float(rmsf_max), float(rmsf_min) + 0.05)
    for result in results:
        _, _, pml = write_structures(result, float(rmsf_min), float(rmsf_max))
        print(f"Rendering {result['spec']['label']} with PyMOL...", flush=True)
        run_pymol(pml)
    save_figure(draw_panel_a(results, float(rmsf_min), float(rmsf_max)), "panel_A_conformational_ensembles_rmsf_putty")
    save_figure(draw_panel_b_dccm(results), "panel_B_dccm_matrices")
    for extension in ("pdf", "png", "svg", "tiff"):
        (OUTPUT_DIR / f"panel_B_dccm_change_vs_P10_0.{extension}").unlink(missing_ok=True)
    summary_rows = write_source_data(results)
    write_documentation(results, float(rmsf_min), float(rmsf_max), summary_rows)
    print(f"Completed: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
