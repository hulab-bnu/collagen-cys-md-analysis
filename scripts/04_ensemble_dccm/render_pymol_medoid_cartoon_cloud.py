#!/usr/bin/env python3
"""Render each MD ensemble with its real C-alpha medoid as the foreground.

The foreground is selected from the 16 sampled frames by minimum mean
pairwise C-alpha RMSD.  It is therefore an observed MD conformation, not an
average, fitted structure, spline, or idealized collagen model.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
ANALYSIS_ROOT = Path(os.environ.get("COLLAGEN_DCCM_SOURCE", REPO_ROOT / "data" / "ensemble_dccm"))
STRUCTURE_ROOT = ANALYSIS_ROOT / "structures"
OUTPUT_ROOT = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "ensemble_dccm"
RENDER_ROOT = OUTPUT_ROOT / "pymol_medoid_cartoon_cloud"
PYMOL = Path(os.environ.get("PYMOL_BIN", shutil.which("pymol") or "pymol"))
SYSTEMS = ("CC", "mC", "NC", "P10_0")

PANEL_RGB = (237, 248, 250)
CANVAS_SIZE = (1800, 540)


def stage_file(source: Path) -> Path:
    destination = RENDER_ROOT / source.name
    if not destination.exists() or destination.stat().st_size != source.stat().st_size:
        shutil.copy2(source, destination)
    return destination


def ca_models(pdb_path: Path) -> list[np.ndarray]:
    """Read every model's C-alpha coordinates from a multi-model PDB."""
    models: list[np.ndarray] = []
    current: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ATOM") and line[12:16].strip() == "CA":
            current.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("ENDMDL"):
            if current:
                models.append(np.asarray(current, dtype=float))
    if len(models) < 2:
        raise ValueError(f"Expected multiple models in {pdb_path}")
    if len({model.shape for model in models}) != 1:
        raise ValueError(f"Inconsistent C-alpha count across frames in {pdb_path}")
    return models


def choose_medoid_state(pdb_path: Path) -> tuple[int, float]:
    """Return the 1-based actual frame nearest to all other aligned frames."""
    coordinates = np.stack(ca_models(pdb_path), axis=0)
    delta = coordinates[:, None, :, :] - coordinates[None, :, :, :]
    rmsd = np.sqrt(np.mean(delta * delta, axis=(2, 3)))
    mean_rmsd = rmsd.mean(axis=1)
    index = int(np.argmin(mean_rmsd))
    return index + 1, float(mean_rmsd[index])


def write_pml(key: str, ensemble: Path, medoid_state: int) -> Path:
    pml_path = RENDER_ROOT / f"{key}_medoid_cartoon_cloud.pml"
    pml_path.write_text(
        "\n".join(
            [
                "reinitialize",
                "set_color panel_blue, [0.929, 0.973, 0.980]",
                "set_color cloud_bluegray, [0.46, 0.59, 0.63]",
                "set_color chain_blue, [0.04, 0.38, 0.60]",
                "set_color chain_teal, [0.05, 0.51, 0.56]",
                "set_color chain_mint, [0.31, 0.66, 0.53]",
                "set_color cys_gold, [0.95, 0.55, 0.03]",
                f"load {ensemble.name}, cloud",
                # This exact state is the C-alpha medoid of the 16-frame cloud.
                f"create front, cloud, {medoid_state}, 1",
                "set all_states, on",
                "hide everything, all",
                "show cartoon, cloud",
                "color cloud_bluegray, cloud",
                "set cartoon_transparency, 0.86, cloud",
                "set cartoon_smooth_loops, on",
                "set cartoon_sampling, 8",
                "set cartoon_loop_radius, 0.24, cloud",
                # Standard tubular PyMOL cartoon of an actual MD snapshot.
                "show cartoon, front",
                "set cartoon_transparency, 0.0, front",
                "set cartoon_smooth_loops, on",
                "set cartoon_sampling, 8",
                "set cartoon_loop_radius, 0.34, front",
                "color chain_blue, front and chain A",
                "color chain_teal, front and chain B",
                "color chain_mint, front and chain C",
                "show sticks, front and resn CYS+CYX+CYM",
                "show spheres, front and resn CYS+CYX+CYM",
                "color cys_gold, front and resn CYS+CYX+CYM",
                "set sphere_scale, 0.28, front and resn CYS+CYX+CYM",
                "set stick_radius, 0.12, front and resn CYS+CYX+CYM",
                "set orthoscopic, on",
                "set depth_cue, 0",
                "set ray_shadows, 0",
                "set antialias, 2",
                "set ray_opaque_background, on",
                "bg_color panel_blue",
                "orient front",
                "zoom front, 0.72",
                f"png {key}_medoid_cartoon_cloud.png, {CANVAS_SIZE[0]}, {CANVAS_SIZE[1]}, 300, 1",
                "quit",
            ]
        ),
        encoding="ascii",
    )
    return pml_path


def crop_to_plate(raw_png: Path, key: str) -> None:
    with Image.open(raw_png) as raw:
        image = raw.convert("RGB")
        background = Image.new("RGB", image.size, PANEL_RGB)
        mask = ImageChops.difference(image, background).convert("L").point(lambda value: 255 if value > 12 else 0)
        bbox = mask.getbbox() or (0, 0, image.width, image.height)
        left, top, right, bottom = bbox
        crop = image.crop((max(0, left - 80), max(0, top - 56), min(image.width, right + 80), min(image.height, bottom + 56)))
        plate = Image.new("RGB", CANVAS_SIZE, PANEL_RGB)
        scale = min((CANVAS_SIZE[0] - 100) / crop.width, (CANVAS_SIZE[1] - 74) / crop.height)
        resized = crop.resize((round(crop.width * scale), round(crop.height * scale)), Image.Resampling.LANCZOS)
        plate.paste(resized, ((CANVAS_SIZE[0] - resized.width) // 2, (CANVAS_SIZE[1] - resized.height) // 2))
        plate.save(OUTPUT_ROOT / f"Fig_docking_MD_{key}_pymol_medoid_cartoon_cloud.png", dpi=(300, 300))
        plate.save(
            OUTPUT_ROOT / f"Fig_docking_MD_{key}_pymol_medoid_cartoon_cloud.tiff",
            dpi=(300, 300),
            compression="tiff_lzw",
        )


def main() -> None:
    if not PYMOL.is_file():
        raise FileNotFoundError(PYMOL)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    RENDER_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for key in SYSTEMS:
        ensemble = stage_file(STRUCTURE_ROOT / f"{key}_ensemble_16_frames.pdb")
        medoid_state, mean_rmsd = choose_medoid_state(ensemble)
        pml = write_pml(key, ensemble, medoid_state)
        result = subprocess.run([str(PYMOL), "-cq", str(pml)], cwd=RENDER_ROOT, capture_output=True, text=True, timeout=180)
        if result.returncode:
            raise RuntimeError(f"PyMOL failed for {key}:\n{result.stdout}\n{result.stderr}")
        crop_to_plate(RENDER_ROOT / f"{key}_medoid_cartoon_cloud.png", key)
        rows.append({"construct": key, "medoid_state_1_based": medoid_state, "mean_CA_RMSD_to_16_frames_A": f"{mean_rmsd:.4f}"})
        print(f"Rendered {key}: medoid state {medoid_state}")

    with (OUTPUT_ROOT / "pymol_medoid_cartoon_cloud_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
