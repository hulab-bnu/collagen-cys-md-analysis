#!/usr/bin/env python3
"""Render a transparent PyMOL cartoon plate for four collagen ensembles."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SOURCE_ROOT = Path(os.environ.get("COLLAGEN_ENSEMBLE_STRUCTURES", REPO_ROOT / "data" / "ensemble_dccm" / "structures"))
OUTPUT_DIR = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "ensemble_dccm"
RENDER_DIR = OUTPUT_DIR / "pymol_transparent_ensemble_renders"
PYMOL = Path(os.environ.get("PYMOL_BIN", shutil.which("pymol") or "pymol"))
STEM = "Fig_conformational_ensembles_Cys_mapping_PyMOL_transparent_300dpi"

SYSTEMS = (
    ("CC", "CC_ensemble_16_frames.pdb", "CC_rmsf_putty.pdb"),
    ("mC", "mC_ensemble_16_frames.pdb", "mC_rmsf_putty.pdb"),
    ("NC", "NC_ensemble_16_frames.pdb", "NC_rmsf_putty.pdb"),
    ("P10", "P10_0_ensemble_16_frames.pdb", "P10_0_rmsf_putty.pdb"),
)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "none",
            "figure.facecolor": "none",
        }
    )


def make_pml(ensemble: Path, putty: Path, output_png: Path) -> str:
    """Use alpha-enabled ray tracing so the plate has a true transparent background."""
    return "\n".join(
        (
            "reinitialize",
            "viewport 2200, 560",
            f'load "{ensemble.as_posix()}", ensemble',
            f'load "{putty.as_posix()}", reference',
            "set all_states, on",
            "set ray_opaque_background, off",
            "set ray_shadows, off",
            "set antialias, 2",
            "set orthoscopic, on",
            "set cartoon_sampling, 14",
            "set_color ensemble_grey, [0.69, 0.78, 0.80]",
            "set_color chain_a, [0.18, 0.49, 0.66]",
            "set_color chain_b, [0.26, 0.66, 0.60]",
            "set_color chain_c, [0.52, 0.77, 0.58]",
            "set_color cysteine_gold, [0.93, 0.60, 0.12]",
            "hide everything, all",
            "show cartoon, ensemble",
            "color ensemble_grey, ensemble",
            "set cartoon_transparency, 0.78, ensemble",
            "show cartoon, reference",
            "cartoon putty, reference",
            "set cartoon_putty_transform, 0, reference",
            "set cartoon_putty_scale_min, 0.55, reference",
            "set cartoon_putty_scale_max, 1.25, reference",
            "color chain_a, reference and chain A",
            "color chain_b, reference and chain B",
            "color chain_c, reference and chain C",
            "show spheres, reference and resn CYS+CYX+CYM",
            "show sticks, reference and resn CYS+CYX+CYM",
            "color cysteine_gold, reference and resn CYS+CYX+CYM",
            "set sphere_scale, 0.25, reference and resn CYS+CYX+CYM",
            "set stick_radius, 0.10, reference and resn CYS+CYX+CYM",
            "orient reference",
            "zoom reference, 1.18",
            f"png {output_png.as_posix()}, 2200, 560, 300, 1",
            "quit",
        )
    )


def render_system(key: str, ensemble_name: str, putty_name: str) -> Path:
    ensemble = SOURCE_ROOT / ensemble_name
    putty = SOURCE_ROOT / putty_name
    if not ensemble.is_file() or not putty.is_file():
        raise FileNotFoundError(f"Missing PDB input for {key}")

    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    cache_dir = RENDER_DIR / "source_pdb"
    cache_dir.mkdir(exist_ok=True)
    cached_ensemble = cache_dir / f"{key}_ensemble.pdb"
    cached_putty = cache_dir / f"{key}_putty.pdb"
    for source, cached in ((ensemble, cached_ensemble), (putty, cached_putty)):
        if not cached.is_file() or cached.stat().st_size != source.stat().st_size:
            shutil.copy2(source, cached)

    pml_path = RENDER_DIR / f"{key}_transparent_ensemble.pml"
    png_path = RENDER_DIR / f"{key}_transparent_ensemble.png"
    pml_path.write_text(make_pml(cached_ensemble, cached_putty, png_path), encoding="ascii")
    completed = subprocess.run(
        [str(PYMOL), "-cq", str(pml_path)],
        cwd=RENDER_DIR,
        capture_output=True,
        text=True,
        timeout=180,
    )
    if completed.returncode or not png_path.is_file():
        raise RuntimeError(f"PyMOL failed for {key}:\n{completed.stdout}\n{completed.stderr}")
    return png_path


def trim_alpha(image: np.ndarray) -> np.ndarray:
    """Remove only transparent PyMOL canvas space, retaining a small margin."""
    if image.ndim != 3 or image.shape[-1] != 4:
        raise ValueError("PyMOL did not return an RGBA image")
    foreground = image[..., 3] > 0.01
    rows, columns = np.where(foreground)
    if not len(rows):
        raise ValueError("Rendered PyMOL image contains no visible atoms")
    margin_y, margin_x = 36, 64
    y0 = max(0, int(rows.min()) - margin_y)
    y1 = min(image.shape[0], int(rows.max()) + margin_y + 1)
    x0 = max(0, int(columns.min()) - margin_x)
    x1 = min(image.shape[1], int(columns.max()) + margin_x + 1)
    return image[y0:y1, x0:x1]


def save_panel(fig: plt.Figure) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(
            OUTPUT_DIR / f"{STEM}.{extension}",
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            pad_inches=0.01,
            transparent=True,
            facecolor="none",
        )
    fig.savefig(
        OUTPUT_DIR / f"{STEM}.tiff",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.01,
        transparent=True,
        facecolor="none",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def save_single_structure(image: np.ndarray, key: str) -> None:
    """Save a cropped PyMOL RGBA render without resampling its atom detail."""
    uint8_image = np.rint(np.clip(image, 0.0, 1.0) * 255).astype(np.uint8)
    output = OUTPUT_DIR / f"{STEM}_{key}.png"
    Image.fromarray(uint8_image, mode="RGBA").save(output, dpi=(300, 300))


def main() -> None:
    configure_style()
    if not PYMOL.is_file():
        raise FileNotFoundError(f"PyMOL executable not found: {PYMOL}")

    rendered: list[tuple[str, np.ndarray]] = []
    for key, ensemble_name, putty_name in SYSTEMS:
        print(f"Rendering {key}...", flush=True)
        image = trim_alpha(plt.imread(render_system(key, ensemble_name, putty_name)))
        save_single_structure(image, key)
        rendered.append((key, image))

    fig, axes = plt.subplots(4, 1, figsize=(4.75, 4.50))
    fig.patch.set_alpha(0)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005, hspace=0.025)
    for axis, (_, image) in zip(axes, rendered):
        axis.set_facecolor("none")
        axis.imshow(image, interpolation="lanczos", aspect="auto")
        axis.set_axis_off()
    save_panel(fig)
    plt.close(fig)
    print(f"Saved {STEM} to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
