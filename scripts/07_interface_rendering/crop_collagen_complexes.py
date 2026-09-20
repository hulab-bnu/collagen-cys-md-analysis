from __future__ import annotations

import json
import os
from pathlib import Path

from analyze_collagen_docking import analyze


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(os.environ.get("COLLAGEN_DOCKING_SOURCE", REPO_ROOT / "data" / "docking_structures"))
OUTPUT_ROOT = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "cropped_collagen"

PLANS = {
    r"CC\CG_P10_GC-20C_P25_53C.pdb.pdb": {
        "strict_equal": (20, 53),
        "interface_safe": (19, 53),
    },
    r"CC\para-P10_GC-P25_70C.pdb.pdb": {
        "strict_equal": (39, 70),
        "interface_safe": (38, 70),
    },
    r"mC\P10_16cys-P25_35cys.pdb.pdb": {
        "strict_equal": (24, 53),
        "interface_safe": (24, 53),
    },
    r"NC\para-P10_GC-P25_70C.pdb.pdb": {
        "strict_equal": (39, 70),
        "interface_safe": (38, 70),
    },
}


def atom_identity(line: str):
    chain = line[21:22].strip() or "_"
    resseq = int(line[22:26])
    return chain, resseq


def crop_pdb(source: Path, output: Path, start: int, end: int, mode: str):
    raw_lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    output.parent.mkdir(parents=True, exist_ok=True)

    kept_by_chain = {chain: [] for chain in "ABCDEF"}
    for line in raw_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain, resseq = atom_identity(line)
        if chain in "ABC":
            keep = start <= resseq <= end
        elif chain in "DEF":
            keep = True
        else:
            keep = False
        if keep:
            kept_by_chain[chain].append(line)

    out = [
        "REMARK 900 COLLAGEN COMPLEX CROPPED WITHOUT COORDINATE TRANSFORMATION",
        f"REMARK 901 SOURCE {source}",
        f"REMARK 902 MODE {mode}",
        f"REMARK 903 CHAINS A-C KEPT RESIDUES {start} TO {end}",
        "REMARK 904 CHAINS D-F KEPT IN FULL; ORIGINAL RESIDUE NUMBERS RETAINED",
        "REMARK 905 REGENERATE TERMINI AND HYDROGENS DURING TOPOLOGY BUILDING",
    ]
    serial = 1
    for chain in "ABCDEF":
        if not kept_by_chain[chain]:
            raise ValueError(f"{source}: no atoms retained for chain {chain}")
        for line in kept_by_chain[chain]:
            out.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        out.append(f"TER   {serial:5d}      {kept_by_chain[chain][-1][17:27]}")
        serial += 1
    out.append("END")
    output.write_text("\n".join(out) + "\n", encoding="ascii")


def main():
    manifest = []
    seen_outputs = set()
    for relative, variants in PLANS.items():
        source = SOURCE_ROOT / relative
        source_result = analyze(source, cutoff=4.5)
        for mode, (start, end) in variants.items():
            output = OUTPUT_ROOT / mode / Path(relative).with_suffix("").with_suffix(".pdb")
            output_key = output.resolve()
            if output_key in seen_outputs:
                raise ValueError(f"duplicate output path: {output}")
            seen_outputs.add(output_key)
            crop_pdb(source, output, start, end, mode)
            result = analyze(output, cutoff=4.5)

            source_contacts = {
                (item["long"], item["short"]): item["distance"] for item in source_result["contacts"]
            }
            output_contacts = {
                (item["long"], item["short"]): item["distance"] for item in result["contacts"]
            }
            lost = sorted(
                [
                    {"long": key[0], "short": key[1], "distance": distance}
                    for key, distance in source_contacts.items()
                    if key not in output_contacts
                ],
                key=lambda item: item["distance"],
            )

            long_counts = [result["chains"][chain]["residue_count"] for chain in "ABC"]
            short_counts = [result["chains"][chain]["residue_count"] for chain in "DEF"]
            manifest.append(
                {
                    "source": str(source),
                    "output": str(output.relative_to(OUTPUT_ROOT)),
                    "mode": mode,
                    "long_window": [start, end],
                    "long_chain_residue_counts": long_counts,
                    "short_chain_residue_counts": short_counts,
                    "equal_length": long_counts == short_counts,
                    "source_contact_pairs_4.5A": len(source_contacts),
                    "retained_contact_pairs_4.5A": len(output_contacts),
                    "retained_fraction": round(len(output_contacts) / len(source_contacts), 4),
                    "lost_contacts": lost,
                    "minimum_inter_molecule_distance_A": result["minimum_inter_molecule_distance"],
                }
            )

    manifest_path = OUTPUT_ROOT / "crop_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    readme = OUTPUT_ROOT / "README.md"
    lines = [
        "# Cropped collagen docking complexes",
        "",
        "Original files were read from the configured COLLAGEN_DOCKING_SOURCE directory and were not modified.",
        "",
        "- `strict_equal`: chains A-C have exactly the same residue count as chains D-F.",
        "- `interface_safe`: retains every original inter-molecule heavy-atom contact within 4.5 Å;",
        "  for the CC/para cases this makes chains A-C one residue longer than D-F.",
        "- Coordinates and original residue numbers are unchanged.",
        "- Chain separators (`TER`) and a final `END` record were standardized.",
        "- These are coordinate inputs, not complete MD systems. Rebuild termini/hydrogens and topology.",
        "",
        "See `crop_manifest.json` for per-file windows and contact-retention validation.",
        "",
    ]
    readme.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
