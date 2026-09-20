#!/usr/bin/env python3
"""Calculate unbiased Cys SG-SG proximity statistics from GROMACS XVG output.

This program is intentionally dependency-free: Python 3.8+ is sufficient.
"""

from __future__ import annotations

import argparse
import csv
import math
from functools import lru_cache
from itertools import combinations
from pathlib import Path


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def read_ndx_indices(path: Path) -> list[int]:
    indices: list[int] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("["):
                continue
            indices.extend(int(value) for value in stripped.split())
    if not indices:
        raise ValueError(f"No atom indices found in {path}")
    return indices


def read_pdb_sg_atoms(path: Path) -> dict[int, dict[str, object]]:
    atoms: dict[int, dict[str, object]] = {}
    valid_resnames = {"CYS", "CYX", "CYM"}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "SG" or line[17:20].strip() not in valid_resnames:
                continue
            serial = int(line[6:11])
            chain = line[21].strip() or "?"
            residue = int(line[22:26])
            insertion = line[26].strip()
            trimer = "ABC" if chain in {"A", "B", "C"} else "DEF" if chain in {"D", "E", "F"} else "other"
            atoms[serial] = {
                "gmx_atom_index": serial,
                "chain": chain,
                "residue_number": residue,
                "insertion_code": insertion,
                "trimer": trimer,
                "label": f"{chain}:C{residue}{insertion}",
            }
    return atoms


def read_xvg_coordinates(path: Path, n_atoms: int) -> tuple[list[float], list[list[tuple[float, float, float]]]]:
    times: list[float] = []
    frames: list[list[tuple[float, float, float]]] = []
    expected = 1 + n_atoms * 3
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@")):
                continue
            values = [float(value) for value in stripped.split()]
            if len(values) != expected:
                raise ValueError(
                    f"Unexpected coordinate layout in {path}: found {len(values)} columns, expected {expected}. "
                    "The GROMACS trajectory command must output x/y/z for every selected SG atom."
                )
            times.append(values[0])
            frames.append([tuple(values[1 + 3 * i : 4 + 3 * i]) for i in range(n_atoms)])
    if not frames:
        raise ValueError(f"No coordinate records found in {path}")
    return times, frames


def distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def pair_class(first: dict[str, object], second: dict[str, object]) -> str:
    if first["chain"] == second["chain"]:
        return "same_chain"
    if first["trimer"] == second["trimer"]:
        return "within_trimer"
    if {first["trimer"], second["trimer"]} == {"ABC", "DEF"}:
        return "between_trimers"
    return "other"


def maximum_matching_count(edges: list[tuple[int, int]], n_atoms: int) -> int:
    """Exact maximum number of non-overlapping contact pairs for a small graph."""
    neighbor_masks = [0] * n_atoms
    for first, second in edges:
        neighbor_masks[first] |= 1 << second
        neighbor_masks[second] |= 1 << first

    @lru_cache(maxsize=None)
    def solve(mask: int) -> int:
        if mask == 0:
            return 0
        first_bit = mask & -mask
        first = first_bit.bit_length() - 1
        best = solve(mask ^ first_bit)
        partners = neighbor_masks[first] & mask
        while partners:
            partner_bit = partners & -partners
            best = max(best, 1 + solve(mask ^ first_bit ^ partner_bit))
            partners ^= partner_bit
        return best

    return solve((1 << n_atoms) - 1)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--reference-pdb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cutoff-nm", type=float, default=0.40)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    indices = read_ndx_indices(args.index)
    pdb_atoms = read_pdb_sg_atoms(args.reference_pdb)
    missing = [index for index in indices if index not in pdb_atoms]
    if missing:
        raise ValueError(
            "Selected SG atom indices do not match the reference PDB atom serials. "
            f"First unmatched indices: {missing[:5]}. Use a snapshot PDB generated from the same repaired trajectory/topology."
        )
    atoms = [pdb_atoms[index] for index in indices]
    times, frames = read_xvg_coordinates(args.coordinates, len(atoms))
    write_csv(
        output_dir / "sg_atom_mapping.csv",
        ["sg_order", "gmx_atom_index", "chain", "residue_number", "insertion_code", "trimer", "label"],
        [{"sg_order": order, **atom} for order, atom in enumerate(atoms)],
    )

    pair_rows: list[dict[str, object]] = []
    time_rows: list[dict[str, object]] = []
    pair_distances: dict[tuple[int, int], list[float]] = {}
    pair_info: dict[tuple[int, int], dict[str, object]] = {}
    for first, second in combinations(range(len(atoms)), 2):
        label = f"{atoms[first]['label']} -- {atoms[second]['label']}"
        info = {
            "pair": label,
            "sg_order_1": first,
            "sg_order_2": second,
            "chain_1": atoms[first]["chain"],
            "residue_1": atoms[first]["residue_number"],
            "chain_2": atoms[second]["chain"],
            "residue_2": atoms[second]["residue_number"],
            "pair_class": pair_class(atoms[first], atoms[second]),
        }
        distances = [distance(frame[first], frame[second]) for frame in frames]
        pair_info[(first, second)] = info
        pair_distances[(first, second)] = distances
        pair_rows.append(
            {
                **info,
                "n_frames": len(distances),
                "median_sg_sg_nm": quantile(distances, 0.50),
                "q1_sg_sg_nm": quantile(distances, 0.25),
                "q3_sg_sg_nm": quantile(distances, 0.75),
                "minimum_sg_sg_nm": min(distances),
                "maximum_sg_sg_nm": max(distances),
                "fraction_lt_0p25_nm": sum(value < 0.25 for value in distances) / len(distances),
                "fraction_lt_0p35_nm": sum(value < 0.35 for value in distances) / len(distances),
                "fraction_lt_0p40_nm": sum(value < 0.40 for value in distances) / len(distances),
                "fraction_lt_0p50_nm": sum(value < 0.50 for value in distances) / len(distances),
            }
        )
        for time_ns, value in zip(times, distances):
            time_rows.append({"time_ns": time_ns, **info, "sg_sg_distance_nm": value})

    pair_rows.sort(key=lambda row: (row["median_sg_sg_nm"], row["minimum_sg_sg_nm"]))
    fields = list(pair_rows[0].keys())
    write_csv(output_dir / "all_sg_sg_pair_summary.csv", fields, pair_rows)
    write_csv(output_dir / "cross_chain_sg_sg_pair_summary.csv", fields, [row for row in pair_rows if row["pair_class"] != "same_chain"])
    write_csv(output_dir / "cross_trimer_sg_sg_pair_summary.csv", fields, [row for row in pair_rows if row["pair_class"] == "between_trimers"])
    write_csv(output_dir / "sg_sg_pair_distance_timeseries.csv", list(time_rows[0].keys()), time_rows)

    count_rows: list[dict[str, object]] = []
    for frame_index, time_ns in enumerate(times):
        cross_chain_edges = []
        cross_trimer_edges = []
        for pair, distances in pair_distances.items():
            if distances[frame_index] >= args.cutoff_nm:
                continue
            kind = pair_info[pair]["pair_class"]
            if kind != "same_chain":
                cross_chain_edges.append(pair)
            if kind == "between_trimers":
                cross_trimer_edges.append(pair)
        count_rows.append(
            {
                "time_ns": time_ns,
                "raw_cross_chain_pairs_lt_cutoff": len(cross_chain_edges),
                "max_nonoverlapping_cross_chain_pairs_lt_cutoff": maximum_matching_count(cross_chain_edges, len(atoms)),
                "raw_cross_trimer_pairs_lt_cutoff": len(cross_trimer_edges),
                "max_nonoverlapping_cross_trimer_pairs_lt_cutoff": maximum_matching_count(cross_trimer_edges, len(atoms)),
                "cutoff_nm": args.cutoff_nm,
            }
        )
    write_csv(output_dir / "potential_disulfide_contact_count.csv", list(count_rows[0].keys()), count_rows)

    candidate_rows = [row for row in pair_rows if row["pair_class"] != "same_chain" and row["fraction_lt_0p40_nm"] > 0]
    write_csv(output_dir / "candidate_sg_sg_pairs_lt_0p40.csv", fields, candidate_rows)
    print(f"Analysed {len(atoms)} SG atoms, {len(pair_rows)} atom pairs and {len(times)} trajectory frames.")
    print("Closest cross-trimer candidates:")
    cross = [row for row in pair_rows if row["pair_class"] == "between_trimers"][:10]
    for row in cross:
        print(
            f"  {row['pair']}: median {row['median_sg_sg_nm']:.3f} nm; "
            f"<0.40 nm occupancy {row['fraction_lt_0p40_nm']:.1%}"
        )


if __name__ == "__main__":
    main()
