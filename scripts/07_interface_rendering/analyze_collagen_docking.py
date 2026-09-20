from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path


def parse_pdb(path: Path):
    atoms = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atom = {
                "record": line[0:6].strip(),
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip() or "_",
                "resseq": int(line[22:26]),
                "icode": line[26:27].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip() or line[12:16].strip()[0],
                "line_no": line_no,
            }
        except (ValueError, IndexError) as exc:
            raise ValueError(f"{path}:{line_no}: malformed coordinate line") from exc
        atoms.append(atom)
    return atoms


def residue_key(atom):
    return atom["chain"], atom["resseq"], atom["icode"], atom["resname"]


def min_distance(a_atoms, b_atoms):
    best = math.inf
    for a in a_atoms:
        if a["element"].upper() == "H":
            continue
        for b in b_atoms:
            if b["element"].upper() == "H":
                continue
            dx = a["x"] - b["x"]
            dy = a["y"] - b["y"]
            dz = a["z"] - b["z"]
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < best:
                best = d2
    return math.sqrt(best)


def analyze(path: Path, cutoff: float):
    atoms = parse_pdb(path)
    residues = defaultdict(list)
    for atom in atoms:
        residues[residue_key(atom)].append(atom)

    chain_residues = defaultdict(list)
    for key in residues:
        chain_residues[key[0]].append(key)
    for chain in chain_residues:
        chain_residues[chain].sort(key=lambda k: (k[1], k[2]))

    chains = sorted(chain_residues)
    if len(chains) != 6:
        raise ValueError(f"{path}: expected six collagen chains, found {chains}")
    long_chains, short_chains = chains[:3], chains[3:]

    chain_summary = {}
    for chain in chains:
        keys = chain_residues[chain]
        chain_atoms = [a for key in keys for a in residues[key]]
        chain_summary[chain] = {
            "residue_count": len(keys),
            "first_residue": keys[0][1],
            "last_residue": keys[-1][1],
            "z_min": round(min(a["z"] for a in chain_atoms), 3),
            "z_max": round(max(a["z"] for a in chain_atoms), 3),
        }

    contacts = []
    contacting_long = set()
    contacting_short = set()
    min_pair = (math.inf, None, None)
    for lchain in long_chains:
        for lkey in chain_residues[lchain]:
            for schain in short_chains:
                for skey in chain_residues[schain]:
                    d = min_distance(residues[lkey], residues[skey])
                    if d < min_pair[0]:
                        min_pair = (d, lkey, skey)
                    if d <= cutoff:
                        contacts.append(
                            {
                                "long": f"{lkey[0]}:{lkey[3]}{lkey[1]}",
                                "short": f"{skey[0]}:{skey[3]}{skey[1]}",
                                "distance": round(d, 3),
                            }
                        )
                        contacting_long.add(lkey)
                        contacting_short.add(skey)

    long_by_chain = {}
    for chain in long_chains:
        nums = sorted(k[1] for k in contacting_long if k[0] == chain)
        long_by_chain[chain] = nums

    short_len = min(chain_summary[c]["residue_count"] for c in short_chains)
    candidate_windows = []
    long_first = max(chain_summary[c]["first_residue"] for c in long_chains)
    long_last = min(chain_summary[c]["last_residue"] for c in long_chains)
    contact_nums = sorted(k[1] for k in contacting_long)
    if contact_nums:
        cmin, cmax = min(contact_nums), max(contact_nums)
        for start in range(long_first, long_last - short_len + 2):
            end = start + short_len - 1
            if start <= cmin and end >= cmax:
                margin = min(cmin - start, end - cmax)
                center_gap = abs((start + end) / 2 - (cmin + cmax) / 2)
                candidate_windows.append((margin, -center_gap, start, end))
        candidate_windows.sort(reverse=True)

    result = {
        "file": str(path),
        "atom_count": len(atoms),
        "chains": chain_summary,
        "assumed_long_molecule": long_chains,
        "assumed_short_molecule": short_chains,
        "contact_cutoff_angstrom": cutoff,
        "contact_pair_count": len(contacts),
        "minimum_inter_molecule_distance": round(min_pair[0], 3),
        "closest_pair": {
            "long": f"{min_pair[1][0]}:{min_pair[1][3]}{min_pair[1][1]}",
            "short": f"{min_pair[2][0]}:{min_pair[2][3]}{min_pair[2][1]}",
        },
        "contacting_long_residues": long_by_chain,
        "contacting_short_residue_count": len(contacting_short),
        "recommended_equal_length_window": (
            {"start": candidate_windows[0][2], "end": candidate_windows[0][3], "length": short_len}
            if candidate_windows
            else None
        ),
        "contacts": sorted(contacts, key=lambda x: x["distance"]),
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = [analyze(path, args.cutoff) for path in args.paths]
    text = json.dumps(results, ensure_ascii=False, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
