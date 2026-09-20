#!/usr/bin/env python3
"""Create an interface-preserving P10:P10 control from the P10:P25 docking PDB."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


P10_LENGTH = 30
P25_INTERFACE_START = 22
P25_INTERFACE_END = P25_INTERFACE_START + P10_LENGTH - 1
FIRST_TRIMER = {"A", "B", "C"}
SECOND_TRIMER = {"D", "E", "F"}


def atom_line(serial: int, source_line: str, new_residue_number: int) -> str:
    """Renumber an ATOM/HETATM line while retaining atom names, coordinates and B factors."""
    return f"{source_line[:6]}{serial:5d}{source_line[11:22]}{new_residue_number:4d}{source_line[26:]}"


def ter_line(serial: int, previous_atom: str) -> str:
    """Write a standard TER record using the final residue identity of a chain."""
    return f"TER   {serial:5d}      {previous_atom[17:20]} {previous_atom[21]}{previous_atom[22:26]}\n"


def unique_residues(lines: list[str]) -> dict[str, list[tuple[int, str]]]:
    residues: dict[str, list[tuple[int, str]]] = defaultdict(list)
    seen: set[tuple[str, int, str]] = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21]
        residue = int(line[22:26])
        name = line[17:20].strip()
        key = (chain, residue, name)
        if key not in seen:
            residues[chain].append((residue, name))
            seen.add(key)
    return residues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    parser.add_argument("--manifest", type=Path, default=None)
    args = parser.parse_args()

    source_lines = args.input_pdb.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    output_lines = [
        "HEADER    P10-P10 INTER-TRIMER DOCKING CONTROL\n",
        "REMARK    Derived from P10-P25 docking coordinates without rigid-body reorientation.\n",
        "REMARK    Chains A-C: retained interface-centred P25 residues 23-52 and renumbered 1-30.\n",
        "REMARK    Chains D-F: original P10 residues 1-30 retained unchanged.\n",
        "REMARK    All six chains are P10 = (PPG)10; coordinates preserve the original interface.\n",
    ]

    serial = 1
    previous_chain: str | None = None
    previous_atom: str | None = None

    for line in source_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21]
        residue = int(line[22:26])
        if chain not in FIRST_TRIMER | SECOND_TRIMER:
            continue
        if chain in FIRST_TRIMER:
            if not (P25_INTERFACE_START <= residue <= P25_INTERFACE_END):
                continue
            new_residue = residue - P25_INTERFACE_START + 1
        else:
            if not (1 <= residue <= P10_LENGTH):
                continue
            new_residue = residue

        if previous_chain is not None and chain != previous_chain and previous_atom is not None:
            output_lines.append(ter_line(serial, previous_atom))
            serial += 1

        output_lines.append(atom_line(serial, line, new_residue))
        serial += 1
        previous_chain = chain
        previous_atom = line

    if previous_atom is None:
        raise RuntimeError("No retained P10 atoms were written; check the input chain and residue numbering.")
    output_lines.append(ter_line(serial, previous_atom))
    output_lines.append("END\n")

    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    args.output_pdb.write_text("".join(output_lines), encoding="ascii")

    residues = unique_residues(output_lines)
    expected_sequence = "PPG" * 10
    for chain in "ABCDEF":
        observed = residues.get(chain, [])
        if len(observed) != P10_LENGTH:
            raise RuntimeError(f"Chain {chain} has {len(observed)} residues, expected {P10_LENGTH}.")
        sequence = "".join({"PRO": "P", "GLY": "G"}.get(name, "X") for _, name in observed)
        if sequence != expected_sequence:
            raise RuntimeError(f"Chain {chain} is not (PPG)10: {sequence}")
        if [residue for residue, _ in observed] != list(range(1, P10_LENGTH + 1)):
            raise RuntimeError(f"Chain {chain} residue numbering is not 1-{P10_LENGTH}.")

    manifest = args.manifest or args.output_pdb.with_suffix(".txt")
    manifest.write_text(
        "P10:P10 docking control\n"
        "=====================\n"
        f"Input: {args.input_pdb.name}\n"
        f"Output: {args.output_pdb.name}\n"
        "Reference complex: P10:P25 docking pose.\n"
        "A-C trimer: P25 residues 23-52 retained, renumbered to 1-30.\n"
        "D-F trimer: original P10 residues 1-30 retained.\n"
        "All chains: P10 = (PPG)10.\n"
        "Rigid-body relationship: unchanged from the source docking pose.\n"
        "Recommended next step: energy minimization followed by the same GROMACS equilibration/production protocol used for CC, mC and NC.\n",
        encoding="ascii",
    )
    print(f"Created: {args.output_pdb}")
    print(f"Manifest: {manifest}")
    print("Validated six 30-residue (PPG)10 chains with preserved source coordinates.")


if __name__ == "__main__":
    main()
