"""Export the production 0-ns mC protein PDB from the archived full-system XTC."""

import argparse
from pathlib import Path
import json
import hashlib
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--xtc", type=Path, required=True, help="Full-system XTC containing the 0-ns frame")
parser.add_argument("--template", type=Path, required=True, help="Protein PDB providing atom, chain and residue labels")
parser.add_argument("--index", type=Path, required=True, help="GROMACS index containing COMPLEX_Protein")
parser.add_argument("--dependency-dir", type=Path)
args = parser.parse_args()
if args.dependency_dir:
    sys.path = [str(args.dependency_dir)] + [p for p in sys.path if "site-packages" not in p.lower()]

import numpy as np
from MDAnalysis.lib.formats.libmdaxdr import XTCFile

xtc = args.xtc.resolve()
template = args.template.resolve()
ndx = args.index.resolve()
groups = {}
for line in ndx.read_text().splitlines():
    line = line.strip()
    if line.startswith("["):
        name = line.strip("[] ")
        groups[name] = []
    elif line and not line.startswith(";"):
        groups[name].extend(map(int, line.split()))
indices = np.array(groups["COMPLEX_Protein"], dtype=int)
records = template.read_text().splitlines()
atoms = [line for line in records if line.startswith(("ATOM  ", "HETATM"))]
assert len(atoms) == len(indices)
assert [int(line[6:11]) for line in atoms] == indices.tolist()
with XTCFile(str(xtc), "r") as reader:
    frame = next(iter(reader))
    assert abs(frame.time) < 1e-6
    coordinates = frame.x[indices - 1].astype(float) * 10.0
    box = frame.box.astype(float) * 10.0
assert np.isfinite(coordinates).all()
lengths = np.linalg.norm(box, axis=1)
angles = [np.degrees(np.arccos(np.clip(np.dot(box[i], box[j]) / (lengths[i] * lengths[j]), -1, 1)))
          for i, j in [(1, 2), (0, 2), (0, 1)]]
output = ["HEADER    M C PRODUCTION TRAJECTORY START", "TITLE     mC protein at 0 ns from repaired full-system XTC",
          "REMARK 950 COORDINATES EXTRACTED FROM THE ARCHIVED 0 NS FRAME.",
          "REMARK 950 ATOM/CHAIN/RESIDUE LABELS FROM THE MATCHING 20 NS PDB.",
          "REMARK 950 THIS IS A PRODUCTION START FRAME, NOT A PRE-MD DOCKING MODEL.",
          f"CRYST1{lengths[0]:9.3f}{lengths[1]:9.3f}{lengths[2]:9.3f}{angles[0]:7.2f}{angles[1]:7.2f}{angles[2]:7.2f} P 1           1"]
index = 0
for line in records:
    if line.startswith(("ATOM  ", "HETATM")):
        xyz = coordinates[index]
        output.append(line[:30] + "".join(f"{value:8.3f}" for value in xyz) + line[54:])
        index += 1
    elif line.startswith("TER"):
        output.append(line)
output.append("END")
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text("\n".join(output) + "\n", encoding="ascii")
roundtrip = np.array([[float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]]
                     for line in output if line.startswith(("ATOM  ", "HETATM"))])
assert np.max(np.abs(roundtrip - coordinates)) <= 0.000501
def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()
metadata = {"time_ns": 0, "coordinate_units": "Angstrom", "atom_count": len(atoms),
            "source_xtc": str(xtc), "source_xtc_sha256": sha(xtc),
            "label_template": str(template), "label_template_sha256": sha(template),
            "index": str(ndx), "index_sha256": sha(ndx),
            "coordinate_roundtrip_max_error_A": float(np.max(np.abs(roundtrip - coordinates)))}
args.output.with_suffix(".provenance.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(json.dumps(metadata, indent=2))
