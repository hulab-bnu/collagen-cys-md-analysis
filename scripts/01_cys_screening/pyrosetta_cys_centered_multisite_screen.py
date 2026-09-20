#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Centered GPP21 Cys-permissibility and construct-aware multi-site screening.

Scientific scope
----------------
This script evaluates local triple-helix compatibility in a 21-triplet window.
For every screened position, the target Gly-X-Y triplet is mapped to template
triplet 11. Real parent-sequence triplets fill the available 61-triplet domain;
GPP padding is used only beyond the real N/C boundaries.

Two screening levels are supported:

1. Single-site screening
   Every structurally permissible X/Y position is mutated to Cys at the center
   of a matched GPP21 window.

2. Construct-aware local multi-site screening
   For NC/NC1/NC3/NC5 (and optionally mC/CC), every experimental Cys is used as
   a window center. Other Cys sites from the same construct are included when
   they fall within the configurable local mutation radius. Multi-mutant energy
   non-additivity is calculated against matched single mutants in the same
   window and replicate.

Important limitation
--------------------
This is a local reduced-Cys permissibility model. It does not place every Cys
from a long construct into one 183-aa pose, form disulfide bonds, or model
D-periodic assembly. Terminal sites evaluated with GPP padding represent an
intrinsic, internalized local context and must be interpreted separately from
their native terminal topology.

Notebook use
------------
Edit NOTEBOOK_CONFIG below, then run the notebook cell containing this script.
No command-line arguments are parsed. When imported as a module, call
main(NOTEBOOK_CONFIG) explicitly.
"""

from __future__ import annotations

import csv
import json
import math
import os
import statistics
import sys
import zlib
from collections import defaultdict
from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Biological input and construct definitions
# ---------------------------------------------------------------------------

FULL_PARENT_COLLAGEN_DOMAIN = (
    "GPPGPPGPPGPPGPPGPPGPPGPPGPPGPP"
    "GPRGEQGPQGLPGKDGEAGAQGPAGPRGPQGPQGLPGPQGPAGPMGPAGFPGER"
    "GEKGEPGTQGAKGDRGETGPVGPRGERGEAGPAGKDGERGPVGPAGPRGPQGPQGLPGPQGPAGAQ"
    "GPPGPPGPPGPPGPPGPPGPPGPPGPPGPP"
    "GGP"
)

DOMAIN_LENGTH = 183
DOMAIN_TRIPLETS = 61
WINDOW_TRIPLETS = 21
WINDOW_CHAIN_LENGTH = 63
CENTER_TRIPLET = 11
N_CHAINS = 3
GPP_PADDING_TRIPLET = "GPP"

CONSTRUCT_CYS_SITES: Dict[str, List[int]] = {
    "mC": [167],
    "CC": [150, 183],
    "NC": [32, 183],
    "NC1": [32, 90, 183],
    "NC3": [32, 60, 90, 117, 183],
    "NC5": [32, 54, 72, 90, 105, 129, 183],
}

PROTECTED_MOTIFS = ("GFPGER",)


# ---------------------------------------------------------------------------
# Protocol defaults
# ---------------------------------------------------------------------------

OUTER_ANCHOR_TRIPLETS = frozenset((1, 2, 3, 19, 20, 21))
TRANSITION_TRIPLETS = frozenset((4, 5, 6, 7, 15, 16, 17, 18))
ANALYSIS_CORE_TRIPLETS = frozenset(range(8, 15))
MOVABLE_TRIPLETS = frozenset(range(4, 19))

DEFAULT_CONSTRUCT_MUTATION_RADIUS_TRIPLETS = 7
DEFAULT_GLY_ANALYSIS_RADIUS_TRIPLETS = 2
DEFAULT_MUTATION_CONTEXT_RADIUS_TRIPLETS = 2

DEFAULT_ANCHOR_COORD_SD = 0.50
DEFAULT_TRANSITION_COORD_SD = 1.50
DEFAULT_GLY_ANCHOR_DISTANCE_SD = 0.75
DEFAULT_ATOM_PAIR_CST_WEIGHT = 1.0
DEFAULT_COORDINATE_CST_WEIGHT = 0.5
DEFAULT_HBOND_ENERGY_CUTOFF = -0.10

DEFAULT_RMSD_FAIL_A = 1.50
DEFAULT_GLY_AREA_FAIL_A2 = 2.50
DEFAULT_FA_REP_FAIL_REU = 10.0

TRACKED_SCORE_TERM_NAMES = (
    "fa_atr",
    "fa_rep",
    "fa_sol",
    "fa_elec",
    "lk_ball_wtd",
    "fa_intra_rep",
    "fa_dun",
    "omega",
    "rama_prepro",
    "p_aa_pp",
    "pro_close",
    "ref",
    "hbond_sr_bb",
    "hbond_lr_bb",
    "hbond_bb_sc",
    "hbond_sc",
    "dslf_fa13",
)


@dataclass
class ScreenConfig:
    """Notebook-editable protocol settings; no command-line parsing is used."""

    template: str = "GPP21.pdb"
    output_dir: str = "centered_gpp21_cys_screen_outputs"
    mode: str = "both"
    constructs: Tuple[str, ...] = ("NC1", "NC3", "NC5")
    scorefxn: str = "ref2015"
    seed: int = 1111
    single_samples: int = 3
    construct_samples: int = 5
    background_relax_repeats: int = 1
    fast_relax_repeats: int = 2
    construct_mutation_radius_triplets: int = (
        DEFAULT_CONSTRUCT_MUTATION_RADIUS_TRIPLETS
    )
    mutation_context_radius_triplets: int = (
        DEFAULT_MUTATION_CONTEXT_RADIUS_TRIPLETS
    )
    gly_analysis_radius_triplets: int = DEFAULT_GLY_ANALYSIS_RADIUS_TRIPLETS
    anchor_coordinate_sd: float = DEFAULT_ANCHOR_COORD_SD
    transition_coordinate_sd: float = DEFAULT_TRANSITION_COORD_SD
    gly_anchor_distance_sd: float = DEFAULT_GLY_ANCHOR_DISTANCE_SD
    atom_pair_cst_weight: float = DEFAULT_ATOM_PAIR_CST_WEIGHT
    coordinate_cst_weight: float = DEFAULT_COORDINATE_CST_WEIGHT
    hbond_energy_cutoff: float = DEFAULT_HBOND_ENERGY_CUTOFF
    rmsd_fail_A: float = DEFAULT_RMSD_FAIL_A
    gly_area_fail_A2: float = DEFAULT_GLY_AREA_FAIL_A2
    fa_rep_fail_REU: float = DEFAULT_FA_REP_FAIL_REU
    include_xy_gly: bool = True
    include_protected_motifs: bool = False
    include_negative_controls: bool = False
    calculate_nonadditivity: bool = True
    dump_pdbs: bool = False
    dry_run: bool = False
    self_test: bool = False


# Edit these values directly in the notebook before running this cell.
NOTEBOOK_CONFIG = ScreenConfig(
    template="GPP21.pdb",
    output_dir="centered_gpp21_cys_screen_outputs",
    mode="both",  # "single", "construct", or "both"
    constructs=("NC1", "NC3", "NC5"),
    single_samples=3,
    construct_samples=5,
)


# PyRosetta symbols are imported lazily so notebook self-tests work without it.
pyrosetta = None
Pose = None
AtomID = None
AtomPairConstraint = None
CoordinateConstraint = None
HarmonicFunc = None
FastRelax = None
MoveMap = None
HBondSet = None
fill_hbond_set = None
SasaCalc = None
mutate_residue = None
pose_from_pdb = None
create_score_function = None
rosetta_random = None
scoring_module = None
SCORE_TYPES: Dict[str, Any] = {}


@dataclass(frozen=True)
class WindowSpec:
    center_global_triplet: int
    parent_sequence: str
    full_to_local: Dict[int, int]
    local_to_full: Dict[int, Optional[int]]
    left_padding_triplets: int
    right_padding_triplets: int


@dataclass(frozen=True)
class ScreenJob:
    kind: str
    label: str
    center_full_position: int
    center_global_triplet: int
    center_register: str
    mutation_full_positions: Tuple[int, ...]
    construct: str
    samples: int
    experimental_site: bool
    control_type: str = ""
    window_center_label: str = ""


# ---------------------------------------------------------------------------
# Pure sequence/window logic
# ---------------------------------------------------------------------------

def triplet_info(full_position: int) -> Tuple[int, str]:
    if not 1 <= full_position <= DOMAIN_LENGTH:
        raise ValueError(f"Full position out of range: {full_position}")
    triplet = (full_position - 1) // 3 + 1
    register = ("G", "X", "Y")[(full_position - 1) % 3]
    return triplet, register


def mutation_label(full_position: int, new_aa: str = "C") -> str:
    return f"{FULL_PARENT_COLLAGEN_DOMAIN[full_position - 1]}{full_position}{new_aa}"


def motif_intervals(sequence: str) -> List[Tuple[str, int, int]]:
    intervals: List[Tuple[str, int, int]] = []
    for motif in PROTECTED_MOTIFS:
        start = sequence.find(motif)
        while start >= 0:
            intervals.append((motif, start + 1, start + len(motif)))
            start = sequence.find(motif, start + 1)
    return intervals


PROTECTED_INTERVALS = motif_intervals(FULL_PARENT_COLLAGEN_DOMAIN)


def protected_motif_at(full_position: int) -> str:
    for motif, start, end in PROTECTED_INTERVALS:
        if start <= full_position <= end:
            return motif
    return ""


def candidate_allowed(
    full_position: int,
    include_xy_gly: bool = True,
    exclude_protected_motifs: bool = True,
) -> Tuple[bool, str]:
    aa = FULL_PARENT_COLLAGEN_DOMAIN[full_position - 1]
    _, register = triplet_info(full_position)
    reasons: List[str] = []

    if register == "G":
        reasons.append("triplet_leading_gly")
    if aa == "C":
        reasons.append("already_cys")
    if aa == "G" and register in ("X", "Y") and not include_xy_gly:
        reasons.append("xy_gly_disabled")
    motif = protected_motif_at(full_position)
    if motif and exclude_protected_motifs:
        reasons.append(f"protected_motif_{motif}")

    return not reasons, ";".join(reasons) if reasons else "allowed"


def build_centered_window(center_global_triplet: int) -> WindowSpec:
    if not 1 <= center_global_triplet <= DOMAIN_TRIPLETS:
        raise ValueError(f"Center triplet out of range: {center_global_triplet}")

    sequence_parts: List[str] = []
    full_to_local: Dict[int, int] = {}
    local_to_full: Dict[int, Optional[int]] = {}
    left_padding = 0
    right_padding = 0

    for local_triplet in range(1, WINDOW_TRIPLETS + 1):
        global_triplet = center_global_triplet + local_triplet - CENTER_TRIPLET
        local_start = (local_triplet - 1) * 3 + 1

        if 1 <= global_triplet <= DOMAIN_TRIPLETS:
            full_start = (global_triplet - 1) * 3 + 1
            triplet_seq = FULL_PARENT_COLLAGEN_DOMAIN[full_start - 1:full_start + 2]
            sequence_parts.append(triplet_seq)
            for offset in range(3):
                full_pos = full_start + offset
                local_pos = local_start + offset
                full_to_local[full_pos] = local_pos
                local_to_full[local_pos] = full_pos
        else:
            sequence_parts.append(GPP_PADDING_TRIPLET)
            if global_triplet < 1:
                left_padding += 1
            else:
                right_padding += 1
            for offset in range(3):
                local_to_full[local_start + offset] = None

    sequence = "".join(sequence_parts)
    if len(sequence) != WINDOW_CHAIN_LENGTH:
        raise AssertionError(f"Window length is {len(sequence)}, expected 63")

    return WindowSpec(
        center_global_triplet=center_global_triplet,
        parent_sequence=sequence,
        full_to_local=full_to_local,
        local_to_full=local_to_full,
        left_padding_triplets=left_padding,
        right_padding_triplets=right_padding,
    )


def local_triplet_for_position(local_position: int) -> int:
    return (local_position - 1) // 3 + 1


def construct_sites_in_local_radius(
    construct: str,
    center_full_position: int,
    radius_triplets: int,
) -> Tuple[int, ...]:
    center_triplet, _ = triplet_info(center_full_position)
    sites = []
    for full_pos in CONSTRUCT_CYS_SITES[construct]:
        site_triplet, _ = triplet_info(full_pos)
        if abs(site_triplet - center_triplet) <= radius_triplets:
            sites.append(full_pos)
    if center_full_position not in sites:
        raise AssertionError("Construct center was excluded from its own local mutation set")
    return tuple(sorted(sites))


def maximal_construct_clusters(
    construct: str,
    radius_triplets: int,
) -> List[Tuple[int, Tuple[int, ...]]]:
    """Return maximal site clusters that fit within center +/- radius triplets."""
    sites = sorted(CONSTRUCT_CYS_SITES[construct])
    feasible: List[Tuple[int, ...]] = []
    for left in range(len(sites)):
        for right in range(left, len(sites)):
            subset = tuple(sites[left:right + 1])
            triplets = [triplet_info(position)[0] for position in subset]
            if max(triplets) - min(triplets) <= 2 * radius_triplets:
                feasible.append(subset)

    maximal = [
        subset
        for subset in feasible
        if not any(set(subset) < set(other) for other in feasible)
    ]
    clusters: List[Tuple[int, Tuple[int, ...]]] = []
    for subset in maximal:
        triplets = [triplet_info(position)[0] for position in subset]
        lower_center = max(triplets) - radius_triplets
        upper_center = min(triplets) + radius_triplets
        midpoint = int(round((min(triplets) + max(triplets)) / 2.0))
        center = max(lower_center, min(upper_center, midpoint))
        center = max(1, min(DOMAIN_TRIPLETS, center))
        clusters.append((center, subset))
    return sorted(clusters, key=lambda item: (item[0], item[1]))


def build_candidate_library(args: ScreenConfig) -> List[Dict[str, Any]]:
    experimental_union = {p for sites in CONSTRUCT_CYS_SITES.values() for p in sites}
    rows: List[Dict[str, Any]] = []
    for full_pos in range(1, DOMAIN_LENGTH + 1):
        triplet, register = triplet_info(full_pos)
        allowed, reason = candidate_allowed(
            full_pos,
            include_xy_gly=args.include_xy_gly,
            exclude_protected_motifs=not args.include_protected_motifs,
        )
        spec = build_centered_window(triplet)
        rows.append({
            "Candidate": mutation_label(full_pos),
            "Full_Position": full_pos,
            "Original_AA": FULL_PARENT_COLLAGEN_DOMAIN[full_pos - 1],
            "Triplet_Index": triplet,
            "Triplet_Register": register,
            "Allowed": allowed,
            "Candidate_Reason": reason,
            "Protected_Motif": protected_motif_at(full_pos),
            "Experimental_Site": "yes" if full_pos in experimental_union else "no",
            "Centered_Local_Position": spec.full_to_local.get(full_pos, ""),
            "Left_GPP_Padding_Triplets": spec.left_padding_triplets,
            "Right_GPP_Padding_Triplets": spec.right_padding_triplets,
            "Context_Class": "padded_intrinsic" if (
                spec.left_padding_triplets or spec.right_padding_triplets
            ) else "fully_native_window",
        })
    return rows


def build_jobs(args: ScreenConfig, candidate_rows: Sequence[Dict[str, Any]]) -> List[ScreenJob]:
    jobs: List[ScreenJob] = []

    if args.mode in ("single", "both"):
        for row in candidate_rows:
            if not row["Allowed"]:
                continue
            full_pos = int(row["Full_Position"])
            triplet, register = triplet_info(full_pos)
            jobs.append(ScreenJob(
                kind="single",
                label=f"single_{mutation_label(full_pos)}",
                center_full_position=full_pos,
                center_global_triplet=triplet,
                center_register=register,
                mutation_full_positions=(full_pos,),
                construct="",
                samples=args.single_samples,
                experimental_site=row["Experimental_Site"] == "yes",
                window_center_label=f"site_{mutation_label(full_pos)}",
            ))

    if args.include_negative_controls and args.mode in ("single", "both"):
        for triplet in (6, 16, 31, 46, 56):
            full_pos = (triplet - 1) * 3 + 1
            jobs.append(ScreenJob(
                kind="control",
                label=f"control_leading_G_{mutation_label(full_pos)}",
                center_full_position=full_pos,
                center_global_triplet=triplet,
                center_register="G",
                mutation_full_positions=(full_pos,),
                construct="",
                samples=args.single_samples,
                experimental_site=False,
                control_type="triplet_leading_G_to_C",
                window_center_label=f"control_triplet_{triplet}",
            ))

    if args.mode in ("construct", "both"):
        for construct in args.constructs:
            if construct not in CONSTRUCT_CYS_SITES:
                raise ValueError(f"Unknown construct: {construct}")
            for center_triplet, local_sites in maximal_construct_clusters(
                construct,
                args.construct_mutation_radius_triplets,
            ):
                center_pos = min(
                    local_sites,
                    key=lambda position: abs(
                        triplet_info(position)[0] - center_triplet
                    ),
                )
                _, register = triplet_info(center_pos)
                mutation_names = "_".join(mutation_label(p) for p in local_sites)
                jobs.append(ScreenJob(
                    kind="construct",
                    label=f"{construct}_cluster_{mutation_names}",
                    center_full_position=center_pos,
                    center_global_triplet=center_triplet,
                    center_register=register,
                    mutation_full_positions=local_sites,
                    construct=construct,
                    samples=args.construct_samples,
                    experimental_site=True,
                    window_center_label=f"cluster_triplet_{center_triplet}",
                ))

    return jobs


# ---------------------------------------------------------------------------
# PyRosetta initialization and pose setup
# ---------------------------------------------------------------------------

def initialize_pyrosetta(args: ScreenConfig) -> str:
    global pyrosetta, Pose, AtomID, AtomPairConstraint, CoordinateConstraint
    global HarmonicFunc, FastRelax, MoveMap, HBondSet, fill_hbond_set
    global SasaCalc, mutate_residue, pose_from_pdb, create_score_function
    global rosetta_random, scoring_module, SCORE_TYPES

    try:
        import pyrosetta as _pyrosetta
        from pyrosetta import create_score_function as _create_score_function
        from pyrosetta import init as _init
        from pyrosetta import pose_from_pdb as _pose_from_pdb
        from pyrosetta.toolbox import mutate_residue as _mutate_residue
        from pyrosetta.rosetta.core.id import AtomID as _AtomID
        from pyrosetta.rosetta.core.kinematics import MoveMap as _MoveMap
        from pyrosetta.rosetta.core.pose import Pose as _Pose
        from pyrosetta.rosetta.core.scoring import constraints as _constraints
        from pyrosetta.rosetta.core.scoring import func as _func
        from pyrosetta.rosetta.core.scoring import hbonds as _hbonds
        from pyrosetta.rosetta.core.scoring import sasa as _sasa
        from pyrosetta.rosetta.core import scoring as _scoring
        from pyrosetta.rosetta.numeric import random as _rosetta_random
        from pyrosetta.rosetta.protocols.relax import FastRelax as _FastRelax
    except ImportError as exc:
        raise RuntimeError(
            "PyRosetta is required for screening. Install/activate PyRosetta, or set "
            "NOTEBOOK_CONFIG.self_test = True to validate sequence/window logic only."
        ) from exc

    pyrosetta = _pyrosetta
    Pose = _Pose
    AtomID = _AtomID
    AtomPairConstraint = _constraints.AtomPairConstraint
    CoordinateConstraint = _constraints.CoordinateConstraint
    HarmonicFunc = _func.HarmonicFunc
    FastRelax = _FastRelax
    MoveMap = _MoveMap
    HBondSet = _hbonds.HBondSet
    fill_hbond_set = _hbonds.fill_hbond_set
    SasaCalc = _sasa.SasaCalc
    mutate_residue = _mutate_residue
    pose_from_pdb = _pose_from_pdb
    create_score_function = _create_score_function
    rosetta_random = _rosetta_random
    scoring_module = _scoring

    flags = (
        f"-ignore_unrecognized_res -ignore_zero_occupancy false -mute all "
        f"-constant_seed -jran {args.seed}"
    )
    _init(flags)

    SCORE_TYPES = {
        name: getattr(scoring_module, name)
        for name in TRACKED_SCORE_TERM_NAMES
        if hasattr(scoring_module, name)
    }
    SCORE_TYPES["atom_pair_constraint"] = scoring_module.atom_pair_constraint
    SCORE_TYPES["coordinate_constraint"] = scoring_module.coordinate_constraint
    return pyrosetta.version()


def clone_pose(pose: Any) -> Any:
    cloned = Pose()
    cloned.assign(pose)
    return cloned


def set_rosetta_seed(seed: int) -> None:
    rosetta_random.rg().set_seed("mt19937", int(seed))


def stable_seed(base_seed: int, label: str, replicate: int) -> int:
    checksum = zlib.crc32(label.encode("utf-8")) & 0x7FFFFFFF
    seed = (base_seed + checksum + replicate * 1009) % 2147483646
    return seed + 1


def clear_pose_constraints(pose: Any) -> None:
    try:
        pose.remove_constraints()
    except Exception as exc:
        raise RuntimeError("Could not remove constraints from pose") from exc


def check_template_pose(pose: Any) -> None:
    expected = WINDOW_CHAIN_LENGTH * N_CHAINS
    if pose.size() != expected:
        raise ValueError(f"Template has {pose.size()} residues; expected {expected}")
    if pose.num_chains() != N_CHAINS:
        raise ValueError(f"Template has {pose.num_chains()} chains; expected 3")
    for chain in range(1, N_CHAINS + 1):
        length = pose.chain_end(chain) - pose.chain_begin(chain) + 1
        if length != WINDOW_CHAIN_LENGTH:
            raise ValueError(f"Template chain {chain} length is {length}; expected 63")
        for local_pos in range(1, WINDOW_CHAIN_LENGTH + 1, 3):
            pose_pos = pose.chain_begin(chain) + local_pos - 1
            if pose.residue(pose_pos).name1() != "G":
                raise ValueError(
                    f"Template chain {chain}, local position {local_pos} is not Gly"
                )


def pose_position(chain: int, local_position: int) -> int:
    return (chain - 1) * WINDOW_CHAIN_LENGTH + local_position


def mutate_pose_to_window_sequence(template_pose: Any, target_sequence: str) -> Any:
    if len(target_sequence) != WINDOW_CHAIN_LENGTH:
        raise ValueError("Target window sequence must contain 63 residues")
    pose = clone_pose(template_pose)
    clear_pose_constraints(pose)
    for local_pos, aa in enumerate(target_sequence, start=1):
        for chain in range(1, N_CHAINS + 1):
            pos = pose_position(chain, local_pos)
            if pose.residue(pos).name1() != aa:
                mutate_residue(pose, pos, aa, pack_radius=0.0)
    return pose


def apply_full_position_mutations(
    pose: Any,
    spec: WindowSpec,
    full_positions: Sequence[int],
) -> None:
    for full_pos in full_positions:
        if full_pos not in spec.full_to_local:
            raise ValueError(
                f"Mutation {mutation_label(full_pos)} is outside centered window "
                f"for triplet {spec.center_global_triplet}"
            )
        local_pos = spec.full_to_local[full_pos]
        for chain in range(1, N_CHAINS + 1):
            pos = pose_position(chain, local_pos)
            if pose.residue(pos).name1() != "C":
                mutate_residue(pose, pos, "C", pack_radius=0.0)


def ca_distance(pose: Any, r1: int, r2: int) -> float:
    return pose.residue(r1).xyz("CA").distance(pose.residue(r2).xyz("CA"))


def add_anchor_constraints(pose: Any, reference_pose: Any, args: ScreenConfig) -> None:
    anchor_atom = AtomID(pose.residue(1).atom_index("CA"), 1)

    for chain in range(1, N_CHAINS + 1):
        for local_pos in range(1, WINDOW_CHAIN_LENGTH + 1):
            triplet = local_triplet_for_position(local_pos)
            if triplet in OUTER_ANCHOR_TRIPLETS:
                sd = args.anchor_coordinate_sd
            elif triplet in TRANSITION_TRIPLETS:
                sd = args.transition_coordinate_sd
            else:
                continue
            pos = pose_position(chain, local_pos)
            atom = AtomID(pose.residue(pos).atom_index("CA"), pos)
            target = reference_pose.residue(pos).xyz("CA")
            if pos != 1:
                pose.add_constraint(
                    CoordinateConstraint(atom, anchor_atom, target, HarmonicFunc(0.0, sd))
                )

    for local_pos in range(1, WINDOW_CHAIN_LENGTH + 1, 3):
        triplet = local_triplet_for_position(local_pos)
        if triplet not in OUTER_ANCHOR_TRIPLETS:
            continue
        chain_positions = [pose_position(chain, local_pos) for chain in range(1, 4)]
        for r1, r2 in (
            (chain_positions[0], chain_positions[1]),
            (chain_positions[1], chain_positions[2]),
            (chain_positions[0], chain_positions[2]),
        ):
            atom1 = AtomID(pose.residue(r1).atom_index("CA"), r1)
            atom2 = AtomID(pose.residue(r2).atom_index("CA"), r2)
            target_distance = ca_distance(reference_pose, r1, r2)
            pose.add_constraint(
                AtomPairConstraint(
                    atom1,
                    atom2,
                    HarmonicFunc(target_distance, args.gly_anchor_distance_sd),
                )
            )


def make_window_movemap() -> Any:
    movemap = MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(True)
    movemap.set_jump(False)
    for chain in range(1, N_CHAINS + 1):
        for local_pos in range(1, WINDOW_CHAIN_LENGTH + 1):
            triplet = local_triplet_for_position(local_pos)
            pos = pose_position(chain, local_pos)
            movemap.set_bb(pos, triplet in MOVABLE_TRIPLETS)
            movemap.set_chi(pos, True)
    return movemap


def relax_window_pose(
    pose: Any,
    reference_pose: Any,
    scorefxn_cst: Any,
    args: ScreenConfig,
    seed: int,
    repeats: int,
) -> Any:
    clear_pose_constraints(pose)
    add_anchor_constraints(pose, reference_pose, args)
    set_rosetta_seed(seed)

    relax = FastRelax(scorefxn_cst, repeats)
    relax.set_movemap(make_window_movemap())
    relax.constrain_relax_to_start_coords(False)
    relax.ramp_down_constraints(False)
    relax.apply(pose)
    return pose


# ---------------------------------------------------------------------------
# Structural and energetic measurements
# ---------------------------------------------------------------------------

def pose_indices_for_local_positions(local_positions: Iterable[int]) -> List[int]:
    indices: List[int] = []
    for chain in range(1, N_CHAINS + 1):
        indices.extend(pose_position(chain, p) for p in sorted(set(local_positions)))
    return indices


def local_positions_for_triplets(triplets: Iterable[int]) -> List[int]:
    positions: List[int] = []
    for triplet in sorted(set(triplets)):
        if 1 <= triplet <= WINDOW_TRIPLETS:
            start = (triplet - 1) * 3 + 1
            positions.extend((start, start + 1, start + 2))
    return positions


def mutation_context_triplets(
    spec: WindowSpec,
    mutation_full_positions: Sequence[int],
    radius: int,
) -> List[int]:
    triplets = set()
    for full_pos in mutation_full_positions:
        local_pos = spec.full_to_local[full_pos]
        center = local_triplet_for_position(local_pos)
        for triplet in range(center - radius, center + radius + 1):
            if 1 <= triplet <= WINDOW_TRIPLETS:
                triplets.add(triplet)
    return sorted(triplets)


def selected_cys_pose_indices(
    spec: WindowSpec,
    mutation_full_positions: Sequence[int],
) -> List[int]:
    local_positions = [spec.full_to_local[p] for p in mutation_full_positions]
    return pose_indices_for_local_positions(local_positions)


def vector_to_tuple(vector: Any) -> Tuple[float, float, float]:
    return float(vector.x), float(vector.y), float(vector.z)


def aligned_ca_rmsd(
    reference_pose: Any,
    mobile_pose: Any,
    measure_indices: Sequence[int],
) -> float:
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("numpy is required for aligned RMSD calculations") from exc

    anchor_local = local_positions_for_triplets(OUTER_ANCHOR_TRIPLETS)
    anchor_indices = pose_indices_for_local_positions(anchor_local)

    ref_anchor = np.array([
        vector_to_tuple(reference_pose.residue(r).xyz("CA")) for r in anchor_indices
    ])
    mob_anchor = np.array([
        vector_to_tuple(mobile_pose.residue(r).xyz("CA")) for r in anchor_indices
    ])
    ref_center = ref_anchor.mean(axis=0)
    mob_center = mob_anchor.mean(axis=0)
    ref_zero = ref_anchor - ref_center
    mob_zero = mob_anchor - mob_center
    covariance = mob_zero.T @ ref_zero
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    ref_measure = np.array([
        vector_to_tuple(reference_pose.residue(r).xyz("CA")) for r in measure_indices
    ])
    mob_measure = np.array([
        vector_to_tuple(mobile_pose.residue(r).xyz("CA")) for r in measure_indices
    ])
    aligned = (mob_measure - mob_center) @ rotation + ref_center
    squared = ((aligned - ref_measure) ** 2).sum(axis=1)
    return float(math.sqrt(float(squared.mean())))


def triangle_area(p1: Any, p2: Any, p3: Any) -> float:
    a = p1.distance(p2)
    b = p2.distance(p3)
    c = p1.distance(p3)
    semiperimeter = 0.5 * (a + b + c)
    area_sq = max(
        semiperimeter
        * (semiperimeter - a)
        * (semiperimeter - b)
        * (semiperimeter - c),
        0.0,
    )
    return math.sqrt(area_sq)


def gly_triangle_changes(
    wt_pose: Any,
    mut_pose: Any,
    spec: WindowSpec,
    mutation_full_positions: Sequence[int],
    radius_triplets: int,
) -> Dict[str, float]:
    triplets = mutation_context_triplets(
        spec, mutation_full_positions, radius_triplets
    )
    deltas: List[float] = []
    distance_deltas: List[float] = []

    for triplet in triplets:
        gly_local = (triplet - 1) * 3 + 1
        residues = [pose_position(chain, gly_local) for chain in range(1, 4)]
        wt_points = [wt_pose.residue(r).xyz("CA") for r in residues]
        mut_points = [mut_pose.residue(r).xyz("CA") for r in residues]
        deltas.append(
            triangle_area(*mut_points) - triangle_area(*wt_points)
        )
        for a, b in ((0, 1), (1, 2), (0, 2)):
            wt_distance = wt_points[a].distance(wt_points[b])
            mut_distance = mut_points[a].distance(mut_points[b])
            distance_deltas.append(abs(mut_distance - wt_distance))

    abs_deltas = [abs(value) for value in deltas]
    return {
        "gly_triangle_count": len(deltas),
        "gly_triangle_mean_delta_A2": statistics.fmean(deltas) if deltas else 0.0,
        "gly_triangle_mean_abs_delta_A2": statistics.fmean(abs_deltas) if abs_deltas else 0.0,
        "gly_triangle_max_abs_delta_A2": max(abs_deltas) if abs_deltas else 0.0,
        "gly_interchain_ca_max_abs_delta_A": max(distance_deltas) if distance_deltas else 0.0,
    }


def score_breakdown(scorefxn: Any, pose: Any) -> Dict[str, float]:
    total = float(scorefxn(pose))
    energies = pose.energies().total_energies()
    result = {"total_score": total}
    for name, score_type in SCORE_TYPES.items():
        if name in ("atom_pair_constraint", "coordinate_constraint"):
            continue
        weight = float(scorefxn.get_weight(score_type))
        result[name] = float(energies[score_type]) * weight
    return result


def constraint_energy(scorefxn_cst: Any, pose: Any) -> Dict[str, float]:
    scorefxn_cst(pose)
    energies = pose.energies().total_energies()
    atom_type = SCORE_TYPES["atom_pair_constraint"]
    coord_type = SCORE_TYPES["coordinate_constraint"]
    atom = float(energies[atom_type]) * float(scorefxn_cst.get_weight(atom_type))
    coord = float(energies[coord_type]) * float(scorefxn_cst.get_weight(coord_type))
    return {
        "atom_pair_constraint_energy": atom,
        "coordinate_constraint_energy": coord,
        "total_constraint_energy": atom + coord,
    }


def selected_residue_total_energy(
    pose: Any,
    scorefxn_raw: Any,
    residue_indices: Sequence[int],
) -> float:
    scorefxn_raw(pose)
    return float(sum(pose.energies().residue_total_energy(r) for r in residue_indices))


def mutation_crosschain_pair_energy(
    pose: Any,
    scorefxn_raw: Any,
    selected_indices: Sequence[int],
) -> Dict[str, float]:
    scorefxn_raw(pose)
    selected = set(selected_indices)
    graph = pose.energies().energy_graph()
    weights = scorefxn_raw.weights()
    total = 0.0
    term_totals = defaultdict(float)

    for r1 in range(1, pose.size() + 1):
        for r2 in range(r1 + 1, pose.size() + 1):
            if r1 not in selected and r2 not in selected:
                continue
            if pose.chain(r1) == pose.chain(r2):
                continue
            edge = graph.find_energy_edge(r1, r2)
            if edge is None:
                continue
            energy_map = edge.fill_energy_map()
            total += float(energy_map.dot(weights))
            for name, score_type in SCORE_TYPES.items():
                if name in ("atom_pair_constraint", "coordinate_constraint"):
                    continue
                term_totals[name] += (
                    float(energy_map[score_type])
                    * float(scorefxn_raw.get_weight(score_type))
                )

    result = {"mutation_crosschain_interaction": total}
    for name, value in term_totals.items():
        result[f"mutation_crosschain_{name}"] = value
    return result


def atom_is_backbone(pose: Any, residue: int, atom_index: int) -> bool:
    try:
        return bool(pose.residue(residue).atom_is_backbone(atom_index))
    except Exception:
        return False


def hbond_category(pose: Any, hbond: Any) -> str:
    donor_backbone = atom_is_backbone(pose, hbond.don_res(), hbond.don_hatm())
    acceptor_backbone = atom_is_backbone(pose, hbond.acc_res(), hbond.acc_atm())
    if donor_backbone and acceptor_backbone:
        return "bb_bb"
    if donor_backbone or acceptor_backbone:
        return "bb_sc"
    return "sc_sc"


def local_interchain_hbonds(
    pose: Any,
    scorefxn_raw: Any,
    local_residue_indices: Sequence[int],
    energy_cutoff: float,
) -> Dict[str, float]:
    scorefxn_raw(pose)
    pose.update_residue_neighbors()
    hbonds = HBondSet()
    fill_hbond_set(pose, False, hbonds)
    local_set = set(local_residue_indices)
    counts = defaultdict(int)
    energies = defaultdict(float)

    for index in range(1, hbonds.nhbonds() + 1):
        hbond = hbonds.hbond(index)
        energy = float(hbond.energy())
        if energy > energy_cutoff:
            continue
        donor = hbond.don_res()
        acceptor = hbond.acc_res()
        if pose.chain(donor) == pose.chain(acceptor):
            continue
        if donor not in local_set and acceptor not in local_set:
            continue
        category = hbond_category(pose, hbond)
        counts[category] += 1
        energies[category] += energy

    result: Dict[str, float] = {}
    for category in ("bb_bb", "bb_sc", "sc_sc"):
        result[f"local_interchain_hbond_{category}_count"] = float(counts[category])
        result[f"local_interchain_hbond_{category}_energy"] = float(energies[category])
    return result


def cys_sasa_metrics(pose: Any, selected_indices: Sequence[int]) -> Dict[str, float]:
    try:
        calculator = SasaCalc()
        calculator.calculate(pose)
        residue_sc_sasa = calculator.get_residue_sasa_sc()
        atom_sasa = calculator.get_atom_sasa()
        sidechain_values: List[float] = []
        sg_values: List[float] = []
        for residue in selected_indices:
            sidechain_values.append(float(residue_sc_sasa[residue]))
            if pose.residue(residue).has("SG"):
                atom_id = AtomID(pose.residue(residue).atom_index("SG"), residue)
                sg_values.append(float(atom_sasa[atom_id]))
        return {
            "mut_cys_sidechain_sasa_mean_A2": statistics.fmean(sidechain_values)
            if sidechain_values else math.nan,
            "mut_cys_sidechain_sasa_min_A2": min(sidechain_values)
            if sidechain_values else math.nan,
            "mut_cys_sg_sasa_mean_A2": statistics.fmean(sg_values)
            if sg_values else math.nan,
            "mut_cys_sg_sasa_min_A2": min(sg_values) if sg_values else math.nan,
        }
    except Exception as exc:
        print(f"WARNING: SASA calculation failed: {exc}", file=sys.stderr)
        return {
            "mut_cys_sidechain_sasa_mean_A2": math.nan,
            "mut_cys_sidechain_sasa_min_A2": math.nan,
            "mut_cys_sg_sasa_mean_A2": math.nan,
            "mut_cys_sg_sasa_min_A2": math.nan,
        }


def cys_geometry_metrics(
    pose: Any,
    spec: WindowSpec,
    mutation_full_positions: Sequence[int],
) -> Dict[str, float]:
    sg_distances: List[float] = []
    cb_distances: List[float] = []
    for full_pos in mutation_full_positions:
        local_pos = spec.full_to_local[full_pos]
        residues = [pose_position(chain, local_pos) for chain in range(1, 4)]
        for a, b in ((0, 1), (1, 2), (0, 2)):
            r1, r2 = residues[a], residues[b]
            if pose.residue(r1).has("SG") and pose.residue(r2).has("SG"):
                sg_distances.append(
                    pose.residue(r1).xyz("SG").distance(pose.residue(r2).xyz("SG"))
                )
            atom1 = "CB" if pose.residue(r1).has("CB") else "CA"
            atom2 = "CB" if pose.residue(r2).has("CB") else "CA"
            cb_distances.append(
                pose.residue(r1).xyz(atom1).distance(pose.residue(r2).xyz(atom2))
            )
    return {
        "mut_cys_sg_distance_mean_A": statistics.fmean(sg_distances)
        if sg_distances else math.nan,
        "mut_cys_sg_distance_min_A": min(sg_distances) if sg_distances else math.nan,
        "mut_cys_sg_distance_max_A": max(sg_distances) if sg_distances else math.nan,
        "mut_cys_cb_distance_mean_A": statistics.fmean(cb_distances)
        if cb_distances else math.nan,
        "mut_cys_cb_distance_min_A": min(cb_distances) if cb_distances else math.nan,
        "mut_cys_pair_count_sg_le_6_5A": float(
            sum(distance <= 6.5 for distance in sg_distances)
        ),
    }


def compare_poses(
    wt_pose: Any,
    mut_pose: Any,
    reference_pose: Any,
    spec: WindowSpec,
    mutation_full_positions: Sequence[int],
    scorefxn_raw: Any,
    scorefxn_cst: Any,
    args: ScreenConfig,
) -> Dict[str, float]:
    wt_scores = score_breakdown(scorefxn_raw, wt_pose)
    mut_scores = score_breakdown(scorefxn_raw, mut_pose)
    wt_constraints = constraint_energy(scorefxn_cst, wt_pose)
    mut_constraints = constraint_energy(scorefxn_cst, mut_pose)

    context_triplets = mutation_context_triplets(
        spec,
        mutation_full_positions,
        args.mutation_context_radius_triplets,
    )
    context_local = local_positions_for_triplets(context_triplets)
    context_indices = pose_indices_for_local_positions(context_local)
    selected_indices = selected_cys_pose_indices(spec, mutation_full_positions)

    wt_local_energy = selected_residue_total_energy(
        wt_pose, scorefxn_raw, context_indices
    )
    mut_local_energy = selected_residue_total_energy(
        mut_pose, scorefxn_raw, context_indices
    )

    wt_pair = mutation_crosschain_pair_energy(
        wt_pose, scorefxn_raw, selected_indices
    )
    mut_pair = mutation_crosschain_pair_energy(
        mut_pose, scorefxn_raw, selected_indices
    )

    wt_hbond = local_interchain_hbonds(
        wt_pose,
        scorefxn_raw,
        context_indices,
        args.hbond_energy_cutoff,
    )
    mut_hbond = local_interchain_hbonds(
        mut_pose,
        scorefxn_raw,
        context_indices,
        args.hbond_energy_cutoff,
    )

    result: Dict[str, float] = {
        "mutation_site_count_in_window": float(len(mutation_full_positions)),
        "delta_total_score": mut_scores["total_score"] - wt_scores["total_score"],
        "delta_total_score_per_site": (
            mut_scores["total_score"] - wt_scores["total_score"]
        ) / len(mutation_full_positions),
        "wt_total_score": wt_scores["total_score"],
        "mut_total_score": mut_scores["total_score"],
        "wt_mutation_context_residue_energy": wt_local_energy,
        "mut_mutation_context_residue_energy": mut_local_energy,
        "delta_mutation_context_residue_energy": mut_local_energy - wt_local_energy,
        "delta_mutation_context_residue_energy_per_site": (
            mut_local_energy - wt_local_energy
        ) / len(mutation_full_positions),
        "mutation_context_ca_rmsd_A": aligned_ca_rmsd(
            wt_pose, mut_pose, context_indices
        ),
        "analysis_core_ca_rmsd_A": aligned_ca_rmsd(
            wt_pose,
            mut_pose,
            pose_indices_for_local_positions(
                local_positions_for_triplets(ANALYSIS_CORE_TRIPLETS)
            ),
        ),
    }

    for name in SCORE_TYPES:
        if name in ("atom_pair_constraint", "coordinate_constraint"):
            continue
        if name in wt_scores and name in mut_scores:
            result[f"wt_{name}"] = wt_scores[name]
            result[f"mut_{name}"] = mut_scores[name]
            result[f"delta_{name}"] = mut_scores[name] - wt_scores[name]

    for name in wt_constraints:
        result[f"wt_{name}"] = wt_constraints[name]
        result[f"mut_{name}"] = mut_constraints[name]
        result[f"delta_{name}"] = mut_constraints[name] - wt_constraints[name]

    result["wt_mutation_crosschain_interaction"] = wt_pair[
        "mutation_crosschain_interaction"
    ]
    result["mut_mutation_crosschain_interaction"] = mut_pair[
        "mutation_crosschain_interaction"
    ]
    result["delta_mutation_crosschain_interaction"] = (
        mut_pair["mutation_crosschain_interaction"]
        - wt_pair["mutation_crosschain_interaction"]
    )
    for key in sorted(set(wt_pair) | set(mut_pair)):
        if key == "mutation_crosschain_interaction":
            continue
        result[f"wt_{key}"] = wt_pair.get(key, 0.0)
        result[f"mut_{key}"] = mut_pair.get(key, 0.0)
        result[f"delta_{key}"] = mut_pair.get(key, 0.0) - wt_pair.get(key, 0.0)

    for key in wt_hbond:
        result[f"wt_{key}"] = wt_hbond[key]
        result[f"mut_{key}"] = mut_hbond[key]
        result[f"delta_{key}"] = mut_hbond[key] - wt_hbond[key]

    result.update(
        gly_triangle_changes(
            wt_pose,
            mut_pose,
            spec,
            mutation_full_positions,
            args.gly_analysis_radius_triplets,
        )
    )
    result.update(cys_sasa_metrics(mut_pose, selected_indices))
    result.update(cys_geometry_metrics(mut_pose, spec, mutation_full_positions))
    return result


# ---------------------------------------------------------------------------
# Replicate execution and aggregation
# ---------------------------------------------------------------------------

def job_metadata(job: ScreenJob, spec: WindowSpec) -> Dict[str, Any]:
    center_local = spec.full_to_local[job.center_full_position]
    reference_triplet, _ = triplet_info(job.center_full_position)
    all_construct_sites = CONSTRUCT_CYS_SITES.get(job.construct, [])
    omitted = sorted(set(all_construct_sites) - set(job.mutation_full_positions))
    return {
        "Job_Type": job.kind,
        "Job_Label": job.label,
        "Construct": job.construct,
        "Control_Type": job.control_type,
        "Window_Center_Label": job.window_center_label,
        "Center_Candidate": mutation_label(job.center_full_position),
        "Center_Full_Position": job.center_full_position,
        "Center_Global_Triplet": job.center_global_triplet,
        "Reference_Site_Global_Triplet": reference_triplet,
        "Center_Register": job.center_register,
        "Center_Local_Position": center_local,
        "Mutation_Full_Positions": ";".join(map(str, job.mutation_full_positions)),
        "Mutation_Labels": ";".join(mutation_label(p) for p in job.mutation_full_positions),
        "Mutation_Count_In_Window": len(job.mutation_full_positions),
        "Construct_Sites_Omitted_From_Local_Window": ";".join(map(str, omitted)),
        "Experimental_Site": "yes" if job.experimental_site else "no",
        "Left_GPP_Padding_Triplets": spec.left_padding_triplets,
        "Right_GPP_Padding_Triplets": spec.right_padding_triplets,
        "Context_Class": "padded_intrinsic" if (
            spec.left_padding_triplets or spec.right_padding_triplets
        ) else "fully_native_window",
        "Parent_Window_Sequence": spec.parent_sequence,
    }


def numeric_values(rows: Sequence[Dict[str, Any]], key: str) -> List[float]:
    values: List[float] = []
    for row in rows:
        value = row.get(key)
        if isinstance(value, (int, float)) and not math.isnan(float(value)):
            values.append(float(value))
    return values


SUMMARY_METRICS = (
    "delta_total_score",
    "delta_total_score_per_site",
    "delta_mutation_context_residue_energy",
    "delta_mutation_context_residue_energy_per_site",
    "delta_mutation_crosschain_interaction",
    "delta_fa_atr",
    "delta_fa_rep",
    "delta_fa_sol",
    "delta_fa_elec",
    "delta_fa_dun",
    "delta_hbond_sr_bb",
    "delta_hbond_lr_bb",
    "delta_hbond_bb_sc",
    "delta_hbond_sc",
    "delta_total_constraint_energy",
    "mutation_context_ca_rmsd_A",
    "analysis_core_ca_rmsd_A",
    "gly_triangle_mean_abs_delta_A2",
    "gly_triangle_max_abs_delta_A2",
    "gly_interchain_ca_max_abs_delta_A",
    "mut_cys_sidechain_sasa_mean_A2",
    "mut_cys_sidechain_sasa_min_A2",
    "mut_cys_sg_sasa_mean_A2",
    "mut_cys_sg_sasa_min_A2",
    "mut_cys_sg_distance_mean_A",
    "mut_cys_sg_distance_min_A",
    "mut_cys_sg_distance_max_A",
    "mut_cys_cb_distance_min_A",
    "mut_cys_pair_count_sg_le_6_5A",
    "delta_local_interchain_hbond_bb_bb_count",
    "delta_local_interchain_hbond_bb_bb_energy",
    "delta_local_interchain_hbond_bb_sc_count",
    "delta_local_interchain_hbond_bb_sc_energy",
    "delta_local_interchain_hbond_sc_sc_count",
    "delta_local_interchain_hbond_sc_sc_energy",
    "nonadditivity_delta_total_score",
    "nonadditivity_delta_context_energy",
    "nonadditivity_delta_crosschain_interaction",
)


def summarize_job_rows(rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    if not rows:
        raise ValueError("Cannot summarize an empty replicate group")
    first = rows[0]
    metadata_keys = (
        "Job_Type",
        "Job_Label",
        "Construct",
        "Control_Type",
        "Window_Center_Label",
        "Center_Candidate",
        "Center_Full_Position",
        "Center_Global_Triplet",
        "Reference_Site_Global_Triplet",
        "Center_Register",
        "Center_Local_Position",
        "Mutation_Full_Positions",
        "Mutation_Labels",
        "Mutation_Count_In_Window",
        "Construct_Sites_Omitted_From_Local_Window",
        "Experimental_Site",
        "Left_GPP_Padding_Triplets",
        "Right_GPP_Padding_Triplets",
        "Context_Class",
        "Parent_Window_Sequence",
    )
    summary = {key: first.get(key, "") for key in metadata_keys}
    summary["Replicate_N"] = len(rows)

    for metric in SUMMARY_METRICS:
        values = numeric_values(rows, metric)
        if not values:
            continue
        summary[f"{metric}_mean"] = statistics.fmean(values)
        summary[f"{metric}_sd"] = statistics.stdev(values) if len(values) > 1 else 0.0
        summary[f"{metric}_min"] = min(values)
        summary[f"{metric}_max"] = max(values)
    return summary


def permission_status(summary: Dict[str, Any], args: ScreenConfig) -> Tuple[str, str]:
    flags: List[str] = []
    if summary.get("mutation_context_ca_rmsd_A_mean", 0.0) > args.rmsd_fail_A:
        flags.append("local_backbone_rmsd")
    if summary.get("gly_triangle_max_abs_delta_A2_mean", 0.0) > args.gly_area_fail_A2:
        flags.append("gly_triangle_geometry")
    if summary.get("delta_fa_rep_mean", 0.0) > args.fa_rep_fail_REU:
        flags.append("steric_repulsion")
    return ("FAIL", ";".join(flags)) if flags else ("PASS", "none")


def percentile_risks(rows: List[Dict[str, Any]], metric: str, output: str) -> None:
    valid = []
    for index, row in enumerate(rows):
        value = row.get(metric)
        if isinstance(value, (int, float)) and not math.isnan(float(value)):
            valid.append((float(value), index))
    if not valid:
        return
    values = sorted(value for value, _ in valid)
    denominator = max(len(values) - 1, 1)
    for value, index in valid:
        lower = sum(v < value for v in values)
        equal = sum(v == value for v in values)
        average_rank = lower + 0.5 * (equal - 1)
        rows[index][output] = average_rank / denominator


def mean_available(row: Dict[str, Any], keys: Sequence[str]) -> float:
    values = [
        float(row[key])
        for key in keys
        if key in row and isinstance(row[key], (int, float)) and not math.isnan(float(row[key]))
    ]
    return statistics.fmean(values) if values else math.nan


def add_family_scores(summaries: List[Dict[str, Any]]) -> None:
    risk_metrics = {
        "Risk_delta_total_per_site": "delta_total_score_per_site_mean",
        "Risk_delta_crosschain": "delta_mutation_crosschain_interaction_mean",
        "Risk_local_rmsd": "mutation_context_ca_rmsd_A_mean",
        "Risk_gly_geometry": "gly_triangle_max_abs_delta_A2_mean",
        "Risk_fa_rep": "delta_fa_rep_mean",
        "Risk_low_sg_sasa": "_negative_sg_sasa",
        "Risk_hbond_bb_energy": "delta_local_interchain_hbond_bb_bb_energy_mean",
        "Risk_nonadditivity": "nonadditivity_delta_total_score_mean",
    }
    for row in summaries:
        sasa = row.get("mut_cys_sg_sasa_mean_A2_mean", math.nan)
        row["_negative_sg_sasa"] = -float(sasa) if isinstance(sasa, (int, float)) else math.nan
    for output, metric in risk_metrics.items():
        percentile_risks(summaries, metric, output)

    for row in summaries:
        row["Stability_Family_Risk"] = mean_available(
            row,
            (
                "Risk_delta_total_per_site",
                "Risk_delta_crosschain",
                *(() if row.get("Job_Type") != "construct" else ("Risk_nonadditivity",)),
            ),
        )
        row["Scaffold_Family_Risk"] = mean_available(
            row, ("Risk_local_rmsd", "Risk_gly_geometry", "Risk_fa_rep")
        )
        row["Accessibility_Family_Risk"] = mean_available(
            row, ("Risk_low_sg_sasa",)
        )
        primary = (
            "Stability_Family_Risk",
            "Scaffold_Family_Risk",
            "Accessibility_Family_Risk",
        )
        row["Overall_Permission_Risk"] = mean_available(row, primary)
        row.pop("_negative_sg_sasa", None)

    ranked = sorted(
        summaries,
        key=lambda row: (
            math.inf if math.isnan(row.get("Overall_Permission_Risk", math.nan))
            else row["Overall_Permission_Risk"]
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["Overall_Rank"] = rank
        row["Overall_N"] = len(ranked)


def run_screening(
    jobs: Sequence[ScreenJob],
    template_pose: Any,
    scorefxn_raw: Any,
    scorefxn_cst: Any,
    args: ScreenConfig,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Tuple[Any, Any]]]:
    jobs_by_triplet: Dict[int, List[ScreenJob]] = defaultdict(list)
    for job in jobs:
        jobs_by_triplet[job.center_global_triplet].append(job)

    replicate_rows: List[Dict[str, Any]] = []
    component_rows: List[Dict[str, Any]] = []
    best_poses: Dict[str, Tuple[Any, Any]] = {}
    best_scores: Dict[str, float] = {}

    for center_triplet in sorted(jobs_by_triplet):
        group = jobs_by_triplet[center_triplet]
        spec = build_centered_window(center_triplet)
        print(
            f"\nCenter triplet {center_triplet:02d}/{DOMAIN_TRIPLETS} | "
            f"jobs={len(group)} | padding={spec.left_padding_triplets}/"
            f"{spec.right_padding_triplets}"
        )

        background = mutate_pose_to_window_sequence(template_pose, spec.parent_sequence)
        background_seed = stable_seed(args.seed, f"background_{center_triplet}", 0)
        background = relax_window_pose(
            background,
            template_pose,
            scorefxn_cst,
            args,
            background_seed,
            args.background_relax_repeats,
        )
        clear_pose_constraints(background)

        max_samples = max(job.samples for job in group)
        for replicate in range(1, max_samples + 1):
            paired_seed = stable_seed(args.seed, f"center_{center_triplet}", replicate)
            wt_pose = clone_pose(background)
            wt_pose = relax_window_pose(
                wt_pose,
                background,
                scorefxn_cst,
                args,
                paired_seed,
                args.fast_relax_repeats,
            )

            single_variant_cache: Dict[int, Any] = {}
            for job in group:
                if replicate > job.samples:
                    continue
                print(
                    f"  replicate {replicate}/{job.samples} | {job.label} | "
                    f"mutations={','.join(mutation_label(p) for p in job.mutation_full_positions)}"
                )
                mut_pose = clone_pose(background)
                apply_full_position_mutations(mut_pose, spec, job.mutation_full_positions)
                mut_pose = relax_window_pose(
                    mut_pose,
                    background,
                    scorefxn_cst,
                    args,
                    paired_seed,
                    args.fast_relax_repeats,
                )

                metrics = compare_poses(
                    wt_pose,
                    mut_pose,
                    background,
                    spec,
                    job.mutation_full_positions,
                    scorefxn_raw,
                    scorefxn_cst,
                    args,
                )

                if (
                    job.kind == "construct"
                    and args.calculate_nonadditivity
                    and len(job.mutation_full_positions) > 1
                ):
                    sum_total = 0.0
                    sum_context = 0.0
                    sum_pair = 0.0
                    multi_context_triplets = mutation_context_triplets(
                        spec,
                        job.mutation_full_positions,
                        args.mutation_context_radius_triplets,
                    )
                    multi_context_indices = pose_indices_for_local_positions(
                        local_positions_for_triplets(multi_context_triplets)
                    )
                    wt_multi_context_energy = selected_residue_total_energy(
                        wt_pose,
                        scorefxn_raw,
                        multi_context_indices,
                    )
                    for full_pos in job.mutation_full_positions:
                        if full_pos not in single_variant_cache:
                            single_pose = clone_pose(background)
                            apply_full_position_mutations(single_pose, spec, (full_pos,))
                            single_pose = relax_window_pose(
                                single_pose,
                                background,
                                scorefxn_cst,
                                args,
                                paired_seed,
                                args.fast_relax_repeats,
                            )
                            single_variant_cache[full_pos] = single_pose
                        single_metrics = compare_poses(
                            wt_pose,
                            single_variant_cache[full_pos],
                            background,
                            spec,
                            (full_pos,),
                            scorefxn_raw,
                            scorefxn_cst,
                            args,
                        )
                        sum_total += single_metrics["delta_total_score"]
                        single_multi_context_energy = selected_residue_total_energy(
                            single_variant_cache[full_pos],
                            scorefxn_raw,
                            multi_context_indices,
                        )
                        single_delta_in_multi_context = (
                            single_multi_context_energy - wt_multi_context_energy
                        )
                        sum_context += single_delta_in_multi_context
                        sum_pair += single_metrics[
                            "delta_mutation_crosschain_interaction"
                        ]
                        component = job_metadata(job, spec)
                        component.update({
                            "Replicate": replicate,
                            "Seed": paired_seed,
                            "Single_Component_Position": full_pos,
                            "Single_Component_Label": mutation_label(full_pos),
                            "Single_delta_total_score": single_metrics["delta_total_score"],
                            "Single_delta_context_energy": single_metrics[
                                "delta_mutation_context_residue_energy"
                            ],
                            "Single_delta_in_multi_context_energy": (
                                single_delta_in_multi_context
                            ),
                            "Single_delta_crosschain_interaction": single_metrics[
                                "delta_mutation_crosschain_interaction"
                            ],
                        })
                        component_rows.append(component)

                    metrics["nonadditivity_delta_total_score"] = (
                        metrics["delta_total_score"] - sum_total
                    )
                    metrics["nonadditivity_delta_context_energy"] = (
                        metrics["delta_mutation_context_residue_energy"] - sum_context
                    )
                    metrics["nonadditivity_delta_crosschain_interaction"] = (
                        metrics["delta_mutation_crosschain_interaction"] - sum_pair
                    )
                else:
                    metrics["nonadditivity_delta_total_score"] = 0.0
                    metrics["nonadditivity_delta_context_energy"] = 0.0
                    metrics["nonadditivity_delta_crosschain_interaction"] = 0.0

                row = job_metadata(job, spec)
                row.update({"Replicate": replicate, "Seed": paired_seed})
                row.update(metrics)
                replicate_rows.append(row)

                score = metrics["delta_total_score"]
                if args.dump_pdbs and (
                    job.label not in best_scores or score < best_scores[job.label]
                ):
                    best_scores[job.label] = score
                    best_poses[job.label] = (clone_pose(wt_pose), clone_pose(mut_pose))

    return replicate_rows, component_rows, best_poses


def summarize_replicates(
    replicate_rows: Sequence[Dict[str, Any]],
    args: ScreenConfig,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in replicate_rows:
        grouped[str(row["Job_Label"])].append(row)

    all_summaries = [summarize_job_rows(rows) for rows in grouped.values()]
    for summary in all_summaries:
        status, flags = permission_status(summary, args)
        summary["Permission_Status"] = status
        summary["Permission_Flags"] = flags

    single = [row for row in all_summaries if row["Job_Type"] == "single"]
    construct = [row for row in all_summaries if row["Job_Type"] == "construct"]
    controls = [row for row in all_summaries if row["Job_Type"] == "control"]
    add_family_scores(single)
    add_family_scores(construct)
    all_summaries = single + construct + controls
    return all_summaries, construct


def construct_level_summary(
    construct_windows: Sequence[Dict[str, Any]],
    single_summaries: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    single_by_position = {
        int(row["Center_Full_Position"]): row for row in single_summaries
    }
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in construct_windows:
        grouped[str(row["Construct"])].append(row)

    output: List[Dict[str, Any]] = []
    for construct, rows in sorted(grouped.items()):
        sites = CONSTRUCT_CYS_SITES[construct]
        statuses = [
            single_by_position[p].get("Permission_Status", "UNKNOWN")
            for p in sites
            if p in single_by_position
        ]
        single_site_coverage = len(statuses)
        risk_rows = [
            row for row in rows
            if isinstance(row.get("Overall_Permission_Risk"), (int, float))
        ]
        worst_row = max(
            risk_rows,
            key=lambda row: row["Overall_Permission_Risk"],
        ) if risk_rows else None

        def aggregate(metric: str, function: Any) -> float:
            values = [
                float(row[metric])
                for row in rows
                if isinstance(row.get(metric), (int, float))
                and not math.isnan(float(row[metric]))
            ]
            return function(values) if values else math.nan

        output.append({
            "Construct": construct,
            "Construct_Cys_Positions": ";".join(map(str, sites)),
            "Construct_Cys_Labels": ";".join(mutation_label(p) for p in sites),
            "Construct_Cys_Count": len(sites),
            "Window_Count": len(rows),
            "Single_Site_Result_Count": single_site_coverage,
            "Single_Site_PASS_Count": statuses.count("PASS") if statuses else "",
            "Single_Site_FAIL_Count": statuses.count("FAIL") if statuses else "",
            "Effective_Cys_Count": statuses.count("PASS") if statuses else "",
            "Mean_Window_Risk": aggregate("Overall_Permission_Risk", statistics.fmean),
            "Worst_Window_Risk": aggregate("Overall_Permission_Risk", max),
            "Worst_Window_Label": worst_row["Job_Label"] if worst_row else "",
            "Mean_Delta_Total_Per_Site": aggregate(
                "delta_total_score_per_site_mean", statistics.fmean
            ),
            "Worst_Delta_Total_Per_Site": aggregate(
                "delta_total_score_per_site_mean", max
            ),
            "Mean_Nonadditivity_Total": aggregate(
                "nonadditivity_delta_total_score_mean", statistics.fmean
            ),
            "Worst_Nonadditivity_Total": aggregate(
                "nonadditivity_delta_total_score_mean", max
            ),
            "Mean_Mutation_Context_RMSD_A": aggregate(
                "mutation_context_ca_rmsd_A_mean", statistics.fmean
            ),
            "Worst_Mutation_Context_RMSD_A": aggregate(
                "mutation_context_ca_rmsd_A_mean", max
            ),
            "Worst_Gly_Triangle_Delta_A2": aggregate(
                "gly_triangle_max_abs_delta_A2_mean", max
            ),
            "Mean_Cys_SG_SASA_A2": aggregate(
                "mut_cys_sg_sasa_mean_A2_mean", statistics.fmean
            ),
            "Minimum_Cys_SG_SASA_A2": aggregate(
                "mut_cys_sg_sasa_min_A2_mean", min
            ),
        })
    return output


# ---------------------------------------------------------------------------
# Output and command-line interface
# ---------------------------------------------------------------------------

def ensure_directory(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_csv(path: str, rows: Sequence[Dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def dump_best_pdbs(
    output_dir: str,
    best_poses: Dict[str, Tuple[Any, Any]],
) -> None:
    pdb_dir = os.path.join(output_dir, "best_pdbs")
    ensure_directory(pdb_dir)
    for label, (wt_pose, mut_pose) in best_poses.items():
        safe_label = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in label)
        wt_pose.dump_pdb(os.path.join(pdb_dir, f"{safe_label}_WT.pdb"))
        mut_pose.dump_pdb(os.path.join(pdb_dir, f"{safe_label}_MUT.pdb"))


def job_plan_rows(jobs: Sequence[ScreenJob]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for job in jobs:
        spec = build_centered_window(job.center_global_triplet)
        row = job_metadata(job, spec)
        row["Samples"] = job.samples
        rows.append(row)
    return rows


def write_manifest(
    output_dir: str,
    args: ScreenConfig,
    pyrosetta_version: str,
    scorefxn_raw: Any,
) -> None:
    weights = {}
    for name, score_type in SCORE_TYPES.items():
        weights[name] = float(scorefxn_raw.get_weight(score_type))
    manifest = {
        "script": (
            os.path.abspath(globals()["__file__"])
            if "__file__" in globals()
            else "<notebook>"
        ),
        "pyrosetta_version": pyrosetta_version,
        "arguments": vars(args),
        "parent_sequence": FULL_PARENT_COLLAGEN_DOMAIN,
        "construct_cys_sites": CONSTRUCT_CYS_SITES,
        "window_protocol": {
            "window_triplets": WINDOW_TRIPLETS,
            "center_triplet": CENTER_TRIPLET,
            "outer_anchor_triplets": sorted(OUTER_ANCHOR_TRIPLETS),
            "transition_triplets": sorted(TRANSITION_TRIPLETS),
            "analysis_core_triplets": sorted(ANALYSIS_CORE_TRIPLETS),
            "movable_triplets": sorted(MOVABLE_TRIPLETS),
        },
        "score_function_weights": weights,
        "scientific_scope": (
            "Local reduced-Cys permissibility and construct-aware overlapping windows; "
            "not a full-length oxidized/disulfide assembly model."
        ),
    }
    with open(
        os.path.join(output_dir, "run_manifest.json"),
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=True)


def run_self_test() -> None:
    assert len(FULL_PARENT_COLLAGEN_DOMAIN) == DOMAIN_LENGTH
    assert DOMAIN_LENGTH // 3 == DOMAIN_TRIPLETS
    assert "C" not in FULL_PARENT_COLLAGEN_DOMAIN

    p32_triplet, p32_register = triplet_info(32)
    p32_window = build_centered_window(p32_triplet)
    assert p32_triplet == 11 and p32_register == "X"
    assert p32_window.full_to_local[32] == 32
    assert p32_window.left_padding_triplets == 0

    p183_triplet, p183_register = triplet_info(183)
    p183_window = build_centered_window(p183_triplet)
    assert p183_triplet == 61 and p183_register == "Y"
    assert p183_window.full_to_local[183] == 33
    assert p183_window.right_padding_triplets == 10

    expected_parent_residues = {
        32: "P",
        54: "A",
        60: "Q",
        72: "A",
        90: "P",
        105: "V",
        117: "A",
        129: "A",
        150: "Q",
        167: "P",
        183: "P",
    }
    for position, aa in expected_parent_residues.items():
        assert FULL_PARENT_COLLAGEN_DOMAIN[position - 1] == aa

    nc5_a72 = construct_sites_in_local_radius("NC5", 72, 7)
    assert nc5_a72 == (54, 72, 90)
    nc3_clusters = maximal_construct_clusters("NC3", 7)
    assert (16, (32, 60)) in nc3_clusters
    assert (25, (60, 90)) in nc3_clusters
    nc5_clusters = maximal_construct_clusters("NC5", 7)
    assert (18, (32, 54, 72)) in nc5_clusters
    assert (36, (90, 105, 129)) in nc5_clusters
    assert build_centered_window(31).parent_sequence[30:33] == (
        FULL_PARENT_COLLAGEN_DOMAIN[90:93]
    )
    print("Self-test passed: sequence, triplet mapping, padding, and constructs are valid.")


def validate_config(args: ScreenConfig) -> None:
    if args.mode not in {"single", "construct", "both"}:
        raise ValueError("mode must be 'single', 'construct', or 'both'")
    unknown_constructs = sorted(set(args.constructs) - set(CONSTRUCT_CYS_SITES))
    if unknown_constructs:
        raise ValueError(f"Unknown constructs: {', '.join(unknown_constructs)}")
    if args.single_samples < 1 or args.construct_samples < 1:
        raise ValueError("Sample counts must be positive")
    if args.fast_relax_repeats < 1 or args.background_relax_repeats < 1:
        raise ValueError("FastRelax repeat counts must be positive")
    if not 0 <= args.construct_mutation_radius_triplets <= 7:
        raise ValueError(
            "Construct mutation radius must be 0-7 so mutated sites remain outside anchors"
        )
    if args.hbond_energy_cutoff >= 0:
        raise ValueError("H-bond energy cutoff must be negative")


def main(config: Optional[ScreenConfig] = None) -> int:
    args = NOTEBOOK_CONFIG if config is None else config
    validate_config(args)
    run_self_test()
    if args.self_test:
        return 0

    ensure_directory(args.output_dir)
    candidate_rows = build_candidate_library(args)
    jobs = build_jobs(args, candidate_rows)
    if not jobs:
        raise RuntimeError("No screening jobs were generated")
    write_csv(os.path.join(args.output_dir, "candidate_library.csv"), candidate_rows)
    write_csv(os.path.join(args.output_dir, "screening_job_plan.csv"), job_plan_rows(jobs))

    if args.dry_run:
        print(f"Dry run completed. Planned jobs: {len(jobs)}")
        print(f"Output directory: {os.path.abspath(args.output_dir)}")
        return 0

    if not os.path.exists(args.template):
        raise FileNotFoundError(f"GPP21 template not found: {args.template}")

    print("Initializing PyRosetta...")
    pyrosetta_version = initialize_pyrosetta(args)
    template_pose = pose_from_pdb(args.template)
    check_template_pose(template_pose)
    clear_pose_constraints(template_pose)

    scorefxn_raw = create_score_function(args.scorefxn)
    scorefxn_cst = scorefxn_raw.clone()
    scorefxn_cst.set_weight(
        SCORE_TYPES["atom_pair_constraint"], args.atom_pair_cst_weight
    )
    scorefxn_cst.set_weight(
        SCORE_TYPES["coordinate_constraint"], args.coordinate_cst_weight
    )

    replicate_rows, component_rows, best_poses = run_screening(
        jobs,
        template_pose,
        scorefxn_raw,
        scorefxn_cst,
        args,
    )
    summaries, construct_windows = summarize_replicates(replicate_rows, args)
    single_summaries = [row for row in summaries if row["Job_Type"] == "single"]
    construct_summary_rows = construct_level_summary(
        construct_windows,
        single_summaries,
    )

    write_csv(
        os.path.join(args.output_dir, "screening_replicates.csv"), replicate_rows
    )
    write_csv(os.path.join(args.output_dir, "screening_summary.csv"), summaries)
    write_csv(
        os.path.join(args.output_dir, "single_candidate_summary.csv"),
        single_summaries,
    )
    write_csv(
        os.path.join(args.output_dir, "construct_window_summary.csv"),
        construct_windows,
    )
    write_csv(
        os.path.join(args.output_dir, "construct_summary.csv"),
        construct_summary_rows,
    )
    write_csv(
        os.path.join(args.output_dir, "nonadditivity_components.csv"),
        component_rows,
    )
    if args.dump_pdbs:
        dump_best_pdbs(args.output_dir, best_poses)
    write_manifest(args.output_dir, args, pyrosetta_version, scorefxn_raw)

    print("\nScreening completed.")
    print(f"Jobs: {len(jobs)}")
    print(f"Replicate rows: {len(replicate_rows)}")
    print(f"Output directory: {os.path.abspath(args.output_dir)}")
    return 0


if __name__ == "__main__":
    main(NOTEBOOK_CONFIG)
