"""Chemical transformation and geometry engine for 3D Functional Group Toolbox."""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

logger = logging.getLogger(__name__)

# Fixed seed so the same replacement always produces the same geometry.
_EMBED_SEED = 0xF00D
# Number of trial rotations about the new bond when placing a group.
_TORSION_STEPS = 36


def _skew(v: np.ndarray) -> np.ndarray:
    """Return the cross-product matrix of a 3-vector."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def _rotation_matrix_from_vectors(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """Compute the 3x3 rotation matrix that rotates the direction of vec1 onto vec2."""
    norm1 = float(np.linalg.norm(vec1))
    norm2 = float(np.linalg.norm(vec2))
    if norm1 < 1e-6 or norm2 < 1e-6:
        return np.eye(3)

    a = vec1 / norm1
    b = vec2 / norm2
    dot = float(np.dot(a, b))

    if dot > 0.999999:
        return np.eye(3)
    if dot < -0.999999:
        # 180 degree turn about any axis perpendicular to a
        ortho = (
            np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        )
        axis = np.cross(a, ortho)
        k = _skew(axis / np.linalg.norm(axis))
        return np.eye(3) + 2 * k @ k

    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    k = _skew(v)
    return np.eye(3) + k + k @ k * ((1.0 - dot) / (s**2))


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotation matrix for `angle` radians about the unit vector `axis` (Rodrigues)."""
    k = _skew(axis)
    return np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * (k @ k)


def _positions(conf: Chem.Conformer) -> np.ndarray:
    return np.asarray(conf.GetPositions(), dtype=float)


def _parse_group(group_smiles: str) -> tuple[Chem.Mol, int, int]:
    """Parse a group SMILES and return (fragment, dummy index, attach-atom index)."""
    fragment = Chem.MolFromSmiles(group_smiles)
    if fragment is None:
        raise ValueError("Invalid functional-group SMILES")

    dummies = [a.GetIdx() for a in fragment.GetAtoms() if a.GetAtomicNum() == 0]
    if not dummies:
        raise ValueError("Functional group has no attachment point")
    if len(dummies) > 1:
        raise ValueError("Functional group must have exactly one attachment point")

    dummy = dummies[0]
    neighbours = fragment.GetAtomWithIdx(dummy).GetNeighbors()
    if len(neighbours) != 1:
        raise ValueError("Functional group attachment point must have one bond")
    return fragment, dummy, neighbours[0].GetIdx()


def _embed_fragment(
    fragment: Chem.Mol, dummy: int, root_atomic_num: int
) -> Chem.Mol | None:
    """Embed the group in 3D, with the dummy standing in for the root atom.

    Atom indices of the heavy atoms are unchanged; hydrogens are appended.
    """
    frag = Chem.RWMol(fragment)
    frag.GetAtomWithIdx(dummy).SetAtomicNum(root_atomic_num or 6)
    frag.GetAtomWithIdx(dummy).SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(frag)
        frag_h = Chem.AddHs(frag.GetMol())
        params = AllChem.ETKDGv3()
        params.randomSeed = _EMBED_SEED
        if AllChem.EmbedMolecule(frag_h, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(frag_h, params) != 0:
                return None
    except (RuntimeError, ValueError) as exc:
        logger.debug("Fragment embedding failed: %s", exc)
        return None
    return frag_h


def _best_torsion(
    rel: np.ndarray,
    origin: np.ndarray,
    axis: np.ndarray,
    others: np.ndarray,
) -> np.ndarray:
    """Rotate `rel` (offsets from origin) about `axis` to maximise clearance from `others`."""
    if len(others) == 0 or len(rel) == 0:
        return rel

    best_rel, best_score = rel, -1.0
    for angle in np.linspace(0.0, 2.0 * np.pi, _TORSION_STEPS, endpoint=False):
        trial = rel @ _axis_rotation(axis, angle).T
        diff = (origin + trial)[:, None, :] - others[None, :, :]
        score = float(np.min(np.einsum("ijk,ijk->ij", diff, diff)))
        if score > best_score:
            best_rel, best_score = trial, score
    return best_rel


def _linear_layout(
    origin: np.ndarray, direction: np.ndarray, count: int
) -> list[np.ndarray]:
    """Fallback zig-zag placement along the bond direction when embedding fails."""
    u = np.array([0.0, 1.0, 0.0])
    if abs(np.dot(direction, u)) > 0.9:
        u = np.array([0.0, 0.0, 1.0])
    u = u - np.dot(u, direction) * direction
    u /= np.linalg.norm(u)
    w = np.cross(direction, u)
    return [
        origin
        + direction * (1.35 * offset)
        + u * (0.35 * (offset % 2))
        + w * (0.2 * ((offset // 2) % 2))
        for offset in range(1, count + 1)
    ]


def replace_atom_with_group(
    mol: Chem.Mol,
    target: int,
    group_smiles: str,
    relax: bool = True,
    max_iters: int = 500,
) -> Chem.Mol:
    """Return a molecule with one atom replaced by a one-attachment functional group.

    The target atom becomes the group's attachment atom and keeps its bonds to
    heavy (or, for a hydrogen target, any) neighbours. Terminal hydrogens on a
    heavy target are dropped, since the group brings its own. The input
    molecule is not modified.

    Args:
        mol: The input RDKit molecule.
        target: Atom index in `mol` to be replaced.
        group_smiles: SMILES string containing a dummy atom [*:1] attachment point.
        relax: If True and 3D conformers exist, performs constrained MMFF/UFF relaxation
               of the newly added atoms while holding parent atoms rigid.
        max_iters: Maximum iterations for force-field minimization.

    Returns:
        New Chem.Mol with the functional group attached and valence-completing hydrogens.
    """
    fragment, dummy, attach = _parse_group(group_smiles)
    if not 0 <= target < mol.GetNumAtoms():
        raise ValueError(f"Atom index {target} is out of range")

    rw = Chem.RWMol(mol)

    # Hydrogens on a heavy target would exceed the new atom's valence; the
    # group's own hydrogens are regenerated by AddHs below.
    if rw.GetAtomWithIdx(target).GetAtomicNum() != 1:
        dropped = sorted(
            (
                n.GetIdx()
                for n in rw.GetAtomWithIdx(target).GetNeighbors()
                if n.GetAtomicNum() == 1 and n.GetDegree() == 1
            ),
            reverse=True,
        )
        for idx in dropped:
            rw.RemoveAtom(idx)
        target -= sum(1 for idx in dropped if idx < target)

    base_num_atoms = rw.GetNumAtoms()
    has_conf = rw.GetNumConformers() > 0
    anchors = [n.GetIdx() for n in rw.GetAtomWithIdx(target).GetNeighbors()]
    root_atomic_num = rw.GetAtomWithIdx(anchors[0]).GetAtomicNum() if anchors else 6

    # Turn the target into the attachment atom, keeping its properties
    # (the host tracks atoms through them) and its existing bonds.
    source = fragment.GetAtomWithIdx(attach)
    target_atom = rw.GetAtomWithIdx(target)
    target_atom.SetAtomicNum(source.GetAtomicNum())
    target_atom.SetFormalCharge(source.GetFormalCharge())
    target_atom.SetIsAromatic(source.GetIsAromatic())
    target_atom.SetIsotope(source.GetIsotope())
    target_atom.SetNumRadicalElectrons(source.GetNumRadicalElectrons())
    target_atom.SetNumExplicitHs(source.GetNumExplicitHs())
    target_atom.SetNoImplicit(source.GetNoImplicit())
    target_atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    mapping = {attach: target}
    for atom in fragment.GetAtoms():
        if atom.GetIdx() not in (dummy, attach):
            mapping[atom.GetIdx()] = rw.AddAtom(Chem.Atom(atom))
    for bond in fragment.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if dummy not in (a, b):
            rw.AddBond(mapping[a], mapping[b], bond.GetBondType())

    Chem.SanitizeMol(rw)
    result = rw.GetMol()

    if has_conf:
        _place_group(
            result,
            target,
            anchors,
            fragment,
            dummy,
            attach,
            mapping,
            root_atomic_num,
            base_num_atoms,
        )

    # Make valence-completing hydrogens explicit for reliable 3D display.
    result = Chem.AddHs(result, addCoords=has_conf)

    if relax and has_conf:
        fixed_atoms = [i for i in range(base_num_atoms) if i != target]
        relax_molecule_with_fixed_atoms(result, fixed_atoms, max_iters=max_iters)

    return result


def _place_group(
    result: Chem.Mol,
    target: int,
    anchors: Sequence[int],
    fragment: Chem.Mol,
    dummy: int,
    attach: int,
    mapping: dict[int, int],
    root_atomic_num: int,
    base_num_atoms: int,
) -> None:
    """Write 3D coordinates for the attachment atom and the new group atoms."""
    conf = result.GetConformer()
    coords = _positions(conf)
    target_pos = coords[target]

    # Point the group away from the atoms the target stays bonded to.
    direction = np.array([1.0, 0.0, 0.0])
    if anchors:
        away = sum(
            (target_pos - coords[a]) / max(np.linalg.norm(target_pos - coords[a]), 1e-6)
            for a in anchors
        )
        if np.linalg.norm(away) < 1e-3:
            away = target_pos - coords[anchors[0]]
        if np.linalg.norm(away) > 1e-3:
            direction = away / np.linalg.norm(away)

    new_atoms = {f: mapping[f] for f in mapping if f != attach}
    others = coords[[i for i in range(base_num_atoms) if i != target]]
    frag_h = _embed_fragment(fragment, dummy, root_atomic_num)

    if frag_h is None:
        placed = _linear_layout(target_pos, direction, len(new_atoms))
        for idx, pos in zip(new_atoms.values(), placed):
            conf.SetAtomPosition(idx, Point3D(*map(float, pos)))
        return

    frag_xyz = _positions(frag_h.GetConformer())
    bond_vec = frag_xyz[attach] - frag_xyz[dummy]
    r_mat = _rotation_matrix_from_vectors(bond_vec, direction)

    # With a single anchor, use the group's own bond length (C-H -> C-C etc.).
    if len(anchors) == 1:
        target_pos = coords[anchors[0]] + direction * float(np.linalg.norm(bond_vec))
        conf.SetAtomPosition(target, Point3D(*map(float, target_pos)))

    # Score clearance with every group atom (hydrogens included) except those on
    # the stand-in root atom, but only write the heavy atoms: AddHs places the rest.
    scored = [
        a.GetIdx()
        for a in frag_h.GetAtoms()
        if a.GetIdx() not in (dummy, attach)
        and not (a.GetAtomicNum() == 1 and a.GetNeighbors()[0].GetIdx() == dummy)
    ]
    rel = (frag_xyz[scored] - frag_xyz[attach]) @ r_mat.T
    rel = _best_torsion(rel, target_pos, direction, others)
    offsets = dict(zip(scored, rel))
    for f_idx, idx in new_atoms.items():
        conf.SetAtomPosition(idx, Point3D(*map(float, target_pos + offsets[f_idx])))


def relax_molecule_with_fixed_atoms(
    mol: Chem.Mol,
    fixed_atom_indices: Sequence[int],
    max_iters: int = 500,
) -> bool:
    """Relax unconstrained atoms using MMFF94 or UFF with specified fixed atom points.

    Returns True if a force field could be set up and minimisation ran, False otherwise.
    """
    if not mol.GetNumConformers():
        return False

    def minimise(ff: object) -> bool:
        if ff is None:
            return False
        for idx in fixed_atom_indices:
            if 0 <= idx < mol.GetNumAtoms():
                ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=max_iters)
        return True

    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None and minimise(
            AllChem.MMFFGetMoleculeForceField(mol, props)
        ):
            return True
    except (RuntimeError, ValueError, AttributeError) as exc:
        logger.debug("MMFF relaxation not applicable: %s", exc)

    try:
        if minimise(AllChem.UFFGetMoleculeForceField(mol)):
            return True
    except (RuntimeError, ValueError, AttributeError) as exc:
        logger.debug("UFF relaxation failed: %s", exc)

    return False
