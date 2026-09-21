"""Chemical transformation and geometry engine for 3D Functional Group Toolbox."""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

logger = logging.getLogger(__name__)


def _rotation_matrix_from_vectors(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """Compute the 3x3 rotation matrix that rotates unit vector vec1 into vec2."""
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
        ortho = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, ortho)
        axis /= np.linalg.norm(axis)
        k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        return np.eye(3) + 2 * k @ k

    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1.0 - dot) / (s ** 2))


def replace_atom_with_group(
    mol: Chem.Mol,
    target: int,
    group_smiles: str,
    relax: bool = True,
    max_iters: int = 500,
) -> Chem.Mol:
    """Return a molecule with one atom replaced by a one-attachment functional group.

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
    fragment = Chem.MolFromSmiles(group_smiles)
    if fragment is None:
        raise ValueError("Invalid functional-group SMILES")

    dummy = next((a.GetIdx() for a in fragment.GetAtoms() if a.GetAtomicNum() == 0), None)
    if dummy is None:
        raise ValueError("Functional group has no attachment point")

    neighbours = list(fragment.GetAtomWithIdx(dummy).GetNeighbors())
    if len(neighbours) != 1:
        raise ValueError("Functional group must have one attachment point")
    attach = neighbours[0].GetIdx()

    orig_num_atoms = mol.GetNumAtoms()
    has_conf = mol.GetNumConformers() > 0
    target_pos = np.array([0.0, 0.0, 0.0])
    dir_vec = np.array([1.0, 0.0, 0.0])
    root_atomic_num = 6

    if has_conf:
        conf = mol.GetConformer()
        t_pt = conf.GetAtomPosition(target)
        target_pos = np.array([t_pt.x, t_pt.y, t_pt.z])
        target_atom_mol = mol.GetAtomWithIdx(target)
        orig_neighbors = target_atom_mol.GetNeighbors()
        if orig_neighbors:
            root_atom = orig_neighbors[0]
            root_atomic_num = root_atom.GetAtomicNum()
            r_pt = conf.GetAtomPosition(root_atom.GetIdx())
            root_pos = np.array([r_pt.x, r_pt.y, r_pt.z])
            v = target_pos - root_pos
            norm_v = np.linalg.norm(v)
            if norm_v > 1e-3:
                dir_vec = v / norm_v

    rw = Chem.RWMol(mol)
    target_atom = rw.GetAtomWithIdx(target)
    source_atom = fragment.GetAtomWithIdx(attach)
    target_atom.SetAtomicNum(source_atom.GetAtomicNum())
    target_atom.SetFormalCharge(source_atom.GetFormalCharge())
    target_atom.SetIsAromatic(source_atom.GetIsAromatic())
    target_atom.SetNoImplicit(False)

    mapping = {attach: target}
    for atom in fragment.GetAtoms():
        if atom.GetIdx() != dummy and atom.GetIdx() != attach:
            mapping[atom.GetIdx()] = rw.AddAtom(Chem.Atom(atom))
    for bond in fragment.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if dummy not in (a, b) and rw.GetBondBetweenAtoms(mapping[a], mapping[b]) is None:
            rw.AddBond(mapping[a], mapping[b], bond.GetBondType())

    Chem.SanitizeMol(rw)
    result = rw.GetMol()

    if result.GetNumConformers():
        conf = result.GetConformer()

        # Embed a 3D copy of the fragment to preserve chemically accurate bond angles and lengths
        frag_copy = Chem.Mol(fragment)
        frag_dummy = next(a.GetIdx() for a in frag_copy.GetAtoms() if a.GetAtomicNum() == 0)
        frag_copy.GetAtomWithIdx(frag_dummy).SetAtomicNum(root_atomic_num if root_atomic_num > 0 else 6)

        embed_success = False
        try:
            embed_res = AllChem.EmbedMolecule(frag_copy, AllChem.ETKDGv3())
            if embed_res != 0:
                embed_res = AllChem.EmbedMolecule(frag_copy)
            embed_success = (embed_res == 0 and frag_copy.GetNumConformers() > 0)
        except (RuntimeError, ValueError, AttributeError):
            embed_success = False

        if embed_success:
            frag_conf = frag_copy.GetConformer()
            d_pos = np.array(frag_conf.GetAtomPosition(frag_dummy))
            a_pos = np.array(frag_conf.GetAtomPosition(attach))
            v_frag = a_pos - d_pos
            r_mat = _rotation_matrix_from_vectors(v_frag, dir_vec)

            for f_atom in fragment.GetAtoms():
                f_idx = f_atom.GetIdx()
                if f_idx in (dummy, attach):
                    continue
                p_f = np.array(frag_conf.GetAtomPosition(f_idx))
                rel_f = p_f - a_pos
                new_pos = target_pos + r_mat @ rel_f
                conf.SetAtomPosition(mapping[f_idx], Point3D(float(new_pos[0]), float(new_pos[1]), float(new_pos[2])))
        else:
            # Fallback linear layout if fragment 3D embedding fails
            u = np.array([0.0, 1.0, 0.0])
            if abs(np.dot(dir_vec, u)) > 0.9:
                u = np.array([0.0, 0.0, 1.0])
            u = u - np.dot(u, dir_vec) * dir_vec
            u_norm = np.linalg.norm(u)
            if u_norm > 1e-4:
                u = u / u_norm
            else:
                u = np.array([0.0, 1.0, 0.0])
            w = np.cross(dir_vec, u)

            added = [idx for idx in mapping.values() if idx != target]
            for offset, idx in enumerate(added, start=1):
                pos = (
                    target_pos
                    + dir_vec * (1.35 * offset)
                    + u * (0.35 * (offset % 2))
                    + w * (0.2 * ((offset // 2) % 2))
                )
                conf.SetAtomPosition(idx, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))

    # Make valence-completing hydrogens explicit for reliable 3D display.
    result = Chem.AddHs(result, addCoords=has_conf)

    if relax and result.GetNumConformers():
        fixed_atoms = [i for i in range(orig_num_atoms) if i != target]
        relax_molecule_with_fixed_atoms(result, fixed_atoms, max_iters=max_iters)

    return result


def relax_molecule_with_fixed_atoms(
    mol: Chem.Mol,
    fixed_atom_indices: Sequence[int],
    max_iters: int = 500,
) -> bool:
    """Relax unconstrained atoms using MMFF94 or UFF with specified fixed atom points.

    Returns True if relaxation converged or completed without error, False otherwise.
    """
    if not mol.GetNumConformers():
        return False

    # Attempt MMFF94 first
    try:
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        if mp is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
            if ff is not None:
                for idx in fixed_atom_indices:
                    if 0 <= idx < mol.GetNumAtoms():
                        ff.AddFixedPoint(idx)
                ff.Minimize(maxIts=max_iters)
                return True
    except (RuntimeError, ValueError, AttributeError) as exc:
        logger.debug("MMFF relaxation not applicable: %s", exc)

    # Fallback to UFF
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
        if ff is not None:
            for idx in fixed_atom_indices:
                if 0 <= idx < mol.GetNumAtoms():
                    ff.AddFixedPoint(idx)
            ff.Minimize(maxIts=max_iters)
            return True
    except (RuntimeError, ValueError, AttributeError) as exc:
        logger.debug("UFF relaxation failed: %s", exc)

    return False
