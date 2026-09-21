from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

Chem = pytest.importorskip("rdkit.Chem")
SOURCE = Path(__file__).resolve().parents[1] / "functional_group_toolbox" / "functional_group_toolbox.py"
spec = importlib.util.spec_from_file_location("functional_group_toolbox_test", SOURCE)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_replacement_preserves_neighbors_and_adds_group_atoms():
    mol = Chem.MolFromSmiles("CC")
    result = module.replace_atom_with_group(mol, 1, "[*:1]CC")
    assert result.GetNumAtoms() == 3
    assert result.GetBondBetweenAtoms(0, 1) is not None
    assert result.GetBondBetweenAtoms(1, 2) is not None


def test_replacement_rejects_missing_attachment_point():
    mol = Chem.MolFromSmiles("CC")
    with pytest.raises(ValueError, match="attachment point"):
        module.replace_atom_with_group(mol, 1, "CC")


def test_group_library_contains_common_single_point_groups():
    assert {"Methyl", "Carboxyl", "Phenyl"}.issubset(module.GROUPS)
    assert all("[*:1]" in smiles for smiles in module.GROUPS.values())
