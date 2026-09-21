"""Comprehensive unit and GUI tests for 3D Functional Group Toolbox."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")
from PyQt6.QtCore import QEvent, QPointF, Qt
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import QApplication

# Ensure repo root is in sys.path
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import functional_group_toolbox_3d as module
from functional_group_toolbox_3d.chemistry import (
    relax_molecule_with_fixed_atoms,
    replace_atom_with_group,
)
from functional_group_toolbox_3d.dialog import FunctionalGroupToolbox, _AtomPickFilter
from functional_group_toolbox_3d.groups import (
    GROUPS,
    get_group_smiles,
    get_groups_by_category,
    search_groups,
)


@pytest.fixture(scope="session")
def qapp():
    """Ensure a QApplication instance exists for GUI tests."""
    app = QApplication.instance()
    if app is None:
        app = QApplication(["--platform", "offscreen"])
    return app


def test_package_metadata_and_public_exports():
    """Verify package public metadata and contract."""
    assert module.PLUGIN_NAME == "3D Functional Group Toolbox"
    assert hasattr(module, "PLUGIN_VERSION")
    assert module.PLUGIN_AUTHOR == "HiroYokoyama"
    assert module.PLUGIN_CATEGORY == "3D Editing"
    assert "rdkit" in module.PLUGIN_DEPENDENCIES
    assert "PyQt6" in module.PLUGIN_DEPENDENCIES


def test_backward_compatible_module_import():
    """Verify functional_group_toolbox_3d.py re-exports match package exports."""
    import functional_group_toolbox_3d.functional_group_toolbox_3d as compat_mod
    assert compat_mod.PLUGIN_NAME == module.PLUGIN_NAME
    assert compat_mod.GROUPS == GROUPS
    assert compat_mod.FunctionalGroupToolbox is FunctionalGroupToolbox
    assert compat_mod.replace_atom_with_group is replace_atom_with_group


def test_group_library_and_categorization():
    """Verify functional group library completeness and categorization."""
    # Must have at least 50 groups
    assert len(GROUPS) >= 50
    assert "Methyl" in GROUPS
    assert "Phenyl" in GROUPS
    assert "Trifluoromethyl" in GROUPS
    assert "Carboxyl" in GROUPS
    assert "Azido" in GROUPS

    # All groups must have a single attachment dummy atom [*:1]
    for name, smiles in GROUPS.items():
        assert "[*:1]" in smiles, f"Group {name} missing [*:1]"
        assert get_group_smiles(name) == smiles

    # Verify categories cover all groups
    all_cat_groups = set(get_groups_by_category("All"))
    assert all_cat_groups == set(GROUPS.keys())
    assert len(get_groups_by_category("Alkyl & Aliphatic")) > 0
    assert len(get_groups_by_category("Aryl & Heteroaryl")) > 0
    assert len(get_groups_by_category("Oxygen & Carbonyl")) > 0
    assert len(get_groups_by_category("Nitrogen & Amine")) > 0
    assert len(get_groups_by_category("Halogen")) > 0
    assert len(get_groups_by_category("Sulfur & Phosphorus")) > 0


def test_search_groups():
    """Test group search filter utility."""
    assert "Phenyl" in search_groups("phenyl")
    assert "Trifluoromethyl" in search_groups("fluoro")
    assert "Dimethylamino" in search_groups("amino")
    assert len(search_groups("")) == len(GROUPS)


def test_replacement_preserves_neighbors_and_adds_group_atoms():
    """Verify chemical replacement on a 2D Mol."""
    mol = Chem.MolFromSmiles("CC")
    result = replace_atom_with_group(mol, 1, "[*:1]CC", relax=False)
    heavy_atoms = [a for a in result.GetAtoms() if a.GetAtomicNum() != 1]
    assert len(heavy_atoms) == 3
    assert result.GetNumAtoms() > len(heavy_atoms)
    assert result.GetBondBetweenAtoms(0, 1) is not None
    assert result.GetBondBetweenAtoms(1, 2) is not None
    assert any(a.GetAtomicNum() == 1 for a in result.GetAtoms())


def test_replacement_rejects_missing_or_multiple_attachments():
    """Verify invalid group SMILES handling."""
    mol = Chem.MolFromSmiles("CC")
    with pytest.raises(ValueError, match="attachment point"):
        replace_atom_with_group(mol, 1, "CC")

    with pytest.raises(ValueError, match="Invalid functional-group"):
        replace_atom_with_group(mol, 1, "invalid_smiles_string")


def test_all_groups_replace_hydrogen_in_3d():
    """Verify that every functional group in GROUPS can replace a hydrogen in a 3D molecule."""
    base_mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(base_mol)
    h_idx = next(a.GetIdx() for a in base_mol.GetAtoms() if a.GetAtomicNum() == 1)

    # Sample key groups across all categories to verify execution and conformer integrity
    sample_groups = [
        "Methyl", "Isopropyl", "tert-Butyl", "Cyclopropyl", "Trifluoromethyl",
        "Vinyl", "Ethynyl", "Phenyl", "4-Pyridyl", "2-Thienyl",
        "Hydroxyl", "Methoxy", "Acetyl", "Carboxyl", "Amino",
        "Cyano", "Nitro", "Azido", "Fluoro", "Thiol", "Methylsulfonyl",
    ]
    for gname in sample_groups:
        smi = GROUPS[gname]
        res = replace_atom_with_group(base_mol, h_idx, smi, relax=True, max_iters=100)
        assert res is not None
        assert res.GetNumConformers() == 1
        assert res.GetNumAtoms() >= base_mol.GetNumAtoms()


def test_3d_relaxation_holds_parent_atoms_rigid():
    """Verify that force-field relaxation keeps parent atoms fixed."""
    benzene = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(benzene)
    conf = benzene.GetConformer()
    orig_positions = {i: conf.GetAtomPosition(i) for i in range(benzene.GetNumAtoms())}

    target_h = 6
    res = replace_atom_with_group(benzene, target_h, GROUPS["Phenyl"], relax=True)
    res_conf = res.GetConformer()

    # Verify all non-target parent atoms stayed in their exact positions
    for i in range(benzene.GetNumAtoms()):
        if i != target_h:
            orig_p = orig_positions[i]
            new_p = res_conf.GetAtomPosition(i)
            dist_sq = (new_p.x - orig_p.x)**2 + (new_p.y - orig_p.y)**2 + (new_p.z - orig_p.z)**2
            assert dist_sq < 1e-4, f"Parent atom {i} moved unexpectedly!"


def test_nitrile_and_second_atom_angle_orientation():
    """Verify that nitrile C-C#N attachment maintains linear 180 degree geometry along the bond vector."""
    benzene = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(benzene)
    conf = benzene.GetConformer()
    target_h = 6
    root_c = benzene.GetAtomWithIdx(target_h).GetNeighbors()[0].GetIdx()
    r_pos = np.array(conf.GetAtomPosition(root_c))

    res = replace_atom_with_group(benzene, target_h, GROUPS["Cyano"], relax=False)
    res_conf = res.GetConformer()

    # Find the nitrogen atom in the attached nitrile group
    nitrile_c = target_h
    nitrile_c_pos = np.array(res_conf.GetAtomPosition(nitrile_c))
    nitrile_n = next(a.GetIdx() for a in res.GetAtomWithIdx(nitrile_c).GetNeighbors() if a.GetAtomicNum() == 7)
    nitrile_n_pos = np.array(res_conf.GetAtomPosition(nitrile_n))

    # Vector 1: root_c -> nitrile_c
    v1 = nitrile_c_pos - r_pos
    # Vector 2: nitrile_c -> nitrile_n
    v2 = nitrile_n_pos - nitrile_c_pos

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_deg = np.arccos(np.clip(cos_angle, -1.0, 1.0)) * 180.0 / np.pi
    # Deviation from collinear (0 degrees between the two sequential bond vectors, i.e. 180 deg bond angle)
    assert angle_deg < 2.0, f"Nitrile is bent! Angle deviation: {angle_deg:.2f} deg"


def test_relax_molecule_fallbacks_and_zero_conformers():
    """Test relaxation when molecule has no conformers or falls back to UFF."""
    mol_2d = Chem.MolFromSmiles("CC")
    assert not relax_molecule_with_fixed_atoms(mol_2d, [0])

    mol_3d = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(mol_3d)

    # Force MMFF to fail by mocking MMFFGetMoleculeProperties to return None
    with patch("rdkit.Chem.AllChem.MMFFGetMoleculeProperties", return_value=None):
        assert relax_molecule_with_fixed_atoms(mol_3d, [0])


def test_pick_filter_click_vs_drag(qapp):
    """Test _AtomPickFilter distinguishes clicks from drags."""
    callback = MagicMock()
    pick_filter = _AtomPickFilter(callback)

    widget = MagicMock()

    # Small movement (click) -> triggers callback
    press_ev = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    release_ev = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(102.0, 101.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    pick_filter.eventFilter(widget, press_ev)
    pick_filter.eventFilter(widget, release_ev)
    assert callback.call_count == 1

    # Large movement (camera rotation/drag) -> does NOT trigger callback
    callback.reset_mock()
    press_drag = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    release_drag = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(160.0, 150.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    pick_filter.eventFilter(widget, press_drag)
    pick_filter.eventFilter(widget, release_drag)
    assert callback.call_count == 0

    assert not pick_filter.eventFilter(widget, None)


def test_dialog_ui_and_search_filter(qapp):
    """Test FunctionalGroupToolbox UI components, category changes, and search filtering."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(mol)

    mock_plotter = MagicMock()
    mock_plotter.interactor = MagicMock()
    mock_plotter.camera_position = [(0, 0, 10), (0, 0, 0), (0, 1, 0)]
    mock_actor = MagicMock()
    mock_plotter.add_point_labels.return_value = mock_actor

    mock_context = MagicMock()
    mock_context.current_mol = mol
    mock_context.plotter = mock_plotter
    mock_context.get_main_window.return_value = None

    dlg = FunctionalGroupToolbox(mock_context)

    # Initial state
    assert dlg.selected_atom_idx is None
    assert dlg.selection_label.text() == "No atom selected"
    assert not dlg.replace_button.isEnabled()
    assert dlg.group_combo.count() == len(GROUPS)

    # Category filter
    dlg.category_combo.setCurrentText("Halogen")
    assert dlg.group_combo.count() == 4
    assert {dlg.group_combo.itemText(i) for i in range(4)} == {"Fluoro", "Chloro", "Bromo", "Iodo"}

    # Search filter
    dlg.category_combo.setCurrentText("All")
    dlg.search_input.setText("phenyl")
    filtered_items = [dlg.group_combo.itemText(i) for i in range(dlg.group_combo.count())]
    assert "Phenyl" in filtered_items
    assert "4-Fluorophenyl" in filtered_items
    assert "Methyl" not in filtered_items

    # Clear search
    dlg.search_input.setText("")
    assert dlg.group_combo.count() == len(GROUPS)

    dlg.close()


def test_dialog_atom_selection_and_3d_labels(qapp):
    """Test selecting an atom, updating UI label, 3D point label creation and clearing."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol)

    mock_plotter = MagicMock()
    mock_plotter.interactor = MagicMock()
    mock_plotter.camera_position = [(0, 0, 10), (0, 0, 0), (0, 1, 0)]
    label_actor = MagicMock()
    mock_plotter.add_point_labels.return_value = label_actor

    mock_context = MagicMock()
    mock_context.current_mol = mol
    mock_context.plotter = mock_plotter
    mock_context.get_main_window.return_value = None

    dlg = FunctionalGroupToolbox(mock_context)

    # Select atom 0
    dlg.selected_atom_idx = 0
    dlg.update_selection_display()

    assert dlg.selection_label.text() == "Selected atom: C0 (index 0)"
    assert dlg.replace_button.isEnabled()
    assert len(dlg.selection_labels) == 1
    assert mock_plotter.add_point_labels.called

    # Clear selection
    dlg.clear_selection()
    assert dlg.selected_atom_idx is None
    assert dlg.selection_label.text() == "No atom selected"
    assert not dlg.replace_button.isEnabled()
    assert len(dlg.selection_labels) == 0
    assert mock_plotter.remove_actor.called

    dlg.close()


def test_dialog_pick_atom_logic(qapp):
    """Test _pick_atom with mocked VTK picker."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    atom0_pos = conf.GetAtomPosition(0)

    mock_plotter = MagicMock()
    mock_plotter.renderer = MagicMock()
    mock_plotter.interactor = MagicMock()

    mock_context = MagicMock()
    mock_context.current_mol = mol
    mock_context.plotter = mock_plotter

    mock_mw = MagicMock()
    mock_view_3d = MagicMock()
    atom_actor = MagicMock()
    mock_view_3d.atom_actor = atom_actor
    mock_mw.view_3d_manager = mock_view_3d
    mock_context.get_main_window.return_value = mock_mw

    dlg = FunctionalGroupToolbox(mock_context)

    # Mock vtk cell picker to return atom0 position
    mock_picker = MagicMock()
    mock_picker.GetActor.return_value = atom_actor
    mock_picker.GetPickPosition.return_value = (atom0_pos.x, atom0_pos.y, atom0_pos.z)
    mock_vtk = MagicMock()
    mock_vtk.vtkCellPicker.return_value = mock_picker

    with patch.dict(sys.modules, {"vtk": mock_vtk}):
        widget = MagicMock()
        widget.devicePixelRatioF.return_value = 1.0
        widget.height.return_value = 600

        # Pick atom 0
        dlg._pick_atom(100, 100, widget)
        assert dlg.selected_atom_idx == 0

        # Pick atom 0 again -> toggles to None
        dlg._pick_atom(100, 100, widget)
        assert dlg.selected_atom_idx is None

        # Pick non-matching actor -> does nothing
        mock_picker.GetActor.return_value = MagicMock()
        dlg._pick_atom(100, 100, widget)
        assert dlg.selected_atom_idx is None

        # Pick too far from any atom
        mock_picker.GetActor.return_value = atom_actor
        mock_picker.GetPickPosition.return_value = (999.0, 999.0, 999.0)
        dlg._pick_atom(100, 100, widget)
        assert dlg.selected_atom_idx is None

    dlg.close()


def test_dialog_replace_atom_action_and_fallbacks(qapp):
    """Test replace_atom execution, warning on no selection, and error handling."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol)

    mock_plotter = MagicMock()
    mock_plotter.interactor = MagicMock()
    mock_plotter.camera_position = [(0, 0, 10), (0, 0, 0), (0, 1, 0)]

    mock_context = MagicMock()
    mock_context.current_mol = mol
    mock_context.plotter = mock_plotter
    mock_context.refresh_3d_view = None

    mock_mw = MagicMock()
    mock_view_3d = MagicMock()
    mock_mw.view_3d_manager = mock_view_3d
    mock_context.get_main_window.return_value = mock_mw

    dlg = FunctionalGroupToolbox(mock_context)

    # Attempt replacement without selection
    with patch("PyQt6.QtWidgets.QMessageBox.warning") as mock_warn:
        dlg.replace_atom()
        assert mock_warn.called

    # Successful replacement using fallback draw_molecule_3d
    dlg.selected_atom_idx = 6  # hydrogen atom
    dlg.update_selection_display()
    dlg.group_combo.setCurrentText("Methyl")
    dlg.replace_atom()

    assert mock_view_3d.draw_molecule_3d.called
    assert mock_context.push_undo_checkpoint.called
    assert dlg.selected_atom_idx is None

    # Error handling branch
    dlg.selected_atom_idx = 0
    with (
        patch("functional_group_toolbox_3d.dialog.replace_atom_with_group", side_effect=RuntimeError("Test error")),
        patch("PyQt6.QtWidgets.QMessageBox.critical") as mock_crit,
    ):
        dlg.replace_atom()
        assert mock_crit.called

    dlg.close()


def test_initialize_and_lifecycle_handlers(qapp):
    """Test initialize(context), _open_toolbox, save, load, and reset handlers."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol)

    mock_plotter = MagicMock()
    mock_plotter.interactor = MagicMock()
    mock_plotter.camera_position = [(0, 0, 10), (0, 0, 0), (0, 1, 0)]

    mock_context = MagicMock()
    mock_context.current_mol = mol
    mock_context.plotter = mock_plotter
    mock_context.get_main_window.return_value = None
    mock_context.get_window.return_value = None

    handlers = {}
    mock_context.register_save_handler.side_effect = lambda h: handlers.setdefault("save", h)
    mock_context.register_load_handler.side_effect = lambda h: handlers.setdefault("load", h)
    mock_context.register_document_reset_handler.side_effect = lambda h: handlers.setdefault("reset", h)

    module.initialize(mock_context)
    mock_context.add_menu_action.assert_called_once_with(
        "3D Edit/3D Functional Group Toolbox...", module._open_toolbox
    )

    # Trigger opening dialog
    module._open_toolbox()
    dlg = mock_context.register_window.call_args[0][1]
    mock_context.get_window.return_value = dlg

    # Opening while existing is visible raises it
    dlg.setVisible(True)
    module._open_toolbox()

    # Save state
    saved = handlers["save"]()
    assert "settings" in saved

    # Load state
    handlers["load"]({"settings": {"last_category": "Halogen", "last_group": "Fluoro", "relax": False}})
    assert dlg.category_combo.currentText() == "Halogen"
    assert dlg.group_combo.currentText() == "Fluoro"
    assert not dlg.relax_checkbox.isChecked()

    # Reset state
    dlg.close()
    dlg.setVisible(False)
    handlers["reset"]()
    assert module._current_settings["last_category"] == "All"
    assert module._current_settings["last_group"] == "Methyl"


def test_dialog_position_near_parent(qapp):
    """Test dialog geometry placement relative to parent widget."""
    from PyQt6.QtWidgets import QWidget
    parent = QWidget()
    parent.resize(800, 600)
    parent.show()

    mock_plotter = MagicMock()
    mock_context = MagicMock()
    mock_context.current_mol = None
    mock_context.plotter = mock_plotter
    mock_context.get_main_window.return_value = parent

    dlg = FunctionalGroupToolbox(mock_context)
    dlg._position_near_parent()
    assert dlg.x() >= 0
    assert dlg.y() >= 0
    dlg.close()
    parent.close()


def test_standalone_module_execution():
    """Test executing functional_group_toolbox_3d.py without an active parent package."""
    import importlib.util
    target_path = REPO_ROOT / "functional_group_toolbox_3d" / "functional_group_toolbox_3d.py"
    spec = importlib.util.spec_from_file_location("__main__", target_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = ""
    spec.loader.exec_module(mod)
    assert hasattr(mod, "FunctionalGroupToolbox")


def test_lifecycle_edge_cases():
    """Test _open_toolbox with no context, save_state before dialog opened, and reset when visible."""
    # Temporarily clear _context
    orig_context = module._context
    orig_opened = module._dialog_opened
    try:
        module._context = None
        module._open_toolbox()  # Early return

        mock_context = MagicMock()
        mock_context.get_window.return_value = None
        module._dialog_opened = False
        handlers = {}
        mock_context.register_save_handler.side_effect = lambda h: handlers.setdefault("save", h)
        mock_context.register_document_reset_handler.side_effect = lambda h: handlers.setdefault("reset", h)
        module.initialize(mock_context)

        # Save before opened returns empty dict
        assert handlers["save"]() == {}

        # Reset when dialog is visible does nothing
        mock_dlg = MagicMock()
        mock_dlg.isVisible.return_value = True
        mock_context.get_window.return_value = mock_dlg
        module._dialog_opened = True
        handlers["reset"]()
        assert module._dialog_opened is True
    finally:
        module._context = orig_context
        module._dialog_opened = orig_opened

