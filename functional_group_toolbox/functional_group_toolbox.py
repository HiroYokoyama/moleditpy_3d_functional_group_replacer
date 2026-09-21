"""Interactive one-attachment functional-group replacement tool."""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QComboBox, QPushButton, QLabel, QMessageBox
from PyQt6.QtCore import Qt, QEvent, QObject
from rdkit import Chem

PLUGIN_NAME = "Functional Group Toolbox"
PLUGIN_VERSION = "2026.09.21"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Replace a selected atom with a common functional group while retaining its neighbours."
PLUGIN_DEPENDENCIES = ["rdkit", "PyQt6"]
GROUPS = {
    "Methyl": "[*:1]C", "Ethyl": "[*:1]CC", "Hydroxyl": "[*:1]O",
    "Amino": "[*:1]N", "Methoxy": "[*:1]OC", "Acetyl": "[*:1]C(=O)C",
    "Carboxyl": "[*:1]C(=O)O", "Cyano": "[*:1]C#N", "Nitro": "[*:1][N+](=O)[O-]",
    "Phenyl": "[*:1]c1ccccc1", "tert-Butyl": "[*:1]C(C)(C)C",
    "Fluoro": "[*:1]F", "Chloro": "[*:1]Cl", "Bromo": "[*:1]Br", "Iodo": "[*:1]I",
}

def replace_atom_with_group(mol, target, group_smiles):
    """Return a molecule with one atom replaced by a one-point group."""
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
    rw = Chem.RWMol(mol)
    target_atom = rw.GetAtomWithIdx(target)
    source_atom = fragment.GetAtomWithIdx(attach)
    target_atom.SetAtomicNum(source_atom.GetAtomicNum())
    target_atom.SetFormalCharge(source_atom.GetFormalCharge())
    target_atom.SetIsAromatic(source_atom.GetIsAromatic())
    mapping = {attach: target}
    for atom in fragment.GetAtoms():
        if atom.GetIdx() != dummy and atom.GetIdx() != attach:
            mapping[atom.GetIdx()] = rw.AddAtom(Chem.Atom(atom))
    for bond in fragment.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if dummy not in (a, b) and rw.GetBondBetweenAtoms(mapping[a], mapping[b]) is None:
            rw.AddBond(mapping[a], mapping[b], bond.GetBondType())
    Chem.SanitizeMol(rw)
    return rw.GetMol()

class _AtomPickFilter(QObject):
    """Observe viewport clicks without consuming host camera gestures."""
    def __init__(self, callback, parent=None):
        super().__init__(parent)
        self.callback = callback
        self.press_pos = None

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
            self.press_pos = event.position().toPoint()
        elif event.type() == QEvent.Type.MouseButtonRelease and self.press_pos is not None:
            pos = event.position().toPoint()
            dx = pos.x() - self.press_pos.x()
            dy = pos.y() - self.press_pos.y()
            if dx * dx + dy * dy <= 25:
                self.callback(pos.x(), pos.y(), obj)
            self.press_pos = None
        return False

class FunctionalGroupToolbox(QWidget):
    def __init__(self, context):
        super().__init__(context.get_main_window())
        self.context = context
        self._pick_filter = None
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(300, 150)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Atom index to replace:"))
        self.atom_combo = QComboBox()
        layout.addWidget(self.atom_combo)
        self.pick_button = QPushButton("Pick atom in 3D")
        self.pick_button.setCheckable(True)
        self.pick_button.setToolTip("Click an atom in the 3D view to select it")
        self.pick_button.toggled.connect(self._toggle_3d_picking)
        layout.addWidget(self.pick_button)
        layout.addWidget(QLabel("Functional group:"))
        self.group_combo = QComboBox()
        self.group_combo.addItems(GROUPS)
        layout.addWidget(self.group_combo)
        button = QPushButton("Replace selected atom")
        button.clicked.connect(self.replace_atom)
        layout.addWidget(button)
        self.refresh_atoms()
        self._install_3d_picking()
        context.register_window("functional_group_toolbox", self)

    def _install_3d_picking(self):
        plotter = getattr(self.context, "plotter", None)
        interactor = getattr(plotter, "interactor", None) if plotter else None
        if interactor is not None:
            self._pick_filter = _AtomPickFilter(self._pick_atom, self)
            interactor.installEventFilter(self._pick_filter)

    def _toggle_3d_picking(self, enabled):
        self.context.show_status_message("3D atom picking enabled." if enabled else "3D atom picking disabled.")

    def _pick_atom(self, x, y, widget):
        if not self.pick_button.isChecked():
            return
        import vtk
        plotter = getattr(self.context, "plotter", None)
        mol = self.context.current_mol
        if plotter is None or mol is None or not mol.GetNumConformers():
            return
        picker = vtk.vtkCellPicker()
        ratio = widget.devicePixelRatioF()
        picker.SetTolerance(0.005)
        picker.Pick(x * ratio, (widget.height() - y) * ratio, 0, plotter.renderer)
        pos = picker.GetPickPosition()
        atom = min(mol.GetAtoms(), key=lambda candidate: self._distance_sq(mol, candidate.GetIdx(), pos))
        self.atom_combo.setCurrentIndex(atom.GetIdx())
        self.context.show_status_message(f"Selected atom {atom.GetIdx()} ({atom.GetSymbol()}).")

    def _distance_sq(self, mol, index, pos):
        point = mol.GetConformer().GetAtomPosition(index)
        return (point.x - pos[0]) ** 2 + (point.y - pos[1]) ** 2 + (point.z - pos[2]) ** 2

    def closeEvent(self, event):
        plotter = getattr(self.context, "plotter", None)
        interactor = getattr(plotter, "interactor", None) if plotter else None
        if interactor is not None and self._pick_filter is not None:
            interactor.removeEventFilter(self._pick_filter)
        super().closeEvent(event)

    def refresh_atoms(self):
        self.atom_combo.clear()
        mol = self.context.current_mol
        if mol:
            self.atom_combo.addItems([f"{a.GetIdx()}: {a.GetSymbol()}" for a in mol.GetAtoms()])

    def replace_atom(self):
        mol = self.context.current_mol
        target = self.atom_combo.currentIndex()
        if mol is None or target < 0:
            return
        try:
            new_mol = replace_atom_with_group(mol, target, GROUPS[self.group_combo.currentText()])
            self.context.current_molecule = new_mol
            self.context.push_undo_checkpoint()
            refresh = getattr(self.context, "refresh_3d_view", None)
            if callable(refresh):
                refresh()
            self.refresh_atoms()
            self.context.show_status_message("Functional group replacement applied.")
        except (RuntimeError, ValueError, AttributeError) as exc:
            QMessageBox.critical(self, "Replacement failed", str(exc))


def initialize(context):
    def show():
        window = context.get_window("functional_group_toolbox")
        if window is None:
            window = FunctionalGroupToolbox(context)
        window.show(); window.raise_(); window.activateWindow()
    context.add_menu_action("Edit/Functional Group Toolbox...", show)
