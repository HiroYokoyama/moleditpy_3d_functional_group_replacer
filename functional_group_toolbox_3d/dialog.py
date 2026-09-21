"""PyQt6 interactive dialog for 3D Functional Group Toolbox."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

from PyQt6.QtCore import QEvent, QObject, Qt
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .chemistry import replace_atom_with_group
from .groups import GROUP_CATEGORIES, GROUPS, get_group_smiles, search_groups

logger = logging.getLogger(__name__)


class _AtomPickFilter(QObject):
    """Observe 3D viewport clicks without consuming host camera rotation/pan gestures."""

    def __init__(self, callback: Any, parent: QObject | None = None):
        super().__init__(parent)
        self.callback = callback
        self.press_pos = None

    def eventFilter(self, obj: QObject | None, event: QEvent | None) -> bool:
        if event is None:
            return False

        if event.type() == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
            self.press_pos = event.position().toPoint()
        elif event.type() == QEvent.Type.MouseButtonRelease and self.press_pos is not None:
            pos = event.position().toPoint()
            dx = pos.x() - self.press_pos.x()
            dy = pos.y() - self.press_pos.y()
            # Threshold to distinguish a click from a camera drag/rotation
            if dx * dx + dy * dy <= 25:
                self.callback(pos.x(), pos.y(), obj)
            self.press_pos = None
        return False


class FunctionalGroupToolbox(QWidget):
    """Interactive functional group replacement tool with real-time 3D picking and labels."""

    def __init__(self, context: Any):
        parent = context.get_main_window() if hasattr(context, "get_main_window") else None
        if not isinstance(parent, QWidget):
            parent = None
        super().__init__(parent)
        self.context = context
        self._pick_filter: _AtomPickFilter | None = None
        self.selected_atom_idx: int | None = None
        self.selection_labels: list[Any] = []

        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setWindowTitle("3D Functional Group Toolbox")
        self.resize(380, 310)

        self._init_ui()
        self._install_3d_picking()
        self._position_near_parent()
        self.context.register_window("functional_group_toolbox", self)

    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # Instructions
        instruction = QLabel("Click an atom in the 3D view to select it for replacement.")
        instruction.setWordWrap(True)
        layout.addWidget(instruction)

        # Selection status display
        self.selection_label = QLabel("No atom selected")
        self.selection_label.setStyleSheet("font-weight: bold; color: #2a7ae2; padding: 2px 0;")
        layout.addWidget(self.selection_label)

        # Category filter
        cat_layout = QHBoxLayout()
        cat_layout.addWidget(QLabel("Category:"))
        self.category_combo = QComboBox()
        self.category_combo.addItems(list(GROUP_CATEGORIES.keys()))
        self.category_combo.currentTextChanged.connect(self._on_category_changed)
        cat_layout.addWidget(self.category_combo, stretch=1)
        layout.addLayout(cat_layout)

        # Search filter
        search_layout = QHBoxLayout()
        search_layout.addWidget(QLabel("Search:"))
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Filter groups by name or SMILES (e.g. phenyl, c1ccccc1, COOH, CF3)...")
        self.search_input.setClearButtonEnabled(True)
        self.search_input.textChanged.connect(self._on_search_changed)
        search_layout.addWidget(self.search_input, stretch=1)
        layout.addLayout(search_layout)

        # Functional group dropdown
        group_layout = QHBoxLayout()
        group_layout.addWidget(QLabel("Group:"))
        self.group_combo = QComboBox()
        self.group_combo.currentTextChanged.connect(self._update_group_preview)
        group_layout.addWidget(self.group_combo, stretch=1)
        layout.addLayout(group_layout)

        # SMILES preview
        self.preview_label = QLabel("")
        self.preview_label.setStyleSheet("color: #666; font-size: 11px;")
        layout.addWidget(self.preview_label)

        # 3D Relaxation option
        self.relax_checkbox = QCheckBox("Relax group 3D geometry (MMFF/UFF force field)")
        self.relax_checkbox.setChecked(True)
        self.relax_checkbox.setToolTip("Optimize the 3D conformation of the added group while holding the parent molecule rigid.")
        layout.addWidget(self.relax_checkbox)

        # Action buttons
        btn_layout = QHBoxLayout()
        self.replace_button = QPushButton("Replace Selected Atom")
        self.replace_button.setEnabled(False)
        self.replace_button.setStyleSheet("font-weight: bold;")
        self.replace_button.clicked.connect(self.replace_atom)
        btn_layout.addWidget(self.replace_button)

        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        btn_layout.addWidget(self.clear_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        btn_layout.addWidget(close_button)

        layout.addLayout(btn_layout)

        self._populate_groups()

    def _populate_groups(self) -> None:
        """Populate the group combo box based on current category and search text."""
        selected_category = self.category_combo.currentText()
        query = self.search_input.text().strip()

        filtered = search_groups(query, category=selected_category)

        self.group_combo.blockSignals(True)
        current = self.group_combo.currentText()
        self.group_combo.clear()
        self.group_combo.addItems(filtered)
        if current in filtered:
            self.group_combo.setCurrentText(current)
        elif filtered:
            self.group_combo.setCurrentIndex(0)
        self.group_combo.blockSignals(False)

        self._update_group_preview(self.group_combo.currentText())


    def _on_category_changed(self, _category: str) -> None:
        self._populate_groups()

    def _on_search_changed(self, _text: str) -> None:
        self._populate_groups()

    def _update_group_preview(self, group_name: str) -> None:
        smi = get_group_smiles(group_name)
        if smi:
            self.preview_label.setText(f"SMILES: {smi}")
        else:
            self.preview_label.setText("")

    def _position_near_parent(self) -> None:
        """Place the tool predictably near the host window and on-screen."""
        parent = self.parentWidget()
        if parent is None:
            return
        parent_rect = parent.frameGeometry()
        screen = parent.screen() or QApplication.primaryScreen()
        if screen is None:
            return
        bounds = screen.availableGeometry()
        x = parent_rect.center().x() - self.width() // 2
        y = parent_rect.center().y() - self.height() // 2
        x = max(bounds.left(), min(x, bounds.right() - self.width()))
        y = max(bounds.top(), min(y, bounds.bottom() - self.height()))
        self.move(x, y)

    def _install_3d_picking(self) -> None:
        """Install interactive pick filter on the 3D plotter interactor (default ON)."""
        plotter = getattr(self.context, "plotter", None)
        interactor = getattr(plotter, "interactor", None) if plotter else None
        if interactor is not None:
            self._pick_filter = _AtomPickFilter(self._pick_atom, self)
            interactor.installEventFilter(self._pick_filter)

    def _pick_atom(self, x: int, y: int, widget: Any) -> None:
        """Handle 3D viewport pick event to select target atom."""
        try:
            import vtk
        except ImportError:
            return

        plotter = getattr(self.context, "plotter", None)
        mol = self.context.current_mol
        if plotter is None or mol is None or not mol.GetNumConformers():
            return

        picker = vtk.vtkCellPicker()
        ratio = widget.devicePixelRatioF() if hasattr(widget, "devicePixelRatioF") else 1.0
        picker.SetTolerance(0.005)
        picker.Pick(x * ratio, (widget.height() - y) * ratio, 0, plotter.renderer)
        picked_actor = picker.GetActor()

        main_window = self.context.get_main_window()
        view_3d = getattr(main_window, "view_3d_manager", None) if main_window else None
        atom_actor = getattr(view_3d, "atom_actor", None) if view_3d else None
        if atom_actor is not None and picked_actor is not atom_actor:
            return

        pos = picker.GetPickPosition()
        closest_atom = min(
            mol.GetAtoms(),
            key=lambda candidate: self._distance_sq(mol, candidate.GetIdx(), pos),
        )
        if self._distance_sq(mol, closest_atom.GetIdx(), pos) > (0.45 ** 2):
            return

        atom_idx = closest_atom.GetIdx()
        # Toggle selection if the same atom is clicked
        if self.selected_atom_idx == atom_idx:
            self.clear_selection()
            self.context.show_status_message("Atom selection cleared.")
            return

        self.selected_atom_idx = atom_idx
        self.update_selection_display()
        self.context.show_status_message(
            f"Selected atom {atom_idx} ({closest_atom.GetSymbol()}) for replacement."
        )

    def _distance_sq(self, mol: Any, index: int, pos: Sequence[float]) -> float:
        point = mol.GetConformer().GetAtomPosition(index)
        return (point.x - pos[0]) ** 2 + (point.y - pos[1]) ** 2 + (point.z - pos[2]) ** 2

    def update_selection_display(self) -> None:
        """Update the dialog UI label and 3D viewport point label."""
        mol = self.context.current_mol
        if self.selected_atom_idx is None or mol is None or self.selected_atom_idx >= mol.GetNumAtoms():
            self.selection_label.setText("No atom selected")
            self.replace_button.setEnabled(False)
            self.clear_selection_labels()
            return

        atom = mol.GetAtomWithIdx(self.selected_atom_idx)
        self.selection_label.setText(
            f"Selected atom: {atom.GetSymbol()}{atom.GetIdx()} (index {atom.GetIdx()})"
        )
        self.replace_button.setEnabled(True)
        self.show_atom_labels()

    def add_selection_label(self, atom_idx: int, label_text: str = "1", color: str = "yellow") -> None:
        """Add a 3D point label in the PyVista viewport matching main app style."""
        plotter = getattr(self.context, "plotter", None)
        if plotter is None or not hasattr(plotter, "add_point_labels"):
            return

        mol = self.context.current_mol
        if mol is None or not mol.GetNumConformers() or atom_idx >= mol.GetNumAtoms():
            return

        cam = None
        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            pass

        conf = mol.GetConformer()
        pt = conf.GetAtomPosition(atom_idx)
        pos = [pt.x, pt.y, pt.z]

        try:
            label_actor = plotter.add_point_labels(
                [pos],
                [label_text],
                point_size=0,
                font_size=12,
                text_color=color,
                always_visible=True,
                show_points=False,
                shape="rect",
                shape_color="gray",
                shape_opacity=0.5,
            )
            self.selection_labels.append(label_actor)

            if cam is not None:
                try:
                    plotter.camera_position = cam
                except (AttributeError, RuntimeError, TypeError):
                    pass

            if hasattr(plotter, "render"):
                plotter.render()
        except (AttributeError, RuntimeError, TypeError, ValueError) as exc:
            logger.debug("Could not add 3D selection label: %s", exc)

    def clear_selection_labels(self) -> None:
        """Remove all 3D label actors from the PyVista plotter."""
        plotter = getattr(self.context, "plotter", None)
        if plotter is not None and hasattr(plotter, "remove_actor"):
            for actor in self.selection_labels:
                try:
                    if actor is not None:
                        plotter.remove_actor(actor)
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    pass
            if hasattr(plotter, "render"):
                plotter.render()
        self.selection_labels.clear()

    def show_atom_labels(self) -> None:
        """Refresh 3D labels for the current selection."""
        self.clear_selection_labels()
        if self.selected_atom_idx is not None:
            self.add_selection_label(self.selected_atom_idx, "1")

    def clear_selection(self) -> None:
        """Clear the current selection in UI and 3D view."""
        self.selected_atom_idx = None
        self.clear_selection_labels()
        self.selection_label.setText("No atom selected")
        self.replace_button.setEnabled(False)

    def replace_atom(self) -> None:
        """Execute replacement of the selected atom with the chosen functional group."""
        mol = self.context.current_mol
        target = self.selected_atom_idx
        group_name = self.group_combo.currentText()
        group_smi = get_group_smiles(group_name)

        if mol is None or target is None or not group_smi:
            QMessageBox.warning(self, "No selection", "Please click an atom in the 3D view first.")
            return

        try:
            relax = self.relax_checkbox.isChecked()
            new_mol = replace_atom_with_group(mol, target, group_smi, relax=relax)
            self.context.current_molecule = new_mol
            self.context.push_undo_checkpoint()

            refresh = getattr(self.context, "refresh_3d_view", None)
            if callable(refresh):
                refresh()
            else:
                main_window = self.context.get_main_window()
                view_3d = getattr(main_window, "view_3d_manager", None) if main_window else None
                if view_3d and hasattr(view_3d, "draw_molecule_3d"):
                    view_3d.draw_molecule_3d(new_mol)

            self.clear_selection()
            self.context.show_status_message(f"Replaced atom {target} with {group_name}.")
        except (RuntimeError, ValueError, AttributeError) as exc:
            QMessageBox.critical(self, "Replacement failed", str(exc))

    def closeEvent(self, event: Any) -> None:
        """Clean up 3D viewport labels and interactor filters on close."""
        self.clear_selection_labels()
        plotter = getattr(self.context, "plotter", None)
        interactor = getattr(plotter, "interactor", None) if plotter else None
        if interactor is not None and self._pick_filter is not None:
            try:
                interactor.removeEventFilter(self._pick_filter)
            except (AttributeError, RuntimeError):
                pass
        super().closeEvent(event)
