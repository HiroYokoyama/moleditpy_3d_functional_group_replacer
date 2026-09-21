"""3D Functional Group Toolbox plugin for MoleditPy."""

from __future__ import annotations

from typing import Any

from .chemistry import relax_molecule_with_fixed_atoms, replace_atom_with_group
from .dialog import FunctionalGroupToolbox
from .groups import (
    GROUP_CATEGORIES,
    GROUPS,
    get_group_smiles,
    get_groups_by_category,
    search_groups,
)

PLUGIN_NAME = "3D Functional Group Toolbox"
PLUGIN_VERSION = "0.4.1"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Replace a selected atom with a common functional group while retaining its neighbours."
PLUGIN_DEPENDENCIES = ["rdkit", "PyQt6"]
PLUGIN_CATEGORY = "3D Editing"
PLUGIN_TAGS = ["3D", "Editing", "Chemistry"]

WINDOW_ID = "functional_group_toolbox"
_context: Any | None = None
_dialog_opened: bool = False
_current_settings: dict[str, Any] = {
    "last_category": "All",
    "last_group": "Methyl",
    "relax": True,
}


def _open_toolbox() -> None:
    """Open or bring to front the 3D Functional Group Toolbox window."""
    global _dialog_opened
    if _context is None:
        return

    existing = _context.get_window(WINDOW_ID)
    if existing is not None and existing.isVisible():
        existing.raise_()
        existing.activateWindow()
        return

    _dialog_opened = True
    window = FunctionalGroupToolbox(_context)
    if _current_settings.get("last_category") in GROUP_CATEGORIES:
        window.category_combo.setCurrentText(_current_settings["last_category"])
    if _current_settings.get("last_group") in GROUPS:
        window.group_combo.setCurrentText(_current_settings["last_group"])
    window.relax_checkbox.setChecked(_current_settings.get("relax", True))
    window.show()
    window.raise_()
    window.activateWindow()


def initialize(context: Any) -> None:
    """Initialize the plugin within MoleditPy host context."""
    global _context
    _context = context

    context.add_menu_action("Edit/3D Functional Group Toolbox...", _open_toolbox)

    def save_state() -> dict[str, Any]:
        if not _dialog_opened:
            return {}
        dlg = context.get_window(WINDOW_ID)
        if dlg is not None and hasattr(dlg, "group_combo"):
            _current_settings["last_category"] = dlg.category_combo.currentText()
            _current_settings["last_group"] = dlg.group_combo.currentText()
            _current_settings["relax"] = dlg.relax_checkbox.isChecked()
        return {"settings": dict(_current_settings)}

    def load_state(data: Any) -> None:
        if isinstance(data, dict):
            saved = data.get("settings")
            if isinstance(saved, dict):
                _current_settings.update(saved)
                dlg = context.get_window(WINDOW_ID)
                if dlg is not None and hasattr(dlg, "group_combo"):
                    if _current_settings.get("last_category") in GROUP_CATEGORIES:
                        dlg.category_combo.setCurrentText(_current_settings["last_category"])
                    if _current_settings.get("last_group") in GROUPS:
                        dlg.group_combo.setCurrentText(_current_settings["last_group"])
                    dlg.relax_checkbox.setChecked(_current_settings.get("relax", True))

    def reset_state() -> None:
        global _dialog_opened
        dlg = context.get_window(WINDOW_ID)
        if dlg is not None and dlg.isVisible():
            return
        _dialog_opened = False
        _current_settings["last_category"] = "All"
        _current_settings["last_group"] = "Methyl"
        _current_settings["relax"] = True

    if hasattr(context, "register_save_handler"):
        context.register_save_handler(save_state)
    if hasattr(context, "register_load_handler"):
        context.register_load_handler(load_state)
    if hasattr(context, "register_document_reset_handler"):
        context.register_document_reset_handler(reset_state)


__all__ = [
    "GROUPS",
    "GROUP_CATEGORIES",
    "PLUGIN_AUTHOR",
    "PLUGIN_CATEGORY",
    "PLUGIN_DEPENDENCIES",
    "PLUGIN_DESCRIPTION",
    "PLUGIN_NAME",
    "PLUGIN_SUPPORTED_MOLEDITPY_VERSION",
    "PLUGIN_SUPPORTED_PYTHON_VERSION",
    "PLUGIN_TAGS",
    "PLUGIN_VERSION",
    "FunctionalGroupToolbox",
    "get_group_smiles",
    "get_groups_by_category",
    "initialize",
    "relax_molecule_with_fixed_atoms",
    "replace_atom_with_group",
    "search_groups",
]
