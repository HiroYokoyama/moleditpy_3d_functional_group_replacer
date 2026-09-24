"""3D Functional Group Replacer plugin for MoleditPy."""

from __future__ import annotations

from typing import Any

from .chemistry import relax_molecule_with_fixed_atoms, replace_atom_with_group
from .dialog import (
    DEFAULT_SETTINGS,
    WINDOW_ID,
    FunctionalGroupReplacer,
    FunctionalGroupToolbox,
)
from .groups import (
    GROUP_ALIASES,
    GROUP_CATEGORIES,
    GROUPS,
    get_group_smiles,
    get_groups_by_category,
    search_groups,
)

PLUGIN_NAME = "3D Functional Group Replacer"
PLUGIN_VERSION = "0.9.0"
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Replace a selected atom with a common functional group while retaining its neighbours."
PLUGIN_DEPENDENCIES = ["rdkit", "PyQt6"]
PLUGIN_CATEGORY = "3D Editing"
PLUGIN_TAGS = ["Utility"]

_context: Any | None = None
_dialog_opened: bool = False
_current_settings: dict[str, Any] = dict(DEFAULT_SETTINGS)


def _is_widget_alive(widget: Any) -> bool:
    """Check whether a Qt widget instance is valid and its C/C++ object has not been deleted."""
    if widget is None:
        return False
    try:
        from PyQt6 import sip

        if sip.isdeleted(widget):
            return False
    except (ImportError, TypeError):
        pass
    try:
        # Calling a simple Qt method confirms whether underlying C++ object is intact
        _ = widget.isVisible()
        return True
    except (RuntimeError, AttributeError):
        return False


def _live_dialog() -> Any | None:
    """Return the open replacer dialog, or None if there is none."""
    if _context is None:
        return None
    dlg = _context.get_window(WINDOW_ID)
    return dlg if _is_widget_alive(dlg) and hasattr(dlg, "apply_settings") else None


def _remember_settings(settings: dict[str, Any]) -> None:
    _current_settings.update(settings)


def _open_replacer() -> None:
    """Open or focus the 3D Functional Group Replacer dialog."""
    global _dialog_opened
    if _context is None:
        return

    _dialog_opened = True
    window = _live_dialog()
    if window is None:
        window = FunctionalGroupReplacer(_context)
        window.apply_settings(_current_settings)
        # Track every change so the choices survive the dialog being closed.
        window.settings_changed.connect(_remember_settings)
    window.show()
    window.raise_()
    window.activateWindow()


def _save_state() -> dict[str, Any]:
    if not _dialog_opened:
        return {}
    return {"settings": dict(_current_settings)}


def _load_state(data: Any) -> None:
    saved = data.get("settings") if isinstance(data, dict) else None
    if not isinstance(saved, dict):
        return
    # Only take known keys with the expected types from the project file.
    for key, default in DEFAULT_SETTINGS.items():
        if isinstance(saved.get(key), type(default)):
            _current_settings[key] = saved[key]
    dlg = _live_dialog()
    if dlg is not None:
        dlg.apply_settings(_current_settings)


def _reset_state() -> None:
    global _dialog_opened
    dlg = _live_dialog()
    if dlg is not None and dlg.isVisible():
        return
    _dialog_opened = False
    _current_settings.clear()
    _current_settings.update(DEFAULT_SETTINGS)


def initialize(context: Any) -> None:
    """Initialize the plugin within MoleditPy host context."""
    global _context
    _context = context

    context.add_menu_action("3D Edit/3D Functional Group Replacer...", _open_replacer)
    context.register_save_handler(_save_state)
    context.register_load_handler(_load_state)
    context.register_document_reset_handler(_reset_state)


__all__ = [
    "GROUPS",
    "GROUP_ALIASES",
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
    "FunctionalGroupReplacer",
    "FunctionalGroupToolbox",
    "get_group_smiles",
    "get_groups_by_category",
    "initialize",
    "relax_molecule_with_fixed_atoms",
    "replace_atom_with_group",
    "search_groups",
]
