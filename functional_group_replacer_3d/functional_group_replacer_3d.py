"""Backward-compatibility module re-exporting the refactored toolbox components."""

from __future__ import annotations

import sys
from pathlib import Path

if not __package__:
    _pkg_root = str(Path(__file__).resolve().parent.parent)
    if _pkg_root not in sys.path:
        sys.path.insert(0, _pkg_root)
    from functional_group_replacer_3d import (
        PLUGIN_AUTHOR,
        PLUGIN_CATEGORY,
        PLUGIN_DEPENDENCIES,
        PLUGIN_DESCRIPTION,
        PLUGIN_NAME,
        PLUGIN_SUPPORTED_MOLEDITPY_VERSION,
        PLUGIN_SUPPORTED_PYTHON_VERSION,
        PLUGIN_TAGS,
        PLUGIN_VERSION,
        initialize,
    )
    from functional_group_replacer_3d.chemistry import (
        relax_molecule_with_fixed_atoms,
        replace_atom_with_group,
    )
    from functional_group_replacer_3d.dialog import (
        FunctionalGroupReplacer,
        FunctionalGroupToolbox,
        _AtomPickFilter,
    )
    from functional_group_replacer_3d.groups import (
        GROUP_CATEGORIES,
        GROUPS,
        get_group_smiles,
        get_groups_by_category,
        search_groups,
    )
else:
    from .__init__ import (
        PLUGIN_AUTHOR,
        PLUGIN_CATEGORY,
        PLUGIN_DEPENDENCIES,
        PLUGIN_DESCRIPTION,
        PLUGIN_NAME,
        PLUGIN_SUPPORTED_MOLEDITPY_VERSION,
        PLUGIN_SUPPORTED_PYTHON_VERSION,
        PLUGIN_TAGS,
        PLUGIN_VERSION,
        initialize,
    )
    from .chemistry import relax_molecule_with_fixed_atoms, replace_atom_with_group
    from .dialog import FunctionalGroupReplacer, FunctionalGroupToolbox, _AtomPickFilter
    from .groups import (
        GROUP_CATEGORIES,
        GROUPS,
        get_group_smiles,
        get_groups_by_category,
        search_groups,
    )


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
    "FunctionalGroupReplacer",
    "FunctionalGroupToolbox",
    "_AtomPickFilter",
    "get_group_smiles",
    "get_groups_by_category",
    "initialize",
    "relax_molecule_with_fixed_atoms",
    "replace_atom_with_group",
    "search_groups",
]
