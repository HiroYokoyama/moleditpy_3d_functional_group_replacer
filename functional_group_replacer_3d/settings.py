"""User preferences in the installer-preserved companion settings.json."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any

from .groups import GROUP_CATEGORIES, GROUPS

logger = logging.getLogger(__name__)
SETTINGS_FILE = Path(__file__).resolve().parent / "settings.json"
DEFAULT_SETTINGS: dict[str, Any] = {
    "last_category": "All",
    "last_group": "Methyl",
    "relax": True,
    "replace_terminal_hydrogen": True,
    "click_to_replace": False,
}


def _validated_settings(data: Any) -> dict[str, Any]:
    settings = dict(DEFAULT_SETTINGS)
    if not isinstance(data, dict):
        return settings
    for key, default in DEFAULT_SETTINGS.items():
        value = data.get(key)
        if not isinstance(value, type(default)):
            continue
        if key == "last_category" and value not in GROUP_CATEGORIES:
            continue
        if key == "last_group" and value not in GROUPS:
            continue
        settings[key] = value
    return settings


def load_settings() -> dict[str, Any]:
    """Load known, valid choices; missing or unreadable files use defaults."""
    try:
        with SETTINGS_FILE.open(encoding="utf-8") as handle:
            return _validated_settings(json.load(handle))
    except FileNotFoundError:
        return dict(DEFAULT_SETTINGS)
    except (OSError, ValueError) as exc:
        logger.warning("Could not read %s; using defaults: %s", SETTINGS_FILE, exc)
        return dict(DEFAULT_SETTINGS)


def save_settings(settings: dict[str, Any]) -> None:
    """Atomically save preferences without truncating a previous valid file."""
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n",
            dir=SETTINGS_FILE.parent, suffix=".tmp", delete=False,
        ) as handle:
            temp_path = Path(handle.name)
            json.dump(_validated_settings(settings), handle, indent=2)
            handle.write("\n")
        os.replace(temp_path, SETTINGS_FILE)
    except OSError as exc:
        logger.warning("Could not write %s: %s", SETTINGS_FILE, exc)
    finally:
        if temp_path is not None:
            try:
                temp_path.unlink(missing_ok=True)
            except OSError:
                logger.debug("Could not remove temporary settings file", exc_info=True)
