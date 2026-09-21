# 3D Functional Group Replacer for MoleditPy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22866924.svg)](https://doi.org/10.5281/zenodo.22866924)
[![CI](https://github.com/HiroYokoyama/moleditpy_3d_functional_group_toolbox/actions/workflows/ci.yml/badge.svg)](https://github.com/HiroYokoyama/moleditpy_3d_functional_group_toolbox/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![MoleditPy](https://img.shields.io/badge/MoleditPy->=4.0.0-orange)](https://github.com/HiroYokoyama/python_molecular_editor)

An interactive, high-precision 3D molecular editing plugin for **[MoleditPy](https://github.com/HiroYokoyama/python_molecular_editor)**. Substitute any selected atom (such as a hydrogen or terminal group) with common functional groups in 3D space with automatic bond alignment and constrained force-field relaxation.

---

## Features

- **Interactive 3D Picking by Default**: Simply click any atom in the 3D viewport to select it as the target. No manual atom index lookups or cumbersome toggle buttons required.
- **Native 3D Selected Atom Label**: Displays an on-screen yellow numbered point label (`"1"`) directly attached to the selected atom in the 3D view (identical to MoleditPy's native 3D geometry editing dialogs). Deselecting, replacing, or closing automatically removes the label.
- **50+ Curated Functional Groups**: Comprehensive library spanning:
  - *Alkyl & Aliphatic*: Methyl, Ethyl, n-Propyl, Isopropyl, n-Butyl, sec-Butyl, Isobutyl, tert-Butyl, Neopentyl, Cyclopropyl, Cyclopentyl, Cyclohexyl, Trifluoromethyl
  - *Alkenyl & Alkynyl*: Vinyl, Allyl, Ethynyl, Propargyl
  - *Aryl & Heteroaryl*: Phenyl, 4-Tolyl, 4-Methoxyphenyl, 4-Fluorophenyl, 4-Chlorophenyl, 4-Nitrophenyl, Benzyl, 2-Pyridyl, 3-Pyridyl, 4-Pyridyl, 2-Thienyl, 2-Furyl
  - *Oxygen & Carbonyl*: Hydroxyl, Methoxy, Ethoxy, Phenoxy, Formyl, Acetyl, Carboxyl, Methoxycarbonyl, Ethoxycarbonyl, Carbamoyl, Acetoxy, Trifluoroacetyl
  - *Nitrogen & Amine*: Amino, Methylamino, Dimethylamino, Acetamido, Cyano, Nitro, Azido, Isocyanato, Isothiocyanato
  - *Halogen*: Fluoro, Chloro, Bromo, Iodo
  - *Sulfur & Phosphorus*: Thiol, Methylsulfanyl, Methylsulfinyl, Methylsulfonyl, Sulfo, Sulfamoyl, Triflyl, Phosphono
- **Real-Time Search by Name or SMILES**: Filter groups instantly by category dropdown or by typing in the search box (e.g. typing `phenyl`, `amino`, or SMILES like `c1ccccc1`, `C(=O)O`).
- **Physically Realistic 3D Conformation**: Calculates the original bond vector from neighbor to target atom, aligns the attachment bond along this vector, and performs constrained **MMFF94 / UFF force field minimization** holding all parent atoms completely rigid.
- **Full Undo/Redo & State Persistence**: Pushes state checkpoints to MoleditPy's undo stack and preserves your last-used category, group, and relaxation options between sessions.

---

## Installation

### Via MoleditPy Plugin Manager (Recommended)
1. Open MoleditPy.
2. Navigate to **Plugins** > **Plugin Manager...**.
3. Locate **3D Functional Group Replacer** and click **Install**.

### Manual Installation
Clone or download this repository into your MoleditPy plugins directory:

```bash
# On Windows:
git clone https://github.com/HiroYokoyama/moleditpy_3d_functional_group_toolbox.git %USERPROFILE%\.moleditpy\plugins\functional_group_replacer_3d

# On Linux/macOS:
git clone https://github.com/HiroYokoyama/moleditpy_3d_functional_group_toolbox.git ~/.moleditpy/plugins/functional_group_replacer_3d
```

Restart MoleditPy or run **Plugins** > **Reload Plugins**.

---

## How to Use

1. Open a molecule in MoleditPy and enter the **3D View**.
2. Open the tool from **3D Edit** > **3D Functional Group Replacer...**.
3. **Select an Atom**: Click any atom in the 3D viewport (e.g., a hydrogen atom or halogen).
   - A yellow `"1"` label immediately pins to the selected atom in the 3D scene.
   - The dialog displays `Selected atom: C0 (index 0)`.
   - Clicking the atom again or pressing **Clear Selection** deselects it.
4. **Choose a Functional Group**:
   - Filter by **Category** or type in the **Search** box (supports group names and SMILES).
   - Select the desired functional group from the dropdown list.
5. **Adjust Options**:
   - Keep **Relax group 3D geometry (MMFF/UFF)** checked for optimized geometry.
6. **Apply**: Click **Replace Selected Atom**.
   - The selected atom is substituted, the group is embedded in 3D, and the view updates automatically with full undo support (`Ctrl+Z`).

---

## Architecture

```text
moleditpy_3d_functional_group_toolbox/
├── functional_group_replacer_3d/
│   ├── __init__.py                # Package metadata, lifecycle, and menu registration
│   ├── chemistry.py               # Core replacement algorithm & constrained 3D relaxation
│   ├── dialog.py                  # PyQt6 UI, 3D pick event filter & PyVista label management
│   ├── groups.py                  # 50+ SMILES definitions and category index
│   └── functional_group_replacer_3d.py # Backward-compatible re-export module
├── tests/
│   ├── test_functional_group_replacer_3d.py # Unit and GUI test suite
│   ├── test_api.py                # Static AST contract check against MoleditPy
│   └── plugin_api_checker.py      # AST API analysis tool
├── .github/
│   └── workflows/
│       ├── ci.yml                 # Multi-platform CI (Linux, Windows)
│       └── release.yml            # Automated GitHub release packaging
├── AGENTS.md                      # AI agent & contributor workflow guidelines
├── CONTRIBUTING.md                # Contribution guidelines
├── CODE_OF_CONDUCT.md             # Code of conduct
├── LICENSE                        # GNU General Public License v3.0
└── pyproject.toml                 # Package configuration
```


---

## Testing

Run the headless test suite:

```bash
python -m pytest tests/ -v
```

Run the API checker against the MoleditPy main app:

```bash
python -m pytest tests/test_api.py -v
```

---

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
