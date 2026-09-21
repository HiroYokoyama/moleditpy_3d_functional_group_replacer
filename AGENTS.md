# Contributor & Agent Guide

## Project Scope

This repository is a MoleditPy plugin for interactive 3D functional group replacement in molecular geometries. It allows users to click an atom in the 3D viewport and substitute it with one of 50+ common functional groups while preserving connectivity, orienting the group along the original bond vector, and performing constrained MMFF94/UFF force field relaxation.

## Development Setup

Install dependencies in your Python environment:

```bash
python -m pip install pytest pytest-cov numpy PyQt6 rdkit
```

## Required Checks

Run the test suite:

```bash
python -m pytest tests/ -v
```

Run API compatibility check against the main app:

```bash
python -m pytest tests/test_api.py -v
```

## Change Expectations

- Keep interactive 3D atom selection active by default; do not introduce manual index comboboxes or selection toggle buttons.
- Ensure 3D viewport point labels (`"1"`) are displayed on selected atoms and properly cleaned up when deselected, replaced, or closed.
- Do not assert hardcoded version numbers in tests.
- Keep `functional_group_toolbox_3d.py` as a backward-compatible shim re-exporting the modular package components.
- Ensure all new functional groups have valid SMILES with a single `[*:1]` attachment point.

## Commits

Make focused commits and use the repository default Git identity. Add this trailer to each commit:

`Assisted-by: Gemini 3.8 Flash`
