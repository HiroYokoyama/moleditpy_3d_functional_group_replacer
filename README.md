# 3D Functional Group Toolbox

A MoleditPy plugin for replacing one selected atom with a common one-attachment functional group while retaining the atom's existing neighbours.

## Features

- Common groups including alkyl, hydroxyl, amino, carbonyl, carboxyl, cyano, nitro, aryl, and halogen groups.
- Atom selection by index in a small toolbox window.
- Undo checkpoint and 3D-view refresh after each replacement.
- Validation of one-point attachment-group definitions.

## Development

```text
python -m pytest tests -q
```

The plugin targets MoleditPy 4.x and Python 3.9–3.14.
