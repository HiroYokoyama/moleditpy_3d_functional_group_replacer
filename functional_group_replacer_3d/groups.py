"""Functional groups library and categorization for 3D Functional Group Toolbox."""

from __future__ import annotations

# Comprehensive library of 1-attachment-point functional groups
# The dummy atom [*:1] indicates the attachment point that connects to the parent molecule.
GROUPS: dict[str, str] = {
    # Alkyl & Aliphatic
    "Methyl": "[*:1]C",
    "Ethyl": "[*:1]CC",
    "n-Propyl": "[*:1]CCC",
    "Isopropyl": "[*:1]C(C)C",
    "n-Butyl": "[*:1]CCCC",
    "sec-Butyl": "[*:1]C(C)CC",
    "Isobutyl": "[*:1]CC(C)C",
    "tert-Butyl": "[*:1]C(C)(C)C",
    "Neopentyl": "[*:1]CC(C)(C)C",
    "Cyclopropyl": "[*:1]C1CC1",
    "Cyclopentyl": "[*:1]C1CCCC1",
    "Cyclohexyl": "[*:1]C1CCCCC1",
    "Trifluoromethyl": "[*:1]C(F)(F)F",
    # Alkenyl & Alkynyl
    "Vinyl": "[*:1]C=C",
    "Allyl": "[*:1]CC=C",
    "Ethynyl": "[*:1]C#C",
    "Propargyl": "[*:1]CC#C",
    # Aryl & Heteroaryl
    "Phenyl": "[*:1]c1ccccc1",
    "4-Tolyl": "[*:1]c1ccc(C)cc1",
    "4-Methoxyphenyl": "[*:1]c1ccc(OC)cc1",
    "4-Fluorophenyl": "[*:1]c1ccc(F)cc1",
    "4-Chlorophenyl": "[*:1]c1ccc(Cl)cc1",
    "4-Nitrophenyl": "[*:1]c1ccc([N+](=O)[O-])cc1",
    "Benzyl": "[*:1]Cc1ccccc1",
    "2-Pyridyl": "[*:1]c1ncccc1",
    "3-Pyridyl": "[*:1]c1cnccc1",
    "4-Pyridyl": "[*:1]c1ccncc1",
    "2-Thienyl": "[*:1]c1sccc1",
    "2-Furyl": "[*:1]c1occc1",
    # Oxygen & Carbonyl
    "Hydroxyl": "[*:1]O",
    "Methoxy": "[*:1]OC",
    "Ethoxy": "[*:1]OCC",
    "Phenoxy": "[*:1]Oc1ccccc1",
    "Formyl": "[*:1]C=O",
    "Acetyl": "[*:1]C(=O)C",
    "Carboxyl": "[*:1]C(=O)O",
    "Methoxycarbonyl": "[*:1]C(=O)OC",
    "Ethoxycarbonyl": "[*:1]C(=O)OCC",
    "Carbamoyl": "[*:1]C(=O)N",
    "Acetoxy": "[*:1]OC(=O)C",
    "Trifluoroacetyl": "[*:1]C(=O)C(F)(F)F",
    "Hydroperoxyl": "[*:1]OO",
    # Nitrogen & Amine
    "Amino": "[*:1]N",
    "Methylamino": "[*:1]NC",
    "Dimethylamino": "[*:1]N(C)C",
    "Acetamido": "[*:1]NC(=O)C",
    "Cyano": "[*:1]C#N",
    "Nitrile": "[*:1]C#N",
    "Nitro": "[*:1][N+](=O)[O-]",
    "Nitroso": "[*:1]N=O",
    "Hydrazinyl": "[*:1]NN",
    "Ureido": "[*:1]NC(=O)N",
    "Guanidino": "[*:1]NC(=N)N",
    "Azido": "[*:1][N]=[N+]=[N-]",
    "Cyanato": "[*:1]OC#N",
    "Isocyanato": "[*:1]N=C=O",
    "Isothiocyanato": "[*:1]N=C=S",
    # Halogen
    "Fluoro": "[*:1]F",
    "Chloro": "[*:1]Cl",
    "Bromo": "[*:1]Br",
    "Iodo": "[*:1]I",
    # Sulfur, Phosphorus & Boron
    "Thiol": "[*:1]S",
    "Methylsulfanyl": "[*:1]SC",
    "Methylsulfinyl": "[*:1]S(=O)C",
    "Methylsulfonyl": "[*:1]S(=O)(=O)C",
    "Sulfo": "[*:1]S(=O)(=O)O",
    "Sulfamoyl": "[*:1]S(=O)(=O)N",
    "Triflyl": "[*:1]S(=O)(=O)C(F)(F)F",
    "Thiocyanato": "[*:1]SC#N",
    "Phosphono": "[*:1]P(=O)(O)O",
    "Boryl": "[*:1]B(O)O",
}

GROUP_CATEGORIES: dict[str, list[str]] = {
    "All": list(GROUPS.keys()),
    "Alkyl & Aliphatic": [
        "Methyl",
        "Ethyl",
        "n-Propyl",
        "Isopropyl",
        "n-Butyl",
        "sec-Butyl",
        "Isobutyl",
        "tert-Butyl",
        "Neopentyl",
        "Cyclopropyl",
        "Cyclopentyl",
        "Cyclohexyl",
        "Trifluoromethyl",
    ],
    "Alkenyl & Alkynyl": [
        "Vinyl",
        "Allyl",
        "Ethynyl",
        "Propargyl",
    ],
    "Aryl & Heteroaryl": [
        "Phenyl",
        "4-Tolyl",
        "4-Methoxyphenyl",
        "4-Fluorophenyl",
        "4-Chlorophenyl",
        "4-Nitrophenyl",
        "Benzyl",
        "2-Pyridyl",
        "3-Pyridyl",
        "4-Pyridyl",
        "2-Thienyl",
        "2-Furyl",
    ],
    "Oxygen & Carbonyl": [
        "Hydroxyl",
        "Methoxy",
        "Ethoxy",
        "Phenoxy",
        "Formyl",
        "Acetyl",
        "Carboxyl",
        "Methoxycarbonyl",
        "Ethoxycarbonyl",
        "Carbamoyl",
        "Acetoxy",
        "Trifluoroacetyl",
        "Hydroperoxyl",
    ],
    "Nitrogen & Amine": [
        "Amino",
        "Methylamino",
        "Dimethylamino",
        "Acetamido",
        "Cyano",
        "Nitrile",
        "Nitro",
        "Nitroso",
        "Hydrazinyl",
        "Ureido",
        "Guanidino",
        "Azido",
        "Cyanato",
        "Isocyanato",
        "Isothiocyanato",
    ],
    "Halogen": [
        "Fluoro",
        "Chloro",
        "Bromo",
        "Iodo",
    ],
    "Sulfur & Phosphorus": [
        "Thiol",
        "Methylsulfanyl",
        "Methylsulfinyl",
        "Methylsulfonyl",
        "Sulfo",
        "Sulfamoyl",
        "Triflyl",
        "Thiocyanato",
        "Phosphono",
        "Boryl",
    ],
}


def get_group_smiles(name: str) -> str | None:
    """Retrieve the SMILES string for a named functional group."""
    return GROUPS.get(name)


def get_groups_by_category(category: str) -> list[str]:
    """Retrieve a list of functional group names in the given category."""
    return GROUP_CATEGORIES.get(category, list(GROUPS.keys()))


def search_groups(query: str, category: str = "All") -> list[str]:
    """Filter group names matching query by name or SMILES."""
    q = query.strip()
    pool = GROUP_CATEGORIES.get(category, list(GROUPS.keys()))
    if not q:
        return pool

    q_lower = q.lower()

    # Try parsing query as SMILES with RDKit to compare canonical SMILES
    can_query = None
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(q) or Chem.MolFromSmiles(f"[*:1]{q}")
        if mol:
            # strip dummy for canonical comparison
            for a in mol.GetAtoms():
                a.SetAtomMapNum(0)
            can_query = Chem.MolToSmiles(mol)
    except (ValueError, TypeError, RuntimeError):
        can_query = None

    exact_matches = []
    sub_matches = []

    for name in pool:
        smi = GROUPS.get(name, "")
        smi_clean = smi.replace("[*:1]", "")

        # 1. Canonical SMILES match via RDKit
        if can_query:
            try:
                from rdkit import Chem

                m = Chem.MolFromSmiles(smi)
                if m:
                    for a in m.GetAtoms():
                        a.SetAtomMapNum(0)
                    if Chem.MolToSmiles(m) == can_query:
                        exact_matches.append(name)
                        continue
                    rw = Chem.RWMol(m)
                    for d in [
                        a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0
                    ]:
                        rw.RemoveAtom(d)
                    if Chem.MolToSmiles(rw.GetMol()) == can_query:
                        exact_matches.append(name)
                        continue
            except (ValueError, TypeError, RuntimeError):
                pass

        # 2. Exact SMILES match (case-sensitive for SMILES)
        if q == smi_clean or q == smi:
            exact_matches.append(name)
            continue

        # 3. Name match
        if q_lower in name.lower():
            sub_matches.append(name)
            continue

        # 4. Case-sensitive substring in SMILES
        if len(q) >= 2 and (q in smi or q in smi_clean):
            sub_matches.append(name)

    # Return exact SMILES matches first, then name/substring matches
    seen = set()
    result = []
    for item in exact_matches + sub_matches:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result
