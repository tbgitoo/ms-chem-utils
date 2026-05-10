from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .formula import counter_to_formula, formula_to_counter
from .fragmentation import _mol_formula_counter, _non_h_counter, _contains_non_h_elements
from .rdkit_query import (
    make_generic_bond_query,
    unsat_units_from_mol,
    heavy_heavy_bond_count,
    submol_from_match,
)

def smiles_fits_into_smiles(
    child_smiles: str,
    parent_smiles: str,
    *,
    h_tolerance_abs: int = 2,
    h_tolerance_fraction_of_parent: float = 0.20,
    bond_tolerance_abs: int = 2,
    bond_tolerance_fraction_of_fragment: float = 0.20,
    max_substruct_matches_to_check: int = 200,
) -> Dict[str, Any]:
    """
    Determine whether child_smiles can fit into parent_smiles as a heavy-atom subgraph,
    allowing:
      - Hydrogen mismatch: +/- max(2, ceil(0.20 * H_parent))
      - Bond-order/unsaturation quirks: child can have additional unsaturation units relative to
        the matched parent region, but:
          * unsaturation must not decrease (child >= parent_fragment)
          * additional unsaturation is limited to max(2, ceil(0.20 * heavy-heavy bonds in fragment))

    Returns a structured dict, never throws for normal bad inputs.
    """

    child_smiles = (child_smiles or "").strip()
    parent_smiles = (parent_smiles or "").strip()

    if not child_smiles:
        return {"ok": False, "fits": False, "reason": "Empty child SMILES"}
    if not parent_smiles:
        return {"ok": False, "fits": False, "reason": "Empty parent SMILES"}

    child = Chem.MolFromSmiles(child_smiles)
    if child is None:
        return {"ok": False, "fits": False, "reason": "Invalid child SMILES"}

    parent = Chem.MolFromSmiles(parent_smiles)
    if parent is None:
        return {"ok": False, "fits": False, "reason": "Invalid parent SMILES"}

    # ---- Stage A: formula containment (non-H) ----
    child_formula = _mol_formula_counter(child)
    parent_formula = _mol_formula_counter(parent)

    child_nonh = _non_h_counter(child_formula)
    parent_nonh = _non_h_counter(parent_formula)

    if not _contains_non_h_elements(parent_nonh, child_nonh):
        return {
            "ok": True,
            "fits": False,
            "reason": "Non-H elemental containment fails (child introduces new element or exceeds counts)",
            "child_formula": counter_to_formula(child_formula),
            "parent_formula": counter_to_formula(parent_formula),
        }

    # Hydrogen tolerance based on parent H count
    h_parent = int(parent_formula.get("H", 0))
    h_tol = max(int(h_tolerance_abs), int(math.ceil(h_tolerance_fraction_of_parent * h_parent)))

    # ---- Stage B: substructure match with generic bonds ----
    query = make_generic_bond_query(child)

    matches = parent.GetSubstructMatches(query, uniquify=True)
    if not matches:
        return {
            "ok": True,
            "fits": False,
            "reason": "No substructure match (even allowing generic bonds)",
            "child_formula": counter_to_formula(child_formula),
            "parent_formula": counter_to_formula(parent_formula),
            "H_tolerance": h_tol,
        }

    # Precompute child unsaturation units (heavy-heavy bonds)
    child_unsat = unsat_units_from_mol(child, atom_set=None)
    child_H = int(child_formula.get("H", 0))

    checked = 0
    passing_examples: List[Dict[str, Any]] = []

    for m in matches:
        checked += 1
        if checked > max_substruct_matches_to_check:
            break

        atom_set = set(m)

        sub = submol_from_match(parent, m)
        if sub is None:
            continue

        sub_h = Chem.AddHs(sub)
        frag_formula = formula_to_counter(rdMolDescriptors.CalcMolFormula(sub_h))
        frag_nonh = _non_h_counter(frag_formula)
        frag_H = int(frag_formula.get("H", 0))

        if frag_nonh != child_nonh:
            continue

        if abs(frag_H - child_H) > h_tol:
            continue

        parent_frag_unsat = unsat_units_from_mol(parent, atom_set=atom_set)
        hh_bonds = heavy_heavy_bond_count(parent, atom_set)
        unsat_tol = max(int(bond_tolerance_abs), int(math.ceil(bond_tolerance_fraction_of_fragment * hh_bonds)))

        if child_unsat < parent_frag_unsat:
            continue

        if (child_unsat - parent_frag_unsat) > unsat_tol:
            continue

        example = {
            "match_atom_indices_in_parent": list(m),
            "fragment_formula_capped": counter_to_formula(frag_formula),
            "child_formula": counter_to_formula(child_formula),
            "parent_formula": counter_to_formula(parent_formula),
            "frag_H": frag_H,
            "child_H": child_H,
            "H_tolerance": h_tol,
            "parent_frag_unsat_units": parent_frag_unsat,
            "child_unsat_units": child_unsat,
            "unsat_tolerance": unsat_tol,
            "fragment_heavy_heavy_bonds": hh_bonds,
        }
        passing_examples.append(example)

        return {
            "ok": True,
            "fits": True,
            "reason": "Found a mapping satisfying element containment, H tolerance, and unsaturation constraints",
            "checked_matches": checked,
            "total_matches_found": len(matches),
            "example": example,
        }

    return {
        "ok": True,
        "fits": False,
        "reason": "Substructure matches exist, but none satisfy H tolerance and/or unsaturation constraints",
        "child_formula": counter_to_formula(child_formula),
        "parent_formula": counter_to_formula(parent_formula),
        "H_tolerance": h_tol,
        "child_unsat_units": child_unsat,
        "checked_matches": checked,
        "total_matches_found": len(matches),
        "first_few_passing_examples": passing_examples[:3],
    }