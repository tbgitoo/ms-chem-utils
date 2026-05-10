from __future__ import annotations

import math
from collections import Counter, deque
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .formula import formula_to_counter, counter_to_formula

def _mol_formula_counter(mol: Chem.Mol) -> Counter:
    """RDKit formula counter using AddHs() so H count is explicit."""
    mol_h = Chem.AddHs(mol)
    f = rdMolDescriptors.CalcMolFormula(mol_h)
    return formula_to_counter(f)

def _non_h_counter(c: Counter) -> Counter:
    """Copy of Counter without 'H'."""
    out = Counter(c)
    out.pop("H", None)
    for k in list(out.keys()):
        if out[k] == 0:
            out.pop(k, None)
    return out

def _contains_non_h_elements(parent_nonh: Counter, child_nonh: Counter) -> bool:
    """Return True if child has no new non-H elements and all counts <= parent."""
    for elem, cnt in child_nonh.items():
        if cnt <= 0:
            continue
        if parent_nonh.get(elem, 0) < cnt:
            return False
    return True

def _heavy_atom_indices(mol: Chem.Mol) -> List[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]

def _adjacency(mol: Chem.Mol) -> Dict[int, List[int]]:
    adj: Dict[int, List[int]] = {a.GetIdx(): [] for a in mol.GetAtoms()}
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)
    return adj

def _is_connected_subset(adj: Dict[int, List[int]], subset: Set[int]) -> bool:
    if not subset:
        return False
    start = next(iter(subset))
    seen = {start}
    q = deque([start])
    while q:
        v = q.popleft()
        for w in adj.get(v, []):
            if w in subset and w not in seen:
                seen.add(w)
                q.append(w)
    return len(seen) == len(subset)

def _enumerate_connected_subsets_of_size(
    mol: Chem.Mol,
    k: int,
    *,
    heavy_only: bool = True,
    limit: int = 100000,
) -> List[Tuple[int, ...]]:
    """
    Enumerate connected subsets of atom indices of size k (by default heavy atoms only).
    Frontier expansion avoids duplicates via a seed index ordering constraint.
    """
    if k <= 0:
        return []

    adj = _adjacency(mol)
    atoms = _heavy_atom_indices(mol) if heavy_only else list(range(mol.GetNumAtoms()))
    atom_set = set(atoms)

    results: Set[Tuple[int, ...]] = set()

    def expand(seed: int) -> None:
        stack = [(set([seed]), set(adj[seed]))]  # (chosen, frontier)
        while stack:
            chosen, frontier = stack.pop()

            if len(chosen) == k:
                results.add(tuple(sorted(chosen)))
                if len(results) >= limit:
                    return
                continue

            for nxt in list(frontier):
                if nxt not in atom_set:
                    continue
                if nxt in chosen:
                    continue
                if nxt < seed:
                    continue

                chosen2 = set(chosen)
                chosen2.add(nxt)

                frontier2 = set(frontier)
                frontier2.remove(nxt)
                frontier2.update(adj[nxt])
                frontier2.difference_update(chosen2)

                stack.append((chosen2, frontier2))
                if len(results) >= limit:
                    return

    for seed in atoms:
        expand(seed)
        if len(results) >= limit:
            break

    return sorted(results)

def _submol_from_atoms(mol: Chem.Mol, atom_indices: List[int]) -> Optional[Chem.Mol]:
    """
    Build submolecule containing only atoms in atom_indices and bonds between them.
    Open valences are left open; AddHs() later caps them.
    """
    atom_set = set(atom_indices)

    em = Chem.EditableMol(Chem.Mol())
    old_to_new: Dict[int, int] = {}

    for old_idx in atom_indices:
        a = mol.GetAtomWithIdx(old_idx)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetIsAromatic(a.GetIsAromatic())
        new_idx = em.AddAtom(na)
        old_to_new[old_idx] = new_idx

    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        if i in atom_set and j in atom_set:
            em.AddBond(old_to_new[i], old_to_new[j], b.GetBondType())

    sub = em.GetMol()
    try:
        Chem.SanitizeMol(sub)
    except Exception:
        return None
    return sub

def can_carve_formula_from_smiles(
    smiles: str,
    formula: str,
    *,
    require_connected: bool = True,
    max_matches_to_return: int = 10,
    enumeration_limit: int = 200000,
    h_tolerance_abs: int = 2,
    h_tolerance_fraction_of_parent: float = 0.20,
) -> Dict[str, Any]:
    """
    Graph-realizability check with relaxed hydrogen matching.

    - Non-H element counts must match the target exactly.
    - H is allowed to differ by +/- max(h_tolerance_abs, ceil(h_tolerance_fraction_of_parent * H_parent)),
      where H_parent is computed from the whole SMILES (after AddHs()).
    """
    target = formula_to_counter(formula)
    if not target:
        return {"ok": False, "reason": "Empty/invalid target formula", "fits": False}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"ok": False, "reason": "Invalid SMILES", "fits": False}

    parent_with_h = Chem.AddHs(mol)
    parent_formula = rdMolDescriptors.CalcMolFormula(parent_with_h)
    parent_counter = formula_to_counter(parent_formula)
    h_parent = int(parent_counter.get("H", 0))

    h_tol = max(int(h_tolerance_abs), int(math.ceil(h_tolerance_fraction_of_parent * h_parent)))

    target_non_h = _non_h_counter(target)
    target_h = int(target.get("H", 0))
    target_heavy = sum(v for k, v in target.items() if k != "H")

    if target_heavy <= 0:
        return {"ok": False, "reason": "Target formula has no heavy atoms", "fits": False}

    mol_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    if target_heavy > mol_heavy:
        return {
            "ok": True,
            "fits": False,
            "reason": "Target heavy atom count exceeds molecule heavy atom count",
            "target_formula": counter_to_formula(target),
            "parent_formula": parent_formula,
            "H_tolerance": h_tol,
            "n_matches": 0,
            "matches": [],
        }

    subsets = _enumerate_connected_subsets_of_size(
        mol, target_heavy, heavy_only=True, limit=enumeration_limit
    )

    adj = _adjacency(mol)
    matches: List[Dict[str, Any]] = []
    n_found = 0

    for subset in subsets:
        atom_idxs = list(subset)

        if require_connected and not _is_connected_subset(adj, set(atom_idxs)):
            continue

        sub = _submol_from_atoms(mol, atom_idxs)
        if sub is None:
            continue

        sub_h = Chem.AddHs(sub)
        frag_formula = rdMolDescriptors.CalcMolFormula(sub_h)
        frag_counter = formula_to_counter(frag_formula)

        frag_non_h = _non_h_counter(frag_counter)
        frag_h = int(frag_counter.get("H", 0))

        if frag_non_h != target_non_h:
            continue
        if abs(frag_h - target_h) > h_tol:
            continue

        n_found += 1
        if len(matches) < max_matches_to_return:
            frag_smiles = Chem.MolToSmiles(Chem.RemoveHs(sub_h), isomericSmiles=True)
            matches.append({
                "atom_indices": atom_idxs,
                "fragment_smiles": frag_smiles,
                "fragment_formula_capped": frag_formula,
                "fragment_H": frag_h,
                "target_H": target_h,
                "H_tolerance": h_tol,
                "delta_H": frag_h - target_h,
            })

    return {
        "ok": True,
        "fits": n_found > 0,
        "target_formula": counter_to_formula(target),
        "target_nonH": counter_to_formula(target_non_h),
        "target_H": target_h,
        "target_heavy_atoms": target_heavy,
        "parent_formula": parent_formula,
        "parent_H": h_parent,
        "H_tolerance": h_tol,
        "n_candidates_enumerated": len(subsets),
        "n_matches": n_found,
        "matches": matches,
    }

def smiles_can_fragment_to_formula(smiles: str, formula: str) -> bool:
    """Convenience boolean wrapper."""
    ret = can_carve_formula_from_smiles(smiles, formula)
    return bool(ret.get("fits"))