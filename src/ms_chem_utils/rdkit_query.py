from __future__ import annotations

from typing import Dict, Optional, Set, Tuple
from rdkit import Chem

def make_generic_bond_query(child: Chem.Mol) -> Chem.Mol:
    """
    Return a query molecule derived from `child` where bond types are made generic
    (single/double/triple/aromatic all match). This avoids using Bond.SetQuery().
    """
    params = Chem.AdjustQueryParameters.NoAdjustments()
    params.makeBondsGeneric = True
    q = Chem.AdjustQueryProperties(child, params)
    return q

def unsat_units_from_mol(mol: Chem.Mol, atom_set: Optional[Set[int]] = None) -> int:
    """
    Unsaturation units based on heavy-heavy bonds:
      single = 0
      aromatic or double = 1
      triple = 2
    If atom_set is provided, only count bonds where both endpoints are in atom_set.
    """
    u = 0
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        if atom_set is not None and (i not in atom_set or j not in atom_set):
            continue
        ai = mol.GetAtomWithIdx(i)
        aj = mol.GetAtomWithIdx(j)
        if ai.GetAtomicNum() <= 1 or aj.GetAtomicNum() <= 1:
            continue

        if b.GetIsAromatic():
            u += 1
        else:
            bt = b.GetBondType()
            if bt == Chem.BondType.DOUBLE:
                u += 1
            elif bt == Chem.BondType.TRIPLE:
                u += 2
    return u




def heavy_heavy_bond_count(mol: Chem.Mol, atom_set=None) -> int:
    """Count heavy-heavy bonds fully inside atom_set."""
    
    if atom_set is None:
        atom_set = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    n = 0
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        if i in atom_set and j in atom_set:
            ai = mol.GetAtomWithIdx(i)
            aj = mol.GetAtomWithIdx(j)
            if ai.GetAtomicNum() > 1 and aj.GetAtomicNum() > 1:
                n += 1
    return n

def submol_from_match(parent: Chem.Mol, match: Tuple[int, ...]) -> Optional[Chem.Mol]:
    """
    Build a submolecule containing the matched atoms and internal bonds.
    Produces open valences where bonds to outside atoms were cut.
    """
    atom_indices = list(match)
    atom_set = set(atom_indices)

    em = Chem.EditableMol(Chem.Mol())
    old_to_new: Dict[int, int] = {}

    for old_idx in atom_indices:
        a = parent.GetAtomWithIdx(old_idx)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetIsAromatic(a.GetIsAromatic())
        new_idx = em.AddAtom(na)
        old_to_new[old_idx] = new_idx

    for b in parent.GetBonds():
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