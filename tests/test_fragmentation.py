# tests/test_fragmentation.py
from tests.conftest import requires_rdkit

@requires_rdkit
def test_fragmentation_smoke():
    from ms_chem_utils.fragmentation import (
        make_generic_bond_query,
        submol_from_match,
    )
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCOC(=O)C")
    assert mol is not None

    q = make_generic_bond_query(mol)
    assert q is not None

    # simplest match: the molecule against itself should match
    matches = mol.GetSubstructMatches(mol)
    assert len(matches) >= 1
    sub = submol_from_match(mol, matches[0])
    assert sub is not None

@requires_rdkit
def test_fragmentation_does_not_throw_on_small_molecule():
    from ms_chem_utils.fragmentation import heavy_heavy_bond_count
    from rdkit import Chem

    mol = Chem.MolFromSmiles("O")
    assert mol is not None
    n = heavy_heavy_bond_count(mol)
    assert isinstance(n, int)
    assert n >= 0