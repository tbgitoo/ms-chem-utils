# tests/test_containment.py
from tests.conftest import requires_rdkit

@requires_rdkit
def test_smiles_fits_basic_true():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    out = smiles_fits_into_smiles("CCO", "CCOC(=O)C")
    assert isinstance(out, dict)
    assert out.get("ok") is True
    assert out.get("fits") is True
    assert "reason" in out
    assert "example" in out

@requires_rdkit
def test_smiles_fits_basic_false():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    # child cannot fit into parent if parent is smaller / unrelated
    out = smiles_fits_into_smiles("CCOC(=O)C", "CCO")
    assert out.get("ok") is True
    assert out.get("fits") is False
    assert "reason" in out

@requires_rdkit
def test_smiles_fits_invalid_smiles_child():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    out = smiles_fits_into_smiles("NOT_A_SMILES", "CCOC(=O)C")
    assert out.get("ok") is False
    assert out.get("fits") is False
    assert "Invalid child SMILES" in out.get("reason", "")

@requires_rdkit
def test_smiles_fits_invalid_smiles_parent():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    out = smiles_fits_into_smiles("CCO", "NOT_A_SMILES")
    assert out.get("ok") is False
    assert out.get("fits") is False
    assert "Invalid parent SMILES" in out.get("reason", "")

@requires_rdkit
def test_smiles_fits_empty_strings():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    out = smiles_fits_into_smiles("", "")
    assert out.get("ok") is False
    assert out.get("fits") is False

@requires_rdkit
def test_smiles_fits_parameter_sanity():
    from ms_chem_utils.containment import smiles_fits_into_smiles

    out = smiles_fits_into_smiles(
        "CCO",
        "CCOC(=O)C",
        h_tolerance_abs=0,
        unsat_tolerance_abs=0,
        bond_tolerance_abs=0,
        max_substruct_matches_to_check=5,
    )
    assert out.get("ok") is True
    # fits may still be true/false depending on your internal logic,
    # but must not throw and must return consistent shape
    assert "fits" in out
    assert "reason" in out