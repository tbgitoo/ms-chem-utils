from ms_chem_utils import smiles_fits_into_smiles

def test_smiles_fits_basic():
    # ethanol fits into ethyl acetate? (as a substructure, generic bonds allow some flexibility)
    out = smiles_fits_into_smiles("CCO", "CCOC(=O)C")
    assert out["ok"] is True
    assert "fits" in out