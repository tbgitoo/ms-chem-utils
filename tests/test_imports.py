# tests/test_imports.py
import importlib

def test_import_package():
    import ms_chem_utils  # noqa: F401

def test_import_submodules():
    # these should import without side effects / NameErrors
    importlib.import_module("ms_chem_utils.containment")
    importlib.import_module("ms_chem_utils.formula")
    importlib.import_module("ms_chem_utils.fragmentation")
    importlib.import_module("ms_chem_utils.spectral_score")
    importlib.import_module("ms_chem_utils.spectra_pubchem")

def test_public_api_smoke():
    from ms_chem_utils.containment import smiles_fits_into_smiles  # noqa: F401