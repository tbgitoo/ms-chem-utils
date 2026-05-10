from ms_chem_utils import formula_to_counter, counter_to_formula

def test_roundtrip_formula():
    c = formula_to_counter("C2H6O")
    assert c["C"] == 2 and c["H"] == 6 and c["O"] == 1
    assert counter_to_formula(c) in ("C2H6O",)  # formatting canonical