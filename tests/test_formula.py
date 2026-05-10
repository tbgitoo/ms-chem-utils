# tests/test_formula.py
def test_formula_roundtrip_simple():
    from ms_chem_utils.formula import formula_to_counter, counter_to_formula

    c = formula_to_counter("C6H6O")
    assert c["C"] == 6
    assert c["H"] == 6
    assert c["O"] == 1
    f = counter_to_formula(c)
    # order may vary; enforce semantic equivalence
    c2 = formula_to_counter(f)
    assert c2 == c

def test_formula_handles_empty():
    from ms_chem_utils.formula import formula_to_counter, counter_to_formula

    c = formula_to_counter("")
    assert isinstance(c, dict)
    f = counter_to_formula(c)
    assert isinstance(f, str)

def test_formula_rejects_nonsense_gracefully():
    from ms_chem_utils.formula import formula_to_counter

    # your function might raise or return empty; pick your intended behavior:
    try:
        c = formula_to_counter("C6H6O+++")
        assert isinstance(c, dict)
    except Exception:
        # acceptable if you explicitly raise; if you prefer no-throw, remove this
        pass