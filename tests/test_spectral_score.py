# tests/test_spectral_score.py
def test_score_monotonicity_smoke():
    """
    Whatever your scoring function is, it should generally score
    a closer match >= a worse match (monotonicity sanity).
    """
    from ms_chem_utils.spectral_score import score_mass_error_ppm

    s_close = score_mass_error_ppm(1.0)   # 1 ppm
    s_far   = score_mass_error_ppm(50.0)  # 50 ppm

    assert isinstance(s_close, (int, float))
    assert isinstance(s_far, (int, float))
    assert s_close >= s_far

def test_score_bounds():
    from ms_chem_utils.spectral_score import score_mass_error_ppm

    s0 = score_mass_error_ppm(0.0)
    s1 = score_mass_error_ppm(1e6)
    assert s0 >= s1