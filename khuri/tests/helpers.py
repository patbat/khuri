"""
===========================================
Facilites that are used in different tests.
===========================================
"""
import numpy as np


def schwarz(func, s, atol=1e-14, rtol=0):
    """Check whether `func` fulfills Schwarz reflection principle."""
    a = func(s)
    b = func(s.conjugate()).conjugate()
    assert np.allclose(a, b, atol=atol, rtol=rtol)


def connected(sheet1, sheet2, s, epsilon=1e-15, rtol=1e-5, atol=1e-8):
    """Check if `sheet1` and `sheet2` are continuously connected."""
    s_plus = s + epsilon * 1j
    s_minus = s_plus.conjugate()

    first_sheet_above = sheet1(s_plus)
    first_sheet_below = sheet1(s_minus)

    second_sheet_above = sheet2(s_plus)
    second_sheet_below = sheet2(s_minus)

    assert np.allclose(first_sheet_above, second_sheet_below, rtol=rtol,
                       atol=atol)
    assert np.allclose(first_sheet_below, second_sheet_above, rtol=rtol,
                       atol=atol)
