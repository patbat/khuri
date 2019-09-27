"""
===========================================
Facilites that are used in different tests.
===========================================
"""
import numpy as np


def schwarz(func, s, atol=1e-14, rtol=0):
    """Check whether `func` fulfills Schwarz reflection principle."""
    a = func(s)
    b = np.conjugate(func(np.conjugate(s)))
    assert np.allclose(a, b, atol=atol, rtol=rtol)
