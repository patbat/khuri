"""
=======================================================
Parametrization of pi+pi -> pi+pi scattering amplitudes
=======================================================

Parametrizations of the pi+pi -> pi+pi scattering amplitudes and phases of the
s0-wave and the p-wave from

    <arXiv:1907.13162 by J. R. Peláez, A. Rodas, J. Ruiz de Elvira>.

Along the real axis the parametrizations are valid up to 2 GeV. The notation
from the aforementioned paper is used in the implementation.
"""
import numpy as np
import numpy.polynomial.chebyshev as cheb # use cheb.chebval or cheb.Chebyshev


def omega(s, s0, alpha):
    """A conformal variable.

    Parameters
    ----------
    s: the Mandelstam variable s
    s0: usually non-negative. Above s0, the conformal variable is complex.
    alpha: the center of the conformal expansion
    """
    sqrt_s = np.sqrt(s)
    diff = alpha * np.sqrt(s_0 - s)
    return (sqrt_s - diff) / (sqrt_s + diff)


def omega_2(s, s_m, constant=2.0):
    """Variable for Chebyshev polynomials used in the high-energy region."""
    sqrt_m = np.sqrt(s_m)
    return 2.0 * (np.sqrt(s) - sqrt_m) / (constant - sqrt_m) - 1.0
