"""
=======================================================
Parametrization of pi+pi -> pi+pi scattering amplitudes
=======================================================

Parametrizations of the pi+pi -> pi+pi scattering amplitudes and phases of the
s0-wave and the p-wave from

    <arXiv:1907.13162 by J. R. Peláez, A. Rodas, J. Ruiz de Elvira>.

Along the real axis the parametrizations are valid up to 2 GeV. The notation
from the aforementioned paper is used in the implementation. The values for
the non-dimensionless fit parameters are given in GeV, so Mandelstam s should
accordingly be given in units of (GeV)^2.
"""
import numpy as np
import numpy.polynomial.chebyshev as cheb # use cheb.chebval or cheb.Chebyshev
from numpy.polynomial.polynomial import Polynomial

from khuri.phase_space import rho



## general ####################################################################


PION_MASS = 0.13957


def amplitude_from_cot(cot_phase):
    def amplitude(s):
        return 1 / (cot_phase(s) - 1j) / rho(PION_MASS, s)
    return amplitude


def amplitude_from_phase(phase):
    def amplitude(s):
        p = phase(s)
        return np.exp(1j * p) * np.sin(p) / rho(PION_MASS, s)
    return amplitude


def second_sheet_from_first(amplitude):
    def second(s):
        first = amplitude(s)
        return first / (1 + 2j * rho(PION_MASS, s) * first)

    return second


def arccot2(x):
    return np.arctan2(1.0, x)


def momentum(s):
    return np.sqrt(s) * rho(PION_MASS, s) / 2.0


def omega(s, s_0, alpha):
    """A conformal variable.

    Parameters
    ----------
    s: the Mandelstam variable s
    s_0: usually non-negative. Above s_0, the conformal variable is complex.
    alpha: the center of the conformal expansion
    """
    sqrt_s = np.sqrt(s)
    diff = alpha * np.sqrt(s_0 - s)
    return (sqrt_s - diff) / (sqrt_s + diff)


def omega_2(s, s_m, constant=2.0):
    """Variable for Chebyshev polynomials used in the high-energy region."""
    sqrt_m = np.sqrt(s_m)
    return 2.0 * (np.sqrt(s) - sqrt_m) / (constant - sqrt_m) - 1.0



## p-wave #####################################################################


P_WAVE_RHO_MASS = 0.7752
P_WAVE_CONFORMAL_COEFFICIENTS_1 = (0.97, 0.12, -0.18, 0.40, 1.65)
P_WAVE_CONFORMAL_COEFFICIENTS_2 = (0.97, 0.11, -0.13, 0.47, 1.36)


def p_wave_generate_cot_phase(rho_mass=P_WAVE_RHO_MASS,
                              conformal_coeff=P_WAVE_CONFORMAL_COEFFICIENTS_1):
    conformal_polynomial = Polynomial(conformal_coeff)
    rho2 = rho_mass**2
    pion = PION_MASS**3
    s_0 = 1.43**2

    def cot_phase(s):
        conf = conformal_polynomial(omega(s, s_0=s_0, alpha=0.3))
        sqrt_s = np.sqrt(s)
        prefactor = sqrt_s * (rho2 - s) / 2 / momentum(s)**3
        remove_spurious = 2 * pion / rho2 / sqrt_s
        return prefactor * (remove_spurious + conf)

    return cot_phase


def p_wave_inelasticity(s, k_0, s_e=1.12**2):
    s = np.asarray(s)

    def above(x):
        return 1.0 - k_0 * (1.0 - x / s_e)**2

    return np.piecewise(s, [s < s_e, s >= s_e], [1, above])



