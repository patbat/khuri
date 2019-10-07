"""
=======================================================
Parametrization of pi+pi -> pi+pi scattering amplitudes
=======================================================

Parametrizations of the pi+pi -> pi+pi scattering amplitudes and phases of the
s0-wave and the p-wave from

    <arXiv:1907.13162 by J. R. Peláez, A. Rodas, J. Ruiz de Elvira>.

Along the real axis the parametrizations are valid up to 2 GeV. The notation
from the aforementioned paper is used in the implementation. The values for
the non-dimensionless fit parameters/default arguments are given in GeV, so
Mandelstam s should accordingly be given in units of (GeV)^2.
"""
import numpy as np
import numpy.polynomial.chebyshev as cheb # use cheb.chebval or cheb.Chebyshev
from numpy.polynomial.polynomial import Polynomial
from scipy.misc import derivative

from khuri.phase_space import rho



## general / helpers ##########################################################


PION_MASS = 0.13957 # in GeV
S_MATCHING = 1.4**2 # in GeV


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


def chebyshev_derivative(number, x):
    """Evaluate the derivative of `number`-th Chebyshev polynomial at `x`."""
    coefficients = np.zeros(number + 1)
    coefficients[number] = 1
    return cheb.chebval(x, cheb.chebder(coefficients))


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


def omega_2(s, s_m=S_MATCHING, constant=2.0):
    """Variable for Chebyshev polynomials used in the high-energy region."""
    sqrt_m = np.sqrt(s_m)
    return 2.0 * (np.sqrt(s) - sqrt_m) / (constant - sqrt_m) - 1.0



## p-wave #####################################################################


P_WAVE_RHO_MASS = 0.7752 # in GeV

P_WAVE_CONFORMAL_COEFFICIENTS_1 = (0.97, 0.12, -0.18, 0.40, 1.65)
P_WAVE_D_0_1 = 11.1
P_WAVE_D_1_1 = 5.6
P_WAVE_K_0_1 = 0.35
P_WAVE_E_0_1 = -0.146
P_WAVE_E_1_1 = 0.337
P_WAVE_E_2_1 = 0.198

P_WAVE_CONFORMAL_COEFFICIENTS_2 = (0.97, 0.11, -0.13, 0.47, 1.36)
P_WAVE_D_0_2 = 3.6
P_WAVE_D_1_2 = 2.4
P_WAVE_K_0_2 = 0.17
P_WAVE_E_0_2 = -0.16
P_WAVE_E_1_2 = 0.043
P_WAVE_E_2_2 = 0.06


def p_wave_generate_cot_phase(conformal_coeff, rho_mass=P_WAVE_RHO_MASS):
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


def p_wave_extend_phase(d_0, d_1, s_m=S_MATCHING):
    """Extend phase to be valid up to 2 GeV."""
    def wrapper(low_energy_phase):
        value = low_energy_phase(s_m)
        d_value = derivative(low_energy_phase, s_m, dx=1e-6)

        delta = ((d_value - d_0 * chebyshev_derivative(2, -1.0)
                  - d_1 * chebyshev_derivative(3, -1.0))
                 / chebyshev_derivative(1, -1.0))

        def high_energy_phase(s):
            chebyshev = cheb.chebval(omega_2(s, s_m=s_m), [0, delta, d_0, d_1])
            return value + chebyshev + delta - d_0 - d_1

        def phase(s):
            s = np.asarray(s)
            return np.piecewise(s,
                                [s <= s_m, s > s_m],
                                [low_energy_phase, high_energy_phase])

        return phase

    return wrapper


def p_wave_generate_phase(conformal_coeff, d_0, d_1,
                          rho_mass=P_WAVE_RHO_MASS,
                          s_m=S_MATCHING):
    cot_phase = p_wave_generate_cot_phase(conformal_coeff, rho_mass=rho_mass)

    def low_energy_phase(s):
        return arccot2(cot_phase(s))

    return p_wave_extend_phase(d_0, d_1, s_m=s_m)(low_energy_phase)


p_wave_phase = p_wave_generate_phase(
                    conformal_coeff=P_WAVE_CONFORMAL_COEFFICIENTS_2,
                    d_0=P_WAVE_D_0_2, d_1=P_WAVE_D_1_2)


def p_wave_inelasticity(s, k_0, s_e=1.12**2):
    s = np.asarray(s)

    def above(x):
        return 1.0 - k_0 * (1.0 - x / s_e)**2

    return np.piecewise(s, [s < s_e, s >= s_e], [1, above])


