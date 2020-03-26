"""
============================================================
The Schenk parametrization pion+pion -> pion+pion scattering
============================================================

The original form of the parametrization is given in
<A. Schenk, Absorption and dispersion of pions at finite temperature.
Nucl. Phys. B363, 97 (1991)>.

More recent values for the parameters are given in
<G. Colangelo, J. Gasser, H.Leutwyler, Nuclear Physics B 603 (2001) 125–179>.
"""
import numpy as np


def tan_phase(mandelstam_s, angular_momentum, peak, coefficients, pion_mass):
    """Return the tangent of the phase.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    angular_momentum: int
        the value of the angular momentum of the partial wave
    peak: float
        the value of Mandelstam s at which the phase equals 90 degree
    coefficients: iterable
        the coefficients of the polynomial in the parametrization
    pion_mass: float
        the value of the pion mass
    """
    threshold = 4.0 * pion_mass**2
    param = mandelstam_s / threshold - 1.0
    phase_space = np.sqrt(1.0 - threshold / mandelstam_s)
    polynomial = sum(coeff * param**i for i, coeff in enumerate(coefficients))
    peak_factor = (threshold - peak) / (mandelstam_s - peak)
    return phase_space * param**angular_momentum * polynomial * peak_factor
