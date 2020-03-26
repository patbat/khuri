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
from khuri.amplitude import from_cot


def tan_phase(mandelstam_s, isospin, pion_mass, peak, coefficients):
    """Return the tangent of the phase.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    peak: float
        the value of Mandelstam s at which the phase equals 90 degree
    coefficients: iterable
        the (dimensionless) coefficients of the polynomial in the
        parametrization
    """
    threshold = 4.0 * pion_mass**2
    param = mandelstam_s / threshold - 1.0
    phase_space = np.sqrt(1.0 - threshold / mandelstam_s)
    polynomial = sum(coeff * param**i for i, coeff in enumerate(coefficients))
    peak_factor = (threshold - peak) / (mandelstam_s - peak)
    return phase_space * param**isospin * polynomial * peak_factor


def partial_wave(mandelstam_s, isospin, pion_mass, peak, coefficients):
    """Return the partial wave amplitude.

    Parameters
    ----------
    mandelstam_s: float or array_like
        the value of the Mandelstam variable s, at which the phase is evaluated
    isospin: int
        the value of the isospin of the partial wave (one of 0, 1, 2)
    pion_mass: float
        the value of the pion mass
    peak: float
        the value of Mandelstam s at which the phase equals 90 degree
    coefficients: iterable
        the (dimensionless) coefficients of the polynomial in the
        parametrization
    """
    @from_cot(pion_mass)
    def amp(s_value):
        return 1.0 / tan_phase(s_value, isospin, pion_mass, peak, coefficients)
    return amp(mandelstam_s)


def partial_wave_lit(mandelstam_s, isospin, pion_mass):
    params = literature_values(isospin, pion_mass)
    return partial_wave(mandelstam_s, isospin, pion_mass, *params)


def literature_values(isospin, pion_mass):
    if isospin == 0:
        return 36.77 * pion_mass**2, (0.22, 0.268, -0.139e-1, -0.139e-2)
    if isospin == 1:
        return 30.72 * pion_mass**2, (0.379e-1, 0.14e-4, -0.673e-4, 0.163e-7)
    if isospin == 2:
        return -21.62 * pion_mass**2, (-0.444e-1, -0.857e-1, -0.221e-2,
                                       -0.129e-3)
