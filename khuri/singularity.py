"""
========================================
Singularities in Khuri Treiman equations
========================================
"""
from collections.abc import Iterable

import numpy as np
from numpy.lib.scimath import sqrt  # handles negative reals correctly

CUBE_ROOT = (-1.0 + sqrt(3) * 1j) / 2.0


def cube_root(value):
    if isinstance(value, Iterable):
        return np.array(tuple(float(x)**(1/3) for x in value))
    return float(value)**(1/3)


def cubic_equation(coeff, number):
    c0, c1, c2, c3 = coeff
    delta0 = c2**2 - 3.0 * c3 * c1
    delta1 = 2.0 * c2**3 - 9.0 * c3 * c2 * c1 + 27.0 * c3**2 * c0
    sqr = sqrt(delta1**2 - 4.0 * delta0**3)
    delta = cube_root((delta1 + sqr) / 2.0)
    if abs(delta) < 1e-20:
        delta = cube_root((delta1 - sqr) / 2.0)
    delta *= CUBE_ROOT**number
    return -1.0 * (c2 + delta + delta0 / delta) / 3.0 / c3


def coefficient_3(cos2):
    return cos2 - 1


def coefficient_2(cos2, decay_mass_2, mass_2):
    return 2.0 * ((decay_mass_2 - 5.0 * mass_2)
                  - cos2 * (decay_mass_2 + 3.0 * mass_2))


def coefficient_1(cos2, decay_mass_2, mass_2):
    return ((cos2 * (decay_mass_2 + 3.0 * mass_2))**2
            - (decay_mass_2 - 5.0 * mass_2)**2)


def coefficient_0(cos2, decay_mass_2, mass_2):
    return - (2.0 * cos2 * mass_2 * (decay_mass_2 - mass_2))**2


def coefficients(cos2, decay_mass_2, mass_2):
    return (coefficient_0(cos2, decay_mass_2, mass_2),
            coefficient_1(cos2, decay_mass_2, mass_2),
            coefficient_2(cos2, decay_mass_2, mass_2),
            coefficient_3(cos2))


def singularity(cos2, decay_mass_2, mass_2, number):
    coeff = coefficients(cos2, decay_mass_2, mass_2)
    return cubic_equation(coeff, number)
