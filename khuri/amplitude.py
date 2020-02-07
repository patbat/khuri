"""
============================================================
General two particle into two particle scattering amplitudes
============================================================

It is assumed that all particles have the same mass.
"""
from functools import wraps

import numpy as np

from khuri.phase_space import rho


def from_cot(mass):
    """Obtain amplitude from cotangent of the phase."""
    def wrapper(cot_phase):
        @wraps(cot_phase)
        def amplitude(*args):
            return 1 / (cot_phase(*args) - 1j) / rho(mass, args[-1])
        return amplitude
    return wrapper


def from_phase(mass, inelasticity=None):
    """Obtain amplitude from the phase."""
    if inelasticity is None:
        inelasticity = lambda _: 1.0
    def wrapper(phase):
        @wraps(phase)
        def amplitude(*args):
            phase_value = phase(*args)
            inelastic = inelasticity(args[-1])
            numerator = inelastic * np.exp(2j * phase_value) - 1
            return numerator / rho(mass, args[-1]) / 2j
        return amplitude
    return wrapper


def second_sheet(mass):
    """Obtain the amplitude on the second sheet from the one on the first."""
    def wrapper(amplitude):
        @wraps(amplitude)
        def second(*args):
            first = amplitude(*args)
            return first / (1 + 2j * rho(mass, args[-1]) * first)

        return second
    return wrapper


def first_sheet(mass):
    """Obtain the amplitude on the first sheet from the one on the second."""
    def wrapper(amplitude):
        @wraps(amplitude)
        def first(*args):
            second = amplitude(*args)
            return second / (1 - 2j * rho(mass, args[-1]) * second)

        return first
    return wrapper
