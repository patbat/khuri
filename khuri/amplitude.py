"""
============================================================
General two particle into two particle scattering amplitudes
============================================================

It is assumed that all particles have the same mass.
"""
import numpy as np

from khuri.phase_space import rho


def from_cot(mass):
    """Obtain amplitude from cotangent of the phase."""
    def wrapper(cot_phase):
        def amplitude(mandelstam_s):
            return 1 / (cot_phase(mandelstam_s) - 1j) / rho(mass, mandelstam_s)
        return amplitude
    return wrapper


def from_phase(mass, inelasticity=None):
    """Obtain amplitude from the phase."""
    if inelasticity is None:
        inelasticity = lambda x: 1.0
    def wrapper(phase):
        def amplitude(mandelstam_s):
            phase_value = phase(mandelstam_s)
            inelastic = inelasticity(mandelstam_s)
            numerator = inelastic * np.exp(2j * phase_value) - 1
            return numerator / rho(mass, mandelstam_s) / 2j
        return amplitude
    return wrapper


def second_sheet(mass):
    """Obtain the amplitude on the second sheet from the one on the first."""
    def wrapper(amplitude):
        def second(mandelstam_s):
            first = amplitude(mandelstam_s)
            return first / (1 + 2j * rho(mass, mandelstam_s) * first)

        return second
    return wrapper


def first_sheet(mass):
    """Obtain the amplitude on the first sheet from the one on the second."""
    def wrapper(amplitude):
        def first(mandelstam_s):
            second = amplitude(mandelstam_s)
            return second / (1 - 2j * rho(mass, mandelstam_s) * second)

        return first
    return wrapper
