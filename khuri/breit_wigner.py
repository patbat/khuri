"""
================================
Breit Wigner amplitude and phase
================================

The expressions are mainly taken from
<D. M. Asner, C. Hanhart, and E. Klempt. “Resonances”.
 In: Review of Particle Physics.  2015. Chap. 48.>.
"""
import numpy as np


def simple_amplitude(mandelstam_s, mass, width, coupling1, coupling2=None):
    """The Breit Wigner amplitude in its simple form."""
    if coupling2 is None:
        coupling2 = coupling1
    denominator = mandelstam_s - mass**2 + np.sqrt(mandelstam_s) * 1.0j * width
    return - coupling1 * coupling2 / denominator


def amplitude(mandelstam_s,
              resonance_mass,
              resonance_width,
              angular_momentum,
              mass,
              coupling1,
              coupling2=None):
    """The Breit Wigner amplitude with and energy dependent width."""
    width = width_energy_dependent(mandelstam_s, resonance_mass,
                                   resonance_width, angular_momentum, mass)
    return simple_amplitude(mandelstam_s, resonance_mass, width, coupling1,
                            coupling2)


def width_energy_dependent(mandelstam_s,
                           resonance_mass,
                           resonance_width,
                           angular_momentum,
                           mass):
    """Energy dependent width."""
    momentum1 = cm_momentum(mandelstam_s, mass)
    momentum2 = cm_momentum(resonance_mass**2, mass)
    momentum = (momentum1 / momentum2) ** (2 * angular_momentum + 1)

    blatt_weiss1 = blatt_weisskopf(angular_momentum, momentum1)
    blatt_weiss2 = blatt_weisskopf(angular_momentum, momentum2)
    blatt_weiss = (blatt_weiss1 / blatt_weiss2)**2

    energy = resonance_mass / np.sqrt(mandelstam_s)

    return resonance_width * momentum * energy * blatt_weiss


def cm_momentum(mandelstam_s, mass):
    """Center of mass momentum."""
    return np.sqrt(mandelstam_s / 4.0 - mass**2)


def blatt_weisskopf(angular_momentum: int, momentum):
    """The Blatt-Weisskopf barrier factors SQUARED.

    Parameters
    ----------
    angular_momentum: the angular momentum quantum number.
    momentum: the magnitude of the momentum, usually expressed in terms of an
        energy scale representing the range of interaction.

    Notes
    -----
    Taken from <Ann. Physik 4 (1995) 404 - 430>.
    """
    if angular_momentum == 0:
        return _zero(momentum)
    if angular_momentum == 1:
        return _one(momentum)
    if angular_momentum == 2:
        return _two(momentum)
    if angular_momentum == 3:
        return _three(momentum)
    if angular_momentum == 4:
        return _four(momentum)
    raise NotImplementedError('Blatt Weisskopf for J > 4 not implemented.')


def _zero(z):
    return 1.0


def _one(z):
    return 2.0 * z / (z + 1.0)


def _two(z):
    return 13.0 * z**2 / ((z - 3.0)**2 + 9.0 * z)


def _three(z):
    return 277.0 * z**3 / (z * (z - 15.0)**2 + 9.0 * (2.0 * z - 5.0))


def _four(z):
    return 12746.0 * z**4 / ((z**2 - 45.0 * z + 105.0)**2
                             + 25.0 * z * (2.0 * z - 21.0)**2)
