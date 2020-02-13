"""
=============================================================================
Poles/residues of partial waves for 2 to 2 scattering of equal mass particles
=============================================================================

Most important functions:
* pole_and_coupling
* pole
* coupling
"""
from typing import Callable, Tuple
import warnings

import numpy as np
from scipy.optimize import root
from scipy.misc import derivative

from khuri.phase_space import rho


# a scattering amplitude/partial wave on the first sheet
# input: mandelstam s
# output: value of the scattering amplitude/partial wave
Amplitude = Callable[[complex], complex]


def lower(func: Callable[..., complex]):
    """The decorated function returns only values in the lower half plane.

    This is useful e.g. when one searches for poles:
    sometimes the pole in the upper plane instead of the one in the lower
    half plane is detected. Since the functions of interest obey the
    Schwarz reflection principle, this will also alter the sign of the
    imaginary part of the residue.
    """
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        return result.real - 1j * abs(result.imag)
    return wrapper


class PoleError(Exception):
    """Raised when pole cannot be found."""

    def __init__(self, message):
        self.message = message


def pole_and_coupling(func: Amplitude,
                      mass: float,
                      guess: complex) -> Tuple[complex, complex]:
    """Determine pole and squared coupling given the amplitude on 1st sheet."""
    pole_position = pole(func, mass, guess)
    coupling_value = coupling(func, mass, pole_position, compute_pole=False)
    return pole_position, coupling_value


def pole(func: Amplitude,
         mass: float,
         guess: complex) -> complex:
    """Determine a pole on the 2nd sheet given the amplitude on the 1st.

    Note
    ----
    The pole position is calculated via finding the root of the denominator of
    the amplitude. In addition to the root associated with the rho
    resonance there is one at threshold, the latter one does not correspond to
    a pole, since the numerator of the IAM is zero at threshold, too. One needs
    to be careful in distinguishing both roots espiacially at high pion masses,
    when they get relatively close to each other.
    """
    result = _pole_from_first_sheet(func, mass, guess)
    if near_threshold(result, mass):
        warnings.warn("The determined pole is close to threshold, "
                      "might be artificial.")
    return result


@lower
def coupling(func: Amplitude,
             mass: float,
             pole_position: complex,
             compute_pole: bool = True) -> complex:
    """Return coupling squared on the 2nd sheet given the amplitude on the 1st.

    Parameters
    ----------
    compute_pole: if True, the pole is computed anew using `pole_position` as
        an initial guess. If False, the pole is set to `pole_position`.

    Note
    ----
    This assumes a p-wave.
    """
    if compute_pole:
        pole_position = pole(func, mass, pole_position)
    denominator = generate_denominator(func, mass)
    numerator = func(pole_position)
    residue = numerator / derivative(denominator, pole_position, dx=1e-3)
    coupling_squared = -48.0 * np.pi / (pole_position - 4.0 * mass**2) * residue
    return coupling_squared


def near_threshold(poles, mass, radius=1.0):
    """Check if any pole lies inside a disk around threshold with given radius.

    Parameters
    ----------
    poles: number or iterable of numbers
    mass: number
    """
    threshold = 4.0 * mass**2
    distances = np.abs(np.asarray(poles) - threshold) / threshold
    return np.any(distances < radius)


def generate_denominator(func: Amplitude,
                         mass: float) -> Callable[[complex], complex]:
    """Generate the denominator of the amplitude on the 2nd sheet."""
    def wrapper(mandelstam_s):
        return 1.0 + 2.0j * rho(mass, mandelstam_s) * func(mandelstam_s)
    return wrapper


@lower
def _pole_from_first_sheet(func: Amplitude,
                           mass: float,
                           guess: complex) -> complex:
    denominator = generate_denominator(func, mass)
    res = root(lambda x: complex_to_array(denominator(complex(*x))),
               [guess.real, guess.imag], method='lm')
    result = complex(*res.x)
    # Sometimes the root finding routine signals success although the
    # denominator does not vanish at the proposed solution, hence the second
    # condition.
    if not res.success or abs(denominator(result) > 1e-10):
        raise PoleError('pole: could not find root')
    return result


def complex_to_array(number):
    """Convert a complex number to an array of two numbers."""
    return np.real(number), np.imag(number)
