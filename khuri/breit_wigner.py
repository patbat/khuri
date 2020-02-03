"""
================================
Breit Wigner amplitude and phase
================================
"""
import numpy as np


def amplitude(mandelstam_s, mass, width, coupling1, coupling2=None):
    """The Breit Wigner amplitude."""
    if coupling2 is None:
        coupling2 = coupling1
    denominator = mandelstam_s - mass**2 + np.sqrt(mandelstam_s) * 1.0j * width
    return - coupling1 * coupling2 / denominator
