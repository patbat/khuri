import numpy as np
import pytest

from khuri.omnes import Omnes
from khuri import madrid
from khuri.tests.helpers import schwarz
from khuri.tests.test_phases import phase


# A small value is added to avoid the evaluation of the parametrization by
# the Madrid group at threshold
THRESHOLD = (2.0 * madrid.PION_MASS)**2 + 1e-6


@pytest.fixture
def omnes_function():
    return Omnes(phase, threshold=THRESHOLD, constant=np.pi, cut=1e10)


@pytest.fixture
def omnes_function_alt():
    return Omnes(phase, threshold=THRESHOLD)


def test_omnes_at_zero(omnes_function, omnes_function_alt):
    assert omnes_function(0) == pytest.approx(1.0)
    assert omnes_function_alt(0) == pytest.approx(1.0)


def test_phase(omnes_function):
    """Check if the phase of the Omnes function agrees with the input phase."""
    mandelstam_s = np.linspace(THRESHOLD, 1400, 20)**2
    omnes_phase = np.angle(omnes_function(mandelstam_s))
    phase_values = phase(mandelstam_s)
    assert np.allclose(omnes_phase, phase_values)


def test_schwarz(omnes_function):
    """Check if the Omnes function fulfills the Schwarz reflection princple."""
    real_part = 760**2
    imaginary_parts = np.linspace(-1e4, 1e4, 20)
    mandelstam_s = real_part + 1j * imaginary_parts
    schwarz(omnes_function, mandelstam_s)
