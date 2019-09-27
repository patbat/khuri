import numpy as np
import pytest

from khuri.omnes import Omnes
from khuri import madrid
from khuri.tests.test_phases import phase


# A small value is added to avoid the evaluation of the parametrization by
# the Madrid group at threshold
THRESHOLD = (2.0 * madrid.PION_MASS)**2 + 1e-6


@pytest.fixture
def omnes_function():
    return Omnes(phase, threshold=THRESHOLD)


def test_omnes_at_zero(omnes_function):
    assert omnes_function(0) == pytest.approx(1.0)


def test_phase(omnes_function):
    """Check if the phase of the Omnes function agrees with the input phase."""
    mandelstam_s = np.linspace(THRESHOLD, 1400, 20)**2
    omnes_phase = np.angle(omnes_function(mandelstam_s))
    phase_values = phase(mandelstam_s)
    assert np.allclose(omnes_phase, phase_values)
