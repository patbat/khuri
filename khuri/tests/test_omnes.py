import numpy as np
import pytest

from khuri.omnes import generate_omnes
from khuri.gsl import IntegrationRoutine
from khuri import madrid
from khuri.tests.helpers import schwarz
from khuri.tests.test_phases import phase


# A small value is added to avoid the evaluation of the parametrization by
# the Madrid group at threshold
THRESHOLD = (2.0 * madrid.PION_MASS)**2 + 1e-6


def get_omnes(**kwargs):
    return generate_omnes(phase, threshold=THRESHOLD, **kwargs)


OMNES_LIST = (
    get_omnes(constant=np.pi, cut=1e10),
    get_omnes(),
    get_omnes(integration_routine=IntegrationRoutine.qag),
    get_omnes(constant=np.pi, cut=1e4,
              integration_routine=IntegrationRoutine.qag),
)


@pytest.mark.parametrize('function', OMNES_LIST)
def test_omnes_at_zero(function):
    assert function(0) == pytest.approx(1.0)


@pytest.mark.parametrize('function', OMNES_LIST)
def test_phase(function):
    """Check if the phase of the Omnes function agrees with the input phase."""
    mandelstam_s = np.linspace(THRESHOLD, 1.4**2, 20)
    omnes_phase = np.angle(function(mandelstam_s))
    phase_values = phase(mandelstam_s)
    assert np.allclose(omnes_phase, phase_values)


@pytest.mark.parametrize('function', OMNES_LIST)
def test_schwarz(function):
    """Check if the Omnes function fulfills the Schwarz reflection princple."""
    real_part = 0.760**2
    imaginary_parts = np.linspace(-1e4, 1e4, 20)
    mandelstam_s = real_part + 1j * imaginary_parts
    schwarz(function, mandelstam_s)
