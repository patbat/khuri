import itertools

import numpy as np
import pytest

from khuri.omnes import generate_omnes
from khuri.gsl import IntegrationRoutine
from khuri import madrid
from khuri.tests.helpers import schwarz
from khuri.tests.test_phases import PHASES, MATCHING_POINTS


# A small value is added to avoid the evaluation of the (non-global)
# parametrization by the Madrid group at threshold
THRESHOLD = (2.0 * madrid.PION_MASS)**2 + 1e-6


def all_omnes_for_phase(phase):
    return (
        generate_omnes(phase, threshold=THRESHOLD, constant=np.pi, cut=1e10),
        generate_omnes(phase, threshold=THRESHOLD),
        generate_omnes(phase, threshold=THRESHOLD, constant=np.pi, cut=5e3,
                       integration_routine=IntegrationRoutine.qag),
        generate_omnes(phase, threshold=THRESHOLD,
                       integration_routine=IntegrationRoutine.qag),
    )


OMNES_LIST = list(itertools.chain.from_iterable(all_omnes_for_phase(phase)
                                                for phase in PHASES))
PHASE_LIST = list(itertools.chain.from_iterable(itertools.repeat(phase, 4)
                                                for phase in PHASES))
OMNES_IDS = [
    'p_wave_cquad_cut',
    'p_wave_cquad',
    'p_wave_qag_cut',
    'p_wave_qag',
    'p_wave_global_cquad_cut',
    'p_wave_global_cquad',
    'p_wave_global_qag_cut',
    'p_wave_global_qag',
    's_wave_global_cquad_cut',
    's_wave_global_cquad',
    's_wave_global_qag_cut',
    's_wave_global_qag',
]


@pytest.mark.parametrize('function', OMNES_LIST, ids=OMNES_IDS)
def test_omnes_at_zero(function):
    assert function(0) == pytest.approx(1.0)


@pytest.mark.parametrize('function, phase', list(zip(OMNES_LIST, PHASE_LIST)),
                         ids=OMNES_IDS)
def test_phase(function, phase):
    """Check if the phase of the Omnes function agrees with the input phase."""
    # The phase of the s-wave exceeds pi close to kaon threshold, thus
    # for simplicity only the region below is checked.
    mandelstam_s = np.linspace(THRESHOLD, min(MATCHING_POINTS) - 0.02, 20)
    omnes_phase = np.angle(function(mandelstam_s))
    phase_values = phase(mandelstam_s)
    assert np.allclose(omnes_phase[1:], phase_values[1:])


@pytest.mark.parametrize('function', OMNES_LIST, ids=OMNES_IDS)
def test_schwarz(function):
    """Check if the Omnes function fulfills the Schwarz reflection princple."""
    real_part = 0.760**2
    imaginary_parts = np.linspace(-1e4, 1e4, 20)
    mandelstam_s = real_part + 1j * imaginary_parts
    schwarz(function, mandelstam_s)
