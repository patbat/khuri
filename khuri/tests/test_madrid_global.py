import numpy as np
import pytest

from khuri import madrid_global as mg
from khuri.tests.helpers import schwarz, connected


THRESHOLD = 4.0 * mg.PION_MASS**2


@pytest.fixture
def vertical_line():
    """Vertical line in the complex s-plane."""
    real_part = 0.76**2
    imaginary_parts = np.linspace(-1e4, 1e4, 100)
    return real_part + 1j * imaginary_parts


@pytest.fixture
def real_axis():
    """Part of the real axis in the complex s-plane."""
    return np.linspace(THRESHOLD, 1.4**2, 100)


@pytest.mark.parametrize('function',
                         [mg.p_wave, mg.p_wave_2, mg.s_wave, mg.s_wave_2],
                         ids=['p_wave', 'p_wave_2', 's_wave', 's_wave_2'])
def test_schwarz(function, vertical_line):
    schwarz(function, vertical_line)


@pytest.mark.parametrize('first, second',
                         [(mg.p_wave, mg.p_wave_2), (mg.s_wave, mg.s_wave_2)],
                         ids=['p_wave', 's_wave'])
def test_connected(first, second, real_axis):
    connected(first, second, real_axis)
