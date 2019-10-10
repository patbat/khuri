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


def test_p_wave_schwarz(vertical_line):
    schwarz(mg.p_wave, vertical_line)


def test_p_wave_second_sheet_schwarz(vertical_line):
    schwarz(mg.p_wave_2, vertical_line)


def test_p_wave_connected(real_axis):
    connected(mg.p_wave, mg.p_wave_2, real_axis)


def test_s_wave_schwarz(vertical_line):
    schwarz(mg.s_wave, vertical_line)


def test_s_wave_second_sheet_schwarz(vertical_line):
    schwarz(mg.s_wave_2, vertical_line)


def test_s_wave_connected(real_axis):
    connected(mg.s_wave, mg.s_wave_2, real_axis)
