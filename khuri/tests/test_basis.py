import numpy as np
import pytest

from khuri import madrid_global, phases, omnes
import khuri.khuri_treiman as kt


THRESHOLD = (2.0 * madrid_global.PION_MASS)**2


@phases.asymptotic1(matching_point=1.12**2)
def phase(s):
    return madrid_global.p_wave_phase(s)


def amplitude(s):
    return madrid_global.p_wave(s)


@pytest.fixture
def omnes_function():
    return omnes.generate_omnes(phase, threshold=THRESHOLD,
                                constant=np.pi, cut=1e10)


@pytest.fixture
def curve():
    return kt.Real(4.0, 100.0)


@pytest.fixture
def grid(curve):
    """Small grid (way too small for real world application)."""
    return kt.GridReal(curve, (5,), 2)


def test_basis(omnes_function, grid):
    """Test if simple construction and evalution work."""
    subtractions = 1
    pion_mass = 1.0
    virtuality = 0.0
    basis = kt.BasisReal(omnes_function,
                         amplitude,
                         subtractions,
                         grid,
                         pion_mass,
                         virtuality)
    assert isinstance(basis(0, 2.0-10.0j), complex)

    with pytest.raises(IndexError):
        basis(1, 10.0)
