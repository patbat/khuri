import pytest
import numpy as np

from khuri.madrid import phase, PION_MASS
import khuri.amplitude as ka


CUT = np.linspace(2.0 * PION_MASS + 1, 4.0 * PION_MASS, 2)**2


@ka.from_phase(PION_MASS)
def amplitude(mandelstam_s):
    return phase(mandelstam_s)

@ka.from_cot(PION_MASS)
def amplitude2(mandelstam_s):
    return 1.0 / np.tan(phase(mandelstam_s))


def test_consistency():
    assert np.allclose(amplitude(CUT), amplitude2(CUT))


class TestSecondSheet:
    second = ka.second_sheet(PION_MASS)(amplitude)
    first = ka.first_sheet(PION_MASS)(second)

    def test_identity(self):
        cls = type(self)
        assert np.allclose(cls.first(CUT), amplitude(CUT))
