import itertools

import pytest
import numpy as np

import khuri.amplitude as ka
from khuri.breit_wigner import partial_wave
from khuri.tests.helpers import connected


PION_MASS = 0.14
CUT = np.linspace(2.0 * PION_MASS + 1, 4.0 * PION_MASS, 50)**2


def amplitude(mandelstam_s):
    return partial_wave(mandelstam_s,
                        resonance_mass=0.77,
                        resonance_width=0.15,
                        angular_momentum=1,
                        mass=PION_MASS)


def phase(mandelstam_s):
    return np.angle(amplitude(mandelstam_s))


@ka.from_phase(PION_MASS)
def amplitude2(mandelstam_s):
    return phase(mandelstam_s)


@ka.from_cot(PION_MASS)
def amplitude3(mandelstam_s):
    return 1.0 / np.tan(phase(mandelstam_s))


class Phase:
    shift = np.pi

    def __init__(self, shift):
        self.shift = shift

    @classmethod
    @ka.from_phase(PION_MASS)
    def static(cls, mandelstam_s):
        return phase(mandelstam_s) + cls.shift

    @ka.from_cot(PION_MASS)
    def non_static(self, mandelstam_s):
        return 1.0 / np.tan(phase(mandelstam_s) + self.shift)


def test_method():
    instance = Phase(shift=np.pi)
    assert instance.non_static(7.0) < np.inf


def test_static_method():
    assert Phase.static(7.0) < np.inf


def test_consistency():
    amplitudes = (amplitude, amplitude2, amplitude3)
    for amp1, amp2 in itertools.combinations(amplitudes, 2):
        assert np.allclose(amp1(CUT), amp2(CUT))


class TestSecondSheet:
    second = ka.second_sheet(PION_MASS)(amplitude)
    first = ka.first_sheet(PION_MASS)(second)

    def test_identity(self):
        cls = type(self)
        assert np.allclose(cls.first(CUT), amplitude(CUT))

    def test_connected(self):
        cls = type(self)
        connected(cls.first, cls.second, CUT)
