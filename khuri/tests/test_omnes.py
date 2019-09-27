import pytest

from khuri.omnes import Omnes


def phase(s):
    return 1.0 + 2.0 / s


def test_omnes_at_zero():
    o = Omnes(phase, threshold=4.0)
    assert o(0) == pytest.approx(1.0)
