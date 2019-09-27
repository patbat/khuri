import numpy as np
import pytest

import khuri.madrid
import khuri.phases


MATCHING_POINT = 1000**2


@khuri.phases.asymptotic1(matching_point=MATCHING_POINT)
def phase(s):
    return khuri.madrid.phase(s)


class TestAsymptotic1:
    def test_limit(self):
        """Test if the phase approaches the desired limit."""
        assert phase(1e12) == pytest.approx(np.pi)

    def test_agreement(self):
        """Test agreement with the original phase."""
        values = np.linspace(300**2, MATCHING_POINT - 10, 500)
        original = khuri.madrid.phase(values)
        wrapper = phase(values)
        assert np.all(original == wrapper)
