import numpy as np
import pytest

from khuri.schenk import partial_wave_lit, tan_phase_lit


ISOSPINS = (0, 1, 2)
ISOSPIN_IDS = tuple(f'isospin {i}' for i in ISOSPINS)


@pytest.mark.parametrize('isospin', ISOSPINS, ids=ISOSPIN_IDS)
def test(isospin):
    mandelstam_s = np.linspace(0.3, 0.5)**2
    pion_mass = 0.135
    tan = tan_phase_lit(mandelstam_s, isospin, pion_mass)
    partial_wave = partial_wave_lit(mandelstam_s, isospin, pion_mass)
    tan2 = np.tan(np.angle(partial_wave))
    assert np.allclose(tan, tan2)


@pytest.mark.parametrize('func', (partial_wave_lit, tan_phase_lit),
                         ids=('partial wave', 'tangent'))
def test_isospin_4(func):
    with pytest.raises(ValueError):
        func(0.3**2, 4, 0.135)
