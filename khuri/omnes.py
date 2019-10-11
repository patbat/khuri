from khuri.gsl import Settings, Integrators
from _khuri_omnes import *


def generate_omnes(*args, integration_routine=Integrators.cquad, **kwargs):
    if integration_routine == Integrators.cquad:
        return OmnesCquad(*args, **kwargs)
    if integration_routine == Integrators.qag:
        return OmnesQag(*args, **kwargs)
    raise ValueError('unknown integration routine')
