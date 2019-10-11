from khuri.gsl import Settings, IntegrationRoutine
from _khuri_omnes import *


def generate_omnes(*args,
                   integration_routine=IntegrationRoutine.cquad,
                   **kwargs):
    """Generate an Omnes function.

    Parameters
    ----------
    phase: callable
        the phase of the Omnes function
    threshold: float
        the threshold in the s-plane
    constant: float, optional
        if specified, the phase is set to constant above `cut`
    cut: float, optional
        if specified, the phase is set to constant above `cut`
    minimal_distance: float, optional
        half the width of a band around the cut. For arguments of the Omnes
        function in this band, a different prescription is used for the
        evaluation of the Omnes function to take care of the singularity in the
        integral.
    conf: Settings, optional
        the settings for the integration routine
    integration_routine: IntegrationRoutine, optional
        specify an adaptive integration routine to be used to compute the
        integral

    Returns
    -------
    A callable that yields the values of the requested Omnes function.
    """
    if integration_routine == IntegrationRoutine.cquad:
        return OmnesCquad(*args, **kwargs)
    if integration_routine == IntegrationRoutine.qag:
        return OmnesQag(*args, **kwargs)
    raise ValueError('unknown integration routine')
