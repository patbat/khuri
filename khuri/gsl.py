from enum import Enum

from _khuri_gsl import *


class IntegrationRoutine(Enum):
    """The different adaptive integration routines.

    `cquad` is more potent in solving slowly converging or otherwise hard
    integrals, while `qag` uses less resources, but is also much less reliable.
    In general, `cquad` is strongly preferred.
    """
    cquad = 1
    qag = 2
