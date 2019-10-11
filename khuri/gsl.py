from enum import Enum

from _khuri_gsl import *


class IntegrationRoutine(Enum):
    """The different adaptive integration routines.

    `cquad` is more potent in solving slowly converging or otherwise hard
    integrals, while `qag` uses less resources.
    """
    cquad = 1
    qag = 2
