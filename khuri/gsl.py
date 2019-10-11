from enum import Enum

from _khuri_gsl import *


class Integrators(Enum):
    cquad = 1
    qag = 2
