"""
===================================
High energy continuation of phases.
===================================
"""
import numpy as np
from scipy.misc import derivative


def asymptotic1(matching_point, limit=np.pi):
    """Return phase that approaches a constant at high energies.

    The decorated phase is used below `matching_point`, above which `limit`
    is approached eventually. The transition is one times differentiable.
    """
    def wrapper(phase):
        value = phase(matching_point)
        d_value = derivative(phase, matching_point, dx=1e-6)
        parameter1 = (limit - value)**2 / matching_point / d_value
        parameter2 = (limit - value) / matching_point / d_value - 1

        if parameter2 <= -1:
            raise ValueError("the continuation of the phase to high energies"
                             " exhibits a pole above the matching point")

        def continued_phase(s):
            if s < matching_point:
                return phase(s)
            return limit - parameter1 / (parameter2 + s / matching_point)

        return np.vectorize(continued_phase, doc=phase.__doc__)

    return wrapper


def asymptotic2(lower, upper, limit=np.pi):
    """Return phase that equals `phase` below `lower` and `limit` above `upper`.
    
    Between `lower` and `upper` there is a in [`lower`,`upper`) one
    times continuously differentiable smooth connection.
    """
    def wrapper(phase):
        def continued_phase(s):
            if s < lower:
                return phase(s)
            if s < upper:
                c = _connect(lower, upper, s)
                return (1.0 - c) * phase(s) + c * limit
            return limit

        return np.vectorize(continued_phase, doc=phase.__doc__)

    return wrapper


def _connect(lower, upper, s):
    a = (s - lower)**2
    b = 3.0*upper - 2.0*s - lower
    c = upper - lower
    return a * b / c**3
