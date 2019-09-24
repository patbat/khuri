"""
==========================================================
The I=J=1 pi pi->pi pi Madrid partial wave parametrization
==========================================================

The phrase "physical units", which is used several times in the documentation,
refers to a consistent choice of units (e.g. all
masses and energies etc. expressed in MeV or everything expressed in powers of
GeV).
"""
import numpy as np

class Parameters:
    """The parameters of the Madrid phase."""
    # constraint fit parameters
    b0_cfd = 1.043
    b1_cfd = 0.19
    lambda1_cfd = 1.39
    lambda2_cfd = -1.70

    # unconstraint fit parameters
    b0_ufd = 1.055
    b1_ufd = 0.15
    lambda1_ufd = 1.57
    lambda2_ufd = -1.96

    s0 = 1050.0**2 # in (MeV)^2

    # The masses are parameters that determine the shape of the phase.
    # They are not necessarily the physical masses of the particles.
    kaon_mass = 496.0 # in MeV
    rho_mass = 773.6 # in MeV
    pion_mass = 139.57 # in MeV

    # The parametrization is not valid for energies below lowest
    # or above highest.
    lowest = 2.0*pion_mass # in MeV
    highest = 1420 # in MeV


def _momentum(s, m):
    """The center of mass momentum.

    Parameters
    ----------
    s: number
        the squared cm energy in physical units
    m: number
        the pion mass in physical units
    """
    return np.sqrt(0.25*s - m**2)


def _w(s, s0=Parameters.s0):
    """Simplify the calculation of the phase.

    Parameters
    ----------
    s: number
        the squared cm energy in physical units
    s0: number, optional
        parameter, needs to have the same units as s
        (default value is given in MeV^2)
    """
    a = np.sqrt(s)
    b = np.sqrt(s0 - s)
    return  (a-b) / (a+b)


def phase_low(s, m_pi, m_rho, b0, b1, s0=Parameters.s0):
    """Calculate the phase for s <= 4m_K^2

    Parameters
    ----------
    s: number
        the squared cm energy in physical units
    m_pi: number
        the pion mass in physical units
    m_rho: number
        the rho mass in physical units
    b0: number
        a parameter determined by fits
    b1: number
        a parameter determined by fits
    s0: number, optional
        parameter, needs to have the same units as s
        (default value is given in MeV^2)
    """
    en = np.sqrt(s)
    coeff = en * (m_rho**2-s) / (2*_momentum(s, m_pi)**3)
    summ = 2*m_pi**3 / (m_rho**2*en) + b0 + b1*_w(s, s0)
    return np.arctan2(1, coeff*summ)


def phase_high(s, m_k, lambda_0, lambda_1, lambda_2):
    """Calculate the phase for (1420 MeV)^2 > s > 4m_K^2

    Parameters
    ----------
    s: number
        the squared cm energy in physical units
    m_k: number
        the kaon mass in physical units
    lambda_0: number
        a parameter determined by fits
    lambda_1: number
        a parameter determined by fits
    lambda_2: number
        a parameter determined by fits
    """
    temp = np.sqrt(s)/(2*m_k) - 1
    return lambda_0 + lambda_1*temp + lambda_2*temp**2


def lambda0(m_pi, m_k, m_rho, b0, b1, s0):
    s = 4.0*m_k**2
    return phase_low(s,m_pi,m_rho,b0,b1,s0)


@np.vectorize
def phase(
        s,
        m_pi=Parameters.pion_mass,
        m_k=Parameters.kaon_mass,
        m_rho=Parameters.rho_mass,
        b0=Parameters.b0_cfd,
        b1=Parameters.b1_cfd,
        lambda_1=Parameters.lambda1_cfd,
        lambda_2=Parameters.lambda2_cfd,
        s0=Parameters.s0):
    """Calculate the Madrid phase for the I=J=1 pi pi-> pi pi partial wave.

    Parameters
    ----------
    s: number
        the squared cm energy in physical units
    m_pi: number
        the pion mass in physical units
    m_k: number
        the kaon mass in physical units
    m_rho: number
        the rho mass in physical units
    b0: number
        a parameter determined by fits
    b1: number
        a parameter determined by fits
    lambda_1: number
        a parameter determined by fits
    lambda_2: number
        a parameter determined by fits
    s0: number, optional
        parameter, needs to have the same units as s
        (default value is given in MeV^2)
    """
    if s < 4*m_pi**2:
        print("below threshold")
    elif s <= 4*m_k**2:
        return phase_low(s, m_pi, m_rho, b0, b1, s0)
    elif s <= 1420**2:
        l0 = lambda0(m_pi, m_k, m_rho, b0, b1, s0)
        return phase_high(s, m_k, l0, lambda_1, lambda_2)
    else:
        print("above region, in which the parametrization of the phase is\
                valid")
