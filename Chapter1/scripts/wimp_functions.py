
import numpy as np
import scipy
from scipy.interpolate import interp1d
from scipy.special import erf


def getWIMPRateIsotope(
        isotope,
        wimp_mass,
        wimp_cross_section,
        recoil_energy,
        earth_velocity,
        wimp_velocity_dispersion,
        escaping_velocity,
        wimp_local_density,
):
    m_fU = 931.494; # unified atomic mass in MeV
    Mt = m_fU * isotope / 1000.
    mu = Mt*wimp_mass / (wimp_mass + Mt)
    m_fVc = 2.9979e5
    m_fK = (3600*24*m_fVc*1e5)*( np.power(m_fVc,2.0) *1e6/(1.6e-10))*1e-6 # conversion factor for the rate to be in unit of dru
    return m_fK * wimp_local_density / ( 2.*np.power(mu,2.) * wimp_mass) * getEta(
        isotope,
        wimp_mass,
        recoil_energy,
        earth_velocity,
        wimp_velocity_dispersion,
        escaping_velocity) * m_fVc * getNuclearSigma(isotope, wimp_mass, wimp_cross_section, recoil_energy)

def getEta(
        isotope,
        wimp_mass,
        recoil_energy,
        earth_velocity,
        wimp_velocity_dispersion,
        escaping_velocity
):
    m_fU = 931.494
    m_fVc = 2.9979e5
    Mt = m_fU * isotope / 1000.
    mu = wimp_mass * Mt / (wimp_mass + Mt)
    v_min = (Mt*recoil_energy / mu) / np.sqrt(2.*Mt*recoil_energy) * m_fVc / 1000.
    x = v_min / wimp_velocity_dispersion
    y = earth_velocity / wimp_velocity_dispersion
    z = escaping_velocity / wimp_velocity_dispersion
    Nesc = scipy.special.erf(z) - 2.*z*np.exp(-z*z) / np.sqrt(np.pi)
    if x<z-y:
        return 1./(2.*Nesc*earth_velocity)*(erf(x+y) - erf(x-y) - 4./np.sqrt(np.pi)*y*np.exp(-z*z))
    elif x>z-y and x<z+y:
        return 1./(2.*Nesc*earth_velocity)*(erf(z) - erf(x-y) - 2./np.sqrt(np.pi)*(y+z-x)*np.exp(-z*z))
    else:
        return 0

def getNuclearSigma(
        isotope,
        wimp_mass,
        wimp_cross_section,
        recoil_energy,
):
    m_fU = 931.494
    Mt = m_fU*isotope/1000.
    mu = Mt*wimp_mass / (wimp_mass + Mt)
    mu_p = 0.938272*wimp_mass / (0.938272 + wimp_mass)
    return np.power(isotope*mu/mu_p, 2.) * wimp_cross_section * np.power(getFormFactor(isotope, wimp_mass, recoil_energy), 2.)

def getFormFactor(
        isotope,
        wimp_mass,
        recoil_energy,
):
    m_fU = 931.494
    Mt = m_fU * isotope / 1000.
    q = np.sqrt(2*Mt*recoil_energy)/ 197.326
    a = 0.52
    s = 0.9
    c = 1.23*np.power(isotope, 1/3.) - 0.60
    r_n = np.sqrt(c*c + 7/3.*np.pi*np.pi*a*a-5.*s*s)
    FF_SI= 3.*np.exp(-q*s*q*s/2.) *(np.sin(q*r_n)-q*r_n*np.cos(q*r_n))/np.power(q*r_n, 3.0)
    #print recoil_energy, c
    return FF_SI