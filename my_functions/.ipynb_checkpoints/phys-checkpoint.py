# physical functions
import numpy as np
from scipy.integrate import quad

def encl_mass(dens, r, *kwargs):
    """
    Return total mass enclosed in r assuming a spherically
    symmetric density profile

    Inputs:
    * dens (function, density profile. Should be of the form dens(r, *kwargs)).
    * r (float, radius at which enclosed mass should be evaluated)
    * + all arguments needed by the dens function

    Outputs:
    * Mass enclosed in r (float)
    """
    return quad(lambda r_: 4*np.pi*r_**2*dens(r_, *kwargs), 0, r)[0]

def Burkert(r, r_s, rho_s):
    """
    Burkert density profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)

    Outputs:
    density (float)
    """
    x = r/r_s
    return rho_s/((1 + x)*(1 + x**2))

def tNFW(r, r_s, rho_s, t):
    """
    Truncated NFW profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)

    Outputs:
    density (float)
    """
    x = r/r_s
    return (rho_s/(x * (1+x)**2)) * (t**2/(x**2 + t**2)) 

def NFW(r, r_s, rho_s):
    """
    NFW profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)
    - or - 
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * M(r)
    

    Outputs:
    density (float)
    """
    x = r/r_s
    return rho_s/(x * (1+x)**2)

def mass_NFW(r, r_s, rho_s):
    """
    Mass enclosed in r for an NFW profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)

    Outputs:
    enclosed mass (float)
    """
    return 4 * np.pi * rho_s * r_s**3 * ( np.log( (r_s + r) /r_s ) - r/(r_s + r))

def Hernquist(r, r_s, rho_s):
    """
    Hernquist profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)

    Outputs:
    density (float)
    """
    x = r/r_s
    return rho_s/(x * (1+x)**3)

def mass_Hernquist(r, r_s, rho_s):
    """
    Mass enclosed in r for a Hernquist profile

    Inputs:
    * r (float, radius variable)
    * r_s (float, scale radius parameter)
    * rho_s (float, scale density parameter)

    Outputs:
    enclosed mass (float)
    """
    x = r/r_s
    return 2 * np.pi * rho_s * r_s**3 * x**2/(x+1)**2

def Beta_model(r, r_c, rho_0, beta):
    """
    Beta model density profile

    Inputs:
    * r (float, radius variable)
    * r_c (float, scale radius parameter)
    * rho_0 (float, central density)
    * beta (float, curvature parameter)

    Outputs:
    density (float)
    """
    return rho_0 * (1 + (r/r_c)**2 )**(-3*beta/2)

def inverse_CPD_Hernquist(x, r_s, rho_s, r_max):
    """
    Inverse commulative probability distribution for a Hernquist profile
    (used for sampling particles)

    Inputs:
    * x (float, number between 0 and 1)
    * r_s (float, scale radius)
    * rho_s (float, scale density)
    * r_max (float, maximum sampling radius in physical units)

    Outputs:
    radius in physical units (float)
    """
    x *= (r_max/r_s)**2/( (r_max/r_s) +1)**2
    return - r_s * np.sqrt(x)/(np.sqrt(x)-1)