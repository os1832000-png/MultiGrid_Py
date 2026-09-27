"""Physical constants matching Fortran MARFA."""
import numpy as np


class Constants:
    """Physical constants - matching Fortran MARFA exactly."""

    c         = 2.99792458e10    # speed of light [cm/s]
    h         = 6.62606957e-27   # Planck constant [erg·s]
    k_B       = 1.380649e-16     # Boltzmann constant [erg/K]
    k_B_SI    = 1.380649e-23     # Boltzmann constant [J/K]
    Na        = 6.02214129e23    # Avogadro constant [1/mol]
    R         = 8.314472e7       # Gas constant [erg/(mol·K)]
    c2        = 1.4388           # second radiation constant [cm·K]
    T_ref     = 296.0            # reference temperature [K]
    P_ref     = 1.0              # reference pressure [atm]
    atm_to_pa = 101325.0         # 1 atm = 101325 Pa
    sqrt_pi   = np.sqrt(np.pi)
    sqrt_ln2  = np.sqrt(np.log(2.0))
    ln2_pi    = np.log(2.0) / np.pi