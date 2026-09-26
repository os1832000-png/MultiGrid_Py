"""Core MARFA calculation engine."""
import logging
from typing import List

import numpy as np

from .constants import Constants
from .molecules import MOLECULAR_MASSES
from .tips import TIPS
from .line_shapes import LineShapes
from .wing_corrections import WingCorrections
from .hitran import SpectralLine

logger = logging.getLogger(__name__)


class MARFA:
    """
    Main MARFA calculation engine.

    Wraps the LBL engine with temperature-dependent intensities,
    pressure broadening, and chi-factors.
    """

    def __init__(self, *, verbose: bool = True):
        self.const = Constants()
        self.tips  = TIPS()
        self.verbose = verbose

    # ------------------------------------------------------------------
    # line parameter helpers
    # ------------------------------------------------------------------

    def calculate_doppler_width(self, nu0: float, T: float, mol_id: int) -> float:
        M = MOLECULAR_MASSES.get(mol_id, 44.0)
        m = M / self.const.Na
        return nu0 * np.sqrt(2.0 * self.const.k_B * T * np.log(2.0) /
                             (m * self.const.c ** 2))

    def calculate_lorentz_width(self,
                                line: SpectralLine,
                                P: float,
                                P_self: float,
                                T: float) -> float:
        pressure_term = line.gamma_air * (P - P_self) + line.gamma_self * P_self
        temperature_factor = (self.const.T_ref / T) ** line.n_air
        return pressure_term * temperature_factor

    def calculate_line_intensity_at_T(self,
                                      line: SpectralLine,
                                      T: float,
                                      P: float) -> float:
        nu_shifted = line.nu + line.delta_air * P
        Q_T        = self.tips.get_Q(line.mol_id, T)
        Q_ref      = self.tips.get_Q(line.mol_id, self.const.T_ref)
        partition_ratio = Q_ref / Q_T if Q_T > 0 else 1.0

        boltzmann_ratio = (np.exp(-self.const.c2 * line.E_low / T) /
                           np.exp(-self.const.c2 * line.E_low / self.const.T_ref))

        c2_nu_T = self.const.c2 * nu_shifted / T
        c2_nu_r = self.const.c2 * nu_shifted / self.const.T_ref
        emission_ratio = (1.0 if c2_nu_T > 50 else
                          (1.0 - np.exp(-c2_nu_T)) /
                          (1.0 - np.exp(-c2_nu_r)))

        return line.S * partition_ratio * boltzmann_ratio * emission_ratio

    # ------------------------------------------------------------------
    # absorption cross-section / coefficient
    # ------------------------------------------------------------------

    def calculate_absorption_cross_section(
        self,
        nu_grid: np.ndarray,
        lines: List[SpectralLine],
        T: float,
        P: float,
        self_broadening_fraction: float = 0.0,
        line_cutoff_cm: float = 25.0,
        wing_correction: str = "none",
    ) -> np.ndarray:
        """
        Monochromatic absorption cross-section sigma(nu) [cm2/molecule].

        Uses a direct line-by-line Voigt summation on ``nu_grid``.  Integral
        of sigma over all nu equals the line intensity S(T).
        """
        if not lines:
            return np.zeros_like(nu_grid)

        sigma  = np.zeros_like(nu_grid, dtype=float)
        P_self = P * self_broadening_fraction
        mol_id = lines[0].mol_id

        chi_func = {
            "tonkov": WingCorrections.tonkov_chi,
            "perrin": WingCorrections.perrin_chi,
        }.get(wing_correction, WingCorrections.no_correction)

        for line in lines:
            S_T = self.calculate_line_intensity_at_T(line, T, P)
            if S_T <= 0.0:
                continue

            gamma_D = self.calculate_doppler_width(line.nu, T, mol_id)
            gamma_L = self.calculate_lorentz_width(line, P, P_self, T)
            if max(gamma_D, gamma_L) <= 0.0:
                continue

            mask = np.abs(nu_grid - line.nu) <= line_cutoff_cm
            if not np.any(mask):
                continue

            x       = nu_grid[mask] - line.nu
            profile = LineShapes.voigt_normalized(x, gamma_L, gamma_D)

            if wing_correction != "none":
                profile = profile * chi_func(x, line.nu)

            sigma[mask] += S_T * profile

        return sigma

    def calculate_absorption_coefficient(
        self,
        nu_grid: np.ndarray,
        lines: List[SpectralLine],
        T: float,
        P: float,
        mole_fraction: float,
        self_broadening_fraction: float = 0.0,
        line_cutoff_cm: float = 25.0,
        wing_correction: str = "none",
    ) -> np.ndarray:
        """
        Monochromatic volume absorption coefficient alpha(nu) [cm-1].

        alpha(nu) = sigma(nu) * n_species
        where n_species = n_total * mole_fraction [molecules / cm3].
        """
        sigma = self.calculate_absorption_cross_section(
            nu_grid, lines, T, P,
            self_broadening_fraction=self_broadening_fraction,
            line_cutoff_cm=line_cutoff_cm,
            wing_correction=wing_correction,
        )

        # Total air number density [molecules / cm3]
        n_total = (P * self.const.atm_to_pa) / (self.const.k_B_SI * T) / 1e6
        return sigma * n_total * mole_fraction