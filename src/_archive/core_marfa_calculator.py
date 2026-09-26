# ==============================================================================
# CORE MARFA CALCULATOR
# ==============================================================================

class MARFA:
    """
    Main MARFA calculation engine.
    Integrates the full multi-grid LBL cascade (LineGridCalc) with
    temperature-dependent intensities, pressure broadening, and chi-factors.
    """

    def __init__(self):
        self.const = Constants()
        self.tips  = TIPS()

    # ------------------------------------------------------------------
    # line parameter helpers
    # ------------------------------------------------------------------

    def calculate_doppler_width(self, nu0: float, T: float, mol_id: int) -> float:
        M = MOLECULAR_MASSES.get(mol_id, 44.0)
        m = M / self.const.Na
        return nu0 * np.sqrt(2.0 * self.const.k_B * T * np.log(2.0) /
                             (m * self.const.c**2))

    def calculate_lorentz_width(self, line: SpectralLine,
                                P: float, P_self: float, T: float) -> float:
        pressure_term      = line.gamma_air * (P - P_self) + line.gamma_self * P_self
        temperature_factor = (self.const.T_ref / T) ** line.n_air
        return pressure_term * temperature_factor

    def calculate_line_intensity_at_T(self, line: SpectralLine,
                                      T: float, P: float) -> float:
        nu_shifted = line.nu + line.delta_air * P
        Q_T        = self.tips.get_Q(line.mol_id, T)
        Q_ref      = self.tips.get_Q(line.mol_id, self.const.T_ref)
        partition_ratio  = Q_ref / Q_T if Q_T > 0 else 1.0
        boltzmann_ratio  = (np.exp(-self.const.c2 * line.E_low / T) /
                            np.exp(-self.const.c2 * line.E_low / self.const.T_ref))
        c2_nu_T = self.const.c2 * nu_shifted / T
        c2_nu_r = self.const.c2 * nu_shifted / self.const.T_ref
        emission_ratio = (1.0 if c2_nu_T > 50 else
                          (1.0 - np.exp(-c2_nu_T)) / (1.0 - np.exp(-c2_nu_r)))
        return line.S * partition_ratio * boltzmann_ratio * emission_ratio
