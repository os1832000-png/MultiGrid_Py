# ------------------------------------------------------------------
    # absorption cross-section -- with full LBL multi-grid cascade
    # ------------------------------------------------------------------

    def calculate_absorption_cross_section(self,
                                           nu_grid: np.ndarray,
                                           lines: List[SpectralLine],
                                           T: float,
                                           P: float,
                                           self_broadening_fraction: float = 0.0,
                                           line_cutoff_cm: float = 25.0,
                                           wing_correction: str = 'none'
                                           ) -> np.ndarray:
        """
        Calculate monochromatic absorption cross-section sigma(nu) [cm2/molecule].

        Uses a direct line-by-line Voigt summation on the supplied nu_grid.
        For each spectral line, the Voigt profile is evaluated at all grid
        points within line_cutoff_cm of the line center and accumulated into
        sigma.

        Parameters
        ----------
        nu_grid   : wavenumber grid [cm-1]
        lines     : list of SpectralLine objects (from HITRANReader)
        T         : temperature [K]
        P         : total pressure [atm]
        self_broadening_fraction : fraction of P that is self-pressure (0-1)
        line_cutoff_cm : wing cutoff [cm-1] -- profile set to zero beyond this
        wing_correction : 'none' | 'tonkov' | 'perrin' (sub-Lorentzian chi-factor)

        Returns
        -------
        sigma : np.ndarray, shape (len(nu_grid),), units cm2/molecule
                Integral of sigma over all nu equals the line intensity S(T).
        """
        if not lines:
            return np.zeros_like(nu_grid)

        sigma  = np.zeros_like(nu_grid, dtype=float)
        P_self = P * self_broadening_fraction
        mol_id = lines[0].mol_id

        chi_func = {
            'tonkov': WingCorrections.tonkov_chi,
            'perrin': WingCorrections.perrin_chi,
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

            if wing_correction != 'none':
                profile = profile * chi_func(x, line.nu)

            sigma[mask] += S_T * profile

        return sigma

    def calculate_absorption_coefficient(self,
                                         nu_grid: np.ndarray,
                                         lines: List[SpectralLine],
                                         T: float,
                                         P: float,
                                         mole_fraction: float,
                                         self_broadening_fraction: float = 0.0,
                                         line_cutoff_cm: float = 25.0,
                                         wing_correction: str = 'none'
                                         ) -> np.ndarray:
        """
        Monochromatic volume absorption coefficient alpha(nu) [cm-1].

        alpha(nu) = sigma(nu) * n_species
        where n_species = n_total * mole_fraction  [molecules / cm3]

        Parameters
        ----------
        nu_grid           : wavenumber grid [cm-1]
        lines             : list of SpectralLine objects
        T                 : temperature [K]
        P                 : total pressure [atm]
        mole_fraction     : volume mixing ratio of absorbing species (0-1)
        self_broadening_fraction : fraction of P for self-broadening (0-1)
        line_cutoff_cm    : wing cutoff [cm-1]
        wing_correction   : 'none' | 'tonkov' | 'perrin'

        Returns
        -------
        alpha : np.ndarray, shape (len(nu_grid),), units cm-1
        """
        sigma   = self.calculate_absorption_cross_section(
                      nu_grid, lines, T, P,
                      self_broadening_fraction=self_broadening_fraction,
                      line_cutoff_cm=line_cutoff_cm,
                      wing_correction=wing_correction)
        # Number density of total air [molecules / cm3]
        # n = P [Pa] / (k_B [J/K] * T [K])  ->  [molecules / m3]  /1e6 -> [/cm3]
        n_total = (P * self.const.atm_to_pa) / (self.const.k_B_SI * T) / 1e6
        return sigma * n_total * mole_fraction

