# ==============================================================================
# PT-TABLE GENERATOR
# ==============================================================================

class PTTableGenerator:
    def __init__(self, marfa: MARFA):
        self.marfa = marfa

    def generate_table(self,
                       nu_grid: np.ndarray,
                       lines: List[SpectralLine],
                       P_grid: np.ndarray,
                       T_grid: np.ndarray,
                       mole_fraction: float,
                       output_file: str):
        print(f"\nGenerating PT table...")
        print(f"  Pressures   : {len(P_grid)} points")
        print(f"  Temperatures: {len(T_grid)} points")
        print(f"  Wavenumbers : {len(nu_grid)} points")

        n_P, n_T, n_nu = len(P_grid), len(T_grid), len(nu_grid)
        table = np.zeros((n_P, n_T, n_nu))

        for i, P in enumerate(P_grid):
            for j, T in enumerate(T_grid):
                print(f"  P={P:.3e} atm  T={T:.1f} K", end='\r')
                table[i, j, :] = self.marfa.calculate_absorption_coefficient(
                    nu_grid, lines, T, P, mole_fraction)

        print(f"\n  Table complete")
        np.save(output_file, table)
        print(f"  Saved: {output_file}")
        return table

