"""PT-table generator for radiative-transfer codes."""
import logging
from typing import List

import numpy as np

from .core import MARFA
from .hitran import SpectralLine

logger = logging.getLogger(__name__)


class PTTableGenerator:
    """Generate a (P, T, nu) lookup table of absorption coefficients."""

    def __init__(self, marfa: MARFA):
        self.marfa = marfa

    def generate_table(
        self,
        nu_grid: np.ndarray,
        lines: List[SpectralLine],
        P_grid: np.ndarray,
        T_grid: np.ndarray,
        mole_fraction: float,
        output_file: str,
    ) -> np.ndarray:
        logger.info("Generating PT table: %d P x %d T x %d nu",
                    len(P_grid), len(T_grid), len(nu_grid))

        n_P, n_T, n_nu = len(P_grid), len(T_grid), len(nu_grid)
        table = np.zeros((n_P, n_T, n_nu))

        for i, P in enumerate(P_grid):
            for j, T in enumerate(T_grid):
                logger.debug("P=%.3e atm  T=%.1f K", P, T)
                table[i, j, :] = self.marfa.calculate_absorption_coefficient(
                    nu_grid, lines, T, P, mole_fraction)

        np.save(output_file, table)
        logger.info("Saved PT table to %s", output_file)
        return table