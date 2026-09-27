"""Total Internal Partition Sums (TIPS)."""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class TIPS:
    """
    Total Internal Partition Sums.

    NOTE
    ----
    The temperature dependence implemented here is a rough T**1.5 scaling
    anchored at 296 K.  It is adequate for qualitative work but is NOT the
    real Gamache et al. (2017) TIPS database.  Replace `_calculate_tips_table`
    with a real table lookup for quantitative accuracy.
    """

    def __init__(self):
        self.Q_ref_296 = {
            1: 178.12,  2: 289.49,  3: 4870.3,  4: 1122.3,
            5: 108.58,  6: 590.43,  7: 216.21,  8: 159.47,
            9: 5792.6,  10: 2379.5, 11: 169.24, 12: 11456.0,
        }
        self.T_grid = np.arange(20, 1002, 2)
        self._calculate_tips_table()

    def _calculate_tips_table(self):
        self.tips_table = {}
        for mol_id in range(1, 13):
            Q_ref = self.Q_ref_296.get(mol_id, 100.0)
            self.tips_table[mol_id] = Q_ref * (self.T_grid / 296.0) ** 1.5

    def get_Q(self, mol_id: int, T: float) -> float:
        if mol_id not in self.tips_table:
            logger.debug("Unknown molecule id %d, falling back to CO2", mol_id)
            mol_id = 2
        T = float(np.clip(T, 20.0, 1000.0))
        idx = int((T - 20.0) / 2.0)
        idx = max(0, min(idx, len(self.T_grid) - 2))
        T1, T2 = self.T_grid[idx], self.T_grid[idx + 1]
        Q1, Q2 = self.tips_table[mol_id][idx], self.tips_table[mol_id][idx + 1]
        return float(Q1 + (Q2 - Q1) * (T - T1) / (T2 - T1))