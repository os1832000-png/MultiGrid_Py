"""Multi-resolution grid state (Fortran GridState / Grids module)."""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class GridState:
    """
    Multi-resolution stencil arrays for the LBL multi-grid algorithm.

    For each level n in [0..9] there are three arrays:
        RKnP[i] -- 'plus'   neighbour contribution at grid point i
        RKn[i]  -- 'center' contribution at grid point i
        RKnL[i] -- 'left'   neighbour contribution at grid point i

    NT0..NT9 are the number of grid points at each level.
    RK[i] is the coarsest flat grid (NT points).
    """

    def __init__(self, delta_wv: float, cut_off: float, H: float):
        self.deltaWV = delta_wv
        self.cutOff  = cut_off
        self.H       = H

        self.H0 = 2.0    * H
        self.H1 = 4.0    * H
        self.H2 = 8.0    * H
        self.H3 = 16.0   * H
        self.H4 = 32.0   * H
        self.H5 = 64.0   * H
        self.H6 = 128.0  * H
        self.H7 = 256.0  * H
        self.H8 = 512.0  * H
        self.H9 = 1024.0 * H

        self.NT  = max(1, int(round(delta_wv / H)))
        self.NT0 = max(1, self.NT // 2)
        self.NT1 = max(1, self.NT // 4)
        self.NT2 = max(1, self.NT // 8)
        self.NT3 = max(1, self.NT // 16)
        self.NT4 = max(1, self.NT // 32)
        self.NT5 = max(1, self.NT // 64)
        self.NT6 = max(1, self.NT // 128)
        self.NT7 = max(1, self.NT // 256)
        self.NT8 = max(1, self.NT // 512)
        self.NT9 = max(1, self.NT // 1024)

        self._allocate()

    def _allocate(self):
        def z(n):
            return np.zeros(n + 2)   # +2 for 1-based index safety

        self.RK = z(self.NT)

        self.RK0P = z(self.NT0);  self.RK0 = z(self.NT0);  self.RK0L = z(self.NT0)
        self.RK1P = z(self.NT1);  self.RK1 = z(self.NT1);  self.RK1L = z(self.NT1)
        self.RK2P = z(self.NT2);  self.RK2 = z(self.NT2);  self.RK2L = z(self.NT2)
        self.RK3P = z(self.NT3);  self.RK3 = z(self.NT3);  self.RK3L = z(self.NT3)
        self.RK4P = z(self.NT4);  self.RK4 = z(self.NT4);  self.RK4L = z(self.NT4)
        self.RK5P = z(self.NT5);  self.RK5 = z(self.NT5);  self.RK5L = z(self.NT5)
        self.RK6P = z(self.NT6);  self.RK6 = z(self.NT6);  self.RK6L = z(self.NT6)
        self.RK7P = z(self.NT7);  self.RK7 = z(self.NT7);  self.RK7L = z(self.NT7)
        self.RK8P = z(self.NT8);  self.RK8 = z(self.NT8);  self.RK8L = z(self.NT8)
        self.RK9P = z(self.NT9);  self.RK9 = z(self.NT9);  self.RK9L = z(self.NT9)

    def reset(self):
        """Zero all arrays (called before each new spectral line)."""
        self._allocate()

    def reconstruct_fine_grid(self) -> np.ndarray:
        """Reconstruct the fine-grid array from the stencil hierarchy."""
        result = self.RK[1:self.NT + 1].copy()

        levels = [
            (self.RK0, self.RK0P, self.RK0L, self.NT0, self.H0),
            (self.RK1, self.RK1P, self.RK1L, self.NT1, self.H1),
            (self.RK2, self.RK2P, self.RK2L, self.NT2, self.H2),
            (self.RK3, self.RK3P, self.RK3L, self.NT3, self.H3),
            (self.RK4, self.RK4P, self.RK4L, self.NT4, self.H4),
            (self.RK5, self.RK5P, self.RK5L, self.NT5, self.H5),
            (self.RK6, self.RK6P, self.RK6L, self.NT6, self.H6),
            (self.RK7, self.RK7P, self.RK7L, self.NT7, self.H7),
            (self.RK8, self.RK8P, self.RK8L, self.NT8, self.H8),
            (self.RK9, self.RK9P, self.RK9L, self.NT9, self.H9),
        ]

        NT = self.NT
        for RKn, RKnP, RKnL, NTn, Hn in levels:
            ratio = max(1, int(round(Hn / self.H)))
            i_arr = np.arange(1, NTn + 1)
            c_idx = i_arr * ratio - 1
            l_idx = c_idx - ratio // 2
            r_idx = c_idx + ratio // 2

            for idx_arr, RK_arr in ((c_idx, RKn),
                                    (l_idx, RKnL),
                                    (r_idx, RKnP)):
                valid = (idx_arr >= 0) & (idx_arr < NT)
                if np.any(valid):
                    np.add.at(result, idx_arr[valid],
                              RK_arr[1:NTn + 1][valid])

        return result