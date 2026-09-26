# ==============================================================================
# MULTI-GRID STATE  (Fortran grid arrays)
# ==============================================================================

class GridState:
    """
    Holds the multi-resolution stencil arrays that correspond to the
    Fortran module-level variables in LineGridCalc / Grids.

    For each level n in [0..9] there are three arrays:
        RKnP[i]  -- 'plus'   neighbor contribution at grid point i
        RKn[i]   -- 'center' contribution at grid point i
        RKnL[i]  -- 'left'   neighbor contribution at grid point i

    NT0..NT9 are the number of grid points at each level.
    RK[i] is the coarsest flat grid (NT points).

    Grid spacing hierarchy (Fortran convention):
        H  = finest spacing  (e.g. deltaWV / NT)
        H0 = 2*H, H1=2*H0, ..., H9=2*H8   (each level doubles)

    deltaWV = width of the central interval (endDeltaWV - startDeltaWV)
    cutOff  = wing cutoff beyond deltaWV
    """

    def __init__(self, delta_wv: float, cut_off: float, H: float):
        """
        Parameters
        ----------
        delta_wv : float
            Width of the central interval [cm-1]  (Fortran deltaWV)
        cut_off  : float
            Wing cutoff half-width [cm-1]           (Fortran cutOff)
        H        : float
            Finest grid spacing [cm-1]              (Fortran H)
        """
        self.deltaWV = delta_wv
        self.cutOff  = cut_off
        self.H       = H

        # Grid spacings (H0 = 2H, H1 = 4H, ..., H9 = 1024H)
        self.H0 = 2.0  * H
        self.H1 = 4.0  * H
        self.H2 = 8.0  * H
        self.H3 = 16.0 * H
        self.H4 = 32.0 * H
        self.H5 = 64.0 * H
        self.H6 = 128.0 * H
        self.H7 = 256.0 * H
        self.H8 = 512.0 * H
        self.H9 = 1024.0 * H

        # Number of grid points at each level (NT = total fine points in deltaWV)
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
        """Allocate all stencil arrays (initialised to zero)."""
        def z(n): return np.zeros(n + 2)   # +2 for 1-based index safety

        self.RK   = z(self.NT)

        self.RK0P = z(self.NT0);  self.RK0  = z(self.NT0);  self.RK0L = z(self.NT0)
        self.RK1P = z(self.NT1);  self.RK1  = z(self.NT1);  self.RK1L = z(self.NT1)
        self.RK2P = z(self.NT2);  self.RK2  = z(self.NT2);  self.RK2L = z(self.NT2)
        self.RK3P = z(self.NT3);  self.RK3  = z(self.NT3);  self.RK3L = z(self.NT3)
        self.RK4P = z(self.NT4);  self.RK4  = z(self.NT4);  self.RK4L = z(self.NT4)
        self.RK5P = z(self.NT5);  self.RK5  = z(self.NT5);  self.RK5L = z(self.NT5)
        self.RK6P = z(self.NT6);  self.RK6  = z(self.NT6);  self.RK6L = z(self.NT6)
        self.RK7P = z(self.NT7);  self.RK7  = z(self.NT7);  self.RK7L = z(self.NT7)
        self.RK8P = z(self.NT8);  self.RK8  = z(self.NT8);  self.RK8L = z(self.NT8)
        self.RK9P = z(self.NT9);  self.RK9  = z(self.NT9);  self.RK9L = z(self.NT9)

    def reset(self):
        """Zero all arrays (called before each new spectral line)."""
        self._allocate()

    def reconstruct_fine_grid(self) -> np.ndarray:
        """
        Reconstruct the full fine-grid absorption array from the stencil
        hierarchy. Fully vectorised -- no Python loops over grid points.

        Each coarse level n contributes its stencil values to the fine grid
        by computing all index positions as integer arrays and using
        np.add.at for scatter-accumulation.

        Returns
        -------
        np.ndarray  shape (NT,)  absorption on the fine grid inside [0, deltaWV]
        """
        result = self.RK[1:self.NT + 1].copy()

        levels = [
            (self.RK0,  self.RK0P,  self.RK0L,  self.NT0,  self.H0),
            (self.RK1,  self.RK1P,  self.RK1L,  self.NT1,  self.H1),
            (self.RK2,  self.RK2P,  self.RK2L,  self.NT2,  self.H2),
            (self.RK3,  self.RK3P,  self.RK3L,  self.NT3,  self.H3),
            (self.RK4,  self.RK4P,  self.RK4L,  self.NT4,  self.H4),
            (self.RK5,  self.RK5P,  self.RK5L,  self.NT5,  self.H5),
            (self.RK6,  self.RK6P,  self.RK6L,  self.NT6,  self.H6),
            (self.RK7,  self.RK7P,  self.RK7L,  self.NT7,  self.H7),
            (self.RK8,  self.RK8P,  self.RK8L,  self.NT8,  self.H8),
            (self.RK9,  self.RK9P,  self.RK9L,  self.NT9,  self.H9),
        ]
        half_step = int(round(1.0))   # stencil half-width in fine-grid indices = Hn/(2H) = Hn/H/2

        for RKn, RKnP, RKnL, NTn, Hn in levels:
            # i runs 1..NTn; centre fine-grid index = round(i*Hn/H) - 1
            ratio = int(round(Hn / self.H))
            i_arr = np.arange(1, NTn + 1)
            c_idx = i_arr * ratio - 1          # centre indices (0-based)
            l_idx = c_idx - ratio // 2         # left  neighbour
            r_idx = c_idx + ratio // 2         # right neighbour

            # clip to valid range [0, NT-1]
            NT = self.NT
            for idx_arr, RK_arr in ((c_idx, RKn), (l_idx, RKnL), (r_idx, RKnP)):
                valid = (idx_arr >= 0) & (idx_arr < NT)
                if np.any(valid):
                    np.add.at(result, idx_arr[valid], RK_arr[1:NTn+1][valid])

        return result

