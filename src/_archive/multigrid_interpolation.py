# ==============================================================================
# MULTI-GRID INTERPOLATION  (Fomin 1995)
# ==============================================================================

class MultiGridInterpolator:
    """
    Eleven-grid interpolation (Fomin 1995).
    Uses GridState + LineGridCalc to accumulate each line's contribution
    across the stencil hierarchy, then reconstructs the fine grid.
    """

    def __init__(self, n_grids: int = 11):
        self.n_grids = n_grids

    def create_grid_hierarchy(self, nu_min: float, nu_max: float,
                              delta_nu_fine: float) -> List[np.ndarray]:
        grids = []
        for i in range(self.n_grids):
            delta = delta_nu_fine * (2.0 ** i)
            grid = np.arange(nu_min, nu_max + delta, delta)
            grids.append(grid)
        return grids

    def accumulate_line(self,
                        gs: 'GridState',
                        lgc: 'LineGridCalc',
                        line_center: float,
                        subinterval_start: float,
                        FSHAPE: Callable[[float], float],
                        EPS: float = 1e-30):
        """
        Call leftLBL_full, centerLBL_full, rightLBL_full for one spectral line.

        Parameters
        ----------
        gs                : GridState  (holds the stencil arrays)
        lgc               : LineGridCalc instance
        line_center       : nu_0  [cm-1]   (Fortran FREQ)
        subinterval_start : UL    [cm-1]   (left edge of current subinterval)
        FSHAPE            : lineshape function  f(offset: float) -> float
        EPS               : convergence threshold
        """
        lgc.leftLBL_full(line_center, subinterval_start, FSHAPE, EPS)
        lgc.centerLBL_full(line_center, subinterval_start, FSHAPE, EPS)
        lgc.rightLBL_full(line_center, subinterval_start, FSHAPE, EPS)
