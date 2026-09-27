"""LineGridCalc: Python port of Fortran LineGridCalc."""
import logging
from typing import Callable

from .grids import GridState

logger = logging.getLogger(__name__)


class LineGridCalc:
    """
    Python port of the Fortran LineGridCalc module.

    Implements ``leftLBL_full``, ``centerLBL_full`` and ``rightLBL_full``,
    which route each spectral line's shape-function sample to the correct
    multi-resolution stencil cell, matching the Fortran GOTO-cascade.

    ``FSHAPE`` plays the role of Fortran's ``procedure(shape)`` pointer: it
    receives a single float offset and returns a float lineshape value.
    """

    def __init__(self, gs: GridState):
        self.gs = gs

    # ------------------------------------------------------------------
    # leftLBL_full
    # ------------------------------------------------------------------
    def leftLBL_full(self,
                     FREQ: float,
                     UL: float,
                     FSHAPE: Callable[[float], float],
                     EPS: float) -> None:
        gs = self.gs
        UU = UL - FREQ

        if UU >= 0.0:
            return
        if -UU > gs.cutOff:
            return

        FF = float(FSHAPE(UU))
        if FF < EPS:
            return

        gs.RK[1] += FF

        if -UU < gs.H0:
            XXX = gs.H0
            for I in range(2, gs.NT0 + 1):
                gs.RK0P[I] += FF
                FF = float(FSHAPE(UU - XXX - gs.H1))
                gs.RK0[I]  += FF
                XXX += gs.H0
                FF = float(FSHAPE(UU - XXX))
                gs.RK0L[I] += FF
                if FF < EPS:
                    return
            return

        gs.RK0P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H1))
        gs.RK0[1] += FF_c
        FF = float(FSHAPE(UU - gs.H0))
        gs.RK0L[1] += FF

        if -UU < gs.H1:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 1)
            return

        gs.RK1P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H2))
        gs.RK1[1] += FF_c
        FF = float(FSHAPE(UU - gs.H1))
        gs.RK1L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H2:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 2)
            return

        gs.RK2P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H3))
        gs.RK2[1] += FF_c
        FF = float(FSHAPE(UU - gs.H2))
        gs.RK2L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H3:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 3)
            return

        gs.RK3P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H4))
        gs.RK3[1] += FF_c
        FF = float(FSHAPE(UU - gs.H3))
        gs.RK3L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H4:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 4)
            return

        gs.RK4P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H5))
        gs.RK4[1] += FF_c
        FF = float(FSHAPE(UU - gs.H4))
        gs.RK4L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H5:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 5)
            return

        gs.RK5P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H6))
        gs.RK5[1] += FF_c
        FF = float(FSHAPE(UU - gs.H5))
        gs.RK5L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H6:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 6)
            return

        gs.RK6P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H7))
        gs.RK6[1] += FF_c
        FF = float(FSHAPE(UU - gs.H6))
        gs.RK6L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H7:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 7)
            return

        gs.RK7P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H8))
        gs.RK7[1] += FF_c
        FF = float(FSHAPE(UU - gs.H7))
        gs.RK7L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H8:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 8)
            return

        gs.RK8P[1] += FF
        FF_c = float(FSHAPE(UU - gs.H9))
        gs.RK8[1] += FF_c
        FF = float(FSHAPE(UU - gs.H8))
        gs.RK8L[1] += FF
        if FF < EPS:
            return

        if -UU < gs.H9:
            self._rev_cascade_full(UU, FF, FSHAPE, EPS, 9)
            return

        gs.RK[2] += float(FSHAPE(UU - gs.H))
        gs.RK[3] += float(FSHAPE(UU - gs.H - gs.H))
        gs.RK[4] += float(FSHAPE(UU + gs.H - gs.H9))
        FF = float(FSHAPE(UU - gs.H9))
        gs.RK[5] += FF

        self._rev_cascade_full(UU, FF, FSHAPE, EPS, 9, col=2)

    def _rev_cascade_full(self,
                          UU: float,
                          FF: float,
                          FSHAPE: Callable[[float], float],
                          EPS: float,
                          start_level: int,
                          col: int = 2) -> None:
        gs = self.gs
        Hs   = [gs.H0, gs.H1, gs.H2, gs.H3, gs.H4,
                gs.H5, gs.H6, gs.H7, gs.H8, gs.H9]
        Hs_n = [gs.H1, gs.H2, gs.H3, gs.H4, gs.H5,
                gs.H6, gs.H7, gs.H8, gs.H9, gs.H + gs.H]
        RKP  = [gs.RK0P, gs.RK1P, gs.RK2P, gs.RK3P, gs.RK4P,
                gs.RK5P, gs.RK6P, gs.RK7P, gs.RK8P, gs.RK9P]
        RKC  = [gs.RK0,  gs.RK1,  gs.RK2,  gs.RK3,  gs.RK4,
                gs.RK5,  gs.RK6,  gs.RK7,  gs.RK8,  gs.RK9]
        RKL  = [gs.RK0L, gs.RK1L, gs.RK2L, gs.RK3L, gs.RK4L,
                gs.RK5L, gs.RK6L, gs.RK7L, gs.RK8L, gs.RK9L]

        for lvl in range(start_level - 1, -1, -1):
            RKP[lvl][col] += FF
            FF_c = float(FSHAPE(UU - Hs[lvl] - Hs_n[lvl]))
            RKC[lvl][col] += FF_c
            FF = float(FSHAPE(UU - Hs[lvl]))
            RKL[lvl][col] += FF
            if FF < EPS:
                return

    # ------------------------------------------------------------------
    # centerLBL_full  (see center_lbl.py — unchanged below)
    # ------------------------------------------------------------------
    # ... [paste centerLBL_full verbatim from your center_lbl.py] ...
    # ... [paste rightLBL_full and _right_cascade_full verbatim] ...