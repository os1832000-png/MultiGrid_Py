# ------------------------------------------------------------------
    # rightLBL -- contributions from [deltaWV, deltaWV+cutOff]
    # ------------------------------------------------------------------
    def rightLBL_full(self, FREQ: float, UL: float,
                      FSHAPE: Callable[[float], float], EPS: float):
        """
        Port of Fortran subroutine rightLBL.
        Mirror of leftLBL but at the right edge of the subinterval.
        """
        gs = self.gs
        UU = UL - FREQ - gs.deltaWV    # offset beyond deltaWV

        if UU >= gs.cutOff:
            return

        FF = float(FSHAPE(UU))
        if FF < EPS:
            return

        NT = gs.NT

        if UU < gs.H0:
            gs.RK0L[gs.NT0] += FF
            FF = float(FSHAPE(UU + gs.H1))
            gs.RK0[gs.NT0]  += FF
            FF = float(FSHAPE(UU + gs.H0))
            gs.RK0P[gs.NT0] += FF
            # fall through to fine loop (label 12)
            XXX = gs.H0
            for I in range(gs.NT0 - 1, 0, -1):
                gs.RK0L[I] += FF
                FF = float(FSHAPE(UU + XXX + gs.H1))
                gs.RK0[I]  += FF
                XXX += gs.H0
                FF = float(FSHAPE(UU + XXX))
                gs.RK0P[I] += FF
                if FF < EPS:
                    return
            gs.RK[1] += FF
            return

        if UU < gs.H1:
            gs.RK1L[gs.NT1] += FF
            FF = float(FSHAPE(UU + gs.H2))
            gs.RK1[gs.NT1]  += FF
            FF = float(FSHAPE(UU + gs.H1))
            gs.RK1P[gs.NT1] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 1)
            return

        if UU < gs.H2:
            gs.RK2L[gs.NT2] += FF
            FF = float(FSHAPE(UU + gs.H3))
            gs.RK2[gs.NT2]  += FF
            FF = float(FSHAPE(UU + gs.H2))
            gs.RK2P[gs.NT2] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 2)
            return

        if UU < gs.H3:
            gs.RK3L[gs.NT3] += FF
            FF = float(FSHAPE(UU + gs.H4))
            gs.RK3[gs.NT3]  += FF
            FF = float(FSHAPE(UU + gs.H3))
            gs.RK3P[gs.NT3] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 3)
            return

        if UU < gs.H4:
            gs.RK4L[gs.NT4] += FF
            FF = float(FSHAPE(UU + gs.H5))
            gs.RK4[gs.NT4]  += FF
            FF = float(FSHAPE(UU + gs.H4))
            gs.RK4P[gs.NT4] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 4)
            return

        if UU < gs.H5:
            gs.RK5L[gs.NT5] += FF
            FF = float(FSHAPE(UU + gs.H6))
            gs.RK5[gs.NT5]  += FF
            FF = float(FSHAPE(UU + gs.H5))
            gs.RK5P[gs.NT5] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 5)
            return

        if UU < gs.H6:
            gs.RK6L[gs.NT6] += FF
            FF = float(FSHAPE(UU + gs.H7))
            gs.RK6[gs.NT6]  += FF
            FF = float(FSHAPE(UU + gs.H6))
            gs.RK6P[gs.NT6] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 6)
            return

        if UU < gs.H7:
            gs.RK7L[gs.NT7] += FF
            FF = float(FSHAPE(UU + gs.H8))
            gs.RK7[gs.NT7]  += FF
            FF = float(FSHAPE(UU + gs.H7))
            gs.RK7P[gs.NT7] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 7)
            return

        if UU < gs.H8:
            gs.RK8L[gs.NT8] += FF
            FF = float(FSHAPE(UU + gs.H9))
            gs.RK8[gs.NT8]  += FF
            FF = float(FSHAPE(UU + gs.H8))
            gs.RK8P[gs.NT8] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 8)
            return

        if UU < gs.H9:
            gs.RK9L[gs.NT9] += FF
            FF = float(FSHAPE(UU + gs.H + gs.H))
            gs.RK9[gs.NT9]  += FF
            FF = float(FSHAPE(UU + gs.H9))
            gs.RK9P[gs.NT9] += FF
            if FF < EPS:
                return
            self._right_cascade_full(UU, FF, FSHAPE, EPS, 9)
            return

        # label 69: beyond all fine grids
        gs.RK[NT]     += float(FSHAPE(UU))
        gs.RK[NT - 1] += float(FSHAPE(UU + gs.H))
        gs.RK[NT - 2] += float(FSHAPE(UU + gs.H + gs.H))
        gs.RK[NT - 3] += float(FSHAPE(UU + gs.H9 - gs.H))

        FF = float(FSHAPE(UU + gs.H9))
        self._right_cascade_full(UU, FF, FSHAPE, EPS, 9)

    def _right_cascade_full(self, UU: float, FF: float,
                             FSHAPE: Callable[[float], float], EPS: float,
                             start_level: int):
        """
        Right-side cascade (Fortran labels 139..131).
        Works from start_level down to level 1 filling NT-1 columns.
        """
        gs = self.gs
        Hs   = [gs.H0, gs.H1, gs.H2, gs.H3, gs.H4,
                gs.H5, gs.H6, gs.H7, gs.H8, gs.H9]
        Hs_n = [gs.H1, gs.H2, gs.H3, gs.H4, gs.H5,
                gs.H6, gs.H7, gs.H8, gs.H9, gs.H + gs.H]
        NTs  = [gs.NT0, gs.NT1, gs.NT2, gs.NT3, gs.NT4,
                gs.NT5, gs.NT6, gs.NT7, gs.NT8, gs.NT9]
        RKP  = [gs.RK0P, gs.RK1P, gs.RK2P, gs.RK3P, gs.RK4P,
                gs.RK5P, gs.RK6P, gs.RK7P, gs.RK8P, gs.RK9P]
        RKC  = [gs.RK0,  gs.RK1,  gs.RK2,  gs.RK3,  gs.RK4,
                gs.RK5,  gs.RK6,  gs.RK7,  gs.RK8,  gs.RK9]
        RKL  = [gs.RK0L, gs.RK1L, gs.RK2L, gs.RK3L, gs.RK4L,
                gs.RK5L, gs.RK6L, gs.RK7L, gs.RK8L, gs.RK9L]

        for lvl in range(start_level - 1, -1, -1):
            N = NTs[lvl] - 1
            RKL[lvl][N] += FF
            FF_c = float(FSHAPE(UU + Hs[lvl] + Hs_n[lvl]))
            RKC[lvl][N] += FF_c
            FF = float(FSHAPE(UU + Hs[lvl]))
            RKP[lvl][N] += FF
            if FF < EPS:
                return

        # label 12: fill finest RK0 loop going downward
        XXX = gs.H0
        for I in range(gs.NT0 - 1, 0, -1):
            gs.RK0L[I] += FF
            FF = float(FSHAPE(UU + XXX + gs.H1))
            gs.RK0[I]  += FF
            XXX += gs.H0
            FF = float(FSHAPE(UU + XXX))
            gs.RK0P[I] += FF
            if FF < EPS:
                return
        gs.RK[1] += FF

