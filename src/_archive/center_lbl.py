# ------------------------------------------------------------------
    # centerLBL -- contributions from [0, deltaWV]
    # ------------------------------------------------------------------
    def centerLBL_full(self, FREQ: float, UL: float,
                       FSHAPE: Callable[[float], float], EPS: float):
        """
        Port of Fortran subroutine centerLBL.
        Handles the central part of the extended subinterval.
        """
        gs = self.gs
        UU = UL - FREQ

        if UU >= gs.deltaWV:
            return

        FF = float(FSHAPE(0.0))
        if FF < EPS:
            return

        NPOINT = 1

        # ---- left-of-centre sweep ------------------------------------
        FA = float(FSHAPE(UU))
        EPS4 = EPS * 0.25
        if FA > EPS4:
            gs.RK[1] += FA

        if UU < gs.H:
            # jump directly to right-of-centre sweep (label 211)
            pass
        else:
            I = 0
            UUU = UU

            # level 0
            if UUU >= gs.H0 + gs.H0:
                done0 = False
                for i in range(1, gs.NT0 + 1):
                    UUU -= gs.H0
                    FF = float(FSHAPE(UUU))
                    if FF < EPS:
                        break
                    gs.RK0P[i] += FA
                    gs.RK0[i]  += float(FSHAPE(UUU + gs.H1))
                    gs.RK0L[i] += FF
                    FA = FF
                    if UUU - gs.H0 < gs.H0:
                        done0 = True
                        break
                I = i * 2 if not done0 else i * 2

            # level 1
            if UUU >= gs.H0:
                IB = I + 1
                for i in range(IB, gs.NT1 + 1):
                    UUU -= gs.H1
                    FF = float(FSHAPE(UUU))
                    if FF < EPS:
                        break
                    gs.RK1P[i] += FA
                    gs.RK1[i]  += float(FSHAPE(UUU + gs.H2))
                    gs.RK1L[i] += FF
                    FA = FF
                    if UUU - gs.H1 < gs.H1:
                        break
                I = i * 2

            # level 2
            if UUU >= gs.H1:
                IB = I + 1
                for i in range(IB, gs.NT2 + 1):
                    UUU -= gs.H2
                    FF = float(FSHAPE(UUU))
                    if FF < EPS:
                        break
                    gs.RK2P[i] += FA
                    gs.RK2[i]  += float(FSHAPE(UUU + gs.H3))
                    gs.RK2L[i] += FF
                    FA = FF
                    if UUU - gs.H2 < gs.H2:
                        break
                I = i * 2

            # levels 3..9 follow same pattern
            level_data = [
                (gs.H2, gs.H3, gs.H4, gs.NT3, gs.RK3P, gs.RK3, gs.RK3L),
                (gs.H3, gs.H4, gs.H5, gs.NT4, gs.RK4P, gs.RK4, gs.RK4L),
                (gs.H4, gs.H5, gs.H6, gs.NT5, gs.RK5P, gs.RK5, gs.RK5L),
                (gs.H5, gs.H6, gs.H7, gs.NT6, gs.RK6P, gs.RK6, gs.RK6L),
                (gs.H6, gs.H7, gs.H8, gs.NT7, gs.RK7P, gs.RK7, gs.RK7L),
                (gs.H7, gs.H8, gs.H9, gs.NT8, gs.RK8P, gs.RK8, gs.RK8L),
                (gs.H8, gs.H9, gs.H + gs.H, gs.NT9,
                 gs.RK9P, gs.RK9, gs.RK9L),
            ]
            for (Hmin, Hn, Hn1, NTn, RKnP, RKn, RKnL) in level_data:
                if UUU >= Hmin:
                    IB = I + 1
                    for i in range(IB, NTn + 1):
                        UUU -= Hn
                        FF = float(FSHAPE(UUU))
                        if FF < EPS:
                            break
                        RKnP[i] += FA
                        RKn[i]  += float(FSHAPE(UUU + Hn1))
                        RKnL[i] += FF
                        FA = FF
                        if UUU - Hn < Hn:
                            break
                    I = i * 2

            # fill coarsest grid
            I = I * 4
            IB = I + 2
            CONSER = UU - (IB - 1) * gs.H
            for ICON in range(IB, gs.NT + 1):
                gs.RK[ICON] += float(FSHAPE(CONSER))
                CONSER -= gs.H
                if CONSER < 0.0:
                    NPOINT = ICON
                    break

        # ---- right-of-centre sweep (label 211) -----------------------
        NPOINT += 1
        UUU = gs.deltaWV - UU
        FA = float(FSHAPE(UUU))

        III = 0
        # level 0 (reverse direction)
        if UUU >= gs.H0 + gs.H0:
            for i in range(gs.NT0, 0, -1):
                III += 1
                UUU -= gs.H0
                FF = float(FSHAPE(UUU))
                if FF < EPS:
                    break
                gs.RK0L[i] += FA
                gs.RK0[i]  += float(FSHAPE(UUU + gs.H1))
                gs.RK0P[i] += FF
                FA = FF
                if UUU - gs.H0 < gs.H0:
                    break

        # level 1 reverse
        if UUU >= gs.H0:
            III = III * 2
            IB = gs.NT1 - III
            for i in range(IB, 0, -1):
                III += 1
                UUU -= gs.H1
                FF = float(FSHAPE(UUU))
                if FF < EPS:
                    break
                gs.RK1L[i] += FA
                gs.RK1[i]  += float(FSHAPE(UUU + gs.H2))
                gs.RK1P[i] += FF
                FA = FF
                if UUU - gs.H1 < gs.H1:
                    break

        # levels 2..9 reverse follow same pattern (abbreviated)
        rev_level_data = [
            (gs.H1, gs.H2, gs.H3, gs.NT2, gs.RK2P, gs.RK2, gs.RK2L),
            (gs.H2, gs.H3, gs.H4, gs.NT3, gs.RK3P, gs.RK3, gs.RK3L),
            (gs.H3, gs.H4, gs.H5, gs.NT4, gs.RK4P, gs.RK4, gs.RK4L),
            (gs.H4, gs.H5, gs.H6, gs.NT5, gs.RK5P, gs.RK5, gs.RK5L),
            (gs.H5, gs.H6, gs.H7, gs.NT6, gs.RK6P, gs.RK6, gs.RK6L),
            (gs.H6, gs.H7, gs.H8, gs.NT7, gs.RK7P, gs.RK7, gs.RK7L),
            (gs.H7, gs.H8, gs.H9, gs.NT8, gs.RK8P, gs.RK8, gs.RK8L),
            (gs.H8, gs.H9, gs.H + gs.H, gs.NT9,
             gs.RK9P, gs.RK9, gs.RK9L),
        ]
        for (Hmin, Hn, Hn1, NTn, RKnP, RKn, RKnL) in rev_level_data:
            if UUU >= Hmin:
                III = III * 2
                IB = NTn - III
                for i in range(IB, 0, -1):
                    III += 1
                    UUU -= Hn
                    FF = float(FSHAPE(UUU))
                    if FF < EPS:
                        break
                    RKnL[i] += FA
                    RKn[i]  += float(FSHAPE(UUU + Hn1))
                    RKnP[i] += FF
                    FA = FF
                    if UUU - Hn < Hn:
                        break

        # fill right side of coarsest grid
        III = III * 4
        I = gs.NT - III
        CONSER_r = gs.deltaWV - UU - (NPOINT - 1) * gs.H  # approximate
        for II in range(NPOINT, I + 1):
            if 1 <= II <= gs.NT:
                gs.RK[II] += float(FSHAPE(CONSER_r))
            CONSER_r -= gs.H
