"""Spectral line shape functions."""
import numpy as np
from scipy.special import wofz

from .constants import Constants


class LineShapes:
    """Spectral line shape functions matching Fortran LineShapes module."""

    @staticmethod
    def voigt_humlicek(x: np.ndarray, a: float) -> np.ndarray:
        """Voigt function via Humlíček (1982) / Faddeeva function."""
        return np.real(wofz(x + 1j * a)) / Constants.sqrt_pi

    @staticmethod
    def voigt_normalized(x: np.ndarray,
                         gamma_L: float,
                         gamma_D: float) -> np.ndarray:
        """
        Voigt profile normalised so that its integral over all x is 1.

        Parameters
        ----------
        x       : offset from line centre [cm-1]
        gamma_L : Lorentz (pressure) half-width [cm-1]
        gamma_D : Doppler half-width [cm-1]
        """
        if gamma_D == 0:
            return LineShapes.lorentz(x, gamma_L)
        if gamma_L == 0:
            return LineShapes.doppler(x, gamma_D)

        x_reduced = x / gamma_D
        a = gamma_L / gamma_D
        K = LineShapes.voigt_humlicek(x_reduced, a)
        # voigt_humlicek already returns Re[w(z)]/sqrt(pi); dividing by
        # gamma_D alone yields a profile whose integral over all x equals 1.
        return K / gamma_D

    @staticmethod
    def lorentz(x: np.ndarray, gamma: float) -> np.ndarray:
        return gamma / (np.pi * (x ** 2 + gamma ** 2))

    @staticmethod
    def doppler(x: np.ndarray, gamma_D: float) -> np.ndarray:
        return (Constants.sqrt_ln2 / (gamma_D * Constants.sqrt_pi)) * \
               np.exp(-Constants.ln2_pi * (x / gamma_D) ** 2)