"""Sub-Lorentzian wing correction (chi-factor) functions."""
import numpy as np


class WingCorrections:
    """Sub-Lorentzian wing correction functions."""

    @staticmethod
    def tonkov_chi(delta_nu: np.ndarray, nu0: float) -> np.ndarray:
        abs_delta = np.abs(delta_nu)
        chi = np.ones_like(delta_nu, dtype=float)
        mask2 = (abs_delta > 25.0) & (abs_delta <= 250.0)
        chi[mask2] = 1.0 - 0.02 * ((abs_delta[mask2] - 25.0) / 225.0) ** 2
        mask3 = abs_delta > 250.0
        chi[mask3] = 0.98 * np.exp(-(abs_delta[mask3] - 250.0) / 200.0)
        return chi

    @staticmethod
    def perrin_chi(delta_nu: np.ndarray, nu0: float) -> np.ndarray:
        abs_delta = np.abs(delta_nu)
        chi = np.ones_like(delta_nu, dtype=float)
        mask = abs_delta > 20.0
        chi[mask] = np.exp(-0.0015 * (abs_delta[mask] - 20.0))
        return chi

    @staticmethod
    def no_correction(delta_nu: np.ndarray, nu0: float) -> np.ndarray:
        return np.ones_like(delta_nu, dtype=float)