"""Tests that line shape integrals equal 1."""
import numpy as np
from marfa.line_shapes import LineShapes


def test_voigt_integral_is_one():
    gamma_L, gamma_D = 0.05, 0.03
    x = np.linspace(-2000, 2000, 2_000_001)
    y = LineShapes.voigt_normalized(x, gamma_L, gamma_D)
    assert abs(np.trapz(y, x) - 1.0) < 1e-3
