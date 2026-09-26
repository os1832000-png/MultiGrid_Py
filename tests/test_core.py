"""Core MARFA engine sanity checks."""
import numpy as np
from marfa import MARFA
from marfa.hitran import SpectralLine


def _one_line(mol=2, nu=2000.0, S=1e-20,
              gamma_air=0.07, gamma_self=0.1,
              E_low=100.0, n_air=0.7, delta_air=0.0):
    return SpectralLine(
        mol_id=mol, iso_id=1, nu=nu, S=S, A=0.0,
        gamma_air=gamma_air, gamma_self=gamma_self,
        E_low=E_low, n_air=n_air, delta_air=delta_air,
    )


def test_empty_lines_returns_zeros():
    m = MARFA()
    nu = np.linspace(1990, 2010, 1001)
    sigma = m.calculate_absorption_cross_section(nu, [], 296.0, 1.0)
    assert sigma.shape == nu.shape
    assert np.all(sigma == 0.0)


def test_cross_section_integral_matches_intensity():
    m = MARFA()
    line = _one_line(S=1e-20, nu=2000.0, E_low=100.0)
    nu = np.linspace(1900.0, 2100.0, 400_001)
    sigma = m.calculate_absorption_cross_section(
        nu, [line], T=296.0, P=1.0, line_cutoff_cm=200.0,
    )
    integral = np.trapz(sigma, nu)
    assert integral > 0
    assert 0.9 * 1e-20 < integral < 1.1 * 1e-20, integral


def test_absorption_coefficient_scales_with_vmr():
    m = MARFA()
    line = _one_line()
    nu = np.linspace(1990.0, 2010.0, 20001)
    a1 = m.calculate_absorption_coefficient(nu, [line], 296.0, 1.0, 1e-4)
    a2 = m.calculate_absorption_coefficient(nu, [line], 296.0, 1.0, 2e-4)
    assert np.allclose(a2, 2.0 * a1)


def test_higher_temperature_widens_line():
    m = MARFA()
    line = _one_line()
    gd_cold = m.calculate_doppler_width(line.nu, 200.0, line.mol_id)
    gd_hot  = m.calculate_doppler_width(line.nu, 400.0, line.mol_id)
    assert gd_hot > gd_cold
    assert abs(gd_hot / gd_cold - np.sqrt(2.0)) < 1e-9
