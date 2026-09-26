"""Example: CO2 absorption calculation with MARFA.

Run with:  python -m marfa.examples.co2_basic
Requires a CO2.par file in the current directory.
"""
import numpy as np
import matplotlib.pyplot as plt

from marfa import MARFA, HITRANReader
from marfa.plotting import plot_spectrum


def main():
    NU_MIN, NU_MAX = 0.0, 4000.0
    NU_RESOLUTION  = 0.05
    T, P           = 323.0, 1.0
    MOLE_FRACTION  = 400e-6

    print("=" * 60)
    print("MARFA v3: CO2 absorption example")
    print("=" * 60)
    print(f"  Range      : {NU_MIN}-{NU_MAX} cm-1")
    print(f"  Resolution : {NU_RESOLUTION} cm-1")
    print(f"  T={T} K, P={P} atm, CO2={MOLE_FRACTION * 1e6:.0f} ppm")

    nu_grid = np.arange(NU_MIN, NU_MAX + NU_RESOLUTION, NU_RESOLUTION)

    print("\nLoading spectral lines ...")
    lines = HITRANReader.read_par_file(
        "CO2.par",
        molecule_id=2,
        wavenumber_min=NU_MIN,
        wavenumber_max=NU_MAX,
        intensity_threshold=1e-30,
    )
    if not lines:
        print("No lines loaded -- check CO2.par")
        return

    marfa = MARFA()
    print("\nComputing absorption coefficient ...")
    alpha = marfa.calculate_absorption_coefficient(
        nu_grid, lines, T, P, MOLE_FRACTION,
        line_cutoff_cm=25.0, wing_correction="none",
    )
    print(f"\n  Max alpha : {alpha.max():.3e} cm-1")

    plot_spectrum(
        nu_grid, alpha,
        title="MARFA v3: CO2 Absorption Coefficient",
        filename="marfa_python_example.png",
        show=True,
    )
    print("\n  Plot saved: marfa_python_example.png")


if __name__ == "__main__":
    main()