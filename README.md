# MARFA

Python port of the Fortran MARFA line-by-line molecular absorption code.

## Install
    pip install --user -e .

## Quick start

    import numpy as np
    from marfa import MARFA, HITRANReader

    lines = HITRANReader.read_par_file("CO2.par", molecule_id=2)
    nu = np.arange(2300.0, 2400.0, 0.01)
    alpha = MARFA().calculate_absorption_coefficient(
        nu, lines, T=296.0, P=1.0, mole_fraction=400e-6)

Reference: arXiv:2411.03418 (2024)
