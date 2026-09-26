# MultiGrid_Py

**Python port of the Fortran MARFA line-by-line molecular absorption code.**

Compute monochromatic absorption cross-sections and volume absorption
coefficients for planetary atmospheres, using HITRAN spectral line data.

Based on the original Fortran code by **Mikhail Razumovskiy**, **Boris Fomin**,
and **Denis Astanin**.

Reference: [MARFA: An Effective Line-by-line Tool for Calculating Molecular
Absorption in Planetary Atmospheres](https://arxiv.org/abs/2411.03418)
(arXiv:2411.03418, 2024).

---

## Features

- **Line-by-line (LBL)** monochromatic absorption from HITRAN `.par` files
- **Voigt, Lorentz, Doppler** line shapes (Voigt via Humlicek / Faddeeva)
- **Temperature-dependent intensities** using TIPS partition sums
- **Pressure and Doppler broadening** with HITRAN air/self parameters
- **Wing corrections** (Tonkov, Perrin) for sub-Lorentzian continua
- **Multi-resolution stencil machinery** (`GridState`, `LineGridCalc`)
  for future PT-table work
- Clean, tested, installable package -- 11 passing pytest tests
- No import-time side effects: lazy matplotlib, no global warnings mutation

---

## Install

Requires Python 3.6+, NumPy, and SciPy.

From source (development):

```bash
git clone https://github.com/os1832000-png/MultiGrid_Py.git
cd MultiGrid_Py
pip install --user -e .
With optional plotting support:

bash
pip install --user -e ".[plot]"
Quick start
python
import numpy as np
from marfa import MARFA, HITRANReader

# 1. Load CO2 lines from a HITRAN .par file
lines = HITRANReader.read_par_file(
    "CO2.par",
    molecule_id=2,
    wavenumber_min=2300.0,
    wavenumber_max=2400.0,
)
print(f"loaded {len(lines)} lines")

# 2. Build a wavenumber grid [cm^-1]
nu = np.arange(2300.0, 2400.0, 0.01)

# 3. Compute volume absorption coefficient [cm^-1]
alpha = MARFA().calculate_absorption_coefficient(
    nu, lines,
    T=296.0,             # temperature [K]
    P=1.0,               # pressure [atm]
    mole_fraction=400e-6,
    line_cutoff_cm=25.0,
    wing_correction="none",   # 'none' | 'tonkov' | 'perrin'
)
print(f"max alpha = {alpha.max():.4e} cm^-1")
For a runnable script see examples/co2_basic.py:

bash
python3 examples/co2_basic.py /path/to/CO2.par
Plotting (optional)
python
from marfa.plotting import plot_spectrum

plot_spectrum(
    nu, alpha,
    title="CO2 2300-2400 cm^-1",
    filename="co2_spectrum.png",
)
matplotlib is imported lazily, inside plot_spectrum. A plain
import marfa never loads matplotlib.

Command line
bash
python3 -m marfa.cli CO2.par --mol 2 --nu-min 2300 --nu-max 2400 \
    --T 296 --P 1.0 --vmr 400e-6 -o co2_spectrum.npz
The output .npz file contains nu and alpha arrays.

Get help with:

bash
python3 -m marfa.cli --help
Public API
Class / function	Purpose
MARFA	Main calculation engine
HITRANReader, SpectralLine	HITRAN .par file reader
TIPS	Total Internal Partition Sums
LineShapes	Voigt, Lorentz, Doppler profiles
WingCorrections	Sub-Lorentzian chi factors (Tonkov, Perrin)
GridState, LineGridCalc	Multi-resolution stencil machinery
AtmosphericProfile	Simple vertical profile container
PTTableGenerator	Pre-computed P-T absorption tables
plotting.plot_spectrum	Optional matplotlib helper
Tests
bash
pip install --user pytest
python3 -m pytest tests/ -v
11 tests covering line-shape normalisation, HITRAN parsing, the core
engine, and import-time safety guarantees.

Notes and limitations
TIPS partition sums use a simplified Q_ref * (T/296)^1.5 model,
not the full Gamache et al. (2017) tables. Adequate for qualitative work.

Wing corrections are simple parametrisations; not validated against
measured continua.

Line mixing is not implemented.

Continuum absorption (H2O, CO2, N2) is not implemented.

The multi-grid stencil classes (GridState, LineGridCalc) are present
but the working spectral path is the flat-grid Voigt sum.

License
MIT. See LICENSE.

Citation
If you use this code in a publication, please cite the original MARFA paper:

Razumovskiy, M., Fomin, B., Astanin, D.
MARFA: An Effective Line-by-line Tool for Calculating Molecular Absorption
in Planetary Atmospheres. arXiv:2411.03418 (2024).

Acknowledgements
This is an independent Python port based on the Fortran MARFA code by
Mikhail Razumovskiy, Boris Fomin, and Denis Astanin. All scientific credit
for the underlying algorithm belongs to them.
