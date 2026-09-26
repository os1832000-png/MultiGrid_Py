# ==============================================================================
# EXAMPLE
# ==============================================================================

def example_basic_calculation():
    print_marfa_banner()
    print("Example: CO2 Absorption Calculation")
    print("=" * 60)

    NU_MIN, NU_MAX = 0.0, 4000.0
    NU_RESOLUTION  = 0.05
    T, P           = 323.0, 1.0
    MOLE_FRACTION  = 400e-6

    print(f"\n  Range      : {NU_MIN}-{NU_MAX} cm-1")
    print(f"  Resolution : {NU_RESOLUTION} cm-1")
    print(f"  T={T} K, P={P} atm,  CO2={MOLE_FRACTION*1e6:.0f} ppm")

    nu_grid = np.arange(NU_MIN, NU_MAX + NU_RESOLUTION, NU_RESOLUTION)

    print("\nLoading spectral lines...")
    lines = HITRANReader.read_par_file(
        "CO2.par",
        molecule_id=2,
        wavenumber_min=NU_MIN,
        wavenumber_max=NU_MAX,
        intensity_threshold=1e-30)

    if not lines:
        print("No lines loaded -- check CO2.par")
        return

    marfa = MARFA()

    print("\nComputing absorption coefficient...")
    alpha = marfa.calculate_absorption_coefficient(
        nu_grid, lines, T, P, MOLE_FRACTION,
        line_cutoff_cm=25.0, wing_correction='none')

    print(f"\n  Max alpha : {alpha.max():.3e} cm-1")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(nu_grid, alpha, 'b-', linewidth=0.6)
    ax.set_xlabel('Wavenumber (cm-1)')
    ax.set_ylabel('Absorption coefficient (cm-1)')
    ax.set_title('MARFA v3: CO2 Absorption Coefficient')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig('marfa_python_example.png', dpi=300)
    print("\n  Plot saved: marfa_python_example.png")
    plt.show()


if __name__ == "__main__":
    example_basic_calculation()