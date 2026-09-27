"""Command-line interface for MARFA."""
import argparse
import logging
import sys

import numpy as np

from .core import MARFA
from .hitran import HITRANReader


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="marfa",
        description="MARFA: LBL molecular absorption calculator",
    )
    parser.add_argument("par_file", help="HITRAN .par file")
    parser.add_argument("--mol", type=int, required=True,
                        help="HITRAN molecule ID (1-12)")
    parser.add_argument("--nu-min", type=float, default=0.0)
    parser.add_argument("--nu-max", type=float, default=4000.0)
    parser.add_argument("--dnu", type=float, default=0.05)
    parser.add_argument("--T", type=float, default=296.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument("--vmr", type=float, default=400e-6)
    parser.add_argument("--cutoff", type=float, default=25.0)
    parser.add_argument("--wing",
                        choices=["none", "tonkov", "perrin"],
                        default="none")
    parser.add_argument("-o", "--output", default=None,
                        help="Save spectrum as .npz")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug logging")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    lines = HITRANReader.read_par_file(
        args.par_file,
        molecule_id=args.mol,
        wavenumber_min=args.nu_min,
        wavenumber_max=args.nu_max,
    )
    if not lines:
        print(f"No lines loaded from {args.par_file}", file=sys.stderr)
        return 1

    nu = np.arange(args.nu_min, args.nu_max + args.dnu, args.dnu)
    marfa = MARFA()

    print(f"Computing alpha(nu) over {len(nu)} points "
          f"(T={args.T} K, P={args.P} atm, VMR={args.vmr:.3e}) ...")
    alpha = marfa.calculate_absorption_coefficient(
        nu, lines, args.T, args.P, args.vmr,
        line_cutoff_cm=args.cutoff,
        wing_correction=args.wing,
    )

    print(f"Max alpha = {alpha.max():.4e} cm^-1")

    if args.output:
        np.savez(args.output, nu=nu, alpha=alpha)
        print(f"Saved: {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())