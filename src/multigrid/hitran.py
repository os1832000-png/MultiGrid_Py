"""HITRAN spectral line reader."""
import logging
import os
from dataclasses import dataclass
from typing import List, Optional

from .molecules import MOLECULE_NAMES

logger = logging.getLogger(__name__)


@dataclass
class SpectralLine:
    """HITRAN spectral line parameters."""
    mol_id:       int
    iso_id:       int
    nu:           float
    S:            float
    A:            float
    gamma_air:    float
    gamma_self:   float
    E_low:        float
    n_air:        float
    delta_air:    float
    g_upper:      int = 0
    g_lower:      int = 0


class HITRANReader:
    """Read HITRAN ``.par`` files into :class:`SpectralLine` objects."""

    @staticmethod
    def read_par_file(
        filename: str,
        molecule_id: Optional[int] = None,
        wavenumber_min: Optional[float] = None,
        wavenumber_max: Optional[float] = None,
        intensity_threshold: Optional[float] = None,
    ) -> List[SpectralLine]:
        lines: List[SpectralLine] = []

        if not os.path.exists(filename):
            logger.warning("HITRAN file not found: %s", filename)
            return lines

        try:
            with open(filename, "r") as f:
                for lineno, line_str in enumerate(f, start=1):
                    if len(line_str) < 160:
                        continue
                    try:
                        mol        = int(line_str[0:2])
                        iso        = int(line_str[2])
                        nu         = float(line_str[3:15])
                        S          = float(line_str[15:25])
                        A          = float(line_str[25:35])
                        gamma_air  = float(line_str[35:40])
                        gamma_self = float(line_str[40:45])
                        E_low      = float(line_str[45:55])
                        n_air      = float(line_str[55:59])
                        delta_air  = float(line_str[59:67])
                    except (ValueError, IndexError):
                        logger.debug("Skipping malformed line %d in %s",
                                     lineno, filename)
                        continue

                    if molecule_id is not None and mol != molecule_id:
                        continue
                    if wavenumber_min is not None and nu < wavenumber_min:
                        continue
                    if wavenumber_max is not None and nu > wavenumber_max:
                        continue
                    if intensity_threshold is not None and S < intensity_threshold:
                        continue

                    lines.append(SpectralLine(
                        mol_id=mol, iso_id=iso, nu=nu, S=S, A=A,
                        gamma_air=gamma_air, gamma_self=gamma_self,
                        E_low=E_low, n_air=n_air, delta_air=delta_air,
                    ))
        except OSError:
            logger.exception("Error reading %s", filename)
            return lines

        if lines:
            name = MOLECULE_NAMES.get(molecule_id, "Unknown")
            logger.info("%s: %d lines loaded from %s", name, len(lines), filename)
        else:
            logger.warning("No lines matched the filters in %s", filename)

        return lines