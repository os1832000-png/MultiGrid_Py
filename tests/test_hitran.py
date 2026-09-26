"""HITRAN reader tests using a synthetic line."""
import tempfile
import os
from marfa.hitran import HITRANReader, SpectralLine


def _make_par_line(mol=2, iso=1, nu=2000.0, S=1e-20,
                   gamma_air=0.07, gamma_self=0.1,
                   E_low=100.0, n_air=0.7, delta_air=-0.001):
    """Build a valid 160+ char HITRAN .par line."""
    parts = []
    parts.append(f"{mol:2d}")              # 0:2
    parts.append(f"{iso:d}")               # 2:3
    parts.append(f"{nu:12.6f}")            # 3:15
    parts.append(f"{S:10.3e}")             # 15:25
    parts.append(f"{0.0:10.3e}")           # 25:35
    parts.append(f"{gamma_air:5.3f}")      # 35:40
    parts.append(f"{gamma_self:5.3f}")     # 40:45
    parts.append(f"{E_low:10.4f}")         # 45:55
    parts.append(f"{n_air:4.2f}")          # 55:59
    parts.append(f"{delta_air:8.6f}")      # 59:67
    line = "".join(parts)
    return line.ljust(160)


def test_read_synthetic_line():
    content = _make_par_line(mol=2, nu=2000.0, S=1e-20)
    with tempfile.NamedTemporaryFile("w", suffix=".par", delete=False) as fh:
        fh.write(content + "\n")
        path = fh.name
    try:
        lines = HITRANReader.read_par_file(path, molecule_id=2)
        assert len(lines) == 1
        ln = lines[0]
        assert isinstance(ln, SpectralLine)
        assert ln.mol_id == 2
        assert abs(ln.nu - 2000.0) < 1e-4
        assert abs(ln.S - 1e-20) / 1e-20 < 1e-6
        assert abs(ln.gamma_air - 0.07) < 1e-4
    finally:
        os.unlink(path)


def test_missing_file_returns_empty():
    lines = HITRANReader.read_par_file("/nonexistent/does_not_exist.par")
    assert lines == []


def test_molecule_id_filter():
    a = _make_par_line(mol=2, nu=2000.0)
    b = _make_par_line(mol=1, nu=2000.0)
    with tempfile.NamedTemporaryFile("w", suffix=".par", delete=False) as fh:
        fh.write(a + "\n" + b + "\n")
        path = fh.name
    try:
        only_co2 = HITRANReader.read_par_file(path, molecule_id=2)
        assert len(only_co2) == 1
        assert only_co2[0].mol_id == 2
        only_h2o = HITRANReader.read_par_file(path, molecule_id=1)
        assert len(only_h2o) == 1
        assert only_h2o[0].mol_id == 1
    finally:
        os.unlink(path)
