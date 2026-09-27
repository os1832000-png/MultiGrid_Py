"""A well-behaved library must not have import-time side effects."""
import subprocess
import sys


def _run(code):
    """Run python -c 'code' in a fresh subprocess, return (stdout, stderr)."""
    r = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True,
    )
    return r.stdout, r.stderr


def test_import_does_not_load_matplotlib():
    code = (
        "import sys\n"
        "import multigrid\n"
        "print('matplotlib' in sys.modules)\n"
    )
    out, err = _run(code)
    assert err == "", f"stderr: {err!r}"
    assert out.strip() == "False", \
        f"matplotlib was loaded on import; stdout={out!r}"


def test_import_adds_no_blanket_ignore_filter():
    code = (
        "import warnings, numpy, scipy\n"
        "from scipy.special import wofz\n"
        "before = list(warnings.filters)\n"
        "import multigrid\n"
        "added = [f for f in warnings.filters if f not in before]\n"
        "blanket = [f for f in added if f[0]=='ignore' and f[2] is Warning]\n"
        "print(len(blanket))\n"
    )
    out, err = _run(code)
    assert err == "", f"stderr: {err!r}"
    assert out.strip() == "0", \
        f"marfa added {out.strip()} blanket-ignore filter(s)"


def test_import_prints_nothing():
    code = "import multigrid\n"
    out, err = _run(code)
    assert out == "", f"stdout on import: {out!r}"
    assert err == "", f"stderr on import: {err!r}"
