"""Optional plotting helpers.

Requires matplotlib.  matplotlib is imported lazily inside each function
so that ``import multigrid`` never pulls it in.
"""
from typing import Optional

import numpy as np


def plot_spectrum(
    nu: np.ndarray,
    alpha: np.ndarray,
    *,
    title: Optional[str] = None,
    filename: Optional[str] = None,
    show: bool = True,
    ax=None,
):
    """
    Plot an absorption spectrum.

    Parameters
    ----------
    nu       : wavenumber grid [cm-1]
    alpha    : absorption coefficient [cm-1]
    title    : optional plot title
    filename : if given, save the figure to this path
    show     : if True, call ``plt.show()``
    ax       : optional existing matplotlib Axes to draw on
    """
    import matplotlib.pyplot as plt   # lazy — only when actually plotting

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(nu, alpha, "b-", linewidth=0.6)
    ax.set_xlabel("Wavenumber (cm⁻¹)")
    ax.set_ylabel("Absorption coefficient (cm⁻¹)")
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    if filename:
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
    if show:
        plt.show()
    return ax