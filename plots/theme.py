"""Repository-wide Matplotlib defaults applied explicitly by a caller.

Importing this module must not change global plotting state. Applications,
notebooks, command-line entry points and rendering workers should call
``apply_theme()`` once before creating figures.

Usage
-----
Pass any Matplotlib ``rcParams`` entry as a mapping; keys keep their dotted
Matplotlib names::

    from plots import theme

    theme.apply_theme({
        "font.size": 16,
        "axes.labelsize": 14,
        "axes.grid": True,
        "legend.frameon": False,
        "figure.figsize": (10, 8),
        "savefig.dpi": 300,
    })

``apply_theme`` merges the defaults below with the supplied mapping. It does
not undo keys that were set earlier. Use ``matplotlib.rc_context`` for a
single figure, or ``matplotlib.rcdefaults()`` followed by ``apply_theme()`` to
return to the repository defaults.

Explicit ``fontsize=``, ``figsize=`` and similar arguments inside a plotting
function still win over ``rcParams``; the theme only controls values that the
function leaves unspecified.
"""

from __future__ import annotations

from collections.abc import Mapping

import matplotlib as mpl

__all__ = ["DEFAULT_RC", "apply_theme"]

# Repository-wide defaults. Keep this list small: product code owns layout,
# while the shared theme only fixes typography and legend chrome. Any other
# valid Matplotlib rcParam can be passed directly to apply_theme() without
# being added here.
DEFAULT_RC = {
    "font.size": 14,
    "legend.frameon": False,
}


def apply_theme(overrides: Mapping[str, object] | None = None) -> dict:
    """Apply the shared plotting defaults and return the resulting rcParams.

    Parameters
    ----------
    overrides : mapping, optional
        Additional ``matplotlib.rcParams`` entries applied after the defaults.
        Keys use Matplotlib's dotted names, for example ``"font.size"`` or
        ``"axes.grid"``. Values are passed through unchanged.

    Returns
    -------
    dict
        A snapshot of the updated ``matplotlib.rcParams``.

    Examples
    --------
    >>> apply_theme({"font.size": 16, "axes.grid": True})
    {...}

    To scope a change to one figure instead of the whole process::

        import matplotlib as mpl
        with mpl.rc_context({"font.size": 12}):
            fig, ax = plot_product(...)

    Notes
    -----
    This function mutates process-wide Matplotlib state. Call it once at the
    start of a rendering worker or notebook, not concurrently per request.
    """
    values = dict(DEFAULT_RC)
    if overrides is not None:
        values.update(dict(overrides))
    mpl.rcParams.update(values)
    return dict(mpl.rcParams)
