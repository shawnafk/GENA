"""Lazy coordinate backend dispatch; all adapters exchange NumPy arrays."""

from importlib import import_module

BACKENDS = ("spacepy", "astropy", "pyspedas")


def validate_backend(backend):
    """Return a normalized backend name or raise for an unknown choice."""
    name = str(backend).lower()
    if name not in BACKENDS:
        raise ValueError(f"backend must be one of {BACKENDS}, got {backend!r}")
    return name


def transform(values, times, *, source, target, unit, origin, backend):
    """Dispatch an explicit backend without fallback or retained backend objects."""
    name = validate_backend(backend)
    adapter = import_module(f"{__name__}.{name}")
    try:
        return adapter.transform(values, times, source=source, target=target,
                                 unit=unit, origin=origin)
    except ModuleNotFoundError as error:
        if error.name and error.name.split('.')[0] == name:
            raise ImportError(f"backend {name!r} requires the {name} package") from error
        raise
