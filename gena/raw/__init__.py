"""Raw packet decoding and raw-to-Parquet conversion."""

from importlib import import_module

__all__ = ['RawConfig', 'RunReport', 'convert_raw']


def __getattr__(name):
    if name in __all__:
        module = import_module(f'{__name__}.decode')
        return getattr(module, name)
    raise AttributeError(name)
