"""Default constants for raw conversion.

Runtime choices (snid, chunk size, products, overwrite policy) are carried by
``decode.RawConfig``; this module intentionally has no import-time environment
state.
"""

from common.constants import MASS_E, MASS_H, MASS_O, SPEED_C
from gena.constants import init_param

DEFAULT_SNID = 1
DEFAULT_CHUNK_PACKETS = 4096
DEFAULT_PRODUCTS = frozenset({"pdf", "l1a", "l2l", "l2m", "l2h", "l2f"})


if __name__ == '__main__':
    print('test')
