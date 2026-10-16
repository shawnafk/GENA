"""
Unified check API with dimension expansion, diagnosis, and classification.

Three orthogonal dimensions for TDC channel checks:
  dir     : X, Y        (spatial axis)
  edge    : 1, 2        (detector edge)
  channel : s, e        (start/end timing)

Core API:
  check(kind, params, prob, dir, edge, channel, re_err, **kwargs)
      -> callable(df) -> boolean mask
  diagnose(df, params, prob, re_err)
      -> DataFrame of per-check booleans
  classify(df, params, prob, re_err)
      -> dict of named category masks
"""

import numpy as np
import pandas as pd
from itertools import product
from functools import reduce, wraps
import operator
from gena.constants import init_param
from gena import constants as gcon
TOL = 0.2
TOL_SUM = 0.2
TOL_DIFF = 0.2
class FilterNode:
    """Deferred boolean-mask calculation.

    Combine masks with ``&``, ``|`` and ``~``; apply them with frame.select().
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, df):
        # 真正执行计算的地方
        return self.func(df)

    def __and__(self, other):
        return FilterNode(lambda df: self.func(df) & other(df))

    def __or__(self, other):
        return FilterNode(lambda df: self.func(df) | other(df))
        
    def __invert__(self):
        return FilterNode(lambda df: ~self.func(df))


def filter_rule(func):
    """Accept df directly, or return a composable callable yielding a mask.

    Select rows with frame.select(rule); every rule returns a boolean mask.
    Positional and keyword df calls have the same behavior.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if (args and isinstance(args[0], pd.DataFrame)) or "df" in kwargs:
            return func(*args, **kwargs)
        return FilterNode(lambda df: func(df, *args, **kwargs))
    return wrapper


# ---------------------------------------------------------------------------
# Dimension definitions
# ---------------------------------------------------------------------------

DIMS = {
    'dir':     ['X', 'Y'],
    'edge':    ['1', '2'],
    'channel': ['s', 'e'],
}



# ---------------------------------------------------------------------------
# Atomic checks (private — use check() as public API)
# ---------------------------------------------------------------------------

def _length_atomic(df, dir, edge, channel, params=init_param(1), prob='a', re_err=TOL):
    """Single-channel TDC value within expected bound."""
    vmax = params[f'{prob.upper()}_{dir}_SUM'] * (1 + re_err)
    #print(f"check T{dir}_{edge}{channel}")    
    return df[f'T{dir}_{edge}{channel}'] < vmax

def _rel_atomic(df, dir, edge, params=init_param(1), prob='a', re_err=TOL):
    #print(f"check T{dir}_{edge}s"<f"check T{dir}_{edge}e")    
    return df[f'T{dir}_{edge}s'] <  df[f'T{dir}_{edge}e'] 


#no edge, always 1-2
def _sum_atomic(df, dir, channel, params=init_param(1), prob='a', re_err=TOL_SUM):
    """Start-channel sum conservation: T{dir}_1s + T{dir}_2s in range."""
    p = prob.upper()
    xu = params[f'{p}_{dir}_SUM'] * (1 + re_err)
    xl = params[f'{p}_{dir}_SUM'] / (1 + re_err)
    if isinstance(channel, (tuple, list)) and len(channel) == 2:
        c1, c2 = channel[0], channel[1]
    else:
        c1, c2 = channel, channel  # 如果只传了 's'，那就默认两边都是 's'
    #print(f"check T{dir}_1{c1} + T{dir}_2{c2}")    
    s = np.array(df[f'T{dir}_1{c1}'] + df[f'T{dir}_2{c2}'])
    return (s < xu) & (s > xl)


def _diff_atomic(df, dir, channel, params = init_param(1), prob='a', re_err=TOL_DIFF):
    """Difference T{dir}_1{ch} - T{dir}_2{ch} within calibrated boundaries."""
    xu = params[f'{prob.upper()}_{dir}R'] * (1 + re_err)
    #xl is negative, so still 1+re_err
    xl = params[f'{prob.upper()}_{dir}L'] * (1 + re_err)
    if isinstance(channel, (tuple, list)) and len(channel) == 2:
        c1, c2 = channel[0], channel[1]
    else:
        c1, c2 = channel, channel
    #print(f"check T{dir}_1{c1} - T{dir}_2{c2}")    
    d = np.array(df[f'T{dir}_1{c1}'] - df[f'T{dir}_2{c2}'])
    return (d < xu) & (d > xl)

# ---------------------------------------------------------------------------
# Unified check() factory
# ---------------------------------------------------------------------------
_DIM_CHECKS = {
    'length': _length_atomic,
    'sum':    _sum_atomic,
    'diff':   _diff_atomic,
    'rel':   _rel_atomic,
}

'''
def _expand(**kw):
    """Expand 'all' to concrete dimension values."""
    dirs = DIMS['dir'] if kw.get('dir') == 'all' else [kw['dir']]
    edges = DIMS['edge'] if kw.get('edge') == 'all' else [kw['edge']]
    channels = DIMS['channel'] if kw.get('channel') == 'all' else [kw['channel']]
    for d, e, c in product(dirs, edges, channels):
        yield d, e, c

def check(kind, dir='all', edge='all',
          channel='all', params=const.init_param(1), prob='a', re_err=0.1, **kwargs):
    """Return a callable(df) -> boolean mask.

    Dimension-aware checks (length, sum, diff):
      dir, edge, channel accept 'all' to expand and AND-combine.
      sum is start-channel-only (no channel expansion).
    Spatial checks (onboard, angle, tofxy, energy, circ, rect):
      kwargs passed through to the underlying function.
    """
    if kind not in _DIM_CHECKS:
        raise ValueError(
            f"Unknown check: {kind}. "
            f"Use: {sorted(_DIM_CHECKS)}"
        )

    atomic = _DIM_CHECKS[kind]
    def dim_fn(df):
        masks = []
        if kind == 'length':
            # length: varies dir x edge x channel
            for d, e, c in _expand(dir=dir, edge=edge, channel=channel):
                masks.append(atomic(df,  d, e, c, params, prob, re_err))
        elif kind == 'sum':
            # sum: varies dir only (start channel, both edges)
            dirs = DIMS['dir'] if dir == 'all' else [dir]
            for d, c in _expand(dir=dir, channel=channel):
                masks.append(atomic(df, d, c, params, prob, re_err))
        else: 
            # diff: varies dir x channel (across edges)
            dirs = DIMS['dir'] if dir == 'all' else [dir]
            channels = DIMS['channel'] if channel == 'all' else [channel]
            for d, c in product(dirs, channels):
                masks.append(atomic(df, d, c, params, prob, re_err))

        result = masks[0]
        for m in masks[1:]:
            result = result & m
        return result

    dim_fn.__name__ = f"{kind}_d{dir}_e{edge}_c{channel}"
    return FilterNode(dim_fn)
'''
def _expand(**kw):
    """
    动态展开传入的维度。
    例如: 
    _expand(dir='all', channel='s') -> 生成 ('X', 's'), ('Y', 's')
    _expand(dir='all', edge='all', channel='e') -> 生成 ('X', '1', 'e'), ('X', '2', 'e')...
    """
    iterables = []
    
    # 遍历你传进来的所有参数 (比如 dir='all', channel='s')
    for key, val in kw.items():
        if val == 'all':
            # 如果是 'all'，就从 DIMS 字典里取出所有的列表去遍历
            iterables.append(DIMS[key])
        else:
            # 如果指定了具体的值，就把它包装成单元素列表
            iterables.append([val])
            
    # 使用 *iterables 动态解包，生成所有维度的笛卡尔积
    for combo in product(*iterables):
        yield combo
        
def check(kind, dir='all', edge='all', channel='all', params=init_param(1), prob='a', re_err=None, **kwargs):
    if kind not in _DIM_CHECKS:
        raise ValueError(
            f"Unknown check: {kind}. "
            f"Use: {sorted(_DIM_CHECKS)}"
        )

    atomic = _DIM_CHECKS[kind]
    
    def dim_fn(df):
        masks = []
        
        if kind == 'length':
            tol = TOL if re_err is None else re_err
            for d, e, c in _expand(dir=dir, edge=edge, channel=channel):
                masks.append(atomic(df, dir=d, edge=e, channel=c, params=params, prob=prob, re_err=tol))
        elif kind == 'rel':
            tol = TOL if re_err is None else re_err
            for d, e in _expand(dir=dir, edge=edge):
                masks.append(atomic(df, dir=d, edge=e, params=params, prob=prob, re_err=tol))
                
        elif kind == 'sum':
            tol = TOL_SUM if re_err is None else re_err
            for d, c in _expand(dir=dir, channel=channel):
                #d_1c + d_2c
                masks.append(atomic(df, dir=d, channel=c, params=params, prob=prob, re_err=tol))
                
        else: 
            tol = TOL_DIFF if re_err is None else re_err
            # diff: varies dir x channel (2个维度)
            for d, c in _expand(dir=dir, channel=channel):
                #d_1c - d_2c
                masks.append(atomic(df, dir=d, channel=c, params=params, prob=prob, re_err=tol))

        # 合并所有的 mask
        result = masks[0]
        for m in masks[1:]:
            result = result & m
        return result

    dim_fn.__name__ = f"{kind}_d{dir}_e{edge}_c{channel}"
    return FilterNode(dim_fn)

# ---------------------------------------------------------------------------
# Trigger / efficiency checks (channel firing, not value-based)
# ---------------------------------------------------------------------------
@filter_rule
def triggered_eff(df):
    """Any of 8 TDC channels fired."""
    s = np.array(
        df['TX_1s'] + df['TX_2s'] + df['TY_1s'] + df['TY_2s']
        + df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e']
    )
    return (s != 0)


@filter_rule
def start_eff(df):
    """All 4 start channels fired."""
    m = np.array(df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'])
    return (m != 0)


@filter_rule
def start_only_eff(df):
    """All 4 start channels fired, no end channels."""
    m = np.array(df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s'])
    s = np.array(df['TX_1e'] + df['TX_2e'] + df['TY_1e'] + df['TY_2e'])
    return (m != 0) & (s == 0)


@filter_rule
def start_stop_eff(df):
    """All 8 channels fired."""
    m = np.array(
        df['TX_1s'] * df['TX_2s'] * df['TY_1s'] * df['TY_2s']
        * df['TX_1e'] * df['TX_2e'] * df['TY_1e'] * df['TY_2e']
    )
    return (m != 0)


# ---------------------------------------------------------------------------
# Spatial / geometric checks
# ---------------------------------------------------------------------------
#level2
@filter_rule
def onboard_eff(df, center=[0, 0], c_size=[90, 30], m_size=[150, 60]):
    """Event within carbon and MCP rectangular regions."""
    cxl = center[0] - c_size[0] / 2
    cxr = center[0] + c_size[0] / 2
    cyd = center[1] - c_size[1] / 2
    cyu = center[1] + c_size[1] / 2
    mxl = center[0] - m_size[0] / 2
    mxr = center[0] + m_size[0] / 2
    myd = center[1] - m_size[1] / 2
    myu = center[1] + m_size[1] / 2
    cond1 = (df['X1'] > cxl) & (df['Y1'] > cyd) & (df['X1'] < cxr) & (df['Y1'] < cyu)
    cond2 = (df['X2'] > mxl) & (df['Y2'] > myd) & (df['X2'] < mxr) & (df['Y2'] < myu)
    return cond1 & cond2


@filter_rule
def angle_eff(df, Ax, Ay):
    """Angle from position differences within bounds."""
    px = np.arctan((df['X1'] - df['X2']) / 50) * 180 / np.pi
    py = np.arctan((df['Y1'] - df['Y2']) / 50) * 180 / np.pi
    return (px > -Ax) & (px < Ax) & (py > -Ay) & (py < Ay)

@filter_rule
def block_eff(df, d, angle_uncert_deg=2.5):
    required_cols = ["X1", "Y1", "PX", "PY"]
    if not all(c in df.columns for c in required_cols):
        raise KeyError("缺失X1/Y1/PX/PY列")
    if d <= 0:
        raise ValueError("d必须>0")

    LX = gcon.X_CF
    LY = gcon.Y_CF
    alpha = np.deg2rad(gcon.THETA_X)
    beta = np.deg2rad(gcon.THETA_Y)

    # 窗口固定边界（你修正后的不对称偏移）
    x_left  = -LX / 2 - d * np.tan(alpha)
    x_right =  LX / 2 + d * np.tan(alpha)
    y_bot   = -LY / 2 - d * np.tan(beta)
    y_top   =  LY / 2 + d * np.tan(beta)

    # 原始角度 + 不确定度上下界（转弧度）
    px0 = np.deg2rad(df["PX"])
    py0 = np.deg2rad(df["PY"])
    d_ang = np.deg2rad(angle_uncert_deg)

    px_low  = px0 - d_ang
    px_high = px0 + d_ang
    py_low  = py0 - d_ang
    py_high = py0 + d_ang

    # 角度不确定带来的落点上下限
    x_hit_min = df["X1"] + d * np.tan(px_low)
    x_hit_max = df["X1"] + d * np.tan(px_high)
    y_hit_min = df["Y1"] + d * np.tan(py_low)
    y_hit_max = df["Y1"] + d * np.tan(py_high)

    # 区间重叠判定：粒子落点范围和窗口范围有交集即有效
    eff_x = (x_hit_min <= x_right) & (x_hit_max >= x_left)
    eff_y = (y_hit_min <= y_top) & (y_hit_max >= y_bot)

    return eff_x & eff_y
@filter_rule
def tofxy_eff(df, p_low=0.937, p_high=1.117):
    #print("check tofx ~ tofy")
    """TOF consistency between X and Y measurements."""
    tofx = ((df['TX_1e'] + df['TX_2e']) - (df['TX_1s'] + df['TX_2s'])) / 2
    tofy = ((df['TY_1e'] + df['TY_2e']) - (df['TY_1s'] + df['TY_2s'])) / 2
    r = tofx / tofy
    #!there are end less than starts
    return (r > p_low) & (r < p_high) & (tofx > 0) & (tofy > 0)


@filter_rule
def energy_eff(df, Emin, Emax):
    """Energy range filter."""
    return (df['Energy'] > Emin) & (df['Energy'] < Emax)


@filter_rule
def select_circ(df, center=[0, 0], radius=10, layer=1):
    """Circular region selection."""
    x, y = {1: (df['X1'], df['Y1']),
            2: (df['X2'], df['Y2']),
            3: (df['PX'], df['PY'])}[layer]
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    return (r < radius)


@filter_rule
def select_rect(df, center=[0, 0], width=10, height=10, layer=1):
    """Rectangular region selection."""
    x, y = {1: (df['X1'], df['Y1']),
            2: (df['X2'], df['Y2']),
            3: (df['PX'], df['PY'])}[layer]
    w, h = width / 2, height / 2
    return ((x > center[0] - w) & (x < center[0] + w)
            & (y > center[1] - h) & (y < center[1] + h))



# ---------------------------------------------------------------------------
# Standard filter factories
# ---------------------------------------------------------------------------
#check X Y, strictly
def golden_std(dir='all',params=init_param(1), prob='a', tofxy=True):
    base = [
        check('length', params=params, prob=prob, dir=dir, edge='all', channel='s'),
        check('rel', params=params, prob=prob, dir=dir,edge='all'),
        check('sum', params=params, prob=prob, dir=dir,channel='s'),
        check('diff', params=params, prob=prob, dir=dir, channel='s'),
        check('diff', params=params, prob=prob, dir=dir, channel='e'),
    ]
    if tofxy and dir == 'all':
        base.append(tofxy_eff)
    return reduce(operator.and_, base)

def soft_std(dir='all',params=init_param(1), prob='a'):
    return golden_std(dir=dir, params=params, prob=prob, tofxy=False)

# ---------------------------------------------------------------------------
# Diagnose & Classify
# ---------------------------------------------------------------------------

_DIAG_COLUMNS = [
    # X direction
    ('X_len_1s', 'length', {'dir': 'X', 'edge': '1', 'channel': 's'}),
    ('X_len_2s', 'length', {'dir': 'X', 'edge': '2', 'channel': 's'}),
    ('X_sum',    'sum',    {'dir': 'X'}),
    ('X_diff_s', 'diff',   {'dir': 'X', 'channel': 's'}),
    ('X_diff_e', 'diff',   {'dir': 'X', 'channel': 'e'}),
    # Y direction
    ('Y_len_1s', 'length', {'dir': 'Y', 'edge': '1', 'channel': 's'}),
    ('Y_len_2s', 'length', {'dir': 'Y', 'edge': '2', 'channel': 's'}),
    ('Y_sum',    'sum',    {'dir': 'Y'}),
    ('Y_diff_s', 'diff',   {'dir': 'Y', 'channel': 's'}),
    ('Y_diff_e', 'diff',   {'dir': 'Y', 'channel': 'e'}),
]


def diagnose(df, params, prob='a', re_err=0.1):
    """Return DataFrame of per-check boolean masks for repair decisions.

    Columns: X_len_1s, X_len_2s, X_sum, X_diff_s, X_diff_e,
             Y_len_1s, Y_len_2s, Y_sum, Y_diff_s, Y_diff_e, tofxy

    Example usage:
        d = diagnose(df, params)
        # X start fully trustable
        x_start_ok = d[['X_len_1s','X_len_2s','X_sum','X_diff_s']].all(axis=1)
        # X good, Y bad -> use X's tof to help repair Y
        can_repair_y_from_x = x_start_ok & ~d['Y_sum']
        # start ok but end diff bad -> use start + external tof for end
        can_repair_x_end = x_start_ok & ~d['X_diff_e']
    """
    result = {}
    for col_name, kind, dim_kw in _DIAG_COLUMNS:
        c = check(kind, params=params, prob=prob, re_err=re_err, **dim_kw)
        result[col_name] = c(df)
    result['tofxy'] = tofxy_eff(df)
    return pd.DataFrame(result, index=df.index)


def classify(df, params, prob='a', re_err=TOL):
    """Classify events by per-direction quality. Convenience over diagnose().

    Returns dict of named boolean masks (mutually exclusive, union = all rows):
      'golden'          - X ok, Y ok, tofxy ok
      'x_good_y_bad'    - X passes all, Y fails
      'x_bad_y_good'    - Y passes all, X fails
      'both_bad'        - neither direction passes
      'xy_good_tof_bad' - both ok, tofxy mismatch
    """
    d = diagnose(df, params, prob, re_err)
    x_start = d[['X_len_1s', 'X_len_2s', 'X_sum', 'X_diff_s']].all(axis=1)
    y_start = d[['Y_len_1s', 'Y_len_2s', 'Y_sum', 'Y_diff_s']].all(axis=1)
    x_ok = x_start & d['X_diff_e']
    y_ok = y_start & d['Y_diff_e']

    return {
        'golden':          x_ok & y_ok & d['tofxy'],
        'x_good_y_bad':    x_ok & ~y_ok,
        'x_bad_y_good':    ~x_ok & y_ok,
        'both_bad':        ~x_ok & ~y_ok,
        'xy_good_tof_bad': x_ok & y_ok & ~d['tofxy'],
    }
