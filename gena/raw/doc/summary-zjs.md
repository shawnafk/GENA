# ZJS Filter Code Analysis Summary

## 1. Overview

The ZJS Python codebase implements a multi-level filtering system for ENA (Energetic Neutral Atom) detector data. The filtering operates on two main data levels:
- **Level 1 (df1)**: Raw TDC (Time-to-Digital Converter) time data
- **Level 2 (df2)**: Converted spatial coordinates and physical parameters

## 2. Data Structure Definitions

### 2.1 Level 1 Data (df1) - Raw TDC Times
Location: `raw/structure.py`, lines 134-155

```python
df1_type = {
    'UUID':          np.dtype('S36'),
    'PID':           np.uint8,
    'BM':            np.uint8,
    'BX':            np.uint16,
    'Timestamp':     np.float64,
    'TimeInternal':  np.float32,
    'EvtLevel':      np.uint8,
    'EP':            np.float32,
    'Tof_0':         np.float32,
    'Tof_1':         np.float32,
    'Tof_2':         np.float32,
    'Tof_3':         np.float32,
    'TX_1s':         np.float32,  # X1 start time
    'TX_1e':         np.float32,  # X1 end time
    'TX_2s':         np.float32,  # X2 start time
    'TX_2e':         np.float32,  # X2 end time
    'TY_1s':         np.float32,  # Y1 start time
    'TY_1e':         np.float32,  # Y1 end time
    'TY_2s':         np.float32,  # Y2 start time
    'TY_2e':         np.float32,  # Y2 end time
}
```

### 2.2 Level 2 Data (df2) - Physical Coordinates
Location: `raw/structure.py`, lines 176-197

```python
df2_type = {
    'UUID':          np.dtype('S36'),
    'PID':           np.uint8,
    'BM':        np.uint8,
    'BX':        np.uint16,
    'Timestamp': np.float64,
    'TimeInternal':  np.float32,
    'EvtLevel':  np.uint8,
    'EP':        np.float32,
    'X1':        np.float32,  # Layer 1 X coordinate
    'X2':        np.float32,  # Layer 2 X coordinate
    'Y1':        np.float32,  # Layer 1 Y coordinate
    'Y2':        np.float32,  # Layer 2 Y coordinate
    'PX':        np.float32,  # Projected X angle
    'PY':        np.float32,  # Projected Y angle
    'Tof_M':     np.float32,
    'Tof_X':     np.float32,
    'Tof_Y':     np.float32,
    'Tof':       np.float32,  # Final time-of-flight
    'Spec':      np.uint8,    # Species (1=H, 2=O)
    'Energy':    np.float32,  # Energy in keV
}
```

## 3. Filter Functions - Level 1 (df1)

Location: `raw/filter.py`, lines 1-82

### 3.1 Signal Presence Filters

| Function | Line | Logic | Description |
|----------|------|-------|-------------|
| `triggered_eff(df)` | 4-6 | `sum != 0` | At least one TDC channel has signal |
| `start_eff(df)` | 9-11 | `TX_1s * TX_2s * TY_1s * TY_2s != 0` | All 4 start channels have signal |
| `start_stop_eff(df)` | 14-16 | Product of all 8 channels != 0 | All start AND end channels have signal |
| `X_start_eff(df)` | 19-21 | `TX_1s * TX_2s != 0` | X start channels have signal |
| `X_start_stop_eff(df)` | 24-26 | Product of 4 X channels != 0 | X start AND end have signal |
| `Y_start_eff(df)` | 29-31 | `TY_1s * TY_2s != 0` | Y start channels have signal |
| `Y_start_stop_eff(df)` | 34-36 | Product of 4 Y channels != 0 | Y start AND end have signal |
| `start_only_eff(df)` | 41-44 | `start_prod != 0 AND end_sum == 0` | Has start but no end signal |

### 3.2 Time Window Filters

| Function | Line | Logic | Description |
|----------|------|-------|-------------|
| `start_ns_eff(df, xlim=50, ylim=40)` | 47-48 | `TX_1s < xlim AND TX_2s < xlim AND TY_1s < ylim AND TY_2s < ylim` | Start signals within absolute limits |
| `start_ns_eff2(df)` | 51-57 | `abs(TX_1s - TX_2s) <= 55 AND abs(TY_1s - TY_2s - 10) <= 35` | Start signals within relative constraints |
| `match_x_eff(df, xl=45, xu=55)` | 64-66 | `(TX_1s + TX_2s) > xl AND < xu` | Sum of X start signals in range |
| `match_y_eff(df, yl=35, yu=45)` | 68-70 | `(TY_1s + TY_2s) > yl AND < yu` | Sum of Y start signals in range |
| `match_eff(df)` | 72-73 | `match_x_eff AND match_y_eff` | Both X and Y sums in range |

### 3.3 Coordinate Conversion

Location: `raw/filter.py`, lines 147-166

```python
def ns_to_xy(df, flag, prob, snid=1):
    """
    Convert TDC time differences to spatial coordinates.

    Args:
        df: DataFrame with TDC times
        flag: 0 for start times, 1 for end times
        prob: 'a' or 'b' for probe selection
        snid: 1 for orbit, 2 for ground

    Returns:
        x, y: Spatial coordinates in mm
    """
    # Uses calibration parameters from constants
    # X = ((TX_2s - TX_1s) - xc) * xcoef
    # Y = ((TY_2s - TY_1s) - yc) * ycoef
```

## 4. Filter Functions - Level 2 (df2)

Location: `raw/filter.py`, lines 91-142

### 4.1 Spatial Filters

| Function | Line | Logic | Description |
|----------|------|-------|-------------|
| `onboard_eff(df, center=[0,0], c_size=[90,30], m_size=[150,60])` | 92-105 | Box containment check | Events within detector boundaries |
| `select_circ(df, center=[0,0], radius=10, layer=1)` | 121-129 | `sqrt((x-cx)^2 + (y-cy)^2) < radius` | Circular region selection |
| `select_rect(df, center=[0,0], width=10, height=10, layer=1)` | 131-142 | Box containment | Rectangular region selection |

### 4.2 Physics Filters

| Function | Line | Logic | Description |
|----------|------|-------|-------------|
| `angle_eff(df, Ax, Ay)` | 107-114 | `abs(arctan((X1-X2)/50)*180/pi) < Ax` | Azimuth/elevation angle limits |
| `energy_eff(df, Emin, Emax)` | 118-119 | `Energy > Emin AND Energy < Emax` | Energy window filter |

## 5. Filter Composition

Location: `raw/filter.py`, lines 76-82

```python
from functools import reduce

def filter_df(df, flts):
    """
    Apply multiple filters to a DataFrame using logical AND.

    Args:
        df: Input DataFrame (df1 or df2)
        flts: List of filter functions

    Returns:
        Filtered ENA_DataFrame
    """
    combined_mask = reduce(lambda x, y: x & y, [flt(df) for flt in flts])
    return load.ENA_DataFrame(df.loc[combined_mask])
```

## 6. Data Conversion Pipeline

Location: `raw/conversion.py`

### 6.1 Level 0 → Level 1 (l0_2_l1)
Lines 28-68

- Converts raw channel values (`ch_tof_0`, `ch_1s`, etc.) to TDC times
- Applies probe-specific channel mapping (Probe A vs Probe B)
- Converts EP (Energy Pulse) values

### 6.2 Level 1 → Level 2 (l1_2_l2)
Lines 168-203

- Calculates TOF from TDC times using `f_tof_mxy()`
- Determines particle species using `ENA_spec()` (1=H, 2=O)
- Converts TDC times to XY coordinates using `ns_to_xy()`
- Calculates energy using `ENA_energy()`
- Computes projection angles PX, PY

## 7. Typical Filter Usage Patterns

### 7.1 Basic Level 1 Filtering
```python
from raw.filter import filter_df, start_eff, start_ns_eff2, match_eff

# Apply basic quality filters on df1
flts = [start_eff, start_ns_eff2, match_eff]
df1_filtered = filter_df(df1, flts)
```

### 7.2 Level 2 Geometric Filtering
```python
from raw.filter import filter_df, onboard_eff, angle_eff, energy_eff

# Apply geometric and physics filters on df2
flts = [
    lambda df: onboard_eff(df, center=[0,0], c_size=[90,30]),
    lambda df: angle_eff(df, Ax=30, Ay=15),
    lambda df: energy_eff(df, Emin=0.5, Emax=250)
]
df2_filtered = filter_df(df2, flts)
```

## 8. Key Design Patterns

1. **Functional Filter Design**: Each filter is a pure function that takes a DataFrame and returns a boolean mask
2. **Composability**: Filters can be combined using `filter_df()` with logical AND
3. **Configurable Parameters**: Most filters accept configurable limits (xlim, ylim, etc.)
4. **Probe-Agnostic**: `ns_to_xy()` handles Probe A and Probe B calibration via parameters
5. **Vectorized Operations**: All filters use NumPy vectorized operations for performance

## 9. Limitations and Notes

1. **No OR Composition**: The current `filter_df()` only supports AND composition
2. **No Built-in OR/NOT**: Filter functions don't provide negation or alternative combinations
3. **Hardcoded Constants**: Some filters have hardcoded limits that may need tuning
4. **Limited Documentation**: Some filter functions lack detailed docstrings
5. **Rollover Handling**: No TDC rollover correction in the filter functions (handled elsewhere)
