# `load.api`

This package is the source access boundary for the loading system. It reads remote or cached source products and returns source-native records together with explicit product layouts. DataFrame construction belongs to the loader layer.

`spedas.py` contains the common pySPEDAS access wrapper and `parse_spedas()`. `read_spedas()` exposes one public read signature while keeping the THEMIS, MMS, and OMNI backend call differences private. `read_spedas()` returns `(data, layouts)`, with each mission reader generating its layouts directly from the product information already available while reading. The public wrapper validates options, dispatches to the selected reader, and applies optional `variable_layouts` overrides. MMS data retain their complete backend variable names here; the shared loader preserves these complete names too. Callers pass the returned layouts to `parse_spedas()` without generating them again. The data records use the pySPEDAS fields `x`, `y`, and optional `v`, `v1`, `v2`, or `v3`.

`indices.py` provides Kyoto/GFZ source access and scalar parsing through `read_indices()` and `parse_indices()`; see [indices API details](indices.md). Product rules, vector declarations, and energy-axis names live in `load.api.spedas`. `OMNI_FIELDS` declares the HRO scalar layout and default variable list alongside `THEMIS_FIELDS` and `MMS_FIELDS`. Each mission reader applies its own rules before returning `(data, layouts)`. A private helper handles only common array representations and auxiliary axes; there is no separate mission-dispatching layout generator. DataFrame construction belongs in `load.loaders`. OMNI local CDF access is handled by its API reader, which returns the same records and layouts as remote access.

"""Parse observation sources data, assemble time-series frames, and expose loaders.

The loading process has three layers:

1. Parsers decode source-specific fields into explicit variable descriptions.
2. :func:`assemble` converts those descriptions into a :class:`common.frame.TimeSeriesDataFrame`.
3. Public ``load_*`` functions resolve files or remote backends, call the appropriate parser, and return the final frame.

Each variable description passed to :func:`assemble` is stored as ``name -> description`` and requires:

    ``time``
        Observation times accepted by :class:`common.time.Time`. Times must not
        contain ``NaT``.
    ``values``
        Values convertible to a regular NumPy array.
    ``time_axis``
        Integer identifying the time dimension of ``values``. The assembler moves
        this dimension to axis 0, whose length must equal the number of observation
        times.
    ``kind`` and ``columns``
        The output representation. A ``scalar`` has shape ``(time,)`` and one
        column. A ``vector`` has shape ``(time, component)`` and one column per
        component. An ``array`` has shape ``(time, ...)`` and one object column;
        each row stores the complete array sample for that time.

    A description may also contain ``units``, source ``attrs``, and ``axes``. Each
    axis description requires ``values``, an output ``column``, and a ``time_axis``.
    ``time_axis=None`` denotes a static coordinate shared by all observations; an
    integer denotes a time-dependent coordinate whose time dimension must match the
    variable. ``role="energy"`` records that axis as the variable's energy
    coordinate.

    ``Datetime`` and ``Timestamp`` are reserved output columns. Data columns must
    be unique across variables. Duplicate times within one variable keep their
    first sample. Variables may use different observation times; their columns are
    outer-joined on the union of times and sorted chronologically. Additional
    provenance is copied to ``frame.attrs`` together with the generated variable
    metadata.

GENA Parquet products are the exception to the parser/assembler path. Particle products contain independent event rows and parameter products already have a table layout, so :func:`load_gena` constructs their specialized frame classes directly.
"""

THEMIS, MMS, and OMNI use the shared loading pipeline in `load.loaders.spedas.load_spedas()`. The public functions in `load.load` apply the shared multi-interval decorator once and dispatch one normalized interval to each source loader. Only `load_spedas()` is exposed as the SPEDAS loader; mission-specific alignment and scientific helpers remain in the same module. GENA and indices retain separate source loaders.

The sole loading entrypoint is `load_spedas(mission, trange, product, ...)`. Mission-specific public wrappers and intermediate product-loading functions have been removed. All missions use `read_spedas()` followed by clipping, `parse_spedas()`, and `assemble()`. Dictionary output retains backend variable names and native axes. Default download policies are configured in `SPEDAS_OPTIONS` and resolved by the read API.

The indices path mirrors SPEDAS: `read_indices(trange, product, ...)` returns records and layouts, `parse_indices(data, layouts)` produces structured variables, and `load.loaders.indices.load_indices()` handles one normalized interval, assembly, and output selection. The public `load_indices()` wrapper handles multi-interval concatenation. `INDICES_FIELDS` declares Kyoto AE/Dst/SYM/ASY and GFZ Kp/ap/Ap/SN/F10.7 fields. All products share offline/local-first/refresh policies and monthly JSON caching.


# Indices API and loader

The indices path follows the same layers as SPEDAS: `read_indices()` returns `(data, layouts)`, `parse_indices(data, layouts)` returns structured variables, and `load_indices()` handles multiple intervals, assembly, concatenation, and the requested return type. Field definitions, URLs, network access, cache handling, fill values, units, and source metadata belong to `load/api/indices.py`; the loader does not implement source-specific downloads or parsers.

```python
from load import read_indices, parse_indices, assemble, load_indices

trange = ["2025-01-01", "2025-01-02"]
data, layouts = read_indices(trange, ["ae", "dst", "f107"])
variables = parse_indices(data, layouts)
frame = assemble(variables, source="indices")

# The complete loading path, including multiple intervals when needed.
frame = load_indices(trange, ["ae", "al", "au", "dst", "kp", "f107"])
```

## Fields

`INDICES_FIELDS` is the source of truth for product names, source datasets, remote fields, and units. Product names are case-insensitive. Omitting `product` loads `ae`, `al`, and `au`.

| Products | Source | Native cadence | Units |
| --- | --- | --- | --- |
| `ae`, `al`, `au`, `ao` | Kyoto HAPI `min_ae` | 1 minute | nT |
| `dst` | Kyoto HAPI `hour_dst` | 1 hour | nT |
| `sym_h`, `sym_d`, `asy_h`, `asy_d` | Kyoto HAPI `min_asysym` | 1 minute | nT |
| `kp` | GFZ `Kp` | 3 hours | Dimensionless |
| `ap` | GFZ `ap` | 3 hours | 2 nT per index unit |
| `ap_daily` | GFZ `Ap` | Daily | 2 nT per index unit |
| `sn` | GFZ `SN` | Daily | Dimensionless |
| `f107` | GFZ `Fobs` | Daily | sfu |
| `f107adj` | GFZ `Fadj` | Daily | sfu |

`f107` is the observed solar radio flux; `f107adj` is adjusted to 1 AU. `f107obs` and `f10.7` select `f107`; `f10.7adj` selects `f107adj`. Aliases return the canonical column name. `ap_daily` explicitly distinguishes daily Ap from three-hourly ap. Values retain the source's numerical scale: ap/Ap are not multiplied by two. One sfu is `10^-22 W m^-2 Hz^-1`.

Kyoto's best-available datasets select the highest available processing level for each period. Source timestamps are retained: Dst is stamped at `HH:29:30`, minute indices on the UTC minute, and GFZ daily products use UTC date labels at `00:00`. Daily labels are not asserted to be the physical observation time of the solar flux measurement. Different cadences are outer-joined without interpolation. Missing source values become `NaN`; network or format errors propagate to the caller instead of returning fabricated data.

## Cache and download policy

`prefix` selects the cache root; otherwise the configured `spedas_base` is used. Source JSON responses are stored under `indices/<source>/<dataset>/YYYYMM.json`. All requested AE components share the same monthly `min_ae` response; SYM/ASY components share `min_asysym`. Kyoto and GFZ requests are limited to individual UTC months, then clipped to the caller's inclusive interval. Month boundaries and overlapping requested intervals are deduplicated.

| `download` | Behavior |
| --- | --- |
| `False`, `"local"`, `"offline"` | Read existing cache only; missing cache raises `FileNotFoundError` |
| `True`, `"download"`, `"local_first"`, `"auto"` | Reuse existing cache and download missing months |
| `"refresh"` | Fetch each month again and replace the cache only after validation succeeds |

The default is `"local_first"`. Existing monthly files are snapshots, including the current month; use `"refresh"` when new samples or processing revisions are needed. Failed requests do not erase existing cache files. This JSON cache replaces the previous ad hoc Kyoto text/GFZ-table cache; `prefix` replaces the old `save_folder` option. The old `dst2` fallback name is removed: the Kyoto `dst` best-available dataset already selects the available processing version.

```python
cached = load_indices(trange, "dst", prefix="/data/space-weather", download=False)
updated = load_indices(trange, "dst", prefix="/data/space-weather", download="refresh")
```

Default output is `TimeSeriesDataFrame`. With `return_type="dict"`, a single selected product returns its `x`/`y` record, including units and source attributes; multiple products return a name-to-record mapping. Native times, units, source field names, citation identifiers, timestamp conventions, and month-specific processing metadata are retained in the structured variable attributes and resulting frame metadata.

## Official references

- [Kyoto HAPI documentation](https://wdc.kugi.kyoto-u.ac.jp/hapi/): datasets, best-available versions, timestamp conventions, and response formats.
- [GFZ data services](https://kp.gfz.de/en/data): JSON API and supported geomagnetic/solar indices.
- [GFZ Kp/ap description](https://doi.org/10.1029/2020SW002641): index definitions and ap/Ap scale.
- [NOAA F10.7 description](https://www.spaceweather.gov/phenomena/f107-cm-radio-emissions): solar radio flux units.


# THEMIS ESA ESD 读取说明

`load_themis()` 支持 THEMIS ESA 的独立 ESD（ESA 3-D distribution）L2 产品。ESD 是角分辨的三维速度空间观测：数据数组的三个维度依次为观测时间、能量通道和固角方向 bin。它与普通 `product="esa"` 返回的 `peir_en_eflux` 不同；后者只有时间和能量维，已经对方向积分，不能用于 pitch-angle distribution。

## 读取

```python
from load import load_themis

trange = ["2025-12-04T14:00:00", "2025-12-04T16:00:00"]
esd = load_themis(trange, "esd", probe="b", datatype="peif", download="local_first")
```

`product="esd"` 只支持 `level="l2"`。`datatype` 传给 pySPEDAS；常用值为 `peif`（full-mode ion）、`peef`（full-mode electron）、`peir`（reduced-mode ion）、`peer`（reduced-mode electron）和 `peib`（burst-mode ion）。默认值为 `peif`。和其他 THEMIS 读取器一样，默认 `download=False` 只读本地缓存；使用 `download=True`、`"local_first"` 或 `"refresh"` 才会访问远程数据。

默认读取五个 ESD CDF 原生变量：`valid`、`data_quality`、`bins`、`eff` 和 `eflux`。与常规 ESA 变量不同，这些名称不带 `tha_`、`thb_` 等探测器前缀；如果传入 `vars`，也必须使用这些裸名称，例如 `vars=["eflux", "bins"]`。

## 返回约定

默认返回 `TimeSeriesDataFrame`。每一行对应 ESD 的一个观测时间；`eflux`、`bins`、`eff` 单元格内均为形状 `(energy-channel, solid-angle-bin)` 的 NumPy 数组，因此整列拼接后的物理形状为 `(time, energy-channel, solid-angle-bin)`。pySPEDAS 的 `v1`、`v2` 坐标轴不会被丢弃，分别保存在 `eflux_v1`、`eflux_v2`（以及相同变量的 companion 列）中。

```python
import numpy as np

cube = np.stack(esd["eflux"].to_numpy())
bin_ok = np.stack(esd["bins"].to_numpy()).astype(bool)
quality_ok = (esd["valid"] == 1) & (esd["data_quality"] == 0)
cube[~quality_ok.to_numpy()] = np.nan
cube = np.where(bin_ok, cube, np.nan)
```

`eflux` 已是 ESD L2 校准后的微分能通量；不要再为了得到物理通量手动除以 `eff`。`bins == 1` 表示该能量–方向 bin 有效，`valid == 1` 和 `data_quality == 0` 是本项目图件使用的保守整条分布筛选。更细的质量位掩码解释应以相应 CDF 的元数据和 THEMIS ESA 数据质量说明为准。

若需要保留原始 pySPEDAS 字典而非 DataFrame，可使用：

```python
raw = load_themis(trange, "esd", probe="b", datatype="peif", return_type="dict")
eflux = raw["eflux"]
time = eflux["x"]
cube = eflux["y"]
energy_channel = eflux["v1"]
solid_angle_bin = eflux["v2"]
```

ESD 的 `v1` 和 `v2` 是 CDF 的通道索引，不是直接可用于 PA 的物理能量或方向。计算 pitch angle 时，还需读取同一 ESD CDF 的 `energy`、`theta`、`phi`、`domega` 及模式索引，并用常规 ESA L2 `peif_magf` 在 ESA DSL 坐标系给出磁场方向。`1.Paper_BBF/plot_themis_esd_pitch_angle_20251204.py` 展示了完整的 ESD-only PA 分布读取和绘图流程。
