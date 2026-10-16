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
