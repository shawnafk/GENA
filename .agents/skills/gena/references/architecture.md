# GENA-CODE-Dev architecture and API routing

For detailed foundation-module behavior, array/time formats, dependencies, and runnable examples, read [common/README.md](../../../../common/README.md). Inspect the live signatures before calling an API.

## Contents

1. Dependency architecture
2. Data contracts
3. Function routing map
4. Loading and configuration
5. GENA processing workflow
6. Cross-mission workflow
7. Important distinctions and hazards

## 1. Dependency architecture

Work from the repository root containing `gena/` and `common/`.

```text
analysis scripts / gena.post.export
    +-- load: gena, themis, mms, omni, indices
    +-- gena: frame, constants, filter, raw, post, utils
    +-- plots: vis, basic, advance/{joint,ani,v3d}
    +-- libs: regions, plasma, spectrum, flux, ephemeris, environment, space
    +-- common: config, frame, time, coord, constants
```

Keep mission modules independent. Put time, coordinates, types and configuration in `common/`, shared plotting in `plots.vis`, and scientific models in `libs/`;
compose preloaded products in `plots.advance` and study-specific analyses in analysis scripts. Mission `vis.py` modules have been removed; product plots live in `plots.basic`, and report orchestration may call `plots`.

The `load` package exports only `load`, `load_pos`, and `load_sw`; the facade lives in `load/load.py`, and mission readers plus `read_*`/`parse_*`/`assemble` are module-level internals. THEMIS, MMS, and OMNI share `load/loaders/spedas.py`, which owns clipping, assembly, and return-type handling; mission-specific loader wrappers have been removed. GENA and indices retain their own loader modules; the shared assembler currently lives in `load/loaders/gena.py`. `libs.regions` owns cross-mission classification; explicit column mappings allow non-THEMIS inputs. Only `load`, `load_pos`, and `load_sw` are re-exported from the package; removed scientific helpers are not re-exported. Loader return types, download policies, configuration and units are unchanged.

Shared loader-only normalization helpers live privately in `load/api/variable.py`. In particular, `_clip_variable_dict()` trims backend `x`/`y` dictionaries before conversion and is not part of the `common` or `load` public API.

Treat `deprecated/` and `gena/raw/filter_od.py` as historical unless the user
explicitly requests legacy behavior. Prefer public functions without a leading
underscore. Inspect private helpers only when diagnosing their public caller.

## 2. Data contracts

### GENA tables

`load("gena", )` returns `gena.frame.ParameterDataFrame` for the PDF engineering product and `gena.frame.ParticleDataFrame` for particle products. `Timestamp` is Unix seconds; `load("gena", )` creates `Datetime` directly with pandas before construction. Both types inherit the shared `common.frame.TimeIndexedDataFrame` behavior, so `Datetime` is a timezone-naive pandas `datetime64[ns]` column whose values are defined as UTC. `ParameterDataFrame` additionally inherits `TimeSeriesDataFrame`; `ParticleDataFrame` does not. Time-series consumers access the normalized time directly with `frame["Datetime"]`; source-field mapping belongs in loaders. `TimeIndexedDataFrame.__init__` requires exactly one `Datetime` column and delegates normalization to `common.time.Time`, which rejects invalid/NaT values. It never infers a missing `Datetime` from `Timestamp` or `Date` and does not automatically add `Timestamp`. Pandas column selections without `Datetime` return ordinary DataFrames; results retaining it preserve the concrete subclass. This construction check does not continuously validate in-place mutations. Plain DataFrames passed to time-series operations must also provide `Datetime`.

Typical product columns:

- PDF/engineering: `Timestamp`, `MJ2000_Star_X/Y/Z`, `MJ2000_q1/q2/q3`.
- L1/L1A: TDC channels such as `TX_1s`, `TX_2s`, `TY_1e`, plus event metadata.
- L2: `Energy`, `Spec`, `X1/Y1`, `X2/Y2`, `PX/PY`, pose columns.
- Statistical weights are external arrays aligned positionally with event rows; science frames are not expected to store `Eff` or `W` columns.
- `counts_nd(..., keepdims=(...))` can retain an axis with only two edges as a length-1 dimension; otherwise such an axis is selected and collapsed.
- `counts_to_flux()` returns `(energy, flux, flux_error)`; flux and error use the same Energy-width, exposure, and geometric-factor denominator. `flux_e()`, `flux_t()`, and `flux_t_e()` are the only public flux wrappers. `flux_to_temperature()` returns temperature, temperature error, slope, slope error, and weighted R² arrays with the Energy axis removed; optional signed drift speed and look-direction cosine arrays apply the known-drift Maxwellian plasma-frame energy correction.

`ParticleDataFrame.select`, `between`, and `split` return particle-frame subsets. Its histogram methods delegate to `gena.post.statistic`. Per-event position and attitude are ordinary columns and have no interpolation methods. `ParameterDataFrame.pos_at(times)` and `rot_at(times)` are the PDF engineering interpolation APIs.

### Other missions

Product conversion belongs to `load`, never to `common.frame`. THEMIS and MMS obtain `(data, layouts)` from `load.api.spedas.read_spedas()`, whose mission-specific readers generate layouts directly using the product information already supplied for reading. Product constants (`THEMIS_FIELDS`, `MMS_FIELDS`, and `OMNI_FIELDS`), default-variable and name-qualification helpers, and layout generation all live in `load.api.spedas`; there is no separate layouts module. OMNI HRO declares scalar output and its default fields there. The public wrapper validates and dispatches, then applies layout overrides; no separate mission-dispatching layout generator is used. The consolidated `load.loaders.spedas` pipeline passes these layouts and decoded backend dictionaries through `parse_spedas()`. MMS backend names remain complete in the read API; the shared loader also preserves complete backend names. OMNI HRO decodes its shared `Epoch` and scalar fields independently; indices use `load.api.indices.read_indices()` to obtain records and layouts from Kyoto HAPI or GFZ JSON, then `parse_indices()` creates scalar descriptions. `INDICES_FIELDS` defines AE/AL/AU/AO, Dst, SYM/ASY, Kp, ap, daily Ap, sunspot number, and observed/adjusted F10.7. Source JSON is cached by UTC month under `prefix/indices`, with one request for shared Kyoto datasets. The index loader only assembles and concatenates intervals. `assemble()` mechanically joins normalized descriptions and prepares `Datetime` and `Timestamp` before constructing the frame. GENA Parquet is converted directly inside `load("gena", )` and does not pass through the pySPEDAS parser or time-series assembler.

Layouts explicitly specify kind, output columns, time dimension, and coordinate columns. Known product vectors may expand into components; unknown multidimensional fields remain array cells. A three-column array is never sufficient evidence for a spatial vector. Only declared energy coordinates receive an `_energy` name. THEMIS/MMS accept `variable_layouts` overrides for additional fields. pySPEDAS static one-dimensional axes and sampled multidimensional axes are handled at the backend boundary, including clipping when the channel count equals the time count.

The removed `dict_to_dataframe`, `dataframe_to_dict`, and `to_variable_dict` are not common APIs. Standard pandas `to_dict(orient=...)` is inherited unchanged. Legacy loader `return_type="dict"` merges native records directly; it does not reconstruct them through a DataFrame. `common.frame.read_variable` reads one declared variable for plotting without converting the entire table or importing loaders. The pandas-backed `common.time.Time` stores scalar Timestamp or 1-D DatetimeIndex values as timezone-naive UTC nanoseconds. Bare numbers are Unix seconds; milliseconds and scientific formats require explicit `format`. Use `.datetime` for tables, `to_timestamp()` for numeric seconds, and `to_tick()` for SpacePy; Astropy supplies scientific type recognition and conversion. Time is a converter, not an array container; indexing, slicing and array metadata belong to its native pandas `.datetime` attribute. The old independent conversion helpers have been removed. Daily expansion uses pandas date_range in callers; day starts use Time.floor("D"). Range consumers pass a two-element list to `Time(trange).datetime`; the loaders do not validate its length or ordering. DataFrame interval grouping lives in `common.frame.split_time_intervals` and directly processes native pandas datetime values without Time conversion; interpolation and miscellaneous helpers live in `libs.utils`.

THEMIS, MMS, OMNI, and index loaders return `TimeSeriesDataFrame` by default. Variables are outer-joined on the union of their native timestamps, retaining each variable's original sample values and cadence. Loaders never resample or interpolate; explicit legacy dictionary output also retains native timestamps. Missing samples at other variables' timestamps remain missing. Three-component vectors use `_x/_y/_z` scalar columns; spectra remain per-row arrays with a companion `_energy` array column. Provenance and the original variable layout are stored in `DataFrame.attrs`.

Pass `return_type="dict"` for the legacy dictionaries:

```python
data[variable]["x"]  # time axis
data[variable]["y"]  # scalar/vector/spectrum values
data[variable].get("v")  # optional energy or secondary axis
```

The shared frame implements `select`, `select_by`, `between`, `split`, and `resample_to`; these operations preserve the concrete DataFrame subclass and metadata. `target.select_by(selector, condition)` synchronizes the selector onto the target's native timestamps, evaluates the condition there, and returns untouched rows from the original target; nearest-neighbor selector synchronization is the default. `frame.resample(tres="4s")` builds a regular UTC grid inside the frame's coverage, while `frame.resample_to(other)` uses the reference frame's timestamps restricted to the intersection of both coverages. Numeric scalar columns are linearly interpolated without extrapolation, while discrete and array-valued columns use nearest samples. The default `max_gap="auto"` limits each column to three times its median native cadence so real outages remain missing; pass an explicit duration to override it or `None` to permit unlimited gap bridging. It also provides `frame.tplot(variable, ...)`, which delegates to `plots.basic.tplot`. The public entry is `from plots.basic import tplot`, implemented in `plots.basic.timeseries`; mission modules use it internally without re-exporting it. Variable-selection helpers live alongside this plotter, not in `plots.vis`. Scalar variables draw one line and metadata-backed three-component variables draw three lines; mission-specific spectrum helpers remain responsible for energy spectrograms.

### Coordinates and attitude

- Particle products store the position and attitude recorded for every event as ordinary columns; projection code reads those columns directly.
- `ParameterDataFrame.pos_at(times)` cleans PDF positions and interpolates them within PDF coverage.
- `ParameterDataFrame.rot_at(times)` cleans PDF quaternion samples and performs SLERP within PDF coverage.
- `gena.post.project.get_spacecraft_pos(l2df, coord)` converts positions already present in an L2 particle frame.
- `common.coord.Coordinates` stores copied NumPy `(N, 3)` components, UTC `datetime64[ns]` times, frame, unit, origin and optional NumPy xyzw attitude. `to(frame, backend="spacepy" | "astropy" | "pyspedas")` returns a new object; lazy adapters in `common._coord_backends` retain no backend objects. QQ/A/B relay through ECI2000 with explicit attitude; AGSM relays through GSM with a fixed 4.2-degree rotation. Conversion never changes origin. Astropy supports native ITRS/GCRS/CIRS/TEME and ICRS/Galactic positions with explicit length units; it does not alias GSE/GSM/ECI2000. Unsupported paths raise without fallback. `resample_at` linearly interpolates native components and SLERPs attitude without extrapolation. `Coordinates.to` delegates to the public `transform_vectors(values, source=..., target=..., times=..., backend=...)` array entry. Fixed local rotations need no time; standard backend conversions require UTC times. Individual frame-to-frame helpers are private. Callers use canonical uppercase names such as ECI2000 and GSM; the public layer performs no frame spelling or alias normalization. FOV, angle, quaternion-array validation, vector rotation and numeric-axis SLERP are also defined in `common.coord`; these functions remain stateless while `Coordinates` assigns the QQ→ECI2000 attitude meaning. The unused legacy `coord.rot` property is removed; the legacy constructor `rot` argument remains in use; old geometry spellings were removed and GENA callers migrated. Legacy `coord(times, values, frame, rot)` remains; `sync` clears derived arrays. Only `Coordinates` is re-exported at the common package level.

## 3. Function routing map

Use `.agents/skills/gena/scripts/api_query.py` before relying on this static map.

| Need | Start in | Principal public API |
|---|---|---|
| Resolve data directories | `common.config` | `resolve_paths`, `configure`, `get_configured_paths`, `resolve_output_path` |
| Normalize times, ranges, or UTC observing days | `common.time` | `Time` |
| Select, split, synchronize, concatenate, or read timestamped frames | `common.frame` | `split_time_intervals`, `TimeSeriesDataFrame.select_by`, `TimeSeriesDataFrame.split`, `TimeSeriesDataFrame.resample`, `TimeSeriesDataFrame.resample_to`, `TimeSeriesDataFrame.resample_at`, `concat`, `read_variable` |
| Interpolate position/rotation/data | `libs.utils` | `i_pos`, `i_rot`, `f_interp_vars` |
| Transform coordinates | `common.coord` | `Coordinates.to`, `Coordinates.resample_at`, `transform_vectors`; legacy `coord` remains available |
| Vector geometry and rotations | `common.coord` | `angle_between`, `angle_to_plane`, `vector_to_fov`, `fov_to_vector`, `quaternion_array`, `rotate_vectors`, `interpolate_quaternions` |
| Load GENA Parquet | `load` | `load("gena", ...)`; types: `gena.frame.ParameterDataFrame`, `gena.frame.ParticleDataFrame` |
| Interpolate PDF position or attitude | `gena.frame.ParameterDataFrame` | `pos_at`, `rot_at` |
| Save Matplotlib figures consistently | `plots.vis` | `savefig` |
| Build event masks | `gena.filter` | `check`, atomic rules, `soft_std`, `golden_std`, `diagnose`, `classify` |
| Apply relative energy efficiency | `gena.constants` | `eff_E`, `eff_H`, `eff_O` |
| Compute GENA science products | `gena.post.statistic` | `CountsNDResult`, `counts_nd`, `counts_to_flux`, `flux_e`, `flux_t`, `flux_t_e`, `flux_to_temperature` |
| Project rays/geometry | `gena.post.project` | `get_spacecraft_pos`, `get_ray`, `get_projection`, `projected_fov_axis_angle`, `axial_mean_angle`, `ray_sphere_intersections` |
| Export L2 geometry, observation summaries, and daily or exact-interval reports | `gena.post.export` | `dump`, `normalize_report_times`, `build_observation_summary`, `observation_summary`, `batch_plot`, `report` |
| Plot solar-wind/index context and QQ pointing geometry | `plots.advance`, `gena.utils` | `env`, `pointing_geometry`, `orbit` |
| Plot GENA 2-D diagnostics | `plots.basic.gena` | `energy`, `layer`, `fov`, `counts`, signal and engineering diagnostics |
| Build and save GENA animations | `plots.ani` | `AnimationDemo3`, `dumpFov` |
| Construct interactive 3-D scenes | `plots.v3d` | `SceneContext`, `prepare_scene`, `make_scene`, `dumpGeo` |
| Decode telemetry | `gena.raw.decode` | `decode` |
| Convert data levels | `gena.raw.conversion` | `l0_2_l1`, `l1_2_l1a`, `l1a_2_l1b`, `l1_2_l2` |
| Load THEMIS | `load` | `load("themis", ...)`, `load_pos("themis", ...)`; use `load("themis", ..., "esd", datatype="peif")` for ESA L2 3-D distributions |
| Compose preloaded mission panels | `plots.advance` | `Panel`, `overview`, `env` |
| Plot THEMIS | `plots.basic.themis` | `themis`, `spectrum`, `flux`, `get_ma_beta`, `velocity_parallel_perp` |
| Classify regions and plasma-sheet flows | `libs.regions` | `classify_liuzzo`, `classify_kiehas`, `classify_runov`, `classify_plasma_sheet_flows`; optional column mappings accept other missions |
| Plot classification strips | `plots.advance` | `plot_region_code` |
| Load MMS | `load` | `load("mms", ...)` |
| Load OMNI | `load` | `load("omni", ...)` |
| Load geomagnetic indices | `load` / `load.api.indices` | `load("indices", ...)`, `read_indices`, `parse_indices`; field definitions in `INDICES_FIELDS` |
| Combine solar-wind context | `load` | `load_sw` |
| Sun/Moon/satellite ephemeris | `libs.ephemeris` | `sun`, `moon`, `satellites`, `update_satellite_orbits` (`satellites` returns GSE positions in `R_E`, rejects out-of-coverage requests by default, and accepts opt-in `update=True`) |
| Plasma physics | `libs.plasma` | `alfven_velocity`, `beta`, `sigma_HH`, `f_maxwellian` |
| Magnetosphere regions/models | `libs.space` | region classifier, bow shock/magnetopause objects, `T96_fieldlines.trace_lines_y0`, `T96NeutralSheet.surface_xyz` |

For the complete current signature list, run:

```bash
python .agents/skills/gena/scripts/api_query.py --list-modules
python .agents/skills/gena/scripts/api_query.py spectrum
python .agents/skills/gena/scripts/api_query.py pos_at
python .agents/skills/gena/scripts/api_query.py --module gena.post map_T
```

Then open the returned `file:line` and inspect the implementation and docstring.
Use `rg -n "symbol_name\\("` to find call sites before changing a public API.

## 4. Loading and configuration

Set `FILE_BASE` or call `configure(file_base=..., output_dir=...)`, which persists explicit runtime choices. GENA range loaders, reports, and 3-D scene loaders always use the configured `pqt_base`; they do not accept a per-call data-root override. Relative artifact paths resolve under the configured result root. `GENA_OUTPUT_DIR` selects a separate root; otherwise `GENA_RES` under the data root is used. Reports default to its reports/ subdirectory. Absolute artifact paths override the output root. `common.config.resolve_paths()` is
read-only; `configure()` validates and creates standard directories.

```python
from common.config import configure
paths = configure(file_base="/data/ENA")
```

GENA supports daily discovery and direct files:

```python
from load import load

l2 = load("gena", trange=["2025-03-12", "2025-03-13"], product="l2_a")
one_file = load("gena", file="/data/20250312_l2_a.pqt")
```

Product suffixes are exact: `pdf`, `l1a_a/b`, `l1b_a/b`, `l2_a/b`, soft
`l2s_a/b`, and strict+spatial golden `l2g_a/b`.
`load(mission, trange, product, *, probe=..., key=...)` dispatches by mission: `product` is a GENA Parquet suffix or a SPEDAS instrument product; `probe` is accepted only by THEMIS/MMS (GENA encodes it in the suffix, OMNI and indices have none); `key` selects a subset and maps to `variables` for THEMIS/MMS/OMNI, `columns` for GENA, and the index product for indices. Every loader accepts one inclusive interval; request several intervals by calling it once per interval.

Mission loaders:

```python
from load import load
from load import load

esa = load("themis", trange, "esa", probe="b")  # strict local read by default
esa_density = load("themis", trange, "esa", probe="b", variables=["peir_density"])
downloaded_esa = load("themis", trange, "esa", probe="b", download=True)
auto_esa = load("themis", trange, "esa", probe="b", download="local_first")
esd = load("themis", trange, "esd", probe="b", datatype="peif", download=True)
fpi = load("mms", trange, "fpi", probe=1)
omni = load("omni", trange)  # BZ/BY GSM, speed, density, pressure
offline_omni = load("omni", trange, download=False)
auroral = load("indices", trange)  # one frame with ae, al, and au columns
dst = load("indices", trange, "dst")
context = load("indices", trange, ["dst", "sym_h", "kp", "ap_daily", "f107"], download="local_first")
legacy_esa = load("themis", trange, "esa", probe="b", return_type="dict")
figure, axis = esa.tplot("peir_avgtemp")
```

THEMIS defaults to strict local operation (`download=False`) and therefore never
contacts the Berkeley remote index unless requested. Use `download=True` to
allow remote access, `download="local_first"` to check local coverage and then
download if incomplete, or `download="refresh"` to force a new download. MMS
continues to default to `download="local_first"`. Its API reader sets the MMS environment and backend data directories on each read; there is no separate configuration entrypoint, import-time auto-configuration, or global wget download patch. Both loaders pass their default
variable lists to the backend, avoiding decoding unused CDF columns. MMS local
reads use the non-SPDF local-file backend and `notplot=True`, so an offline call
neither queries the SPDF index nor fills the pytplot registry.

THEMIS defaults cover reduced ion/electron ESA (`peir/peer`), full ion/electron
SST (`psif/psef`), ion/electron moments (`peim/peem`), spin-resolution FGM, and
the L1 `state` product. `state` is the position/velocity source behind
`load_pos("themis", ...)`, and its default level is `l1`
rather than `l2`; pass an explicit `level=` only for products that offer several.
Pass `variables=[...]` to replace the relevant default list; probe prefixes are
optional in these user-supplied names. Use the shared `variables=` argument for explicit field selection.
THEMIS ESD is a separate L2-only ESA distribution product. Call `load("themis", ..., "esd", datatype="peif")` for full-mode ion distributions; it defaults to the CDF-native `valid`, `data_quality`, `bins`, `eff`, and `eflux` variables. ESD arrays retain their three-dimensional `(time, energy-channel, solid-angle-bin)` values in `TimeSeriesDataFrame` cells, while pySPEDAS channel axes `v1`/`v2` are preserved as `<variable>_v1`/`<variable>_v2` companion array columns. Pass `return_type="dict"` to retain the direct pySPEDAS-style `x`/`y`/`v1`/`v2` dictionaries.
MMS retains its standard FPI and HPCA moment defaults. OMNI defaults to
`BZ_GSM`, `BY_GSM`, `flow_speed`, `proton_density`, and `Pressure`.
For NASA SPDF OMNI downloads, a loopback system proxy is bypassed only for the
NASA archive hosts; if another configured proxy returns no data, the loader
retries those hosts once with a direct connection.
Calling the index loader without ``product`` loads AE, AL, and AU into one
outer-joined frame; explicit strings still select individual index products.
Multi-interval requests merge dictionary output into one flat `x`/`y` record
for a single product, matching the single-interval shape. MMS requests forward
their `level` to the backend and derive FPI `datatype` from the requested
variables (`energyspectr_*` selects `dis-dist`, otherwise `dis-moms`).

THEMIS, MMS, OMNI, and index calls can download data. Do not trigger remote I/O
during a code-only review. Use existing local data or mocks for tests unless the
user requested data retrieval.

`plots.basic.themis.time_averaged_spectra()` returns one ``(energy, flux)`` pair.
Non-finite averaged channels are removed pairwise, and the combined ESA+SST
points are sorted by energy without interpolation across the instrument gap.
GENA ``flux_e`` and ``flux_t_e`` accept ``normalized=True, E_norm=(low, high)`` to divide each spectrum by its mean flux over that reference band; GENA ``E_norm`` is in keV, whereas the corresponding THEMIS visualization functions use eV. ``flux_t`` has already collapsed Energy and therefore cannot apply an Energy-band normalization.

## 5. GENA processing workflow

### L2 analysis

```python
from gena import filter as gf
from load import load
from gena.post import statistic as stat

l2 = load("gena", trange=trange, product="l2_a")
l2 = l2.select(gf.energy_eff(Emin=4, Emax=100))
weights = build_event_weights(l2)
centers, flux, error = stat.flux_e(l2, energy_bins, exposure=exposure_seconds, geometric_factor=geometric_factor, weights=weights)
```

Filter rules have two modes:

```python
mask = gf.energy_eff(l2, Emin=4, Emax=100)       # immediate boolean mask
rule = gf.energy_eff(Emin=4, Emax=100) & gf.onboard_eff()  # deferred
subset = l2.select(rule)
```

All filter rules return masks. Use `frame.select(rule)` for rows; `gena.filter.select` and `SelectorNode` were removed.

### Projection

```python
from gena.post.project import get_projection
from gena.post.export import dump

xyz = get_projection(l2, prob="A", coord="GSM")
table = dump(l2, prob="A", coord="GSM")
```

Projection results preserve event order and row count. Events with missing or invalid attitude have NaN directions and intersections.

### Available models and historical modules

`libs.flux` contains `j_maxwell`, `j_kappa`, `j_dual_model`, `fit_dual_moment_constrained`, and `qk_estimate_ena_counts`. Its `integral_ena_counts` is an unfinished draft that raises `NotImplementedError`. `gena.post.ena_drift`, `libs.distribution`, `libs.energy_transport`, `libs.bbf_ena`, and `libs.ena_response` are absent from the current source tree. Do not offer them as callable APIs. `libs.spectrum.normalize_flux_by_energy_band` is shared by GENA statistics and THEMIS spectra.

### Observation reports

```python
from gena.post import export

exact = export.report("A", ["2026-04-07T02:15:00Z", "2026-04-07T03:40:00Z"])
continuous = export.report("A", ["2026-04-07", "2026-04-08"])
discrete = export.report("B", ["2026-04-07", "2026-04-12", "2026-04-27"], output_dir="reports")
complete = export.report("A", "2026-04-07", include_animation=True, include_geometry=True)
catalog = export.observation_summary("A", ["2026-04-07", "2026-04-12", "2026-04-27"])
```

`report()` has one required `time` input: a scalar selects one UTC day, two date-only endpoints expand to inclusive daily reports, two endpoints containing clock precision select one exact interval that is not split at midnight, and longer date-only sequences select discrete days. `observation_summary()` uses exactly the same time normalization and writes one workbook spanning all selected intervals. Target classification uses the mean probe source-looking direction relative to GSM +X: Earth below 10 degrees, mid-tail from 15 to 40 degrees, far-tail above 120 degrees; other angles are unknown. There is no targets mapping parameter in the current API. Each report interval is clipped during product loading, passed through one `batch_plot()` call, and saved as independent 300 dpi JPEG figures; exact intervals use a `START_END` output key and daily reports keep `YYYYMMDD`. FOV MP4 and 3-D HTML export are opt-in through `include_animation=True` and `include_geometry=True`. `plots.advance.env()` plots caller-supplied OMNI/index frames and never loads or downloads data; `report(..., environment_sources=[...])` may include that preloaded context. `plots.basic.gena` owns detector/event diagnostics, `plots.advance.ani` owns Matplotlib animations, and `plots.advance.v3d` owns all interactive Plotly 3-D implementation. For one scene, `plots.advance.v3d.make_scene(time)` automatically loads that UTC day. For repeated times or an animation, build `context = plots.advance.v3d.prepare_scene(trange)` once, then call `context.make_scene(time)` or `context.make_timerange(...)`; pose, satellite-orbit, and space-weather files are loaded only during context construction. Reusable cross-mission plotting lives in `plots.advance`; study-specific science calculations remain in analysis scripts.

The Plotly scene renders all spatial coordinates in Earth radii. PDF spacecraft positions are Moon-relative MJ2000 metres, JPL Moon positions are kilometres, and `libs.ephemeris.satellites()` positions are GSE `R_E`; `plots.advance.v3d.scene()` performs the first two conversions and does not rescale satellite orbits. GSM `Z=0` and `Y=0` are rendered as translucent coordinate-reference planes by default and can be disabled with `include_reference_planes=False`; neither plane should be interpreted as the physical plasma sheet. By default, `SceneContext` calls `libs.space.T96_fieldlines.trace_lines_y0()` to trace 21 T96 field lines from northern nightside SM footpoints sampled between 70 and 82 degrees, after transforming each footpoint to GSM with the time-dependent dipole tilt from `geopack.recalc`; pressure, Dst, and IMF By/Bz specify the external model, and nonzero IMF By can carry the resulting three-dimensional lines away from GSM `Y=0`. It also constructs `libs.space.T96NeutralSheet`, a numerical total T96+IGRF `Bx_GSM=0` surface on a tail `(X,Y)` grid, with upper and lower translucent bounds separated by the user-configurable `neutral_sheet_half_thickness` visual half-thickness; this band is not a density-derived plasma-sheet boundary. The scene annotation updates with the field geometry and reports UTC, dipole tilt, Pdyn, Dst, IMF By/Bz, the A-FOV crossing of `Z=0`, and the modeled `Bx=0` position there. Single scenes use the exact requested time; animations cache field lines and neutral-sheet surfaces at `field_line_cadence="1h"` by default. Pass `include_field_lines=False` or `include_neutral_sheet=False` to omit either model product. Satellite current-position markers are plotted over separately cached trajectories spanning one hour before and after the instantaneous scene time, or one hour outside each end of an animated time range; configure this with `satellite_trajectory_padding` and `satellite_marker_size`. Satellite traces are omitted when a requested date lies outside the stored orbit coverage rather than extrapolating nonphysical coordinates. `satellites(..., update=True)` explicitly enables an SSCWeb network update through the latest requested timestamp. `update_satellite_orbits()` updates multiple keys atomically, preserves the local cadence, prefers definitive `themisa/themisd/themise/artemisp1/artemisp2` data, and uses their SSC predictive counterparts only beyond definitive coverage; MMS1 uses `mms1` directly. The companion `orbits_all.meta.json` records source observatories, coverage, coordinate system, and units.

### Raw decode

Use `gena.raw.decode.decode(raw_file, output_folder, by_day=True, snid=1)` only
for explicit conversion tasks; it writes parquet products and log files. Inspect
`gena/raw/constants.py` and environment-controlled pipeline flags before running.
Do not edit `gena/raw/jinfa_ref.py` as the active production path.

## 6. Cross-mission workflow

All mission frames already share `Datetime`/`Timestamp`; load each mission independently, inspect variable names/shapes, then explicitly synchronize with `TimeSeriesDataFrame.resample_to` only when a common cadence is scientifically required. All products, including THEMIS ESA, retain native outer-joined timestamps; no loader implicitly aligns or interpolates variables.
Keep provenance (probe, coordinate system, product level, units) in names or
metadata. Plotting helpers expect already-loaded dictionaries/DataFrames; they
should not own download or analysis policy.

For environment context at selected times, `load.load_sw()` combines
OMNI and Dst and falls back to defaults on load failure. Use lower-level loaders
when failure must be explicit rather than silently replaced.

## 7. Important distinctions and hazards

- Inspect source signatures: the repository is under active development and the
  README or this reference may lag.
- `gena.post.project.get_spacecraft_pos(l2df, coord)` reads the discrete position columns already stored on particle rows; PDF engineering evaluation uses `ParameterDataFrame.pos_at(times)` and `rot_at(times)`.
- Pass the correct probe (`a`/`b`, `A`/`B` as documented) through filtering,
  conversion, projection, and calibration.
- Treat flux defaults of area, solid angle, and exposure equal to `1` as proxies,
  not automatically calibrated physical flux.
- `Spec == 2` means oxygen in current L2 products; other numeric species values
  are treated as hydrogen by efficiency helpers.
- Prefer `l2_a/b` for standard L2, `l2s_a/b` for soft quality, and `l2g_a/b` only
  when strict plus spatial selection is intended.
- Preserve user data and output files. Loaders and raw decoding may perform I/O;
  plotting/statistics should normally work on in-memory objects.

TimeIndexedDataFrame owns UTC normalization, selection, interval splitting, plotting, and subclass preservation. ParticleDataFrame directly inherits it, retains pandas bin aggregation for event counts and histograms, and does not expose resample_to or interpolate event rows. TimeSeriesDataFrame implements interpolating resample(tres), resample_to(other), and resample_at(times); ParameterDataFrame inherits it and adds PDF-specific pos_at(times) and rot_at(times). select_by may sample a physical selector at event timestamps while keeping all original matching events. TimeSeriesDataFrame.resample returns a frame rather than a pandas Resampler; use pd.DataFrame(frame).resample(rule, on="Datetime") explicitly for aggregation.

`frame.resample_at(times)` accepts an explicit one-dimensional UTC time sequence (including Unix seconds), preserves its order and duplicates, and returns missing physical values outside source coverage. It is the shared public implementation behind resample and resample_to; ParticleDataFrame does not expose it.


## Index source boundary

`load/api/indices.py` owns Kyoto HAPI AE/Dst/SYM/ASY and GFZ Kp/ap/Ap/SN/F10.7 access, field definitions, cache and download policy, scalar decoding, and source metadata. `read_indices()` returns `(data, layouts)` and `parse_indices()` creates structured variables. `load/loaders/indices.py` owns one `load("indices", )` pipeline; it does not contain per-product loading wrappers. Default products are AE/AL/AU; `prefix` specifies the cache root and `download=False` is strictly offline, `"local_first"` reuses cache and fetches missing months, and `"refresh"` replaces validated monthly JSON. Source times are retained: Kyoto Dst at HH:29:30, GFZ daily labels at midnight rather than an invented noon; ap and daily Ap retain index units of 2 nT. Single-product dictionary output remains a flat x/y record. The old dst2 alias and text-cache/save_folder interface are replaced by best-available dst and the shared prefix/download interface. See `load/api/indices.md` for products and examples.

The current `load/load.py` facade imports the source loaders from `load/loaders` and read/parse APIs from `load/api`; it now also re-exports `assemble` through `load/__init__`. The position helpers behind `load_pos()` are private to the facade and accept a scalar or one-dimensional time sequence. The facade and package initializer no longer import or export the removed `esa_mode2_angles`, `themis_decompress`, or `read_themis_cdf` helpers.
