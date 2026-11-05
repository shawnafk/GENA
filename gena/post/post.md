# GENA L2 Post Processing

Current-source audit (2026-09-17): only `statistic.py`, `project.py`, and `export.py` are present here. The linked historical science notes and `ena_drift.py` sections below describe previous work, not available modules. See [the current design review](../../docs/design_review_20260917.md) for validation and missing-file details.

This folder contains post-processing functions for GENA Level-2 particle rows. At this stage each row is treated as one ENA event with time, direction, position/projection, energy, and species; optional statistical weights are supplied as external arrays aligned positionally with the rows.

The end-to-end handoff for separating ion leakage from remote ENA candidates, distinguishing lobe/plasma-sheet regions, constructing matched backgrounds, and grading events is documented in [`ion_leakage_region_event_selection_handoff.md`](ion_leakage_region_event_selection_handoff.md).

The detailed archive for the 2025-11-12/13 storm-time Earth-pointing observation, including energy-resolved FOV morphology, source-ray geometry, substorm timing, unbinned de-dispersion, THEMIS-B data-quality corrections, and the quantitative ion-leakage budget, is documented in [`earth_pointing_observation_20251112_analysis.md`](earth_pointing_observation_20251112_analysis.md).

The statistical assumptions, formulas, and limitations behind `hist()` and `ave()` are documented in [`statistics.md`](statistics.md).

The Scime temperature derivation, `counts_nd` processing chain, API examples, map conventions, and interpretation limits are documented in [`temperature.md`](temperature.md).

The boundary is:

- Methods that return a `ParticleDataFrame` subset belong on `ParticleDataFrame` itself, for example `select`, `between`, and `split`.
- Functions that return derived statistics, spectra, fluxes, maps, or fitted physical quantities belong in the post layer.

## Data Contract

Expected columns by feature:

| Feature | Required columns |
|---------|------------------|
| Time binning | `Timestamp` or `Date` |
| Weighted statistics | external weight array aligned positionally with the event rows |
| Energy spectra | `Energy`; optional weights are passed as an external array |
| Weight calibration | `Energy` and `Spec`; keep the resulting weights outside the event DataFrame |
| FOV maps | `PX`, `PY`, `Energy`; optional weights are passed separately |
| Projection maps | external projected `x/y` arrays plus `Energy`; optional weights are passed separately |
| Ray/projection geometry | `Timestamp`, `MJ2000_Star_X/Y/Z`, `MJ2000_q1/q2/q3`, `X1/Y1/X2/Y2` |
| GSM two-band drift inversion | Ray/projection geometry and `Energy`; explicit energy and GSM longitude/latitude edges |

## Implemented Functions

### `statistic.py`

- `CountsNDResult`, `counts_nd`: N-dimensional raw and optionally weighted event counts without flux normalization; weights are external arrays, two-edge axes are selected and collapsed by default, and `keepdims` can retain a one-bin axis.
- `counts_to_flux`: selects the matching raw or weighted counts and errors, divides both by Energy width, exposure, and geometric factor, and optionally multiplies both by the Energy center when `en=True`.
- `flux_e`: Energy differential flux from explicit Energy edges; exposure may be a scalar or an array broadcastable to the spectrum.
- `flux_t`: time-dependent band-averaged differential flux from explicit Date edges and one two-edge Energy range.
- `flux_t_e`: time–Energy differential flux from explicit Date and Energy edges.
- `df_resample`: time-window grouping for derived statistics.
- `window_hist0d`, `window_hist1d`, `window_hist2d`, `window_ave`: generic binned statistics.
- `flux_to_temperature`: performs inverse-variance weighted zero-drift or known-drift Maxwellian fits along any Energy axis and returns temperature, temperature error, slope, slope error, and weighted R² arrays with that axis removed.

Flux proxy convention:

```text
J(E, t) ~= sum_i W_i / (Delta E * Delta t * G)
```

The default one-sigma count error is:

```text
sigma_J ~= sqrt(sum_i W_i^2) / (Delta E * Delta t * G)
```

where `Delta E` must come from an explicit input energy window or explicit energy-bin edges. `W` denotes an external per-event weight array aligned positionally with the input rows; it is not stored in the event DataFrame. `G` is the user-provided geometric factor. If it is left as `1`, the result is a normalized count-rate proxy rather than an absolute physical flux.

This error term is the propagated counting/statistical uncertainty after event weighting. It does not include calibration uncertainty from the efficiency model.

## External Event Weight Convention

Construct weights outside the event DataFrame and pass the resulting one-dimensional array directly to statistical functions:

```python
weights = 1.0 / efficiency_func(l2["Energy"].to_numpy())
result = counts_nd(l2, axes={"Energy": energy_edges}, weights=weights)
```

For each event:

```text
efficiency = fitted efficiency(Energy, Spec) / fitted efficiency(10 keV, Spec)
weight     = 1 / efficiency
```

Thus the relative efficiency and weight are both 1 at 10 keV. The external weight array must have the same length and row order as the event DataFrame.

### `project.py`

本模块按“基础数组几何 → 事件几何 → 目标定义 → 通用视场投影 → 月球遮挡”排列。目标定义与观察者分离，不加载工程数据、太阳风或星历，不自动插值，不保留旧接口兼容包装。

| 分区 | 公开函数 | 输入与结果 |
| --- | --- | --- |
| 基础数组几何 | `plane_intersections`、`ray_sphere_intersections`、`sphere_limb`、`axial_mean_angle` | 数组求交、球体精确切线轮廓及轴方向平均 |
| 位置转换 | `position_as_coord` | 位置数组加可选偏移后转换坐标系，返回坐标对象 |
| 事件几何 | `ray_from_L2`、`event_axes` | 逐行读取姿态，返回射线或三根坐标轴的坐标对象 |
| 事件投影 | `project_events`、`projected_fov_axis_angle` | 返回交点数组，或每行角度及平均角度 |
| 目标定义 | `dipole_fieldlines` | 生成命名的 SM/R_E 参考偶极磁力线，不依赖观察者或姿态 |
| 通用视场 | `obj_in_view` | 输入观察者、姿态、目标字典，返回各目标的角度、视场掩码、裁剪线段和可见性标记 |
| 月球遮挡 | `moon_blocking_mask`、`moon_blocking_intervals` | 从 GENA 表提取同行位置和姿态，定义月球球体并复用 `obj_in_view` |

空间位置使用 R_E；GENA 表的 `pos()` 返回月心 MJ2000 米位置，进入位置转换或视场函数前除以 `RE`。`rot()` 返回逐行 QQ→ECI2000 的 xyzw 数组，无效姿态保留 NaN。`position_as_coord` 先添加与位置使用相同源坐标系和单位的 `pos_offset`，再从 `ref_coord` 转到 `res_coord`。L2 的 `project_events`、`projected_fov_axis_angle` 和导出 `dump` 需显式传入地心 ECI2000 月球位置偏移。`gena.utils.velocity` 复用 `ray_from_L2`，生成单位粒子飞行方向和 km/s 速度；飞行方向为层 1 → 层 2，仪器观测方向为探头 -Z。

`obj_in_view(time, position, attitude, objects, ...)` 每次处理一个 UTC 时刻和观察者位置。`position` 使用 `ref_coord`（默认 ECI2000）；各目标可指定自己的 `coord`，但所有目标与观察者必须共用原点和 R_E 单位。视场朝 `probe` 的 -Z，`ax_deg`、`ay_deg` 为矩形半视角，默认 22.5°、10°。目标名称没有特殊含义，地球、月球或其他球体都使用同一种定义。

| `kind` | 几何字段 | 用途 |
| --- | --- | --- |
| `sphere` | `center` 三维位置、`radius` | 地球、月球等球体；投影时根据观察者计算切线轮廓 |
| `curve` | `points` 为 `(N, 3)` 数组 | 磁力线、磁尾边界、轨迹、预先准备的 limb 曲线 |
| `points` | `points` 为 `(N, 3)` 数组或一个 `(3,)` 位置 | 卫星或其他离散目标 |

返回字典保留目标名称和元数据。`px_py` 为前半球的角度，后方或无效点为 NaN；`in_view` 为逐点视场掩码，`visible_px_py` 额外将视场外点置为 NaN。曲线和球体轮廓的 `segments_px_py` 为 `(M, 2, 2)` 裁剪线段，两端都在视场外但穿过视场的线段也会保留。`any_in_view` 表示目标与视场相交，`all_in_view` 表示整个目标在视场内；球体另有 `center_px_py` 和 `covers_view`（球体视盘覆盖整个视场）。球体的这三个布尔判定使用解析几何，`samples` 只控制轮廓绘图精度；其他目标按给定点或折线处理。各目标独立投影，不计算目标之间的遮挡。

```python
from load import load_data
from libs import ephemeris
from common.constants import RE, RM
from common.time import Time
from gena.post.project import project_events, dipole_fieldlines, obj_in_view
from gena.utils import velocity
import pandas as pd

# 逐事件计算：加载与月球星历准备放在调用端。
events = load_data("gena", trange, "l2_a")
vectors = velocity(events, probe="A", coord="GSM")
pos_offset = ephemeris.moon(events["Datetime"], geocentric=True) * 1000 / RE
xyz = project_events(events, pos_offset=pos_offset, probe="A", coord="GSM",
                     plane_normal=(1, 0, 0), plane_point=(-65, 0, 0))

# 定义目标：与观察者姿态和视场无关。
stamp = Time(time).datetime
moon_position = ephemeris.moon(stamp, geocentric=True) * 1000 / RE
objects = {
    "earth": {"kind": "sphere", "center": (0, 0, 0), "radius": 1., "coord": "ECI2000"},
    "moon": {"kind": "sphere", "center": moon_position, "radius": RM / RE, "coord": "ECI2000"},
    **dipole_fieldlines(l_shells=(4, 8)),
}

# 准备观察者，再把同一组目标交给通用视场函数。
pad = pd.Timedelta("1h")
pose = load_data("gena", [stamp - pad, stamp + pad], "pdf")
position = pose.pos_at([stamp])[0] / RE + moon_position
attitude = pose.rot_at([stamp])
image = obj_in_view(stamp, position, attitude, objects, probe="A")
```

磁尾模型在调用端通过 `load_sw` 和 `libs.space.T96_fieldlines` 生成，然后将每条 GSM/R_E 曲线定义为 `{"kind": "curve", "points": line, "coord": "GSM"}` 加入 `objects`。卫星当前位置使用 `{"kind": "points", "points": positions, "coord": "GSE"}`，位置应提前加载并对齐到观测时刻。`plots.gena.fov(events, projection=image)` 可直接绘制这些通用结果；目标的可选 `style` 字典设置颜色、线宽和图例文字。

月球遮挡使用月心位置，因此月球目标中心为原点，不需要读取月球星历。先联合筛选有效位置和姿态，再按 `gap` 抽样；返回的 `partial_block` 表示月球视盘与视场相交，`full_block` 表示整个月球视盘位于视场内，`covers_view` 才表示月球覆盖整个视场。`moon_blocking_intervals(full_disk=True)` 按 `full_block` 分组，默认按 `partial_block` 分组。

旧的 `instrument_projection` 已移除，目标定义与 `obj_in_view` 直接替代其功能，不保留别名。

### Historical `ena_drift.py` (not present)

二维特定能段漂移拟合的完整物理推导、调用示例和结果解释见 [`ena_drift.md`](ena_drift.md)。

- `ena_velocity_unit_vectors_gsm`: returns the physical layer-1-to-layer-2 ENA particle-flight direction in GSM, using the same convention as `ray_from_L2` and `gena.utils.velocity`.
- `gsm_angles_from_unit_vectors`: converts physical ENA velocity directions to GSM longitude and latitude.
- `bin_two_band_gsm_counts`: makes raw `(band, longitude, latitude)` counts from two non-overlapping energy bands. Longitude is unwrapped from the first supplied edge, so a grid may cross -180/180 degrees.
- `vb_from_ena_2d`: fits one speed along a supplied GSM direction from the conditional two-band angular-count likelihood, with `eff_H(E)` and `sigma_HH(E)` folded into the forward model.

`libs.ena_drift.rest_flux` provides the isotropic zero-bulk-velocity Maxwellian-core + Kappa-tail proton flux. `libs.ena_drift.rest_2_lab` transforms that spectrum using one scalar bulk speed and averages over all direction cosines for comparison with an omnidirectional laboratory-frame spectrum. `gena.post.ena_drift.vb_from_ena_2d` is the directional two-band fit when the GSM drift direction is supplied externally.

### `export.py`

- `dump`: export event geometry, spacecraft frame axes, and normalized ENA rays to a tabular product.
- `batch_plot`: call the reusable 2-D functions in `plots.basic.gena` once to build all available static report figures; data loading and file saving remain separate.
- `build_observation_summary`: split already-loaded PDF engineering samples into large and continuous small observation intervals, calculate the QQ orbit and probe-boresight projection, classify Earth/mid-tail/far-tail pointing from the probe `-Z` angle to the GSM `+X` Sun-Earth axis, and count matching ENA events.
- `observation_summary`: use the same automatic time semantics as `report`, load PDF plus the requested ENA product, and write one `observation_summary.xlsx` workbook for a day, inclusive date range, exact interval, or sequence of discrete days.
- `report`: automatically normalize one time input, load and clip daily products, build and save independent 300 dpi JPEG figures including QQ/FOV-plane geometry and the L1A joint distribution, and optionally invoke `plots.ani.dumpFov(...)` and `plots.v3d.prepare_scene(...)` for MP4 and 3-D HTML output. Preloaded OMNI/index frames may be passed through `environment_sources`; the report never downloads those products implicitly. The 3-D scene uses GSM coordinates in `R_E`; unavailable satellite ephemerides are omitted instead of extrapolated.

```python
from gena.post import export

export.report("A", ["2026-04-07", "2026-04-08"])
export.report("B", ["2026-04-07", "2026-04-12", "2026-04-27"], output_dir="reports")
export.report("A", "2026-04-07", include_animation=True, include_geometry=True)
summary = export.observation_summary("A", ["2026-04-07", "2026-04-08"], output_dir="reports")
```

Observation targets use the mean probe `-Z` look direction within each segment. Its angle to the GSM `+X` Sun-Earth axis is classified as `earth` below 10 degrees, `mid_tail` from 15 through 40 degrees, and `far_tail` above 120 degrees. Uncovered angle ranges are labeled `unknown`; the measured angle is retained as `sun_earth_angle_deg` in the workbook.

### `gena.utils`

- `plots.advance.env(*sources, cols=..., f=None, ax=None)`: plot selected columns from caller-supplied solar-wind and geomagnetic-index frames. Each `cols` entry creates one dynamic panel; a tuple groups several columns on one panel.
- `pointing_geometry(pdf, probe=..., coord=...)`: return the Earth-centred QQ orbit, nominal source-looking FOV direction, and intersections with the XY, XZ, and YZ coordinate planes.
- `orbit(...)`: plot the reusable `pointing_geometry` result with a dynamic panel layout.

## Additional ENA Post Products Worth Keeping

For ENA L2 analysis, the post layer should cover these product families:

- Exposure/quality masks: good-time intervals, dead-time handling, and invalid detector-state exclusion.
- Calibrated response correction: energy efficiency, spatial efficiency, geometric factor, and FOV solid angle.
- Spectral products: all-species and species-separated energy spectra, time-energy spectra, and background-subtracted spectra.
- Imaging products: detector FOV maps, projected-plane maps, and line/profile cuts.
- Uncertainty products: weighted-count errors, propagated background subtraction errors, and fit quality flags.
- Source diagnostics: temperature maps/profiles, spectral slopes, and species composition ratios.

The current implementation covers the basic statistical, spectral, imaging, background-subtraction, and temperature-slope products.  Exposure/dead-time handling and fully calibrated response constants still require mission calibration inputs.

## Temperature Inference Derivation

The current temperature estimate follows Scime et al. (2002) and assumes the observed ENA differential flux in one spatial/time bin can be approximated as:

```text
J_ENA(E) = C sigma(E) E exp(-E / T_i)
```

Definitions:

- `E`: ENA energy in keV.
- `J_ENA(E)`: observed differential ENA flux or flux proxy.
- `sigma(E)`: charge-exchange cross section.
- `T_i`: source ion temperature in keV.
- `C`: slowly varying line-of-sight/geometric/density factor treated as
  constant within one bin.

Taking the logarithm gives:

```text
ln[J_ENA(E) / (sigma(E) E)] = ln(C) - E / T_i
```

Therefore, fitting:

```text
y = ln[J_ENA(E) / (sigma(E) E)]
x = E
y = a + b x
```

gives:

```text
T_i = -1 / b
```

For a known signed bulk drift speed `V_d` along a chosen positive reference axis and a look-direction cosine `mu`, the plasma-frame energy is:

```text
E_d     = m V_d^2 / 2
E_prime = E + E_d - 2 sqrt(E E_d) mu sign(V_d)
```

The drift-aware fit replaces `x = E` by `x = E_prime` while retaining the measured Energy in `sigma(E) E`. Pass `drift_velocity_kms` and `look_direction_cosine` to `flux_to_temperature`; both can be scalars or arrays broadcastable to the temperature output shape. Their defaults give the original zero-drift result.

Implementation details:

- Energy is histogrammed by `counts_nd` and divided by each explicit bin width through the same normalization core used by `flux_e`, `flux_t`, and `flux_t_e`.
- Optional external event weights are positionally aligned with events; all efficiency and detector-response correction remains the caller's responsibility.
- The caller chooses the Energy edges used to construct the flux; the default fit uses centers from 8 through 25 keV.
- Bins with non-positive flux, non-positive energy, or non-positive cross section are excluded.
- The fit propagates `flux_err / flux` into log space, uses inverse-variance weights, and reports the matching weighted coefficient of determination.
- A non-negative fitted slope is treated as invalid and returns `NaN` temperature.
- Drift velocity is an external known or scanned parameter; its uncertainty is not included in the returned statistical temperature error.

Important limitation:

This is a local spectral-slope estimate. It is not a full line-of-sight inversion, and a finite one-dimensional spectrum generally cannot uniquely separate temperature, Kappa-tail hardness, and bulk drift. A single direction cosine is exact for one look direction; a broad FOV direction distribution requires response-weighted forward modeling.
