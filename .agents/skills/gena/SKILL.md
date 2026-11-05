---
name: gena
description: Navigate, query, call, diagnose, and extend the GENA-CODE-Dev scientific Python architecture for Queqiao-2 GENA ENA data, raw-to-L2 conversion, filtering, calibration, coordinate projection, statistics, plots, observation reports, and coordinated THEMIS, MMS, OMNI, geomagnetic-index, ephemeris, plasma, and magnetosphere analysis. Use whenever work targets this repository or mentions GENA functions, ENA_DataFrame, product levels/suffixes, detector filters, spacecraft pose/FOV, cross-mission loaders, space-environment context, or selecting the correct existing API before implementing new code.
---

# GENA

Use the repository as the source of truth. Route the task to the owning module,
query current signatures without importing side-effectful modules, inspect the
implementation, then call or modify the narrowest public API.

## Establish context

1. Work from the repository root containing `gena/` and `common/`. If the
   current directory is elsewhere, locate the root before resolving relative
   paths.
2. Read [references/architecture.md](references/architecture.md) sections 1–3
   for dependency direction, data contracts, and function routing.
3. Read the task-specific section of that reference:
   - Section 4 for configuration and loaders.
   - Section 5 for GENA analysis, projection, or raw decoding.
   - Section 6 for cross-mission work.
   - Section 7 before changing APIs or interpreting calibrated products.
4. Treat `deprecated/`, `TODO/`, generated artifacts, and paper-specific scripts
   as consumers or historical context, not the canonical API, unless explicitly
   placed in scope.

## Discover functions before calling them

Use the bundled AST query tool. It reads source without importing project modules
or triggering configuration/download side effects.

```bash
python .agents/skills/gena/scripts/api_query.py --list-modules
python .agents/skills/gena/scripts/api_query.py load
python .agents/skills/gena/scripts/api_query.py --module gena.post spectrum
python .agents/skills/gena/scripts/api_query.py 'pos|rot' --regex
```

Add `--legacy` only when the task explicitly requires historical raw modules.

For every candidate:

1. Open the reported `file:line` and read the complete function/class plus its
   directly used helpers.
2. Search call sites with `rg -n "symbol_name\\("` to confirm established usage,
   expected shapes, and compatibility impact.
3. Prefer functions exported through `__all__` or documented public entrypoints.
   Use private helpers only for diagnosis or when the public wrapper lacks the
   required capability.
4. Check the live signature again immediately before writing a call. Do not
   reproduce signatures solely from the reference file.

## Route by layer

Follow this ownership order:

- `common`: configuration, data types, time normalization/interpolation, coordinates, and physical constants. `plots.vis` owns shared plotting primitives.
- `gena.frame`, `gena.filter`, `gena.constants`: GENA table types, detector masks, and calibration/efficiency. `load` owns local GENA product loading.
- `gena.post`: derived science statistics, ray/projection geometry, and reports.
- `gena.raw`: telemetry packet parsing and raw/L0/L1/L2 conversion with writes.
- `load`: GENA and geomagnetic-index loading through `load_gena`/`load_indices`, and THEMIS/MMS/GOES/OMNI through the sole `load_spedas(mission, ...)` entrypoint; external observations are loaded into `TimeSeriesDataFrame`. `libs.regions` owns cross-mission region classification.
- `plots.basic`: GENA, THEMIS, MMS, and OMNI product plots. Old mission `vis`
  modules have been removed; `plots.vis` owns primitives.
- `plots.advance`: reusable cross-mission composition of preloaded data
  (`Panel`, `overview`, `env`); no implicit downloading or synchronization.
- `libs`: ephemeris, plasma physics, space-environment composition, and boundary
  models.
- analysis directories such as `1.Paper_BBF`: study-specific orchestration only.

Put new behavior at the lowest reusable owner. Do not duplicate a shared helper
inside an analysis script or create sideways mission dependencies.

## Work with data safely

- Normalize time inputs through `common.time.Time`; use `.datetime` for native pandas UTC nanoseconds and explicit `format` for numeric inputs other than Unix seconds.
- Inspect columns, dictionary keys, array shapes, units, probe, level, and
  coordinate frame before calculation.
- Expect GENA loaders to return `ENA_DataFrame`; expect THEMIS/MMS/OMNI/index loaders to return `TimeSeriesDataFrame` by default, with per-variable `x`/`y` dictionaries available only through the explicit legacy return option.
- Keep filtering as boolean masks until row selection. `FilterNode` supports
  `&`, `|`, and `~`; Row selection is owned by `TimeSeriesDataFrame.select`; no separate selector class is retained.
- Read stored GENA pose with `frame.pos()` (Moon-relative MJ2000 metres) and `frame.rot()` (QQ→ECI2000 xyzw arrays, invalid rows NaN); both preserve row order without interpolation. Only PDF `ParameterDataFrame.pos_at(times)` and `rot_at(times)` interpolate engineering data.
- Define FOV targets separately from observers: sphere center/radius, curve points, or discrete points. Use `gena.post.project.obj_in_view` for viewing and `dipole_fieldlines` for reference target definitions; lunar blocking shares this projection. Prepare ephemeris, satellite positions and model curves before calling it.
- Pass probe identity consistently through filters, conversion, calibration,
  projection, and plotting.
- Treat unit-area, unit-solid-angle, or unit-exposure flux defaults as proxies
  unless mission calibration values are supplied.

## Control I/O boundaries

Inspect/review tasks do not authorize downloads or output writes.

- `resolve_paths()` only resolves configuration; `configure()` may validate and
  create directories.
- THEMIS, MMS, OMNI, and index loaders may access networks.
- `gena.raw.decode.decode()` writes parquet and logs.
- Report Excel export writes a workbook.

For tests, prefer explicit temporary paths, small in-memory frames, local files,
and mocked loader boundaries. Perform real downloads or raw conversion only when
the user requests them or they are an inherent implementation step.

## Implement and verify

When changing code:

1. Preserve public return contracts and product semantics unless the requested
   change explicitly alters them.
2. Update docstrings and [references/architecture.md](references/architecture.md)
   when module ownership, public names, data contracts, or loader/product
   behavior changes.
3. Add the smallest test that exercises the changed boundary. Include both
   probe A and B when probe-specific calibration or filtering is involved.
4. Run targeted tests, `python -m py_compile` on edited modules, and import smoke
   checks. Avoid network-backed tests by default.
5. State clearly when validation used synthetic data only and when real
   telemetry/calibration validation remains outstanding.

## Common call sequence

Use this as a routing example, not as a substitute for checking live signatures:

```python
from gena import filter as gf
from load import load_data
from gena.post import statistic as stat

l2 = load_data("gena", trange=trange, product="l2_a")
l2 = l2.select(gf.energy_eff(Emin=4, Emax=100) & gf.onboard_eff())
weights = build_event_weights(l2)
centers, flux, flux_err = stat.flux_e(l2, energy_bins, exposure=exposure_seconds, geometric_factor=geometric_factor, weights=weights)
```

For loader, projection, raw decode, and cross-mission examples, read the
corresponding section of [references/architecture.md](references/architecture.md).
