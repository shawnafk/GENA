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

THEMIS, MMS, and OMNI use the shared loading pipeline in `load.loaders.spedas.load_spedas()`. Shared interval normalization lives in `load.loaders.common`; the mission loaders own clipping, assembly, concatenation, and return-type handling. Only `load_spedas()` is exposed as the SPEDAS loader; mission-specific alignment and scientific helpers remain in the same module. GENA and indices retain separate source loaders.

The sole loading entrypoint is `load_spedas(mission, trange, product, ...)`. Mission-specific public wrappers and intermediate product-loading functions have been removed. All missions use `read_spedas()` followed by clipping, `parse_spedas()`, and `assemble()`. Dictionary output retains backend variable names and native axes. Default download policies are configured in `SPEDAS_OPTIONS` and resolved by the read API.

The indices path mirrors SPEDAS: `read_indices(trange, product, ...)` returns records and layouts, `parse_indices(data, layouts)` produces structured variables, and `load.loaders.indices.load_indices()` only handles interval iteration, assembly, concatenation, and output selection. `INDICES_FIELDS` declares Kyoto AE/Dst/SYM/ASY and GFZ Kp/ap/Ap/SN/F10.7 fields. All products share offline/local-first/refresh policies and monthly JSON caching.
