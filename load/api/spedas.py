"""Read pySPEDAS products and normalize their backend variable dictionaries.

This module owns the boundary between remote or cached CDF products and the
repository's loader layer.  ``read_spedas`` keeps the mission-specific
pySPEDAS call signatures in one place; ``parse_spedas`` only converts the
returned ``x``/``y``/``v`` mapping into the explicit variable descriptions
consumed by an assembler.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import os
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common.config import get_configured_paths
from common.time import Time


__all__ = ["parse_spedas", "read_spedas"]

"""Shared mission and product rules for backend variable conversion."""
THEMIS_FIELDS = {
    "esa": {
        "default": lambda p: [f"th{p}_{m}_{f}" for m in ("peir", "peer") for f in ("avgtemp", "en_eflux", "data_quality", "density", "velocity_gsm", "flux")],
        "modes": ("peeb", "peef", "peer", "peib", "peif", "peir"),
        "vectors": ("flux", "eflux", "magf", "velocity_dsl", "velocity_gse", "velocity_gsm"),
        "energy": ("en_eflux",),
    },
    "sst": {
        "default": lambda p: [f"th{p}_{m}_{f}" for m in ("psif", "psef") for f in ("avgtemp", "en_eflux", "data_quality", "density", "velocity_gsm", "flux")],
        "modes": ("psif", "psef", "psir", "pser"),
        "vectors": ("flux", "eflux", "magf", "velocity_dsl", "velocity_gse", "velocity_gsm"),
        "energy": ("en_eflux",),
    },
    "mom": {
        "default": lambda p: [f"th{p}_{m}_{f}" for m in ("peim", "peem") for f in ("data_quality", "density", "flux", "eflux", "velocity_gsm", "ptot")],
        "modes": ("peim", "peem"),
        "vectors": ("flux", "eflux", "magf", "velocity_dsl", "velocity_gse", "velocity_gsm"),
        "energy": (),
    },
    "fgm": {
        "default": lambda p: [f"th{p}_fgm_fgs_quality", f"th{p}_fgs_btotal", f"th{p}_fgs_gsm"],
        "modes": (), "vectors": (), "energy": (),
    },
    "esd": {"default": lambda p: ["valid", "data_quality", "bins", "eff", "eflux"], "modes": (), "vectors": (), "energy": ()},
}

MMS_FIELDS = {
    "fpi": {
        "default": lambda p: ["temppara_fast", "tempperp_fast", "numberdensity_fast", "bulkv_gse_fast"],
        "vectors": {f"{f}_{c}_{r}" for f in ("bulkv", "bulkv_spintone", "heatq") for c in ("dbcs", "gse", "err") for r in ("fast", "brst")},
        "energy": {f"energyspectr_{d}_{r}" for d in ("px", "mx", "py", "my", "pz", "mz", "omni") for r in ("fast", "brst")},
        "stem": "dis",
    },
    "hpca": {
        "default": lambda p: [f"{s}_{n}" for s in ("hplus", "heplus", "oplus", "heplusplus") for n in ("number_density", "ion_bulk_velocity_GSM", "tpara", "tperp")] + ["hplus_ion_bulk_velocity", "hplus_scalar_temperature"],
        "vectors": {*(f"{s}_ion_bulk_velocity{suffix}" for s in ("hplus", "heplus", "heplusplus", "oplus") for suffix in ("", "_GSM")), "B_GSE_spin_avg", "B_GSM_spin_avg"},
        "energy": set(), "stem": "hpca",
    },
}

OMNI_FIELDS = {
    "hro": {
        "default": ["BZ_GSM", "BY_GSM", "flow_speed", "proton_density", "Pressure"],
        "kind": "scalar",
        "time_axis": 0,
    },
}

SPEDAS_OPTIONS = {
    "themis": {
        "products": list(THEMIS_FIELDS),
        "probes": ["a", "b", "c", "d", "e"],
        "default_probe": "a",
        "download": False,
    },
    "mms": {
        "products": list(MMS_FIELDS),
        "probes": ["1", "2", "3", "4"],
        "default_probe": "1",
        "download": "local_first",
    },
    "omni": {
        "products": list(OMNI_FIELDS),
        "probes": [],
        "default_probe": None,
        "download": True,
    },
}


def default_variables(mission, product, probe):
    table = {"themis": THEMIS_FIELDS, "mms": MMS_FIELDS, "omni": OMNI_FIELDS}[mission]
    try:
        defaults = table[product]["default"]
        return defaults(probe) if callable(defaults) else list(defaults)
    except KeyError as exc:
        raise ValueError(f"unsupported {mission} product {product!r}") from exc


def qualify_variables(mission, product, variables, probe):
    if variables is None:
        return None
    if mission == "omni":
        return list(variables)
    if mission == "themis":
        if product == "esd":
            return [str(name).lstrip("_") for name in variables]
        prefix = f"th{probe}_"
        return [name if str(name).startswith("th") else prefix + str(name).lstrip("_") for name in variables]
    stem = f"mms{probe}_{MMS_FIELDS[product]['stem']}_"
    return [name if str(name).startswith("mms") else stem + str(name) for name in variables]


def sample_is_present(value):
    if isinstance(value, (np.ndarray, list, tuple)):
        return True
    missing = pd.isna(value)
    return not bool(missing) if np.ndim(missing) == 0 else True


def _variable_layouts(data, *, vectors=(), energy=()):
    """Describe arrays and auxiliary axes using explicit variable name sets."""
    layouts = {}
    for name, variable in data.items():
        kind = "vector" if name in vectors else "scalar" if np.ndim(variable["y"]) == 1 else "array"
        axes = {}
        for key, value in variable.items():
            if re.fullmatch(r"v\d*", key):
                is_energy = key == "v" and name in energy
                axes[key] = {
                    "column": f"{name}_energy" if is_energy else f"{name}_{key}",
                    "time_axis": None if np.ndim(value) <= 1 else 0,
                    "role": "energy" if is_energy else None,
                }
        layouts[name] = {"kind": kind, "time_axis": 0, "columns": [f"{name}_{c}" for c in ("x", "y", "z")] if kind == "vector" else [name], "axes": axes}
    return layouts


def _configured_spedas_dir(prefix):
    if prefix is not None:
        return str(Path(prefix).expanduser())
    return get_configured_paths()["spedas_base"]


def _read_themis(trange, *, product, variables, probe, level, prefix, download, **kwargs):
    from pyspedas.projects.themis.config import CONFIG
    from pyspedas.projects.themis.load import load as backend_load

    CONFIG["local_data_dir"] = _configured_spedas_dir(prefix) + "/themis_data/"
    call_kwargs = dict(kwargs)
    call_kwargs.update(
        instrument=product,
        probe=probe,
        level=level,
        notplot=True,
        time_clip=True,
    )
    if variables is not None:
        call_kwargs["varnames"] = variables
    call_kwargs.pop("no_update", None)
    if download in (False, "offline", "local"):
        data = backend_load(trange=trange, no_update=True, **call_kwargs) or {}
    elif download in (True, "download", "local_first", "auto", "refresh"):
        data = backend_load(
            trange=trange,
            no_update=False,
            force_download=download == "refresh",
            **call_kwargs,
        ) or {}
    else:
        raise ValueError(
            "download must be False/'local'/'offline', True/'download', "
            "'local_first'/'auto', or 'refresh'"
        )
    spec = THEMIS_FIELDS[product]
    root = f"th{probe}_"
    vectors = {root + mode + "_" + field for mode in spec["modes"] for field in spec["vectors"]}
    energy = {root + mode + "_" + field for mode in spec["modes"] for field in spec["energy"]}
    if product == "fgm":
        vectors.update(root + rate + "_" + coord for rate in ("fgs", "fgl", "fgh") for coord in ("dsl", "gse", "gsm"))
    return data, _variable_layouts(data, vectors=vectors, energy=energy)


def _read_mms(trange, *, product, variables, probe, prefix, download, **kwargs):
    from pyspedas.projects import mms as backend

    directory = _configured_spedas_dir(prefix) + "/mms_data/"
    os.environ["MMS_DATA_DIR"] = directory
    from pyspedas.projects.mms.mms_config import CONFIG

    CONFIG["local_data_dir"] = directory
    if product not in ("fpi", "hpca"):
        raise ValueError("MMS product must be 'fpi' or 'hpca'")
    names = list(variables or [])
    full_names = []
    stem = f"mms{probe}_dis_" if product == "fpi" else f"mms{probe}_hpca_"
    for name in names:
        full_names.append(name if str(name).startswith("mms") else stem + str(name))

    def call(no_update):
        common = dict(
            trange=trange,
            probe=str(probe),
            spdf=False,
            no_update=no_update,
            notplot=True,
            time_clip=True,
            varnames=full_names,
            **kwargs,
        )
        if product == "fpi":
            return backend.mms_load_fpi(datatype="dis-moms", **common) or {}
        return backend.mms_load_hpca(**common) or {}

    if download in (False, "offline"):
        raw = call(True)
    elif download in (True, "download", "local_first", "auto", "refresh"):
        raw = call(False)
    else:
        raise ValueError("download must be local/offline, download/refresh, or local_first/auto")
    data = {name: value for name, value in raw.items() if name in full_names}
    spec = MMS_FIELDS[product]
    vectors = {stem + name for name in spec["vectors"]}
    energy = {stem + name for name in spec["energy"]}
    return data, _variable_layouts(data, vectors=vectors, energy=energy)


def _read_omni(trange, *, variables, prefix, download, **kwargs):
    if download in (False, "offline", "local"):
        from spacepy import pycdf

        dates = Time(trange).datetime
        start, end = dates[0], dates[1]
        directory = Path(_configured_spedas_dir(prefix)) / "omni_data" / "hro_1min"
        parts = {name: [] for name in variables}
        for month in pd.period_range(start.to_period("M"), end.to_period("M"), freq="M"):
            date = month.start_time
            filename = directory / str(date.year) / f'omni_hro_1min_{date:%Y%m%d}_v01.cdf'
            if not filename.exists():
                logging.getLogger(__name__).warning("Missing OMNI file %s", filename)
                continue
            with pycdf.CDF(str(filename)) as cdf:
                times = Time(np.asarray(cdf["Epoch"])).datetime
                mask = np.asarray((times >= start) & (times <= end))
                for name in variables:
                    field = cdf[name]
                    values = np.asarray(field, dtype=float)
                    if values.ndim != 1 or len(values) != len(times):
                        raise ValueError(f"OMNI HRO field {name!r} must be scalar on Epoch; got {values.shape}")
                    attrs = dict(field.attrs)
                    values = values[mask].copy()
                    if "FILLVAL" in attrs:
                        values[values == attrs["FILLVAL"]] = np.nan
                    parts[name].append({"x": times[mask].to_numpy(), "y": values, "attrs": attrs, "units": attrs.get("UNITS")})
        raw = {
            name: {**items[0], "x": np.concatenate([item["x"] for item in items]), "y": np.concatenate([item["y"] for item in items])}
            for name, items in parts.items() if items
        }
    else:
        import requests
        from pyspedas.projects.omni.config import CONFIG
        from pyspedas.projects.omni.load import load as backend_load

        CONFIG["local_data_dir"] = _configured_spedas_dir(prefix) + "/omni_data/"
        kwargs.setdefault("varnames", list(variables or []))
        remote_url = CONFIG["remote_data_dir"]
        proxies = requests.utils.get_environ_proxies(remote_url)
        raw = backend_load(trange=trange, level="hro", notplot=True, **kwargs) or {}
        if not raw and proxies:
            raw = backend_load(trange=trange, level="hro", notplot=True, **kwargs) or {}
        for variable in raw.values():
            times = np.asarray(variable["x"])
            if np.issubdtype(times.dtype, np.number) and times.size and np.nanmedian(np.abs(times)) > 1e16:
                variable["x"] = Time(pd.to_datetime(times, unit="ns")).datetime
    spec = OMNI_FIELDS["hro"]
    layouts = {name: {"kind": spec["kind"], "time_axis": spec["time_axis"], "columns": [name]} for name in raw}
    return raw, layouts


def read_spedas(
    mission,
    trange,
    *,
    product=None,
    variables=None,
    probe=None,
    level="l2",
    prefix=None,
    download=None,
    variable_layouts=None,
    **kwargs,
):
    """Read one mission product and return ``(data, layouts)``.

    Validate options and dispatch to a mission reader, which returns its data
    and product layouts together. Data records retain the pySPEDAS x/y/v
    format. Apply variable_layouts overrides to the returned layouts.
    """
    key = str(mission).lower()
    if key not in SPEDAS_OPTIONS:
        choices = ", ".join(repr(name) for name in SPEDAS_OPTIONS)
        raise ValueError(f"mission must be one of {choices}")
    options = SPEDAS_OPTIONS[key]
    if download is None:
        download = options["download"]
    if product is None:
        if len(options["products"]) != 1:
            choices = ", ".join(repr(name) for name in options["products"])
            raise ValueError(f"product is required for {key}; choose one of {choices}")
        product = options["products"][0]
    product = str(product).lower()
    if product not in options["products"]:
        choices = ", ".join(repr(name) for name in options["products"])
        raise ValueError(f"unsupported {key} product {product!r}; choose one of {choices}")
    if options["probes"]:
        probe = options["default_probe"] if probe is None else str(probe).lower()
        if probe not in options["probes"]:
            choices = ", ".join(repr(name) for name in options["probes"])
            raise ValueError(f"unsupported {key} probe {probe!r}; choose one of {choices}")
    elif probe is not None:
        raise ValueError(f"{key} does not accept a probe")
    variables = qualify_variables(
        key, product,
        default_variables(key, product, probe) if variables is None else variables,
        probe,
    )
    if key == "themis":
        data, layouts = _read_themis(
            trange, product=product, variables=variables, probe=probe,
            level=level, prefix=prefix, download=download, **kwargs,
        )
    elif key == "mms":
        data, layouts = _read_mms(
            trange, product=product, variables=variables, probe=probe,
            prefix=prefix, download=download, **kwargs,
        )
    elif key == "omni":
        data, layouts = _read_omni(
            trange, variables=variables, prefix=prefix, download=download, **kwargs,
        )

    layouts.update(variable_layouts or {})
    return data, layouts


def parse_spedas(data, layouts):
    """Convert pySPEDAS variable mappings into explicit loader descriptions."""
    variables = {}
    for name, variable in data.items():
        if not isinstance(variable, Mapping) or not {"x", "y"}.issubset(variable):
            raise ValueError(f"Variable {name!r} must contain backend time and values")
        if name not in layouts:
            raise ValueError(f"Variable {name!r} has no product layout")
        layout = layouts[name]
        attrs = variable.get("attrs", {})
        axes = {}
        for key, spec in layout.get("axes", {}).items():
            if key not in variable:
                raise ValueError(f"Variable {name!r} is missing backend axis {key!r}")
            axes[key] = {**spec, "values": variable[key]}
        variables[name] = {
            "time": variable["x"],
            "values": variable["y"],
            "time_axis": layout["time_axis"],
            "kind": layout["kind"],
            "columns": layout["columns"],
            "axes": axes,
            "units": variable.get(
                "units",
                attrs.get("UNITS") if isinstance(attrs, Mapping) else None,
            ),
            "attrs": deepcopy(attrs),
        }
    return variables

