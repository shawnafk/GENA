"""Variable labels, units and display scales shared by plotting modules.

The registry is deliberately a fallback, not a constraint: loader metadata in
``DataFrame.attrs["variables"]`` wins when it provides a label or unit, and any
unknown variable still plots under its resolved column name.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "SUFFIX_INFO",
    "VARIABLE_INFO",
    "VariableInfo",
    "axis_label",
    "describe",
]


@dataclass(frozen=True)
class VariableInfo:
    """Display metadata for one scalar or vector variable."""

    label: str | None = None
    unit: str | None = None
    scale: str = "linear"


VARIABLE_INFO: dict[str, VariableInfo] = {
    # OMNI / solar wind
    "BZ_GSM": VariableInfo(r"$B_z$", "nT"),
    "BY_GSM": VariableInfo(r"$B_y$", "nT"),
    "flow_speed": VariableInfo(r"$V_{sw}$", "km/s"),
    "proton_density": VariableInfo(r"$n_p$", "cm$^{-3}$", "log"),
    "Pressure": VariableInfo(r"$P_{dyn}$", "nPa"),
    "Beta": VariableInfo(r"$\beta$", ""),
    "Mgs_mach_num": VariableInfo(r"$M_{ms}$", ""),
    "Vx": VariableInfo(r"$V_x$", "km/s"),
    # Kyoto / geomagnetic indices
    "AE_INDEX": VariableInfo("AE", "nT"),
    "SYM_H": VariableInfo("SYM-H", "nT"),
    "dst": VariableInfo("Dst", "nT"),
    "ae": VariableInfo("AE", "nT"),
    "al": VariableInfo("AL", "nT"),
    "au": VariableInfo("AU", "nT"),
    "kp": VariableInfo("Kp", ""),
    # Common THEMIS / MMS product names
    "thb_peir_density": VariableInfo(r"$n_i$", "cm$^{-3}$", "log"),
    "thb_peir_avgtemp": VariableInfo(r"$T_i$", "eV", "log"),
    "thb_peir_en_eflux": VariableInfo("Eflux", "eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$", "log"),
    "thb_fgs_btotal": VariableInfo(r"$|B|$", "nT"),
    "thb_fgs_gsm": VariableInfo(r"$B_{GSM}$", "nT"),
    "numberdensity_fast": VariableInfo(r"$n_i$", "cm$^{-3}$", "log"),
    "tempperp_fast": VariableInfo(r"$T_\perp$", "eV", "log"),
    "temppara_fast": VariableInfo(r"$T_\parallel$", "eV", "log"),
    "bulkv_gse_fast": VariableInfo(r"$V_i$", "km/s"),
    "hplus_number_density": VariableInfo(r"$n_{H+}$", "cm$^{-3}$", "log"),
    "heplus_number_density": VariableInfo(r"$n_{He+}$", "cm$^{-3}$", "log"),
    "oplus_number_density": VariableInfo(r"$n_{O+}$", "cm$^{-3}$", "log"),
    "hplus_tperp": VariableInfo(r"$T_{\perp,H+}$", "eV", "log"),
    "heplus_tperp": VariableInfo(r"$T_{\perp,He+}$", "eV", "log"),
    "oplus_tperp": VariableInfo(r"$T_{\perp,O+}$", "eV", "log"),
}


# Suffix rules cover per-probe THEMIS names (thb_*, thc_*, ...) and MMS names
# that carry a product prefix. Exact matches above take precedence.
SUFFIX_INFO: tuple[tuple[str, VariableInfo], ...] = (
    ("_density", VariableInfo(r"$n_i$", "cm$^{-3}$", "log")),
    ("_avgtemp", VariableInfo(r"$T_i$", "eV", "log")),
    ("_en_eflux", VariableInfo("Eflux", "eV cm$^{-2}$ s$^{-1}$ sr$^{-1}$ eV$^{-1}$", "log")),
    ("_velocity_gsm", VariableInfo(r"$V_i$", "km/s")),
    ("_velocity_gse", VariableInfo(r"$V_i$", "km/s")),
    ("_fgs_gsm", VariableInfo(r"$B_{GSM}$", "nT")),
    ("_fgs_btotal", VariableInfo(r"$|B|$", "nT")),
    ("_ion_bulk_velocity_gsm", VariableInfo(r"$V_i$", "km/s")),
    ("_number_density", VariableInfo(r"$n_i$", "cm$^{-3}$", "log")),
    ("_tempperp_fast", VariableInfo(r"$T_\perp$", "eV", "log")),
    ("_temppara_fast", VariableInfo(r"$T_\parallel$", "eV", "log")),
    ("_bulkv_gse_fast", VariableInfo(r"$V_i$", "km/s")),
    ("_tperp", VariableInfo(r"$T_\perp$", "eV", "log")),
    ("_tpara", VariableInfo(r"$T_\parallel$", "eV", "log")),
)


def describe(name: str) -> VariableInfo:
    """Return display metadata for a resolved variable name.

    Lookup order: exact registry entry, then the longest matching suffix rule,
    then an empty :class:`VariableInfo` so callers can fall back to the raw
    column name.
    """
    if not name:
        return VariableInfo()
    exact = VARIABLE_INFO.get(name)
    if exact is not None:
        return exact
    lowered = name.lower()
    for suffix, info in SUFFIX_INFO:
        if lowered.endswith(suffix):
            return info
    return VariableInfo()


def axis_label(info: VariableInfo, fallback: str | None = None) -> str:
    """Format ``label (unit)`` for an axis, falling back to a column name."""
    label = info.label or fallback
    if label and info.unit:
        return f"{label} ({info.unit})"
    if label:
        return label
    return info.unit or ""
