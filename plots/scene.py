"""Stable 3-D scene contract for reports and the quick-look service.

This module adds a small, explicit layer over :mod:`plots.v3d`:

* :class:`SceneSpec` is the serializable request (time range, satellites,
  target mode, layers, frame cap);
* :func:`prepare` returns a cached :class:`SceneBundle` and exposes omitted
  satellites plus space-weather availability;
* :func:`add_target` draws the Earth/mid-tail/far-tail target and plane;
* :func:`to_html` emits an iframe-friendly document that can share one
  ``plotly.min.js`` instead of embedding several megabytes per scene.

The existing ``plots.v3d.prepare_scene`` / ``make_scene`` / ``dumpGeo`` API is
unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

import plotly.graph_objects as go
import plotly.io as pio

from common.config import resolve_output_path
from plots import v3d

__all__ = [
    "SceneBundle",
    "SceneSpec",
    "add_target",
    "prepare",
    "target_point",
    "to_html",
    "write_html",
]

_TARGETS = {
    "earth": (0.0, 0.0, 0.0),
    "mid_tail": (-30.0, 0.0, 0.0),
    "far_tail": (-60.0, 0.0, 0.0),
}


def target_point(mode: str) -> tuple[float, float, float] | None:
    """Return the GSM target point for an observation mode.

    ``other`` and ``unknown`` have no target and return ``None``.
    """
    normalized = str(mode).strip().lower()
    return _TARGETS.get(normalized)


@dataclass(frozen=True)
class SceneSpec:
    """Serializable request for one cached 3-D scene."""

    trange: tuple
    satellites: tuple[str, ...] = ()
    target_mode: str = "earth"
    include_field_lines: bool = False
    include_neutral_sheet: bool = False
    include_reference_planes: bool = True
    use_space_weather: bool = True
    space_weather: dict | None = None
    max_frames: int = 24
    satellite_marker_size: float = 6.0
    quality: str = "quick"

    def __post_init__(self):
        if self.target_mode not in {*_TARGETS, "other", "unknown"}:
            raise ValueError(f"unsupported target_mode: {self.target_mode!r}")
        if self.quality not in {"quick", "standard", "full"}:
            raise ValueError(f"unsupported quality: {self.quality!r}")


@dataclass
class SceneBundle:
    """A prepared :class:`v3d.SceneContext` and its request."""

    context: object
    spec: SceneSpec

    @property
    def omitted_satellites(self) -> tuple[str, ...]:
        return tuple(getattr(self.context, "omitted_satellites", ()))

    @property
    def space_weather_status(self) -> str:
        """Describe where the 3-D boundary model gets solar-wind input.

        ``provided`` means the caller supplied arrays; ``load_sw`` means the
        context used the loader's own fallback-capable path; ``none`` means
        the model will use its internal defaults.
        """
        source = getattr(self.context, "space_weather_source", None)
        if source in {"provided", "load_sw", "none"}:
            return source
        weather = getattr(self.context, "space_weather", {}) or {}
        if not weather:
            return "none"
        expected = len(getattr(self.context, "cache_timestamps", ()))
        usable = sum(
            1 for values in weather.values()
            if expected and len(values) == expected
        )
        if usable == len(weather):
            return "provided"
        if usable:
            return "partial"
        return "none"

    def make_scene(self, time):
        figure = self.context.make_scene(time)
        add_target(figure, self.spec.target_mode)
        return figure

    def make_timerange(self, trange=None, *, max_frames=None):
        figure = self.context.make_timerange(
            trange if trange is not None else self.spec.trange,
            max_frames=max_frames if max_frames is not None else self.spec.max_frames,
        )
        add_target(figure, self.spec.target_mode)
        return figure

    def to_html(self, **kwargs):
        """Render the animated scene and return an HTML document."""
        return to_html(self.make_timerange(), **kwargs)

    def write_html(self, output_path, **kwargs):
        """Render the animated scene and write it atomically."""
        return write_html(self.make_timerange(), output_path, **kwargs)


def prepare(spec: SceneSpec) -> SceneBundle:
    """Preload and cache all data needed by ``spec``.

    ``quality`` selects a conservative layer budget: ``quick`` omits T96 field
    lines and the neutral sheet, ``standard`` adds the neutral sheet, and
    ``full`` honours the explicit ``include_*`` flags.
    """
    include_field_lines = spec.include_field_lines
    include_neutral_sheet = spec.include_neutral_sheet
    if spec.quality == "quick":
        include_field_lines = False
        include_neutral_sheet = False
    elif spec.quality == "standard":
        include_field_lines = False
        include_neutral_sheet = True

    context = v3d.prepare_scene(
        list(spec.trange),
        satellites=spec.satellites,
        use_space_weather=spec.use_space_weather,
        space_weather=spec.space_weather,
        include_field_lines=include_field_lines,
        include_neutral_sheet=include_neutral_sheet,
        include_reference_planes=spec.include_reference_planes,
        satellite_marker_size=spec.satellite_marker_size,
    )
    return SceneBundle(context=context, spec=spec)


def add_target(figure: go.Figure, mode: str, *, plane: bool = True) -> go.Figure:
    """Add the target marker and, for tail modes, its GSM plane."""
    target = target_point(mode)
    if target is None:
        return figure
    colors = {"earth": "#2e8b57", "mid_tail": "#c47f17", "far_tail": "#6b5bb5"}
    color = colors.get(str(mode).strip().lower(), "#6b7280")
    figure.add_trace(
        go.Scatter3d(
            x=[target[0]],
            y=[target[1]],
            z=[target[2]],
            mode="markers+text",
            marker=dict(size=8, color=color, symbol="diamond"),
            text=[str(mode)],
            textposition="top center",
            name=f"{mode} target",
            showlegend=True,
        )
    )
    if plane and target[0] != 0.0:
        figure.add_traces(
            v3d.plane(
                (1.0, 0.0, 0.0, -target[0]),
                x_range=(-80.0, 20.0),
                y_range=(-40.0, 40.0),
                z_range=(-40.0, 40.0),
                color=color,
                opacity=0.06,
                name=f"X={target[0]:g} R_E",
            )
        )
    return figure


def to_html(
    figure: go.Figure,
    *,
    include_plotlyjs: bool = False,
    plotlyjs_url: str = "/static/vendor/plotly.min.js",
    div_id: str = "gena-scene",
) -> str:
    """Return a standalone or shared-Plotly HTML document for ``figure``.

    With ``include_plotlyjs=False`` the document references ``plotlyjs_url``
    once, which is the right choice for a local web service. Set
    ``include_plotlyjs=True`` for a self-contained file.
    """
    body = pio.to_html(
        figure,
        full_html=bool(include_plotlyjs),
        include_plotlyjs=bool(include_plotlyjs),
        div_id=div_id,
        config={"responsive": True, "displaylogo": False},
    )
    if include_plotlyjs:
        return body
    script = (
        f'<script src="{escape(plotlyjs_url, quote=True)}"></script>'
        if plotlyjs_url else ""
    )
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<style>html,body,#{div_id}{{height:100%;margin:0}}</style>"
        f"</head><body>{script}{body}</body></html>"
    )


def write_html(
    figure: go.Figure,
    output_path,
    *,
    include_plotlyjs: bool = False,
    plotlyjs_url: str = "/static/vendor/plotly.min.js",
    div_id: str = "gena-scene",
) -> Path:
    """Atomically write :func:`to_html` output to ``output_path``."""
    destination = resolve_output_path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        temporary.write_text(
            to_html(
                figure,
                include_plotlyjs=include_plotlyjs,
                plotlyjs_url=plotlyjs_url,
                div_id=div_id,
            ),
            encoding="utf-8",
        )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
