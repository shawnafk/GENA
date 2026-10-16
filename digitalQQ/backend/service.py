"""Worker entrypoints, historical weather adapter and reproducible exports."""

from datetime import timedelta
from io import BytesIO
import csv
import io
import json
import zipfile
import numpy as np
import pandas as pd
from .models import Scene, ExportRequest, Weather
from .physics import attitude, environment, json_safe


def compute(kind, payload):
    scene = Scene.model_validate(payload)
    if kind == "attitude":
        return attitude(scene)
    if kind == "solve":
        return attitude(scene, force_solve=True)
    if kind == "boundaries":
        return json_safe(environment(scene, False))
    if kind == "environment":
        return json_safe(environment(scene))
    raise ValueError("Unknown computation")


def load_weather(iso, download=False):
    """Nearest complete OMNI sample within five minutes, never fill missing values."""
    from load import load_indices, load_spedas

    time = pd.Timestamp(iso)
    if time.tz is None:
        raise ValueError("时间必须包含时区")
    if time > pd.Timestamp.now(tz="UTC"):
        raise ValueError("未来时刻没有 OMNI 观测数据，请使用手动场景")
    start = (time - pd.Timedelta("5min")).isoformat()
    stop = (time + pd.Timedelta("5min")).isoformat()
    mapping = {
        "Pressure": "pressure",
        "BY_GSM": "by",
        "BZ_GSM": "bz",
        "Beta": "beta",
        "Mgs_mach_num": "mach",
    }
    raw = load_spedas("omni", [start, stop], "hro", variables=list(mapping), download=download)
    if raw is None or raw.empty:
        raise ValueError("没有可用的 OMNI 历史数据；可检查数据目录或选择允许下载")
    if any(k not in raw for k in mapping):
        raise ValueError("OMNI 缺少动压、IMF、beta 或磁声马赫数，不能填充完整场景")
    valid = raw.replace([np.inf, -np.inf], np.nan).dropna(subset=list(mapping)).copy()
    # OMNI fill values are positive sentinels well outside physically usable ranges.
    valid = valid[
        (valid.Pressure > 0)
        & (valid.Pressure <= 100)
        & (valid.Beta >= 0)
        & (valid.Beta <= 100)
        & (valid.Mgs_mach_num > 1)
        & (valid.Mgs_mach_num <= 100)
        & (valid.BY_GSM.abs() <= 100)
        & (valid.BZ_GSM.abs() <= 100)
    ]
    if valid.empty:
        raise ValueError("OMNI 样本缺失或含无效填充值")
    offsets = (pd.to_datetime(valid.Datetime, utc=True) - time).abs()
    row = valid.iloc[int(np.argmin(offsets.to_numpy()))]
    observed = (
        pd.Timestamp(row.Datetime).tz_localize("UTC")
        if pd.Timestamp(row.Datetime).tz is None
        else pd.Timestamp(row.Datetime)
    )
    if abs((observed - time).total_seconds()) > 300:
        raise ValueError("最近完整 OMNI 样本超过五分钟，未自动填充")
    # The existing Kyoto Dst loader has no download switch. Preflight its exact
    # monthly cache files so the UI's local-only choice cannot initiate a fetch.
    if not download:
        from common.config import get_configured_paths
        from pathlib import Path

        folder = Path(get_configured_paths()["spedas_base"])
        days = [time - pd.Timedelta("1h"), time + pd.Timedelta("1h")]
        missing = [
            str(folder / f'dst{t.strftime("%y%m")}.for.request')
            for t in days
            if not (folder / f'dst{t.strftime("%y%m")}.for.request').exists()
        ]
        if missing:
            raise ValueError("缺少本地 Dst 文件；允许下载或手动输入 Dst：" + missing[0])
    dst = load_indices(
        [
            (time - pd.Timedelta("1h")).isoformat(),
            (time + pd.Timedelta("1h")).isoformat(),
        ],
        product="dst",
    )
    if dst is None or dst.empty or "dst" not in dst:
        raise ValueError("缺少 Dst 历史观测，未自动填充")
    valid_dst = dst[np.isfinite(dst.dst) & (dst.dst.abs() < 1000)]
    if valid_dst.empty:
        raise ValueError("Dst 数据无有效样本")
    delta = (pd.to_datetime(valid_dst.Datetime, utc=True) - time).abs()
    i = int(np.argmin(delta.to_numpy()))
    if delta.iloc[i] > pd.Timedelta("1h"):
        raise ValueError("Dst 最近样本超过一小时")
    data = {v: float(row[k]) for k, v in mapping.items()}
    data.update(
        dst=float(valid_dst.iloc[i].dst),
        source="OMNI HRO + Kyoto Dst",
        observed_at=observed.isoformat(),
    )
    return Weather(**data).model_dump()


def export_report(payload):
    req = ExportRequest.model_validate(payload)
    from plotly import graph_objects as go

    scene = req.scene
    midpoint = req.start + (req.stop - req.start) / 2
    times = pd.date_range(req.start, req.stop, freq=f"{req.cadence_minutes}min").union(
        pd.DatetimeIndex([midpoint, req.stop])
    )
    rows = []
    midpoint = req.start + (req.stop - req.start) / 2
    for time in times:
        current = scene.model_copy(update={"time": time.to_pydatetime()})
        a = attitude(current)
        row = {
            "beijing_time": time.tz_convert("Asia/Shanghai").isoformat(),
            "utc_time": time.isoformat(),
            "theta_roll_deg": a["roll"],
            "theta_pitch_deg": a["pitch"],
            "pointing_error_deg": a["pointing_error"],
            "sun_angle_after_deg": a["sun_angle"],
            "middle_target_ray_miss_RE": a["ray_miss"],
            "sun_avoidance_gt45": a["sun_ok"],
        }
        row.update(
            {f"qq_to_eci2000_quat_{k}": v for k, v in zip("xyzw", a["quaternion"])}
        )
        row.update({f"position_gsm_{k}_RE": v for k, v in zip("xyz", a["position"])})
        rows.append(row)
    table = pd.DataFrame(rows)
    a = attitude(scene)
    env = environment(scene)
    fig = go.Figure()
    colors = {"magnetopause": "#33c9c1", "bowshock": "#7b8bff", "sheet": "#f5ab55"}
    for name in colors:
        grid = env.get(name)
        if grid:
            fig.add_trace(
                go.Surface(
                    x=grid[0],
                    y=grid[1],
                    z=grid[2],
                    opacity=0.22,
                    showscale=False,
                    name=name,
                    colorscale=[[0, colors[name]], [1, colors[name]]],
                )
            )
    p = np.array(a["position"])
    for vec, name, color in [
        (-np.array(a["axes"])[2], "QQ -Z", "#ffb454"),
        (np.array(a["target_direction"]), "Target", "#e85eab"),
    ]:
        line = np.array([p, p + 70 * vec])
        fig.add_trace(
            go.Scatter3d(
                x=line[:, 0],
                y=line[:, 1],
                z=line[:, 2],
                name=name,
                line={"color": color},
            )
        )
    for i, (name, color) in enumerate([("A FOV", "#49dec1"), ("B FOV", "#8b9bff")]):
        corners = np.array(a["fov"])[4 * i : 4 * i + 4] * 70 + p
        points = []
        for j in range(4):
            points.extend([p, corners[j], corners[(j + 1) % 4], p])
        v = np.array(points)
        fig.add_trace(
            go.Scatter3d(
                x=v[:, 0],
                y=v[:, 1],
                z=v[:, 2],
                mode="lines",
                name=name,
                line={"color": color},
            )
        )
    for l in env["field_lines"]:
        l = np.array(l)
        fig.add_trace(
            go.Scatter3d(
                x=l[:, 0],
                y=l[:, 1],
                z=l[:, 2],
                mode="lines",
                showlegend=False,
                line={"color": "#8bb5d5", "width": 1},
            )
        )
    # Earth and Moon spheres retain physical radius.
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 32), np.linspace(0, np.pi, 20))
    for center, radius, name, color in [
        (np.zeros(3), 1, "Earth", "#3179ba"),
        (np.array(a["moon"]), 1737.4 * 1000 / 6378100, "Moon", "#b5bfcc"),
    ]:
        fig.add_trace(
            go.Surface(
                x=center[0] + radius * np.cos(u) * np.sin(v),
                y=center[1] + radius * np.sin(u) * np.sin(v),
                z=center[2] + radius * np.cos(v),
                showscale=False,
                name=name,
                colorscale=[[0, color], [1, color]],
            )
        )
    fig.update_layout(
        title=f"digitalQQ | {scene.time.isoformat()}",
        scene={
            "aspectmode": "data",
            "xaxis_title": "GSM X (RE)",
            "yaxis_title": "GSM Y (RE)",
            "zaxis_title": "GSM Z (RE)",
        },
    )
    report = BytesIO()
    metadata = {
        "scene": scene.model_dump(mode="json"),
        "interval": {
            "start": req.start.isoformat(),
            "stop": req.stop.isoformat(),
            "cadence_minutes": req.cadence_minutes,
            "stop_is_boundary": True,
        },
        "coordinate": "GSM / RE",
        "quaternion": "QQ -> ECI2000 xyzw",
        "ephemeris": "DE432s geometric, TDB",
        "igrf": "IGRF14",
        "models": env["models"],
        "model_errors": env["errors"],
        "orbit_source": a["source"],
    }
    from .render import figures

    images = figures(table, a, env)
    with zipfile.ZipFile(report, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in images.items():
            z.writestr(name, data)
        z.writestr("scene.json", scene.model_dump_json(indent=2))
        z.writestr(
            "calculation_metadata.json",
            json.dumps(metadata, ensure_ascii=False, indent=2),
        )
        z.writestr(
            "parameters_10min.csv" if req.cadence_minutes == 10 else "parameters.csv",
            table.to_csv(index=False),
        )
        z.writestr(
            "parameters_midpoint.csv",
            pd.DataFrame(
                [
                    table.iloc[
                        int(np.argmin(abs(pd.to_datetime(table.utc_time) - midpoint)))
                    ]
                ]
            ).to_csv(index=False),
        )
        z.writestr(
            "parameters_daily_noon.csv",
            table[table.beijing_time.str.contains("T12:00:00")].to_csv(index=False),
        )
        z.writestr("middle_3d.html", fig.to_html(include_plotlyjs=True))
        z.writestr(
            "观测姿态计算说明.md",
            "# digitalQQ 观测规划\n\n"
            + f'轨道来源：{a["source"]}。场景时刻：{scene.time.isoformat()}。参考时刻：{scene.reference_time.isoformat()}。保持方式：{scene.hold}。\n\n固定 GSM Z → X，角度遵循右手定则；四元数为 QQ → ECI2000 xyzw。等离子体片厚度是可视化假设，T96 远尾结果为模型外推。太阳夹角严格大于 45° 作为本工具的检查条件。此结果用于规划，不是正式飞控指令编码。\n\nHTML 几何图取 scene.time；参数表覆盖导出时段并包含终点边界。PNG 同时保存在报告内，也可从网页单独导出。\n',
        )
    return report.getvalue()
