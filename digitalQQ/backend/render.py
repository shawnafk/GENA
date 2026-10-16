"""Standalone scientific PNGs for the downloadable observation report."""

from io import BytesIO
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from .physics import FOV


def png(fig):
    data = BytesIO()
    fig.savefig(data, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return data.getvalue()


def figures(table, a, env):
    dates = (
        __import__("pandas")
        .to_datetime(table.beijing_time, utc=True)
        .dt.tz_convert("Asia/Shanghai")
        .dt.tz_localize(None)
    )
    history, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, key, label in zip(
        axes,
        ["theta_roll_deg", "theta_pitch_deg", "pointing_error_deg"],
        ["theta_roll (deg)", "theta_pitch (deg)", "Target error (deg)"],
    ):
        ax.plot(dates, table[key], color="#2c719b")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    axes[-1].set_xlabel("Beijing time (UTC+8)")
    history.suptitle("digitalQQ | GSM Z then X | Moon proxy or imported orbit")
    history.tight_layout()
    sun, ax = plt.subplots(figsize=(10, 4))
    ax.plot(
        dates,
        table.sun_angle_after_deg,
        color="#97639b",
        label="QQ -Z to spacecraft -> Sun",
    )
    ax.axhline(45, color="red", ls="--", label="45 deg limit")
    ax.set_ylabel("Sun angle (deg)")
    ax.set_xlabel("Beijing time (UTC+8)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    ax.grid(alpha=0.25)
    ax.legend()
    p = np.array(a["position"])
    target = np.array(a["target_direction"])
    length = min(100, np.linalg.norm(np.array(a["target"]) - p)) if a["target"] else 75
    flat, ax2 = plt.subplots(2, 3, figsize=(15, 8))
    static, ax3 = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "3d"})
    for j, basis in enumerate([np.array(a["reference_axes"]), np.array(a["axes"])]):
        lines = [
            (np.array([p, p - length * basis[2]]), "purple", "-"),
            (np.array([p, p + length * target]), "orange", "--"),
            (np.array(a["track"]), "black", "-"),
        ]
        for v, c in zip(basis, ["red", "green", "blue"]):
            lines.append((np.array([p, p + 12 * v]), c, "-"))
        corners = FOV @ basis * length + p
        for k, c in [(0, "green"), (4, "blue")]:
            for v in corners[k : k + 4]:
                lines.append((np.array([p, v]), c, "-"))
            lines.append((corners[[k, k + 1, k + 2, k + 3, k]], c, "-"))
        for idx, axis in [
            (list(range(3)), ax3[j]),
            ([0, 1], ax2[j, 0]),
            ([0, 2], ax2[j, 1]),
            ([1, 2], ax2[j, 2]),
        ]:
            for values, color, style in lines:
                axis.plot(*values[:, idx].T, color=color, ls=style, lw=1.3, alpha=0.75)
            axis.scatter(*np.zeros(len(idx)), color="royalblue", s=25)
            if a["target"]:
                axis.scatter(
                    *np.array(a["target"])[idx], color="orange", marker="*", s=60
                )
            axis.set_xlabel(f'GSM {"XYZ"[idx[0]]} (RE)')
            axis.set_ylabel(f'GSM {"XYZ"[idx[1]]} (RE)')
            axis.set_title(
                "Before: Earth pointing" if j == 0 else "After: candidate attitude"
            )
            if len(idx) == 3:
                axis.set_zlabel("GSM Z (RE)")
                axis.set_box_aspect([1, 1, 1])
                axis.view_init(elev=22, azim=125)
            else:
                axis.set_aspect("equal", adjustable="box")
                axis.grid(alpha=0.25)
    for f in [flat, static]:
        f.suptitle(f'digitalQQ | {a["time"]} | purple: QQ -Z; orange: target')
        f.tight_layout()
    return {
        "parameter_history.png": png(history),
        "sun_avoidance.png": png(sun),
        "middle_2d.png": png(flat),
        "middle_3d.png": png(static),
    }
