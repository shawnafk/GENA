import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist-min";
import type { Attitude, Environment, Scene, Vec, Grid } from "./types";
import { add, scale } from "./geometry";
export type PlotHandle = HTMLDivElement;
export default function SceneView({
  attitude: a,
  environment: e,
  scene,
  layers,
  view,
  onReady,
}: {
  attitude: Attitude | null;
  environment: Environment | null;
  scene: Scene;
  layers: Record<string, boolean>;
  view: string;
  onReady: (v: PlotHandle) => void;
}) {
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (ref.current) onReady(ref.current);
    return () => {
      if (ref.current) Plotly.purge(ref.current);
    };
  }, []);
  useEffect(() => {
    if (!ref.current) return;
    const traces: Plotly.Data[] = [];
    const line = (
      points: Vec[],
      name: string,
      color: string,
      width = 3,
      dash = "solid",
    ) =>
      traces.push({
        type: "scatter3d",
        mode: "lines",
        x: points.map((v) => v[0]),
        y: points.map((v) => v[1]),
        z: points.map((v) => v[2]),
        name,
        line: { color, width, dash: dash as "solid" },
        hoverinfo: "name",
      });
    const surface = (g: Grid, name: string, color: string, opacity: number) =>
      traces.push({
        type: "surface",
        x: g[0] as number[][],
        y: g[1] as number[][],
        z: g[2] as number[][],
        name,
        opacity,
        colorscale: [
          [0, color],
          [1, color],
        ],
        showscale: false,
        showlegend: true,
        hoverinfo: "name",
        connectgaps: false,
      });
    const sphere = (c: Vec, r: number, name: string, color: string) => {
      const g: Grid = [[], [], []];
      for (let i = 0; i < 25; i++) {
        g.forEach((d) => d.push([]));
        for (let j = 0; j < 40; j++) {
          const t = (i * Math.PI) / 24,
            p = (j * 2 * Math.PI) / 39;
          const v = add(c, [
            r * Math.sin(t) * Math.cos(p),
            r * Math.sin(t) * Math.sin(p),
            r * Math.cos(t),
          ]);
          v.forEach((x, k) => g[k][i].push(x));
        }
      }
      surface(g, name, color, 1);
    };
    sphere([0, 0, 0], 1, "地球", "#458bc9");
    if (e) {
      if (layers.mp) surface(e.magnetopause, "磁层顶", "#30d6c5", 0.16);
      if (layers.bs) surface(e.bowshock, "弓激波", "#8799ff", 0.12);
      if (layers.sheet && e.sheet) {
        surface(e.sheet, "T96 中性片", "#e9a45b", 0.38);
        for (const sign of [-1, 1]) {
          const g: Grid = [
            e.sheet[0],
            e.sheet[1],
            e.sheet[2].map((row) =>
              row.map((v) =>
                v === null ? null : v + sign * scene.half_thickness,
              ),
            ),
          ];
          surface(g, "假设片层厚度", "#d99150", 0.09);
        }
      }
      if (layers.field)
        e.field_lines.forEach((l, i) =>
          line(l, "T96 磁力线 " + (i + 1), "#7299bd", 1),
        );
    }
    if (a) {
      sphere(a.moon, 0.2724, "月球", "#bcc7d5");
      if (layers.orbit) line(a.track, "鹊桥轨道", "#b0becf", 2);
      traces.push({
        type: "scatter3d",
        mode: "markers",
        x: [a.position[0]],
        y: [a.position[1]],
        z: [a.position[2]],
        name: "鹊桥",
        marker: { size: 5, color: "#fff" },
      });
      const rayLength = a.target
        ? Math.min(
            100,
            Math.hypot(...a.target.map((v, i) => v - a.position[i])),
          )
        : 75;
      if (layers.reference)
        line(
          [a.position, add(a.position, scale(a.reference_axes[2], -rayLength))],
          "参考 −Z",
          "#7a8a9e",
          3,
          "dash",
        );
      if (layers.measured && a.measured_axes)
        line(
          [a.position, add(a.position, scale(a.measured_axes[2], -rayLength))],
          "实测 −Z",
          "#fb80b1",
          4,
        );
      a.axes.forEach((v, i) =>
        line(
          [a.position, add(a.position, scale(v, 12))],
          ["QQ +X", "QQ +Y", "QQ +Z"][i],
          ["#ff6b77", "#51dbad", "#719aff"][i],
          4,
        ),
      );
      line(
        [a.position, add(a.position, scale(a.axes[2], -rayLength))],
        "候选 −Z",
        "#ffbc6c",
        5,
      );
      line(
        [a.position, add(a.position, scale(a.target_direction, rayLength))],
        "目标方向",
        "#ff70b5",
        3,
        "dash",
      );
      if (layers.sun)
        line(
          [a.position, add(a.position, scale(a.sun_direction, 25))],
          "太阳方向",
          "#f7e18a",
          3,
        );
      if (a.target)
        traces.push({
          type: "scatter3d",
          mode: "markers",
          x: [a.target[0]],
          y: [a.target[1]],
          z: [a.target[2]],
          name: "目标",
          marker: { size: 5, color: "#ff70b5", symbol: "diamond" },
        });
      if (layers.fov)
        for (let k = 0; k < 2; k++) {
          const corners = a.fov
            .slice(k * 4, k * 4 + 4)
            .map((v) => add(a.position, scale(v, rayLength)));
          const points: Vec[] = [];
          for (let i = 0; i < 4; i++)
            points.push(
              a.position,
              corners[i],
              corners[(i + 1) % 4],
              a.position,
            );
          line(points, k ? "B 视场" : "A 视场", k ? "#93a0ff" : "#5be0bc", 2);
        }
    }
    const cameras: Record<string, Partial<Plotly.Camera>> = {
      perspective: { eye: { x: 1.3, y: 1.4, z: 0.9 } },
      xy: { eye: { x: 0, y: 0, z: 2.2 }, up: { x: 0, y: 1, z: 0 } },
      xz: { eye: { x: 0, y: -2.2, z: 0 }, up: { x: 0, y: 0, z: 1 } },
      yz: { eye: { x: 2.2, y: 0, z: 0 }, up: { x: 0, y: 0, z: 1 } },
    };
    let camera = cameras[view];
    if (view === "satellite" && a) {
      const center = scale(a.position, 1 / 80),
        look = scale(a.axes[2], -0.25);
      camera = {
        eye: { x: center[0], y: center[1], z: center[2] },
        center: {
          x: center[0] + look[0],
          y: center[1] + look[1],
          z: center[2] + look[2],
        },
        up: { x: a.axes[1][0], y: a.axes[1][1], z: a.axes[1][2] },
      };
    }
    const axis = (title: string) => ({
      title: { text: title, font: { color: "#879bb5", size: 11 } },
      range: [-80, 80],
      gridcolor: "#1b2b3d",
      zerolinecolor: "#34435b",
      showbackground: false,
      color: "#7c90aa",
    });
    Plotly.react(
      ref.current,
      traces,
      {
        paper_bgcolor: "#0b1420",
        plot_bgcolor: "#0b1420",
        font: { color: "#bdcce0", family: "Inter, sans-serif" },
        margin: { l: 0, r: 0, t: 0, b: 0 },
        showlegend: false,
        uirevision: view,
        scene: {
          aspectmode: "cube",
          xaxis: axis("GSM X · Rᴇ"),
          yaxis: axis("GSM Y · Rᴇ"),
          zaxis: axis("GSM Z · Rᴇ"),
          camera,
        },
        annotations: [
          {
            text: "GSM / Rᴇ",
            xref: "paper",
            yref: "paper",
            x: 0.02,
            y: 0.98,
            showarrow: false,
            font: { color: "#6f849f", size: 11 },
          },
        ],
      },
      {
        responsive: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: "png",
          filename: "middle_20_21_3d",
          width: 1600,
          height: 1000,
          scale: 1,
        },
      },
    );
  }, [a, e, scene.half_thickness, layers, view]);
  useEffect(() => {
    const observer = new ResizeObserver(() => {
      if (ref.current) Plotly.Plots.resize(ref.current);
    });
    if (ref.current) observer.observe(ref.current);
    return () => observer.disconnect();
  }, []);
  return <div className="plot" ref={ref} data-testid="scene-plot" />;
}
