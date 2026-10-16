import { useEffect, useMemo, useRef, useState } from "react";
import Plotly from "plotly.js-dist-min";
import SceneView from "./SceneView";
import {
  initial,
  type Scene,
  type Attitude,
  type Environment,
  type Orbit,
  type Vec,
} from "./types";
import { request, save, Latest, localTime, utc } from "./api";
import { preview, rayHit, scale } from "./geometry";

function NumberField({
  label,
  value,
  onChange,
  step = 0.1,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
  step?: number;
}) {
  return (
    <label className="field">
      <span>{label}</span>
      <input
        aria-label={label}
        type="number"
        value={value}
        step={step}
        onChange={(e) => {
          if (e.target.value !== "" && Number.isFinite(+e.target.value))
            onChange(+e.target.value);
        }}
      />
    </label>
  );
}
function Section({
  title,
  children,
  index,
}: {
  title: string;
  children: React.ReactNode;
  index: string;
}) {
  return (
    <section>
      <h2>
        <span>{index}</span>
        {title}
      </h2>
      {children}
    </section>
  );
}
export default function App() {
  const [scene, setScene] = useState<Scene>(initial),
    [a, setA] = useState<Attitude | null>(null),
    [env, setEnv] = useState<Environment | null>(null);
  const [error, setError] = useState(""),
    [envError, setEnvError] = useState(""),
    [actionError, setActionError] = useState(""),
    [status, setStatus] = useState("正在初始化"),
    [busy, setBusy] = useState(false),
    [envBusy, setEnvBusy] = useState(false),
    [actionBusy, setActionBusy] = useState(false);
  const [view, setView] = useState("perspective"),
    [playing, setPlaying] = useState(false),
    [start, setStart] = useState("2026-10-20T00:00"),
    [stop, setStop] = useState("2026-10-22T00:00"),
    [cadence, setCadence] = useState(10);
  const [layers, setLayers] = useState<Record<string, boolean>>({
    mp: true,
    bs: true,
    field: true,
    sheet: true,
    orbit: true,
    fov: true,
    reference: true,
    measured: true,
    sun: false,
  });
  const [download, setDownload] = useState(false),
    [importFrame, setImportFrame] = useState<Orbit["frame"]>("ECI2000"),
    [importUnit, setImportUnit] = useState<Orbit["unit"]>("km");
  const plot = useRef<HTMLDivElement | null>(null),
    gate = useRef(new Latest()),
    envGate = useRef(new Latest()),
    actionGate = useRef(new Latest());
  const update = (patch: Partial<Scene>) =>
    setScene((s) => ({ ...s, ...patch }));
  const weather = (key: string, v: number) =>
    setScene((s) => ({
      ...s,
      weather: { ...s.weather, [key]: v, source: "manual", observed_at: null },
    }));
  // Increment on every scene change; stale successes and errors cannot replace a newer scene.
  useEffect(() => {
    const id = gate.current.next(),
      controller = new AbortController();
    setBusy(true);
    setError("");
    const timer = setTimeout(
      () =>
        request<Attitude>("/api/attitude", scene, controller.signal)
          .then((r) => {
            if (gate.current.current(id)) {
              setA(r);
              setStatus("精确计算已更新");
              setBusy(false);
            }
          })
          .catch((e) => {
            if (gate.current.current(id) && e.name !== "AbortError") {
              setError(e.message);
              setBusy(false);
            }
          }),
      180,
    );
    return () => {
      clearTimeout(timer);
      controller.abort();
      gate.current.next();
    };
  }, [scene]);
  const environmentKey = JSON.stringify({
    time: scene.time,
    weather: scene.weather,
    magnetic: scene.magnetic && !playing,
    quality: scene.quality,
  });
  useEffect(() => {
    const id = envGate.current.next(),
      controller = new AbortController();
    const environmentScene = { ...scene, magnetic: scene.magnetic && !playing };
    setEnvBusy(true);
    setEnvError("");
    const timer = setTimeout(async () => {
      try {
        const fast = await request<Environment>(
          "/api/environment?part=boundaries",
          environmentScene,
          controller.signal,
        );
        if (!envGate.current.current(id)) return;
        setEnv((prev) => ({
          ...fast,
          field_lines: scene.magnetic ? (prev?.field_lines ?? []) : [],
          sheet: scene.magnetic ? (prev?.sheet ?? null) : null,
          field_time: prev?.field_time,
          sheet_time: prev?.sheet_time,
        }));
        if (environmentScene.magnetic) {
          const full = await request<Environment>(
            "/api/environment",
            environmentScene,
            controller.signal,
          );
          if (!envGate.current.current(id)) return;
          setEnv((prev) => ({
            ...full,
            field_lines: full.errors.some((e) => e.startsWith("T96 磁力线"))
              ? (prev?.field_lines ?? [])
              : full.field_lines,
            sheet: full.errors.some((e) => e.startsWith("T96 中性片"))
              ? (prev?.sheet ?? null)
              : full.sheet,
            field_time: full.errors.some((e) => e.startsWith("T96 磁力线"))
              ? prev?.field_time
              : full.time,
            sheet_time: full.errors.some((e) => e.startsWith("T96 中性片"))
              ? prev?.sheet_time
              : full.time,
          }));
          setEnvError(full.errors.join("；"));
        }
        setEnvBusy(false);
      } catch (e) {
        if (envGate.current.current(id) && (e as Error).name !== "AbortError") {
          setEnvError((e as Error).message);
          setEnvBusy(false);
        }
      }
    }, 300);
    return () => {
      clearTimeout(timer);
      controller.abort();
      envGate.current.next();
    };
  }, [environmentKey]);
  useEffect(() => {
    if (!playing) return;
    const timer = setInterval(() => {
      if (busy) return;
      setScene((s) => {
        const next = new Date(s.time).getTime() + cadence * 60000;
        if (next > new Date(utc(stop)).getTime()) {
          setPlaying(false);
          return s;
        }
        return { ...s, time: new Date(next).toISOString() };
      });
    }, 600);
    return () => clearInterval(timer);
  }, [playing, busy, stop, cadence]);
  const displayed = useMemo(
    () =>
      a && scene.hold === "earth" ? preview(a, scene.roll, scene.pitch) : a,
    [a, scene.roll, scene.pitch, scene.hold],
  );
  const current =
    !!a &&
    new Date(a.time).getTime() === new Date(scene.time).getTime() &&
    !busy &&
    !error;
  const hits = useMemo(() => {
    if (!displayed || !env) return [];
    return [
      ["magnetopause", "磁层顶"],
      ["bowshock", "弓激波"],
      ["sheet", "中性片"],
    ].map(([key, label]) => {
      const grid = env[key as "sheet"];
      return {
        label,
        distance: grid
          ? rayHit(displayed.position, scale(displayed.axes[2], -1), grid)
          : null,
      };
    });
  }, [displayed, env]);
  async function solve() {
    const id = actionGate.current.next(),
      snapshot = JSON.stringify(scene);
    setActionBusy(true);
    setActionError("");
    try {
      const r = await request<Attitude>("/api/solve-pointing", scene);
      if (!actionGate.current.current(id)) return;
      setScene((s) =>
        JSON.stringify(s) === snapshot
          ? {
              ...s,
              roll: r.solution.roll,
              pitch: r.solution.pitch,
              reference_time: s.time,
            }
          : s,
      );
      setStatus(
        r.solution.reachable ? "已求解精确指向" : "目标不可达，显示最小误差解",
      );
    } catch (e) {
      if (actionGate.current.current(id)) setActionError((e as Error).message);
    } finally {
      if (actionGate.current.current(id)) setActionBusy(false);
    }
  }
  async function omni() {
    const id = actionGate.current.next(),
      time = scene.time;
    setActionBusy(true);
    setActionError("");
    try {
      const w = await request<Scene["weather"]>(
        `/api/space-weather?time=${encodeURIComponent(time)}&download=${download}`,
      );
      if (actionGate.current.current(id))
        setScene((s) => (s.time === time ? { ...s, weather: w } : s));
    } catch (e) {
      if (actionGate.current.current(id)) setActionError((e as Error).message);
    } finally {
      if (actionGate.current.current(id)) setActionBusy(false);
    }
  }
  async function importOrbit(file: File) {
    setPlaying(false);
    setActionError("");
    try {
      const text = await file.text();
      let orbit: Orbit;
      if (file.name.endsWith(".json")) orbit = JSON.parse(text);
      else {
        const lines = text
            .trim()
            .replace(/^\uFEFF/, "")
            .split(/\r?\n/),
          header = lines
            .shift()!
            .split(",")
            .map((s) => s.trim());
        const samples = lines.map((line) => {
          const v = line.split(",").map((s) => s.trim()),
            get = (key: string) => v[header.indexOf(key)];
          const q = ["qx", "qy", "qz", "qw"].map((k) => get(k));
          if (
            q.some((v) => v !== undefined && v !== "") &&
            !q.every((v) => v !== undefined && v !== "")
          )
            throw new Error("每行必须提供完整 qx,qy,qz,qw 或全部省略");
          return {
            time: get("time"),
            position: ["x", "y", "z"].map((k) => Number(get(k))) as Vec,
            ...(q.every((v) => v !== undefined && v !== "")
              ? { quaternion: q.map(Number) }
              : {}),
          };
        });
        if (!["time", "x", "y", "z"].every((k) => header.includes(k)))
          throw new Error("CSV 必须含 time,x,y,z；四元数可选 qx,qy,qz,qw");
        orbit = {
          source: "imported",
          frame: importFrame,
          unit: importUnit,
          samples,
        };
      }
      const valid = await request<Orbit>("/api/orbit/validate", orbit);
      if (valid.source !== "imported")
        throw new Error("请导入 source=imported 的轨道");
      update({
        orbit: valid,
        time: valid.samples[0].time,
        reference_time: valid.samples[0].time,
      });
      setStatus(`已导入 ${valid.samples.length} 个轨道样本`);
    } catch (e) {
      setActionError((e as Error).message);
    }
  }
  async function loadScene(file: File) {
    setPlaying(false);
    setActionError("");
    try {
      const candidate = JSON.parse(await file.text());
      const valid = await request<Scene>("/api/scene/validate", candidate);
      setScene(valid);
      setStatus("场景已恢复");
    } catch (e) {
      setActionError((e as Error).message);
    }
  }
  async function exportReport() {
    setActionBusy(true);
    setActionError("");
    try {
      const r = await fetch("/api/export", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          scene,
          start: utc(start),
          stop: utc(stop),
          cadence_minutes: cadence,
        }),
      });
      if (!r.ok) {
        const e = await r.json();
        throw new Error(
          typeof e.detail === "string" ? e.detail : JSON.stringify(e.detail),
        );
      }
      save("digitalQQ_observation.zip", await r.blob(), "application/zip");
    } catch (e) {
      setActionError((e as Error).message);
    } finally {
      setActionBusy(false);
    }
  }
  const timestamp = new Date(scene.time).getTime(),
    begin = new Date(utc(start)).getTime(),
    end = new Date(utc(stop)).getTime();
  return (
    <div className="app">
      <header>
        <div className="brand">
          <span className="brand-icon">◈</span>
          <div>
            <strong>
              digital<span>QQ</span>
            </strong>
            <small>鹊桥 · 空间环境与观测规划</small>
          </div>
        </div>
        <div className="header-meta">
          <span className="tag">LOCAL SIMULATION</span>
          <span className={"signal " + (busy || envBusy ? "pending" : "")} />
          {busy
            ? "姿态计算中"
            : envBusy
              ? "磁场计算中"
              : error || envError
                ? "计算失败"
                : "场景就绪"}
        </div>
        <button
          onClick={() => save("scene.json", JSON.stringify(scene, null, 2))}
        >
          保存场景
        </button>
        <label className="button">
          打开场景
          <input
            hidden
            type="file"
            accept=".json"
            onChange={(e) => {
              if (e.target.files?.[0]) void loadScene(e.target.files[0]);
              e.target.value = "";
            }}
          />
        </label>
      </header>
      <main>
        <aside className="left">
          <Section index="01" title="空间环境">
            <label className="field">
              <span>北京时间 · UTC+8</span>
              <input
                aria-label="北京时间"
                type="datetime-local"
                value={localTime(scene.time)}
                onChange={(e) => {
                  if (e.target.value) update({ time: utc(e.target.value) });
                }}
              />
            </label>
            <div className="fields">
              {[
                ["pressure", "动压 · nPa"],
                ["by", "IMF By · nT"],
                ["bz", "IMF Bz · nT"],
                ["dst", "Dst · nT"],
                ["beta", "等离子体 β"],
                ["mach", "磁声马赫数"],
              ].map(([key, label]) => (
                <NumberField
                  key={key}
                  label={label}
                  value={scene.weather[key as "pressure"]}
                  onChange={(v) => weather(key, v)}
                />
              ))}
            </div>
            <div className="source">
              {scene.weather.source === "manual"
                ? "手动场景"
                : scene.weather.source}
              {scene.weather.observed_at && (
                <small>观测 {localTime(scene.weather.observed_at)}</small>
              )}
            </div>
            <button className="wide" disabled={actionBusy} onClick={omni}>
              读取 OMNI 历史参数 ↗
            </button>
            <label className="check">
              <input
                type="checkbox"
                checked={download}
                onChange={(e) => setDownload(e.target.checked)}
              />
              允许下载缺失历史数据
            </label>
          </Section>
          <Section index="02" title="轨道与参考">
            <label className="field">
              <span>轨道来源</span>
              <select
                aria-label="轨道来源"
                value={scene.orbit.source}
                onChange={(e) => {
                  if (e.target.value === "moon")
                    update({ orbit: { ...scene.orbit, source: "moon" } });
                }}
              >
                <option value="moon">月心近似</option>
                <option
                  value="imported"
                  disabled={scene.orbit.source !== "imported"}
                >
                  导入的真实 / 预测轨道
                </option>
              </select>
            </label>
            <div className="fields">
              <label className="field">
                <span>CSV 坐标系</span>
                <select
                  value={importFrame}
                  onChange={(e) =>
                    setImportFrame(e.target.value as Orbit["frame"])
                  }
                >
                  {["ECI2000", "GSM", "GSE"].map((v) => (
                    <option key={v}>{v}</option>
                  ))}
                </select>
              </label>
              <label className="field">
                <span>CSV 位置单位</span>
                <select
                  value={importUnit}
                  onChange={(e) =>
                    setImportUnit(e.target.value as Orbit["unit"])
                  }
                >
                  {["km", "m", "RE"].map((v) => (
                    <option key={v}>{v}</option>
                  ))}
                </select>
              </label>
            </div>
            <label className="button wide">
              导入轨道 JSON / CSV
              <input
                hidden
                type="file"
                accept=".json,.csv"
                onChange={(e) => {
                  if (e.target.files?.[0]) void importOrbit(e.target.files[0]);
                  e.target.value = "";
                }}
              />
            </label>
            <p className="hint">
              CSV: time,x,y,z；可选 qx,qy,qz,qw。时间必须含时区，四元数为 QQ →
              ECI2000。
            </p>
            <div className="notice">
              {scene.orbit.source === "moon"
                ? "当前使用月心位置近似鹊桥轨道"
                : "真实轨道范围内插值，超出覆盖不外推"}
            </div>
          </Section>
          <Section index="03" title="图层与模型">
            <div className="layer-grid">
              {Object.entries({
                mp: "磁层顶",
                bs: "弓激波",
                field: "T96 磁力线",
                sheet: "中性片",
                orbit: "卫星轨道",
                fov: "A / B 视场",
                reference: "参考姿态",
                measured: "实测姿态",
                sun: "太阳方向",
              }).map(([key, label]) => (
                <label className="check" key={key}>
                  <input
                    type="checkbox"
                    checked={layers[key]}
                    onChange={(e) =>
                      setLayers((v) => ({ ...v, [key]: e.target.checked }))
                    }
                  />
                  {label}
                </label>
              ))}
            </div>
            <NumberField
              label="片层假设半厚度 · Rᴇ"
              value={scene.half_thickness}
              onChange={(v) => update({ half_thickness: v })}
            />
            <label className="check">
              <input
                type="checkbox"
                checked={scene.magnetic}
                onChange={(e) => update({ magnetic: e.target.checked })}
              />
              计算 T96 磁场结构
            </label>
            <label className="field">
              <span>网格精度</span>
              <select
                value={scene.quality}
                onChange={(e) =>
                  update({ quality: e.target.value as Scene["quality"] })
                }
              >
                <option value="preview">预览 · 快速</option>
                <option value="fine">精细 · 后台计算</option>
              </select>
            </label>
            <p className="hint">
              中性片为 T96 + IGRF 的 Bx=0 曲面；厚度为人为假设。远尾为模型外推。
            </p>
          </Section>
        </aside>
        <div className="center">
          <div className="scene-toolbar">
            <div>
              <span className="eyebrow">OBSERVATION GEOMETRY</span>
              <h1>磁层空间 · 姿态预演</h1>
            </div>
            <div className="view-buttons">
              {[
                ["perspective", "3D"],
                ["xy", "XY"],
                ["xz", "XZ"],
                ["yz", "YZ"],
                ["satellite", "卫星"],
              ].map(([k, label]) => (
                <button
                  key={k}
                  className={view === k ? "active" : ""}
                  onClick={() => setView(k)}
                >
                  {label}
                </button>
              ))}
            </div>
          </div>
          <div className="scene-wrap">
            <SceneView
              attitude={displayed}
              environment={env}
              scene={scene}
              layers={layers}
              view={view}
              onReady={(v) => (plot.current = v)}
            />
            <div className="scene-badge">
              GSM <span>•</span> {a?.source || "月心近似"}
              <span>•</span>几何{" "}
              {localTime(a?.time || scene.time).replace("T", " ")} 北京时间
            </div>
            <div className="legend">
              <span style={{ color: "#30d6c5" }}>● 磁层顶</span>
              <span style={{ color: "#8799ff" }}>● 弓激波</span>
              <span style={{ color: "#e9a45b" }}>● 中性片</span>
              <span style={{ color: "#ffbc6c" }}>━ −Z</span>
            </div>
          </div>
          <div className="diagnostics">
            <div>
              <small>目标方向误差</small>
              <strong>
                {current ? a!.pointing_error.toFixed(4) : "—"}
                <em>°</em>
              </strong>
            </div>
            <div>
              <small>太阳夹角</small>
              <strong className={current && !a!.sun_ok ? "bad" : "good"}>
                {current ? a!.sun_angle.toFixed(3) : "—"}
                <em>°</em>
              </strong>
            </div>
            <div>
              <small>目标射线距离</small>
              <strong>
                {current
                  ? a!.ray_miss === null
                    ? "方向目标"
                    : a!.ray_miss.toExponential(2)
                  : "—"}
                <em>{a?.ray_miss === null ? "" : "Rᴇ"}</em>
              </strong>
            </div>
            <div>
              <small>求解状态</small>
              <strong className="compact">
                {current
                  ? a!.solution.reachable
                    ? "目标可达"
                    : "存在残差"
                  : "等待计算"}
              </strong>
            </div>
          </div>
          <div className="messages" aria-live="polite">
            {[error, envError, actionError].filter(Boolean).map((v, i) => (
              <div role="alert" className="error" key={i}>
                {v}
              </div>
            ))}
            {!error && !envError && !actionError && (
              <p>
                {status}
                {envBusy ? " · 磁场结构计算中，磁场图层保留上次结果" : ""}
                {env &&
                new Date(env.time).getTime() !== new Date(scene.time).getTime()
                  ? " · 环境图层来自上一时刻"
                  : ""}
                {busy ? " · 角度拖动为即时预览，诊断等待精确计算" : ""}
              </p>
            )}
          </div>
          <div className="timeline">
            <div className="timeline-top">
              <h2>
                时间回放 <span>BEIJING TIME</span>
              </h2>
              <button
                className="play"
                disabled={!!error || !(end > begin)}
                onClick={() => {
                  if (!playing && (timestamp < begin || timestamp >= end))
                    update({ time: utc(start) });
                  setPlaying(!playing);
                }}
              >
                {playing ? "Ⅱ 暂停" : "▶ 播放"}
              </button>
            </div>
            <input
              aria-label="时间轴"
              type="range"
              min={begin}
              max={Math.max(begin + 1, end)}
              step={cadence * 60000}
              value={Math.max(begin, Math.min(end, timestamp))}
              onChange={(e) => {
                setPlaying(false);
                update({ time: new Date(+e.target.value).toISOString() });
              }}
            />
            <div className="timeline-fields">
              <label>
                起点
                <input
                  aria-label="回放起点"
                  type="datetime-local"
                  value={start}
                  onChange={(e) => {
                    if (e.target.value) setStart(e.target.value);
                  }}
                />
              </label>
              <label>
                终点
                <input
                  aria-label="回放终点"
                  type="datetime-local"
                  value={stop}
                  onChange={(e) => {
                    if (e.target.value) setStop(e.target.value);
                  }}
                />
              </label>
              <label>
                步长 / 分钟
                <input
                  type="number"
                  min="1"
                  max="1440"
                  value={cadence}
                  onChange={(e) =>
                    setCadence(Math.max(1, Math.min(1440, +e.target.value)))
                  }
                />
              </label>
            </div>
          </div>
        </div>
        <aside className="right">
          <Section index="04" title="姿态调整">
            <span className="tag amber">FIXED GSM Z → X</span>
            {(["roll", "pitch"] as const).map((key, i) => (
              <div className="angle-control" key={key}>
                <NumberField
                  label={`theta_${key} · °`}
                  value={scene[key]}
                  onChange={(v) => update({ [key]: v })}
                  step={0.01}
                />
                <input
                  aria-label={`${key} 滑块`}
                  type="range"
                  min="-180"
                  max="180"
                  step=".01"
                  disabled={scene.hold === "tracking"}
                  value={scene[key]}
                  onChange={(e) => update({ [key]: +e.target.value })}
                />
                <small>绕固定 GSM +{i ? "X" : "Z"} 轴 · 右手正方向</small>
              </div>
            ))}
            <label className="field">
              <span>姿态保持方式</span>
              <select
                aria-label="姿态保持方式"
                value={scene.hold}
                onChange={(e) =>
                  update({
                    hold: e.target.value as Scene["hold"],
                    reference_time: scene.time,
                  })
                }
              >
                <option value="earth">固定角度 · 地心跟踪</option>
                <option value="inertial">固定 ECI2000 惯性姿态</option>
                <option value="tracking">持续跟踪目标</option>
              </select>
            </label>
            <small className="hint">
              参考时刻：{localTime(scene.reference_time).replace("T", " ")}{" "}
              北京时间
            </small>
            <button
              className="wide"
              onClick={() => update({ reference_time: scene.time })}
            >
              以当前时刻锁定参考
            </button>
            {scene.hold === "tracking" && a && (
              <p className="hint">
                当前求解角：{a.roll.toFixed(4)}° / {a.pitch.toFixed(4)}°
              </p>
            )}
          </Section>
          <Section index="05" title="观测目标">
            <label className="field">
              <span>目标定义</span>
              <select
                aria-label="目标定义"
                value={scene.target_mode}
                onChange={(e) =>
                  update({
                    target_mode: e.target.value as Scene["target_mode"],
                  })
                }
              >
                <option value="point">指定 GSM 位置</option>
                <option value="far">远磁尾 · GSM −X</option>
                <option value="sheet">T96 中性片 XY 位置</option>
              </select>
            </label>
            {scene.target_mode !== "far" && (
              <div className="target-fields">
                {["X", "Y", "Z"].map((label, i) => (
                  <NumberField
                    key={label}
                    label={`${label} · Rᴇ`}
                    value={scene.target[i]}
                    onChange={(v) => {
                      const target = [...scene.target] as Vec;
                      target[i] = v;
                      update({ target });
                    }}
                  />
                ))}
              </div>
            )}
            {scene.target_mode === "sheet" && (
              <p className="hint">
                Z 由 T96 Bx=0 求解，输入的 Z 不参与目标定义。
              </p>
            )}
            <button
              className="primary wide"
              disabled={actionBusy}
              onClick={solve}
            >
              反解 roll / pitch →
            </button>
            <p className="hint">
              远磁尾固定 pitch=0；不可达时返回最小残差，明确报告误差。
            </p>
          </Section>
          <Section index="06" title="交会与导出">
            <p className="hint">
              候选 −Z 与已绘制网格的首个前向交点；视场边界用于几何预览。
            </p>
            {hits.map((h) => (
              <div className="hit" key={h.label}>
                <span>{h.label}</span>
                <b>
                  {h.distance === null
                    ? "无交点"
                    : h.distance.toFixed(2) + " Rᴇ"}
                </b>
              </div>
            ))}
            {a?.measured_error !== null && a?.measured_error !== undefined && (
              <div className="hit">
                <span>候选与实测 −Z 偏差</span>
                <b>{a.measured_error.toFixed(3)}°</b>
              </div>
            )}
            <button
              className="wide"
              disabled={!current}
              onClick={async () => {
                try {
                  await navigator.clipboard.writeText(
                    JSON.stringify(
                      {
                        time: a!.time,
                        reference_time: scene.reference_time,
                        hold: scene.hold,
                        frame: "GSM Z → X",
                        theta_roll_deg: a!.roll,
                        theta_pitch_deg: a!.pitch,
                        qq_to_eci2000_xyzw: a!.quaternion,
                      },
                      null,
                      2,
                    ),
                  );
                  setStatus("姿态参数已复制");
                } catch {
                  setActionError("剪贴板不可用，请保存场景或导出报告");
                }
              }}
            >
              复制姿态参数
            </button>
            <button
              className="wide"
              disabled={busy || envBusy || !!error || !a}
              onClick={() => {
                if (plot.current)
                  void Plotly.downloadImage(plot.current, {
                    format: "png",
                    filename: "middle_3d",
                    width: 1600,
                    height: 1000,
                  });
              }}
            >
              导出当前 PNG
            </button>
            <button
              className="primary wide"
              disabled={actionBusy}
              onClick={exportReport}
            >
              {actionBusy ? "正在处理…" : "导出观测报告 ZIP"}
            </button>
            <p className="hint">
              ZIP 含场景、逐时 CSV、四元数、交互
              HTML、元数据和说明。规划参数用于预演。
            </p>
          </Section>
        </aside>
      </main>
      <footer>
        <span>
          digitalQQ / 0.1 <span className="muted">· DE432s / IGRF14 / T96</span>
        </span>
        <span>本机观测规划 · UTC 计算 / UTC+8 显示 · 45° 避日检查</span>
      </footer>
    </div>
  );
}
