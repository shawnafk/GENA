# digitalQQ · 鹊桥观测规划

FastAPI + React/TypeScript + Plotly WebGL 的本机观测几何模拟。支持手动空间环境参数、OMNI 历史填充、经验磁层顶/弓激波、T96 磁力线与中性片、鹊桥轨道、A/B 视场、姿态反解和时间回放。全部应用代码位于本目录，复用上层仓库科学模型，不依赖外部 `prepare_obs`。

## 启动

在 GENA-CODE-Dev 根目录，用已经安装仓库科学依赖的 Python 环境运行：

```bash
python -m pip install -r digitalQQ/requirements.txt
cd digitalQQ/frontend
npm ci
npm run build
cd ../..
./digitalQQ/start.sh
```

浏览器打开 <http://127.0.0.1:8765>。服务仅绑定本机。启动脚本在缺少构建产物时自动安装并构建前端；修改前端后需要重新 `npm run build`，修改 Python 后重启服务。

开发时，在两个终端分别运行后端 `python -m uvicorn digitalQQ.backend.app:app --host 127.0.0.1 --port 8765 --reload` 和前端 `cd digitalQQ/frontend && npm run dev`；Vite 将 `/api` 代理到后端。接口文档位于 `/docs`。

T96 使用已安装 geopack 的本地 IGRF14 系数，禁止其导入时联网更新。若首次安装 geopack 后系数尚未下载，错误信息会给出目标路径；将随应用提供的 `data/igrf14coeffs.txt` 复制到该路径即可。SpacePy 的 IGRF14 独立读取应用内系数文件。IGRF14 系数来源为 NOAA/IAGA，覆盖 2025–2030 预测区间。

## 场景与数据格式

默认场景是北京时间 2026-10-21 00:00 的月心近似、中磁尾 GSM `(-30,0,0) R_E` 目标。`examples/scene_20261020_21.json` 可以通过“打开场景”恢复。时间始终带时区，后端统一为 UTC，网页固定显示 UTC+8，不受浏览器时区影响。

- 位置：地心 ECI2000、GSM 或 GSE，单位 km、m 或 RE；RE=6378100 m。网页所有几何统一为地心 GSM / RE。
- 姿态：单位四元数 `[qx,qy,qz,qw]`，含义 QQ → ECI2000。参考 −Z 指向地心，参考 X 位于 ECI2000 XY 平面；实测姿态单独显示，不作为隐式参考姿态。
- 角度：度，主动右手旋转，先固定 GSM +Z 的 roll，再固定 GSM +X 的 pitch。网页即时几何使用相同顺序，诊断和导出由 Python 精确计算。
- 环境：pressure 为 nPa；By/Bz 为 GSM nT；Dst 为 nT；beta 和磁声马赫数无量纲。手动修改后来源变为 manual。

轨道 CSV 必须含 `time,x,y,z`，可选完整的 `qx,qy,qz,qw`。CSV 坐标系和单位通过页面选择。时间必须严格递增且至少有两行；真实预测或遥测均可导入。四元数要么每行提供，要么每行都省略。

```csv
time,x,y,z,qx,qy,qz,qw
2026-10-20T00:00:00+08:00,-320000,180000,10000,0,0,0,1
2026-10-22T00:00:00+08:00,-360000,120000,16000,0,0,0,1
```

轨道 JSON 使用 `{"source":"imported","frame":"ECI2000","unit":"km","samples":[{"time":"...Z","position":[x,y,z],"quaternion":[qx,qy,qz,qw]}]}`，四元数字段可省略。`examples/synthetic_orbit.json` 是导入格式演示，里面的坐标是合成数据，不是真实鹊桥轨道。位置先逐样本转换到 ECI2000 再线性插值，姿态采用 SLERP；超出覆盖时拒绝计算，需要显式改选月心近似。

## 使用流程

1. 选择时间与轨道来源，输入环境参数，或通过“读取 OMNI 历史参数”填充；允许下载默认关闭，数据目录沿用仓库 `FILE_BASE` 配置。
2. 勾选图层、选择预览或精细网格、修改假设片层半厚度。3D 视角可自由旋转缩放，也可选择 XY/XZ/YZ 或卫星视角。
3. 手动拖动 roll/pitch，或选择指定 GSM 点、远尾方向、中性片 XY 点并反解。中性片目标的 Z 由 Bx=0 求解，不使用输入 Z。不可达时显示最小残余误差，不伪称精确指向。
4. 在时间轴上回放，并选择固定角度地心跟踪、固定参考时刻 ECI2000 惯性姿态、持续跟踪目标。惯性保持使用 reference_time，按钮可以将当前时刻设为参考。
5. 保存/加载自包含场景 JSON、复制带时刻和四元数的姿态参数、下载 PNG 或导出观测报告 ZIP。JSON 内保留输入轨道样本，便于复现。

太阳角使用卫星到太阳的真实几何射线，检查严格 `>45°`。交会距离是候选 −Z 与当前绘制的有限曲面网格的首个前向交点，是网格近似结果；A/B 视场显示安装方向与边界，不代表密度积分或 ENA 辐射模拟。

## 物理范围及历史数据

磁层顶和弓激波直接调用仓库 `libs.space` 的经验公式；磁层顶在 Bz=0 时采用正侧极限以避免原公式该点的分段空隙。曲面在 100 RE 内裁剪，并避开弓激波公式渐近线。边界模型输出有限网格，无效部分显示为空洞。

中性片为 T96 外场 + IGRF 内场的 `Bx=0` 根；片层半厚度是用户假设，不是密度或温度模型推导的真实等离子体片边界。磁力线和片层的远尾结果属于模型外推。模型错误会显示在界面，磁场图层保留上次有效结果；更新中会显示其尚未更新的状态。

OMNI 采用 HRO 历史数据，要求距离选定时间不超过五分钟的完整动压、IMF、beta、马赫数样本；Dst 来自 Kyoto，最近样本不超过一小时。填充值、缺列、过期和未来数据明确报错，不使用共享 `sw_keys` 的默认值回退。Dst loader 没有 download 参数，本应用在禁止下载时预先检查对应本地月度文件。最新上游太阳风数据及传播时延不在此版本内。

## API 与后台计算

| 接口 | 功能 |
|---|---|
| `GET /api/health` | 本机服务状态 |
| `POST /api/scene/validate` | 验证及规范化场景 |
| `POST /api/orbit/validate` | 验证导入轨道 |
| `POST /api/attitude` | 按保持方式计算姿态、视场、太阳角、目标误差 |
| `POST /api/solve-pointing` | 在当前时刻反解目标方向 |
| `POST /api/environment?part=boundaries` | 快速边界网格 |
| `POST /api/environment` | 边界、T96 磁力线与中性片 |
| `GET /api/space-weather?time=...&download=false` | 历史参数填充 |
| `POST /api/export` | 场景与时段报告 ZIP |

姿态与边界使用一个串行计算进程，T96 与导出使用另一个串行进程，保持 geopack 全局状态隔离。环境按时间、参数与精度缓存，最多保留 24 个结果；角度改变不触发环境重算。回放期间更新快速边界并保留上次磁场，暂停后计算当前时刻 T96，避免慢任务积压。前端参数停止变化 300 ms 后请求环境，姿态请求防抖 180 ms；拖动角度期间立即更新几何。请求 generation gate 和 AbortController 同时屏蔽过期成功结果及错误。实际诊断在最新精确结果到达前显示等待状态。

## 报告结构与验证

报告包含 `scene.json`、`calculation_metadata.json`、`parameters_10min.csv`（非十分钟步长则为 parameters.csv）、`parameters_midpoint.csv`、`parameters_daily_noon.csv`、`middle_3d.html`、`middle_3d.png`、`middle_2d.png`、`parameter_history.png`、`sun_avoidance.png` 和 `观测姿态计算说明.md`。逐时表包含终点边界和精确区间中点；HTML/几何 PNG 取导出场景时刻。导出限制七天和最多 2000 个常规采样点，以控制本机计算量。

```bash
python -m pytest digitalQQ/tests -q
python -m compileall -q digitalQQ/backend
cd digitalQQ/frontend
npm test
npm run build
npx playwright install chromium
node browser-check.mjs
```

浏览器检查需要先启动本机服务，截图与检查下载保存在忽略的 `artifacts/`。回归基准为北京时间 2026-10-21 00:00：roll=-30.932073968°、pitch=-1.643111289°，月心近似下目标误差约机器精度。测试使用本地星历和合成轨道，OMNI 缺失与完整样本通过模拟 loader 验证；真实遥测和外部历史下载需在配置实际数据后核验。复制输出是规划参数 JSON，不是已编码的飞控指令。

本机预览基准：姿态约 0.08 秒，边界约 0.04 秒，默认 T96 首次计算约 29 秒，缓存返回小于 0.001 秒。实际用时受参数、模型根和网格精度影响；T96 超时会明确显示。
