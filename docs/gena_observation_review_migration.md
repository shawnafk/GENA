# GENA observation review 功能迁移说明

> 本文记录历史迁移过程。现行入口已于 2026-09-17 核对更新；文末测试次数属于当时记录，本次验证见 [design 修改报告](design_review_20260917.md)。

更新日期：2026-09-14

## 1. 变更目的

原 `gena/post/review.py` 同时承担数据加载、观测分段、指向计算、环境参数绘图、快速浏览图生成和 Excel 输出，并且依赖已经不存在的 `merge` 包和旧版 `load_gena_from()` 接口。该脚本在导入时还会立即加载固定的长时间范围，因此不适合作为当前 GENA 包中的可复用模块。

本次修改删除了 `gena/post/review.py`，将报告和观测汇总功能迁移到 `gena.post.export`，将环境参数与 QQ/FOV 几何的可复用计算和绘图功能迁移到 `gena.utils`。数据加载、几何计算、绘图和文件输出现在具有明确边界。

## 2. 迁移结果

| 原功能 | 新位置 | 新接口 |
| --- | --- | --- |
| 根据 PDF 时间间断划分观测 | `gena.post.export` | `build_observation_summary()` |
| 自动加载指定日期或时间段并生成 Excel | `gena.post.export` | `observation_summary()` |
| QQ 地心轨道和视场中心方向计算 | `gena.utils` | `pointing_geometry()` |
| QQ 轨道及 XY/XZ/YZ 平面交点绘图 | `gena.utils` | `orbit()` |
| OMNI、Dst、AE/AL/AU 等环境参数绘图 | `gena.utils` | `env()` |
| L1A joint signal 图 | `gena.post.export.report()` | 自动保存 `joint_signal_distribution` 图 |
| QQ/FOV 几何报告图 | `gena.post.export.report()` | 自动保存 `orbit_pointing` 图 |

已有的能量分布、FOV、Layer 1/2、start signal 和 end signal 绘图继续复用 `plots.basic.gena`，没有在新模块中重复实现。

## 3. 观测汇总表接口

### 3.1 自动加载并输出 Excel

```python
from common.config import configure
from gena.post import export

configure(file_base="/path/to/data")

summary = export.observation_summary(
    "A",
    ["2026-04-07", "2026-04-09"],
    output_dir="reports",
)

print(summary.attrs["output_path"])
# /absolute/path/to/reports/observation_summary.xlsx
```

函数签名：

```python
observation_summary(
    prob,
    time,
    *,
    output_dir="reports",
    large_gap="1h",
    small_gap="5min",
    coord="GSM",
    targets=None,
    ena_product="l2g",
)
```

`prob` 接受 `A` 或 `B`。默认使用对应探头的严格加空间筛选产品 `l2g_a` 或 `l2g_b` 统计 ENA 事件；如需其他产品，可以显式传入 `ena_product="l2"` 或 `ena_product="l2s"`。

### 3.2 时间输入语义

`observation_summary()` 与 `export.report()` 共用 `normalize_report_times()`，时间输入含义完全一致：

```python
# 一个完整 UTC 日
export.observation_summary("A", "2026-04-07")

# 日期端点：包含首尾日期的连续日序列
export.observation_summary("A", ["2026-04-07", "2026-04-09"])

# 带时钟的两个端点：一个精确时间段，不按午夜拆分
export.observation_summary(
    "A",
    ["2026-04-07T23:50:00Z", "2026-04-08T00:10:00Z"],
)

# 三个或更多日期：离散 UTC 日序列
export.observation_summary(
    "B",
    ["2026-04-07", "2026-04-12", "2026-04-27"],
)
```

所有选定时间的观测记录写入同一个 `observation_summary.xlsx`，工作表名称为 `observations`。

### 3.3 已加载数据接口

如果调用者已经加载 PDF 和 ENA 数据，可以直接构建 DataFrame 而不写文件：

```python
from load import load_gena
from gena.post.export import build_observation_summary

trange = ["2026-04-07T00:00:00Z", "2026-04-08T00:00:00Z"]
pdf = load_gena(trange=trange, product="pdf")
l2g = load_gena(trange=trange, product="l2g_a")

summary = build_observation_summary(pdf, l2g, "A")
```

## 4. 观测分段规则

观测分段分为两层：

1. PDF 相邻样本时间差大于 `large_gap` 时开始一个新的大观测，默认阈值为 1 小时。
2. 每个大观测内部，相邻样本时间差大于 `small_gap` 时开始一个新的连续小观测，默认阈值为 5 分钟。

Excel 中每一行对应一个连续小观测。`parent_l_level_index` 标识所属大观测，`sub_observation_index` 标识该大观测内部的小观测序号。相同大观测的行使用相同背景色，但大观测关系同时保存在数值列中，不依赖颜色表达。

`duration_seconds` 是该小观测中第一条和最后一条 PDF 样本之间的时间跨度，不等同于经过死时间、开关机状态或 good-time interval 修正后的有效曝光时间。`ena_counts_per_second` 是事件数除以该跨度得到的浏览量，不是绝对物理通量。

## 5. QQ 位置、视场方向和观测模式

`gena.utils.pointing_geometry()` 使用 PDF 中的 `MJ2000_Star_X/Y/Z` 和 `MJ2000_q1/q2/q3`：

1. 将月心相对的 QQ ECI2000 位置与 JPL 月球位置相加，得到地心 ECI2000 位置。
2. 将位置转换到指定坐标系，默认是 GSM，空间单位为地球半径 `R_E`。
3. 使用对应探头坐标系中沿 `-Z` 的标称 source-looking 视场中心方向。
4. 只保留沿视场前向延伸的交点，计算中心线与 XY、XZ 和 YZ 三个坐标面的交点。

```python
from gena import utils as gut
from load import load_gena

pdf = load_gena(trange=["2026-04-07", "2026-04-08"], product="pdf")
geometry = gut.pointing_geometry(pdf, probe="A", coord="GSM")
```

返回表包含：

- `spacecraft_x/y/z`：QQ 地心位置。
- `look_x/y/z`：归一化 source-looking 视场中心方向。
- `xy_x/y/z`：视场中心线与 `Z=0` 平面的前向交点。
- `xz_x/y/z`：视场中心线与 `Y=0` 平面的前向交点。
- `yz_x/y/z`：视场中心线与 `X=0` 平面的前向交点。

观测模式使用每个小观测内 `Y=0` 平面交点的平均 `(X, Z)`，选择距离最近的目标中心。默认目标中心为：

| 模式 | GSM `Y=0` 平面目标中心 `(X, Z)` |
| --- | --- |
| `earth` | `(0, 0) R_E` |
| `mid_tail` | `(-30, 0) R_E` |
| `far_tail` | `(-60, 0) R_E` |

默认值是任务级标称目标，不是磁层物理边界。汇总表保留 `target_distance_re`，便于检查实际指向与标称中心的偏离。任务定义变化时应显式覆盖：

```python
targets = {
    "earth": (0.0, 0.0),
    "mid_tail": (-25.0, 0.0),
    "far_tail": (-65.0, 0.0),
}

summary = export.observation_summary(
    "A",
    ["2026-04-07", "2026-04-09"],
    targets=targets,
)
```

如果视场中心线与 GSM `Y=0` 平面没有有效前向交点，`target_mode` 为 `unknown`，投影及距离字段为缺失值。

## 6. observation_summary.xlsx 字段

| 字段 | 含义 |
| --- | --- |
| `observation_index` | 所有小观测的全局顺序编号 |
| `parent_l_level_index` | 所属大观测编号 |
| `sub_observation_index` | 大观测内部的小观测编号 |
| `probe` | 探头 A 或 B |
| `start_time`、`end_time` | 小观测首末 PDF 样本时间，UTC |
| `duration_seconds` | 首末样本的时间跨度，秒 |
| `spacecraft_x/y/z_mean` | 小观测内 QQ GSM 平均位置，`R_E` |
| `pointing_x_min/max/mean` | 视场中心线与 GSM `Y=0` 平面交点的 X 统计，`R_E` |
| `pointing_z_min/max/mean` | 视场中心线与 GSM `Y=0` 平面交点的 Z 统计，`R_E` |
| `target_mode` | 最近的标称观测目标，或 `unknown` |
| `target_distance_re` | 平均投影点到所选目标中心的距离，`R_E` |
| `earth_centered` | `target_mode == "earth"` 的兼容布尔列 |
| `ena_product` | 实际统计的 ENA 产品，例如 `l2g_a` |
| `ena_total_counts` | 小观测首末时间内包含端点的 ENA 事件数 |
| `ena_counts_per_second` | `ena_total_counts / duration_seconds`，仅用于快速浏览 |

工作簿将首行设为筛选表头并冻结，时间保持为可排序的 Excel 日期值，数值字段保持数值类型。

## 7. 环境参数绘图

`plots.advance.env()` 是纯绘图函数，不负责加载或下载 OMNI、Dst、AE、AL、AU 等数据。调用者将已加载的 `TimeSeriesDataFrame` 传入；函数按参数顺序查找包含目标列的第一个数据源。

```python
from gena import utils as gut
from indices.load import load as load_indices
from omni.load import load as load_omni

trange = ["2025-12-04T00:00:00Z", "2025-12-05T03:00:00Z"]
omni = load_omni(trange, download=False)
indices = load_indices(trange)
dst = load_indices(trange, "dst")

fig, axes = gut.env(
    omni,
    indices,
    dst,
    cols=[
        "dst",
        ("ae", "al", "au"),
        "flow_speed",
        "proton_density",
        "BZ_GSM",
        "BY_GSM",
    ],
)
```

`cols` 决定 axes 的形状：一个字符串产生一个面板，一个字符串元组或列表将多条曲线画到同一面板。传入外部 `f` 和 `ax` 时，axes 数量必须与 `cols` 生成的面板数量一致。函数返回 `(fig, ax)`；单面板返回一个 `Axes`，多面板返回 axes 列表。

`report()` 可以接收预加载环境数据并将环境图保存到报告目录：

```python
artifacts = export.report(
    "A",
    trange,
    environment_sources=[omni, indices, dst],
    environment_cols=[
        "dst",
        ("ae", "al", "au"),
        "flow_speed",
        "BZ_GSM",
        "BY_GSM",
    ],
)
```

该接口不会隐式调用 OMNI 或指数下载器。离线报告是否包含环境图完全由调用者提供的数据决定。

## 8. QQ 轨道和坐标面交点图

`gena.utils.orbit()` 默认生成四个时间序列面板：QQ 的 GSM X/Y/Z 轨道、FOV 与 XY 平面交点、FOV 与 XZ 平面交点、FOV 与 YZ 平面交点。

```python
fig, axes = gut.orbit(pdf, probe="A", coord="GSM")
```

已经计算过几何结果时，可以避免重复坐标转换：

```python
geometry = gut.pointing_geometry(pdf, probe="A", coord="GSM")
fig, axes = gut.orbit(geometry=geometry)
```

与 `env()` 一样，`orbit()` 支持通过 `cols` 自定义面板和曲线组合，也支持调用者传入 `f` 和 `ax`。

## 9. report 新增输出

当数据列可用时，`export.report()` 新增以下静态 JPEG：

- `{prob}_orbit_pointing.jpg`：QQ 位置及视场中心线与三个坐标面的交点。
- `{prob}_joint_signal_distribution.jpg`：L1A 的 X/Y start/end 两端联合分布。
- `{prob}_environment.jpg`：仅在传入 `environment_sources` 时生成。

原有 summary、数据级计数、工程参数、start/end signal、能谱、Layer 1/2 和 FOV 图保持不变。轨道图只有在 PDF 包含完整位置和姿态列时才生成，缺少几何列不会阻断其他报告图。

## 10. 删除和兼容性说明

- `gena/post/review.py` 已删除。
- 源码、测试和 GENA 文档中不存在对 `gena.post.review`、`post.review` 或 `export_observation_excel` 的残留引用。
- 原脚本误将 `l2s` 称为 strict；新汇总默认使用当前产品语义中的严格加空间筛选产品 `l2g`。
- 新接口不在模块导入阶段加载数据，也不使用固定的 2017–2026 时间范围。
- `build_observation_summary()` 保留纯内存计算边界，`observation_summary()` 才负责加载和写文件。

## 11. 验证结果

本次迁移的定向验证包括：

- A/B 两个探头的观测分段和 earth/mid-tail/far-tail 分类。
- 单日、连续日期、精确时间段和离散日期时间语义。
- `l2g_a/b` 事件计数和计数率。
- QQ/FOV 的 XY、XZ、YZ 前向交点计算。
- `env()` 的动态 axes、同面板多变量和外部 axes 数量检查。
- `orbit()` 对预计算几何数据的复用。
- `report()` 的 orbit 和 joint JPEG 输出。
- Excel 表头、日期类型、冻结窗格、筛选和渲染检查。
- 修改模块的 `py_compile` 和导入 smoke test。

直接相关的 31 项测试通过，另有 29 项时间序列接口测试通过。完整仓库测试仍存在与本次迁移无关的既有收集问题，包括缺失的 `1.Paper_BBF` 历史脚本、缺失的旧 `ena_drift`/`libs.distribution` 模块、`libs/flux.py` 的既有缩进错误，以及 `plots.basic.gena` 未重导出 `plots.vis.tplot` 的既有测试失败。

本次几何与工作簿验证使用合成数据完成。默认目标中心需要结合实际任务规划定义确认；真实遥测上的指向模式、姿态异常处理和有效曝光时间仍应使用任务数据单独核验。
