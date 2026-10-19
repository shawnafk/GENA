# 绘图模块

仓库级绘图入口为 `plots`。`basic` 按数据源集中组织产品绘图，`advance` 负责跨数据源联合绘图，`plots.vis` 保留画布、坐标轴、colorbar 和保存图片等周边工具。

```text
plots/
├── __init__.py
├── basic/
│   ├── __init__.py
│   ├── timeseries.py  # tplot 与变量读取辅助函数
│   ├── gena.py
│   ├── themis.py
│   ├── mms.py
│   └── omni.py
├── vis.py
└── advance/
    ├── __init__.py
    ├── joint.py
    ├── ani.py
    └── v3d.py
```

## 单独绘图

固定布局的产品函数使用 `plots.vis.pass_axes`：不提供 `f/ax` 或显式传 `None` 时自动创建；只传 `f` 时在该画布内创建子图；传 `ax` 时复用原轴，并检查轴数量及所属画布。单子图函数收到一个 Axes，多子图函数收到按行展开的列表。调用时推荐使用关键字 `f=`、`ax=`，现有位置参数调用仍可使用。动态行数布局继续使用 `axes_for_columns()`。

```python
from plots.basic import gena, themis, mms, omni

# l2、esa、fpi_data 均由调用者提前加载。
figure, axis = gena.energy(l2)
figure, axis = themis.spectrum(esa, "peir_en_eflux")
figure, axes = mms.fpi(fpi_data)
```

各函数保持原来的参数与返回值；名字相同但产品语义不同的函数通过模块名区分。通用时间序列入口是 `from plots.basic import tplot`，实现位于 `plots/basic/timeseries.py`；变量查找与数据提取辅助函数与它放在一起。各卫星模块只在内部复用，不再重导出 `tplot`。此次迁移保留 THEMIS 绘图中配套的统计计算接口，没有更改单位、归一化和数值算法。动画与交互式三维场景分别使用 `plots.ani`、`plots.v3d`。

## 联合绘图

```python
from plots import Panel, overview
from plots.basic import gena, themis, mms, omni

figure, axes = overview(
    [
        Panel(gena.counts, ([l2],), {"window": "1min"}),
        Panel(themis.fgm_bfield, (fgm,)),
        Panel(themis.spectrum, (esa, "peir_en_eflux")),
        Panel(mms.fpi, (fpi_data,), nrows=3),
        Panel(omni.plot_panel, (omni_data, "BZ_GSM")),
    ],
    trange=["2025-12-04T14:00:00", "2025-12-04T16:00:00"],
    figsize=(12, 16),
)
```

`Panel` 表示一次绘图调用；`kwargs` 应为字典，省略时自动创建空字典。`args` 为位置参数元组，`kwargs` 为绘图参数字典，`nrows` 指定产品占用的行数。单个时间序列或能谱 panel 使用默认一行；FPI 使用三行，HPCA 使用五行，THEMIS 概览的行数应与其 `panels` 参数一致。调度器统一提供 `ax`，函数显式接受 `f` 时也会提供 figure；请勿在 `args` 或 `kwargs` 中再次传入它们。

`overview` 返回 `(figure, axes)`，其中 `axes` 始终为数据坐标轴的扁平列表，不包含 colorbar。支持传入已有的 `f`、`ax`；所有坐标轴必须属于同一张图，数量必须匹配总行数。绘图函数应复用收到的坐标轴。

该入口用于共享 UTC 时间轴的联合图；能量分布直方图、FOV、空间投影等横轴不是时间的图应单独绘制。各数据源保留原始时间采样，`trange` 仅设置最终显示范围，不自动下载、截取、插值或同步数据。需要跨数据源运算时，在绘图前显式完成科学计算与时间同步。自定义 plotter 的数据加载及其他行为由调用者负责。

## 太阳风与地磁指数

```python
from plots import env

figure, axes = env(
    omni_data,
    index_data,
    cols=["dst", ("ae", "al", "au"), "flow_speed", "BZ_GSM"],
)
```

`env` 从传入的数据中按顺序寻找变量，支持多个变量共用一个 panel；不负责加载数据。它保留原返回约定：单行返回一个 axis，多行返回 axes 列表。

## 导入路径迁移

| 原入口 | 新入口 |
|---|---|
| `gena.vis` | `plots.basic.gena` |
| `themis.vis` | `plots.basic.themis` |
| `mms.vis` | `plots.basic.mms` |
| `omni.vis` | `plots.basic.omni` |
| `gena.utils.env` | `plots.advance.env` |

旧 mission `vis.py` 文件和 `gena.utils.env` 入口已移除，不保留兼容转接。历史分析脚本需自行更新导入路径。`indices.vis` 原本只有公共 `tplot` 重导出，也已移除；通用指数时序图直接调用 `plots.basic.tplot`，指数联合图使用 `plots.advance.env`，Dst panel 使用 `plots.basic.omni.dst`。

通用动态布局和列分组工具已迁到 `plots.vis.axes_for_columns`、`plots.vis.column_groups`，供环境图与 GENA 轨道图共享。`gena.utils.orbit`、`pointing_geometry` 及速度计算保持原接口。

## 运行时结果目录

```python
from common.config import configure
from plots import savefig

configure(output_dir="/data/results/run_20260917", mkdir=False, validate=False)
savefig(figure, "figures/overview.png")
```

相对输出路径解析到本次结果目录；绝对路径保持显式覆盖。配置也可通过 `GENA_OUTPUT_DIR` 指定，未指定时使用数据根下的 `GENA_RES`。报告默认写入该目录的 `reports/` 子目录；`dumpGeo()` 默认按日期建子目录。直接调用第三方对象自身的保存方法不受此配置管理。

## 通用时间序列入口（2026-09-20）

```python
from plots.basic import tplot

figure, axis = tplot(data, "density")
tplot(data, "velocity_gsm", ax=existing_axis)
```

`DataFrame.tplot()` 与顶层便捷入口 `plots.tplot` 均委托同一实现。旧 `plots.vis.tplot` 和各卫星模块的 `tplot` 导出已移除，不保留转接。该批迁移保持变量选择、单位、时间裁剪和标量/矢量画线逻辑不变；相关回归为 77 项通过，保留既有 IGRF 有效期和 tight_layout 警告。
