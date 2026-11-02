# 绘图模块

仓库级绘图入口为 `plots`。产品绘图按数据源放在扁平模块中，跨数据源组合放在 `plots.joint`，画布/坐标轴/保存等基础设施放在 `plots.tools`。旧的 `plots.basic`、`plots.advance`、`plots.vis` 路径已删除，不再提供转接。

```text
plots/
├── __init__.py       # 公共入口，按需加载子模块
├── theme.py          # 显式 apply_theme()，import 不修改 rcParams
├── labels.py         # 变量标签、单位、显示尺度
├── panels.py         # PanelSpec / plot_panel：声明式单轴面板
├── tools.py          # savefig、pass_axes、axes_for_columns、colorbar 工具
├── timeseries.py     # tplot、变量解析
├── joint.py          # Panel、overview、env、plot_region_code
├── gena.py           # GENA 产品
├── themis.py         # THEMIS 产品
├── mms.py            # MMS FPI/HPCA 产品
├── omni.py           # OMNI/地磁指数面板
├── ani.py            # 动画 frames / encoder / preset
├── scene.py          # 3D 场景契约、目标点、共享 Plotly HTML
└── v3d.py            # Plotly 三维场景实现
```

## 公共入口

```python
import plots
from plots import gena, themis, mms, omni
from plots import Panel, overview, env, plot_region_code, savefig, tplot

plots.theme.apply_theme()   # 应用字体/legend 默认值；import 本身不改全局状态
```

- `plots.tplot`：通用标量/三分量时间序列绘图。
- `plots.Panel` / `plots.overview`：把已加载数据组合到共享 UTC 时间轴。
- `plots.env`：太阳风与地磁指数面板。
- `plots.savefig`：统一保存并关闭 figure。
- `plots.gena` / `themis` / `mms` / `omni`：任务产品绘图。
- `plots.panels`：`PanelSpec` 声明式面板，替代大量“tplot + ylabel”包装函数。
- `plots.labels`：变量标签、单位和建议尺度。
- `plots.ani`：预计算帧 + MP4 编码。
- `plots.scene`：3D 场景请求、目标点/平面、共享 Plotly HTML；`plots.v3d` 保留底层实现。

## 单独绘图

固定布局的产品函数使用 `plots.tools.pass_axes`：不提供 `f/ax` 或显式传 `None` 时自动创建；只传 `f` 时在该画布内创建子图；传 `ax` 时复用原轴，并检查轴数量及所属画布。单子图函数收到一个 Axes，多子图函数收到按行展开的列表。调用时推荐使用关键字 `f=`、`ax=`，现有位置参数调用仍可使用。动态行数布局使用 `plots.tools.axes_for_columns()`。

```python
from plots import gena, themis, mms

# l2、esa、fpi_data 均由调用者提前加载。
figure, axis = gena.energy(l2)
figure, axis = themis.spectrum(esa, "peir_en_eflux")
figure, axes = mms.fpi(fpi_data)

# 只传入当前实际存在的产品；缺失的产品不会被强制要求。
figure, axes = themis.overview(esa=esa, fgm=fgm)
figure, axes = mms.overview(fpi_data=fpi_data, hpca_data=hpca_data)
```

绘图函数只接收已加载数据，不读文件、不下载、不同步数据，也不调用 `plt.show()`。函数返回 `(figure, axes)`；单子图返回一个 Axes，多子图返回扁平列表。缺变量时抛出带可用变量示例的 `KeyError`，不静默画空线。

## 声明式单轴面板

只差一个 ylabel、log 轴、零线或 legend 的图，不再单独写函数：

```python
from plots.panels import PanelSpec, plot_panel

plot_panel(omni_data, PanelSpec(("BZ_GSM",), zero=True))
plot_panel(
    fgm,
    PanelSpec(
        ("fgs_btotal", "fgs_gsm"),
        ylabel=r"$B_{GSM}$ (nT)",
        zero=True,
        legend=True,
        key_kwargs=({"color": "k"}, {}),
    ),
)
```

`PanelSpec` 的字段：`keys`、`ylabel`、`yscale`、`zero`、`legend`、`key_kwargs`、`label`。`plot_panels` 把多个 spec 放到共享 UTC 轴上。

## 联合绘图

```python
from plots import Panel, overview, panels
from plots import gena, themis, mms, omni

figure, axes = overview(
    [
        Panel(gena.counts, ([l2],), {"window": "1min"}),
        Panel(
            panels.plot_panel,
            (fgm, panels.PanelSpec(("fgs_btotal", "fgs_gsm"), zero=True, legend=True)),
        ),
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
| `plots.basic.gena` | `plots.gena` |
| `plots.basic.themis` | `plots.themis` |
| `plots.basic.mms` | `plots.mms` |
| `plots.basic.omni` | `plots.omni` |
| `plots.basic.tplot` | `plots.tplot` |
| `plots.advance.env` | `plots.env` |
| `plots.advance.Panel` / `overview` | `plots.Panel` / `plots.overview` |
| `plots.vis.*` | `plots.tools.*` |

历史分析脚本需要自行更新导入路径。`gena.utils.orbit`、`pointing_geometry` 及速度计算保持原接口。

## 变量标签与单位

`plots.labels` 提供变量名到标签、单位和建议尺度的回退表。加载器 `DataFrame.attrs["variables"]` 中的标签优先；未知变量仍按列名绘制。

```python
from plots.labels import axis_label, describe

info = describe("thb_peir_density")
axis_label(info)        # '$n_i$ (cm$^{-3}$)'
```

`plots.tplot` 在 axes 还没有 ylabel 时自动使用该表，不会覆盖调用者显式设置的标签。

## 动画

`plots.ani` 分为三层：

```python
from plots.ani import AnimationSpec, gena_fov_frames, render_fov

frames = gena_fov_frames(l2, window="5min", max_frames=36)
render_fov(
    frames,
    "fov.mp4",
    spec=AnimationSpec(fps=5, dpi=120, poster=True),
    progress=lambda i, n: print(i, n),
    cancel=lambda: False,
)
```

- `gena_fov_frames` 只做窗口直方图，不创建 Matplotlib figure；
- `render_fov` 负责 MP4、poster、进度、取消和原子替换；
- `dump_fov` / `dumpFov` 保留一键入口，`AnimationDemo3` 保留旧构造方式。

## 三维场景

`plots.scene` 在 `plots.v3d` 之上提供稳定契约：

```python
from plots.scene import SceneSpec, prepare, write_html, target_point

spec = SceneSpec(
    trange=("2025-12-04", "2025-12-05"),
    target_mode="mid_tail",      # earth / mid_tail / far_tail / other
    quality="quick",
    max_frames=24,
    # 可选：显式传入 OMNI 数组，避免 load_sw 的默认值兜底
    # space_weather={"BZ_GSM": ..., "BY_GSM": ..., ...},
)
bundle = prepare(spec)
figure = bundle.make_timerange()
bundle.write_html("scene.html", include_plotlyjs=False,
                  plotlyjs_url="/static/vendor/plotly.min.js")
```

- `SceneSpec` 是可序列化请求；`quality="quick"` 默认不画 T96 场线/中性片；
- `space_weather` 传入时，`bundle.space_weather_status == "provided"`；用 `load_sw` 时为 `"load_sw"`，关闭时为 `"none"`；
- `add_target` 画目标点和 `X=-30/-60 R_E` 平面；
- `bundle.omitted_satellites` 与 `bundle.space_weather_status` 明确缺失/默认数据；
- `to_html` / `write_html` 支持共享 Plotly JS，避免每个场景内嵌数 MB；`bundle.to_html()` / `bundle.write_html()` 是便捷方法；
- `plots.v3d.prepare_scene/make_scene/dumpGeo` 仍保留为底层入口；`dumpGeo` 的 `include_plotlyjs=False` 会使用共享脚本 URL。

## 运行时结果目录

```python
from common.config import configure
from plots import savefig

configure(output_dir="/data/results/run_20260917", mkdir=False, validate=False)
savefig(figure, "figures/overview.png")
```

相对输出路径解析到本次结果目录；绝对路径保持显式覆盖。配置也可通过 `GENA_OUTPUT_DIR` 指定，未指定时使用数据根下的 `GENA_RES`。报告默认写入该目录的 `reports/` 子目录；`dumpGeo()` 默认按日期建子目录。直接调用第三方对象自身的保存方法不受此配置管理。

## 通用时间序列入口

```python
from plots import tplot

figure, axis = tplot(data, "density")
tplot(data, "velocity_gsm", ax=existing_axis)
```

`DataFrame.tplot()` 与顶层便捷入口 `plots.tplot` 均委托同一实现，实现在 `plots.timeseries`。变量查找支持完整名、带探头前缀和唯一后缀；歧义后缀会报错，不会静默选择。

## 全局样式

`plots` 在 import 时不修改 `matplotlib.rcParams`。需要统一字体和 legend 样式时，在应用入口显式调用：

```python
from plots import theme
theme.apply_theme()
```

渲染 worker、notebook 和命令行脚本应各自在创建 figure 前调用一次。`theme.apply_theme(overrides)` 可覆盖单项，并返回更新后的 rcParams 快照。
