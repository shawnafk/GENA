# design 文档核查与项目修改报告

日期：2026-09-17。核查对象：本地工作目录 `/Users/shawn/syncthing/GENA-CODE-Dev` 中的 `design.md` 和实际源码。第 11 节是章节编号，其中实际列出 10 项问题。本报告逐项记录核查结果、实施内容、验证和剩余限制。

## 1. 修改结果

已完成加载入口、类型定义、时间工具、基础绘图和高级绘图入口的整理，统一项目输出路径配置，修复报告的缺失导入与 `libs.flux` 的语法错误，并同步修订 design、README、相关说明和 GENA 技能的架构指引。修改以移动现有实现和调整调用为主。

本次主要回归集合为 **153 项通过**。另外检查 `tests/test_flux.py` 时，现有模型的 7 项测试通过，4 项测试因源码中已经缺失的旧 API 失败。全量测试收集仍有 11 个错误，原因是已经缺失的研究脚本或模型模块，详见第 6 节。因此，本报告不宣称整个历史测试集已通过。

## 2. design 是否符合实际

原文关于混合分层、GENA 专用处理、默认表格式、THEMIS ESA 产品时间适配、区域判据列映射、报告具有 I/O 等描述基本符合源码。主要问题在于现状、拟议结构与历史功能混在一起。

| 原文位置或说法 | 修改前源码证据 | 核查判断及修正 |
| --- | --- | --- |
| 第 3 节把 GENA load 留在 gena，第 8 节拟议迁移 | `gena/load.py` 同时包含类、读取和姿态辅助函数 | 事实正确；本次落实拆分并更新现状 |
| `common` 中既有表类型又有基础绘图 | `common/timeseries.py` 与 `common/vis.py` 存在 | 事实正确；类型和绘图分别迁入新归属 |
| 第 6 节标题写“待实现功能” | `libs.regions` 已提供 `columns` 映射；现有跨列名测试可运行 | 标题不准确；改成“已有能力及适用边界” |
| ani/v3d “缺失” | 文件实际在 `plots/advance/ani.py` 和 `v3d.py`，export 仍从 gena 导入 | 缺的是正确引用；没有重新编写或恢复一套动画实现 |
| `plots.advance` 提供 Panel 等入口 | 实现已在 `joint.py`，但 advance 目录没有 `__init__.py` 导出入口 | 原说明遗漏迁移后的入口断裂；本次补齐 |
| `libs.flux` 第 679 行后为空函数体 | `integral_ena_counts(...)` 后直接开始顶层 `__all__` | 属实，阻止整个模块导入 |
| libs 中的能量输运、BBF 正演和响应推断 | 当前 `.py` 文件不在 libs，部分文件在 `undefined.zip` | 不能列成现有 API；已明确历史/存档状态 |
| `libs.environment` 与 OMNI load 重复 | `sw_keys()` 组合 OMNI、Dst，插值并提供失败回退 | 职责不同，保留该组合服务 |
| 本表 `pos()/rot()` 也进行清理和插值 | 源码确实调用位置插值和 SLERP | 原 design 此处正确，保留；修正技能参考中“本表旋转不插值”的旧说法 |
| common 完全不依赖绘图层 | 共享表 `tplot()` 方法在调用时导入绘图 | 应说明这一便利方法的延迟依赖；不能声称严格无反向调用 |
| 配置只是在拼路径 | 实际还涉及环境变量、校验、建目录与缓存 | 这些边界有必要；简化实现，保留职责分离 |

额外发现：`common.spectrum` 已缺失，但 GENA 统计和 THEMIS 能谱仍导入它，导致报告与产品绘图无法导入。项目现有 `undefined.zip` 包含原始 `spectrum.py`，本次将它复用到 `libs/spectrum.py`。归一化和误差传播算法沿用存档，只将 NumPy 2 的异常入口改为 `np.exceptions.AxisError`。

## 3. 第 11 节逐项修改

### 3.1 问题 1：GENA 加载入口统一

2026-10-15 更新：所有 parser、`assemble()` 和具名任务加载器集中在 `load/load.py`，GENA 入口为 `load_gena()`。GENA 的表类型和位置姿态方法位于 `gena.frame`。

依赖顺序为 `load → gena.frame → common.frame`。类型定义不反向导入加载器；`gena.post.project` 直接引用实际归属模块中的类型和处理函数。

加载器的文件发现、时间裁剪、排序、事件去重、Parquet 列投影、损坏文件恢复和姿态插值算法均保留。已同步现行代码、测试和使用示例。旧 `gena/load.py` 被移除，没有增加转接文件。

### 3.2 问题 2：过滤接口统一且简洁

保留已广泛使用的 `FilterNode`，只承担延迟掩码计算和 `&`、`|`、`~` 组合。移除没有现行调用方的 `SelectorNode` 和 `gena.filter.select()`，选行统一复用共享表的 `select()`，因此子类和 `attrs` 也由同一实现维护。

`filter_rule()` 改为简短的两种调用分派，增加 `functools.wraps`，保留函数名称、说明和可检查签名；同时支持 `energy_eff(df=frame, ...)` 的关键字形式。删除文件末尾注释字符串中的重复过滤实现。文件由 549 行减少到 474 行。

```python
from gena import filter as gf

# 直接获得掩码。
mask = gf.energy_eff(l2, Emin=4, Emax=100)

# 组合条件，选行复用数据表方法。
rule = gf.energy_eff(Emin=4, Emax=100) & gf.onboard_eff()
selected = l2.select(rule)
```

通道长度、起止关系、求和/差分、TOF、空间与能量阈值、golden/soft 规则及 A/B 标定参数均保留。没有引入新的过滤框架。

### 3.3 问题 3：配置为何需要这些层，以及如何简化

`common/config.py` 从 160 行整理为 80 行。保留三个已有职责：`resolve_paths()` 只解析；`configure()` 显式校验和建目录；`get_configured_paths()` 供加载器复用已校验结果。新增的 `resolve_output_path()` 只负责解析结果路径。

原 `configure(file_base=...)` 只返回字典，后续无参数加载并不会自动使用这次传入的目录；旧缓存键也不包含实际环境变量解析结果。这会让“运行开始配置一次”的预期落空。现在显式配置会保存到当前进程，缓存按实际解析后的路径和校验选项建立，环境变化不会继续命中旧目录。返回缓存的字典副本，防止调用者误改缓存内容。

导入 config 不读取文件、不建目录。单纯保存已加载的数据不必验证数据盘。环境变量原有的 `GENA_RAW`、`GENA_PQT`、`SPEDAS`、`GENA_RES` 仍作为数据根下的子路径使用。

### 3.4 问题 4：类型单独管理

共享类型及配套表操作、同步和变量读取实现统一放在 `common/frame.py`，GENA 专用类型放在 `gena/type.py`。没有把 GENA 仪器姿态处理放进 common，也没有为了拆文件新增基类或多层代理。

2026-10-15 更新：产品和存储格式转换统一位于 `load/load.py`。`parse_spedas()` 和来源 parser 生成明确变量说明，`assemble()` 拼成时间序列表，具名 `load_*` 函数组织完整加载过程。`common.frame` 保存表类型、实例操作、`concat()` 和 `read_variable()`；构造时检查 `Datetime` 并调用 `common.time.Time(...).datetime`，`to_dict()` 保持 pandas 行为，`tplot()` 按需导入 `plots.basic`。

旧 `common/timeseries.py` 和 `common/type.py` 已被移除；现行使用方统一引用 `common.frame`。

### 3.5 问题 5：时间操作职责集中

时间格式转换统一通过 `common.time.Time`：构造表时调用 `Time(values).datetime` 得到不携带时区信息的 UTC 纳秒时间。旧独立函数 `normalize_datetime()` 和 `normalize_time()` 已删除；裸数字表示 Unix 秒，毫秒、儒略日等输入需显式指定 `format`。

`split_time_intervals()` 位于 `common.frame`，直接处理 pandas datetime Series 或 DatetimeIndex，不执行格式转换。表重采样保留在 `TimeSeriesDataFrame` 中，因为它依赖列、类型、缺失值和元数据；粒子事件类型 `ENA_DataFrame` 不继承插值重采样能力。

另外移除 `utils.py` 中被后一个同名函数覆盖的旧 `is_multi_layer()` 定义；实际运行的后一个实现保持不变。`Ticktock` 的导入移到使用它的 `tick()` 内，使基础表类型不因时间工具导入而立即加载 SpacePy。

### 3.6 问题 6：基础绘图移到 plots

`common/vis.py` 整体移到 `plots/vis.py`。坐标轴装饰器、通用 `tplot`、axes 分配、colorbar、刻度和卫星位置标签等继续复用原实现。产品绘图、联合图、GENA 轨道图与报告同步引用新位置。

`plots` 顶层提供 `savefig`、`tplot` 和 `vis`；common 顶层移除这两个绘图函数的重导出。旧 `common/vis.py` 不保留转接。原来该模块对 `libs.ephemeris` 的依赖现在位于绘图层，common 基础绘图反向依赖科学服务的问题随归属调整消除。

此次没有扩展到全局样式改造。`plots.vis` 和部分产品模块仍会设置少量 `rcParams`，design 中明确保留这一现状。

### 3.7 问题 7：接通报告、动画、三维与联合图

`gena.post.export` 的动画、三维导入分别改为 `plots.advance.ani` 和 `plots.advance.v3d`。增加 `plots/advance/__init__.py`，从实际实现 `joint.py` 导出 `Panel`、`overview`、`env`、`plot_region_code`，使原有公共联合绘图入口可用。

复用用户已经移动的两个文件，没有复制另一份动画或三维代码。`dumpGeo()` 与其他 GENA 数据入口统一读取 `configure()` 生效后的 `pqt_base`，不再接受独立的数据目录覆盖参数。

实际验证包括：A/B 探头报告的 JPEG 文件、观测汇总 Excel、带动画控件的三维 HTML、报告的可选动画/几何调度；另用 30 行合成事件调用 `dumpFov()`，成功生成一个 33,786 字节的 MP4，文件头包含 `ftyp`。这次 MP4 使用了真实 ffmpeg 编码器，不仅检查函数是否被调用。

### 3.8 问题 8：flux 语法错误及未完成函数

为原空函数 `integral_ena_counts()` 增加说明并明确抛出 `NotImplementedError`，恢复整个模块的解析与导入。它的原参数只有能量、通量和 density，没有路径长度、面积、视场、效率和积累时间，不能据此擅自定义探测器计数算法。已有计算应使用 `qk_estimate_ena_counts()`，并显式指定科学参数。

把文件末尾的积分演示放入 `if __name__ == "__main__"`，避免导入时执行演示积分和打印。Maxwellian/Kappa 分布、共享漂移和约束拟合实现未改动。

本次没有把该草稿包装成看似可用的计数函数，也没有恢复测试中引用的 `dual_moment_partition()`、`fit_dual_fixed_partition()`。修改前备份中同样不存在这两个定义；4 项失败详见第 6 节。

### 3.9 问题 9：文档入口和实际能力同步

ESD 和 L1 角分布文档现在使用统一入口 `load_themis()`。历史 observation review 中的绘图入口更新到 `plots.basic.gena` 和 `plots.advance.env`，同时增加历史记录提示，原测试次数明确不作为本次验证。

同步中文/英文 README、绘图说明和 GENA 后处理说明中的现行导入、目录与过滤用法。后处理说明保留历史科学背景，但标出当前没有的文件和模型。design 重写了当前目录、依赖关系、输出配置及问题处理状态。

依据 GENA 技能的实现核查要求，同步更新 `.agents/skills/gena/SKILL.md` 与 `references/architecture.md`，以便后续工具查询使用新路径。修正了技能参考里过时的本表姿态、区域目标分类和缺失模型说明。

### 3.10 问题 10：统一运行时结果目录

新增 `configure(output_dir=...)` 和环境变量 `GENA_OUTPUT_DIR`。所有相对输出名经 `resolve_output_path()` 解析到运行时结果目录；绝对路径是显式覆盖。未指定独立输出目录时，继续使用 `FILE_BASE` 下的 `GENA_RES`。

| 输出入口 | 新行为 |
| --- | --- |
| `plots.savefig` / `plots.vis.savefig` | 相对图片名位于结果根目录下 |
| `export.report` | 默认位于 `<结果根>/reports/<日期或精确时间段>/` |
| `export.observation_summary` | 默认 `<结果根>/reports/observation_summary.xlsx` |
| `export.dump(..., output_path=...)` | 相对 CSV 路径位于结果根下，自动创建父目录 |
| `ani.dumpFov` / `AnimationDemo3.save` | 相对动画路径位于结果根下 |
| `v3d.dumpGeo` | 默认 `<结果根>/YYYYMMDD/3d.html`，可显式指定子目录 |

以上入口复用同一个路径解析函数，不各自读取环境变量。输入数据目录、原始遥测转换输出和星历缓存没有混入图件结果目录。直接使用第三方 `figure.savefig()`、`DataFrame.to_csv()` 或 `figure.write_html()` 的用户代码仍由调用者传路径，项目没有修改第三方对象的行为。

## 4. 新调用方式

```python
from common.config import configure
from common.frame import TimeSeriesDataFrame
from gena.data import ENA_DataFrame
from load import load_gena
from gena import filter as gf
from gena.post import export
from plots import savefig
from plots.basic import gena as gena_plot

# 使用实际存在的数据根目录；显式配置会被后续无参数加载复用。
configure(file_base="/data/ENA", output_dir="/data/results/run_20260917")

l2 = load_gena(trange=["2026-04-07", "2026-04-08"], product="l2_a")
selected = l2.select(gf.energy_eff(Emin=4, Emax=100) & gf.onboard_eff())
figure, axis = gena_plot.energy(selected)
savefig(figure, "figures/energy.png")

# 写入 /data/results/run_20260917/reports/20260407/。
export.report("A", "2026-04-07")
```

仅处理已经加载的数据时，可用 `configure(output_dir="/data/results/run_20260917", mkdir=False, validate=False)`；保存函数会按需创建实际输出的父目录。

| 旧入口 | 现行入口 |
| --- | --- |
| `gena.load.load` | `load.load_gena` |
| `gena.load.ENA_DataFrame` | `gena.frame.ParticleDataFrame` 或 `gena.frame.ParameterDataFrame` |
| `common.timeseries.TimeSeriesDataFrame` 及表适配器 | `common.frame` |
| `common.timeseries.normalize_datetime/split_time_intervals` | `common.time` |
| `common.vis` | `plots.vis` |
| `common.savefig/tplot` | `plots.savefig/tplot` |
| `gena.ani` / `gena.v3d` | `plots.advance.ani` / `plots.advance.v3d` |
| `gena.filter.select(rule)(frame)` | `frame.select(rule)` |

旧路径不保留兼容层。未位于当前工作目录的个人脚本或 notebook 需要按表更新引用。

## 5. 修改范围与复用证据

当前加载实现集中在 `load/load.py`，GENA 表类型位于 `gena/frame.py`；共享表操作位于 `common/frame.py`，共享绘图位于 `plots/vis.py`。旧 `gena/load.py`、按任务拆分的加载模块、`common/timeseries.py`、`common/type.py` 和 `common/vis.py` 均已移除。

实际逻辑调整集中在 config、filter、输出路径、缺失入口和 flux 导入修复。其他四个外部加载器、`libs.regions`、各产品图、联合图和多数既有测试仅更新模块路径或类型导入。没有修改遥测解码、级别转换、仪器参数、坐标变换或区域判据公式。

通过修改前备份与当前源码的 AST 比较，以下 9 个代码块完全一致（忽略源码行号和空行）：

| 迁移 | 保持一致的代码块 |
| --- | --- |
| `gena/load.py → load/load.py` | GENA Parquet 发现、读取和产品类型选择 |
| `gena/load.py → gena/frame.py` | `ParameterDataFrame`、`ParticleDataFrame` 及 PDF 工程量插值 |
| `common/timeseries.py → common/utils.py` | `normalize_datetime`、`split_time_intervals` |

新增边界测试覆盖运行时目录持久使用、环境变量变更后的缓存、缓存字典隔离、无数据盘条件下保存图片、报告默认目录、A/B 探头规则、元数据保留、关键字过滤调用、函数签名、UTC 和时间数字含义。已有科学测试复用原测试数据与断言。

本目录没有 `.git`，因此未使用 Git diff 或提交。修改前已创建完整范围备份：`/tmp/gena-design-review.smId8B/before.tar.gz`。其中包括本次涉及的源码、测试、文档和技能文件；`undefined.zip` 原存档未修改。被移除的旧模块和选行层可从该备份恢复。该路径属于系统临时目录，需要长期保留时请另存。

## 6. 本次验证记录

环境：Python 3.14.6、NumPy 2.4.6、pandas 3.0.3，绘图使用 Agg 后端。验证使用合成事件、小型本地 Parquet、临时输出及模拟加载边界，没有下载真实观测数据来充当科学验证。

### 6.1 主要回归：153 项通过

```bash
MPLBACKEND=Agg python -m pytest -q \
  tests/test_design_migration.py tests/test_gena_load.py \
  tests/test_timeseries_io.py tests/test_common_utils.py \
  tests/test_common_vis.py tests/test_plots.py \
  tests/test_gena_export_report.py tests/test_gena_v3d_output.py \
  tests/test_gena_counts_nd.py tests/test_gena_statistic_windows.py \
  tests/test_gena_temperature.py tests/test_gena_project.py \
  tests/test_gena_utils_observation.py tests/test_themis_vis.py \
  tests/test_themis_region_classifier.py tests/test_gena_raw_constants.py
```

结果：`153 passed, 18 warnings`。这些测试覆盖类型与加载、时间处理、统计、温度、投影、区域列映射、产品图、实际报告文件和三维输出。

最后一次路径归一化和动画保存入口调整后，重新执行 `test_design_migration.py`、`test_gena_export_report.py`、`test_gena_v3d_output.py`，结果为 `29 passed, 18 warnings`。此数是相关集合的复验，不与 153 重复累加。

### 6.2 flux：7 项通过，4 项失败

`tests/test_flux.py` 中 Maxwellian/Kappa 叠加、漂移、总矩约束、固定 R、纯 Maxwellian、纯 Kappa 及非法 R 的现有测试通过。以下 4 项失败仍存在：

- `test_manual_partition_preserves_total_moments`：缺 `dual_moment_partition`。
- `test_fixed_partition_jointly_recovers_kappa_and_drift_from_energy_flux`：缺 `dual_moment_partition`。
- `test_fixed_partition_handles_a_pure_maxwellian_without_a_kappa_fit`：缺 `fit_dual_fixed_partition`。
- `test_fixed_partition_handles_a_pure_kappa_population`：缺 `fit_dual_fixed_partition`。

这两个函数在修改前的 `libs/flux.py` 中也没有定义；本次语法修复让测试能够进入这些断言，并非迁移删除了它们。没有以新增兼容函数或改变测试断言的方式掩盖失败。

### 6.3 全量收集：11 个既有缺失项

执行 `MPLBACKEND=Agg python -m pytest --collect-only -q`，成功收集 179 项测试，同时出现 11 个收集错误：

| 受影响测试 | 缺失对象 |
| --- | --- |
| `test_bbf_ena.py` | `libs.bbf_ena` |
| `test_ena_response.py` | `libs.ena_response` |
| `test_energy_transport.py` | `libs.energy_transport` |
| `test_ena_drift.py` | `libs.distribution` |
| `test_ena_drift_rest_model.py` | `gena.post.ena_drift` |
| `test_cg_fov_fit.py` | `1.Paper_BBF/ena_cg_fov_peak/fit_cg_alpha.py` |
| `test_conditional_cg_duration_sweep.py` | 同研究目录下 `run_conditional_cg_duration_sweep.py` |
| `test_conditional_cg_short_windows.py` | 同研究目录下 `fit_conditional_cg_short_windows.py` |
| `test_local_fov_peak.py` | 同研究目录下 `analyze_local_fov_peak.py` |
| `test_dec04_vparallel.py` | `1.Paper_BBF/ena_Vbulk/fit_dec04_secondary_vparallel.py` |
| `test_ena_pitch_angle_estimation.py` | `1.Paper_BBF/ENA_pitch_angle_estimation.py` |

未删除这些测试或把它们自动改成跳过。只有报告和产品绘图直接需要的 spectrum 从项目存档复用，其余已经移走的研究模块没有批量恢复。

### 6.4 语法、入口、输出与警告

- 对当前五个主包中 51 个非 deprecated Python 文件进行编译检查，无错误。
- AST API 查询现在能列出 `libs.flux`，没有此前空函数体的语法警告。
- 检查 14 个核心模块导入，全部成功；`plots.config` 在当前环境也可导入。
- 在独立 Python 进程中导入两个表类型，确认没有加载 `matplotlib.pyplot` 或统一加载模块。
- 在现行源码和 tests 中搜索旧 `gena.load`、`common.timeseries`、`common.vis`、旧动画/三维导入及 `common.spectrum`，没有活动引用残留；deprecated 历史代码不在本次改造范围。
- 已执行真实 ffmpeg MP4 导出；JPEG、Excel、HTML 文件由现有测试读取验证。
- SpacePy 对 2026 年时间发出 IGRF 有效期警告，并使用 2025 年有效日期；这限制了科学解释，软件测试通过不表示该日期的磁场计算已获得正确系数。
- 三维依赖 geopack 在导入时尝试列出 IGRF 系数，打印获取列表失败提示；本次没有更改第三方库的导入行为。
- JPEG 报告有 Matplotlib `tight_layout` 的 axes 兼容性警告；文件格式与 300 dpi 验证通过，但未逐页人工审核全部图件布局。

## 7. 仍需明确的范围

本次完成的是第 11 节提出的架构整理和直接阻断修复。`integral_ena_counts` 仍是明确不可用的草稿；缺失的历史科学 API 和研究目录未恢复。真实遥测重解码、A/B 绝对标定和跨仪器判据的物理适用性不属于这批合成测试的证明范围。

`plots.config` 的 `myplot` 来自当前机器的外部路径，未列入 `pyproject.toml`；当前可导入不代表全新安装环境自动具备该依赖。绘图全局样式、`gena.utils` 的混合职责、THEMIS 产品绘图中的部分统计便利计算仍按 design 中的“当前例外”保留。

后续调用请优先使用本报告第 4 节的新入口和统一保存函数。科学计算参数、探头、单位、时间范围和输出目录都保留在调用处，便于复查分析过程。
