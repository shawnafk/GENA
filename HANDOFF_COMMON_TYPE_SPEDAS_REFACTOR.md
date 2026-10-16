# Handoff：拆分 `common.type` 与 pySPEDAS/CDF 适配逻辑

> 这是一份历史交接记录，正文中的 `common.type` 和独立 `common.frame` 描述对应当时的中间状态。当前表类型、表操作、`concat` 与 `read_variable` 已统一到 `common.frame`，现行代码不要再导入 `common.type`。

日期：2026-10-10

## 本轮落实结果（2026-10-10）

以下是本轮完成的实现；后文保留原始审查记录，其中“当前源码”等描述指修改前的状态。

- 后续共识更新：`TimeIndexedDataFrame.__init__`（由 `TimeSeriesDataFrame` 继承）强制要求唯一的 `Datetime` 列，仅调用 `common.time.Time(...).datetime` 完成时间归一化；不会从 `Date` 或 `Timestamp` 推断时间。无效时间由公共时间函数拒绝，`Timestamp` 仍由 loader 按需提供。`to_dict`、文本和 HTML 显示继续使用 pandas 实现；不含 `Datetime` 的选列结果返回普通 DataFrame。
- `common.type` 保留共享表类型、选行、分段和重采样；`concat` 与 `read_variable` 位于 `common.frame`，已移除专用字典转换与数组装箱函数。
- THEMIS/MMS 在统一的 `load/load.py` 中声明产品字段布局，`parse_spedas()` 适配已解码的后端变量。已知向量明确展开；未知多维数据保留为数组；额外字段可通过 `variable_layouts` 指定类别、列名、时间维度和附加轴。
- OMNI HRO 独立读取共享 `Epoch`、标量字段、填充值和单位；indices 从各自来源生成标量说明；GENA Parquet 由 `load_gena()` 生成时间列。它们不需要经过 pySPEDAS 表转换。
- `assemble()` 只执行明确变量说明的外连接与列字典组装，再调用类型构造函数。静态轴与随时间变化的轴明确区分，避免轴长度恰好等于样本数时被误裁剪。
- `return_type="dict"` 直接合并原始变量，保留原生时间、无效样本和后端属性；私有 `frame_to_spedas_dict` 仅供显式导出具有变量元数据的表，不占用 pandas 的 `to_dict`。
- 绘图通过 `common.frame.read_variable` 读取一个目标变量，不再把整表转换为旧式字典，也不依赖私有 loader。
- 已迁移缺失 `common.utils` 的相关导入：时间函数位于 `common.time`，插值和其他工具位于 `libs.utils`；现有类型路径统一使用 `gena.type`。修正了 `common.coord` 中越过顶层包的工具导入。
- README、设计文档和 gena 架构参考已同步。

离线验证覆盖显式向量、三能道数组、高维分布、静态与动态轴、非第零维时间、独立时间外连接、单位/属性、重复时间、Parquet、OMNI CDF 模拟字段、绘图、同步、统计与 A/B 探针投影。真实任务 CDF 尚未逐产品验证，也未触发下载。扩展报告测试仍受已有的 `plots.advance.ani` 缺失阻断。

### 手动构造方式

构造函数只接收已经整理好的表格数据，必须提供 `Datetime`；该列的类型由公共时间函数统一。例如：

```python
import pandas as pd
from common.type import TimeSeriesDataFrame

seconds = [1.0, 2.0]
columns = {
    "Datetime": seconds,  # Unix 秒由 common.time.Time(...).datetime 统一转换
    "Timestamp": seconds,
    "density": [10.0, 20.0],
}
frame = TimeSeriesDataFrame(columns)
frame.to_dict(orient="list")  # 标准 pandas 行为
```

## 目标

`common.type.TimeSeriesDataFrame` 应保持为跨任务通用的时间序列表类型，只负责标准 `Date`/`Timestamp`、选行、分段、同步和表拼接。当前 `dict_to_dataframe()`、`dataframe_to_dict()`、`to_variable_dict()` 与 `_array_cells()` 实际处理的是 pySPEDAS 从 CDF 解码后返回的 `x/y/v/v1/v2/v3` 变量字典，以及少数主动模仿该格式的加载结果，不应继续被描述为通用 DataFrame 转换能力。

本 handoff 记录本轮审查已经达成的判断、当前源码状态、严重风险和建议迁移顺序。此文档只记录后续工作，没有继续修改转换代码。

## 本轮已经达成的判断

- `_clip_variable_dict()` 只服务外部观测加载器，位于统一加载模块 `load.load` 内部，不属于公共接口。
- 多个时间序列表的组合应保留为模块函数，不应改成实例方法或类方法；重复时间应通过稳定排序保证较早输入优先。
- `sync()` 的 linear 和 nearest 路径应共用一个带 `method` 参数的内部函数，不引入 SciPy；linear 使用 `np.interp`，nearest 使用 `np.searchsorted`，外围代码负责缺失值、禁止外推和 `max_gap`。
- 同步前应先清理每列无效样本，并同时过滤对应时间；相关命名使用 `valid`，不使用 `present`。当前实现已有 `_value_is_valid()` 和 `_clean_samples()`。
- `TimeSeriesDataFrame.to_variable_dict()` 作为单对象转换入口在接口形式上是合理的，但其语义只适用于带 pySPEDAS 变量元信息的 frame，因此不能被当作普通 DataFrame 的通用转换方法。
- `_array_cells()` 本身并非只支持二维数组：它始终沿第 0 维拆分，输入为 `(time, a, b, ...)` 时，每个单元格保留形状 `(a, b, ...)`。真正严重的问题是第 0 维被隐式假定为时间，以及变量类别由数组形状和列名猜测。

## `dict_to_dataframe()` 实际支持的格式

当前输入不是任意字典，也不是任意 CDF 内容，而是以下 pySPEDAS 风格变量字典：

```python
{
    "variable_name": {
        "x": time_values,
        "y": observations,
        "v": optional_axis,
        "v1": optional_axis_1,
        "v2": optional_axis_2,
        "v3": optional_axis_3,
        "units": optional_units,
        "attrs": optional_backend_metadata,
    }
}
```

THEMIS、MMS 和 OMNI 加载器直接处理这种后端结果，indices 加载器则先把自身结果包装成相同形状。这个事实说明转换代码属于加载边界的 pySPEDAS 适配层，而不是 `common.type` 的基础能力。

## 当前实现中的严重问题

### 1. 用形状和变量名猜测变量类别

当前代码把任意 `(time, 3)` 数组识别成三分量向量，除非变量名包含 `"en_eflux"`：

```python
elif values.ndim == 2 and values.shape[1] == 3 and "en_eflux" not in name:
```

这会导致三能道能谱、三类别标志、三元素矩阵切片等数据被错误展开为 `<name>_x/<name>_y/<name>_z`。反过来，只有名称命中 `"en_eflux"` 的三能道数据才能逃过错误分类。变量意义不应由字符串片段决定。

### 2. 隐式假定第 0 维是时间

`len(dates) == len(values)` 只验证 `values.shape[0]` 与时间数量一致，随后 `_array_cells()` 固定沿第 0 维拆分。对于 pySPEDAS 常见输出，这通常成立；对于其他 CDF 解码器、转置数组或显式指定了其他时间轴的数据，这个假设会产生错误结果或含糊报错。

### 3. 坐标轴名称被写死

代码只认识 `v`、`v1`、`v2`、`v3`，并把 `v` 强制命名为 `<name>_energy`。`v` 不一定表示能量，其他数据源也可能使用完全不同的轴名称。这个命名规则是 pySPEDAS/THEMIS 使用习惯，不是通用时间序列表规则。

### 4. 输出列名被写死

三分量向量固定展开为 `_x/_y/_z`，`v` 固定变成 `_energy`，其余轴固定变成 `_<axis_key>`。这要求下游绘图和逆转换理解同一套命名规则，使 `common.type`、加载器和绘图层产生隐式耦合。

### 5. 逆转换依赖专用 `attrs`

`dataframe_to_dict()` 只有在 `frame.attrs["variables"]` 完整存在时才知道哪些列组成一个变量、变量属于 scalar/vector/array，以及哪些列是附加轴。普通 DataFrame 即使列和值完全合法，也不能可靠转换。经过不保留或修改 `attrs` 的 Pandas 操作后，逆转换也可能静默返回空字典或不完整结果。

### 6. 往返转换不是完全保真的

`dict_to_dataframe()` 把单位写入 `attrs["variables"][name]["units"]`，但当前 `dataframe_to_dict()` 没有把单位写回变量字典。其他后端属性也没有完整保存。数组逆转换使用 `np.stack()`，因此不规则形状的单元格无法还原。

### 7. `_array_cells()` 的名称隐藏了关键假设

函数可以处理二维、三维和更高维数组，也可以在直接收到序列时容纳不同形状的单元格；它的真实语义是“把第 0 维的每个时间样本复制成一个 object 单元格”。建议后续改名为 `_time_sample_cells()`，并由调用者显式传入或验证时间轴，而不是继续使用含糊的 `_array_cells()`。

### 8. 不应把 pySPEDAS 转换塞进 `TimeSeriesDataFrame.__init__`

当前源码有注释提出把 `dict_to_dataframe()` 放入构造函数。不要这样做。构造函数需要保持与 Pandas DataFrame 的正常构造方式兼容；让它同时猜测普通 mapping、pySPEDAS 变量字典、CDF 路径和 Parquet 路径会让输入语义更不明确，并增加子类构造、Pandas 内部操作和 `_constructor` 调用时的风险。CDF/pySPEDAS 与 Parquet 应由各自加载器解析，解析完成后再构造表。

## 推荐的模块边界

pySPEDAS 适配代码与命名加载器集中在 `load/load.py`：

```text
common/type.py
    TimeSeriesDataFrame
    concat(...) 或命名一致的通用拼接函数
    通用时间同步内部函数

load/load.py
    parse_spedas(...)
    assemble(...)
    frame_to_spedas_dict(...)
    load_gena/load_themis/load_mms/load_omni/load_indices
```

THEMIS 和 MMS 的加载器先调用 `parse_spedas()`，再把返回的明确变量说明交给 `assemble()`。OMNI CDF parser 和 indices loader 直接生成同一种变量说明。

Parquet 数据不经过 pySPEDAS parser。`load_gena()` 继续读取表结构并构造对应的 GENA frame，CDF 与 Parquet 的入口不需要统一到 `TimeSeriesDataFrame.__init__`。

## 变量类别必须显式确定

短期最低限度可以先用附加轴判断数组类型，避免三能道能谱被识别成向量：

```python
axis_keys = ("v", "v1", "v2", "v3")
has_axis = any(key in variable for key in axis_keys)

if values.ndim == 1:
    kind = "scalar"
elif values.ndim == 2 and values.shape[1] == 3 and not has_axis:
    kind = "vector"
else:
    kind = "array"
```

这只是过渡方案，因为没有附加轴的 `(time, 3)` 数组仍然存在歧义。最终方案应由加载器或后端元信息显式给出 `kind`、时间轴、分量名称、附加轴名称和单位；通用层只验证并执行，不再根据形状或变量名推断含义。

## 绘图层需要一起解耦

`plots.basic.timeseries` 当前直接导入 `dataframe_to_dict()`，`plots.basic.themis` 也会调用 frame 的旧式字典转换方法。绘图不应为了读取一个变量而把整个表还原成 pySPEDAS 字典。

面向 frame 元信息的只读变量访问由 `common.frame.read_variable()` 提供；通用绘图使用这个函数读取目标变量，不依赖 `load` 的内部 parser。

## 当前工作区的接口不一致

开始后续重构前必须先处理以下当前状态，否则测试和加载器不能可靠运行：

- `common/type.py` 当前定义 `concat()`，但 `common.__init__`、加载器和测试仍引用 `concat`。
- `TimeSeriesDataFrame` 当前定义 `to_dict()`，但加载器和测试仍调用 `to_variable_dict()`；直接覆盖 Pandas 的 `DataFrame.to_dict()` 还会改变用户熟悉的公共行为，因此不应保留这个无参数覆盖版本。
- 仓库当前缺少 `common/utils.py`，而 `common.type`、加载器、绘图和测试仍从 `common.utils` 导入；pytest 因 `ModuleNotFoundError` 在收集阶段停止。
- `common/type.py` 中存在审查期间加入的临时注释和未完成命名，应在正式重构时统一清理。
- `.agents/skills/gena/references/architecture.md` 和设计文档仍把 `dict_to_dataframe()`、`dataframe_to_dict()` 描述为 `common.type` 的共享能力，迁移完成后必须同步修订。

## 建议迁移顺序

1. 先恢复一个可导入、命名一致的基线：解决 `common.utils` 缺失，统一 `concat`/`concat`，恢复 `to_variable_dict()` 或一次性迁移全部调用方，禁止覆盖 Pandas `to_dict()` 的既有语义。
2. 为当前 pySPEDAS 字典行为建立独立测试夹具，覆盖 scalar、三分量 vector、三能道 spectrum、二维能谱、三维 ESD、静态轴、随时间变化的轴、单位和独立时间轴。
3. 在 `load/load.py` 中集中实现 parser、`assemble()` 和命名加载器，并用显式变量说明替代按数组形状猜测类别。
4. 更新 THEMIS、MMS、OMNI 和 indices 加载器的导入及 `return_type="dict"` 路径，然后从 `common.__all__` 和 `common.__init__` 移除 pySPEDAS 转换函数。
5. 用显式变量说明替代 `shape[1] == 3` 和 `"en_eflux"` 判断，并明确第 0 维时间假设；对不满足要求的输入给出包含变量名和实际 shape 的错误。
6. 将绘图变量访问从“整表还原为旧式字典”改为直接读取 frame 列和 `attrs["variables"]`，保持绘图层不依赖加载层。
7. 决定旧式字典转换是否需要面向用户公开；如果只是 `return_type="dict"` 的内部实现，保持私有。如果确有外部用户，公开入口应明确命名为 pySPEDAS 转换，不能继续使用泛化名称。
8. 更新 README、架构参考和设计记录，并运行离线回归测试。

## 必须增加或保留的测试

- `(time,)` 标量变量转换与往返。
- `(time, 3)` 三分量向量由显式说明识别并展开。
- `(time, 3)` 三能道能谱即使变量名不含 `en_eflux` 也不能被识别成向量。
- `(time, energy)` 能谱及静态 `v` 轴。
- `(time, energy, angle)` ESD 数据及 `v1/v2` 轴。
- `(time, a, b, c)` 更高维数组按第 0 维形成单元格。
- 不同变量使用不同原生时间轴时的 outer join。
- 时间轴不在第 0 维时明确拒绝，除非调用方显式指定轴。
- 不规则形状样本明确拒绝或明确支持，不能在 `np.asanyarray()`/`np.stack()` 中得到难以理解的异常。
- 变量名完全不包含任务或仪器关键词时仍按显式说明正确转换。
- 单位和需要保留的后端属性能够往返。
- 重复时间的稳定优先级、空输入和缺少必需键的错误行为。
- linear/nearest 同步、同步前无效样本清理、数组单元格复制和 `max_gap` 行为保持不变。

## 完成标准

- `common/type.py` 中不再出现 `x/y/v/v1/v2/v3`、`en_eflux` 或 pySPEDAS 专用列名推断。
- `TimeSeriesDataFrame` 不覆盖 Pandas `to_dict()` 的标准签名与行为。
- 所有加载器仍保持当前面向用户的 DataFrame 默认返回和显式旧式字典返回能力。
- 三能道能谱与三分量向量有确定且可测试的区分方式。
- CDF/pySPEDAS 与 Parquet 分别由加载层处理，`TimeSeriesDataFrame.__init__` 不解析文件或猜测专用字典。
- `common`、`load`、`plots` 和测试中的函数名称、导入与文档完全一致。
- 目标测试、语法检查和导入检查通过；不触发网络下载。
