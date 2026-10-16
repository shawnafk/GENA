# 坐标数组与可选转换后端

`common.coord.Coordinates` 保存 NumPy 三分量数组与 UTC 时间轴，以及坐标系、单位、参考原点和可选姿态。它不继承、不保存任何后端的坐标对象。`common._coord_backends` 负责适配 SpacePy、Astropy 和 PySPEDAS，依赖只在调用对应转换时加载。

坐标系名称由调用方直接使用标准大写形式，例如 `ECI2000`、`GSM`、`GSE`、`QQ`。公共层不自动修改大小写或将 `J2000` 改为 `ECI2000`；PySPEDAS 适配器仍将标准 `ECI2000` 翻译为其内部要求的 `J2000`。

## 基本用法

```python
from common.coord import Coordinates

position = Coordinates(
    [[10.0, 1.0, -2.0], [12.0, 2.0, -1.0]],
    ["2024-01-01T00:00:00", "2024-01-01T00:00:10"],
    frame="ECI2000",
    unit="Re",
    origin="earth",
)
gsm = position.to("GSM", backend="spacepy")
gsm_spedas = position.to("GSM", backend="pyspedas")
midpoint = position.resample_at(["2024-01-01T00:00:05"])
```

`values` 始终为 `(N, 3)` 浮点数组，单个 `(3,)` 输入扩展为一行；`times` 为 `(N,)` 的 `datetime64[ns]` 数组，按项目约定表示 UTC。标量时间广播到全部样本，时间序列必须与样本数一致。输入数组在构造时复制；公开数组可以修改，后续转换始终读取当前数据，没有派生缓存。修改数组时需要自行保持时间与样本对齐。

数值时间经 `common.time.Time` 解释为 Unix 秒；其他数值时间格式先显式构造 `Time(value, format=...)`。`to()` 和 `resample_at()` 返回新的 `Coordinates`，原对象不会改变。

## 直接转换数组

所有坐标系之间的数组转换统一使用 `transform_vectors()`；`Coordinates.to()` 调用同一实现。结果保持输入的 `(3,)` 或 `(N, 3)` 形状，类中的数据仍统一保存为 `(N, 3)`。

```python
from common.coord import transform_vectors

qq = transform_vectors([0.0, 0.0, -1.0], source="A", target="QQ")
agsm = transform_vectors([[1.0, 2.0, 3.0]], source="GSM", target="AGSM")
gsm = transform_vectors(position.values, times=position.times,
                        source="ECI2000", target="GSM", backend="pyspedas")
```

A/B/QQ 内部与 GSM/AGSM 之间的固定旋转无需时间；跨 QQ 与 ECI2000 的姿态旋转需要姿态；调用标准后端的转换需要时间。提供时间或姿态时，入口统一检查它们与样本数是否一致。源和目标坐标系之间的小函数为私有实现，不再分别公开。

`common` 顶层仅重新导出 `Coordinates`；数组转换和几何运算都从 `common.coord` 显式导入。

## 后端范围

| 后端 | 支持范围 | 约束 |
|---|---|---|
| `spacepy` | SpacePy 的笛卡尔标准坐标系，如 ECI2000、GSE、GSM、GEO、SM | 使用 SpacePy 自身转换后端，参考原点必须为 earth |
| `pyspedas` | cotrans 的标准坐标系，如 J2000、GSE、GSM、GEO、SM | 输入直接使用数组，不创建 tplot 变量；ECI2000 映射为 J2000 |
| `astropy` | ITRS、GCRS、CIRS、TEME，以及 ICRS、GALACTIC；GEO 映射为 ITRS | 使用原生位置变换，要求长度单位；地心与太阳系质心之间的原点变换不在本接口范围内 |

Astropy 不将 GCRS 或 ICRS 冒充 ECI2000，也不提供项目自定义的 GSE/GSM 扩展。指定不支持的转换会抛出 `NotImplementedError`，不会自动改用其他后端。不存在的 backend 抛出 `ValueError`，缺少所选依赖时给出 `ImportError`。无需标准后端的本地旋转只检查 backend 名称，无需安装对应库。

```python
terrestrial = Coordinates([[6378.1, 0.0, 0.0]], "2024-01-01", "ITRS", unit="km")
geocentric = terrestrial.to("GCRS", backend="astropy")
```

Astropy 使用随安装提供的 IERS 表，关闭转换过程中的自动下载。`Re` 按项目 `common.constants.RE` 换算。不同后端的天文与磁场模型可能不同，不要求跨后端结果逐位一致。

## 项目坐标系与姿态

A/B 与 QQ 之间使用固定安装旋转；QQ 与 ECI2000 之间使用姿态；AGSM 与 GSM 之间使用固定 4.2° 旋转。这些步骤由公共层完成，不在适配器中重复实现。

```python
probe_vectors = Coordinates(
    [[0.0, 0.0, -1.0]],
    "2024-01-01",
    "A",
    attitude=[0.0, 0.0, 0.0, 1.0],
)
look_gsm = probe_vectors.to("GSM", backend="spacepy")
```

姿态四元数为 SciPy 的 xyzw 顺序，定义 QQ → ECI2000 的 `R @ v`，构造后仅保存为 NumPy 数组。四元数格式检查、应用与 SLERP 分别由 `common.coord.quaternion_array`、`rotate_vectors`、`interpolate_quaternions` 完成；`Coordinates` 负责定义姿态的 QQ → ECI2000 含义。接受单四元数、`(1, 4)` 或 `(N, 4)`，也接受 SciPy `Rotation` 并提取其数值。跨越 QQ/A/B 与标准坐标系时必须提供姿态；A/B/QQ 内部转换无需姿态。可以在 `to(..., attitude=...)` 中显式指定。

转换不会给向量添加位置偏移，也不会把月心位置自动改为地心位置。SpacePy/PySPEDAS 路径保持数值单位，Astropy 原生路径仅用于明确长度单位的位置。加减运算要求时间轴、坐标系、单位和参考原点一致；本地坐标系还要求姿态一致。

## 时间重采样与几何函数

`resample_at()` 在原始坐标系内逐分量线性插值，姿态使用 SLERP。源时间排序、重复时间保留第一条，目标顺序与重复项保留；禁止超出源覆盖范围，不进行外推。方向向量线性插值后不保证单位长度，需要调用方按用途归一化。

`common.coord` 提供无时间轴、无坐标系依赖的 `angle_between`、`angle_to_plane`、`vector_to_fov` 和 `fov_to_vector`。FOV 函数由 GENA 投影与视场判断使用；角度函数目前没有生产调用，作为基础几何操作与坐标接口集中在同一文件，但不属于坐标转换入口。角度默认为度，向量夹角与平面角可通过 `degrees=False` 返回弧度。平面角默认有符号，正号由法向量决定。

FOV 函数通过 `direction="propagation"` 使用 +Z，或 `direction="look"` 使用 -Z；两者的横向正方向均为 +X/+Y。零角指向对应的中心轴。反变换限制角度严格处于 `(-90°, 90°)`；正向变换可返回后半球角度，因此仅在前半球保证往返。自定义基向量按行存储，形状为 `(3, 3)` 或逐样本 `(N, 3, 3)`，必须正交归一。

## 旧接口

`coord(times, values, frame, rot, backend="spacepy")`、`convert()`、`obj["GSM"]` 保留。`convert()` 存储旧接口的数组结果，`sync()` 更新原始数据并清除这些结果；不再保留 SpacePy `temp` 对象。未使用的旧 `coord.rot` 属性已移除；旧构造中的 `rot` 参数仍用于传入姿态。旧构造中省略姿态或传入 `0` 仍表示单位姿态，新类则要求显式提供跨本地坐标系边界的姿态。

逐坐标系函数及重复的旧几何函数名已移除。坐标转换统一迁移到 `transform_vectors()`；`fov_2_vec`/`vec_2_fov` 迁移到 `common.coord.fov_to_vector`/`vector_to_fov`。标量 FOV 角度的新返回为 `(3,)`；需要一行数组时显式使用 `np.atleast_2d`。基向量统一使用 `(3, 3)` 或逐样本 `(N, 3, 3)`。旧 `coord` 加减也要求原始坐标系和时间轴一致；需要先显式转换、对齐后再运算。
