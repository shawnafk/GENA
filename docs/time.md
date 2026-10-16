# 通用时间类 `Time`

`common.time.Time` 为标量和一维时间序列提供统一的输入及转换接口。内部保存无时区标记、含义为 UTC 的 pandas `Timestamp` 或 `DatetimeIndex`，精度固定为纳秒。DataFrame 的时间列仍是原生 `datetime64[ns]`，可以继续使用 `.dt`、时间筛选、排序和重采样。

## 构造和访问原生时间

```python
from common import Time

single = Time("2026-10-12T08:00:00.123456789+08:00")
single.datetime       # Timestamp('2026-10-12 00:00:00.123456789')
single.isscalar       # True

sequence = Time([0, 1, 2])
sequence.datetime     # DatetimeIndex，1970-01-01 UTC 的第 0、1、2 秒
sequence.datetime[0]  # pandas Timestamp
sequence.datetime[:1] # 长度为 1 的 DatetimeIndex

df["Datetime"] = sequence.datetime
```

构造接口为 `Time(value, *, format=None, scale=None)`。输入可以是 Python `datetime`、`date`、字符串、NumPy `datetime64`、pandas `Timestamp`、`DatetimeIndex`、Series、数字、已有 `Time`、Astropy `Time` 或 SpacePy `Ticktock`。一维序列可以混合不同类型的标量。Series 的索引和名称不保存到时间对象中；`DatetimeIndex` 表示时间值本身。

无时区输入解释为 UTC。带时区输入转换到 UTC 后去掉时区标记。数组保持数组形状，包括单元素和空数组；`Time` 不实现数组容器接口；索引、切片、迭代、长度和形状直接使用 `.datetime`。多维数组不接受，需要调用方先明确选择或展开时间轴。

## 数字格式和时间尺度

裸数字只表示 Unix 秒，不根据数值大小猜测。毫秒、纳秒和科学时间格式通过 `format` 指定：

```python
Time(2451545).datetime                          # 1970-01-29 08:59:05
Time(2451545, format="jd").datetime             # 2000-01-01 12:00:00 UTC
Time(51544.5, format="mjd").datetime            # 同一时刻
Time(946728000000, format="unix_ms").datetime   # 同一时刻
Time(946728000123456789, format="unix_ns")       # 保留整数纳秒
Time(2451545, format="jd", scale="tt")           # 输入为 JD(TT)，保存为 UTC
```

普通日期时间、pandas、NumPy 和 Unix 数字直接通过 pandas 处理。JD、MJD、GPS 等 Astropy 支持的科学格式、非 UTC 时间尺度和外部科学时间对象通过 Astropy 转换到 UTC，再保存成 pandas 时间。已有 Astropy `Time` 的时间尺度由对象自身定义；不会把 TAI 钟面时间直接标记为 UTC。原始时间尺度不保存在新对象中。

## 输出方法

| 方法 | 标量输出 | 数组输出 |
|---|---|---|
| `.datetime` | `pd.Timestamp` | `pd.DatetimeIndex` |
| `to_timestamp()` | 浮点 Unix 秒 | 浮点数组 |
| `to_iso()` | UTC ISO 字符串，9 位小数和 `Z` 后缀 | 字符串数组 |
| `to_str(fmt="%Y-%m-%d %H:%M:%S")` | 格式化字符串 | 字符串数组 |
| `to_numpy(unit="ns")` | `np.datetime64` | datetime64 数组 |
| `to_astropy()` | UTC Astropy Time 标量 | UTC Astropy Time 数组 |
| `to_tick()` | 一个时间元素的 SpacePy Ticktock | SpacePy Ticktock |
| `convert(format, *, scale="utc")` | Astropy 格式对应的标量值 | 相应数组 |

```python
t = Time(2451545, format="jd")
t.to_timestamp()                 # 946728000.0
t.to_iso()                       # '2000-01-01T12:00:00.000000000Z'
t.to_str("%Y/%m/%d %H:%M")        # '2000/01/01 12:00'
t.to_numpy("s")                  # 降为秒精度
t.convert("mjd")                      # 51544.5
t.convert("jd", scale="tt")           # JD(TT)
t.to_astropy().tai                # 按需使用 Astropy 的其他操作
t.to_tick()                      # SpacePy 坐标转换所需的 ticks
```

Unix 秒从 pandas 保存的纳秒整数计算，标量与数组使用相同规则。浮点秒受浮点数表示精度限制。`to_str()` 使用 pandas 的 `strftime`，其中 `%f` 表示微秒；需要完整纳秒文本时使用 `to_iso()`。普通日期运算直接通过 `.datetime` 使用 pandas 的 `Timedelta`、日期偏移和时间方法。

## 范围、观测日期与连续区间

`Time` 将未定义的属性和方法转交给内部 pandas 时间对象，可以直接调用 `t.floor("D")`、`t.ceil("D")`、`t.strftime(...)` 或读取 `t.year`。返回 pandas 原生结果，不改变实例；范围展开由调用方使用 `pd.date_range`。

DataFrame 分段由 `common.frame.TimeIndexedDataFrame.split()` 负责，调用同模块的 `split_time_intervals()`，只操作已有的 pandas datetime Series 或 DatetimeIndex，不经过 `Time`，也不接受数字或字符串时间。分段前排序，只有相邻时间差严格大于阈值才开启新区间，返回的端点保留输入时区。

## 边界和旧接口迁移

`None`、`NaT`、非有限数字、布尔数、复数、掩码时间、无效日期和超出 pandas 纳秒范围的日期均拒绝。pandas 无法保存 UTC 闰秒，因此来自 Astropy 或 SpacePy 的闰秒输入会报错，不替换成相邻时刻。Astropy 提供科学类型识别和转换，SpacePy 在需要识别其输入或输出 Ticktock 时导入。普通输入的时间计算直接使用 pandas。

`normalize_days` 已删除，取天使用 `t.floor("D")`，范围展开使用 `pd.date_range`。旧的 `normalize_time`、`normalize_datetime`、`convert_time`、`time_to_timestamp`、`tick`、`dt_2_str`、`is_scalar` 和 `normalize_trange` 已删除，`libs.utils` 中的相应旧导出也已移除。迁移方式如下：

| 旧调用 | 新调用 |
|---|---|
| `normalize_time(value)` / `normalize_datetime(values)` | `Time(value).datetime` |
| `normalize_time(value, numeric="julian")` | `Time(value, format="jd").datetime` |
| `time_to_timestamp(value)` | `Time(value).to_timestamp()` |
| `tick(value)` | `Time(value).to_tick()` |
| `dt_2_str(value)` | `Time(value).to_str()` |
| `normalize_days(...)` | `t.floor("D")` |
| `normalize_trange(trange)` | `Time(trange).datetime`，范围入口检查两个端点 |
| 时间标量判断 | `Time(value).isscalar` |

旧自动识别的毫秒或 JD 输入必须补上对应 `format`。所有新入口都统一 UTC，已不再保留来源时区。

使用 `.datetime` 取得原生 pandas 对象进行比较和索引；普通 pandas 属性和方法也可以直接通过 Time 实例访问：

```python
t = Time("2025-10-10 00:00:00").datetime
tl = Time(["2025-10-10", "2025-10-11"]).datetime
mask = t < tl
```

