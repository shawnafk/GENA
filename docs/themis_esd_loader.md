# THEMIS ESA ESD 读取说明

`load_themis()` 支持 THEMIS ESA 的独立 ESD（ESA 3-D distribution）L2 产品。ESD 是角分辨的三维速度空间观测：数据数组的三个维度依次为观测时间、能量通道和固角方向 bin。它与普通 `product="esa"` 返回的 `peir_en_eflux` 不同；后者只有时间和能量维，已经对方向积分，不能用于 pitch-angle distribution。

## 读取

```python
from load import load_themis

trange = ["2025-12-04T14:00:00", "2025-12-04T16:00:00"]
esd = load_themis(trange, "esd", probe="b", datatype="peif", download="local_first")
```

`product="esd"` 只支持 `level="l2"`。`datatype` 传给 pySPEDAS；常用值为 `peif`（full-mode ion）、`peef`（full-mode electron）、`peir`（reduced-mode ion）、`peer`（reduced-mode electron）和 `peib`（burst-mode ion）。默认值为 `peif`。和其他 THEMIS 读取器一样，默认 `download=False` 只读本地缓存；使用 `download=True`、`"local_first"` 或 `"refresh"` 才会访问远程数据。

默认读取五个 ESD CDF 原生变量：`valid`、`data_quality`、`bins`、`eff` 和 `eflux`。与常规 ESA 变量不同，这些名称不带 `tha_`、`thb_` 等探测器前缀；如果传入 `vars`，也必须使用这些裸名称，例如 `vars=["eflux", "bins"]`。

## 返回约定

默认返回 `TimeSeriesDataFrame`。每一行对应 ESD 的一个观测时间；`eflux`、`bins`、`eff` 单元格内均为形状 `(energy-channel, solid-angle-bin)` 的 NumPy 数组，因此整列拼接后的物理形状为 `(time, energy-channel, solid-angle-bin)`。pySPEDAS 的 `v1`、`v2` 坐标轴不会被丢弃，分别保存在 `eflux_v1`、`eflux_v2`（以及相同变量的 companion 列）中。

```python
import numpy as np

cube = np.stack(esd["eflux"].to_numpy())
bin_ok = np.stack(esd["bins"].to_numpy()).astype(bool)
quality_ok = (esd["valid"] == 1) & (esd["data_quality"] == 0)
cube[~quality_ok.to_numpy()] = np.nan
cube = np.where(bin_ok, cube, np.nan)
```

`eflux` 已是 ESD L2 校准后的微分能通量；不要再为了得到物理通量手动除以 `eff`。`bins == 1` 表示该能量–方向 bin 有效，`valid == 1` 和 `data_quality == 0` 是本项目图件使用的保守整条分布筛选。更细的质量位掩码解释应以相应 CDF 的元数据和 THEMIS ESA 数据质量说明为准。

若需要保留原始 pySPEDAS 字典而非 DataFrame，可使用：

```python
raw = load_themis(trange, "esd", probe="b", datatype="peif", return_type="dict")
eflux = raw["eflux"]
time = eflux["x"]
cube = eflux["y"]
energy_channel = eflux["v1"]
solid_angle_bin = eflux["v2"]
```

ESD 的 `v1` 和 `v2` 是 CDF 的通道索引，不是直接可用于 PA 的物理能量或方向。计算 pitch angle 时，还需读取同一 ESD CDF 的 `energy`、`theta`、`phi`、`domega` 及模式索引，并用常规 ESA L2 `peif_magf` 在 ESA DSL 坐标系给出磁场方向。`1.Paper_BBF/plot_themis_esd_pitch_angle_20251204.py` 展示了完整的 ESD-only PA 分布读取和绘图流程。
