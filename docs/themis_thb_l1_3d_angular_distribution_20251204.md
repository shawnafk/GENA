# 读取 THB L1 三维角分布并估计 ENA 方向性增强

本文记录 2025-12-04 14:00--16:00 UTC 的 THB ESA 离子数据如何被读取为能量--角度分布，并说明本次用于 GENA 快速 ENA 计数估计的方向性增强计算。这里的“3D”指速度空间中的能量维度加两个方向维度，并非空间三维粒子分布。

## 数据源与变量

本次使用的本地 CDF 文件如下：

| 产品 | 文件 | 变量 | 作用 |
| --- | --- | --- | --- |
| THB ESA L1 | `1.Paper_BBF/P3-forward-modeling/ena_simulation/SPEDAS/themis_data/thb/l1/esa/2025/thb_l1_esa_20251204_v02.cdf` | `thb_eir_050x24_time`, `thb_eir_050x24` | 角分辨离子分布 |
| THB ESA L2 | `/Volumes/Rsync/ENA/SPEDAS/themis_data/thb/l2/esa/2025/thb_l2_esa_20251204_v01.cdf` | `thb_peir_time`, `thb_peir_en_eflux`, `thb_peir_en_eflux_yaxis` | 绝对能谱归一化 |
| THB SST L2 | `/Volumes/Rsync/ENA/SPEDAS/themis_data/thb/l2/sst/2025/thb_l2_sst_20251204_v01.cdf` | `thb_psif_time`, `thb_psif_en_eflux`, `thb_psif_en_eflux_yaxis` | 24.6 keV 以上的全向谱 |

L1 数组 `thb_eir_050x24` 的原始形状为 `(time, 50, 24)`：每个记录有 50 个方向 bin 和 24 个 ESA 能道。该文件在本时段的采样间隔约为 4.304 s。L1 和 L2 `peir` 时间轴逐样本匹配；本次检查到的最大时间差约为 70 微秒。

## 读取与解压

L1 数据采用 THEMIS 的 8-bit 平方根压缩表，不能将原始 `uint8` 直接当作计数使用。统一加载入口导出的 `themis_decompress()` 实现了官方反查表，`esa_mode2_angles()` 返回 reduced-mode-2 的官方视向与每个方向 bin 的固角。

```python
import numpy as np
from cdflib import CDF

from load import esa_mode2_angles, themis_decompress

l1_path = (
    "1.Paper_BBF/P3-forward-modeling/ena_simulation/SPEDAS/"
    "themis_data/thb/l1/esa/2025/thb_l1_esa_20251204_v02.cdf"
)

with CDF(l1_path) as cdf:
    l1_time = cdf.varget("thb_eir_050x24_time")
    compressed = cdf.varget("thb_eir_050x24")

# 单个时刻：原始形状为 (angle, energy)，解压后转为 (energy, angle)。
counts_ea = themis_decompress(compressed[0]).T.astype(float)
look_ea, omega_ea = esa_mode2_angles()

assert counts_ea.shape == (24, 50)
assert look_ea.shape == (24, 50, 3)
assert omega_ea.shape == (24, 50)
```

`omega_ea` 的每一个能道沿 50 个 bin 求和为 (4\pi) sr。三个不同大小的单 bin 固角是 0.15028、0.42505 和 0.47828 sr。相比之下，GENA 20° x 10° 的矩形 FOV 为 0.06054 sr，因此 GENA FOV 小于所有 THB ESA 角 bin。

## 用 L2 能谱给 L1 角分布定标

L1 提供角向相对形状，L2 `peir_en_eflux` 提供绝对差分能量通量。对每一个能道 (E)，先计算 L1 的固角加权平均：

\[
\bar C(E)=\frac{\sum_a C(E,a)\,\Omega(E,a)}{\sum_a\Omega(E,a)}.
\]

定义无量纲角向对比度：

\[
R(E,a)=\frac{C(E,a)}{\bar C(E)}.
\]

则对应方向的绝对粒子差分通量可写为：

\[
j(E,a)=\frac{F_{\rm L2}(E)}{E}\,R(E,a),
\]

其中 (F_{\rm L2}) 的单位是 `eV/(cm2 s sr eV)`，能量 (E) 的单位是 eV。故 (j) 的单位为 `cm-2 s-1 sr-1 eV-1`，可直接进入 `libs.plasma.sigma_HH(E / 1000)` 与 ENA 计数积分。

```python
mean_counts = (counts_ea * omega_ea).sum(axis=1) / omega_ea.sum(axis=1)
relative_angular_flux = counts_ea / mean_counts[:, None]

# eflux_l2_e 的单位为 eV/(cm2 s sr eV)，energy_e 为 eV。
particle_flux_ea = (
    eflux_l2_e[:, None] / energy_e[:, None] * relative_angular_flux
)
```

实际计算时只使用 4--24.6 keV 的 ESA 角分辨能道；L1 没有提供 SST 能段的方向信息。

## 为什么按 1 分钟累积 L1 计数

单个 4.3 s L1 记录在 4--25 keV 的部分高能角 bin 中只有零到几颗解压计数。例如，若某能道只有一颗计数，`max(C)/mean(C)` 会因泊松涨落出现二十至八十倍的伪“束流”。这不能代表一个 60 s GENA 积分。

因此对每个整分钟先累积解压后的 L1 计数：

```python
counts_ea_1min = themis_decompress(compressed[index_in_minute]).sum(axis=0).T
```

本次每个完整分钟包含 11--14 个 L1 样本。随后用该分钟内平均的 ESA/SST L2 能谱作权重。这样得到的是与 `qk_estimate_ena_counts(..., accumulation_time_s=60)` 时间定义一致的方向性估计。

## 单一视线的增强计算

不能让每个能道各取自己的最亮方向，因为 GENA 在一个积分内只有一个 FOV 指向。具体做法是：

1. 以官方 `look_ea` 中的方向为候选 GENA 视线；
2. 对每个能道，选择与候选视线夹角最小的 ESA 方向 bin；
3. 用同一个候选视线跨全部能道计算加权增强，并在候选视线中取最大值。

ESA 覆盖能段的最佳一致视线增强为：

\[
K_{\rm ESA}=\max_{\hat{\ell}}
\frac{\int \sigma_{\rm HH}(E)j_{\rm omni}(E)
R[E,a(E,\hat{\ell})],dE}
{\int \sigma_{\rm HH}(E)j_{\rm omni}(E),dE}.
\]

其中 (a(E,\hat{\ell})) 是该能道最接近候选视线 \(\hat{\ell}\) 的 ESA bin，且 `sigma_HH` 调用时能量由 eV 转为 keV。

对于总的 4--200 keV 计数，本次将没有 L1 角信息的 SST 部分保守地设为各向同性。若 (f_{\rm ESA}) 是 σj 积分中 4--24.6 keV 的比例，则：

\[
K_{4-200}=1+f_{\rm ESA}(K_{\rm ESA}-1).
\]

## 本时段结果

以下结果均是每分钟累积后得到的“ESA 可分辨尺度上的最佳指向”结果；它们不是 GENA-A 实际指向的预测。

| 量 | 5% | 中位数 | 95% | 最大值 |
| --- | ---: | ---: | ---: | ---: |
| (K_{\rm ESA})，4--24.6 keV | 3.48 | 4.86 | 7.67 | 8.27 |
| (K_{4-200})，SST 各向同性 | 1.75 | 3.87 | 7.05 | 7.68 |
| 每能道独立取最亮 bin 的宽松上限 (K_{4-200}) | 2.04 | 4.29 | 7.14 | 8.36 |

15:24 UTC 的一分钟具有最大的全向谱 ENA 基线预测：0.986 counts/min。L1 角分布给出的该分钟增强为 4.31，因此定向预测为 4.25 counts/min。全时段联合最大的分钟预测出现在 14:44 UTC：基线 0.719 counts/min，增强 6.66，得到 4.79 counts/min。

## 适用范围与限制

- 这一步只使用了 THB 的局地离子 VDF；它检验的是“若源区离子 VDF 具有 THB 所见的角向集中度，方向性可将快估提高多少”。它不证明 GENA-A 的 LOS 上离子分布与 THB 局地分布相同。
- 结果已经要求所有能道采用同一个视线，避免能道间不物理的最大值拼接。
- THB ESA 的最小角 bin 仍约为 GENA FOV 的 2.5 倍。因此，若存在远小于 0.15 sr 的窄束流，L1 会把它在 bin 内平均掉；本分析无法测量或排除这种子-bin 束流。
- 尚未将 GENA-A 的姿态、探测器 FOV 和 THB DSL 坐标做同一时刻的坐标转换。要得到实际 GENA-A LOS 的增强，应先从 `l2g_a.pqt` 的姿态构造 GENA 视线，再转换至可与 THB 角分布比较的共同坐标系，而不是在 THB 全天空中取最大方向。
- 结果沿用了原先 `tvis.flux` 的“直接使用 L2 谱”比较方式；最终科学分析应明确 `thb_peir_data_quality` 的位掩码策略，并检查质量筛选对峰值分钟的影响。
