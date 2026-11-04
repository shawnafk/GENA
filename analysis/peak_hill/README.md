# GENA ENA peak/hill 与 THEMIS 离子投掷角对比

一次性的对比分析，用来检查 2025-11-03/05 两个 ENA 时段的离子投掷角分布是否不同，以及这种不同能否解释 ENA 计数率的差异。

## 运行

```bash
# 主分析：GENA-A + THEMIS-B（ESA/ESD/FGM），输出 out/fig1..fig4 与 out/summary.json
python analysis/peak_hill/ena_pa_comparison.py

# 近地交叉检查：THEMIS-A（ESD 只覆盖到 2025-11-04），输出 out/tha_crosscheck.json
python analysis/peak_hill/tha_crosscheck.py
```

只算数值不画图：`python analysis/peak_hill/ena_pa_comparison.py --no-figures`。

## 依赖与数据

仓库内的 `load`、`gena`、`common`、`libs` 需要能从仓库根目录导入；GENA pqt 与 THEMIS 数据目录由 `FILE_BASE` / `common.config` 解析。THEMIS ESD 若本地缺失，脚本会从 `themis.ssl.berkeley.edu` 下载到 `data/`（THEMIS-A 到 `data_tha/`）。

如果当前环境中 pySPEDAS 因 `~/.spacepy` 不可写而导入失败，可在运行时把 `HOME` 指向可写目录，例如：

```bash
HOME=$PWD/analysis/peak_hill/_home python analysis/peak_hill/ena_pa_comparison.py
```

## 口径

- GENA 计数率是 `l2f_a` 全部事件的墙钟计数，ENA 谱用几何因子 1，只能作相对比较。
- ESD 只用 `valid==1`、`data_quality==0`、与 `peif_magf` 匹配 10 s 以内的记录；`look` 按粒子速度方向处理（与 `peif_velocity_dsl` 做过余弦校验）。
- 投掷角增益是 20°×10° 视场内的通量除以同能量全向平均，按每条 ESD 记录自己的 GENA 视线方向计算后再做窗口平均。
- 用户给定的 T1/T2 窗口在 11-04 09:00–23:00 重叠，`FOCUS` 里另给了两个不重叠的聚焦窗口，结论对两组窗口一致。

结论与数字见 `summary.md`。
