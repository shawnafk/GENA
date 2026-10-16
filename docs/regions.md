# 跨任务区域判定

外部数据加载放在 `load/`，区域与流结构判定放在 `libs.regions`，区域色带绘图放在 `plots.advance.plot_region_code`。原 `themis.regions` 已移除，无兼容入口。

`classify_liuzzo`、`classify_kiehas`、`classify_runov` 接收 `TimeSeriesDataFrame`，保留原时间同步规则与判据。默认识别 THEMIS 列名；其他数据源显式传入 `columns` 映射，不需要把数据伪装成 THEMIS。结果来源取自输入的 `attrs['source']`，未提供时为 `None`。

例如，提前加载 MMS FPI 离子矩之后，可构造标量离子温度并调用 Kiehas 判据。以下列名为示例，请以实际加载结果为准：

```python
from libs.regions import classify_kiehas

moments = fpi_data.copy()
moments["ion_temperature_ev"] = (
    moments["mms1_dis_temppara_fast"]
    + 2 * moments["mms1_dis_tempperp_fast"]
) / 3
result = classify_kiehas(moments, columns={
    "density": "mms1_dis_numberdensity_fast",
    "temperature": "ion_temperature_ev",
})
```

这一温度组合要求输入平行和垂直温度具有相同单位和物理定义。密度单位为 cm⁻³，温度为 eV；Liuzzo 判据的磁场必须为 GSM XYZ 分量（nT），位置回调接收 Unix 秒并返回 GSM 位置（R_E）。`columns={'density': name, 'magnetic': [bx, by, bz]}` 显式指定相应字段。数据缺失或时间间隔过大时仍遵循原同步规则。

Runov 判据需要第一张表中的离子密度、离子温度、电子密度和电子温度，以及第二张表中的高能离子密度和温度；映射键分别为 `ion_density`、`ion_temperature`、`electron_density`、`electron_temperature`、`hot_density`、`hot_temperature`。缺少高能矩时不能把普通 FPI 矩重复传入来替代；仪器能量覆盖与原 ESA/SST 组合的可比性需要科学判断。

`classify_plasma_sheet_flows` 原本使用 `beta_i`、`bx_nt`、`vx_kms` 等物理量列，可用于已准备好这些量的不同任务数据。可选的 `esa_l1` / `esa_l2` 方向分布检查只适用于 THEMIS reduced ESA，不是 MMS VDF 检查。

区域色带接受含 `region_code` 的分类结果：

```python
from plots import Panel, overview, plot_region_code

figure, axes = overview([
    Panel(plot_region_code, kwargs={"result": region_result, "label": "MMS1"}),
])
```

Liuzzo 返回 `region_code`，Kiehas 和 Runov 返回布尔筛选结果，不直接等同于完整的区域分类。此次重组与合成数据测试验证的是接口复用和数值行为保持一致，没有重新验证原文判据或完成真实 MMS 数据的科学验证。
