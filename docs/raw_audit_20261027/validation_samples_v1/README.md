# v1 验证样本包

`sample.raw` 是确定性的 6296 字节输入：4 字节前缀、一个含 A/B 事件的合格包、同一包重复一次、一个校验和损坏的包，以及一个只有 100 字节的尾部截断包。每个合格包仍包含 62 个槽位，其中 60 个是结构性全零槽。

运行：

```python
from gena.raw.pipeline import process_raw
process_raw("sample.raw", "package", chunk_packets=1, by_day=False)
```

将输出 `manifest.json` 中的包统计、处理漏斗和产品行数与 `expected.json` 比较。`chunk_packets=1` 特意覆盖跨读块包头；把它改为 `2` 或 `8` 时，包 UUID、产品行数和逐列数值应保持一致。Parquet 行组边界会随处理块变化，因此不能用文件字节哈希比较不同块大小的输出。样本只验证程序的可重复性和拒绝路径，不模拟真实异常事件的发生率或真值。
