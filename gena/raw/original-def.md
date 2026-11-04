# GENA raw 原始字段与转换定义

本文按当前工作区中的 `gena/raw/packet.py`、`conversion.py`、`decode.py`、`structure.py`、`utils.py` 和 `gena/constants.py` 逐字段记录从原始 2064 字节包到 PDF、L0（`df0`）、L1、L1A、L2 的定义与类型转换。文中 `>` 表示大端（big-endian），`u1/i2/u4/i4` 分别表示无符号 8 位、有符号 16 位、无符号 32 位、有符号 32 位；`声明类型` 指 `structure.py` 中的 `*_type`，`实际类型` 指当前代码在 pandas 中实际生成的 dtype。

转换主链是：

```text
raw 2064-byte packet
  ├─ packets_2_pdf()  -> PDF 工程参数/姿态/位置
  └─ packets_2_df0()  -> L0 事件槽（df0）
        └─ l0_2_l1()  -> L1（TDC ns、EP）
              └─ l1_2_l1a() -> L1A（八通道事件 + PDF 插值）
                    └─ l1_2_l2() -> L2（TOF、位置、种类、能量）
```

当前工作区有两个阅读前置问题，检查字段时请先留意：`gena/raw/__init__.py:3` 仍然 `from .pipeline import process_raw`，但当前目录没有 `pipeline.py`，所以 `import gena.raw` 会先失败；`decode.py:275-276` 调用 `rep.l1a_2_l1b`，而 `repair.py:487` 定义的是 `repair_l1a`，并且 `repair_l1a` 会直接 `raise NotImplementedError`。这不影响直接阅读源码，但如果要运行字段检查，需要先处理这两处。

## 1. 2064 字节包的二进制布局

包结构定义在 [structure.py:42-53](structure.py#L42-L53)，总长 2064 字节：

| 字段 | 相对偏移 | 字节数 | numpy dtype | 含义 | 定义 |
| --- | ---: | ---: | --- | --- | --- |
| `header` | 0 | 2 | `u1[2]` | 同步字，通常为 `E2 25` | [structure.py:43](structure.py#L43) |
| `seq1` | 2 | 2 | `>u2` | APID，正常值 `0x0C4E` | [structure.py:44](structure.py#L44) |
| `seq2` | 4 | 2 | `>u2` | 高 2 位 `BM`，低 14 位 `BX` | [structure.py:45](structure.py#L45) |
| `dataLen` | 6 | 2 | `>u2` | 数据长度，正常值 2055 | [structure.py:46](structure.py#L46) |
| `time_s` | 8 | 4 | `>u4` | 包 UTC 秒计数（基准 1609430400） | [structure.py:47](structure.py#L47) |
| `time_ms` | 12 | 2 | `>u2` | 包内毫秒，正常 `< 1000` | [structure.py:48](structure.py#L48) |
| `ep1` | 14 | 32 | `ep1_struct` | 高压、监视量、温度、事件计数 | [structure.py:49](structure.py#L49) |
| `ep2` | 46 | 32 | `ep2_struct` | 平台组、姿态/位置、指向、参数表 | [structure.py:50](structure.py#L50) |
| `events` | 78 | 1984 | `event_struct[62]` | 62 个 32 字节事件槽 | [structure.py:51](structure.py#L51) |
| `verification_sum` | 2062 | 2 | `>u2` | 字节和校验字 | [structure.py:52](structure.py#L52) |

包级常量和校验条件：

| 名称 | 值/条件 | 定义 |
| --- | --- | --- |
| `PACKET_MAGIC` | `b'\xe2\x25\x0c\x4e'` | [packet.py:7](packet.py#L7) |
| `PACKET_SIZE` | `packet_struct.itemsize` = 2064 | [packet.py:8](packet.py#L8) |
| `PACKET_DATA_LENGTH` | `2064 - 8 - 1` = 2055 | [packet.py:9](packet.py#L9) |
| `sync` | `header == [0xE2, 0x25]` | [packet.py:485](packet.py#L485) |
| `apid` | `seq1 == 0x0C4E` | [packet.py:486](packet.py#L486) |
| `length` | `dataLen == 2055` | [packet.py:487](packet.py#L487) |
| `time_ms` | `< 1000` | [packet.py:488](packet.py#L488) |
| `checksum` | `sum(raw[8:2062]) & 0xffff == verification_sum` | [packet.py:483-489](packet.py#L483-L489) |

`ep1_struct` 定义在 [structure.py:5-19](structure.py#L5-L19)，相对 `ep1` 基址 14 的偏移如下：

| 字段 | 相对偏移 | 字节数 | dtype | 含义 | 定义 |
| --- | ---: | ---: | --- | --- | --- |
| `HV_N` | 0 | 2 | `>i2` | 负高压原始码 | [structure.py:6](structure.py#L6) |
| `HV_P` | 2 | 2 | `>u2` | 正高压原始码 | [structure.py:7](structure.py#L7) |
| `HV_C` | 4 | 2 | `>i2` | 加速/偏转高压原始码 | [structure.py:8](structure.py#L8) |
| `HV_G` | 6 | 2 | `>i2` | guard 高压原始码 | [structure.py:9](structure.py#L9) |
| `HV_A` | 8 | 2 | `>i2` | A 高压原始码 | [structure.py:10](structure.py#L10) |
| `HV_B` | 10 | 2 | `>i2` | B 高压原始码 | [structure.py:11](structure.py#L11) |
| `monitors` | 12 | 8 | `int8[8]` | 8 路模拟监视量原始码 | [structure.py:12](structure.py#L12) |
| `temps` | 20 | 3 | `u1[3]` | A/B/PCB 温度原始码 | [structure.py:13](structure.py#L13) |
| `stat` | 23 | 1 | `u1` | 状态字节 | [structure.py:14](structure.py#L14) |
| `Evt_Num_Tot` | 24 | 2 | `>u2` | 总事件计数 | [structure.py:15](structure.py#L15) |
| `Evt_Num_Eff` | 26 | 2 | `>u2` | 有效事件计数 | [structure.py:16](structure.py#L16) |
| `Evt_Num_Aux` | 28 | 2 | `>u2` | 辅助事件计数 | [structure.py:17](structure.py#L17) |
| `instruction_counts` | 30 | 2 | `u1[2]` | 指令计数，当前未进入产品 | [structure.py:18](structure.py#L18) |

`ep2_struct` 定义在 [structure.py:21-28](structure.py#L21-L28)，相对 `ep2` 基址 46 的偏移如下：

| 字段 | 相对偏移 | 字节数 | dtype | 含义 | 定义 |
| --- | ---: | ---: | --- | --- | --- |
| `version` | 0 | 1 | `u1` | 版本 | [structure.py:22](structure.py#L22) |
| `err` | 1 | 1 | `u1` | 错误字 | [structure.py:23](structure.py#L23) |
| `stat` | 2 | 1 | `u1` | 平台组别：0=时钟/天线，1=四元数/月心位置 | [structure.py:24](structure.py#L24) |
| `flex` | 3 | 24 | `u1[24]` | 复用字节区，按 `stat` 解释 | [structure.py:25](structure.py#L25) |
| `args_addr` | 27 | 1 | `u1` | 参数表地址 | [structure.py:26](structure.py#L26) |
| `args_value` | 28 | 4 | `u1[4]` | 参数值原始字节 | [structure.py:27](structure.py#L27) |

`event_struct` 定义在 [structure.py:30-40](structure.py#L30-L40)，第一个事件绝对偏移 78，槽间步长 32 字节：

| 字段 | 相对偏移 | 字节数 | dtype | 含义 | 定义 |
| --- | ---: | ---: | --- | --- | --- |
| `flag` | 0 | 1 | `u1` | bit7=A/B，bit6=有效/辅助，低 4 位子模式 | [structure.py:31](structure.py#L31) |
| `time_s` | 1 | 2 | `>u2` | 事件内部秒计数 | [structure.py:32](structure.py#L32) |
| `time_us` | 3 | 2 | `>u2` | 事件内部微秒/计数码，代码按 20 µs 解释 | [structure.py:33](structure.py#L33) |
| `se_amp` | 5 | 3 | `u1[3]` | 电子幅度三字节 | [structure.py:34](structure.py#L34) |
| `c0` | 8 | 8 | `>u2[4]` | 4 个主 TOF 测量字 | [structure.py:35](structure.py#L35) |
| `c1` | 16 | 4 | `>u2[2]` | 通道 1 的 start/end 字 | [structure.py:36](structure.py#L36) |
| `c3` | 20 | 4 | `>u2[2]` | 通道 3 的 start/end 字 | [structure.py:37](structure.py#L37) |
| `c5` | 24 | 4 | `>u2[2]` | 通道 5 的 start/end 字 | [structure.py:38](structure.py#L38) |
| `c7` | 28 | 4 | `>u2[2]` | 通道 7 的 start/end 字 | [structure.py:39](structure.py#L39) |

每个 TDC 字都是 16 位：高 3 位是通道号，低 13 位是测量值。当前代码只取低 13 位，没有保存也没有校验高 3 位，见 [packet.py:554-571](packet.py#L554-L571)。

## 2. PDF：包级工程参数、姿态、位置与参数表

`packets_2_pdf()` 在 [packet.py:328-460](packet.py#L328-L460)。一行对应一个包。`pdf_type` 的声明类型在 [structure.py:56-109](structure.py#L56-L109)。

### 2.1 包标识、时间与平台状态

| PDF 列 | 来源与表达式 | 单位/比例 | 实际 dtype | `pdf_type` 声明 | 定义 |
| --- | --- | --- | --- | --- | --- |
| `BM` | `seq2 >> 14` | 0..3 | `uint8` | `uint8` | [packet.py:335](packet.py#L335) |
| `BX` | `seq2 & 0x3fff` | 0..16383 | `uint16` | `uint16` | [packet.py:337](packet.py#L337) |
| `Timestamp` | `1609430400 + time_s + time_ms * 0.001` | 秒（UTC） | `float64` | `float64` | [packet.py:346](packet.py#L346) |
| `PlatformGroup` | `ep2.stat` | 0/1/其他 | `uint8` | 未声明 | [packet.py:396](packet.py#L396) |
| `PosePresent` | `ep2.stat == 1` | bool | `bool` | 未声明 | [packet.py:397](packet.py#L397) |

`Timestamp` 的基准 `1609430400` 对应 2020-12-31 16:00:00 UTC，即 2021-01-01 00:00:00 北京时间，见 [packet.py:339-346](packet.py#L339-L346)。

### 2.2 高压与监视量

| PDF 列 | 来源与表达式 | 比例 | 实际 dtype | `pdf_type` 声明 | 定义 |
| --- | --- | --- | --- | --- | --- |
| `HV_N` | `ep1.HV_N * 50000/32768` | V | `float64` | `float32` | [packet.py:350-351](packet.py#L350-L351) |
| `HV_P` | `ep1.HV_P * 50000/32768` | V | `float64` | `float32` | [packet.py:350-352](packet.py#L350-L352) |
| `HV_C` | `ep1.HV_C * 10000/32768` | V | `float64` | `float32` | [packet.py:354-355](packet.py#L354-L355) |
| `HV_G` | `ep1.HV_G * 10000/32768` | V | `float64` | `float32` | [packet.py:354-356](packet.py#L354-L356) |
| `HV_A` | `ep1.HV_A * 10000/32768` | V | `float64` | `float32` | [packet.py:354-357](packet.py#L354-L357) |
| `HV_B` | `ep1.HV_B * 10000/32768` | V | `float64` | `float32` | [packet.py:354-358](packet.py#L354-L358) |
| `A_N12` | `monitors[0] * 50/128` | V | `float64` | `float32` | [packet.py:360-362](packet.py#L360-L362) |
| `A_N5` | `monitors[1] * 20/128` | V | `float64` | `float32` | [packet.py:364-365](packet.py#L364-L365) |
| `A_3_3` | `monitors[2] * 10/128` | V | `float64` | `float32` | [packet.py:367-368](packet.py#L367-L368) |
| `D_3_3` | `monitors[3] * 10/128` | V | `float64` | `float32` | [packet.py:367-369](packet.py#L367-L369) |
| `HV_C_I` | `monitors[4] * 10/128` | 电流监视 | `float64` | `float32` | [packet.py:367-371](packet.py#L367-L371) |
| `HV_N_I` | `monitors[5] * 10/128` | 电流监视 | `float64` | `float32` | [packet.py:367-372](packet.py#L367-L372) |

`monitors` 在 [structure.py:12](structure.py#L12) 声明为 `np.int8`，但 `50/128`、`20/128`、`10/128` 都是按 0..255 满量程设计的正系数。当前实现会把最高位为 1 的原始字节解释成负数；例如原始字节 200 时 `A_N12 = -21.875`，按无符号应为 78.125。这一项需要优先核对。

### 2.3 温度与事件计数

| PDF 列 | 来源与表达式 | 实际 dtype | `pdf_type` 声明 | 定义 |
| --- | --- | --- | --- | --- |
| `Temp_A` | `temps[0]` → 见下式 | `float32` | `float32` | [packet.py:376-380](packet.py#L376-L380) |
| `Temp_B` | `temps[1]` → 见下式 | `float32` | `float32` | [packet.py:376-381](packet.py#L376-L381) |
| `Temp_PCB` | `temps[2]` → 见下式 | `float32` | `float32` | [packet.py:376-382](packet.py#L376-L382) |
| `Evt_Num_Tot` | `Evt_Num_Tot.byteswap().view(byteorder)` | `uint16` | `uint32` | [packet.py:383-386](packet.py#L383-L386) |
| `Evt_Num_Eff` | 同上 | `uint16` | `uint32` | [packet.py:383-387](packet.py#L383-L387) |
| `Evt_Num_Aux` | 同上 | `uint16` | `uint32` | [packet.py:383-388](packet.py#L383-L388) |

温度换算在 [packet.py:376-379](packet.py#L376-L379)：

```text
temps_raw = temps * 0.078125
log_arg   = temps_raw * 10000 / (5000 * (4.97 - temps_raw))
log_arg   = log_arg if log_arg > 0 else NaN
Temp      = 4100 * 298.15 / (298.15 * ln(log_arg) + 4100) - 273.15
```

`Temp_A/B/PCB` 最后分别 `astype(np.float32)`；`log_arg <= 0` 或非有限值产生 NaN。

`Evt_Num_*` 的 `byteswap().view(dtype.newbyteorder())` 是为了兼容 `np.concatenate` 可能把大端结构化数组转成原生字节序的情况；实际输出类型是 `uint16`，不是 `pdf_type` 声明的 `uint32`。

### 2.4 姿态、位置与指向

| PDF 列 | 来源与表达式 | 单位/比例 | 实际 dtype | `pdf_type` 声明 | 定义 |
| --- | --- | --- | --- | --- | --- |
| `MJ2000_q1` | `pose_present ? int32(flex[0:4]) * 2^-30 : 0` | 四元数 x | `float32` | `float32` | [packet.py:402-403](packet.py#L402-L403) |
| `MJ2000_q2` | `pose_present ? int32(flex[4:8]) * 2^-30 : 0` | 四元数 y | `float32` | `float32` | [packet.py:402-404](packet.py#L402-L404) |
| `MJ2000_q3` | `pose_present ? int32(flex[8:12]) * 2^-30 : 0` | 四元数 z | `float32` | `float32` | [packet.py:402-405](packet.py#L402-L405) |
| `MJ2000_Star_X` | `pose_present ? int32(flex[12:16]) : 0` | m，月心 J2000 | `int32` | `int32` | [packet.py:407-408](packet.py#L407-L408) |
| `MJ2000_Star_Y` | `pose_present ? int32(flex[16:20]) : 0` | m | `int32` | `int32` | [packet.py:407-409](packet.py#L407-L409) |
| `MJ2000_Star_Z` | `pose_present ? int32(flex[20:24]) : 0` | m | `int32` | `int32` | [packet.py:407-410](packet.py#L407-L410) |
| `Ant_Ang` | `pointing_present ? int32(flex[6:10]) * 1e-5 : 0` | deg | `float64` | `float32` | [packet.py:426-428](packet.py#L426-L428) |
| `Ant_Azi` | `pointing_present ? int32(flex[10:14]) * 1e-5 : 0` | deg | `float64` | `float32` | [packet.py:426-429](packet.py#L426-L429) |
| `Star_L` | `pointing_present ? int32(flex[14:18]) * 1e-5 : 0` | 视场角 | `float64` | `float32` | [packet.py:431-432](packet.py#L431-L432) |
| `Star_NZ` | `pointing_present ? int32(flex[18:22]) * 1e-5 : 0` | 视场角 | `float64` | `float32` | [packet.py:431-433](packet.py#L431-L433) |
| `Shadow` | `pointing_present ? flex[22] : 0` | 标志 | `uint8` | `float32` | [packet.py:434](packet.py#L434) |

这里的 `pose_present = (ep2.stat == 1)`，`pointing_present = (ep2.stat == 0)`，见 [packet.py:393-397](packet.py#L393-L397)。`pose_present` 为假时位置和四元数全部写 0；`pointing_present` 为假时指向字段全部写 0。

`view_as(..., '>i4')` 在 [packet.py:319-320](packet.py#L319-L320)，它把连续 4 字节按大端有符号 32 位读取；`q_dl = 2**-30` 在 [packet.py:402](packet.py#L402)，`a_dl = m_dl = 1e-5` 在 [packet.py:427-431](packet.py#L427-L431)。

### 2.5 参数表地址、参数值与阈值

| PDF 列 | 来源与表达式 | 实际 dtype | `pdf_type` 声明 | 定义 |
| --- | --- | --- | --- | --- |
| `ArgsAddress` | `ep2.args_addr` | `uint8` | 未声明 | [packet.py:437-440](packet.py#L437-L440) |
| `ArgsValue` | `uint32(ep2.args_value, big-endian)` | `uint32` | 未声明 | [packet.py:437-441](packet.py#L437-L441) |
| `Th_A_MCP` | `args_addr == 53`：`u32(args_value[::-1]) & 0xFFF` | `int16` | `int16` | [packet.py:448](packet.py#L448) |
| `Th_B_MCP` | `args_addr == 57`：同上 | `int16` | `int16` | [packet.py:449](packet.py#L449) |
| `Th_A_X` | `args_addr == 55`：同上 | `int16` | `int16` | [packet.py:451](packet.py#L451) |
| `Th_A_Y` | `args_addr == 59`：同上 | `int16` | `int16` | [packet.py:452](packet.py#L452) |
| `Th_B_X1` | `args_addr == 142`：`signed_7bits(args_value[1])` | `int16` | `int16` | [packet.py:454](packet.py#L454) |
| `Th_B_X2` | `args_addr == 142`：`signed_7bits(args_value[3])` | `int16` | `int16` | [packet.py:455](packet.py#L455) |
| `Th_B_Y1` | `args_addr == 146`：`signed_7bits(args_value[1])` | `int16` | `int16` | [packet.py:456](packet.py#L456) |
| `Th_B_Y2` | `args_addr == 146`：`signed_7bits(args_value[3])` | `int16` | `int16` | [packet.py:457](packet.py#L457) |
| `Th_B_G` | `args_addr == 150`：`int8(args_value[0])` | `int8` | `int8` | [packet.py:458](packet.py#L458) |

`ArgsValue` 使用 `view_as(paras, '>u4')`；`Th_A_*` 使用 `view_as(paras[..., ::-1], np.uint32)`，先把 4 字节反序再按原生 `uint32` 读取，最后 `& 0xFFF`。`signed_7bits` 在 [packet.py:304-314](packet.py#L304-L314)：

```text
bit6 = 符号位
value = raw & 0x3F
if bit6: value = -value
```

它是符号-幅值表示，不是补码：`0x40 -> 0`、`0x41 -> -1`、`0x7F -> -63`。如果协议中的 7 位阈值是补码，这里应为 `0x40 -> -64`、`0x7F -> -1`，需要按原始协议确认。

## 3. L0 事件槽（`df0`）

`packets_2_df0()` 在 [packet.py:514-582](packet.py#L514-L582)。一个包展开 62 行，`PID=1..62`；`UUID` 是整包字节的 SHA-256 再做 UUID5，相同内容包在不同 chunk 中保持同一 UUID，见 [packet.py:502-511](packet.py#L502-L511)。

| `df0` 列 | 来源与表达式 | 声明 dtype | 实际 dtype | 定义 |
| --- | --- | --- | --- | --- |
| `UUID` | `uuid5(namespace, sha256(packet_bytes))`，每包重复 62 次 | `S36` | `str` | [packet.py:524-525](packet.py#L524-L525) |
| `PID` | `tile(arange(1, 63), num_packets)` | `uint8` | `uint8` | [packet.py:527](packet.py#L527) |
| `BM` | `seq2 >> 14`，每包重复 62 次 | `uint8` | `uint8` | [packet.py:529](packet.py#L529) |
| `BX` | `seq2 & 0x3fff`，每包重复 62 次 | `uint16` | `uint16` | [packet.py:530](packet.py#L530) |
| `Timestamp` | `1609430400 + time_s + time_ms * 0.001`，每包重复 62 次 | `float64` | `float64` | [packet.py:531](packet.py#L531) |
| `is_empty_slot` | 32 字节事件槽全 0 | 未声明 | `bool` | [packet.py:534-535](packet.py#L534-L535) |
| `prob_flag` | `flag >> 7` | `uint8` | `uint8` | [packet.py:537](packet.py#L537) |
| `EvtLevel` | `(flag >> 6) & 1` | `uint8` | `uint8` | [packet.py:539](packet.py#L539) |
| `TimeInternal` | `time_s + time_us * 0.00002`，再 `astype(float32)` | `float32` | `float32` | [packet.py:545-549](packet.py#L545-L549) |
| `EP` | `((se_amp[1] << 8) \| se_amp[2]) & 0xFFF` | `uint16` | `uint16` | [packet.py:551-552](packet.py#L551-L552) |
| `ch_tof_0` | `c0[0] & 0x1FFF` | `uint16` | `uint16` | [packet.py:555-558](packet.py#L555-L558) |
| `ch_tof_1` | `c0[1] & 0x1FFF` | `uint16` | `uint16` | [packet.py:555-558](packet.py#L555-L558) |
| `ch_tof_2` | `c0[2] & 0x1FFF` | `uint16` | `uint16` | [packet.py:555-558](packet.py#L555-L558) |
| `ch_tof_3` | `c0[3] & 0x1FFF` | `uint16` | `uint16` | [packet.py:555-558](packet.py#L555-L558) |
| `ch_1s` / `ch_1e` | `c1[0]/c1[1] & 0x1FFF` | `uint16` | `uint16` | [packet.py:562-574](packet.py#L562-L574) |
| `ch_3s` / `ch_3e` | `c3[0]/c3[1] & 0x1FFF` | `uint16` | `uint16` | [packet.py:565-576](packet.py#L565-L576) |
| `ch_5s` / `ch_5e` | `c5[0]/c5[1] & 0x1FFF` | `uint16` | `uint16` | [packet.py:568-578](packet.py#L568-L578) |
| `ch_7s` / `ch_7e` | `c7[0]/c7[1] & 0x1FFF` | `uint16` | `uint16` | [packet.py:571-580](packet.py#L571-L580) |

`ch_*` 只保留低 13 位；`c0`、`c1`、`c3`、`c5`、`c7` 的高 3 位通道号没有进入 `df0`。`TimeInternal` 的 `0.00002` 表示 20 µs 计数，和字段名/注释中的“微秒”并不完全一致，而且结果被压到 float32。

`is_empty_slot` 只表示 32 字节事件槽全为 0；这种槽的 `flag` 也是 0，因此 `prob_flag=0`，会被 `l0_2_l1(df0, 'a')` 当作 A 探头行选中。当前 `l0_2_l1()` 不读取 `is_empty_slot`，所以 L1 中会保留这些空槽；L1A 的八通道掩码随后会把它们排除，但 L1 行数和 L1 产品本身仍包含空槽。

## 4. L0 → L1

`l0_2_l1()` 在 [conversion.py:42-79](conversion.py#L42-L79)。它先按 `prob_flag` 选择 A/B 探头，再把 TDC 原码转换成 ns，把 EP 原码转换成幅度值。

### 4.1 TDC 转换

每个 TDC 原码通过 `tdc_by_index()` 转换，见 [conversion.py:13-20](conversion.py#L13-L20)：

```text
code <= 94          -> 0.0
code > 94           -> (code - 94) * 0.081   [ns]
输出 dtype          -> float32
```

`0.081` 是 ns/code，`94` 是无效/无测量门限；这两个常量硬编码在 [conversion.py:14-19](conversion.py#L14-L19)。

### 4.2 EP 转换

`ep_by_index()` 在 [conversion.py:22-34](conversion.py#L22-L34)：

```text
code == 0 或 4095   -> 0.0
0 < code < 4095     -> (2044 - code) * 100/123
其他/非有限         -> NaN
输出 dtype          -> float32
```

`4095` 和 `0` 是遥测约定的无测量码；`2044`、`100/123` 是幅度标定常数，见 [conversion.py:29-33](conversion.py#L29-L33)。

### 4.3 探头与通道映射

`comb_xy_channel` 在 [conversion.py:50-58](conversion.py#L50-L58)：

| 条件 | `TX_1` | `TX_2` | `TY_1` | `TY_2` | 定义 |
| --- | --- | --- | --- | --- | --- |
| probe A（`prob_flag==0`） | `ch_5` | `ch_7` | `ch_1` | `ch_3` | [conversion.py:50-52](conversion.py#L50-L52) |
| probe B，`snid=1`（在轨） | `ch_1` | `ch_3` | `ch_5` | `ch_7` | [conversion.py:55-56](conversion.py#L55-L56) |
| probe B，`snid=2`（地面） | `ch_7` | `ch_5` | `ch_3` | `ch_1` | [conversion.py:57-58](conversion.py#L57-L58) |

所有 `*_s`/`*_e` 分别从同名字段的 `s`/`e` 版本转换。`snid=2` 的 B 探头 X/Y 通道相对 `snid=1` 是对调的；这一映射需要在轨/地面配置确认。

### 4.4 L1 输出字段

| L1 列 | 来源 | 公式/转换 | dtype | 定义 |
| --- | --- | --- | --- | --- |
| `UUID`、`PID`、`BM`、`BX`、`Timestamp`、`TimeInternal`、`EvtLevel` | `df0` 同名列 | 直接复制 | 同 `df0` | [conversion.py:61-63](conversion.py#L61-L63) |
| `Tof_0`..`Tof_3` | `ch_tof_0`..`ch_tof_3` | `tdc_by_index` | `float32` | [conversion.py:65-69](conversion.py#L65-L69) |
| `TX_1s`、`TX_1e`、`TX_2s`、`TX_2e` | `ch_5s/5e/7s/7e`（A）等 | `tdc_by_index` | `float32` | [conversion.py:73-75](conversion.py#L73-L75) |
| `TY_1s`、`TY_1e`、`TY_2s`、`TY_2e` | `ch_1s/1e/3s/3e`（A）等 | `tdc_by_index` | `float32` | [conversion.py:73-75](conversion.py#L73-L75) |
| `EP` | `df0.EP` | `ep_by_index` | `float32` | [conversion.py:77](conversion.py#L77) |

`df1_type` 的声明类型在 [structure.py:136-157](structure.py#L136-L157)；除 `UUID` 实际是 `str` 外，其余列基本按声明生成。

## 5. L1 → L1A

`l1_2_l1a()` 在 [conversion.py:82-130](conversion.py#L82-L130)。L1A 只保留八个位置计时通道（`TX/Y` 的 `1/2` 与 `s/e`）全部非零的事件，再用 PDF 插值 HV、姿态和位置。

| L1A 列 | 来源 | 转换 | 声明 dtype | 实际 dtype | 定义 |
| --- | --- | --- | --- | --- | --- |
| 基础列 `UUID`..`EP`、`Tof_0..3`、`TX_*`、`TY_*` | L1 同名列 | 八通道掩码后直接复制 | 同 `df1_type` | 同 L1 | [conversion.py:95-108](conversion.py#L95-L108) |
| `HV_C` | PDF `HV_C` | `get_para` 线性插值 | `float32` | `float64` | [conversion.py:114](conversion.py#L114) |
| `HV_G` | PDF `HV_G` | 同上 | `float32` | `float64` | [conversion.py:115](conversion.py#L115) |
| `MJ2000_q1` | PDF 四元数 | SLERP 后 `as_quat()[:,0]` | `float32` | `float64` | [conversion.py:118-120](conversion.py#L118-L120) |
| `MJ2000_q2` | 同上 | `as_quat()[:,1]` | `float32` | `float64` | [conversion.py:118-121](conversion.py#L118-L121) |
| `MJ2000_q3` | 同上 | `as_quat()[:,2]` | `float32` | `float64` | [conversion.py:118-122](conversion.py#L118-L122) |
| `MJ2000_Star_X` | PDF `MJ2000_Star_X` | `CubicSpline` 插值 | `int32` | `float64` | [conversion.py:125-126](conversion.py#L125-L126) |
| `MJ2000_Star_Y` | PDF `MJ2000_Star_Y` | 同上 | `int32` | `float64` | [conversion.py:125-127](conversion.py#L125-L127) |
| `MJ2000_Star_Z` | PDF `MJ2000_Star_Z` | 同上 | `int32` | `float64` | [conversion.py:125-128](conversion.py#L125-L128) |

八通道掩码在 [conversion.py:105](conversion.py#L105)：`(df1[check_cols] != 0).all(axis=1)`。注意 `NaN != 0` 为真，因此如果调用方传入带 NaN 的 TDC 列，NaN 通道会被当成“非零”保留；当前正常 L0→L1 路径不会产生 NaN TDC。L1A 的基础列直接来自 `df1.loc[mask, base_cols]`，dtype 跟随 L1；在正常 L0→L1 路径中 TDC/EP 是 `float32`，标识列是 `uint8/uint16`，`Timestamp` 是 `float64`。

HV 插值由 [utils.py:23-34](utils.py#L23-L34)：

```text
idx = (y != 0)
interp1d(x[idx], y[idx], bounds_error=False, fill_value="extrapolate")
```

也就是说，0 被当作缺失值，超出 PDF 覆盖范围会线性外推；插值样本来自当前 `packets_2_pdf` 的 PDF 行（`decode.py:196`），不是整个文件。

四元数插值由 [utils.py:59-168](utils.py#L59-L168)：先按 `(q1 != 0) & (q2 != 0) & (q3 != 0)` 选样本，再用 `w = sqrt(1 - q1^2 - q2^2 - q3^2)` 补实部，`Rot.from_quat` 后 `Slerp`。这一步会丢弃任何含零分量的合法四元数（例如 `(0.707, 0, 0)`）。

位置插值由 [utils.py:256-313](utils.py#L256-L313)：有效点 `>= 4` 时对 X/Y/Z 分别做 `CubicSpline`；不足 4 点时退回 `func_interp` 的三个独立线性插值。时间轴是包级 `Timestamp`。

`Timestamp` 在整个链路中都是包级时间：`packets_2_df0()` 把一个包的同一个 `Timestamp` 重复给 62 个事件（[packet.py:531](packet.py#L531)），`TimeInternal` 只是内部计数，未加到 `Timestamp` 上。因此同一包内所有事件共享一个时间戳，事件级时间分辨率等于包级时间分辨率（通常 1 s），而不是 `TimeInternal` 的 20 µs 计数分辨率。

`df1a_type` 的声明类型在 [structure.py:184-215](structure.py#L184-L215)；其中 `MJ2000_Star_*` 声明为 `int32`，实际插值结果是浮点。

## 6. L1A → L2

`l1_2_l2()` 在 [conversion.py:228-266](conversion.py#L228-L266)。L2 的派生量如下。

### 6.1 TOF

| L2 列 | 公式 | 单位 | dtype | 定义 |
| --- | --- | --- | --- | --- |
| `Tof_X` | `((TX_1e + TX_2e) - (TX_1s + TX_2s)) / 2` | ns | `float32` | [conversion.py:133-140](conversion.py#L133-L140) |
| `Tof_Y` | `((TY_1e + TY_2e) - (TY_1s + TY_2s)) / 2` | ns | `float32` | [conversion.py:138-140](conversion.py#L138-L140) |
| `Tof_M` | `Tof_1` | ns | `float32` | [conversion.py:141](conversion.py#L141) |
| `Tof` | `ENA_flytime(df2) + f_tof_ep(df2)` | ns | `float32` | [conversion.py:250](conversion.py#L250) |

`ENA_flytime()` 在 [conversion.py:180-201](conversion.py#L180-L201)：主 TOF 与 X/Y 都在 5 ns 内时取 `Tof_M`；否则只有在 X/Y 绝对差 `< 5 ns` 或比值在 `0.9..1.1` 时才取 `min(Tof_X, Tof_Y)`；其余为 NaN。阈值 `5`、`0.9`、`1.1` 硬编码在 [conversion.py:194-196](conversion.py#L194-L196)。

`f_tof_ep()` 在 [conversion.py:144-154](conversion.py#L144-L154)：

```text
E  = |HV_C - HV_G|                         [按数值当作 eV]
ve = sqrt(2E / MASS_E) * SPEED_C           [m/s]
te = (2 * 3.1 mm + 47 mm) / ve * 1e9       [ns]
```

`3.1e-3`、`47e-3`、`MASS_E`、`SPEED_C` 见 [conversion.py:146-153](conversion.py#L146-L153) 和 [../../common/constants.py:29-49](../../common/constants.py#L29-L49)。

### 6.2 粒子种类

`ENA_spec()` 在 [conversion.py:157-177](conversion.py#L157-L177)：

| `Tof` 条件 | `EP` 条件 | `Spec` 值 | 含义 |
| --- | --- | --- | --- |
| 非有限或 `<= 0` | 任意 | 3 | 无效 TOF |
| `< 28.79` | 任意 | 1 | 氢 |
| `> 161.55` | 任意 | 2 | 氧 |
| 28.79..161.55 | `EP <= 1580` | 1 | 氢（模糊带低幅度） |
| 28.79..161.55 | `EP > 1580` | 2 | 氧（模糊带高幅度） |
| 28.79..161.55 | `EP` 非有限 | 0 | 未定 |

阈值 `28.79`、`161.55`、`1580` 在 [conversion.py:165-176](conversion.py#L165-L176)。`Spec` 输出 `uint8`。

### 6.3 位置

`X1/Y1/X2/Y2` 由 `params['ns_2_xy']` 计算，见 [conversion.py:253-259](conversion.py#L253-L259) 和 [../constants.py:218-231](../constants.py#L218-L231)：

```text
start: X1 = KX * (TX_1s - TX_2s) + BX
       Y1 = KY * (TY_1s - TY_2s) + BY
end:   X2 = KX * (TX_1e - TX_2e) + BX
       Y2 = KY * (TY_1e - TY_2e) + BY
```

对探头 `p` 和方向 `d`，斜率/截距由 [../constants.py:148-155](../constants.py#L148-L155)：

```text
K = D / (R - L)
B = -D * (R + L) / (R - L) / 2
```

其中 `D_X = 145 mm`、`D_Y = 55 mm`，`L/R` 是该探头该方向的左右边界时间差，见 [../constants.py:179-202](../constants.py#L179-L202)。输出 `float32`。

### 6.4 能量与角度

| L2 列 | 公式 | 单位 | dtype | 定义 |
| --- | --- | --- | --- | --- |
| `Energy` | `0.5 * mass * beta^2 / 1000` | keV | `float32` | [conversion.py:204-219](conversion.py#L204-L219) |
| `PX` | `arctan2(X1 - X2, 50) * 180/pi` | deg | `float32` | [conversion.py:263](conversion.py#L263) |
| `PY` | `arctan2(Y1 - Y2, 50) * 180/pi` | deg | `float32` | [conversion.py:264](conversion.py#L264) |

`Energy` 的中间量：

```text
distance = sqrt((X1-X2)^2 + (Y1-Y2)^2 + 50^2)   [mm]
beta     = distance*1e-3 / (Tof*1e-9) / SPEED_C
mass     = MASS_O if Spec == 2 else MASS_H       [eV]
Energy   = 0.5 * mass * beta^2 / 1000            [keV]
```

只有 `Tof` 有限且 `>0`、`Spec` 为 1/2、`X1/X2/Y1/Y2` 全部有限时才计算；否则 NaN，见 [conversion.py:209-218](conversion.py#L209-L218)。

### 6.5 L2 中未赋值的列

`df2_type` 在 [structure.py:311-344](structure.py#L311-L344) 声明了 `GRP: uint8`，但 `l1_2_l2()` 的赋值列表 [conversion.py:242-264](conversion.py#L242-L264) 没有 `GRP`。因此只要有数据行，`GRP` 实际会变成 `float64` 全 NaN；空的 L2 才显示声明时的 `uint8`。检查 L2 时不要把它当作有效分组列。

### 6.6 L2S / L2G 选择使用哪些字段

L2S/L2G 不是新的字段转换，而是在 L2 上按 L1A 的 TDC 和 L2 的位置做布尔选择。`soft_std()`/`golden_std()` 在 [gena/filter.py:421-434](../filter.py#L421-L434)，具体的 `check()` 在 [gena/filter.py:230-275](../filter.py#L230-L275)。

| 产品 | 使用的字段/条件 | 参数 | 定义 |
| --- | --- | --- | --- |
| L2S `soft_std` | `T{dir}_1s < SUM*(1+0.2)`、`T{dir}_2s < SUM*(1+0.2)` | `X_SUM`/`Y_SUM` | [gena/filter.py:421-428](../filter.py#L421-L428) |
| L2S | `T{dir}_1s < T{dir}_1e`、`T{dir}_2s < T{dir}_2e` | 无 | [gena/filter.py:113-116](../filter.py#L113-L116) |
| L2S | `SUM/1.2 < T{dir}_1s + T{dir}_2s < SUM*1.2` | `X_SUM`/`Y_SUM` | [gena/filter.py:119-131](../filter.py#L119-L131) |
| L2S | `XL*1.2 < T{dir}_1s - T{dir}_2s < XR*1.2`，终点同理 | `XL/XR/YL/YR` | [gena/filter.py:133-151](../filter.py#L133-L151) |
| L2G | L2S 全部条件 | 同上 | [gena/filter.py:421-431](../filter.py#L421-L431) |
| L2G | `tofx/tofy` 在 `0.937..1.117` 且两者均 `>0` | `p_low/p_high` | [gena/filter.py:379-387](../filter.py#L379-L387) |
| L2G | `X1/Y1`、`X2/Y2` 在 `c_size=[150,60]` 的矩形内 | 150×60 mm | [gena/filter.py:316-329](../filter.py#L316-L329)、[decode.py:346-347](decode.py#L346-L347) |

L2S/L2G 的列与 L2 完全相同，只是行被上述掩码筛选；`decode` 中 L2S/L2G 分别由 `df1a_a.loc[func_soft...]`、`df1a_b.loc[func_gold...]` 重新进入 `l1_2_l2`，见 [decode.py:326-347](decode.py#L326-L347)。

## 7. 声明类型与实际类型差异汇总

pandas 给现有列赋值时，实际 dtype 不一定保持 `structure.py` 的声明；以下差异在当前代码中可以直接观察到。

| 对象/列 | `structure.py` 声明 | 当前实际类型 | 原因/位置 |
| --- | --- | --- | --- |
| PDF `HV_N`、`HV_P`、`HV_C`、`HV_G`、`HV_A`、`HV_B` | `float32` | `float64` | 整数原始码乘 Python `float` 常数，见 [packet.py:350-358](packet.py#L350-L358) |
| PDF `A_N12`、`A_N5`、`A_3_3`、`D_3_3`、`HV_C_I`、`HV_N_I` | `float32` | `float64` | 同上，见 [packet.py:360-372](packet.py#L360-L372) |
| PDF `Ant_Ang`、`Ant_Azi`、`Star_L`、`Star_NZ` | `float32` | `float64` | 整数乘 `1e-5`，未 `astype`，见 [packet.py:426-433](packet.py#L426-L433) |
| PDF `Shadow` | `float32` | `uint8` | `np.where(..., flex[22], 0)`，见 [packet.py:434](packet.py#L434) |
| PDF `Evt_Num_Tot/Eff/Aux` | `uint32` | `uint16` | `byteswap().view(...)`，见 [packet.py:383-388](packet.py#L383-L388) |
| PDF `PlatformGroup`、`PosePresent`、`ArgsAddress`、`ArgsValue` | 未声明 | `uint8`/`bool`/`uint8`/`uint32` | 新列未加入 `pdf_type`，见 [packet.py:396-397](packet.py#L396-L397)、[packet.py:440-441](packet.py#L440-L441) |
| DF0 `UUID` | `S36` | `str` | `np.repeat(uuids, 62)` 为 Python 字符串，见 [packet.py:524-525](packet.py#L524-L525) |
| DF0 `is_empty_slot` | 未声明 | `bool` | 新增列，见 [packet.py:534-535](packet.py#L534-L535) |
| L1A `HV_C`/`HV_G` | `float32` | `float64` | `interp1d` 输出，见 [conversion.py:114-115](conversion.py#L114-L115) |
| L1A `MJ2000_q1..q3` | `float32` | `float64` | `Rotation.as_quat()` 输出，见 [conversion.py:118-122](conversion.py#L118-L122) |
| L1A `MJ2000_Star_X/Y/Z` | `int32` | `float64` | `CubicSpline`/线性插值输出，见 [conversion.py:125-128](conversion.py#L125-L128) |
| L2 `GRP` | `uint8` | `float64` 全 NaN | 从未赋值，见 [structure.py:317](structure.py#L317) 和 [conversion.py:242-264](conversion.py#L242-L264) |
| L2 基础列 | `float32`/`uint8` 等 | 可能跟随输入 DataFrame 的 dtype | `df2[k] = df1[k]` 直接赋值，见 [conversion.py:242-243](conversion.py#L242-L243) |

## 8. `astype` / `view` / 显式 dtype 清单

### 8.1 `astype` 调用

| 文件:行 | 表达式 | 源类型 | 目标类型 | 说明 |
| --- | --- | --- | --- | --- |
| [packet.py:335](packet.py#L335) | `(seq2 >> 14).astype(...)` | 整数 | `uint8` | PDF `BM` |
| [packet.py:337](packet.py#L337) | `(seq2 & 0x3fff).astype(...)` | 整数 | `uint16` | PDF `BX` |
| [packet.py:380](packet.py#L380) | `Temps[:,0].astype(...)` | `float64` | `float32` | PDF `Temp_A` |
| [packet.py:381](packet.py#L381) | `Temps[:,1].astype(...)` | `float64` | `float32` | PDF `Temp_B` |
| [packet.py:382](packet.py#L382) | `Temps[:,2].astype(...)` | `float64` | `float32` | PDF `Temp_PCB` |
| [packet.py:396](packet.py#L396) | `group.astype(...)` | `uint8` | `uint8` | PDF `PlatformGroup` |
| [packet.py:403](packet.py#L403) | `np.where(...).astype(...)` | `float64` | `float32` | PDF `MJ2000_q1` |
| [packet.py:404](packet.py#L404) | 同上 | `float64` | `float32` | PDF `MJ2000_q2` |
| [packet.py:405](packet.py#L405) | 同上 | `float64` | `float32` | PDF `MJ2000_q3` |
| [packet.py:408](packet.py#L408) | `np.where(...).astype(...)` | 整数 | `int32` | PDF `MJ2000_Star_X` |
| [packet.py:409](packet.py#L409) | 同上 | 整数 | `int32` | PDF `MJ2000_Star_Y` |
| [packet.py:410](packet.py#L410) | 同上 | 整数 | `int32` | PDF `MJ2000_Star_Z` |
| [packet.py:440](packet.py#L440) | `para_addr.astype(...)` | `uint8` | `uint8` | PDF `ArgsAddress` |
| [packet.py:441](packet.py#L441) | `view_as(...).astype(...)` | `uint32` | `uint32` | PDF `ArgsValue` |
| [packet.py:448](packet.py#L448) | `(... & 0xFFF).astype(...)` | `uint32` | `int16` | PDF `Th_A_MCP` |
| [packet.py:449](packet.py#L449) | 同上 | `uint32` | `int16` | PDF `Th_B_MCP` |
| [packet.py:451](packet.py#L451) | 同上 | `uint32` | `int16` | PDF `Th_A_X` |
| [packet.py:452](packet.py#L452) | 同上 | `uint32` | `int16` | PDF `Th_A_Y` |
| [packet.py:454](packet.py#L454) | `signed_7bits(...).astype(...)` | `int8` | `int16` | PDF `Th_B_X1` |
| [packet.py:455](packet.py#L455) | 同上 | `int8` | `int16` | PDF `Th_B_X2` |
| [packet.py:456](packet.py#L456) | 同上 | `int8` | `int16` | PDF `Th_B_Y1` |
| [packet.py:457](packet.py#L457) | 同上 | `int8` | `int16` | PDF `Th_B_Y2` |
| [packet.py:529](packet.py#L529) | `np.repeat(...).astype(...)` | 整数 | `uint8` | DF0 `BM` |
| [packet.py:530](packet.py#L530) | `np.repeat(...).astype(...)` | 整数 | `uint16` | DF0 `BX` |
| [packet.py:531](packet.py#L531) | `np.repeat(...).astype(...)` | 浮点 | `float64` | DF0 `Timestamp` |
| [packet.py:549](packet.py#L549) | `(...).astype(...)` | `float64` | `float32` | DF0 `TimeInternal` |
| [packet.py:552](packet.py#L552) | `se_amp[...,1].astype(...)`、`se_amp[...,2].astype(...)` | `uint8` | `uint16` | 构造 `EP` 中间值 |
| [packet.py:573](packet.py#L573) | `ch01_tail[...,0].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_1s` |
| [packet.py:574](packet.py#L574) | `ch01_tail[...,1].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_1e` |
| [packet.py:575](packet.py#L575) | `ch03_tail[...,0].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_3s` |
| [packet.py:576](packet.py#L576) | `ch03_tail[...,1].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_3e` |
| [packet.py:577](packet.py#L577) | `ch05_tail[...,0].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_5s` |
| [packet.py:578](packet.py#L578) | `ch05_tail[...,1].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_5e` |
| [packet.py:579](packet.py#L579) | `ch07_tail[...,0].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_7s` |
| [packet.py:580](packet.py#L580) | `ch07_tail[...,1].flatten().astype(...)` | `uint16` | `uint16` | DF0 `ch_7e` |
| [conversion.py:16](conversion.py#L16) | `np.empty(len(tdc), dtype=...)` | 新建 | `float32` | L1 TDC 输出 |
| [conversion.py:19](conversion.py#L19) | `(tdc - 94) * np.float32(0.081)` | `float32` | `float32` | L1 TDC 比例 |
| [conversion.py:31](conversion.py#L31) | `np.full(..., dtype=...)` | 新建 | `float32` | L1 EP 输出 |
| [conversion.py:151](conversion.py#L151) | `np.full(..., dtype=...)` | 新建 | `float32` | 电子修正输出 |
| [conversion.py:168](conversion.py#L168) | `np.zeros(..., dtype=...)` | 新建 | `uint8` | L2 `Spec` |
| [conversion.py:198](conversion.py#L198) | `np.full(..., dtype=...)` | 新建 | `float32` | L2 `Tof` |
| [conversion.py:210](conversion.py#L210) | `np.full(..., dtype=...)` | 新建 | `float32` | L2 `Energy` |
| [conversion.py:263](conversion.py#L263) | `arctan2(...).astype(...)` | `float64` | `float32` | L2 `PX` |
| [conversion.py:264](conversion.py#L264) | `arctan2(...).astype(...)` | `float64` | `float32` | L2 `PY` |

### 8.2 `view` 调用

| 文件:行 | 表达式 | 作用 |
| --- | --- | --- |
| [packet.py:319-320](packet.py#L319-L320) | `np.ascontiguousarray(bytes_array).view(dtype)` | 4 字节大端转 `>i4`/`>u4`/`int8` |
| [packet.py:386-388](packet.py#L386-L388) | `byteswap().view(dtype.newbyteorder())` | `Evt_Num_*` 大小端兼容 |
| [packet.py:441](packet.py#L441) | `view_as(paras, '>u4')` | 参数表值转 `uint32` |
| [packet.py:448-452](packet.py#L448-L452) | `view_as(paras[..., ::-1], np.uint32)` | 参数值反序后取低 12 位 |
| [packet.py:458](packet.py#L458) | `view_as(paras[...][:, :1], np.int8)` | `Th_B_G` |
| [packet.py:482](packet.py#L482) | `p.view(np.uint8)` | 整包字节和校验 |
| [packet.py:534](packet.py#L534) | `np.ascontiguousarray(p_dict_list).view(np.uint8)` | 判断 32 字节槽是否全零 |

## 9. 标定参数与物理常数

### 9.1 位置标定

`init_param()` 在 [../constants.py:167-233](../constants.py#L167-L233) 定义：

| 参数 | 值 | 含义 | 定义 |
| --- | ---: | --- | --- |
| `DX` | 145 | X 方向 MCP 长度，mm | [../constants.py:179](../constants.py#L179) |
| `DY` | 55 | Y 方向 MCP 长度，mm | [../constants.py:180](../constants.py#L180) |
| `A_XL` | -54.91726 | A 探头 X 左边界时间差，ns | [../constants.py:182](../constants.py#L182) |
| `A_XR` | 50.44742 | A 探头 X 右边界时间差，ns | [../constants.py:183](../constants.py#L183) |
| `A_YL` | -16.68598 | A 探头 Y 左边界时间差，ns | [../constants.py:184](../constants.py#L184) |
| `A_YR` | 41.62719 | A 探头 Y 右边界时间差，ns | [../constants.py:185](../constants.py#L185) |
| `B_XL` | -51.18360 | B 探头 X 左边界时间差，ns | [../constants.py:188](../constants.py#L188) |
| `B_XR` | 53.40479 | B 探头 X 右边界时间差，ns | [../constants.py:189](../constants.py#L189) |
| `B_YL`（`snid=1`） | -16.68598 | B 探头在轨 Y 左边界 | [../constants.py:191-193](../constants.py#L191-L193) |
| `B_YR`（`snid=1`） | 41.62719 | B 探头在轨 Y 右边界 | [../constants.py:191-193](../constants.py#L191-L193) |
| `B_YL`（`snid=2`） | -34.71570 | B 探头地面 Y 左边界 | [../constants.py:194-196](../constants.py#L194-L196) |
| `B_YR`（`snid=2`） | 23.59747 | B 探头地面 Y 右边界 | [../constants.py:194-196](../constants.py#L194-L196) |
| `A_X_SUM` | 50.0 | A 探头 X 起止和辅助值 | [../constants.py:206](../constants.py#L206) |
| `A_Y_SUM` | 38.5 | A 探头 Y 起止和辅助值 | [../constants.py:207](../constants.py#L207) |
| `B_X_SUM` | 55.0 | B 探头 X 起止和辅助值 | [../constants.py:210](../constants.py#L210) |
| `B_Y_SUM` | 43.0 | B 探头 Y 起止和辅助值 | [../constants.py:211](../constants.py#L211) |
| `KX/KY/BX/BY` | 由 `K_B` 计算 | 时间差→mm 的斜率/截距 | [../constants.py:148-155](../constants.py#L148-L155) |
| `ns_2_xy` | 闭包 | `K * (T1 - T2) + B` | [../constants.py:218-231](../constants.py#L218-L231) |

`init_param` 最后用 `cast_dict_to_float32` 把数值参数转成 `np.float32`，见 [../constants.py:157-165](../constants.py#L157-L165) 和 [../constants.py:216](../constants.py#L216)。

### 9.2 物理常数

| 常数 | 值 | 单位 | 定义 |
| --- | ---: | --- | --- |
| `SPEED_C` | 299792458 | m/s | [../../common/constants.py:29](../../common/constants.py#L29) |
| `MASS_E` | 510998.918 | eV | [../../common/constants.py:33](../../common/constants.py#L33) |
| `MASS_H` | `MASS_E * 1837.15267343` | eV | [../../common/constants.py:37](../../common/constants.py#L37) |
| `MASS_O` | `MASS_H * 15.8695` | eV | [../../common/constants.py:49](../../common/constants.py#L49) |
| `RE` | 6378100.0 | m | [../../common/constants.py:56](../../common/constants.py#L56) |
| `RM` | `0.27 * RE` | m | [../../common/constants.py:60](../../common/constants.py#L60) |

### 9.3 硬编码的转换/筛选常数

| 常数 | 值 | 用途 | 定义 |
| --- | ---: | --- | --- |
| 包时间基准 | 1609430400 | `Timestamp` 的 UTC 秒基准 | [packet.py:339-346](packet.py#L339-L346) |
| `50000/32768` | 1.52587890625 | `HV_N`/`HV_P` V/code | [packet.py:350](packet.py#L350) |
| `10000/32768` | 0.30517578125 | `HV_C/G/A/B` V/code | [packet.py:354](packet.py#L354) |
| `50/128` | 0.390625 | `A_N12` 满量程系数 | [packet.py:360](packet.py#L360) |
| `20/128` | 0.15625 | `A_N5` 满量程系数 | [packet.py:364](packet.py#L364) |
| `10/128` | 0.078125 | `A_3_3`/`D_3_3`/电流监视系数 | [packet.py:367](packet.py#L367) |
| 温度系数 | 0.078125 | `temps` 原始码→电压 | [packet.py:376](packet.py#L376) |
| `q_dl` | `2**-30` | 四元数定点当量 | [packet.py:402](packet.py#L402) |
| `a_dl`/`m_dl` | `1e-5` | 指向角定点当量，deg | [packet.py:427-431](packet.py#L427-L431) |
| TDC 门限 | 94 | `code <= 94` 记为无测量 | [conversion.py:14](conversion.py#L14) |
| TDC 比例 | 0.081 | ns/code | [conversion.py:19](conversion.py#L19) |
| EP 上限码 | 2044 | `(2044 - code)` 的中心 | [conversion.py:33](conversion.py#L33) |
| EP 比例 | 100/123 | 幅度/code | [conversion.py:33](conversion.py#L33) |
| 电子加速区 | `3.1e-3` m | `f_tof_ep` | [conversion.py:146](conversion.py#L146) |
| 电子飞行区 | `47e-3` m | `f_tof_ep` | [conversion.py:147](conversion.py#L147) |
| 主 TOF 绝对窗 | 5 ns | `ENA_flytime` | [conversion.py:194](conversion.py#L194) |
| X/Y TOF 比值窗 | 0.9..1.1 | `ENA_flytime` | [conversion.py:196](conversion.py#L196) |
| H/O 分界 | 28.79 / 161.55 ns | `ENA_spec` | [conversion.py:165-166](conversion.py#L165-L166) |
| EP 种类阈值 | 1580 | `ENA_spec` 模糊带 | [conversion.py:175-176](conversion.py#L175-L176) |
| 能量几何常数 | 50 mm | `ENA_energy` 与 `PX/PY` | [conversion.py:213](conversion.py#L213)、[conversion.py:263-264](conversion.py#L263-L264) |
| `TimeInternal` 计数当量 | 0.00002 s | `time_us`→秒 | [packet.py:549](packet.py#L549) |

## 10. 当前检查重点

以下是字段层面最值得先核对的点，均可在上表中定位：

1. `monitors` 声明为 `int8`，但缩放系数按无符号满量程设计；最高位为 1 时 `A_N12`、`A_N5`、`A_3_3`、`D_3_3`、`HV_C_I`、`HV_N_I` 会变负，见 [structure.py:12](structure.py#L12) 与 [packet.py:360-372](packet.py#L360-L372)。
2. `pdf_type` 没有 `PlatformGroup`、`PosePresent`、`ArgsAddress`、`ArgsValue`，且 `HV_*`、`Ant_*`、`Shadow`、`Evt_Num_*` 的声明类型与实际类型不一致，见 [structure.py:56-109](structure.py#L56-L109) 与第 7 节。
3. `df1a_type` 声明 `MJ2000_Star_*` 为 `int32`、四元数为 `float32`，实际插值给的是 `float64`，见 [structure.py:184-215](structure.py#L184-L215) 与 [conversion.py:114-128](conversion.py#L114-L128)。
4. `df2_type` 的 `GRP` 没有赋值，实际为全 NaN `float64`，见 [structure.py:317](structure.py#L317) 与 [conversion.py:242-264](conversion.py#L242-L264)。
5. `TimeInternal` 用 20 µs 当量并强制 `float32`，大数值时 20 µs 分辨率会被浮点精度吞掉；同时它没有加到包级 `Timestamp` 上，见 [packet.py:549](packet.py#L549) 与 [packet.py:531](packet.py#L531)。
6. `signed_7bits` 采用符号-幅值，不是补码；`Th_B_X1/X2/Y1/Y2` 的负值映射需要按协议确认，见 [packet.py:304-314](packet.py#L304-L314)。
7. `l1_2_l1a` 的八通道掩码是 `!= 0`，NaN 会被当作非零保留；正常 L1 路径不产生 NaN TDC，但外部传入 L1 时需要额外校验，见 [conversion.py:105](conversion.py#L105)。
8. `is_empty_slot` 只表示 32 字节全零；全零槽的 `prob_flag=0`，会被 `l0_2_l1(...,'a')` 选中，当前 `l0_2_l1()` 没有读取 `is_empty_slot`，见 [packet.py:534-535](packet.py#L534-L535) 与 [conversion.py:50-57](conversion.py#L50-L57)。
9. 当前 `gena/raw/__init__.py` 引用不存在的 `pipeline.py`，`decode.py` 的 repair 分支引用不存在的 `rep.l1a_2_l1b`；阅读字段时以现存 `packet.py`/`conversion.py` 为准，运行前需先修复这两处导入。

## 11. 附：L1B / repair 字段（默认不进入有效测量值）

`repair.py` 保留了历史修复路径，`df1b_type` 在 [structure.py:243-301](structure.py#L243-L301) 声明。默认 `GENA_USE_REPAIR=0`，且 `repair_l1a()` 在 [repair.py:487-502](repair.py#L487-L502) 直接 `raise NotImplementedError`，所以下列字段不属于当前“有效测量值”链路。

| L1B 字段 | 声明 dtype | 含义 | 定义 |
| --- | --- | --- | --- |
| `UUID`、`PID`、`BM`、`BX` | `S36`/`uint8`/`uint8`/`uint16` | 基础标识 | [structure.py:244-247](structure.py#L244-L247) |
| `HV_C`、`HV_G` | `float32` | 插值后的高压 | [structure.py:249-250](structure.py#L249-L250) |
| `MJ2000_q1..q3` | `float32` | 姿态 | [structure.py:251-253](structure.py#L251-L253) |
| `MJ2000_Star_X/Y/Z` | `int32` | 月心 J2000 位置，m | [structure.py:254-256](structure.py#L254-L256) |
| `Timestamp`、`TimeInternal`、`EvtLevel`、`EP` | `float64`/`float32`/`uint8`/`float32` | 时间与幅度 | [structure.py:259-262](structure.py#L259-L262) |
| `Tof_0`、`Tof_1` | `float32` | 主 TOF 槽 | [structure.py:263-264](structure.py#L263-L264) |
| `TX_*`、`TY_*` | `float32` | 修复后的 TDC 时间 | [structure.py:265-272](structure.py#L265-L272) |
| `Type` | `uint8` | 修复等级（0/1..4/10/20/30/31/40/41/255） | [structure.py:291](structure.py#L291) |
| `Process_Mask` | `uint32` | 修复步骤位掩码 | [structure.py:292](structure.py#L292) |
| `Tofs`、`Tofe`、`Tofx0`、`Tofy0`、`Tofx1`、`Tofy1` | `float32` | 起止/轴 TOF 中间量 | [structure.py:277-282](structure.py#L277-L282) |
| `tmpX1`、`tmpX2`、`tmpY1`、`tmpY2` | `float32` | 位置中间量 | [structure.py:283-286](structure.py#L283-L286) |
| `pop_TX_1e`、`pop_TX_2e`、`pop_TY_1e`、`pop_TY_2e` | `float32` | 被替换下来的候选截止信号 | [structure.py:287-290](structure.py#L287-L290) |
| `Xs_valid`、`Ys_valid`、`Xe_valid`、`Ye_valid` | `uint8` | 起止时间差是否在边界内 | [structure.py:295-298](structure.py#L295-L298) |
| `start_valid`、`end_valid` | `uint8` | 起止位置是否在板内 | [structure.py:299-300](structure.py#L299-L300) |

`Type` 的语义注释在 [structure.py:217-241](structure.py#L217-L241)。这些中间量会出现在 `repair.py` 的 `_repair_start`、`_repair_end`、`_repair_method_1/2` 中；当前生产路径不写 L1B，也不把推定值写入黄金样本。
