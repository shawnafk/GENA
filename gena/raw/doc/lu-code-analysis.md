# LU Fortran 代码深度分析文档

## 1. 变量定义与物理意义

### 1.1 输入数据变量（从 TDC 读取）

| 变量名 | 类型 | 物理意义 | 备注 |
|--------|------|----------|------|
| `date_str` | CHARACTER*10 | 日期字符串 | 格式: YYYY-MM-DD |
| `time_str` | CHARACTER*14 | 时间字符串 | 格式: HH:MM:SS |
| `TOF0`, `TOF1` | REAL | 飞行时间原始值 | ns，TOF0通常更小 |
| `X1s`, `X1e` | REAL | X方向左通道时间 | s=start(起始), e=end(终止) |
| `X2s`, `X2e` | REAL | X方向右通道时间 | 单位: ns |
| `Y1s`, `Y1e` | REAL | Y方向下通道时间 | 探测器读出顺序 |
| `Y2s`, `Y2e` | REAL | Y方向上通道时间 | 单位: ns |
| `AMP` | REAL | 信号幅度 | 用于粒子鉴别 |

### 1.2 校准参数变量（从 ITERANCE.NML 读取）

#### Probe A 参数:
| 参数 | 数值 | 物理意义 |
|------|------|----------|
| `AXL` | -54.917 ns | X方向左边界时间 |
| `AXR` | 50.447 ns | X方向右边界时间 |
| `AYD` | -16.686 ns | Y方向下边界时间 |
| `AYU` | 41.627 ns | Y方向上边界时间 |
| `AETof` | 2.537 ns | Probe A 飞行时间偏移 |

#### Probe B 参数:
| 参数 | 数值 | 物理意义 |
|------|------|----------|
| `BXL` | -53.405 ns | X方向左边界时间 |
| `BXR` | 51.184 ns | X方向右边界时间 |
| `BYD` | -16.686 ns | Y方向下边界时间 |
| `BYU` | 41.627 ns | Y方向上边界时间 |
| `BETof` | 2.441 ns | Probe B 飞行时间偏移 |

### 1.3 计算派生变量

| 变量 | 计算公式 | 物理意义 | 单位 |
|------|----------|----------|------|
| `CofX` | `145 / (XR - XL)` | X方向转换系数 | mm/ns |
| `CofY` | `55 / (YU - YD)` | Y方向转换系数 | mm/ns |
| `Z0X` | `(XR + XL) / 2` | X方向时间零点 | ns |
| `Z0Y` | `(YU + YD) / 2` | Y方向时间零点 | ns |
| `CX` | `MAX(ABS(XR), ABS(XL))` | X方向最大时间偏差 | ns |
| `CY` | `MAX(ABS(YU), ABS(YD))` | Y方向最大时间偏差 | ns |

### 1.4 中间计算变量

| 变量 | 物理意义 | 说明 |
|------|----------|------|
| `Tofx0` | X方向起始时间 | 两个X通道的平均值修正 |
| `Tofy0` | Y方向起始时间 | 两个Y通道的平均值修正 |
| `Tofs` | 最终起始时间 | MAX(Tofx0, Tofy0) |
| `Tofe` | 最终终止时间 | MIN(Tofy1, Tofx1) 或其他逻辑 |
| `Xs`, `Ys` | 起始位置坐标 | 单位: mm |
| `Xe`, `Ye` | 终止位置坐标 | 单位: mm |
| `KX`, `KY` | 起始信号有效性标志 | 0=有效, 1=缺失/越界 |
| `KKX`, `KKY` | 终止信号有效性标志 | 0=有效, 1=缺失/越界 |

---

## 2. 数据修复方法详解

### 2.1 时间数据预处理（Line 150-157）

```fortran
IF (Tof0.GT.2.0.AND.Tof1.LT.Tof0) THEN
    Tof1=Tof0
    Tof0=0.0
END IF
```
**目的**: 确保 TOF 数据顺序正确，将较大的值放在 Tof1。

```fortran
IF (X1e.NE.0.0.AND.X1e.LT.X1s) X1e=X1s+663.47
IF (X2e.NE.0.0.AND.X2e.LT.X2s) X2e=X2s+663.47
IF (Y1s.NE.0.0.AND.Y1e.LT.Y1s) Y1e=Y1s+663.47
IF (Y2e.NE.0.0.AND.Y2s.LT.Y1s) Y2e=Y1s+663.47
```
**目的**: 处理 TDC 计数器翻转（rollover）问题。

**物理背景**: TDC 是一个循环计数器，当计数值从最大值回到 0 时，如果终止时间小于起始时间，需要加上一个周期（663.47 ns）来修正。

**⚠️ BUG 警告**: 最后一行 `Y2e=Y1s+663.47` 应该是 `Y2e=Y2s+663.47`（用了 Y1s 而不是 Y2s）。

### 2.2 零信号检测（Line 159-164）

```fortran
IF((X1e.EQ.0.0.AND.X2E.EQ.0.0).OR.(Y1e.EQ.0.0.AND.Y2E.EQ.0.0)) THEN
    WRITE(6,111)...  ! 写入 VOID 文件
    NVO=NVO+1
    GOTO 10
END IF
```
**逻辑**: 如果 X 方向或 Y 方向的两个通道都为零，判定为无效信号，丢弃到 VOID 文件。

---

## 3. 坐标计算方法

### 3.1 标准坐标计算（双信号有效）

当两个通道都有效时（在边界范围内）：

```fortran
Tofx0 = (X2s + X1s - Cx) / 2
Xs = (X2s - X1s - Z0X) * CofX

Tofy0 = (Y2s + Y1s - Cy) / 2
Ys = (Y2s - Y1s - Z0Y) * CofY
```

**物理意义**:
- `X2s - X1s`: 时间差，正比于 X 位置
- `Z0X`: 零点校正，使中心位置对应 X=0
- `CofX`: 将时间差转换为空间距离（mm）

### 3.2 单信号缺失时的坐标重建

当 Y 信号缺失但 X 信号有效时（Line 192-199）：

```fortran
IF(KY.EQ.1.AND.KX.EQ.0) Tofs=Tofx0  ! 使用 X 时间作为参考

IF(ABS(Y2S-TofS).LT.ABS(Y1s-TofS)) THEN
    ! Y2s 更接近 Tofs，使用 Y2s 重建
    Ys=(2*(TofS-Y2s)-CY-Z0Y)*CofY    ! Line 194
    if((Y1s-Tofs).GT.CY) Y1e=Y1s      ! 检查另一通道是否翻转
ELSE
    ! Y1s 更接近 Tofs，使用 Y1s 重建
    Ys=(2*(Tofs-Y1s)+CY-Z0Y)*CofY    ! Line 197
    if((Y2s-Tofs).GT.CY) Y2e=Y2s
END IF
```

**Line 194 详解**:
```fortran
Ys=(2*(TofS-Y2s)-CY-Z0Y)*CofY
```
推导过程：
- 假设 Y2s 是 "正确" 的通道（更接近 Tofs）
- 重建公式: `Ys = (Δt_Y - Z0Y) * CofY`
- 其中 `Δt_Y = 2*(TofS - Y2s) - CY`
- 这相当于假设 `TofS ≈ (Y1s + Y2s - CY) / 2`，反解出 Y1s 再计算差值

当 X 信号缺失但 Y 信号有效时（Line 202-209）：

```fortran
IF(KY.EQ.0.AND.KX.EQ.1) Tofs=Tofy0  ! 使用 Y 时间作为参考

IF(ABS(X2S-TofS).LT.ABS(X1s-TofS)) THEN
    Xs=(2*(X2s-TofS)-Cx-Z0X)*CofX    ! Line 204
    IF((X1s-Tofs).GT.CX) X1e=X1s
ELSE
    Xs=(2*(Tofs-X1s)+CX-Z0X)*Cofy    ! Line 207 ⚠️ BUG!
    IF((X2s-Tofs).GT.CX) X2e=X2s
END IF
```

**Line 207 详解**:
```fortran
Xs=(2*(Tofs-X1s)+CX-Z0X)*Cofy    ! 注意: 这里用的是 Cofy 而不是 CofX!
```

**⚠️ 关键 BUG**: 这里使用了 `Cofy` 而不是 `CofX`，这会导致 X 坐标被错误地缩放！正确的代码应该是:
```fortran
Xs=(2*(Tofs-X1s)+CX-Z0X)*CofX    ! 正确的写法
```

### 3.3 终止位置坐标计算（Line 216-229）

与起始位置计算类似：

```fortran
20  KKX=0
    KKY=0
    IF((X2e-X1e).GT.XL.AND.(X2e-X1e).LT.XR) THEN 
        Tofx1=(X2e+X1e-Cx)/2
        Xe=(X2e-X1e-Z0X)*CofX        ! 标准计算
    else 
        KKX=1
    end if
    IF((Y2e-Y1e).GT.YD.AND.(Y2e-Y1e).LT.YU) THEN 
        Tofy1=(Y2e+Y1e-Cy)/2
        Ye=(Y2e-Y1e-Z0Y)*CofY        ! 标准计算
    else
        KKY=1
    end if
```

### 3.4 终止位置单信号重建（Line 247-259）

```fortran
IF(KY.EQ.1.AND.KX.EQ.0) Tofe=Tofx1
IF(ABS(Y2e-Tofe).LT.ABS(Y1e-Tofe)) THEN
    Ye=(2*(Y2e-Tofe)-Cy-Z0Y)*CofY    ! Line 249
ELSE
    Ye=(2*(Tofe-Y1e)+Cy-Z0Y)*CofY    ! Line 251
END IF

IF(KY.EQ.0.AND.KX.EQ.1) Tofe=Tofy1
IF(ABS(X2e-Tofe).LT.ABS(X1e-Tofe)) THEN
    Xe=(2*(X2e-Tofe)-Cx-Z0X)*CofX    ! Line 256
ELSE
    Xe=(2*(Tofe-X1e)+Cx-Z0X)*Cofy    ! Line 258 ⚠️ BUG!
END IF
```

**Line 249 详解**:
```fortran
Ye=(2*(Y2e-Tofe)-Cy-Z0Y)*CofY
```
使用 Y2e 通道重建 Y 坐标，假设 Y2e 更接近 Tofe。

**Line 256 详解**:
```fortran
Xe=(2*(X2e-Tofe)-Cx-Z0X)*CofX
```
这是正确的 X 坐标计算，使用 X2e 通道。

**⚠️ Line 258 BUG**:
```fortran
Xe=(2*(Tofe-X1e)+Cx-Z0X)*Cofy    ! 又是 Cofy 而不是 CofX!
```

---

## 4. 坐标计算公式总结

### 4.1 标准公式（双信号有效）

```
Xs = (Δt_Xs - Z0X) * CofX
Ys = (Δt_Ys - Z0Y) * CofY
Xe = (Δt_Xe - Z0X) * CofX
Ye = (Δt_Ye - Z0Y) * CofY

其中:
Δt_Xs = X2s - X1s
Δt_Ys = Y2s - Y1s
Δt_Xe = X2e - X1e
Δt_Ye = Y2e - Y1e
```

### 4.2 单信号重建公式

当使用参考时间 `T`（从另一方向获得）重建坐标时：

**对于 Y 坐标**（Line 194 和 197）:
```
如果使用 Y2（Y2 更接近 T）:
    Ys = (2*(T - Y2) - CY - Z0Y) * CofY

如果使用 Y1（Y1 更接近 T）:
    Ys = (2*(T - Y1) + CY - Z0Y) * CofY
```

**对于 X 坐标**（Line 204 和 207）:
```
如果使用 X2（X2 更接近 T）:
    Xs = (2*(X2 - T) - CX - Z0X) * CofX     ! Line 204，正确

如果使用 X1（X1 更接近 T）:
    Xs = (2*(T - X1) + CX - Z0X) * CofX     ! Line 207，应该是 CofX 不是 Cofy
                                            ! 当前代码错误地使用了 Cofy
```

### 4.3 物理意义解释

单信号重建的物理原理：

假设我们只有 Y2 的时间值，但知道参考时间 `Tofs`（来自 X 方向）。因为:
```
Tofs ≈ (Y1s + Y2s - CY) / 2
```

我们可以反推出 Y1s:
```
Y1s ≈ 2*Tofs - Y2s + CY
```

然后计算时间差:
```
Δt_Y = Y2s - Y1s = Y2s - (2*Tofs - Y2s + CY) = 2*Y2s - 2*Tofs - CY = 2*(Y2s - Tofs) - CY
```

由于 Y2s 可能大于或小于 Tofs，需要调整符号:
```
当 Y2s < Tofs 时: Δt_Y = 2*(Tofs - Y2s) - CY
当 Y1s < Tofs 时: Δt_Y = 2*(Tofs - Y1s) + CY
```

这就是 Line 194 和 197 的公式来源。

---

## 5. 后续处理与粒子鉴别

### 5.1 时间差计算

```fortran
DTof = Tofe - Tofs + ETof
```
其中 `ETof` 是探针特定的偏移量（Probe A: 2.537 ns, Probe B: 2.441 ns）。

### 5.2 有效范围检查

| 检查项 | 条件 | 失败处理 |
|--------|------|----------|
| TOF 范围 | `6 ≤ DTof ≤ 187` ns | 写入 INITIAL |
| X 位移 | `ABS(Xe - Xs) ≤ 45` mm | 写入 INITIAL |
| Y 位移 | `ABS(Ye - Ys) ≤ 15` mm | 写入 INITIAL |
| 终止 X | `ABS(Xe) ≤ 75` mm | 写入 INITIAL |
| 终止 Y | `ABS(Ye) ≤ 27.5` mm | 写入 INITIAL |

### 5.3 粒子鉴别

```fortran
V = DL / DTof * 1000.0  ! 速度 (mm/ns)

IF(DTof < 40.0) THEN
    NPar = 1          ! 氢 (H)
    Eng = HE * V^2    ! HE = 5.2247
ELSE IF(DTof > 66.3) THEN
    NPar = 0          ! 氧 (O)
    Eng = OE * V^2    ! OE = 82.927
ELSE
    ! 40 ≤ DTof ≤ 66.3，使用幅度鉴别
    IF(AMP > 1580) THEN
        NPar = 0      ! 氧（高幅度）
        Eng = OE * V^2
    ELSE
        NPar = 1      ! 氢（低幅度）
        Eng = HE * V^2
    END IF
END IF
```

### 5.4 角度计算

```fortran
Azi = ATAN(DX / 50.0) * 180.0 / PI        ! 方位角
Elv = ATAN(DY / 50.0) * 180.0 / PI ± 5.0  ! 仰角 (+5 for A, -5 for B)
```

---

## 6. BUG 汇总与建议

### 6.1 已确认的 BUG

| 位置 | 问题 | 影响 |
|------|------|------|
| Line 157 | `Y2e=Y1s+663.47` 应为 `Y2e=Y2s+663.47` | Y2 rollover 修正错误 |
| Line 207 | `Xs=...*Cofy` 应为 `Xs=...*CofX` | X 坐标缩放错误 |
| Line 258 | `Xe=...*Cofy` 应为 `Xe=...*CofX` | X 坐标缩放错误 |

### 6.2 影响分析

**CofX vs Cofy BUG**:
- `CofX = 145 / (XR - XL) ≈ 145 / 105 = 1.38 mm/ns`
- `CofY = 55 / (YU - YD) ≈ 55 / 58 = 0.95 mm/ns`
- 错误使用 `Cofy` 会导致 X 坐标被低估约 31%

### 6.3 修复建议

1. **立即修复**: 将 Line 207 和 Line 258 的 `Cofy` 改为 `CofX`
2. **修复 Line 157**: 将 `Y1s` 改为 `Y2s`
3. **添加注释**: 解释单信号重建的物理原理
4. **使用参数**: 将硬编码值（145, 55, 663.47 等）改为命名参数

---

## 7. 关键代码行详细注释

### Line 194: Y 坐标单信号重建（使用 Y2s）
```fortran
Ys=(2*(TofS-Y2s)-CY-Z0Y)*CofY
```
- **场景**: Y1s 缺失或越界，但 Y2s 有效
- **逻辑**: 使用 X 方向得到的 Tofs 作为参考，重建 Y 坐标
- **公式推导**: 假设 `Tofs ≈ (Y1s + Y2s - CY)/2`，反推 Y1s 后计算位置

### Line 207: X 坐标单信号重建（使用 X1s）⚠️ BUG
```fortran
Xs=(2*(Tofs-X1s)+CX-Z0X)*Cofy    ! BUG: 应该是 CofX
```
- **场景**: X2s 缺失或越界，但 X1s 有效
- **逻辑**: 使用 Y 方向得到的 Tofs 作为参考，重建 X 坐标
- **问题**: 使用了 Y 方向的转换系数 `Cofy` 而不是 X 方向的 `CofX`

### Line 249: Y 坐标终止位置重建（使用 Y2e）
```fortran
Ye=(2*(Y2e-Tofe)-Cy-Z0Y)*CofY
```
- **场景**: Y1e 缺失，使用 Y2e 重建终止 Y 坐标
- **注意**: 这里用的是 `Y2e - Tofe`（与起始位置的 `Tofs - Y2s` 符号相反）
- **原因**: 重建公式的对称性

### Line 256: X 坐标终止位置重建（使用 X2e）
```fortran
Xe=(2*(X2e-Tofe)-Cx-Z0X)*CofX
```
- **场景**: X1e 缺失，使用 X2e 重建终止 X 坐标
- **注意**: 这是**正确**的写法，使用了 `CofX`
- **对比**: 与 Line 258 形成对比，后者有 BUG

---

## 附录：变量命名对照表

| 代码变量 | 物理含义 | 维度 |
|----------|----------|------|
| X1, Y1 | X左/Y下通道 | - |
| X2, Y2 | X右/Y上通道 | - |
| s (suffix) | 起始 (start) | - |
| e (suffix) | 终止 (end) | - |
| AXL/AXR | Probe A X 左/右边界 | ns |
| AYU/AYD | Probe A Y 上/下边界 | ns |
| BXL/BXR | Probe B X 左/右边界 | ns |
| BYU/BYD | Probe B Y 上/下边界 | ns |
| CofX/CofY | 转换系数 | mm/ns |
| Z0X/Z0Y | 时间零点 | ns |
| CX/CY | 最大时间偏差 | ns |
