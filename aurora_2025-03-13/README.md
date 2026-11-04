# 2025-03-13 极光图像（THEMIS ASI）

**结论先说**：这一天的极光图像能从 THEMIS 全天空相机（ASI）公开档案里拿到，已经下载并整理成下面这套 PNG。DMSP 这边拿不到这一天的图：SSUSI 的公开档案在 2025-02-18 之后就没有了（JHU/APL 于 2025-02-14 停止支持），CDAWeb 上的 DMSP 粒子数据停得更早，所以 2025-03-13 的 DMSP 极光图在公开渠道不存在（详见下文「DMSP 部分」）。

## 1. 当晚空间天气（用来判断该看哪几个小时）

OMNI 每小时值（CDAWeb HAPI，`OMNI2_H0_MRG1HR`）：

| 时间 (UT) | Kp | Dst (nT) | AE (nT) |
|---|---|---|---|
| 03-13 03:30 | 4.7 | -46 | 580 |
| 03-13 05:30 | 4.7 | -52 | 316 |
| 03-13 07:30 | 5.0 | -54 | 609 |
| 03-13 09:30 | 4.3 | -66 | 546 |
| 03-13 11:30 | 4.3 | -55 | 512 |
| 03-13 12:30 | 5.0 | -53 | 888 |

Kp 4–5、Dst 最低约 −66 nT，属于中等扰动，极光椭圆扩展到 55°–60° 磁纬附近。北美的夜间对应 **2025-03-13 02:00–13:00 UT**，这一整套图就取自这个窗口。

有两个对成像不利的条件，图上都能直接看到：

- **满月**（2025-03-14 满月），月光把天空背景抬高到 3000–13000 counts，画面里那个亮斑就是月亮；
- **云**：很多站点当晚是多云/有霾，云被月光照亮后比极光还亮。相对最干净的是 **GAKO（Gakona, 62.4°N）04–07 UT** 和 **INUV（Inuvik, 68.4°N）05–13 UT**。

## 2. 目录内容

```
frames/<site>/<site>_20250313_HHMMSSut.png   单帧全天空图（12 站 × 10–12 个时刻，共 135 张）
mosaics/<site>_2025-03-13_sequence.png       每站当晚的时间序列拼图（4 列）
overview/stations_2025-03-13_0800ut.png      12 个站点同一分钟（08:00 UT）的横向对比
overview/keograms_6stations_2025-03-13.png    6 个站 02–13 UT 的 keogram 叠在一起（最能看出极光活动）
overview/gako_2025-03-13_0530-0533ut.gif      GAKO 05:30–05:33 UT 的 3 秒分辨率动画（60 帧）
keograms/<site>_2025-03-13_keogram.png       每站 02–13 UT 的整夜 keogram
quicklook/<day>_<hh>_<site>_<kind>.jpg        档案自带的每小时快视图（keogram / montage / average）
frames_index.csv                             每张图的站点、经纬度、UTC 时刻、源文件名
_raw/                                        原始存档（每分钟一个 .pgm.gz，20 帧 @3 s）
code/                                        下载与绘图脚本
```

`keograms/` 里的横轴是 UT，纵轴是沿天顶方向切的那一条天空；斜条纹是极光射线扫过，横条纹多是云和月光。

## 3. 数据来源

THEMIS ASI 原始图像由 UCalgary 公开档案分发（Berkeley 的 `themis.ssl.berkeley.edu` 也在，但当天目录是 `thg_l1_asf` 每小时 ~100 MB 的 CDF，取单帧不划算）：

```
https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream0/
    <YYYY>/<MM>/<DD>/<site>_themisNN/ut<HH>/<YYYYMMDD>_<HHMM>_<site>_themisNN_full.pgm.gz
https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream2/
    <YYYY>/<MM>/<DD>/<site>_themisNN/ut<HH>/<YYYYMMDD>_<HH>_<site>_themisNN_full-{keogram,montage,average}.pgm.jpg
```

每个 `_full.pgm.gz` 是一分钟的数据：20 帧 256×256、16 bit，3 秒曝光；文件头里带站点经纬度、光学投影（`a1=1.5, x0=382, y0=150`）、曝光时刻。注意只有一帧的头里带完整站点信息，后面 19 帧的头只有时间。

## 4. 显示方法

为了 135 张图彼此可比，统一用同一套拉伸（`code/render_asi.py`）：

1. 3×3 中值滤波去掉星点和宇宙线；
2. 天空本底 = 视场圆内（r = 0.97 R）像素中值；
3. 线性拉伸到圆内第 99.7 百分位，gamma 0.7；
4. 视场圆外全黑。

这是**展示用**拉伸，不是辐射定标；要定量用流量，得走 `stream1`/`thg_l2_asi_cal` 的定标产品。

## 5. DMSP 部分：为什么没有 2025-03-13 的图

- `ssusi.jhuapl.edu` 首页现在写着「由于经费不足，SSUSI 网站已不可用」；
- NASA SPDF 上的 SSUSI 存档（`/pub/data/dmsp/dmspf16|f17|f18/ssusi/data/edr-aurora/2025/`）最后一天是 **DOY 049 = 2025-02-18**，之后没有文件；
- CDAWeb 的 DMSP SSJ/SSIES 数据集更早：`DMSP-F18_SSJ_PRECIPITATING-ELECTRONS-IONS` 覆盖 2009–2014，`pub/data/dmsp/*` 目录停在 2019 年。

也就是说这一天既没有 SSUSI 极光图像，也没有能替代的 DMSP 粒子沉降数据；要 DMSP 原图只能走 NOAA CLASS（需要账号）。如果之后要补，可以用同一天的 DMSP SSJ 沉降能通量画 MLT–MLAT 极光沉降图，但那不是图像产品。

## 6. 想换日期/站点/时刻

`code/download_themis_asi.py` 里改 `DAY`、`SITES`、`MINUTES`，`code/download_quicklooks.py` 里改 `HOURS`，跑完再跑 `code/render_asi.py` 即可。缺文件的时刻，服务器会回 301，脚本按「该分钟不在档案里」处理。

## 7. 引用

THEMIS ASI 数据版权归 UC Berkeley / UCalgary；使用时请按 <https://data.phys.ucalgary.ca/data/how_to_cite> 致谢。
