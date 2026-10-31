# 测试规范

这份文档说明本仓库写测试时的约定和原则，面向 `tests/` 目录。目标是让测试**可读、可重复、离线、失败时能直接指出问题**，而不是追求覆盖率数字。

## 1. 目录与发现规则

pytest 的默认发现规则是：文件名匹配 `test_*.py` 或 `*_test.py`，函数名以 `test_` 开头，类名以 `Test` 开头且不能有 `__init__`。本仓库在 `pyproject.toml` 里固定了：

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-ra"
```

`-ra` 会在结尾汇总所有非通过项（skipped/xfailed/failed）的原因，所以平时只需要：

```bash
python -m pytest -q            # 全量
python -m pytest tests/test_load_real_cases.py -v
python -m pytest -q -k "themis or mms"
```

目录按用途分层，文件名直接说明它测什么：

```text
tests/
  test_raw_conversion.py      # 单元：解码/转换公式
  test_raw_pipeline.py        # 单元：分包、清单、幂等
  test_load_audit.py          # 缺陷复现（mock 边界，不依赖数据）
  test_load_real_cases.py     # 缺陷复现（真实本地数据）
  test_load_expected.py       # 期望行为（xfail，待修清单）
```

新增文件时先问自己：它属于"单元"、"真实数据集成"还是"已知缺陷清单"？三者混在一个文件里会让失败信息难以解读。

## 2. 一条测试只验证一件事

采用 arrange / act / assert 三段式，中间不要夹带无关断言：

```python
def test_gena_three_endpoints_cover_both_days():
    # arrange
    _require(PQT_BASE / "20251205_l2s_a.pqt")

    # act
    frame = load("gena", ["20251204", "20251206"], "l2s_a")

    # assert
    assert frame["Datetime"].max() >= pd.Timestamp("2025-12-05")
```

一条测试里出现多个互不相关的断言时，第一个失败会掩盖后面的信息。确实需要验证一组相关性质时，用 `assert` 之后加一句说明性消息：

```python
assert len(frame) > 0, "GENA l2s_a 应该至少有一行事件"
```

## 3. 命名：描述行为，不描述实现

推荐 `test_<对象>_<条件>_<期望>`，读起来像一句话：

```python
test_gena_two_date_endpoints_cover_only_the_first_day
test_dst_loads_from_local_cache
test_mms_fpi_without_local_files_is_labelled
```

本仓库额外使用三个前缀，含义固定，看到名字就知道测试的立场：

| 前缀 | 含义 | 期望结果 |
|---|---|---|
| `test_control_*` | 正向对照：正常路径本来就该工作 | 通过 |
| `test_claim_*` | 复现当前缺陷：断言"坏行为确实发生" | 通过（修好后应改为 `test_*` 并断言正确行为） |
| `test_fixed_*` | 修复后的期望：证明某个最小改动能让路径工作 | 通过 |
| 普通 `test_*` | 长期回归测试：断言正确行为 | 通过 |
| `test_*` + `xfail` | 已知缺陷，待修 | xfailed |

不要用 `test_1`、`test_ok`、`test_load` 这类名字；失败时日志里只有名字，名字必须自带信息。

## 4. 断言：用 `assert`，但要用对

- 布尔、存在性、形状直接用 `assert`，不要 `assert x == True` 或 `assert bool(x)`。
- 浮点用 `pytest.approx`，不要 `==`：

```python
assert values["Pressure"][0] == pytest.approx(1.0)
```

- numpy 数组用 `np.testing.assert_allclose`（会给最大误差）或 `np.array_equal`（要求严格相等）：

```python
np.testing.assert_allclose(positions, expected, rtol=1e-6)
assert np.isfinite(positions).all()
```

- 预期异常必须用 `pytest.raises`，并尽量带 `match` 锁定原因；`match` 是**正则**，括号、方括号要转义：

```python
with pytest.raises(ValueError, match="unsupported themis product 'state'"):
    load("themis", TRANGE, "state", probe="b", download=False)
```

只写 `pytest.raises(ValueError)` 也能过，但下一个人不知道你为什么期待它失败。

- 不要在测试里捕获异常再 `print`；让失败自然发生，pytest 会给出 traceback。
- 不要断言日志文本、DataFrame 的 `repr`、字典顺序这类易变细节。

## 5. 测试替身：fixture、monkeypatch、tmp_path

### 5.1 fixture 管理共享资源

```python
@pytest.fixture(scope="module")
def pqt_base():
    return Path(resolve_paths()["pqt_base"])
```

- `scope="function"`（默认）：每个测试独立，优先用它。
- `scope="module"`：只在昂贵的准备（读大文件）时用，且 fixture 必须是只读的。
- fixture 不要返回可变全局对象，否则测试之间会互相污染。

### 5.2 monkeypatch 替换边界，而不是替换实现

原则是**在哪里查找，就在哪里替换**。`load_spedas` 内部调用的名字来自 `load.loaders.spedas.read_spedas`，所以：

```python
monkeypatch.setattr(spedas_loader, "read_spedas", lambda *a, **k: ({}, {}))
```

替换 `load.api.spedas.read_spedas` 不会影响已经导入到 `load.loaders.spedas` 的引用。这个区别是本仓库 `test_load_audit.py` 里真实踩过的坑。

`monkeypatch.setitem` 适合字典和 `sys.modules`：

```python
monkeypatch.setitem(spedas_api.THEMIS_FIELDS, "state", {...})
monkeypatch.setitem(sys.modules, "pyspedas.projects", fake_pkg)
```

`monkeypatch` 会在测试结束时自动还原，不要手动 `del`/赋值回去；也不要用 `unittest.mock.patch` 和 `monkeypatch` 混用。

### 5.3 数据用 fake，不要用 MagicMock

测试解析逻辑时，返回**形状正确的普通字典/数组**比 `MagicMock` 好得多：`MagicMock` 对任何属性都返回 mock，错误会被推迟到很远的地方才暴露。

```python
def _sw_variable(value: float) -> dict:
    times = pd.DatetimeIndex([...])
    return {"x": times, "y": np.full(len(times), value, dtype=float)}
```

### 5.4 文件用 `tmp_path`

```python
def test_gena_frame_type_follows_data(tmp_path):
    path = tmp_path / "qq_eng.pqt"
    frame.to_parquet(path, engine="pyarrow")
```

不要在仓库根目录、`tests/` 或用户 home 下写文件。需要临时改 `HOME`（例如 spacepy 首次导入要写 `~/.spacepy`）时，用环境变量在调用前设置，并在文档里说明；不要提交依赖特定 home 的测试。

## 6. 测试必须是确定的、离线的

- **不联网**：所有 loader 调用显式传 `download=False`，或者 mock 掉后端。THEMIS 默认离线，MMS 默认 `local_first` 会联网，OMNI 默认会联网，写测试时一定要注意默认值。
- **不依赖当前时间**：需要"现在"时用参数注入或 `freezegun`，不要在断言里用 `datetime.now()`。
- **不依赖执行顺序**：任何测试单独跑、倒序跑都必须通过；共享状态只能通过只读 fixture 传递。
- **不依赖 CWD**：路径一律从 `resolve_paths()` 或 `tmp_path` 推导，不要 `Path("data")`。
- **不 sleep**：等待异步/线程用事件或超时断言，不要 `time.sleep` 碰运气。
- **随机数固定种子**：`np.random.default_rng(0)`，并且断言统计性质而不是具体随机值。
- **不改全局配置**：测试里不要 `configure()` 改 `FILE_BASE`；需要隔离时用 `monkeypatch.setenv` 或传 `prefix=`。

## 7. 参数化：用数据表代替复制粘贴

```python
@pytest.mark.parametrize(("probe", "snid"), [("a", 1), ("b", 1), ("b", 2)])
def test_calibration_uses_probe_and_site(probe, snid):
    ...
```

- 参数组合要覆盖**容易错的组合**（A/B 探头、轨道/地面、日期边界），而不是所有笛卡尔积。
- 用 `ids=` 或参数元组让失败信息可读：`ids=["orbit-a", "orbit-b", "ground-b"]`。
- 一个参数化用例失败时，pytest 会显示具体参数值，这比十个复制函数更好维护。

## 8. 该覆盖什么：边界优先于行数

本仓库的科学计算代码，优先测这些点，而不是追求覆盖率：

- **空输入**：`len(df) == 0`、空 time range、空 dict。
- **缺失/无效值**：NaN、`NaT`、inf、全填充值的通道、`q=0` 的合法四元数分量。
- **单元素**：只有一行事件、只有一个时间点（本仓库的 `project()` 和插值在这类输入上出过问题）。
- **形状边界**：长度 2 的 bin 边缘、1 维 vs 2 维数组、`(N, 3)` 与 `(3,)`。
- **探头与站点**：A/B、`snid=1/2`、轨道/地面标定必须成对覆盖。
- **单位与坐标**：km/m/`R_E`、GSM/GSE/ECI2000、度/弧度。
- **时间语义**：闭区间端点、日期字符串等于 00:00、时区、UTC 无时区标记。
- **幂等与副作用**：重复调用、重复写文件、`attrs` 是否丢失。

不要为了覆盖率去测 getter、常量值、`__repr__`；这类测试只会增加维护成本。

## 9. 真实数据测试：数据缺失时 `skip`，不要静默通过

本仓库有很多依赖本地数据的功能。写集成测试时：

```python
def _require(path: Path) -> Path:
    if not Path(path).exists():
        pytest.skip(f"local data missing: {path}")
    return Path(path)
```

- 数据不存在时 `skip`，并在消息里写清缺哪个文件；不要 `return`（那会让"没测"看起来像"通过"）。
- 数据存在时断言**具体数值或范围**，不要只断言 `len > 0`；例如 Dst 的采样时刻应满足 `minute == 29 and second == 30`。
- 测试用的日期窗口写在模块顶部常量里，便于整体替换。
- 真实数据测试要控制运行时间（本仓库目前整套约 7 秒）；超过几秒的用例考虑单独目录或 marker。

## 10. 已知缺陷怎么表达：`xfail`、`skip`、`--runxfail`

三种标记语义完全不同：

| 写法 | 含义 | 修好之后 |
|---|---|---|
| `skip` | 当前环境无法测试（缺数据、缺依赖） | 仍 skip，直到环境具备 |
| `xfail(strict=False)` | 已知会失败，但不关注 | XPASS，不影响结果 |
| `xfail(strict=True)` | 已知会失败，且**修好后必须回来删标记** | XPASS，被报成失败 |
| 普通测试 | 断言正确行为 | 一直通过 |

已知缺陷用 `strict=True`，这样修好时 pytest 会提醒你删除标记、把测试转成回归测试：

```python
@pytest.mark.xfail(
    strict=True,
    reason="load_pos('themis') always raises ValueError despite local state CDFs",
)
def test_load_pos_themis_returns_positions():
    ...
```

默认运行只显示 `xfailed`，想看真实错误用 `--runxfail`（pytest 会忽略 xfail 标记，真正执行并报错）：

```bash
python -m pytest tests/test_load_expected.py -q -rx    # 列出待修清单和原因
python -m pytest tests/test_load_expected.py -q --runxfail   # 看红色错误
```

注意：断言"当前坏行为"的测试（`test_claim_*` 里用 `pytest.raises` 或 `assert 坏结果`）是**故意写反**的，它们通过不代表代码正确；修复后要同步翻转断言。仓库里同时保留"断言坏行为"和"断言好行为"两份测试时，必须用命名前缀和 docstring 说清楚，否则一定会误读。

## 11. 反模式清单

- 断言实现细节（私有函数、内部字段、日志字符串），重构就红。
- 测试里 `try/except` 吞掉异常，或 `assert True` 占位。
- 用 `MagicMock` 模拟领域对象，错误被推迟到下游。
- 一个测试覆盖多个不相关分支，第一个失败掩盖其余。
- 测试之间靠执行顺序共享变量。
- 依赖网络、当前时间、当前用户目录、当前工作目录。
- 数据缺失时 `return`，让测试"通过"。
- 断言 `approx` 却给 `rel=1` 这种宽到没有意义的容差。
- 把 `xfail` 当"暂时不想修"的垃圾桶，且不写 reason。

## 12. 新增测试的步骤

1. 先写一个会失败的测试（断言期望行为），确认它因为**预期原因**失败，而不是 ImportError 或拼写错误。
2. 再写实现，直到测试通过；如果需要固定当前缺陷，另写 `test_claim_*` 并说明修复后要翻转。
3. 补一个正向对照（`test_control_*`），证明正常路径没有被误伤。
4. 跑 `python -m pytest -q`，确认没有拖慢整体、没有污染其它测试。
5. 检查测试是否独立：单独跑一次、加 `-p no:randomly`（如装了插件）或倒序跑一次。
6. 把测试名、失败信息当成给未来自己的文档来写。

## 13. 模板

```python
"""一句话说明这个文件测什么，以及依赖哪些本地数据。"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from common.config import resolve_paths
from load import load


TRANGE = ["20251204", "20251205"]
_PATHS = resolve_paths()


def _require(path: Path) -> Path:
    if not Path(path).exists():
        pytest.skip(f"local data missing: {path}")
    return Path(path)


def test_control_gena_l2s_a_loads():
    _require(Path(_PATHS["pqt_base"]) / "20251204_l2s_a.pqt")
    frame = load("gena", TRANGE, "l2s_a")
    assert len(frame) > 0
    assert {"Energy", "Tof", "Datetime"} <= set(frame.columns)
```

## 14. 提交前检查清单

- [ ] 文件名 `test_*.py`，函数名 `test_*`，名字能说明行为。
- [ ] 一条测试只验证一件事，失败信息足够定位。
- [ ] 浮点用 `pytest.approx`/`assert_allclose`，异常用 `pytest.raises(..., match=...)`。
- [ ] 没有联网、没有写仓库目录、没有依赖 CWD/HOME/当前时间。
- [ ] 数据缺失时 `skip` 并写清原因，不是 `return`。
- [ ] 已知缺陷用 `xfail(strict=True, reason=...)`，修复后删除标记。
- [ ] `python -m pytest -q` 全量通过，运行时间没有明显变长。
- [ ] 单独运行该文件也能通过。

## 15. 自检：这条测试真的有效吗

最有效的检查是**反向验证**：临时把实现改坏（改错一个符号、删掉一个判断、把 `return` 换成 `pass`），测试应该变红；如果仍然绿，说明测试没有真正覆盖那行逻辑。

另一种检查是问自己："如果这段代码被同事重写，但行为不变，这条测试还会通过吗？" 会，说明测的是行为；不会，说明测的是实现，应该改写。

## 16. 参考：本仓库现有的测试文件

| 文件 | 类型 | 默认结果 | 说明 |
|---|---|---|---|
| `tests/test_raw_conversion.py` | 单元 | 通过 | 解码公式、探头/站点标定、无效输入 |
| `tests/test_raw_pipeline.py` | 单元 | 通过 | 分包、清单、拒绝覆盖 |
| `tests/test_load_audit.py` | 缺陷复现（mock） | 通过 | 断言当前缺陷，修复后需翻转 |
| `tests/test_load_real_cases.py` | 缺陷复现（真实数据） | 通过 | 用 2025-12-04/05 的真实调用复现 |
| `tests/test_load_expected.py` | 待修清单 | xfailed | 断言正确行为，`--runxfail` 看红 |
