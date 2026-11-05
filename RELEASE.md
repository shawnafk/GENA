# GENA-CODE 发行说明

当前首个发行候选版本为 `0.2.0rc1`，基于当前源码接口。版本号由 `gena/_version.py` 定义，Python 安装元数据和 `gena.__version__` 从同一来源读取。完整 ZIP 另带生成的 `VERSION` 和 `release-manifest.json`，记录版本、源码提交及每个文件的 SHA-256。

## 下载与安装

从 GitHub Release 下载 `gena-code-0.2.0rc1-full.zip` 和同名 `.zip.sha256`。先核对 ZIP 的 SHA-256，再解压到独立版本目录。GitHub 自动生成的 Source code 压缩包只含 Git 文件，完整运行资源包含在我们构建的 full ZIP 中。

Windows 使用 Python 3.12–3.14，在解压目录运行以下命令。安装脚本先验证文件，在目录下创建 `.venv`，安装基础依赖及 Excel、3D 可选依赖，再检查依赖和导入。

```powershell
Get-FileHash .\gena-code-0.2.0rc1-full.zip -Algorithm SHA256
# 解压后进入 gena-code-0.2.0rc1
powershell -ExecutionPolicy Bypass -File .\install.ps1
.\.venv\Scripts\python.exe -c "import gena; print(gena.__version__)"
```

macOS 和 Linux 可使用以下步骤。安装需要可访问的 Python 包源，本发行物包含运行资源，Python 第三方依赖仍通过 pip 安装。依赖范围由 `pyproject.toml` 声明，首版尚未提供各平台的锁定环境和离线 wheel 集合。

```bash
python scripts/verify_release.py
python -m venv .venv
.venv/bin/python -m pip install '.[excel,3d]'
.venv/bin/python -m pip check
.venv/bin/python -c 'import gena; from load import load_data; from libs import ephemeris; print(gena.__version__)'
```

从其他目录运行分析脚本时，使用这个虚拟环境的 Python。用户观测数据通过 `FILE_BASE` 配置。升级时解压新版本并安装到它自己的虚拟环境，回滚时使用旧目录和旧环境。

## 资源与构建

完整包包含 `de432s.bsp`、`orbits_all.npz` 及轨道元数据。轨道快照截止 2026-08-27 左右，各卫星的精确覆盖见 `libs/orbits_all.meta.json`。超出覆盖的请求遵循现有报错或场景省略行为，显式轨道更新会改变本地资源，此后需要重新解压原始发行包才能通过原始指纹校验。资源来源见 `NOTICE`。

构建器只提取 Git 提交中五个运行包内的 Python 文件及明确列出的发行文件，工程 PDF、分析结果、缓存和用户数据不进入完整包。大型资源保留在 Git 外，构建时严格校验已提交的 `release-resources.json`。更新资源需要先更新并审阅这个清单，再提交和发布新版本。

```bash
# 在 Git 仓库中，工作目录干净且运行资源指纹匹配后执行
python scripts/build_release.py --ref HEAD
# 验证发行包后，对同一提交建立版本标签
git tag -a v0.2.0rc1 -m 'GENA-CODE 0.2.0rc1'
```

正式发布使用固定 tag 构建，上传完整 ZIP 和外部校验文件作为 GitHub Release 附件。构建器会拒绝已有同名输出，同一版本的内容应保持固定。验证器用于发现文件缺失和损坏，下载时应先通过可信 Release 页提供的外部指纹核对整个 ZIP。

科学结果、列名、单位、坐标系、文件级别及函数输入输出语义的变化都需要记录兼容性影响。当前 `0.x` 阶段以 minor 版本区分破坏性变更，兼容修复用 patch，发行候选用 `rc`。正式 `1.0` 后按语义化版本规则发布。
