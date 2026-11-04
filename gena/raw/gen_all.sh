set -euo pipefail

RAW_DIR=/Users/shawn/ENA_BASE/GENA_RAW
OUT_DIR=/Users/shawn/ENA_BASE/GENA_PQT
CLI=/Users/shawn/syncthing/GENA-CODE-Dev/gena/raw/raw_2_pqt.py
STAMP=$(date +%Y%m%d_%H%M%S)
PARTS_DIR="${OUT_DIR}.parts_${STAMP}"

# 1) 备份旧输出，避免和新结果混在一起
if [ -d "$OUT_DIR" ]; then
    mv "$OUT_DIR" "${OUT_DIR}.backup_${STAMP}"
fi
mkdir -p "$OUT_DIR" "$PARTS_DIR/logs"
export OUT_DIR PARTS_DIR CLI

# 2) 每个 raw 独立输出到 PARTS_DIR/<raw basename>，互不覆盖
find "$RAW_DIR" -maxdepth 1 -type f -name '*.raw' -print0 |
xargs -0 -P 4 -n 1 sh -c '
    f="$1"
    b=$(basename "$f" .raw)
    log="$PARTS_DIR/logs/${b}.log"
    printf "START %s\n" "$b"
    "$CLI" -q \
        -o "$PARTS_DIR/$b" \
        --products pdf,l1,l1a,l2,l2s,l2g \
        --chunk-packets 4096 \
        "$f" >"$log" 2>&1
    rc=$?
    printf "DONE  %s rc=%d\n" "$b" "$rc"
    exit "$rc"
' _

# 3) 按“日期 + 产品”合并临时目录中的结果
python - "$PARTS_DIR" "$OUT_DIR" <<'PY'
import sys
from pathlib import Path
import pandas as pd

parts = Path(sys.argv[1])
out = Path(sys.argv[2])
groups = {}

for f in sorted(parts.glob('*/*.pqt')):
    groups.setdefault(f.name, []).append(f)

for name, files in sorted(groups.items()):
    frames = [pd.read_parquet(f) for f in files]
    df = pd.concat(frames, ignore_index=True)

    if {'UUID', 'PID'}.issubset(df.columns):
        # 事件身份：(UUID, PID)。UUID 是稳定 packet UUID。
        df = df.drop_duplicates(subset=['UUID', 'PID'], keep='first')
    elif 'Timestamp' in df.columns:
        # PDF 没有 UUID，用包时间去重。
        df = df.drop_duplicates(subset=['Timestamp'], keep='first')

    target = out / name
    df.to_parquet(target, index=False, engine='pyarrow', compression='zstd')
    print(f'{name}: {len(df)} rows')
PY

echo "DONE: $OUT_DIR"
echo "PARTS: $PARTS_DIR"
