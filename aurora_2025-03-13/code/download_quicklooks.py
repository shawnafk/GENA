"""Download THEMIS ASI hour quicklooks (keogram + montage + hourly average) from UCalgary.

stream2 holds ready-made per-hour views:
  <YYYYMMDD>_<HH>_<site>_<imager>_full-keogram.pgm.jpg   600x256, time ->, sky row ->
  <YYYYMMDD>_<HH>_<site>_<imager>_full-montage.pgm.jpg   60 frames of the hour, 1/min
  <YYYYMMDD>_<HH>_<site>_<imager>_full-average.pgm.jpg   hour-average frame
"""
from __future__ import annotations

import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = "https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream2"
DAY = "20250313"
OUT = Path(__file__).resolve().parents[1] / "quicklook"

SITES = {
    "atha": "atha_themis02", "fsim": "fsim_themis08", "fsmi": "fsmi_themis10",
    "gako": "gako_themis20", "gill": "gill_themis19", "inuv": "inuv_themis17",
    "kapu": "kapu_themis21", "pina": "pina_themis18", "rank": "rank_themis12",
    "snkq": "snkq_themis10", "talo": "talo_themis15", "tpas": "tpas_themis05",
}
HOURS = [f"{h:02d}" for h in range(2, 14)]
KINDS = ["keogram", "montage", "average"]


def url_for(imager: str, hour: str, kind: str) -> str:
    return (f"{BASE}/{DAY[:4]}/{DAY[4:6]}/{DAY[6:]}/{imager}/ut{hour}/"
            f"{DAY}_{hour}_{imager}_full-{kind}.pgm.jpg")


def get(url: str, dest: Path, retries: int = 2) -> str:
    if dest.exists() and dest.stat().st_size > 5_000:
        return "cached"
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                data = resp.read()
            dest.write_bytes(data)
            return "ok"
        except urllib.error.HTTPError as exc:
            if exc.code in (301, 302, 404):
                return f"missing({exc.code})"
        except Exception:  # noqa: BLE001
            pass
    return "fail"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [(s, h, k) for s in SITES for h in HOURS for k in KINDS]
    stats: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {
            pool.submit(get, url_for(SITES[s], h, k), OUT / f"{DAY}_{h}_{SITES[s]}_{k}.jpg"): (s, h, k)
            for s, h, k in tasks
        }
        for fut, (s, h, k) in futures.items():
            status = fut.result()
            stats[status] = stats.get(status, 0) + 1
            if status not in ("ok", "cached"):
                print(f"{status:12s} {s} {h} {k}", flush=True)
    print("summary:", stats, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
