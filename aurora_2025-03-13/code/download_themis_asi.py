"""Download THEMIS ASI all-sky frames for 2025-03-13 (UTC) from the UCalgary open archive.

Archive layout (one file per minute, ~2.3 MB gzip, 20 x 256x256 16-bit frames at 3 s):
  https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream0/
      <YYYY>/<MM>/<DD>/<site>_themisNN/ut<HH>/<YYYYMMDD>_<HHMM>_<site>_themisNN_full.pgm.gz

Missing minutes are served as an nginx 301 with no Location header, so a 301 means
"not in the archive" and is not retried.
"""
from __future__ import annotations

import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = "https://data.phys.ucalgary.ca/sort_by_project/THEMIS/asi/stream0"
DAY = "20250313"
RAW = Path(__file__).resolve().parents[1] / "_raw"

SITES = {
    "gill": "gill_themis19",
    "kapu": "kapu_themis21",
    "rank": "rank_themis12",
    "fsmi": "fsmi_themis10",
    "talo": "talo_themis15",
    "inuv": "inuv_themis17",
    "gako": "gako_themis20",
    "snkq": "snkq_themis10",
    "atha": "atha_themis02",
    "pina": "pina_themis18",
    "tpas": "tpas_themis05",
    "fsim": "fsim_themis08",
}

MINUTES = ["0400", "0445", "0530", "0615", "0700", "0745", "0830", "0915",
           "1000", "1045", "1130", "0800"]


def url_for(site: str, minute: str) -> str:
    imager = SITES[site]
    return f"{BASE}/{DAY[:4]}/{DAY[4:6]}/{DAY[6:]}/{imager}/ut{minute[:2]}/{DAY}_{minute}_{imager}_full.pgm.gz"


def get(url: str, dest: Path, retries: int = 2) -> str:
    if dest.exists() and dest.stat().st_size > 50_000:
        return "cached"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                data = resp.read()
            dest.write_bytes(data)
            return "ok"
        except urllib.error.HTTPError as exc:
            if exc.code in (301, 302, 404):
                return f"missing({exc.code})"
            last = f"http{exc.code}"
        except Exception as exc:  # noqa: BLE001
            last = type(exc).__name__
    return "fail:" + last


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    tasks = [(s, m) for s in SITES for m in MINUTES]
    stats: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {
            pool.submit(get, url_for(s, m), RAW / f"{DAY}_{m}_{SITES[s]}_full.pgm.gz"): (s, m)
            for s, m in tasks
        }
        for fut in futures:
            site, minute = futures[fut]
            status = fut.result()
            stats[status.split("(")[0]] = stats.get(status.split("(")[0], 0) + 1
            if status != "ok":
                print(f"{status:12s} {site:5s} {minute}", flush=True)
    print("\nsummary:", stats, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
