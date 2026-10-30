#!/usr/bin/env python3
"""Ambiguity of every in-house dataset (29 populations) with the in-house settings.

Settings: global SOAP basis Co Cr Li Na O P V (identical to the per-structure basis after
normalisation), K = 100, alpha 1.25, m 5, same-frame candidates, signed moments, T8 with the
20 eV/A frame-force filter, no |m| filter.  Each dataset is an independent neighbour pool.
Resumable: datasets with a complete metadata.json of the same settings version are skipped.

    python scripts/run_inhouse_all.py --out results/inhouse [--datasets cro2_u licoo2_u] [--n-jobs 12]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from magnetic_ambiguity.cli import add_setting_args, settings_from_args  # noqa: E402
from magnetic_ambiguity.inhouse import BY_KEY, DATASETS  # noqa: E402
from magnetic_ambiguity.pipeline import run_file  # noqa: E402
from magnetic_ambiguity.settings import INHOUSE_SOAP_SPECIES  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "inhouse")
    ap.add_argument("--datasets", nargs="*", default=None, help="keys; default all 29")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--query-chunk", type=int, default=2048)
    add_setting_args(ap)
    args = ap.parse_args()
    if args.species is None:
        args.species = list(INHOUSE_SOAP_SPECIES)
    settings = settings_from_args(args)
    version = settings.version()
    keys = args.datasets or [d.key for d in DATASETS]
    run_dir = args.out / version
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "settings.json").write_text(json.dumps(settings.to_dict(), indent=2, default=str) + "\n")
    rows = []
    for key in keys:
        d = BY_KEY[key]
        out = run_dir / key
        meta = out / "metadata.json"
        if meta.exists() and json.loads(meta.read_text()).get("settings_version") == version:
            print(f"{key}: complete, skipping", flush=True)
            summary = pd.read_csv(out / "summary.csv").iloc[0].to_dict()
        else:
            summary = run_file(d.path, d.element, settings, out, n_jobs=args.n_jobs, query_chunk=args.query_chunk)
        rows.append({"dataset": key, "label": d.label, "material": d.material, "analysis_group": d.analysis_group,
                     "hubbard_u": d.hubbard_u, "temperature_K": d.temperature_K, **summary})
        pd.DataFrame(rows).to_csv(run_dir / "summary_all.csv", index=False, float_format="%.10g")
    t = pd.DataFrame(rows)
    tcol = f"T{round(settings.dataset.tail_fraction * 100)}_muB"
    pd.set_option("display.width", 200)
    print(t[["dataset", "element", "n_sites", tcol, "A_muB", "median_muB", "k_eff_median"]].to_string(index=False))
    print(f"\nwritten {run_dir / 'summary_all.csv'}")


if __name__ == "__main__":
    main()
