"""Command line: `python -m magnetic_ambiguity <command> ...`

  run        one extxyz pool, one element -> out_dir/{sites.parquet, neighbors.npz, summary.csv, metadata.json}
  rescore    recompute a_i and T_p from stored neighbours with another alpha / m / fixed k (no new kNN)
  summarize  dataset statistics from an existing sites.parquet (other tail fraction / force filter)
  settings   print the default (or modified) settings and their version hash

Every scientific parameter is an option; the defaults are the 2026-10-28 definition.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .settings import (AmbiguitySettings, DatasetSettings, NeighborSettings, PopulationSettings,
                       Settings, SoapSettings, INHOUSE_SOAP_SPECIES)


def _none_or_float(text: str):
    return None if text.lower() in ("none", "off", "nan", "") else float(text)


def add_setting_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("SOAP")
    g.add_argument("--r-cut", type=float, default=6.0)
    g.add_argument("--n-max", type=int, default=8)
    g.add_argument("--l-max", type=int, default=3)
    g.add_argument("--species", nargs="*", default=None,
                   help="global SOAP species basis; omit for the per-structure basis (identical after "
                        "normalisation); 'inhouse' = Co Cr Li Na O P V")
    g = p.add_argument_group("neighbours")
    g.add_argument("--k", type=int, default=100, help="exact neighbours stored / adaptive cap (default 100)")
    g.add_argument("--alpha", type=float, default=1.25, help="twins: d_ij <= alpha * d_i,m (default 1.25)")
    g.add_argument("--m-ref", type=int, default=5, help="reference neighbour and minimum k_eff (default 5)")
    g.add_argument("--exclude-same-frame", action="store_true",
                   help="previous rule: exclude every site of the query frame (with --k 20 = cross-frame K20)")
    g = p.add_argument_group("ambiguity / dataset")
    g.add_argument("--moments", choices=("signed", "unsigned"), default="signed")
    g.add_argument("--tail", type=float, default=0.08, help="T_p tail fraction p (default 0.08 -> T8)")
    g.add_argument("--extra-tails", type=float, nargs="*", default=[0.20])
    g.add_argument("--fmax", type=_none_or_float, default=20.0,
                   help="frame-force filter on query sites, eV/A (default 20; 'none' disables)")
    g.add_argument("--tau", type=float, default=0.5, help="threshold of P(a_i > tau), muB")
    g = p.add_argument_group("population")
    g.add_argument("--max-abs-moment", type=_none_or_float, default=None,
                   help="drop whole frames with any |m| above this (muB); MatPES/MP-MM used 8")
    g.add_argument("--magmom-key", default=None)
    g.add_argument("--force-key", default=None)


def settings_from_args(a: argparse.Namespace) -> Settings:
    species = a.species
    if species is not None:
        species = tuple(INHOUSE_SOAP_SPECIES) if species == ["inhouse"] else tuple(species)
    return Settings(
        soap=SoapSettings(r_cut=a.r_cut, n_max=a.n_max, l_max=a.l_max, species=species),
        neighbors=NeighborSettings(k_exact=a.k, alpha=a.alpha, m_ref=a.m_ref, exclude_same_frame=a.exclude_same_frame),
        ambiguity=AmbiguitySettings(moments=a.moments),
        dataset=DatasetSettings(tail_fraction=a.tail, extra_tail_fractions=tuple(a.extra_tails),
                                fmax_filter_eV_A=a.fmax, tau_muB=a.tau),
        population=PopulationSettings(abs_moment_max_muB=a.max_abs_moment, magmom_key=a.magmom_key,
                                      force_key=a.force_key),
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="magnetic_ambiguity", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("run", help="SOAP + exact kNN + adaptive a_i + T_p for one pool")
    p.add_argument("--input", required=True, type=Path, help="extxyz with per-atom moments (and DFT forces)")
    p.add_argument("--element", required=True)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--n-jobs", type=int, default=1, help="dscribe threads")
    p.add_argument("--query-chunk", type=int, default=4096, help="query sites per distance block (memory knob)")
    p.add_argument("--no-neighbors", action="store_true", help="do not store neighbors.npz")
    add_setting_args(p)

    p = sub.add_parser("rescore", help="other band rule from stored neighbours")
    p.add_argument("--dir", required=True, type=Path, help="output directory of a previous run")
    p.add_argument("--fixed-k", type=int, default=None)
    p.add_argument("--out", type=Path, default=None, help="write sites_rescored.parquet/summary_rescored.csv here")
    add_setting_args(p)

    p = sub.add_parser("summarize", help="dataset statistics from sites.parquet")
    p.add_argument("--sites", required=True, type=Path)
    p.add_argument("--column", default="ambiguity_muB")
    add_setting_args(p)

    p = sub.add_parser("settings", help="print settings and version hash")
    add_setting_args(p)

    a = ap.parse_args(argv)
    settings = settings_from_args(a)

    if a.command == "settings":
        print(json.dumps({"settings_version": settings.version(), "soap_version": settings.soap_version(),
                          **settings.to_dict()}, indent=2, default=str))
        return 0
    if a.command == "run":
        from .pipeline import run_file
        summary = run_file(a.input, a.element, settings, a.out, n_jobs=a.n_jobs, query_chunk=a.query_chunk,
                           save_neighbors=not a.no_neighbors)
        print(json.dumps(summary, indent=2, default=str))
        return 0
    if a.command == "rescore":
        from .pipeline import rescore
        sites, summary = rescore(a.dir, settings, alpha=a.alpha, m_ref=a.m_ref, fixed_k=a.fixed_k)
        if a.out:
            a.out.mkdir(parents=True, exist_ok=True)
            sites.to_parquet(a.out / "sites_rescored.parquet", index=False)
            pd.DataFrame([summary]).to_csv(a.out / "summary_rescored.csv", index=False)
        print(json.dumps(summary, indent=2, default=str))
        return 0
    if a.command == "summarize":
        from .pipeline import summarize_sites
        sites = pd.read_parquet(a.sites)
        print(json.dumps(summarize_sites(sites, settings, column=a.column), indent=2, default=str))
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
