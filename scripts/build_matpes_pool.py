#!/usr/bin/env python3
"""Stream a MatPES JSON release into one extxyz pool file per target element.

Population rules (source project matpes/README.md, "Population rules"):
  * one pool per element: every frame containing >= 1 atom of the element, all such atoms are sites;
  * drop a WHOLE frame if any atomic |magmom| > 8 muB (--max-abs-moment);
  * do not append the separate isolated-atom reference file;
  * duplicates retained; signed site moments retained (source final collinear moments);
  * abort on missing / non-finite moments or partial occupancies;
  * parent id = MatPES id without its final underscore-delimited step (``..._999484_73`` -> ``..._999484``).
Per-atom arrays written: dft_magmom (muB), REF_forces (eV/A).  Frame info: matpes_id, parent_id,
source_record, functional, REF_energy, formula.  The output is consumed by
``python -m magnetic_ambiguity run --input <pool>.extxyz --element <El> --max-abs-moment 8``
(the |m| filter is already applied here; passing it again is harmless).

    python scripts/build_matpes_pool.py --source /path/to/MatPES-PBE-2025.2-charges.json \
        --out data/matpes/pools/pbe --elements V Co Fe Cr Mn Ni Cu Ti Nb Mo Tc Ru Rh Pd
Streams the 2 GB JSON once with ijson.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter
from pathlib import Path

import ijson
import numpy as np
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ELEMENTS = ("V", "Co", "Fe", "Cr", "Mn", "Ni", "Cu", "Ti", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd")


def parent_id(matpes_id: str | None, record: int) -> str:
    if not matpes_id:
        return f"source-record-{record}"
    head, sep, tail = matpes_id.rpartition("_")
    return head if sep and tail.isdigit() else matpes_id


def frame_from_record(obj: dict, record: int) -> tuple[Atoms, np.ndarray]:
    sites = obj["structure"]["sites"]
    if len(sites) != int(obj["nsites"]):
        raise RuntimeError(f"record {record}: nsites mismatch")
    symbols, moments, xyz = [], [], []
    for i, site in enumerate(sites):
        species = site.get("species", [])
        if len(species) != 1 or float(species[0].get("occu", 0.0)) != 1.0:
            raise RuntimeError(f"record {record} atom {i}: partial/multiple occupancy")
        m = site.get("properties", {}).get("magmom")
        if m is None or not math.isfinite(float(m)):
            raise RuntimeError(f"record {record} atom {i}: missing/non-finite magmom")
        symbols.append(str(species[0]["element"]))
        moments.append(float(m))
        xyz.append([float(v) for v in site["xyz"]])
    cell = np.asarray(obj["structure"]["lattice"]["matrix"], float)
    atoms = Atoms(symbols=symbols, positions=np.asarray(xyz), cell=cell, pbc=True)
    moments = np.asarray(moments, float)
    mid = obj.get("matpes_id")
    atoms.info.update({"matpes_id": mid or f"source-record-{record}", "parent_id": parent_id(mid, record),
                       "source_record": record, "functional": str(obj.get("functional", "")),
                       "REF_energy": float(obj["energy"]), "formula": str(obj.get("formula_pretty", ""))})
    atoms.arrays["dft_magmom"] = moments
    atoms.arrays["REF_forces"] = np.asarray(obj["forces"], float)
    return atoms, moments


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--elements", nargs="+", default=list(DEFAULT_ELEMENTS))
    ap.add_argument("--max-abs-moment", type=float, default=8.0)
    ap.add_argument("--max-records", type=int, default=None, help="bounded test run")
    ap.add_argument("--flush-every", type=int, default=2000, help="frames per append to each pool file")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    elements = list(dict.fromkeys(args.elements))
    buffers = {el: [] for el in elements}
    counts = Counter()
    sites = Counter()
    paths = {el: args.out / f"{el}.extxyz" for el in elements}
    for p in paths.values():
        if p.exists():
            raise FileExistsError(f"{p} exists; remove it or choose another --out")

    def flush(el):
        if buffers[el]:
            write(paths[el], buffers[el], format="extxyz", append=True)
            buffers[el].clear()

    with open(args.source, "rb") as fh:
        for record, obj in enumerate(ijson.items(fh, "item", use_float=True)):
            if args.max_records is not None and record >= args.max_records:
                break
            counts["records_seen"] += 1
            atoms, moments = frame_from_record(obj, record)
            if np.any(np.abs(moments) > args.max_abs_moment):
                counts["frames_abs_moment_dropped"] += 1
                continue
            present = set(atoms.get_chemical_symbols()) & set(elements)
            if not present:
                counts["frames_without_target_elements"] += 1
                continue
            counts["frames_selected"] += 1
            for el in present:
                buffers[el].append(atoms)
                sites[el] += int(np.sum(np.asarray(atoms.get_chemical_symbols()) == el))
                if len(buffers[el]) >= args.flush_every:
                    flush(el)
            if record % 20000 == 0 and record:
                print(f"  record {record}: {dict(counts)}", flush=True)
    for el in elements:
        flush(el)
    manifest = {"source": str(args.source.resolve()), "source_size_bytes": args.source.stat().st_size,
                "elements": elements, "counts": dict(counts), "target_sites": dict(sites),
                "population_rules": {"whole_frame_abs_moment_max_muB": args.max_abs_moment, "duplicates": "retained",
                                     "separate_isolated_atoms": "not appended",
                                     "missing_or_nonfinite_magmom": "fatal", "partial_occupancy": "fatal"},
                "arrays": {"dft_magmom": "muB, source final collinear moment", "REF_forces": "eV/A"},
                "pools": {el: str(paths[el]) for el in elements}}
    (args.out / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["counts"], indent=1)); print(json.dumps(manifest["target_sites"], indent=1))


if __name__ == "__main__":
    main()
