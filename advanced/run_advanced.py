#!/usr/bin/env python3
"""Run the advanced-method inputs and write outputs under this directory.

Each case writes <system>/<method>_<sym>.out next to its input.  A machine
readable summary is written to results.json.  Runs are executed by a small
process pool so two cases can share the two GPUs.
"""
import argparse
import json
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def parse_output(text, returncode, elapsed):
    def last(pattern):
        m = re.findall(pattern, text, re.M)
        return float(m[-1]) if m else None

    aborted = "[ABORT]" in text
    error = re.search(r"\[ABORT\]\s*\n?\s*\*?\s*(.+)", text)
    return {
        "converged": "SCF run converged" in text,
        "aborted": aborted,
        "abort_message": error.group(1).strip().strip("* ") if error else None,
        "program_ended": "PROGRAM ENDED AT" in text,
        "energy_hartree": last(r"ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$"),
        "nkpoints": last(r"BRILLOUIN\| List of Kpoints.*?\]\s+(\d+)"),
        "returncode": returncode,
        "wall_seconds": elapsed,
    }


def run_case(args):
    binary, env, system, method, symmetry, timeout, mpiranks = args
    name = f"{method}_{'sym' if symmetry else 'nosym'}"
    directory = ROOT / system.lower()
    inp = directory / f"{name}.inp"
    out = directory / f"{name}.out"
    launcher = directory / f"{name}.launcher.log"
    if mpiranks > 1:
        cmd = ["mpiexec", "-n", str(mpiranks), "--bind-to", "none",
               binary, "-i", inp.name, "-o", out.name]
    else:
        cmd = [binary, "-i", inp.name, "-o", out.name]
    started = time.monotonic()
    try:
        with launcher.open("w") as handle:
            proc = subprocess.run(cmd, cwd=directory, env=env, stdout=handle,
                                  stderr=subprocess.STDOUT, timeout=timeout)
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        rc = "timeout"
    elapsed = round(time.monotonic() - started, 1)
    status = parse_output(out.read_text(errors="replace"), rc, elapsed)
    status.update({"case": name, "system": system, "method": method,
                   "symmetry": symmetry})
    return status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="/workspace/install/bin/cp2k.psmp")
    parser.add_argument("--data-dir", default="/workspace/data")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--mpiranks", type=int, default=2)
    parser.add_argument("--only", default=None, help="comma list e.g. wannier,admm")
    args = parser.parse_args()

    env = os.environ | {"OMP_NUM_THREADS": str(args.threads),
                        "OPENBLAS_NUM_THREADS": "1",
                        "CP2K_DATA_DIR": args.data_dir}
    methods = args.only.split(",") if args.only else None
    cases = []
    for system in ("al", "alas"):
        present = {p.stem.rsplit("_", 1)[0] for p in (ROOT / system).glob("*.inp")}
        for method in ("gw", "rpa", "wannier", "admm", "hfx", "mp2"):
            if method not in present:
                continue
            if methods and method not in methods:
                continue
            for symmetry in (True, False):
                cases.append((args.binary, env, system, method, symmetry,
                              args.timeout, args.mpiranks))

    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_case, c): c for c in cases}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps(result), flush=True)

    results.sort(key=lambda r: (r["system"], r["method"], r["symmetry"]))
    (ROOT / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    ok = sum(1 for r in results if r["converged"] and r["program_ended"])
    print(f"\n{ok}/{len(results)} converged and ended cleanly", flush=True)


if __name__ == "__main__":
    main()
