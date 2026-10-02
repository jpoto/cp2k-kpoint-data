#!/usr/bin/env python3
"""Snap the lattice metric of every structure to its detected point group.

The Materials-Project cells shipped in `regtest_inputs.tar.gz` differ from ideal
symmetry by ~1e-6 (in the lattice parameters).  CP2K's symmetry detection with
`EPS_SYMMETRY 1e-5` treats them as symmetric, so a symmetry-reduced run
integrates a *symmetrized* cell while a `FULL_GRID` run integrates the actual,
very slightly distorted cell.  On cells where the k-set changes a lot (primitive
cubic Al/W, body-centred BaFe2As2) this shows up as a ~1e-7 Ha energy and
sub-bar-to-few-bar stress difference between the reduced and full runs.

This script removes that methodological discrepancy at the source: for each
structure it builds the spglib symmetry dataset, averages the metric tensor

    g' = (1/|G|) sum_R R^T g R

over the space-group rotations `R`, and writes the resulting lattice parameters
back to `structures.json`.  `g'` is exactly invariant under the detected point
group, so the full and reduced runs now integrate the *same* cell.  Atomic
coordinates are left untouched: spglib detects the intended space group for 11
of the 12 cells directly (LiMg_imma is genuinely lower-symmetry at this
tolerance and is not a large-error case).

The original (as-given) lattice is preserved under `lattice_asgiven` for
traceability.  Run with `--check` to only report what would change.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import spglib

ROOT = Path(__file__).resolve().parent
DEG = 180.0 / np.pi

ELEMENTS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,
    "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51,
    "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
}


def element_of(kind):
    if kind.get("element"):
        return kind["element"]
    return kind["name"].split("_")[0]


def lattice_matrix(a, b, c, alpha, beta, gamma):
    """Rows are lattice vectors (spglib convention)."""
    al, be, ga = np.radians([alpha, beta, gamma])
    va = [a, 0.0, 0.0]
    vb = [b * np.cos(ga), b * np.sin(ga), 0.0]
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([va, vb, [cx, cy, cz]])


def metric_params(g):
    a, b, c = np.sqrt(np.diag(g))
    alpha = np.degrees(np.arccos(np.clip(g[1, 2] / (b * c), -1.0, 1.0)))
    beta = np.degrees(np.arccos(np.clip(g[0, 2] / (a * c), -1.0, 1.0)))
    gamma = np.degrees(np.arccos(np.clip(g[0, 1] / (a * b), -1.0, 1.0)))
    return a, b, c, alpha, beta, gamma


def symmetrized_metric(lattice, frac, numbers, symprec):
    dataset = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=symprec)
    if dataset is None:
        raise RuntimeError("spglib could not determine a symmetry dataset")
    g = lattice @ lattice.T
    gs = np.zeros((3, 3))
    for rotation in dataset.rotations.astype(float):
        gs += rotation.T @ g @ rotation
    gs /= len(dataset.rotations)
    return (gs + gs.T) / 2, dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--structures", type=Path, default=ROOT / "structures.json")
    parser.add_argument("--symprec", type=float, default=1.0e-5,
                        help="spglib symprec used to detect the point group")
    parser.add_argument("--check", action="store_true",
                        help="report changes without writing")
    args = parser.parse_args()

    doc = json.loads(args.structures.read_text())
    for system in doc["systems"]:
        lat = system["lattice"]
        lattice = lattice_matrix(lat["a"], lat["b"], lat["c"],
                                 lat["alpha"], lat["beta"], lat["gamma"])
        frac = [[x, y, z] for _, x, y, z in system["coords"]]
        # map every atom to its kind's atomic number
        kind_num = {kind["name"]: ELEMENTS[element_of(kind)]
                    for kind in system["kinds"]}
        numbers = [kind_num[name] for name, _, _, _ in system["coords"]]

        gs, dataset = symmetrized_metric(lattice, frac, numbers, args.symprec)
        a, b, c, alpha, beta, gamma = metric_params(gs)
        new = {"a": round(a, 10), "b": round(b, 10), "c": round(c, 10),
               "alpha": round(alpha, 8), "beta": round(beta, 8),
               "gamma": round(gamma, 8)}
        old = {k: lat[k] for k in ("a", "b", "c", "alpha", "beta", "gamma")}
        maxdiff = max(abs(new[k] - old[k]) for k in old)
        status = "OK" if dataset.number == system["space_group"] else "DIFF"
        print(f"{status:4s} {system['label']:13s} sg={dataset.number:3d} "
              f"{dataset.international:8s} max|dlat|={maxdiff:.2e} "
              f"a={new['a']:.7f} b={new['b']:.7f} c={new['c']:.7f} "
              f"al={new['alpha']:.6f} be={new['beta']:.6f} ga={new['gamma']:.6f}")
        if not args.check:
            system.setdefault("lattice_asgiven", old)
            system["lattice"] = new

    if args.check:
        print("(check only; no file written)")
        return
    doc["_comment"] = (
        "Structures from the pinned regtest_*.inp inputs (Materials Project), "
        "expressed as canonical ABC/ALPHA_BETA_GAMMA lattice parameters. The "
        "lattice metric has been snapped to the spglib-detected point group "
        "(see symmetrize_structures.py) so that symmetry-reduced and FULL_GRID "
        "runs integrate the same cell; the as-given parameters are kept under "
        "'lattice_asgiven'. Coordinates are fractional (SCALED T). Expected IBZ "
        "counts come from expected_ibz_counts.json.")
    args.structures.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"Wrote {args.structures}")


if __name__ == "__main__":
    main()
