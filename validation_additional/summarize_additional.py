#!/usr/bin/env python3
"""Summarise the additional k-point IBZ cases as Markdown and LaTeX tables.

Reads `structures.json` (space group, sites, expected irreducible counts) and
`results_additional.json` (measured counts), then writes:

  * `additional_tables.md`   — human-readable tables
  * `additional_tables.tex`  — matching LaTeX tables

The central quantity is the fraction of the full Monkhorst-Pack grid that the
point-group plus time-reversal reduction avoids:

    avoided = 1 - N_IBZ / N_full

For the 4x4x4 mesh N_full = 64 and for 6x6x6 N_full = 216.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FULL = {"4x4x4-shifted": 64, "4x4x4-gamma": 64, "6x6x6-shifted": 216, "6x6x6-gamma": 216}
MESH_ORDER = ["4x4x4-shifted", "4x4x4-gamma", "6x6x6-shifted", "6x6x6-gamma"]


def load():
    doc = json.loads((ROOT / "structures.json").read_text())
    measured = {}
    results_path = ROOT / "results_additional.json"
    if results_path.exists():
        for r in json.loads(results_path.read_text()):
            measured[(r["label"], r["mesh"])] = r
    return doc, measured


def expected_counts(system):
    """Kind-resolved counts (magnetic KIND splits lower the symmetry)."""
    return system.get("expected_ibz_kind_resolved", system["expected_ibz"])


def fmt_avoided(nk, full):
    if nk is None:
        return "--"
    return f"{100.0 * (1.0 - nk / full):.1f}%"


def markdown(doc, measured):
    lines = ["# Additional k-point IBZ reduction cases", "",
             "Structures from `regtest_inputs.tar.gz` (Materials Project). "
             "`N_full` is the unreduced Monkhorst-Pack/MacDonald grid; `N_IBZ` is the "
             "point-group + time-reversal reduced count. `avoided` is the fraction of "
             "k-points that symmetry removes, `1 - N_IBZ/N_full`.", "",
             "## Systems", "",
             "| System | Formula | Space group | Crystal system | Sites |",
             "| --- | --- | --- | --- | --- |"]
    for s in doc["systems"]:
        lines.append(f"| {s['label']} | {s['formula']} | #{s['space_group']} "
                     f"{s['space_group_symbol']} | {s['crystal_system']} | {s['nsites']} |")

    lines += ["", "## Irreducible k-points per mesh", "",
              "| System | Space group | " + " | ".join(MESH_ORDER) + " |",
              "| --- | --- | " + " | ".join("---" for _ in MESH_ORDER) + " |"]
    for s in doc["systems"]:
        cells = []
        for mesh in MESH_ORDER:
            nk = expected_counts(s)[mesh.replace("-", "_")]
            meas = measured.get((s["label"], mesh))
            if meas and meas.get("nk") is not None:
                cells.append(f"{int(meas['nk'])}/{FULL[mesh]}")
            else:
                cells.append(f"{nk}/{FULL[mesh]}")
        lines.append(f"| {s['label']} | #{s['space_group']} {s['space_group_symbol']} | "
                     + " | ".join(cells) + " |")

    lines += ["", "## Fraction of the full grid avoided by symmetry", "",
              "| System | Space group | " + " | ".join(MESH_ORDER) + " |",
              "| --- | --- | " + " | ".join("---" for _ in MESH_ORDER) + " |"]
    for s in doc["systems"]:
        cells = [fmt_avoided(expected_counts(s)[m.replace("-", "_")], FULL[m])
                 for m in MESH_ORDER]
        lines.append(f"| {s['label']} | #{s['space_group']} {s['space_group_symbol']} | "
                     + " | ".join(cells) + " |")

    lines += ["", "## Points avoided per mesh", "",
              "Absolute number of k-points removed from the full grid, `N_full - N_IBZ`.", "",
              "| System | Space group | " + " | ".join(MESH_ORDER) + " |",
              "| --- | --- | " + " | ".join("---" for _ in MESH_ORDER) + " |"]
    for s in doc["systems"]:
        cells = []
        for m in MESH_ORDER:
            nk = expected_counts(s)[m.replace("-", "_")]
            cells.append(f"{FULL[m] - nk} / {FULL[m]}")
        lines.append(f"| {s['label']} | #{s['space_group']} {s['space_group_symbol']} | "
                     + " | ".join(cells) + " |")

    if measured:
        lines += ["", "## Measured CP2K results", "",
                  "`N_IBZ` is the k-point list CP2K printed; `expected` is the spglib count. "
                  "For the anti-ferromagnetic CoF2/BaFe2As2 cases the expected count is the "
                  "*kind-resolved* one, because splitting an element into two magnetic KINDs "
                  "intentionally lowers the point group.", "",
                  "| System | Mesh | N_IBZ | expected | count | SCF |",
                  "| --- | --- | --- | --- | --- | --- |"]
        for s in doc["systems"]:
            for mesh in MESH_ORDER:
                r = measured.get((s["label"], mesh))
                if not r:
                    continue
                exp = r.get("expected_nk")
                if exp is None:
                    exp = s.get("expected_ibz_kind_resolved", s["expected_ibz"])[mesh.replace("-", "_")]
                nk = "--" if r.get("nk") is None else str(int(r["nk"]))
                count = "yes" if r.get("nk") is not None and int(r["nk"]) == exp else "NO"
                if r.get("converged"):
                    scf = "converged"
                elif r.get("nk") is not None:
                    scf = "not converged"
                else:
                    scf = "aborted"
                lines.append(f"| {s['label']} | {mesh} | {nk} | {exp} | {count} | {scf} |")

        hard = [(s["label"], m) for s in doc["systems"] for m in MESH_ORDER
                if (r := measured.get((s["label"], m))) and not r.get("converged")]
        if hard:
            lines += ["", "## Runs that did not converge", "",
                      "These inputs need SCF/mixing improvements before they can be used as "
                      "regression cases; the k-point list itself is still valid.", ""]
            for label, mesh in hard:
                lines.append(f"- `inputs/{label}-{mesh}-sym/{label}-{mesh}-sym.inp`")

    settings = ["", "## SCF settings used", "",
                "| System | smearing | mixer | alpha | Nbuffer | EPS_SCF | MAX_SCF |",
                "| --- | --- | --- | --- | --- | --- | --- |"]
    for s in doc["systems"]:
        scf = s["scf"]
        smear = scf.get("smearing")
        smear_txt = f"FERMI_DIRAC {smear['temperature_k']} K" if smear else "none"
        mix = scf.get("mixing")
        if mix:
            mixer = mix["method"].replace("_MIXING", "")
            alpha, nbuf = str(mix["alpha"]), str(mix["nbuffer"])
        else:
            mixer, alpha, nbuf = "default DIIS", "--", "--"
        settings.append(f"| {s['label']} | {smear_txt} | {mixer} | {alpha} | {nbuf} | "
                        f"{scf['eps_scf']:.0e} | {scf['max_scf']} |")
    lines += settings
    return "\n".join(lines) + "\n"


def latex(doc, measured):
    out = []
    out.append(r"\begin{table}[htbp]")
    out.append(r"\centering\small\setlength{\tabcolsep}{4pt}")
    out.append(r"\caption{Structures used for the additional k-point irreducible-Brillouin-zone "
               r"tests, with the fraction of a full $4^3$ and $6^3$ Monkhorst-Pack/MacDonald "
               r"grid removed by point-group and time-reversal symmetry. Shifts are "
               r"MacDonald $(1/2,1/2,1/2)$; ``$\Gamma$'' denotes the Gamma-centred mesh.}")
    out.append(r"\label{tab:additional-ibz}")
    out.append(r"\begin{tabular}{lllrrrr}")
    out.append(r"\toprule")
    out.append(r"System & Space group & Crystal & \multicolumn{2}{c}{$4^3$} & "
               r"\multicolumn{2}{c}{$6^3$} \\")
    out.append(r"\cmidrule(lr){4-5}\cmidrule(lr){6-7}")
    out.append(r" & & system & shifted & $\Gamma$ & shifted & $\Gamma$ \\")
    out.append(r"\midrule")
    for s in doc["systems"]:
        cells = [fmt_avoided(expected_counts(s)[m.replace("-", "_")], FULL[m])
                 for m in MESH_ORDER]
        cells = [c.replace("%", r"\%") for c in cells]
        sg = s["space_group_symbol"].replace("_", r"\_")
        out.append(f"{s['label'].replace('_', r'\_')} & {sg} & {s['crystal_system']} & "
                   + " & ".join(cells) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    out.append(r"\end{table}")
    return "\n".join(out) + "\n"


def main():
    doc, measured = load()
    md = markdown(doc, measured)
    tex = latex(doc, measured)
    (ROOT / "additional_tables.md").write_text(md)
    (ROOT / "additional_tables.tex").write_text(tex)
    print(md)


if __name__ == "__main__":
    main()
