#!/usr/bin/env python3
"""Generate advanced-method CP2K inputs for a metal (Al, fcc) and a
semiconductor (AlAs, zincblende) at 4x4x4 Monkhorst-Pack, each with and
without k-point symmetry.

Systems and geometry:

  Al    primitive fcc cell, a = 4.05 A (conventional cubic), 4 Al atoms
  AlAs  primitive zincblende cell from validation_additional/structures.json
        (a = 4.013398 A, alpha = beta = gamma = 60.0 deg), 2 atoms

Only the basis sets and the method sections vary; those follow the
corresponding QS regtests:

  gw       tests/QS/regtest-gw-kpoints/G0W0_kpoints_in_self_energy.inp
  rpa      tests/QS/regtest-rpa-cubic-scaling/RPA_kpoints_H2O_kpsym.inp
  wannier  tests/QS/regtest-kp-6/h_fcc_wannier90_scf_mp.inp
  admm     tests/QS/regtest-kp-hfx-ri-admm-2/diamond_admmp.inp
  hfx      tests/QS/regtest-kp-hfx-ri/hBN_gpw_pbe0_kpsym.inp (ccGRB + AUTO_BASIS;
           As has no cc-TZ/RI_TZ basis, so the diamond template is not reusable)
  mp2      tests/QS/regtest-ri-laplace-mp2/RI_laplace_MP2_H2O.inp
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# systems
# ---------------------------------------------------------------------------

METHODS = ("gw", "rpa", "wannier", "admm", "hfx", "mp2")

# primitive fcc cell, from validation_additional/structures.json, 1 atom
AL_FCC = {
    "name": "al",
    "lattice": {"a": 2.8559542459167657, "b": 2.855954291634759, "c": 2.855955,
                "abc": "60.0000033779252 60.00000284838631 60.000002428463056"},
    "sites": [("Al", 0, 0, 0)],
    "smearing": True,
    "wannier_functions": 4,
    "methods": ("rpa", "wannier", "admm", "hfx"),
}

# primitive zincblende cell, from validation_additional/structures.json
ALAS_ZB = {
    "name": "alas",
    "lattice": {"a": 4.013398250535005, "b": 4.013397970922498, "c": 4.01339746,
                "abc": "59.9999959688638 59.99999827352137 59.9999965443567"},
    "sites": [("Al", 0, 0, 0), ("As", .25, .25, .25)],
    "smearing": False,
    "wannier_functions": 4,
    "methods": METHODS,
}

SYSTEMS = (AL_FCC, ALAS_ZB)

# element data: charge, MOLOPT RI_HFX basis, MOLOPT_UZH DZVP, ccGRB-D
ELEMENT = {
    "Al": {
        "q": 3,
        "ri_molopt": "RI_45_4_4_3_2_0_0_0_1.7e-06",
        "molopt_uzh": "DZVP-MOLOPT-PBE-GTH-q3",
        "ccgrb": "ccGRB-D-q3",
        "admm": "admm-dz-q3",
    },
    "As": {
        "q": 5,
        "ri_molopt": "RI_46_5_4_3_2_0_0_0_8.1e-06",
        "molopt_uzh": "DZVP-MOLOPT-PBE-GTH-q5",
        "ccgrb": "ccGRB-D-q5",
        "admm": "admm-dz-q5",
    },
}


def elements(system):
    result = []
    for el, *_ in system["sites"]:
        if el not in result:
            result.append(el)
    return result


def coords(system):
    return "\n".join(f"      {el} {x:g} {y:g} {z:g}" for el, x, y, z in system["sites"])


def cell(system):
    lat = system["lattice"]
    lines = [f"      ABC {lat['a']} {lat['b']} {lat['c']}"]
    if lat["abc"]:
        lines.append(f"      ALPHA_BETA_GAMMA {lat['abc']}")
    lines.append("      PERIODIC XYZ")
    return "    &CELL\n" + "\n".join(lines) + "\n    &END CELL"


def kind(el, basis, potential, ri_hfx=None, aux_fit=None):
    lines = [f"    &KIND {el}",
             f"      BASIS_SET {basis}",
             f"      POTENTIAL {potential}"]
    if ri_hfx:
        lines.append(f"      BASIS_SET RI_HFX {ri_hfx}")
    if aux_fit:
        lines.append(f"      BASIS_SET AUX_FIT {aux_fit}")
    lines.append("    &END KIND")
    return "\n".join(lines)


def kinds(system, builder, **kwargs):
    return [builder(el, **kwargs.get(el, {})) for el in elements(system)]


def dft_head(project, basis_files, potential_file, auto_basis=None, cutoff=600,
             rel_cutoff=60):
    bf = "\n".join(f"    BASIS_SET_FILE_NAME {b}" for b in basis_files)
    auto = f"    AUTO_BASIS RI_HFX {auto_basis}\n" if auto_basis else ""
    return f"""&GLOBAL
  PROJECT {project}
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
{bf}
    POTENTIAL_FILE_NAME {potential_file}
    SORT_BASIS EXP
{auto}    &QS
      EPS_DEFAULT 1.0E-12
      EPS_PGF_ORB 1.0E-5
      METHOD GPW
    &END QS
    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF {rel_cutoff}
    &END MGRID
"""


def kpoints(symmetry, verbose=False, eps_symmetry=False):
    lines = ["    &KPOINTS",
             "      SCHEME MONKHORST-PACK 4 4 4",
             f"      SYMMETRY {'ON' if symmetry else 'OFF'}"]
    if not symmetry:
        lines.append("      FULL_GRID ON")
    if eps_symmetry:
        lines.append("      EPS_SYMMETRY 1.0E-8")
    if verbose:
        lines.append("      VERBOSE T")
    lines.append("      WAVEFUNCTIONS COMPLEX")
    lines.append("    &END KPOINTS")
    return "\n".join(lines) + "\n"


def scf(system, eps="1.0E-9", max_scf=200, added_mos="-1", mixing=True, ignore=False,
        cholesky_off=False):
    lines = ["    &SCF",
             f"      EPS_SCF {eps}",
             f"      MAX_SCF {max_scf}",
             "      SCF_GUESS ATOMIC"]
    if cholesky_off:
        lines += ["      CHOLESKY OFF",
                  "      EPS_EIGVAL 1.0E-8"]
    lines.append(f"      ADDED_MOS {added_mos}")
    lines += ["      &DIAGONALIZATION",
              "        ALGORITHM STANDARD",
              "      &END DIAGONALIZATION"]
    if mixing:
        lines += ["      &MIXING",
                  "        METHOD BROYDEN_MIXING",
                  "        ALPHA 0.15",
                  "      &END MIXING"]
    lines += ["      &PRINT",
              "        &RESTART OFF",
              "        &END RESTART",
              "      &END PRINT"]
    if system["smearing"]:
        lines += ["      &SMEAR",
                  "        METHOD FERMI_DIRAC",
                  "        ELECTRONIC_TEMPERATURE 300",
                  "      &END SMEAR"]
    lines.append("    &END SCF")
    return "\n".join(lines) + "\n"


def tail(system, kind_blocks, extra_print=""):
    kb = "\n".join(kind_blocks)
    return f"""    &PRINT
      &KPOINTS ON
        FILENAME __STD_OUT__
      &END KPOINTS
{extra_print}    &END PRINT
  &END DFT
  &SUBSYS
{cell(system)}
    &COORD
      SCALED T
{coords(system)}
    &END COORD
{kb}
  &END SUBSYS
&END FORCE_EVAL
"""


def gw(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_gw_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_MOLOPT", "BASIS_RI_MOLOPT"], "GTH_POTENTIALS")
    kp = kpoints(symmetry)
    xc = """    &XC
      &WF_CORRELATION
        MEMORY 2000
        &INTEGRALS
          SIZE_LATTICE_SUM 3
          &WFC_GPW
            CUTOFF 200
            REL_CUTOFF 20
          &END WFC_GPW
        &END INTEGRALS
        &LOW_SCALING
          KPOINTS 4 4 4
          REGULARIZATION_RI 1.0E-3
        &END LOW_SCALING
        &RI_RPA
          RPA_NUM_QUAD_POINTS 6
          &GW
            CORR_OCC 1
            CORR_VIRT 1
            KPOINTS_SELF_ENERGY 1 1 1
            RI_SIGMA_X
          &END GW
        &END RI_RPA
      &END WF_CORRELATION
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, "DZVP-MOLOPT-SR-GTH", f"GTH-PBE-q{ELEMENT[el]['q']}",
               ri_hfx=ELEMENT[el]["ri_molopt"]) for el in elements(system)]
    return head + kp + scf(system, eps="1.0E-5", added_mos="10000000") + xc + tail(system, kb)


def rpa(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_rpa_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_MOLOPT", "BASIS_RI_MOLOPT"], "GTH_POTENTIALS")
    kp = kpoints(symmetry)
    xc = """    &XC
      &WF_CORRELATION
        MEMORY 2000
        NUMBER_PROC 1
        &INTEGRALS
          SIZE_LATTICE_SUM 3
          &WFC_GPW
            CUTOFF 200
            REL_CUTOFF 20
          &END WFC_GPW
        &END INTEGRALS
        &LOW_SCALING
          DO_EXTRAPOLATE_KPOINTS FALSE
          DO_KPOINTS
          MEMORY_CUT 1
        &END LOW_SCALING
        &RI_RPA
          RPA_NUM_QUAD_POINTS 6
        &END RI_RPA
      &END WF_CORRELATION
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, "DZVP-MOLOPT-SR-GTH", f"GTH-PBE-q{ELEMENT[el]['q']}",
               ri_hfx=ELEMENT[el]["ri_molopt"]) for el in elements(system)]
    return head + kp + scf(system, eps="1.0E-5", added_mos="10000000") + xc + tail(system, kb)


def wannier(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_wannier_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_MOLOPT_UZH"], "POTENTIAL_UZH")
    kp = kpoints(symmetry, verbose=True, eps_symmetry=True)
    w90 = f"""      &WANNIER90 ON
        KPOINTS_SOURCE SCF
        SEED_NAME adv_wannier_{name}
        WANNIER_FUNCTIONS {system['wannier_functions']}
      &END WANNIER90
"""
    xc = """    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, ELEMENT[el]["molopt_uzh"], f"GTH-PBE-q{ELEMENT[el]['q']}")
          for el in elements(system)]
    return head + kp + scf(system, eps="1.0E-9") + xc + tail(system, kb, extra_print=w90)


def admm(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_admm_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_ccGRB_UZH", "BASIS_ADMM_UZH"], "POTENTIAL_UZH",
                    auto_basis="SMALL")
    admm = """    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_PURIFICATION_METHOD NONE
      EXCH_CORRECTION_FUNC PBEX
      EXCH_SCALING_MODEL MERLOT
      METHOD BASIS_PROJECTION
    &END AUXILIARY_DENSITY_MATRIX_METHOD
"""
    kp = kpoints(symmetry)
    xc = """    &XC
      &HF
        FRACTION 1.0
        &INTERACTION_POTENTIAL
          CUTOFF_RADIUS 0.5
          POTENTIAL_TYPE TRUNCATED
        &END INTERACTION_POTENTIAL
        &RI
          MEMORY_CUT 2
          RI_METRIC IDENTITY
        &END RI
      &END HF
      &XC_FUNCTIONAL NONE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, ELEMENT[el]["ccgrb"], f"GTH-PBE-q{ELEMENT[el]['q']}",
               aux_fit=ELEMENT[el]["admm"]) for el in elements(system)]
    return head + admm + kp + scf(system, eps="1.0E-6", max_scf=100) + xc + tail(system, kb)


def hfx(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_hfx_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_ccGRB_UZH"], "POTENTIAL_UZH", auto_basis="SMALL")
    kp = kpoints(symmetry, verbose=True)
    xc = """    &XC
      &HF
        FRACTION 1.0
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE IDENTITY
        &END INTERACTION_POTENTIAL
        &RI
          EPS_FILTER 1.0E-10
          MEMORY_CUT 2
          RI_METRIC IDENTITY
        &END RI
      &END HF
      &XC_FUNCTIONAL NONE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, ELEMENT[el]["ccgrb"], f"GTH-PBE0-q{ELEMENT[el]['q']}")
          for el in elements(system)]
    return head + kp + scf(system, eps="1.0E-6", max_scf=100) + xc + tail(system, kb)


def mp2(system, symmetry):
    name = system["name"]
    head = dft_head(f"adv_mp2_{name}_{'sym' if symmetry else 'nosym'}",
                    ["BASIS_MOLOPT", "BASIS_RI_MOLOPT"], "GTH_POTENTIALS")
    kp = kpoints(symmetry)
    xc = """    &XC
      &HF
        FRACTION 1.0
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE IDENTITY
        &END INTERACTION_POTENTIAL
        &RI
          MEMORY_CUT 2
          RI_METRIC IDENTITY
        &END RI
        &SCREENING
          EPS_SCHWARZ 1.0E-8
          SCREEN_ON_INITIAL_P FALSE
        &END SCREENING
      &END HF
      &WF_CORRELATION
        MEMORY 2000
        NUMBER_PROC 1
        &INTEGRALS
          ERI_METHOD GPW
          &WFC_GPW
            CUTOFF 200
            REL_CUTOFF 20
            EPS_FILTER 1.0E-6
            EPS_GRID 1.0E-6
          &END WFC_GPW
        &END INTEGRALS
        &RI_SOS_MP2
          QUADRATURE_POINTS 5
        &END RI_SOS_MP2
      &END WF_CORRELATION
      &XC_FUNCTIONAL NONE
      &END XC_FUNCTIONAL
    &END XC
"""
    kb = [kind(el, "DZVP-MOLOPT-SR-GTH", f"GTH-PBE-q{ELEMENT[el]['q']}",
               ri_hfx=ELEMENT[el]["ri_molopt"]) for el in elements(system)]
    return head + kp + scf(system, eps="1.0E-7", added_mos="10000000") + xc + tail(system, kb)


BUILDERS = {"gw": gw, "rpa": rpa, "wannier": wannier,
            "admm": admm, "hfx": hfx, "mp2": mp2}


def main():
    count = 0
    for system in SYSTEMS:
        outdir = ROOT / system["name"]
        outdir.mkdir(parents=True, exist_ok=True)
        wanted = set()
        for method in system["methods"]:
            for symmetry in (True, False):
                name = f"{method}_{'sym' if symmetry else 'nosym'}.inp"
                wanted.add(name)
                (outdir / name).write_text(BUILDERS[method](system, symmetry))
                count += 1
        for stale in outdir.glob("*.inp"):
            if stale.name not in wanted:
                stale.unlink()
                stale.with_suffix(".out").unlink(missing_ok=True)
                stale.with_suffix(".launcher.log").unlink(missing_ok=True)
    print(f"wrote {count} input files under {ROOT}")


if __name__ == "__main__":
    main()
