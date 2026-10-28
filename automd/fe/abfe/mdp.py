"""Per-window free-energy MDP generation."""
from __future__ import annotations

from pathlib import Path

from .rundir import write_mdp
from .schedule import Schedule, validate

_COMMON = ["cutoff-scheme", "nstlist", "pbc", "coulombtype", "rcoulomb", "vdwtype",
           "vdw-modifier", "rvdw-switch", "rvdw", "dispcorr", "fourierspacing",
           "pme-order", "ewald-rtol"]
_TC = ["tc-grps", "tau-t", "ref-t"]
_CONSTR = ["constraints", "constraint-algorithm"]
_BARO = ["tau-p", "compressibility", "ref-p", "pcoupltype"]

_NVT_PS, _NPT_PS = 100.0, 200.0


def _lam(v) -> str:
    return " ".join(f"{x:.4f}" for x in v)


def _copy(base: dict, keys) -> dict:
    return {k: base[k] for k in keys if k in base}


def fe_mdps(base: dict, s: Schedule, ligand: str, state: int, prod_ns: float,
            seed: int) -> dict[str, dict]:
    validate(s)
    if not 0 <= state < s.n_states:
        raise ValueError(f"state {state} out of range 0..{s.n_states - 1}")
    dt = float(base.get("dt", "0.002"))
    fe = {
        "free-energy": "yes",
        "couple-moltype": ligand,
        "couple-lambda0": "vdw-q",
        "couple-lambda1": "none",
        "couple-intramol": "no",
        "sc-alpha": "0.5",
        "sc-power": "1",
        "sc-r-power": "6",
        "nstdhdl": "100",
        "calc-lambda-neighbors": "-1",
        "dhdl-print-energy": "total",
        "separate-dhdl-file": "yes",
        "restraint-lambdas": _lam(s.restraint),
        "coul-lambdas": _lam(s.coul),
        "vdw-lambdas": _lam(s.vdw),
        "init-lambda-state": str(state),
    }

    def steps(ps):
        return str(int(round(ps / dt)))

    def md(integrator_nsteps, extra, barostat):
        d = {"integrator": "sd", "dt": base.get("dt", "0.002"),
             "nsteps": integrator_nsteps, "nstcomm": "100", "nstlog": "5000",
             "nstenergy": "5000"}
        d.update(_copy(base, _COMMON + _TC + _CONSTR))
        if barostat:
            d.update({"pcoupl": base.get("pcoupl", "C-rescale")})
            d.update(_copy(base, _BARO))
        else:
            d["pcoupl"] = "no"
        d.update(extra)
        d.update(fe)
        return d

    em = {"integrator": "steep", "nsteps": "5000", "emtol": "1000", "emstep": "0.001",
          "nstlog": "100", "nstenergy": "100"}
    em.update(_copy(base, _COMMON + _TC))
    em.update(fe)

    nvt = md(steps(_NVT_PS), {"gen-vel": "yes", "gen-temp": base.get("ref-t", "300").split()[0],
                              "gen-seed": str(seed), "ld-seed": str(seed),
                              "continuation": "no"}, False)
    npt = md(steps(_NPT_PS), {"continuation": "yes", "gen-vel": "no",
                              "ld-seed": str(seed)}, True)
    prod = md(str(int(round(prod_ns * 1000.0 / dt))),
              {"continuation": "yes", "gen-vel": "no", "ld-seed": str(seed),
               "nstxout-compressed": "5000"}, True)
    return {"em": em, "nvt": nvt, "npt": npt, "prod": prod}


def write_window_mdps(outdir, base: dict, s: Schedule, ligand: str, state: int,
                      prod_ns: float, seed: int) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, d in fe_mdps(base, s, ligand, state, prod_ns, seed).items():
        p = outdir / f"{name}.mdp"
        write_mdp(d, p)
        paths.append(p)
    return paths
