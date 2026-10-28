"""abfe_analyze: free energies from the per-window GROMACS dhdl files of an abfe_setup folder.

    abfe_analyze.sh ABFE_DIR [--equil-ns 0.5] [--temperature T] [--allow-missing]

Per leg and repeat (``<ABFE>/<leg>/rep<r>/lam_NN/prod.xvg``) it discards the first
``--equil-ns``, subsamples every window by its statistical inefficiency
(alchemlyb ``decorrelate_u_nk`` on ``dE`` / ``decorrelate_dhdl``), and runs MBAR
(main estimator) and TI (trapezoid; all dH/dlambda components summed, so the
Boresch restraint term in the bonded-lambda column is included). It also computes
the MBAR overlap matrix and a forward/backward convergence analysis (10 blocks).

The legs are combined with the analytic Boresch correction from boresch.json:

    dG_bind = dG_solvent,decouple - dG_complex,decouple - dG_restr,off
    dG_hyd  = -dG_solvent,decouple                       (no complex leg)

where dG_<leg>,decouple = G(last state, decoupled) - G(state 0, coupled) and
dG_restr,off (negative) is the free energy of releasing the restraint of the
decoupled ligand into V0. Results go to <ABFE>/results/: abfe_result.txt,
abfe_result.csv, overlap_<leg>_rep<r>.png, convergence_<leg>_rep<r>.png.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .restraints import KB, KCAL

# warning thresholds (Global Constraints / spec "Analysis"); kcal/mol except overlap
OVERLAP_MIN = 0.03
FWD_BWD_MAX = 0.5
MBAR_TI_MAX = 0.5
SPREAD_MAX = 1.0
MIN_SAMPLES = 10           # fewer after subsampling -> use the unsubsampled window
CONV_OVERLAP_FLOOR = 1e-4  # below this MBAR is meaningless: skip the 20-fit convergence
CONV_BLOCKS = 10
SHORT_TOL = 0.99           # window shorter than 99 % of prod_ns -> warning


class AnalysisError(RuntimeError):
    pass


@dataclass
class LegResult:
    """One leg x repeat; energies in kJ/mol."""
    dg_mbar: float
    err_mbar: float
    dg_ti: float
    err_ti: float
    min_overlap: float
    fwd_bwd_diff: float
    n_states: int
    notes: list = field(default_factory=list)
    overlap: object = None          # MBAR overlap matrix (np.ndarray) or None
    convergence: object = None      # alchemlyb forward/backward DataFrame (kT) or None
    missing: list = field(default_factory=list)


def log(msg: str) -> None:
    print(f"[abfe_analyze] {msg}", flush=True)


def _xvg_temperature(path) -> float | None:
    from alchemlyb.parsing.util import anyopen
    with anyopen(str(path), "r") as f:
        for line in f:
            if isinstance(line, bytes):
                line = line.decode()
            if not line.startswith(("#", "@")):
                return None
            m = re.search(r"T = ([0-9.eE+-]+)", line)
            if line.startswith("@ subtitle") and m:
                return float(m.group(1))
    return None


def _last_time(df) -> float:
    return float(df.index.get_level_values("time").max())


def _window_lambda(df):
    vals = df.index[0][1:]
    return vals[0] if len(vals) == 1 else tuple(vals)


def _subsample(df, func, label: str, notes: list, **kw):
    try:
        out = func(df, **kw)
    except Exception as e:  # pymbar timeseries on very short/constant series
        notes.append(f"{label}: subsampling failed ({type(e).__name__}: {e}); "
                     f"using all {len(df)} samples")
        return df
    if len(out) < MIN_SAMPLES:
        notes.append(f"{label}: subsampling left {len(out)} < {MIN_SAMPLES} samples; "
                     f"using all {len(df)} (correlated) samples")
        return df
    return out


def analyze_files(files, T: float, equil_ns: float, *, prod_ns: float | None = None,
                  label: str = "") -> LegResult:
    """MBAR + TI + overlap + forward/backward on the given window dhdl files."""
    import alchemlyb
    from alchemlyb.convergence import forward_backward_convergence
    from alchemlyb.estimators import MBAR, TI
    from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
    from alchemlyb.preprocessing import decorrelate_dhdl, decorrelate_u_nk, slicing

    files = [Path(f) for f in files]
    if not files:
        raise AnalysisError(f"{label}: no dhdl files")
    notes: list[str] = []
    pre = f"{label} " if label else ""
    kT = KB * T
    us, hs = [], []
    for f in files:
        name = f.parent.name if f.parent.name.startswith("lam_") else f.name
        tf = _xvg_temperature(f)
        if tf is not None and abs(tf - T) > 0.5:
            raise AnalysisError(f"{f}: dhdl file says T = {tf:g} K but the analysis uses "
                                f"{T:g} K (setup.json / --temperature)")
        try:
            u = extract_u_nk(str(f), T=T)
            h = extract_dHdl(str(f), T=T)
        except Exception as e:
            raise AnalysisError(f"{f}: cannot parse ({type(e).__name__}: {e}); "
                                "truncated or not a GROMACS dhdl file?") from e
        tend = _last_time(u)
        if prod_ns is not None and tend < SHORT_TOL * prod_ns * 1000.0:
            notes.append(f"{pre}{name}: short window ({tend / 1000:g} of {prod_ns:g} ns)")
        if equil_ns > 0:
            u = slicing(u, lower=equil_ns * 1000.0)
            h = slicing(h, lower=equil_ns * 1000.0)
        if len(u) == 0:
            raise AnalysisError(f"{f}: no data after discarding --equil-ns {equil_ns:g} "
                                f"(window ends at {tend / 1000:g} ns)")
        us.append(_subsample(u, decorrelate_u_nk, f"{pre}{name} u_nk", notes, method="dE"))
        hs.append(_subsample(h, decorrelate_dhdl, f"{pre}{name} dH/dl", notes))

    columns = list(us[0].columns)
    n_states = len(columns)
    sampled = sorted(columns.index(_window_lambda(u)) for u in us)

    dg_ti = err_ti = math.nan
    try:
        ti = TI().fit(alchemlyb.concat(hs))
        dg_ti = float(ti.delta_f_.iloc[0, -1]) * kT
        err_ti = float(ti.d_delta_f_.iloc[0, -1]) * kT
    except Exception as e:
        notes.append(f"{pre}TI failed ({type(e).__name__}: {e})")

    dg_mbar = err_mbar = min_ov = math.nan
    overlap = None
    try:
        mbar = MBAR().fit(alchemlyb.concat(us))
        dg_mbar = float(mbar.delta_f_.iloc[0, -1]) * kT
        err_mbar = float(mbar.d_delta_f_.iloc[0, -1]) * kT
        overlap = np.asarray(mbar.overlap_matrix)
        pairs = [overlap[a, b] for a, b in zip(sampled[:-1], sampled[1:])]
        min_ov = float(min(pairs)) if pairs else math.nan
    except Exception as e:
        notes.append(f"{pre}MBAR failed ({type(e).__name__}: {e})")

    conv = None
    fwd_bwd = math.nan
    if not (min_ov >= CONV_OVERLAP_FLOOR):
        notes.append(f"{pre}forward/backward convergence skipped (MBAR neighbour overlap "
                     f"{min_ov:.2g} < {CONV_OVERLAP_FLOOR:g})")
    else:
        try:
            conv = forward_backward_convergence(us, "MBAR", num=CONV_BLOCKS)
            half = conv.iloc[CONV_BLOCKS // 2 - 1]  # 50 % of the data from each end
            fwd_bwd = abs(float(half["Forward"]) - float(half["Backward"])) * kT
        except Exception as e:
            notes.append(f"{pre}forward/backward convergence failed "
                         f"({type(e).__name__}: {e})")
            conv = None
    return LegResult(dg_mbar, err_mbar, dg_ti, err_ti, min_ov, fwd_bwd, n_states,
                     notes=notes, overlap=overlap, convergence=conv)


def _window_file(rep_dir: Path, i: int) -> Path:
    return Path(rep_dir) / f"lam_{i:02d}" / "prod.xvg"


def missing_windows(rep_dir, n_states: int) -> list[int]:
    """State indices whose lam_NN/prod.xvg is absent or empty."""
    out = []
    for i in range(n_states):
        f = _window_file(rep_dir, i)
        if not f.is_file() or f.stat().st_size == 0:
            out.append(i)
    return out


def analyze_leg(rep_dir, T: float, equil_ns: float, *, n_states: int | None = None,
                prod_ns: float | None = None, allow_missing: bool = False,
                label: str = "") -> LegResult:
    rep_dir = Path(rep_dir)
    if n_states is None:
        n_states = len(sorted(rep_dir.glob("lam_*/")))
    miss = missing_windows(rep_dir, n_states)
    if miss and not allow_missing:
        raise AnalysisError(f"{label or rep_dir}: missing windows {miss}")
    files = [_window_file(rep_dir, i) for i in range(n_states) if i not in miss]
    res = analyze_files(files, T, equil_ns, prod_ns=prod_ns, label=label)
    if res.n_states != n_states:
        raise AnalysisError(f"{label or rep_dir}: dhdl files have {res.n_states} states, "
                            f"expected {n_states}")
    res.missing = miss
    if miss:
        res.notes.insert(0, f"{label}: missing windows {miss} analysed without them "
                            "(--allow-missing)" + ("; TI end point missing, TI is wrong"
                                                    if 0 in miss or n_states - 1 in miss
                                                    else ""))
    return res


def combine(setup: dict, legs: dict, boresch_dg) -> dict:
    """Combine legs per repeat; returns kcal/mol values (and kJ/mol copies).

    legs: {(leg, repeat): LegResult} (kJ/mol); boresch_dg: dG_restr,off in kJ/mol
    (negative) or None for hydration.
    """
    binding = "complex" in setup["legs"]
    reps = sorted({r for (_, r) in legs})
    per, per_kj, skipped = {}, {}, []
    for r in reps:
        s = legs.get(("solvent", r))
        c = legs.get(("complex", r))
        if s is None or (binding and c is None):
            skipped.append(r)
            continue
        per[r], per_kj[r] = {}, {}
        for est, (v, e) in {"MBAR": ("dg_mbar", "err_mbar"), "TI": ("dg_ti", "err_ti")}.items():
            if binding:
                val = getattr(s, v) - getattr(c, v) - boresch_dg
                err = math.hypot(getattr(s, e), getattr(c, e))
            else:
                val, err = -getattr(s, v), getattr(s, e)
            per_kj[r][est] = (val, err)
            per[r][est] = (val / KCAL, err / KCAL)
    mean = {}
    for est in ("MBAR", "TI"):
        vals = np.array([per[r][est][0] for r in per])
        if len(vals) == 0:
            mean[est] = (math.nan, math.nan)
        elif len(vals) == 1:
            mean[est] = per[next(iter(per))][est]
        else:
            mean[est] = (float(vals.mean()), float(vals.std(ddof=1) / math.sqrt(len(vals))))
    return {
        "quantity": "dG_bind" if binding else "dG_hyd",
        "formula": ("dG_bind = dG_solvent,decouple - dG_complex,decouple - dG_restr,off"
                    if binding else "dG_hyd = -dG_solvent,decouple"),
        "restr_kJ": boresch_dg if binding else None,
        "per_repeat": per, "per_repeat_kJ": per_kj, "mean": mean,
        "n_repeats": len(per),
        "error_kind": "SEM over repeats" if len(per) > 1 else
                      "statistical error of the single repeat",
        "skipped_repeats": skipped,
    }


def warnings(legs: dict, combined: dict | None = None) -> list[str]:
    """Warnings with the spec thresholds (they never fail the analysis)."""
    out = []
    for (leg, r), res in sorted(legs.items()):
        tag = f"{leg} rep{r}"
        if not math.isfinite(res.dg_mbar):
            out.append(f"{tag}: MBAR estimate not available")
        if not math.isfinite(res.min_overlap) or res.min_overlap < OVERLAP_MIN:
            out.append(f"{tag}: minimum neighbour overlap {res.min_overlap:.3g} < "
                       f"{OVERLAP_MIN} (add lambda states or sample longer)")
        if math.isfinite(res.fwd_bwd_diff) and res.fwd_bwd_diff / KCAL > FWD_BWD_MAX:
            out.append(f"{tag}: forward/backward (50 % of the data each) differ by "
                       f"{res.fwd_bwd_diff / KCAL:.2f} > {FWD_BWD_MAX} kcal/mol "
                       "(not converged)")
        d = abs(res.dg_mbar - res.dg_ti) / KCAL
        if math.isfinite(d) and d > MBAR_TI_MAX:
            out.append(f"{tag}: MBAR and TI differ by {d:.2f} > {MBAR_TI_MAX} kcal/mol")
        out.extend(f"{tag}: {n}" if not n.startswith(leg) else n for n in res.notes)
    if combined:
        q = combined["quantity"]
        for r in combined["skipped_repeats"]:
            out.append(f"rep{r}: a leg is missing, repeat left out of {q}")
        vals = [v["MBAR"][0] for v in combined["per_repeat"].values()]
        vals = [v for v in vals if math.isfinite(v)]
        if len(vals) > 1 and max(vals) - min(vals) > SPREAD_MAX:
            out.append(f"{q} (MBAR) repeat spread {max(vals) - min(vals):.2f} > "
                       f"{SPREAD_MAX} kcal/mol")
    return out


# ---------------------------------------------------------------- output files

def _placeholder(path: Path, title: str, msg: str) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.5, 0.5, msg, ha="center", va="center", wrap=True, fontsize=9)
    fig.savefig(path, dpi=100)
    plt.close(fig)


def _plots(res: LegResult, leg: str, r: int, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from alchemlyb.visualisation import plot_convergence, plot_mbar_overlap_matrix

    p = outdir / f"overlap_{leg}_rep{r}.png"
    if res.overlap is not None:
        try:
            ax = plot_mbar_overlap_matrix(res.overlap)
            ax.set_title(f"MBAR overlap {leg} rep{r}")
            ax.figure.savefig(p, dpi=100, bbox_inches="tight")
            plt.close(ax.figure)
        except Exception as e:
            _placeholder(p, f"overlap {leg} rep{r}", f"plot failed: {e}")
    else:
        _placeholder(p, f"overlap {leg} rep{r}", "MBAR failed; no overlap matrix")
    p = outdir / f"convergence_{leg}_rep{r}.png"
    if res.convergence is not None:
        try:
            ax = plot_convergence(res.convergence, units="kcal/mol")
            ax.set_title(f"forward/backward {leg} rep{r}")
            ax.figure.savefig(p, dpi=100, bbox_inches="tight")
            plt.close(ax.figure)
        except Exception as e:
            _placeholder(p, f"convergence {leg} rep{r}", f"plot failed: {e}")
    else:
        why = next((n for n in res.notes if "convergence" in n), "not computed")
        _placeholder(p, f"convergence {leg} rep{r}", why)


def _f(v: float, e: float | None = None, scale: float = 1.0) -> str:
    if not math.isfinite(v):
        return "nan"
    s = f"{v / scale:9.3f}"
    if e is not None:
        s += f" +- {e / scale:.3f}" if math.isfinite(e) else " +- nan"
    return s


def write_results(outdir: Path, setup: dict, legs: dict, comb: dict, warns: list,
                  *, T: float, equil_ns: float, boresch: dict | None, abfe: Path,
                  args_line: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    q = comb["quantity"]
    binding = q == "dG_bind"
    L = []
    L.append(f"ABFE analysis ({args_line}; {datetime.datetime.now():%Y-%m-%d %H:%M})")
    L.append(f"folder: {abfe}")
    L.append(f"ligand: {setup.get('ligand', '?')} (net charge {setup.get('ligand_charge', '?')})"
             f"   T = {T:g} K   equilibration discarded: {equil_ns:g} ns per window")
    L.append("estimators: alchemlyb MBAR (main) and TI (trapezoid, all dH/dl components); "
             "windows subsampled by statistical inefficiency")
    L.append("legs: " + ", ".join(f"{l} ({setup['n_states'][l]} states)"
                                  for l in setup["legs"])
             + f"; repeats analysed: {sorted({r for _, r in legs})}")
    L.append("")
    L.append("Sign convention")
    L.append("  dG_<leg>,decouple = G(last lambda state, ligand decoupled) - G(state 0, coupled)")
    if binding:
        L.append("  dG_restr,off = analytic free energy (Boresch 2003, eq. 32) of releasing the "
                 "restraint on the decoupled ligand into V0 = 1.66054 nm^3; negative")
        L.append(f"  {comb['formula']}")
        L.append("  (more negative dG_bind = stronger binding)")
    else:
        L.append(f"  {comb['formula']}   (hydration free energy; no complex leg)")
    L.append("")
    if binding:
        rk = comb["restr_kJ"]
        L.append(f"Restraint correction (boresch.json, T = {boresch.get('T', T):g} K)")
        L.append(f"  dG_restr,off = {rk:.3f} kJ/mol = {rk / KCAL:.3f} kcal/mol"
                 f"   (enters dG_bind as -dG_restr,off = {-rk / KCAL:+.3f} kcal/mol)")
        L.append("")
    L.append("Per leg and repeat, dG_decouple in kcal/mol (+- statistical error)")
    L.append(f"  {'leg':8s} {'rep':>3s}  {'MBAR':>18s}  {'TI':>18s}  {'min overlap':>11s}"
             f"  {'|fwd-bwd| 50%':>13s}  states")
    for (leg, r), res in sorted(legs.items()):
        fb = f"{res.fwd_bwd_diff / KCAL:13.3f}" if math.isfinite(res.fwd_bwd_diff) else \
             f"{'n/a':>13s}"
        L.append(f"  {leg:8s} {r:3d}  {_f(res.dg_mbar, res.err_mbar, KCAL):>18s}  "
                 f"{_f(res.dg_ti, res.err_ti, KCAL):>18s}  {res.min_overlap:11.3g}  {fb}  "
                 f"{res.n_states}" + (f" (missing {res.missing})" if res.missing else ""))
    L.append("")
    L.append(f"{q} per repeat")
    for r, d in comb["per_repeat"].items():
        kj = comb["per_repeat_kJ"][r]
        L.append(f"  rep{r}: MBAR {_f(*d['MBAR'])} kcal/mol ({_f(*kj['MBAR'])} kJ/mol)"
                 f"   TI {_f(*d['TI'])} kcal/mol")
    L.append("")
    L.append(f"Result ({comb['n_repeats']} repeat(s); error = {comb['error_kind']})")
    for est in ("MBAR", "TI"):
        v, e = comb["mean"][est]
        L.append(f"  {q} ({est:4s}) = {_f(v, e)} kcal/mol   = {_f(v, e, 1 / KCAL)} kJ/mol")
    L.append("")
    L.append(f"Warnings ({len(warns)})")
    L.extend(f"  WARNING: {w}" for w in warns)
    if not warns:
        L.append("  none")
    (outdir / "abfe_result.txt").write_text("\n".join(L) + "\n")

    with (outdir / "abfe_result.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "repeat", "estimator", "value_kcal", "error_kcal"])
        for (leg, r), res in sorted(legs.items()):
            w.writerow([f"dG_{leg}_decouple", r, "MBAR", f"{res.dg_mbar / KCAL:.4f}",
                        f"{res.err_mbar / KCAL:.4f}"])
            w.writerow([f"dG_{leg}_decouple", r, "TI", f"{res.dg_ti / KCAL:.4f}",
                        f"{res.err_ti / KCAL:.4f}"])
        if binding:
            w.writerow(["dG_restr_off", "all", "analytic", f"{comb['restr_kJ'] / KCAL:.4f}",
                        "0.0000"])
        for r, d in comb["per_repeat"].items():
            for est in ("MBAR", "TI"):
                w.writerow([q, r, est, f"{d[est][0]:.4f}", f"{d[est][1]:.4f}"])
        for est in ("MBAR", "TI"):
            v, e = comb["mean"][est]
            w.writerow([q, "mean", est, f"{v:.4f}", f"{e:.4f}"])


# ---------------------------------------------------------------- CLI

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="abfe_analyze.sh", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("abfe", type=Path, help="folder made by abfe_setup.sh")
    p.add_argument("--equil-ns", type=float, default=0.5,
                   help="discard the first N ns of every production window (default 0.5)")
    p.add_argument("--temperature", type=float, default=None,
                   help="K; default: setup.json (the effective run temperature)")
    p.add_argument("--allow-missing", action="store_true",
                   help="analyse legs with missing windows instead of refusing")
    p.add_argument("--outdir", type=Path, default=None, help="default ABFE/results")
    return p


def _repeats(setup: dict, leg_dir: Path) -> list[int]:
    n = (setup.get("options") or {}).get("repeats")
    if n:
        return list(range(int(n)))
    return sorted(int(p.name[3:]) for p in leg_dir.glob("rep*") if p.name[3:].isdigit())


def run(args) -> Path:
    abfe = args.abfe.resolve()
    sj = abfe / "setup.json"
    if not sj.is_file():
        raise AnalysisError(f"{abfe}: no setup.json (not an abfe_setup folder)")
    setup = json.loads(sj.read_text())
    if args.equil_ns < 0:
        raise AnalysisError("--equil-ns must be >= 0")
    legs_present = setup["legs"]
    if "solvent" not in legs_present:
        raise AnalysisError("setup.json has no solvent leg")
    temps = {leg: float(args.temperature if args.temperature is not None
                        else setup["temperature"][leg]) for leg in legs_present}
    T = temps["solvent"]
    boresch, restr = None, None
    pre_warn = []
    if "complex" in legs_present:
        bj = abfe / "boresch.json"
        if not bj.is_file():
            raise AnalysisError(f"{bj} missing (complex leg needs the restraint correction)")
        boresch = json.loads(bj.read_text())
        restr = float(boresch["dG_restr_off_kJ"])
        if abs(float(boresch.get("T", temps["complex"])) - temps["complex"]) > 0.5:
            pre_warn.append(f"boresch.json correction is for T = {boresch['T']:g} K but the "
                            f"analysis uses {temps['complex']:g} K")

    # refuse missing windows before any heavy work
    plan, missing = [], []
    for leg in legs_present:
        n = int(setup["n_states"][leg])
        for r in _repeats(setup, abfe / leg):
            m = missing_windows(abfe / leg / f"rep{r}", n)
            if m:
                missing.append(f"{leg} rep{r}: missing windows {m}")
            if len(m) < n:
                plan.append((leg, r, n))
            else:
                missing[-1] += " (no data at all)"
    if missing and not args.allow_missing:
        raise AnalysisError("; ".join(missing)
                            + " -- rerun those windows or use --allow-missing")

    outdir = (args.outdir or abfe / "results").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    legs = {}
    for leg, r, n in plan:
        log(f"{leg} rep{r}: {n} states, T = {temps[leg]:g} K")
        res = analyze_leg(abfe / leg / f"rep{r}", temps[leg], args.equil_ns, n_states=n,
                          prod_ns=(setup.get("prod_ns") or {}).get(leg),
                          allow_missing=args.allow_missing, label=f"{leg} rep{r}")
        legs[(leg, r)] = res
        _plots(res, leg, r, outdir)
    comb = combine(setup, legs, restr)
    warns = pre_warn + [m for m in missing if "no data" in m] + warnings(legs, comb)
    write_results(outdir, setup, legs, comb, warns, T=T, equil_ns=args.equil_ns,
                  boresch=boresch, abfe=abfe,
                  args_line="abfe_analyze.sh " + " ".join(sys.argv[1:]
                                                          if args._argv is None
                                                          else args._argv))
    for est in ("MBAR", "TI"):
        v, e = comb["mean"][est]
        log(f"{comb['quantity']} ({est}) = {_f(v, e).strip()} kcal/mol")
    for w in warns:
        log(f"WARNING: {w}")
    log(f"results: {outdir / 'abfe_result.txt'}")
    return outdir


def _quiet_alchemlyb() -> None:
    """alchemlyb logs every subsampling step at DEBUG level through loguru."""
    try:
        from loguru import logger
        logger.remove()
        logger.add(sys.stderr, level="WARNING")
    except ImportError:
        pass


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    _quiet_alchemlyb()
    args._argv = argv
    try:
        run(args)
    except (AnalysisError, KeyError, FileNotFoundError, ValueError) as e:
        print(f"abfe_analyze: error: {e}", file=sys.stderr)
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
