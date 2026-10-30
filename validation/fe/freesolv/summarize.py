#!/usr/bin/env python
"""Summarise the FreeSolv hydration validation (bcc vs abcg2).

usage: summarize.py [--work DIR] [--analyze] [--equil-ns N] [--out results.md] [--methods bcc abcg2]

Reads <work>/<RES>_<method>/ABFE/results/abfe_result.{csv,txt} (dG_hyd, MBAR = main estimate) and
molecules.csv (experiment). --analyze runs abfe_analyze.sh first on systems that have no results yet.
Prints per-molecule table, per-method MUE / RMSE / Pearson R / mean signed error, per-molecule
warnings, and the pass criteria (bcc MUE <= 1.5 kcal/mol; abcg2 MUE <= bcc MUE + 0.2). With --out
the same text is written as markdown. Missing systems are listed and left out of the statistics.
"""
import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def read_molecules(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_result(d):
    """(dG_hyd MBAR, err, dG TI, warnings) of one ABFE folder, or None if not analysed."""
    csv_p, txt_p = d / "ABFE" / "results" / "abfe_result.csv", d / "ABFE" / "results" / "abfe_result.txt"
    if not csv_p.is_file():
        return None
    mbar = ti = err = None
    with open(csv_p) as fh:
        for row in csv.DictReader(fh):
            if row["quantity"] == "dG_hyd" and row["repeat"] == "mean":
                if row["estimator"] == "MBAR":
                    mbar, err = float(row["value_kcal"]), float(row["error_kcal"])
                elif row["estimator"] == "TI":
                    ti = float(row["value_kcal"])
    if mbar is None:
        return None
    warns = []
    if txt_p.is_file():
        warns = [ln.split("WARNING:", 1)[1].strip() for ln in txt_p.read_text().splitlines()
                 if ln.strip().startswith("WARNING:")]
    return mbar, err, ti, warns


def stats(pairs):
    """pairs: [(calc, expt)] -> dict(n, mue, rmse, r, mse)."""
    n = len(pairs)
    if n == 0:
        return dict(n=0, mue=math.nan, rmse=math.nan, r=math.nan, mse=math.nan)
    diff = [c - e for c, e in pairs]
    mue = sum(abs(x) for x in diff) / n
    rmse = math.sqrt(sum(x * x for x in diff) / n)
    mse = sum(diff) / n
    r = math.nan
    if n > 2:
        mc, me = sum(c for c, _ in pairs) / n, sum(e for _, e in pairs) / n
        sc = math.sqrt(sum((c - mc) ** 2 for c, _ in pairs))
        se = math.sqrt(sum((e - me) ** 2 for _, e in pairs))
        if sc > 0 and se > 0:
            r = sum((c - mc) * (e - me) for c, e in pairs) / (sc * se)
    return dict(n=n, mue=mue, rmse=rmse, r=r, mse=mse)


def f(x, fmt="{:.2f}"):
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else fmt.format(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, default=Path("/cta/users/akocak/maint/abfe_test/freesolv"))
    ap.add_argument("--molecules", type=Path, default=HERE / "molecules.csv")
    ap.add_argument("--methods", nargs="+", default=["bcc", "abcg2"])
    ap.add_argument("--analyze", action="store_true", help="run abfe_analyze.sh where results are missing")
    ap.add_argument("--equil-ns", type=float, default=0.5)
    ap.add_argument("--out", type=Path, default=None, help="write the report as markdown")
    a = ap.parse_args()

    mols = read_molecules(a.molecules)
    data, missing = {}, []   # data[(res, method)] = result
    for m in mols:
        for meth in a.methods:
            d = a.work / f"{m['res']}_{meth}"
            if a.analyze and (d / "ABFE" / "setup.json").is_file() and read_result(d) is None:
                r = subprocess.run(["bash", str(REPO / "bin" / "abfe_analyze.sh"), str(d / "ABFE"),
                                    "--equil-ns", str(a.equil_ns)], capture_output=True, text=True)
                if r.returncode:
                    print(f"{m['res']}_{meth}: analysis failed: {r.stderr.strip().splitlines()[-1:]}", file=sys.stderr)
            res = read_result(d)
            if res is None:
                missing.append(f"{m['res']}_{meth}")
            else:
                data[(m["res"], meth)] = res

    L = ["# FreeSolv hydration free energies: bcc vs abcg2", "",
         "dG_hyd = -dG_solvent,decouple (MBAR, 1 repeat, error = statistical error of the single repeat); "
         "kcal/mol; experiment from FreeSolv database.json.", ""]
    head = "| res | molecule | FreeSolv id | expt |" + "".join(
        f" {m} calc | {m} err |" for m in a.methods)
    L += [head, "|" + "---|" * (4 + 2 * len(a.methods))]
    for m in mols:
        cells = []
        for meth in a.methods:
            r = data.get((m["res"], meth))
            if r is None:
                cells += ["n/a", "n/a"]
            else:
                e = r[0] - float(m["expt_kcal"])
                cells += [f"{r[0]:.2f} +- {r[1]:.2f}", f"{e:+.2f}"]
        L.append(f"| {m['res']} | {m['name']} | {m['freesolv_id']} | {float(m['expt_kcal']):.2f} +- "
                 f"{float(m['d_expt_kcal']):.2f} | " + " | ".join(cells) + " |")
    L += ["", "## Statistics vs experiment", "",
          "| charge model | n | MUE | RMSE | R | mean signed error |", "|---|---|---|---|---|---|"]
    st = {}
    for meth in a.methods:
        pairs = [(data[(m["res"], meth)][0], float(m["expt_kcal"])) for m in mols if (m["res"], meth) in data]
        st[meth] = stats(pairs)
        s = st[meth]
        L.append(f"| {meth} | {s['n']} | {f(s['mue'])} | {f(s['rmse'])} | {f(s['r'], '{:.3f}')} | {f(s['mse'], '{:+.2f}')} |")
    L += ["", "## Pass criteria", ""]
    if "bcc" in st and "abcg2" in st and st["bcc"]["n"] and st["abcg2"]["n"]:
        b, c = st["bcc"]["mue"], st["abcg2"]["mue"]
        full = st["bcc"]["n"] == st["abcg2"]["n"] == len(mols)
        L.append(f"- bcc MUE <= 1.5 kcal/mol: {b:.2f} -> {'PASS' if b <= 1.5 else 'FAIL'}")
        L.append(f"- abcg2 MUE <= bcc MUE + 0.2: {c:.2f} vs {b + 0.2:.2f} -> {'PASS' if c <= b + 0.2 else 'FAIL'}")
        if not full:
            L.append(f"- (provisional: only {st['bcc']['n']} bcc and {st['abcg2']['n']} abcg2 of {len(mols)} molecules done)")
    else:
        L.append("- not evaluable yet (need both methods)")
    L += ["", "## Warnings per system", ""]
    nw = 0
    for m in mols:
        for meth in a.methods:
            r = data.get((m["res"], meth))
            if r is None:
                continue
            ws = list(r[3])
            if r[2] is not None and abs(r[0] - r[2]) > 0.5 and not any("MBAR and TI" in w for w in ws):
                ws.append(f"MBAR and TI differ by {abs(r[0] - r[2]):.2f} kcal/mol")
            e = r[0] - float(m["expt_kcal"])
            if abs(e) > 2.0:
                ws.append(f"|calc - expt| = {abs(e):.2f} > 2 kcal/mol")
            for w in ws:
                L.append(f"- {m['res']}_{meth} ({m['name']}): {w}")
                nw += 1
    if not nw:
        L.append("- none")
    if missing:
        L += ["", "## Missing results", "", ", ".join(missing)]
    text = "\n".join(L) + "\n"
    print(text)
    if a.out:
        a.out.write_text(text)


if __name__ == "__main__":
    main()
