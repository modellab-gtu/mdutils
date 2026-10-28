"""abfe_setup: build the ABFE (or hydration) folder from mdutils MD run folders.

    abfe_setup.sh --complex DIR --solvent DIR --ligand LIG [options]
    abfe_setup.sh --solvent DIR --ligand LIG            (hydration free energy)

Per leg it writes topol_fe.top (define-free, co-alchemical ions, Boresch
restraint for the complex), start.gro (= md_<r>.gro, ions moved to the end),
posre_ref.gro, index.ndx (rebuilt for start.gro), rep<r>/lam_NN/{em,nvt,npt,
prod}.mdp and rep<r>/run_abfe.sbatch; plus setup.json, boresch.json and
submit_abfe.sh at the top. Nothing is written into the run folders.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

from .mdp import write_window_mdps
from .rundir import RunDir, load_rundir
from .schedule import SCHEDULE_DIR, Schedule, load_schedule, validate
from .slurm import TRUBA_DEFAULTS, write_leg_sbatch, write_submit
from .topology import Topology, ligand_charge, preprocess

REPO_ROOT = Path(__file__).resolve().parents[3]
LEGS = ("complex", "solvent")


class SetupError(RuntimeError):
    pass


def log(msg: str) -> None:
    print(f"[abfe_setup] {msg}", flush=True)


# --------------------------------------------------------------------------- CLI
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="abfe_setup.sh",
        description="Set up GROMACS ABFE windows (complex + solvent legs) from mdutils MD "
                    "run folders; with only --solvent, a hydration free energy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--complex", type=Path, help="protein-ligand MD run folder (prepare_prolig_md.sh)")
    p.add_argument("--solvent", type=Path, help="ligand-in-water MD run folder (prepare_lig_md.sh)")
    p.add_argument("--ligand", required=True,
                   help="ligand [ moleculetype ] name (= residue name), e.g. JZ4")
    p.add_argument("--replica", type=int, default=0,
                   help="MD replica r whose md_<r>.gro/tpr/xtc are used")
    p.add_argument("--repeats", type=int, default=3, help="independent repeats per leg")
    p.add_argument("--complex-ns", type=float, default=5.0, help="production ns per complex window")
    p.add_argument("--solvent-ns", type=float, default=3.0, help="production ns per solvent window")
    p.add_argument("--schedule", default=None,
                   help="lambda schedule YAML, or a name in automd/fe/abfe/schedules "
                        "(e.g. test_short); default: the built-in 36/31-state schedule")
    p.add_argument("--cluster", choices=("gtu", "truba"), default="gtu")
    p.add_argument("--complex-node", default="cws01", help="gtu node for the complex leg")
    p.add_argument("--solvent-node", default="cws02", help="gtu node for the solvent leg")
    p.add_argument("--truba-account", default=os.environ.get("USER", ""))
    p.add_argument("--truba-partition", default=TRUBA_DEFAULTS["partition"])
    p.add_argument("--truba-ntasks", type=int, default=TRUBA_DEFAULTS["ntasks"],
                   help="MPI ranks per window")
    p.add_argument("--truba-omp", type=int, default=TRUBA_DEFAULTS["omp"],
                   help="OpenMP threads per rank")
    p.add_argument("--truba-gmx-module", default=TRUBA_DEFAULTS["gmx_module"])
    p.add_argument("--truba-time", default=TRUBA_DEFAULTS["time"])
    p.add_argument("--outdir", type=Path, default=Path("ABFE"), help="output folder (new)")
    p.add_argument("--overwrite", action="store_true",
                   help="replace an existing ABFE output folder")
    p.add_argument("--force-restraint", action="store_true",
                   help="accept a Boresch anchor angle within 10 deg of 0/180")
    p.add_argument("--temperature", type=float, default=None,
                   help="run every FE stage of both legs at T (K): overrides ref-t and "
                        "gen-temp; required when the legs' md.mdp temperatures differ "
                        "by > 0.5 K")
    p.add_argument("--submit", action="store_true",
                   help="run submit_abfe.sh at the end (gtu only)")
    return p


def resolve_paths(args) -> None:
    """Make every user-supplied path absolute (relative to the invocation cwd) before
    anything changes directory or spawns a subprocess (grompp runs elsewhere).
    A --schedule that is not an existing file is a schedule NAME and stays as is."""
    for k in ("complex", "solvent", "outdir"):
        v = getattr(args, k)
        if v is not None:
            setattr(args, k, Path(v).expanduser().resolve())
    if args.schedule is not None:
        p = Path(args.schedule).expanduser()
        if p.is_file():
            args.schedule = str(p.resolve())


def resolve_schedule(spec: str | None) -> Path:
    if spec is None:
        return SCHEDULE_DIR / "default.yaml"
    p = Path(spec).expanduser()
    if p.is_file():
        return p.resolve()
    named = SCHEDULE_DIR / f"{spec}.yaml"
    if named.is_file():
        return named
    avail = ", ".join(sorted(x.stem for x in SCHEDULE_DIR.glob("*.yaml")))
    raise SetupError(f"--schedule {spec}: no such file or built-in schedule ({avail})")


def _inside(child: Path, parent: Path) -> bool:
    child, parent = child.resolve(), parent.resolve()
    return child == parent or parent in child.parents


# ------------------------------------------------------------------- helpers
def _gmx() -> str:
    exe = shutil.which("gmx")
    if exe is None:
        raise SetupError("gmx not found: module load gromacs/2026.1-tmpi")
    return exe


def _run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True,
                       env={**os.environ, "GMX_MAXBACKUP": "-1"}, **kw)
    if r.returncode != 0:
        raise SetupError(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-2000:]}")
    return r


def _ndx_groups(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    cur = None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("["):
            cur = s.strip("[] \t")
            groups[cur] = []
        elif s and cur is not None:
            groups[cur].extend(int(x) for x in s.split())
    return groups


def build_index(gro: Path, out: Path, tcgrps: list[str], lig_res: str,
                must_contain: dict[str, list[int]] | None = None) -> dict[str, list[int]]:
    """index.ndx for ``gro`` with every tc group (make_ndx defaults + named groups).

    ``must_contain`` maps a group name to 1-based atom numbers it must contain
    (the co-alchemical ions in Water_and_ions).
    """
    gmx = _gmx()
    _run([gmx, "make_ndx", "-f", str(gro), "-o", str(out)], input="q\n")
    for g in tcgrps:
        groups = _ndx_groups(out)
        if g in groups:
            continue
        if g == lig_res:
            sel = f"r {lig_res}"
        elif g == f"Protein_{lig_res}" and "Protein" in groups:
            sel = f'"Protein" | r {lig_res}'
        elif g == "Water_and_ions":
            w = "Water" if "Water" in groups else "SOL"
            sel = f'"{w}" | "Ion"' if "Ion" in groups else f'"{w}"'
        else:
            raise SetupError(f"tc group [ {g} ] (from MDP/md.mdp) cannot be rebuilt for "
                             f"{gro.name}; known: {', '.join(groups)}")
        n = len(groups)
        _run([gmx, "make_ndx", "-f", str(gro), "-n", str(out), "-o", str(out)],
             input=f"{sel}\nname {n} {g}\nq\n")
    groups = _ndx_groups(out)
    for g in tcgrps:
        if not groups.get(g):
            raise SetupError(f"tc group [ {g} ] missing or empty in {out}")
    for g, atoms in (must_contain or {}).items():
        miss = sorted(set(atoms) - set(groups.get(g, [])))
        if miss:
            raise SetupError(f"index group [ {g} ] lacks atoms {miss}")
    return groups


def _gro_atoms(gro: Path):
    """(resid, resname, atomname) per atom; the gro resid field is 5 wide (wraps)."""
    lines = gro.read_text().splitlines()
    n = int(lines[1])
    return [(int(l[0:5]), l[5:10].strip(), l[10:15].strip()) for l in lines[2:2 + n]]


def check_boresch_names(b, gro: Path) -> None:
    """The six anchors' names (``RESNAME+RESID:ATOM`` from md_<r>.tpr) must match
    residue name, residue number (modulo 100000, the gro field width) and atom
    name of start.gro at b.atoms."""
    atoms = _gro_atoms(gro)
    for idx, name in zip(b.atoms, b.names):
        resid_part, aname = name.rsplit(":", 1)
        if idx >= len(atoms):
            raise SetupError(f"Boresch atom {idx} ({name}) beyond {gro} ({len(atoms)} atoms)")
        resid, res, at = atoms[idx]
        tail = resid_part[len(res):] if resid_part.startswith(res) else ""
        ok = (at == aname[:5] and tail.lstrip("-").isdigit()
              and int(tail) % 100000 == resid % 100000
              and len(resid_part) - len(tail) <= 5)
        if not ok:
            raise SetupError(f"Boresch atom {idx}: boresch.json says {name}, "
                             f"{gro.name} has {res}{resid}:{at}")


def _ligand_signature(top: Topology, ligand: str):
    """Per atom (name, type, mass, charge) of the ligand moleculetype."""
    mt = top.require_moleculetype(ligand)
    return [(a.name, a.type, a.mass, a.q) for a in mt.atoms]


def compare_ligands(ligand: str, a: list, b: list, qtol: float = 1e-3) -> None:
    """The ligand must be the same molecule in both legs; raise naming the atoms."""
    if len(a) != len(b):
        raise SetupError(f"ligand {ligand}: {len(a)} atoms in --complex, {len(b)} in "
                         "--solvent: not the same ligand")
    bad = []
    for i, (x, y) in enumerate(zip(a, b), 1):
        why = []
        if x[0] != y[0]:
            why.append(f"name {x[0]}/{y[0]}")
        if x[1] != y[1]:
            why.append(f"type {x[1]}/{y[1]}")
        if abs(x[2] - y[2]) > 1e-3:
            why.append(f"mass {x[2]:g}/{y[2]:g}")
        if abs(x[3] - y[3]) > qtol:
            why.append(f"charge {x[3]:+.4f}/{y[3]:+.4f}")
        if why:
            bad.append(f"#{i} {x[0]} ({', '.join(why)})")
    if bad:
        raise SetupError(f"ligand {ligand} differs between --complex and --solvent "
                         f"(complex/solvent): {'; '.join(bad[:10])}"
                         + (f"; ... {len(bad) - 10} more" if len(bad) > 10 else ""))


def _versions() -> dict:
    from importlib import metadata
    out = {"python": sys.version.split()[0]}
    for pkg in ("MDAnalysis", "MDRestraintsGenerator", "alchemlyb", "pymbar", "numpy",
                "PyYAML"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    try:
        r = subprocess.run([_gmx(), "--version"], capture_output=True, text=True)
        out["gromacs"] = next((l.split(":", 1)[1].strip() for l in r.stdout.splitlines()
                               if l.startswith("GROMACS version")), None)
    except SetupError:
        out["gromacs"] = None
    return out


def _git_describe() -> str:
    r = subprocess.run(["git", "-C", str(REPO_ROOT), "describe", "--always", "--dirty",
                        "--tags"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else "unknown"


def seed_for(repeat: int, state: int) -> int:
    return 1000 * repeat + state + 1


# ---------------------------------------------------------------------- legs
def build_leg(leg: str, run: RunDir, leg_dir: Path, outdir: Path, args, sched: Schedule,
              meta: dict, warnings: list[str]) -> None:
    from .ions import apply_coalchemical, check_charges, choose_ions

    import MDAnalysis as mda

    leg_dir.mkdir(parents=True)
    lig = args.ligand
    with tempfile.TemporaryDirectory(dir=outdir, prefix=f".pp_{leg}_") as td:
        top = Topology.parse(preprocess(run, Path(td)))
    mt = top.require_moleculetype(lig)
    n_lig = sum(c for n, c in top.molecules if n == lig)
    if n_lig != 1:
        raise SetupError(f"{leg}: ligand {lig} appears {n_lig} times in [ molecules ] (need 1)")
    lig_res = mt.atoms[0].res
    q = ligand_charge(top, lig)
    T = args.temperature if args.temperature is not None else run.temperature
    base = dict(run.md_mdp)
    if args.temperature is not None:
        n_tc = len(base["ref-t"].split())
        base["ref-t"] = " ".join([f"{T:g}"] * n_tc)   # gen-temp follows ref-t
    start, posre = leg_dir / "start.gro", leg_dir / "posre_ref.gro"

    ions = []
    if q != 0:
        u = mda.Universe(str(run.md_gro))
        idx = choose_ions(u, top, lig, q)
        apply_coalchemical(top, run.md_gro, start, idx, lig)
        natoms = len(u.atoms)
        for j, k in enumerate(idx):
            a = u.atoms[k]
            ions.append({"md_index0": int(k), "start_index0": natoms - len(idx) + j,
                         "resname": a.resname, "name": a.name,
                         "moleculetype": f"{top_name_at(top, natoms - len(idx) + j)}"})
        log(f"{leg}: ligand charge {q:+d}; co-alchemical ion(s) "
            + ", ".join(f"{i['name']} md#{i['md_index0']} -> start#{i['start_index0']}"
                        for i in ions))
    else:
        shutil.copyfile(run.md_gro, start)
        shutil.copyfile(run.md_gro, posre)
        check_charges(top, lig)
    if not posre.is_file():
        raise SetupError(f"{posre} was not written")

    if leg == "complex":
        from .restraints import check_collinear, find_boresch, to_json
        log("complex: searching Boresch anchors (MDRestraintsGenerator) on "
            f"{run.md_tpr.name}/{run.md_xtc.name}")
        # MDAnalysis writes .<xtc>_offsets.npz/.lock next to the trajectory it
        # opens: read symlinks in a scratch folder so the run folder stays untouched
        with tempfile.TemporaryDirectory(dir=outdir, prefix=".mdrg_") as td:
            tpr, xtc = Path(td) / run.md_tpr.name, Path(td) / run.md_xtc.name
            tpr.symlink_to(run.md_tpr.resolve())
            xtc.symlink_to(run.md_xtc.resolve())
            b_md = find_boresch(tpr, xtc, lig_res, T, force=args.force_restraint)
        if args.force_restraint:
            try:
                check_collinear(b_md)
            except ValueError as e:
                warnings.append(f"--force-restraint: {e}")
        moved = sorted(i["md_index0"] for i in ions)
        b = replace(b_md, atoms=tuple(a - sum(m < a for m in moved) for a in b_md.atoms))
        check_boresch_names(b, start)
        from .restraints import intermolecular_block
        top.append_section(intermolecular_block(b))
        bj = to_json(b, T, outdir / "boresch.json")
        bj["atoms_0based_md_tpr"] = list(b_md.atoms)
        bj["note"] = ("atoms_0based/atoms_1based follow complex/start.gro and topol_fe.top; "
                      "atoms_0based_md_tpr follow md_<r>.tpr")
        (outdir / "boresch.json").write_text(json.dumps(bj, indent=2) + "\n")
        meta["boresch"] = {"atoms_1based": bj["atoms_1based"], "names": bj["names"],
                           "dG_restr_off_kJ": bj["dG_restr_off_kJ"],
                           "dG_restr_off_kcal": bj["dG_restr_off_kcal"]}
        log(f"complex: Boresch {', '.join(b.names)}; dG_restr,off = "
            f"{bj['dG_restr_off_kcal']:.3f} kcal/mol")

    top.write(leg_dir / "topol_fe.top")

    tc = run.md_mdp.get("tc-grps", "").split()
    if not tc:
        raise SetupError(f"{leg}: no tc-grps in {run.path / 'MDP' / 'md.mdp'}")
    natoms = len(_gro_atoms(start))
    must = {"Water_and_ions": [i["start_index0"] + 1 for i in ions]} \
        if ions and "Water_and_ions" in tc else None
    if ions and must is None:
        raise SetupError(f"{leg}: tc-grps {tc} has no Water_and_ions group for the "
                         "co-alchemical ion(s)")
    build_index(start, leg_dir / "index.ndx", tc, lig_res, must)

    prod_ns = args.complex_ns if leg == "complex" else args.solvent_ns
    seeds = {}
    for r in range(args.repeats):
        for s in range(sched.n_states):
            w = leg_dir / f"rep{r}" / f"lam_{s:02d}"
            sd = seed_for(r, s)
            write_window_mdps(w, base, sched, lig, s, prod_ns, sd)
            seeds[f"rep{r}/lam_{s:02d}"] = sd
    meta["ligand_charge"][leg] = q
    meta["alch_ions"][leg] = ions
    meta["n_states"][leg] = sched.n_states
    meta["temperature"][leg] = T
    meta["md_temperature"][leg] = run.temperature
    meta["prod_ns"][leg] = prod_ns
    meta["seeds"][leg] = seeds
    meta["schedule"][leg] = {"restraint": sched.restraint, "coul": sched.coul,
                             "vdw": sched.vdw}
    meta["sources"][leg] = {"dir": str(run.path.resolve()), "replica": run.replica,
                            "top": str(run.top), "gro": str(run.md_gro),
                            "tpr": str(run.md_tpr), "xtc": str(run.md_xtc),
                            "natoms": natoms, "tc_grps": tc}
    meta["_signature"][leg] = _ligand_signature(top, lig)
    log(f"{leg}: {sched.n_states} states x {args.repeats} repeat(s), {prod_ns:g} ns prod, "
        f"T = {T:g} K" + (f" (--temperature; md.mdp {run.temperature:g} K)"
                          if args.temperature is not None else ""))


def top_name_at(top: Topology, index0: int) -> str:
    i = 0
    for name, cnt in top.molecules:
        n = len(top.require_moleculetype(name).atoms)
        if index0 < i + n * cnt:
            return name
        i += n * cnt
    raise IndexError(index0)


def node_opts(args, leg: str) -> dict:
    if args.cluster == "truba":
        return {"account": args.truba_account, "partition": args.truba_partition,
                "ntasks": args.truba_ntasks, "omp": args.truba_omp,
                "gmx_module": args.truba_gmx_module, "time": args.truba_time}
    return {"node": args.complex_node if leg == "complex" else args.solvent_node}


def run_setup(args, argv=None) -> Path:
    resolve_paths(args)
    legs = [leg for leg in LEGS if getattr(args, leg) is not None]
    if not legs:
        raise SetupError("give --solvent (hydration) or --complex and --solvent (ABFE)")
    if args.complex is not None and args.solvent is None:
        raise SetupError("--complex needs --solvent (the same ligand in water) for the "
                         "solvent leg")
    if args.repeats < 1:
        raise SetupError("--repeats must be >= 1")
    if args.cluster == "truba" and not args.truba_account:
        raise SetupError("--cluster truba needs --truba-account")
    outdir = args.outdir.resolve()
    for leg in legs:
        src = getattr(args, leg)
        if not src.is_dir():
            raise SetupError(f"--{leg} {src}: not a directory")
        if _inside(outdir, src):
            raise SetupError(f"--outdir {outdir} is inside the {leg} run folder {src}; "
                             "ABFE output is never written into a run folder")
    if outdir.exists():
        if not args.overwrite:
            raise SetupError(f"{outdir} exists; use --overwrite to replace it")
        if any(outdir.iterdir()) and not (outdir / "setup.json").is_file():
            raise SetupError(f"{outdir} is not an abfe_setup folder (no setup.json); "
                             "refusing to --overwrite it")
    sched_path = resolve_schedule(args.schedule)
    scheds = {leg: load_schedule(sched_path, leg) for leg in legs}
    for s in scheds.values():
        validate(s)
    runs = {leg: load_rundir(getattr(args, leg), args.replica) for leg in legs}
    temps = {leg: runs[leg].temperature for leg in legs}
    if args.temperature is not None:
        if args.temperature <= 0:
            raise SetupError("--temperature must be > 0 K")
    elif max(temps.values()) - min(temps.values()) > 0.5:
        raise SetupError("md.mdp temperatures differ between legs ("
                         + ", ".join(f"{l} {t:g} K" for l, t in temps.items())
                         + "): give --temperature T to run both legs at T")
    _gmx()

    outdir.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.partial_", dir=outdir.parent))
    warnings: list[str] = []
    meta = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": "abfe_setup.sh " + " ".join(
            map(str, sys.argv[1:] if argv is None else argv)),
        "options": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "outdir": str(outdir), "ligand": args.ligand, "legs": legs,
        "schedule_file": str(sched_path),
        "ligand_charge": {}, "alch_ions": {}, "n_states": {}, "temperature": {},
        "md_temperature": {}, "temperature_override": args.temperature,
        "prod_ns": {}, "seeds": {}, "schedule": {}, "sources": {}, "_signature": {},
        "seed_rule": "seed = 1000*repeat + state + 1, same rule in both legs "
                     "(gen-seed = ld-seed = seed)",
    }
    try:
        for leg in legs:
            build_leg(leg, runs[leg], partial / leg, partial, args, scheds[leg], meta,
                      warnings)
        qs = set(meta["ligand_charge"].values())
        if len(qs) != 1:
            raise SetupError(f"ligand charge differs between legs: {meta['ligand_charge']}")
        meta["ligand_charge_per_leg"] = meta["ligand_charge"]
        meta["ligand_charge"] = qs.pop()
        sig = meta.pop("_signature")
        if len(legs) == 2:
            compare_ligands(args.ligand, sig["complex"], sig["solvent"])
        sbatches = []
        for leg in legs:
            for r in range(args.repeats):
                p = write_leg_sbatch(partial / leg, leg, meta["n_states"][leg], r,
                                     args.cluster, node_opts(args, leg))
                sbatches.append(str(p.relative_to(partial)))
        write_submit(partial, sbatches)
        meta["sbatch"] = sbatches
        meta["cluster"] = {"name": args.cluster,
                           **{leg: node_opts(args, leg) for leg in legs}}
        meta["versions"] = _versions()
        meta["git_describe"] = _git_describe()
        meta["warnings"] = warnings
        (partial / "setup.json").write_text(json.dumps(meta, indent=2) + "\n")
        if outdir.exists():
            shutil.rmtree(outdir)
        partial.rename(outdir)
    except BaseException:
        shutil.rmtree(partial, ignore_errors=True)
        raise
    for w in warnings:
        log(f"WARNING: {w}")
    log(f"done: {outdir} ({', '.join(f'{l} {meta['n_states'][l]} states' for l in legs)}, "
        f"{args.repeats} repeat(s)); submit with: bash {outdir / 'submit_abfe.sh'}")
    return outdir


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        outdir = run_setup(args, argv)
    except (SetupError, ValueError, RuntimeError, FileNotFoundError) as e:
        print(f"abfe_setup: error: {e}", file=sys.stderr)
        raise SystemExit(1)
    if args.submit:
        if args.cluster == "truba":
            log("TRUBA: copy the ABFE folder to TRUBA and run 'bash submit_abfe.sh' there")
        else:
            r = subprocess.run(["bash", str(outdir / "submit_abfe.sh")], text=True)
            if r.returncode != 0:
                print("abfe_setup: error: submit_abfe.sh failed", file=sys.stderr)
                return r.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
