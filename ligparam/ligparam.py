#!/usr/bin/env python3
"""ligparam: small-molecule force-field parameters for Amber, GROMACS and OpenMM in one call (GTU ModelLab).

GAFF2 route : SDF/MOL2 -> antechamber (gaff2 types; bcc | abcg2 | nagl | user charges) -> parmchk2 -s 2 -> tleap
              -> LIG.prmtop/inpcrd/lib/frcmod/mol2, ParmEd -> LIG.top/LIG.itp/LIG.gro, LIG.xml (OpenMM).
OpenFF route: SDF -> OpenFF Sage (openff-2.2.1) + NAGL or AM1-BCC charges -> Interchange -> same output set.
Net charge  : taken from the SDF formal charges (--charge auto), never re-perceived from a PDB.
Geometry    : --opt am1 (default) = antechamber's own AM1 minimisation, as always (abcg2/bcc). Any other --opt
              relaxes the ligand first (hf = Gaussian HF/6-31G*; aimnet2|uma|ani2x|mace-off = NNP via covparam.opt;
              none = as given) and runs AM1 with maxcyc=0 on that geometry. The geometry is used for the charges only:
              every output keeps the input coordinates.
"""
import argparse, json, os, shlex, shutil, subprocess as sp, sys, datetime

def die(msg): sys.exit(f"ligparam: ERROR: {msg}")
def run(cmd, log):
    with open(log, "a") as fh:
        fh.write(f"\n$ {cmd}\n"); fh.flush()
        r = sp.run(cmd, shell=True, stdout=fh, stderr=sp.STDOUT)
    if r.returncode: die(f"command failed (see {log}): {cmd}")

def read_mol(path):
    from rdkit import Chem
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"): m = Chem.MolFromMolFile(path, removeHs=False)
    elif ext == ".mol2": m = Chem.MolFromMol2File(path, removeHs=False)
    else: die("input must be .sdf/.mol or .mol2 (bond orders and formal charges are needed; PDB is not accepted)")
    if m is None: die(f"RDKit could not read {path}")
    if m.GetNumConformers() == 0 or not m.GetConformer().Is3D(): die("input has no 3D coordinates")
    return m

def check_hydrogens(m):
    """Warn if heavy atoms look under-protonated (implicit H left by the SDF)."""
    missing = sum(a.GetNumImplicitHs() for a in m.GetAtoms())
    if missing: die(f"{missing} implicit hydrogen(s): add explicit H first (e.g. obabel -h, or protonate with --protonate-ph)")

def main():
    ap = argparse.ArgumentParser(prog="ligparam", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="ligand .sdf/.mol2 with explicit H and 3D coordinates")
    ap.add_argument("-o", "--outdir", default=None, help="output folder (default: <input stem>_<ff>_<charges>)")
    ap.add_argument("--ff", choices=["gaff2", "openff"], default="gaff2")
    ap.add_argument("--charge-model", choices=["abcg2", "bcc", "nagl", "user"], default=None,
                    help="gaff2: abcg2 (default) | bcc | nagl | user (charges from the MOL2); openff: nagl (default) | bcc")
    ap.add_argument("--charge", default="auto", help="net charge: auto = sum of SDF formal charges (default), or an integer to cross-check")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--openff-ff", default="openff-2.2.1.offxml")
    ap.add_argument("--protonate-ph", type=float, default=None, help="re-protonate with Open Babel at this pH first (changes the net charge!)")
    ap.add_argument("--opt", choices=OPT_CHOICES, default=None,
                    help="geometry for the AM1 step of abcg2/bcc charges: am1 (default; antechamber's AM1 minimisation) | "
                         "hf (Gaussian HF/6-31G* Opt on cws02; writes opt_hf.gjf + run_opt.sh, then rerun with --geom) | "
                         "aimnet2 | uma | ani2x | mace-off (NNP; needs covparam) | none (input geometry). Not am1: AM1 runs "
                         "with maxcyc=0 on that geometry. Charges only: all outputs keep the input coordinates")
    ap.add_argument("--geom", default=None, help="--opt hf: the optimised geometry (Gaussian .log from run_opt.sh, or .xyz "
                                                 "in the input's atom order); implies --opt hf")
    ap.add_argument("--local", action="store_true", help="--opt hf: run Gaussian here in the foreground (no Slurm job)")
    a = ap.parse_args()

    cm = a.charge_model or ("abcg2" if a.ff == "gaff2" else "nagl")
    if a.ff == "openff" and cm not in ("nagl", "bcc"): die("--ff openff takes --charge-model nagl or bcc")
    if cm == "user" and not a.input.lower().endswith(".mol2"): die("--charge-model user needs a .mol2 with charges")
    opt = a.opt or ("hf" if a.geom else "am1")
    if a.geom and opt != "hf": die(f"--geom is the result of the --opt hf job; it can't be combined with --opt {opt}")
    if opt != "am1" and a.ff != "gaff2": die("--opt applies to --ff gaff2 (OpenFF assigns its own charges); drop --opt")
    if opt != "am1" and cm not in ("abcg2", "bcc"):
        die(f"--opt sets the geometry of the AM1 step of abcg2/bcc charges; --charge-model {cm} has no AM1 step "
            "(and the outputs keep the input coordinates), so --opt would change nothing: drop it")
    if a.local and opt != "hf": die("--local applies to --opt hf")
    if a.geom: a.geom = os.path.abspath(a.geom)
    a.input = os.path.abspath(a.input); inp = a.input; stem = os.path.splitext(os.path.basename(inp))[0]
    out = os.path.abspath(a.outdir or f"{stem}_{a.ff}_{cm}")
    os.makedirs(out, exist_ok=True); log = os.path.join(out, "ligparam.log"); open(log, "w").close()
    R = a.resname[:3].upper()

    from rdkit import Chem
    if a.protonate_ph is not None:
        src = os.path.join(out, "input_pH.sdf")
        run(f"obabel '{inp}' -O '{src}' -p {a.protonate_ph}", log); inp = src
    m = read_mol(inp); check_hydrogens(m)
    q_formal = Chem.GetFormalCharge(m)
    if a.charge != "auto":
        if int(a.charge) != q_formal:
            die(f"--charge {a.charge} disagrees with the SDF formal charges (sum {q_formal}); fix the input or its protonation")
    q = q_formal
    nH = sum(1 for x in m.GetAtoms() if x.GetAtomicNum() == 1)
    opt_info = {}
    if opt != "am1":
        mq, opt_info = pre_optimise(m, q, opt, a, out, log)
        opt = opt_info["opt"]
        mq.SetProp("_Name", R); w = Chem.SDWriter(os.path.join(out, "charge_geom.sdf")); w.write(mq); w.close()
    sdf = os.path.join(out, "input.sdf")
    m.SetProp("_Name", R); w = Chem.SDWriter(sdf); w.write(m); w.close()
    os.chdir(out)

    if a.ff == "gaff2":
        if cm in ("abcg2", "bcc"):
            if opt == "am1":
                run(f"antechamber -i input.sdf -fi sdf -o {R}.mol2 -fo mol2 -c {cm} -nc {q} -at gaff2 -rn {R} -pf y -dr no", log)
            else:  # charges on the --opt geometry (AM1 single point), then the input coordinates back into the MOL2
                run(f"antechamber -i charge_geom.sdf -fi sdf -o {R}.mol2 -fo mol2 -c {cm} -nc {q} -at gaff2 -rn {R} -pf y "
                    f"-dr no -ek \"{AM1_SP}\"", log)
                set_mol2_coords(f"{R}.mol2", m.GetConformer().GetPositions())
        elif cm == "user":
            run(f"antechamber -i '{os.path.abspath(a.input)}' -fi mol2 -o {R}.mol2 -fo mol2 -at gaff2 -rn {R} -pf y -dr no", log)
        else:  # nagl charges on gaff2 types
            run(f"antechamber -i input.sdf -fi sdf -o {R}_types.mol2 -fo mol2 -at gaff2 -rn {R} -pf y -dr no", log)
            from openff.toolkit import Molecule
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
            om = Molecule.from_file("input.sdf"); om.assign_partial_charges(NAGL_MODEL, toolkit_registry=NAGLToolkitWrapper())
            with open(f"{R}.chg", "w") as fh: fh.write("\n".join(f"{c:.6f}" for c in om.partial_charges.m) + "\n")
            run(f"antechamber -i {R}_types.mol2 -fi mol2 -o {R}.mol2 -fo mol2 -c rc -cf {R}.chg -at gaff2 -rn {R} -pf y -dr no", log)
        fix_mol2_charge(f"{R}.mol2", q)
        # -s 2 = GAFF2. Default is GAFF1, and the word form "-s gaff2" is unreliable in AT25 parmchk2.
        run(f"parmchk2 -i {R}.mol2 -f mol2 -o {R}.frcmod -s 2", log)
        with open("tleap.in", "w") as fh:
            fh.write(f"source leaprc.gaff2\n{R} = loadmol2 {R}.mol2\nloadamberparams {R}.frcmod\ncheck {R}\n"
                     f"saveoff {R} {R}.lib\nsaveamberparm {R} {R}.prmtop {R}.inpcrd\nquit\n")
        run("tleap -f tleap.in", log)
        if not os.path.exists(f"{R}.prmtop"): die("tleap did not write the prmtop (see leap.log)")
        import parmed as pmd
        s = pmd.load_file(f"{R}.prmtop", f"{R}.inpcrd")
        s.save(f"{R}.top", overwrite=True); s.save(f"{R}.gro", overwrite=True)
        write_itp(s, R)
        charges = [x.charge for x in s.atoms]; ff_desc = "GAFF2 (gaff2.dat + parmchk2 -s 2)"
        try:
            import openmm.app as app, openmm.unit as u
            sysm = s.createSystem(nonbondedMethod=app.NoCutoff)
            import openmm; open(f"{R}_system.xml", "w").write(openmm.XmlSerializer.serialize(sysm))
        except Exception as e: print(f"ligparam: note: OpenMM system XML not written ({e})")
    else:
        from openff.toolkit import Molecule, ForceField
        from openff.interchange import Interchange
        from openff.units import unit
        om = Molecule.from_file("input.sdf")
        for at in om.atoms: at.metadata["residue_name"] = R
        if cm == "nagl":
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
            om.assign_partial_charges(NAGL_MODEL, toolkit_registry=NAGLToolkitWrapper())
        else:
            from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
            om.assign_partial_charges("am1bcc", use_conformers=om.conformers, toolkit_registry=AmberToolsToolkitWrapper())
        ic = Interchange.from_smirnoff(ForceField(a.openff_ff), [om], charge_from_molecules=[om])
        ic.box = [5, 5, 5] * unit.nanometer
        ic.to_gromacs(prefix=R); ic.to_prmtop(f"{R}.prmtop"); ic.to_inpcrd(f"{R}.inpcrd")
        import parmed as pmd
        write_itp(pmd.load_file(f"{R}.top", xyz=f"{R}.gro"), R)
        import openmm; open(f"{R}_system.xml", "w").write(openmm.XmlSerializer.serialize(ic.to_openmm()))
        charges = list(om.partial_charges.m); ff_desc = a.openff_ff

    qsum = sum(charges)
    if abs(qsum - q) > 0.01: die(f"partial charges sum to {qsum:.4f}, formal charge is {q}")
    rec = {"input": a.input, "resname": R, "net_charge": q, "n_atoms": m.GetNumAtoms(), "n_H": nH,
           "forcefield": ff_desc, "charge_model": cm, "charge_sum": round(qsum, 4),
           "opt": opt if (a.ff == "gaff2" and cm in ("abcg2", "bcc")) else None,
           "versions": versions(), "date": datetime.datetime.now().isoformat(timespec="seconds")}
    rec.update({k: v for k, v in opt_info.items() if k != "opt" and v is not None})
    json.dump(rec, open("ligparam.json", "w"), indent=1)
    print(f"ligparam: {R}  net charge {q:+d}  charges sum {qsum:+.4f}  ({ff_desc}, {cm})  ->  {out}")
    if opt != "am1":
        print(f"  charge geometry: --opt {opt} ({opt_info.get('opt_note')}); AM1 maxcyc=0 on it. Outputs keep the input coordinates")
    print(f"  Amber: {R}.prmtop {R}.inpcrd" + (f" {R}.mol2 {R}.frcmod {R}.lib" if a.ff == "gaff2" else ""))
    print(f"  GROMACS: {R}.itp (moleculetype + atomtypes, #include in your topology) {R}.top {R}.gro   OpenMM: {R}_system.xml")

NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
OPT_CHOICES = ["am1", "hf", "aimnet2", "uma", "ani2x", "mace-off", "none"]
# antechamber's own sqm keywords plus maxcyc=0. -ek REPLACES antechamber's keyword list: "maxcyc=0" alone makes
# sqm fall back to its PM3 default, and BCC/ABCG2 corrections on PM3 charges are wrong.
AM1_SP = "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, maxcyc=0"
G16_MODULE, G16_NODE, G16_CPUS, G16_MEM = "gaussian/g16-avx2", "cws02", 16, "32G"

def need_covparam(opt):
    """covparam (NNP relaxation) is a separate package; a clean message instead of a traceback when it is missing or unreadable."""
    try:
        import covparam.opt, covparam.errors, covparam.model  # noqa: F401
    except Exception as e:
        die(f"--opt {opt} needs covparam (module load covparam/1.0, or use --opt hf / none): cannot import it ({type(e).__name__}: {e})")

def pre_optimise(m, q, opt, a, out, log):
    """Geometry for the charges (--opt other than am1): (mol with that geometry, provenance dict for ligparam.json).
    hf without --geom/--local writes the Gaussian job (+ opt_hf.json with the requested model) and exits."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    info = {"opt": opt, "opt_requested": opt, "opt_note": None, "geom": None}
    if a.geom:
        side = os.path.join(os.path.dirname(a.geom), "opt_hf.json")
        if os.path.isfile(side):  # written with the job: keeps a requested NNP model and the fallback reason
            try: info.update({k: v for k, v in json.load(open(side)).items() if k in ("opt_requested", "opt_note")})
            except (OSError, ValueError): pass
        x = read_geom(a.geom, m); info["geom"] = a.geom
        info["opt_note"] = "; ".join(t for t in (info["opt_note"], f"HF/6-31G* geometry from {a.geom}") if t)
    elif opt == "none":
        info["opt_note"] = "input geometry, no optimisation"
        return Chem.Mol(m), info
    else:
        if opt != "hf":
            need_covparam(opt)
            from covparam.opt import optimise
            from covparam.errors import CovparamCheckError, CovparamInputError
            try:
                res = optimise(nnp_fragment(m, q), opt, os.path.join(out, f"opt_{opt}"), log)
            except (CovparamCheckError, CovparamInputError) as e: die(str(e))
            if res.model_used == "hf":
                print(f"ligparam: note: {res.note}"); info.update(opt="hf", opt_note=res.note)
            else:
                x = res.mol.GetConformer().GetPositions(); info["opt_note"] = res.note
        if info["opt"] == "hf":
            write_hf_job(m, q, out)
            json.dump({k: info[k] for k in ("opt_requested", "opt_note")}, open(os.path.join(out, "opt_hf.json"), "w"), indent=1)
            if not a.local:
                print(f"ligparam: wrote {out}/opt_hf.gjf and {out}/run_opt.sh (HF/6-31G* optimisation, {G16_NODE}, "
                      f"{G16_CPUS} cores). Submit it, then rerun with the optimised geometry:\n"
                      f"  sbatch {out}/run_opt.sh\n  ligparam {rerun_args(out)}")
                sys.exit(0)
            run("bash -lc " + shlex.quote(f"module load {G16_MODULE} && export GAUSS_SCRDIR=$(mktemp -d /tmp/ligparam_g16.XXXXXX)"
                                         f" && trap 'rm -rf \"$GAUSS_SCRDIR\"' EXIT && cd {shlex.quote(out)} && "
                                         f'"$g16root/g16/g16" opt_hf.gjf'), log)
            info["geom"] = os.path.join(out, "opt_hf.log"); x = read_geom(info["geom"], m)
            info["opt_note"] = "; ".join(t for t in (info["opt_note"], "Gaussian HF/6-31G* (--local)") if t)
    mq = Chem.Mol(m); conf = mq.GetConformer()
    for i, p in enumerate(x): conf.SetAtomPosition(i, Point3D(*map(float, p)))
    return mq, info

def nnp_fragment(m, q):
    from covparam.model import Fragment
    return Fragment(mol=m, cap_atom_idx=set(), backbone_idx={}, residue_idx=list(range(m.GetNumAtoms())), net_charge=q)

def hf_opt_route(natoms):
    """covparam's Gaussian Opt keyword for a free ligand (no frozen coordinates; kept here so hf needs no covparam):
    MaxCycles=500, CalcFC (analytic HF Hessian at the start) up to 50 atoms, MaxStep=10 (smaller trust radius)."""
    return "Opt=(" + ",".join(["MaxCycles=500"] + (["CalcFC"] if natoms <= 50 else []) + ["MaxStep=10"]) + ")"

def write_hf_job(m, q, out):
    """opt_hf.gjf (HF/6-31G* Opt, no frozen coordinates) + run_opt.sh (cws02, node-local scratch)."""
    conf = m.GetConformer()
    xyz = [f"{at.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}"
           for at, p in ((at, conf.GetAtomPosition(at.GetIdx())) for at in m.GetAtoms())]
    open(os.path.join(out, "opt_hf.gjf"), "w").write("\n".join(
        ["%chk=opt_hf.chk", f"%nprocshared={G16_CPUS}", "%mem=30GB", f"#P HF/6-31G* {hf_opt_route(m.GetNumAtoms())} NoSymm",
         "", "ligparam HF/6-31G* optimisation", "", f"{q} 1"] + xyz) + "\n\n")
    sh = os.path.join(out, "run_opt.sh")
    open(sh, "w").write(f"""#!/bin/bash -l
#SBATCH -J lp-opt-hf
#SBATCH -w {G16_NODE}
#SBATCH --export=NONE
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {G16_CPUS}
#SBATCH --mem={G16_MEM}
#SBATCH -o {out}/run_opt.slurm-%j.out
# ligparam --opt hf: HF/6-31G* optimisation (written by ligparam). Afterwards: ligparam ... --opt hf --geom opt_hf.log
module purge
module load {G16_MODULE}
export GAUSS_SCRDIR=/tmp/$SLURM_JOB_ID
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT
# the module APPENDS to PATH, so a personal g16 earlier on PATH would win: call the module's binary
cd {shlex.quote(out)}
"$g16root/g16/g16" opt_hf.gjf
""")
    os.chmod(sh, 0o755)

def rerun_args(out):
    """This command line with --opt/--local replaced by --opt hf --geom OUT/opt_hf.log (opt_hf.json keeps the requested model)."""
    args, skip = [], False
    for x in sys.argv[1:]:
        if skip: skip = False; continue
        if x == "--opt": skip = True; continue
        if x.startswith("--opt=") or x == "--local": continue
        args.append(x)
    return " ".join(shlex.quote(x) for x in args + ["--opt", "hf", "--geom", f"{out}/opt_hf.log"])

def read_geom(path, m):
    """Coordinates (atom order of m) from a Gaussian optimisation log or an .xyz file."""
    import re
    from rdkit import Chem
    pt = Chem.GetPeriodicTable(); n = m.GetNumAtoms(); want = [at.GetAtomicNum() for at in m.GetAtoms()]
    if not os.path.isfile(path): die(f"--geom {path}: no such file")
    txt = open(path, errors="replace").read()
    if path.lower().endswith((".log", ".out")):
        lines = txt.rstrip().splitlines()
        if not lines or "Normal termination" not in lines[-1]:
            die(f"--geom {path}: Gaussian did not finish normally (still running, or failed: see the end of the log)")
        if "Optimization completed" not in txt:
            die(f"--geom {path}: no converged optimisation in this log")
        blocks = (re.findall(r"Input orientation:.*?-{20,}\n.*?-{20,}\n(.*?)\n *-{20,}", txt, re.S)
                  or re.findall(r"Standard orientation:.*?-{20,}\n.*?-{20,}\n(.*?)\n *-{20,}", txt, re.S))
        if not blocks: die(f"--geom {path}: no geometry found (not a Gaussian log?)")
        rows = [r.split() for r in blocks[-1].splitlines() if r.strip()]
        z = [int(r[1]) for r in rows]; x = [[float(v) for v in r[3:6]] for r in rows]
    else:
        rows = [r.split() for r in txt.splitlines()[2:] if r.strip()]
        try:
            z = [int(r[0]) if r[0].isdigit() else pt.GetAtomicNumber(r[0].capitalize()) for r in rows]
            x = [[float(v) for v in r[1:4]] for r in rows]
        except (ValueError, IndexError, RuntimeError): die(f"--geom {path}: not an .xyz file (or a Gaussian .log)")
    if len(z) != n: die(f"--geom {path}: {len(z)} atoms, the input has {n} atoms")
    if z != want: die(f"--geom {path}: atoms are not in the input's order (elements differ at atom {next(i for i in range(n) if z[i] != want[i]) + 1})")
    return x

def set_mol2_coords(path, xyz):
    """Replace the x, y, z columns of a MOL2's atom rows (--opt: charges from the optimised geometry, coordinates of the input)."""
    txt = open(path).read()
    head, rest = txt.split("@<TRIPOS>ATOM\n", 1)
    atoms, tail = rest.split("@<TRIPOS>", 1)
    import re
    rows = atoms.splitlines(); out = []
    for r, p in zip(rows, xyz):  # keep the row's layout: each new number right-aligned in the old field's width
        g = re.match(r"(\s*\S+\s+\S+)(\s+\S+)(\s+\S+)(\s+\S+)(.*)$", r).groups()
        out.append(g[0] + "".join(f"{v:.4f}".rjust(len(w)) for v, w in zip(p, g[1:4])) + g[4])
    open(path, "w").write(head + "@<TRIPOS>ATOM\n" + "\n".join(out) + "\n@<TRIPOS>" + tail)

def fix_mol2_charge(path, q):
    """antechamber rounds charges to 6 decimals, so sums like 0.999 or -2.001 occur: spread the residual over the
    heavy atoms in 1e-6 steps so the total is exactly the formal charge (as GROMACS/tleap expect)."""
    txt = open(path).read()
    head, rest = txt.split("@<TRIPOS>ATOM\n", 1)
    atoms, tail = rest.split("@<TRIPOS>", 1)
    rows = atoms.splitlines(); ch = [round(float(r.split()[8]) * 1e6) for r in rows]
    resid = round(q * 1e6) - sum(ch)
    if resid == 0: return
    heavy = [i for i, r in enumerate(rows) if not r.split()[5].lower().startswith("h")] or list(range(len(rows)))
    step = 1 if resid > 0 else -1
    for k in range(abs(resid)): ch[heavy[k % len(heavy)]] += step
    out = []
    for r, c in zip(rows, ch):
        f = r.split(); idx = r.rfind(f[8])
        out.append(f"{r[:idx]}{c / 1e6:.6f}{r[idx + len(f[8]):]}")
    open(path, "w").write(head + "@<TRIPOS>ATOM\n" + "\n".join(out) + "\n@<TRIPOS>" + tail)

def write_itp(s, R):
    """GROMACS .itp for inclusion in a protein topology: [ atomtypes ] then the [ moleculetype ] block."""
    import parmed as pmd, io, re
    top = pmd.gromacs.GromacsTopologyFile.from_structure(s)
    buf = f"{R}._tmp.top"; top.write(buf)
    txt = open(buf).read(); os.remove(buf)
    sec = re.split(r"(?m)^(?=\[ )", txt)
    keep = [b for b in sec if b.startswith("[ atomtypes ]")]
    mol = []; on = False
    for b in sec:
        if b.startswith("[ moleculetype ]"): on = True
        if b.startswith("[ system ]"): on = False
        if on: mol.append(b)
    open(f"{R}.itp", "w").write(f"; {R} from ligparam (GTU ModelLab). [ atomtypes ] must come before any [ moleculetype ]\n"
                                f"; in your topology: merge these atomtypes with the protein force field's if both are present.\n\n"
                                + "".join(keep) + "\n" + "".join(mol))

def versions():
    v = {}
    import glob
    for f in glob.glob(os.path.join(os.environ.get("AMBERHOME", sys.prefix), "conda-meta", "ambertools-dac-*.json")):
        v["ambertools (antechamber)"] = json.load(open(f))["version"]
    for mod in ("rdkit", "parmed", "openff.toolkit", "openff.interchange", "openff.nagl", "openmm"):
        try:
            x = __import__(mod, fromlist=["x"]); v[mod] = getattr(x, "__version__", getattr(getattr(x, "version", None), "version", "?"))
        except Exception: pass
    return v

if __name__ == "__main__":
    main()
