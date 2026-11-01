#!/usr/bin/env python3
"""ligparam: small-molecule force-field parameters for Amber, GROMACS and OpenMM in one call (GTU ModelLab).

GAFF2 route : SDF/MOL2 -> antechamber (gaff2 types; bcc | abcg2 | nagl | user charges) -> parmchk2 -s 2 -> tleap
              -> LIG.prmtop/inpcrd/lib/frcmod/mol2, ParmEd -> LIG.top/LIG.itp/LIG.gro, LIG.xml (OpenMM).
OpenFF route: SDF -> OpenFF Sage (openff-2.2.1) + NAGL or AM1-BCC charges -> Interchange -> same output set.
Net charge  : taken from the SDF formal charges (--charge auto), never re-perceived from a PDB.
"""
import argparse, json, os, shutil, subprocess as sp, sys, datetime

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
    a = ap.parse_args()

    cm = a.charge_model or ("abcg2" if a.ff == "gaff2" else "nagl")
    if a.ff == "openff" and cm not in ("nagl", "bcc"): die("--ff openff takes --charge-model nagl or bcc")
    if cm == "user" and not a.input.lower().endswith(".mol2"): die("--charge-model user needs a .mol2 with charges")
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
    sdf = os.path.join(out, "input.sdf")
    m.SetProp("_Name", R); w = Chem.SDWriter(sdf); w.write(m); w.close()
    os.chdir(out)

    if a.ff == "gaff2":
        if cm in ("abcg2", "bcc"):
            run(f"antechamber -i input.sdf -fi sdf -o {R}.mol2 -fo mol2 -c {cm} -nc {q} -at gaff2 -rn {R} -pf y -dr no", log)
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
           "versions": versions(), "date": datetime.datetime.now().isoformat(timespec="seconds")}
    json.dump(rec, open("ligparam.json", "w"), indent=1)
    print(f"ligparam: {R}  net charge {q:+d}  charges sum {qsum:+.4f}  ({ff_desc}, {cm})  ->  {out}")
    print(f"  Amber: {R}.prmtop {R}.inpcrd" + (f" {R}.mol2 {R}.frcmod {R}.lib" if a.ff == "gaff2" else ""))
    print(f"  GROMACS: {R}.itp (moleculetype + atomtypes, #include in your topology) {R}.top {R}.gro   OpenMM: {R}_system.xml")

NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"

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
