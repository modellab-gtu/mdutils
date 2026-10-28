#!/usr/bin/env python3
"""Build the protonated ligand (explicit H, correct bond orders and formal charges) used by the mdutils prepare scripts.

Sources, best first:
  --sdf / --mol2 FILE   user-supplied ligand: used as is (must already carry explicit hydrogens).
  --pdb FILE --ccd ID   crystal ligand from RCSB: bond orders and formal charges from the Chemical Component
                        Dictionary ideal SDF (RDKit AssignBondOrdersFromTemplate), then hydrogens added on the
                        crystal heavy-atom coordinates. The CCD form is usually neutral.
  --pdb FILE            plain PDB: Open Babel guesses bond orders and adds H from the geometry (least reliable).
Optional --ph X re-protonates with Open Babel (obabel -p X).

Writes --out (MOL2 with unique atom names, residue name --resname) and the same molecule as <out stem>.sdf.
Exits non-zero if the formal charges do not sum to --charge.
"""
import argparse, os, subprocess, sys, tempfile, urllib.request

def die(msg):
    sys.exit(f"ligand_input: ERROR: {msg}")

def note(msg):
    print(f"ligand_input: {msg}", flush=True)

def obabel(args):
    r = subprocess.run(["obabel"] + args, capture_output=True, text=True)
    if r.returncode:
        die("obabel " + " ".join(args) + "\n" + r.stderr)

def from_file(path, sanitize=True):
    from rdkit import Chem
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"):
        m = Chem.MolFromMolFile(path, removeHs=False, sanitize=sanitize)
    elif ext == ".mol2":  # RDKit's MOL2 reader is strict; go through Open Babel's SDF writer
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False).name
        obabel([path, "-O", tmp]); m = Chem.MolFromMolFile(tmp, removeHs=False, sanitize=sanitize); os.remove(tmp)
    else:
        die(f"unsupported ligand format: {path}")
    if m is None:
        die(f"RDKit could not read {path}")
    return m

def protonate(mol, ph):
    from rdkit import Chem
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, "in.sdf"), os.path.join(d, "out.sdf")
        Chem.MolToMolFile(mol, a); obabel([a, "-O", b, "-p", str(ph)])
        m = Chem.MolFromMolFile(b, removeHs=False)
    if m is None:
        die("re-protonation with Open Babel failed")
    note(f"re-protonated with Open Babel at pH {ph}")
    return m

def from_ccd(pdb, ccd):
    """Crystal heavy atoms + CCD template bond orders/charges; H added with coordinates."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    url = f"https://files.rcsb.org/ligands/download/{ccd.upper()}_ideal.sdf"
    tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False).name
    try:
        urllib.request.urlretrieve(url, tmp)
        tmpl = Chem.MolFromMolFile(tmp, removeHs=True)
    except Exception as e:
        note(f"could not fetch the CCD template {url} ({e})"); return None
    finally:
        if os.path.exists(tmp): os.remove(tmp)
    xtal = Chem.MolFromPDBFile(pdb, removeHs=True, sanitize=False)
    if tmpl is None or xtal is None:
        return None
    if xtal.GetNumAtoms() != tmpl.GetNumAtoms():
        note(f"crystal ligand has {xtal.GetNumAtoms()} heavy atoms, CCD {ccd} has {tmpl.GetNumAtoms()} (missing atoms?)")
        return None
    try:
        m = AllChem.AssignBondOrdersFromTemplate(tmpl, xtal)
        Chem.SanitizeMol(m)
    except Exception as e:
        note(f"CCD template did not match the crystal ligand ({e})"); return None
    note(f"bond orders and formal charges from CCD {ccd.upper()} (formal charge {Chem.GetFormalCharge(m):+d})")
    return m

def write_outputs(m, out, resname):
    from rdkit import Chem
    stem = os.path.splitext(out)[0]
    sdf = stem + ".sdf"
    m.SetProp("_Name", resname)
    Chem.MolToMolFile(m, sdf)
    obabel([sdf, "-O", out])
    # unique atom names (element + running number per element) and the residue name, as tleap/antechamber expect
    lines, out_lines, sect, count = open(out).read().splitlines(), [], None, {}
    for ln in lines:
        if ln.startswith("@<TRIPOS>"):
            sect = ln.strip()
            if sect != "@<TRIPOS>UNITY_ATOM_ATTR": out_lines.append(ln)
            continue
        if sect == "@<TRIPOS>ATOM" and ln.strip():
            f = ln.split()
            el = f[5].split(".")[0]
            el = el[0].upper() + el[1:].lower()
            count[el] = count.get(el, 0) + 1
            name = f"{el}{count[el]}"
            out_lines.append(f"{int(f[0]):7d} {name:<8s} {float(f[2]):10.4f} {float(f[3]):10.4f} {float(f[4]):10.4f} "
                             f"{f[5]:<9s} {1:>4d} {resname:<8s} {float(f[8]):10.4f}")
            continue
        if sect == "@<TRIPOS>UNITY_ATOM_ATTR":   # Open Babel's charge records: antechamber only warns about them
            continue
        if sect == "@<TRIPOS>MOLECULE" and ln.strip() and out_lines[-1].startswith("@<TRIPOS>MOLECULE"):
            out_lines.append(resname); continue
        out_lines.append(ln)
    open(out, "w").write("\n".join(out_lines) + "\n")
    return sdf

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sdf"); src.add_argument("--mol2"); src.add_argument("--pdb")
    ap.add_argument("--ccd", help="CCD id for --pdb (RCSB ligands)")
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--ph", type=float, default=None)
    ap.add_argument("--resname", required=True)
    ap.add_argument("--out", required=True, help="output MOL2 with hydrogens")
    a = ap.parse_args()
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.warning")

    if a.sdf or a.mol2:
        m = from_file(a.sdf or a.mol2)
        route = "user file"
    else:
        m = from_ccd(a.pdb, a.ccd) if a.ccd else None
        route = "CCD template"
        if m is None:
            note("WARNING: bond orders are guessed by Open Babel from the PDB geometry; check the result, or give "
                 "--local-ligand lig.sdf with explicit H and formal charges")
            with tempfile.TemporaryDirectory() as d:
                t = os.path.join(d, "h.sdf"); obabel([a.pdb, "-O", t, "-h"]); m = from_file(t)
            route = "Open Babel guess"
    if a.ph is not None:
        m = protonate(m, a.ph)
    elif route == "CCD template":
        m = Chem.AddHs(m, addCoords=True)
    missing = sum(at.GetNumImplicitHs() for at in m.GetAtoms())
    if missing:
        die(f"{missing} hydrogen(s) are implicit in the input: supply the ligand with explicit hydrogens "
            "(e.g. obabel lig.sdf -O lig_h.sdf -h, after checking the protonation) or use --lig-ph")
    q = Chem.GetFormalCharge(m)
    nH = sum(1 for at in m.GetAtoms() if at.GetAtomicNum() == 1)
    sdf = write_outputs(m, a.out, a.resname)
    note(f"{a.resname}: {m.GetNumAtoms()} atoms ({nH} H), formal charge {q:+d}, source: {route} -> {a.out}, {sdf}")
    if q != a.charge:
        die(f"the ligand's formal charges sum to {q:+d} but --charge is {a.charge:+d}. "
            + ("The CCD form is usually neutral: for an ionised ligand give a protonated --local-ligand SDF "
               "(or --lig-ph 7.4) together with the matching --charge." if route == "CCD template" else
               "Fix the protonation of the ligand file or --charge."))

if __name__ == "__main__":
    main()
