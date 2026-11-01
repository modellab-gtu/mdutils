"""ligparam (repo copy) on the aspirin example: outputs written, charges sum to the formal charge."""
import json
import subprocess

OUTPUTS = ["LIG.mol2", "LIG.frcmod", "LIG.lib", "LIG.prmtop", "LIG.inpcrd", "LIG.itp", "LIG.top", "LIG.gro",
           "LIG_system.xml", "ligparam.json", "ligparam.log"]


def mol2_charges(path):
    txt = path.read_text()
    atoms = txt.split("@<TRIPOS>ATOM\n", 1)[1].split("@<TRIPOS>", 1)[0]
    return [float(r.split()[8]) for r in atoms.splitlines() if r.strip()]


def test_aspirin_default(tmp_path, example, ligparam_bin):
    out = tmp_path / "asp"
    r = subprocess.run([str(ligparam_bin), "-i", str(example("aspirin.sdf")), "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    for f in OUTPUTS:
        assert (out / f).is_file(), f
    rec = json.loads((out / "ligparam.json").read_text())
    assert rec["net_charge"] == 0 and rec["charge_model"] == "abcg2"
    assert abs(rec["charge_sum"]) < 1e-4
    assert abs(sum(mol2_charges(out / "LIG.mol2"))) < 1e-5
