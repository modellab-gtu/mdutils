import subprocess

import pytest

from automd.fe.abfe.rundir import load_rundir
from automd.fe.abfe.topology import Topology, ligand_charge, preprocess


@pytest.fixture
def pp(tmp_path, gmx):
    def _pp(rd):
        return preprocess(load_rundir(rd), tmp_path / rd.name)
    return _pp


def test_preprocess_both_routes(pp, lig_rundir, prolig_rundir, parmed_prolig_rundir):
    for rd in (lig_rundir, prolig_rundir, parmed_prolig_rundir):
        p = pp(rd)
        assert p.is_file()
        assert "#include" not in p.read_text()


def test_charge_neutral(pp, lig_rundir, prolig_rundir, parmed_prolig_rundir,
                        charged_lig_rundir):
    for rd in (lig_rundir, prolig_rundir, parmed_prolig_rundir, charged_lig_rundir):
        t = Topology.parse(pp(rd))
        assert abs(t.charge("A")) < 1e-3
        assert abs(t.charge("B")) < 1e-3


def test_ligand_charge(pp, lig_rundir, prolig_rundir, charged_lig_rundir):
    assert ligand_charge(Topology.parse(pp(lig_rundir)), "ACE") == -1
    assert ligand_charge(Topology.parse(pp(prolig_rundir)), "JZ4") == 0
    assert ligand_charge(Topology.parse(pp(charged_lig_rundir)), "LID") == 1


def test_parmed_route_names_and_repeated_molecules(pp, parmed_prolig_rundir):
    t = Topology.parse(pp(parmed_prolig_rundir))
    assert "system1" in t.moleculetypes and "JZ4" in t.moleculetypes
    names = [n for n, _ in t.molecules]
    assert names.count("NA") > 1 and isinstance(t.molecules, list)


def test_missing_ligand_moleculetype(pp, lig_rundir):
    t = Topology.parse(pp(lig_rundir))
    with pytest.raises(ValueError) as e:
        ligand_charge(t, "XYZ")
    assert "XYZ" in str(e.value) and "ACE" in str(e.value)


def test_b_state_defaults_to_a(pp, lig_rundir):
    a = Topology.parse(pp(lig_rundir)).moleculetypes["ACE"].atoms[0]
    assert (a.typeB, a.qB, a.massB) == (a.type, a.q, a.mass)


def test_write_verbatim_when_unedited(pp, lig_rundir, tmp_path):
    src = pp(lig_rundir)
    Topology.parse(src).write(tmp_path / "o.top")
    assert (tmp_path / "o.top").read_text() == src.read_text()


def test_edit_roundtrip(pp, lig_rundir, tmp_path):
    t = Topology.parse(pp(lig_rundir))
    t.moleculetypes["ACE"].atoms[0].qB = t.moleculetypes["ACE"].atoms[0].q + 0.5
    t.write(tmp_path / "e.top")
    t2 = Topology.parse(tmp_path / "e.top")
    assert abs(t2.charge("B") - t2.charge("A") - 0.5) < 1e-4


def _mol_block(top, name):
    """Lines of the processed top belonging to moleculetype `name`."""
    out, inside = [], False
    for ln in top.read_text().splitlines():
        s = ln.split(";")[0].strip()
        if s.startswith("[") and s.strip("[] ").lower() == "moleculetype":
            inside = False
        out.append(ln) if inside else None
        if s.startswith("[") and s.strip("[] ").lower() == "moleculetype":
            inside = "pending"
        elif inside == "pending" and s:
            inside = s.split()[0] == name
            if inside:
                out.append(ln)
    return "\n".join(out)


def test_processed_has_no_define_restraints(pp, lig_rundir, prolig_rundir):
    lig = pp(lig_rundir)
    assert "#ifdef" not in lig.read_text()
    assert "[ position_restraints ]" not in _mol_block(lig, "ACE")
    assert "[ settles ]" in _mol_block(lig, "SOL")
    pl = pp(prolig_rundir)
    assert "#ifdef" not in pl.read_text()
    names = list(Topology.parse(pl).moleculetypes)
    assert "[ position_restraints ]" not in _mol_block(pl, names[0])
    assert "[ position_restraints ]" not in pl.read_text()


def _grompp(gmx, tmp_path, rd, top):
    run = load_rundir(rd)
    mdp = tmp_path / "x.mdp"
    mdp.write_text("integrator = steep\nnsteps = 0\n")
    return subprocess.run(
        [gmx, "grompp", "-f", str(mdp), "-c", str(run.md_gro), "-p", str(top),
         "-o", str(tmp_path / "x.tpr"), "-maxwarn", "2"],
        cwd=tmp_path, capture_output=True, text=True)


@pytest.mark.parametrize("which", ["lig_rundir", "parmed_prolig_rundir"])
def test_roundtrip_grompp(gmx, pp, tmp_path, request, which):
    rd = request.getfixturevalue(which)
    t = Topology.parse(pp(rd))
    t.write(tmp_path / "rt.top")
    r = _grompp(gmx, tmp_path, rd, tmp_path / "rt.top")
    assert r.returncode == 0, r.stderr[-1500:]


def test_structural_edits(gmx, pp, parmed_prolig_rundir, tmp_path):
    rd = parmed_prolig_rundir
    t = Topology.parse(pp(rd))
    na = t.moleculetypes["NA"]
    q_before = t.charge("A")
    block = ("[ moleculetype ]\nNA_ALCH 3\n\n[ atoms ]\n"
             "1 " + na.atoms[0].type + " 1 NA NA 1 1.0 22.99 "
             + na.atoms[0].type + " 0.0 22.99\n\n[ position_restraints ]\n"
             "1 1 1000 1000 1000\n")
    t.add_moleculetype("NA_ALCH", block)
    assert "NA_ALCH" in t.moleculetypes
    mols = list(t.molecules)
    i = [k for k, (n, _) in enumerate(mols) if n == "NA"][0]
    assert mols[i][1] == 1
    mols.pop(i)
    mols.append(("NA_ALCH", 1))
    t.set_molecules(mols)
    t.append_section("[ intermolecular_interactions ]\n[ distance_restraints ]\n"
                     "; dummy\n")
    t.write(tmp_path / "s.top")
    t2 = Topology.parse(tmp_path / "s.top")
    assert "NA_ALCH" in t2.moleculetypes
    assert t2.molecules[-1] == ("NA_ALCH", 1)
    assert abs(t2.charge("A") - q_before) < 1e-3
    assert abs(t2.charge("B") - (t2.charge("A") - 1)) < 1e-3
    text = (tmp_path / "s.top").read_text()
    assert text.rstrip().endswith("; dummy")
    assert "[ position_restraints ]" in text
    # atom edits made before and after structural ops both survive
    t2.moleculetypes["NA_ALCH"].atoms[0].qB = 0.5
    t2.append_section("; end")
    t2.moleculetypes["NA_ALCH"].atoms[0].qB = 0.25
    t2.write(tmp_path / "s2.top")
    assert abs(Topology.parse(tmp_path / "s2.top").moleculetypes["NA_ALCH"].atoms[0].qB - 0.25) < 1e-6


def test_add_after_and_duplicates(pp, lig_rundir, tmp_path):
    t = Topology.parse(pp(lig_rundir))
    with pytest.raises(ValueError):
        t.add_moleculetype("ACE", "[ moleculetype ]\nACE 3\n[ atoms ]\n")
    t.add_moleculetype("X1", "[ moleculetype ]\nX1 3\n[ atoms ]\n"
                       "1 c3 1 X1 C1 1 0.0 12.01\n", after="ACE")
    names = list(t.moleculetypes)
    assert names.index("X1") == names.index("ACE") + 1
    t.moleculetypes["Y"] = type(t.moleculetypes["X1"])("Y")
    with pytest.raises(ValueError):
        t.write(tmp_path / "bad.top")


def test_duplicate_moleculetype_in_file(tmp_path):
    f = tmp_path / "d.top"
    f.write_text("[ moleculetype ]\nA 3\n[ atoms ]\n[ moleculetype ]\nA 3\n")
    with pytest.raises(ValueError):
        Topology.parse(f)
