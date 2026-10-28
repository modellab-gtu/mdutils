"""Tests need LD_LIBRARY_PATH=/cta/capps/envs/ligparam/lib (system libstdc++ lacks CXXABI_1.3.15 for scipy)."""

def test_versions():
    import alchemlyb, pymbar, MDRestraintsGenerator, MDAnalysis
    assert alchemlyb.__version__.startswith("2.5")
    assert pymbar.__version__.startswith("4.0")
    assert MDAnalysis.__version__.startswith("2.11")


def test_package_importable():
    import automd.fe.abfe as a
    assert a.__version__ == "0.1.0"


def test_mdrg_runtime(prolig_rundir):
    """MDRestraintsGenerator 0.2.1 (pinned scipy<1.8 upstream) works with the
    installed scipy: ligand/host atom search + Boresch restraint on JZ4."""
    import MDAnalysis as mda
    from MDRestraintsGenerator.search import find_ligand_atoms, find_host_atoms
    from MDRestraintsGenerator.restraints import FindBoreschRestraint

    u = mda.Universe(str(prolig_rundir / "md_0.tpr"), str(prolig_rundir / "md_0.xtc"))
    u.trajectory[0]
    l_sel = "resname JZ4 and not name H*"
    l_atoms = find_ligand_atoms(u, l_selection=l_sel, p_align="protein and name CA")
    atom_set = []
    for la in l_atoms:
        for pa in find_host_atoms(u, la[0]):
            atom_set.append((la, pa))
    assert atom_set
    r = FindBoreschRestraint(u, atom_set)
    r.run(start=0, stop=10)
    rest = r.restraint
    atoms = set(rest.bond.atomgroup.indices)
    for obj in list(rest.angles) + list(rest.dihedrals):
        atoms |= set(obj.atomgroup.indices)
    assert len(atoms) == 6
