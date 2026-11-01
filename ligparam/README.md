# ligparam: small-molecule parameters for Amber, GROMACS and OpenMM

```
module load ligparam/2026.10
ligparam -i lig.sdf                          # GAFF2 + ABCG2 charges (default)
ligparam -i lig.sdf --charge-model bcc       # GAFF2 + AM1-BCC
ligparam -i lig.sdf --charge-model nagl      # GAFF2 types + NAGL charges (graph-network AM1-BCC, no QM, instant)
ligparam -i lig.sdf --ff openff              # OpenFF Sage 2.2.1 + NAGL (or --charge-model bcc)
```
`run_examples.sh` here runs four ligands: aspirin (neutral), protonated lidocaine (+1), acetate (−1) and methyl
phosphate (−2).

## Input
- An **SDF or MOL2 with explicit hydrogens, 3D coordinates and correct formal charges**: the net charge is the sum of the
  formal charges. PDB files are refused because they carry no bond orders or charges.
- `--charge N` only cross-checks: ligparam stops if N disagrees with the file.
- `--protonate-ph 7.4` re-protonates with Open Babel first (this changes the charge; check the result).

## Geometry before the charges (`--opt`, GAFF2 route)
```
ligparam -i lig.sdf                          # --opt am1 (default): antechamber's own AM1 minimisation, as before
ligparam -i lig.sdf --opt aimnet2            # AIMNet2 relaxation first (also uma | ani2x | mace-off), then AM1 maxcyc=0
ligparam -i lig.sdf --opt hf                 # writes opt_hf.gjf + run_opt.sh (HF/6-31G*, cws02, 16 cores) and stops:
sbatch <outdir>/run_opt.sh                   #   then rerun with the optimised geometry:
ligparam -i lig.sdf --opt hf --geom <outdir>/opt_hf.log
ligparam -i lig.sdf --opt hf --local         # the same, Gaussian in the foreground
ligparam -i lig.sdf --opt none               # charges on the input geometry (AM1 maxcyc=0)
```
- With any `--opt` other than `am1`, AM1 runs with `maxcyc=0` on the given geometry, and the output coordinates are that
  geometry fitted (best rigid fit) onto the input pose. `ligparam.json` records `opt`.
- ANI-2x and MACE-OFF are neutral-only, ANI-2x covers H C N O S F Cl: a ligand they can't treat falls back to `hf` (said
  on screen and in `ligparam.json`).
- `--opt` matters for `abcg2`/`bcc`. For `nagl` and `user` charges the geometry only changes the output coordinates.
  `--ff openff` does not take `--opt`.

## Output (folder `<name>_<ff>_<charges>/`)
| File | Use |
|---|---|
| `LIG.prmtop`, `LIG.inpcrd` | Amber (and gmx_MMPBSA / MMPBSA.py) |
| `LIG.mol2`, `LIG.frcmod`, `LIG.lib` | tleap: `loadamberparams LIG.frcmod`, `loadoff LIG.lib` in your complex build (GAFF2 route) |
| `LIG.itp` | GROMACS: `[ atomtypes ]` + `[ moleculetype ]`, to `#include` in a protein topology |
| `LIG.top`, `LIG.gro` | GROMACS, ligand alone |
| `LIG_system.xml` | OpenMM System |
| `ligparam.json`, `ligparam.log` | method, charges, versions, every command run |

## Notes
- **Same AmberTools as `md/ambertools-25`** (25.3): antechamber/sqm/tleap and the GAFF2 files are byte-identical, so
  charges match the rest of the cluster.
- **parmchk2 runs with `-s 2` (GAFF2).** Its default is GAFF1. With GAFF2 atom types, the default writes GAFF1 values for
  GAFF2-only types (e.g. `nx`/`ny`/`nz`, protonated amines) into the frcmod, and these override gaff2.dat. Example:
  methylammonium C–N becomes 283.3 kcal/mol/Å², 1.511 Å instead of 238.8, 1.484 Å. If you run parmchk2 yourself, use
  `-s 2`: AmberTools 25 parmchk2 parses the word form `-s gaff2` unreliably (GAFF1 in about one run in four).
- Charge sums are made exact: antechamber's 6-decimal rounding (e.g. 0.999) is spread over the heavy atoms.
- The toolbox is also on PATH: `acpype`, `antechamber`, `parmchk2`, `tleap`, OpenFF toolkit/Interchange/NAGL,
  openmmforcefields, ParmEd, RDKit, Open Babel, MDAnalysis.
- Not included: RESP from Gaussian (use metalprep/MCPB or a manual G16 `Pop=MK` + `antechamber -c resp` run), CGenFF
  (licence), espaloma.

## Validation (GTU, 2026-10-28)
- Four ligands × four routes (GAFF2/ABCG2, GAFF2/AM1-BCC, GAFF2/NAGL, OpenFF/NAGL): every charge sum equals the formal
  charge, and the H count is unchanged.
- **GROMACS energies, two independent conversions:** acpype and ParmEd topologies of the same ligand give identical
  single-point energies, term by term (≤ 0.001 kJ/mol; 0.02 kJ/mol for the dianion, from acpype's own charge rounding).
- `LIG.itp` inside a hand-written topology gives the same energy as the full `LIG.top`.
