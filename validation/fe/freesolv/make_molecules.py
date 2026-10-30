#!/usr/bin/env python
"""Build molecules.csv from FreeSolv database.json (MobleyLab/FreeSolv).

usage: make_molecules.py database.json molecules.csv

Matches the 16 validation molecules by exact (case-insensitive) iupac or nickname;
SMILES, experimental dG_hyd and its uncertainty are copied from the database, never typed.
The database file is untrusted data: it is only parsed with json.
"""
import csv
import json
import sys

# (residue name, label, accepted FreeSolv names)
MOLECULES = [
    ("F01", "methanol", ["methanol"]),
    ("F02", "ethanol", ["ethanol"]),
    ("F03", "1-butanol", ["butan-1-ol", "1-butanol"]),
    ("F04", "acetone", ["acetone", "propan-2-one"]),
    ("F05", "ethyl acetate", ["ethyl acetate"]),
    ("F06", "dimethyl ether", ["dimethyl ether", "methoxymethane"]),
    ("F07", "benzene", ["benzene"]),
    ("F08", "toluene", ["toluene"]),
    ("F09", "phenol", ["phenol"]),
    ("F10", "anisole", ["anisole"]),
    ("F11", "chlorobenzene", ["chlorobenzene"]),
    ("F12", "nitrobenzene", ["nitrobenzene"]),
    ("F13", "pyridine", ["pyridine"]),
    ("F14", "acetamide", ["acetamide"]),
    ("F15", "cyclohexane", ["cyclohexane"]),
    ("F16", "imidazole", ["imidazole", "1h-imidazole"]),
]


def main(db_path, out_path):
    with open(db_path) as fh:
        db = json.load(fh)
    rows, missing = [], []
    for res, label, names in MOLECULES:
        hits = [(k, v) for k, v in db.items()
                if any(str(v.get(f, "")).strip().lower() in names for f in ("iupac", "nickname"))]
        if len(hits) != 1:
            missing.append((label, [k for k, _ in hits]))
            continue
        k, v = hits[0]
        rows.append([res, label, k, v["iupac"].strip(), v["smiles"], v["expt"], v["d_expt"]])
    if missing:
        sys.exit(f"unmatched or ambiguous molecules: {missing}")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["res", "name", "freesolv_id", "iupac", "smiles", "expt_kcal", "d_expt_kcal"])
        w.writerows(rows)
    print(f"{len(rows)} molecules -> {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
