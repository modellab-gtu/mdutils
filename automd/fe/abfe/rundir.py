"""Read an mdutils MD run folder (topology, index, MDP, replica outputs)."""
from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from pathlib import Path


def read_mdp(path) -> dict[str, str]:
    """Parse an MDP file: lowercase keys, '_' -> '-', ';' comments stripped."""
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        key = key.strip().lower().replace("_", "-")
        if key:
            out[key] = val.strip()
    return out


def write_mdp(d: dict, path) -> None:
    """Write an MDP file with aligned ``key = value`` lines."""
    width = max((len(str(k)) for k in d), default=0)
    lines = [f"{str(k):<{width}} = {v}" for k, v in d.items()]
    Path(path).write_text("\n".join(lines) + "\n")


@dataclass
class RunDir:
    path: Path
    top: Path
    ndx: Path | None
    md_mdp: dict
    replica: int
    md_tpr: Path
    md_xtc: Path
    md_gro: Path
    temperature: float


def _find_top(path: Path) -> Path:
    sb = path / "run_md.sbatch"
    if sb.is_file():
        for line in sb.read_text().splitlines():
            if "grompp" in line and re.search(r"-o\s+md_", line):
                toks = shlex.split(line)
                if "-p" in toks and toks.index("-p") + 1 < len(toks):
                    return path / toks[toks.index("-p") + 1]
    return path / "topol.top"


def load_rundir(path, replica: int = 0) -> RunDir:
    path = Path(path)
    files = {ext: path / f"md_{replica}.{ext}" for ext in ("tpr", "xtc", "gro")}
    for p in files.values():
        if not p.is_file():
            raise FileNotFoundError(
                f"{p} not found (md_{replica}.* missing): run the MD replica first"
            )
    mdp = read_mdp(path / "MDP" / "md.mdp")
    ref_t = mdp.get("ref-t", "").split()
    if not ref_t:
        raise ValueError(f"no ref-t in {path / 'MDP' / 'md.mdp'}")
    ndx = path / "index.ndx"
    return RunDir(
        path=path,
        top=_find_top(path),
        ndx=ndx if ndx.is_file() else None,
        md_mdp=mdp,
        replica=replica,
        md_tpr=files["tpr"],
        md_xtc=files["xtc"],
        md_gro=files["gro"],
        temperature=float(ref_t[0]),
    )
