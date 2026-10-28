"""Lambda schedules for ABFE legs (restraint, Coulomb, vdW components)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

SCHEDULE_DIR = Path(__file__).parent / "schedules"


@dataclass
class Schedule:
    restraint: list[float]
    coul: list[float]
    vdw: list[float]

    @property
    def n_states(self) -> int:
        return len(self.vdw)


def validate(s: Schedule) -> None:
    n = len(s.restraint)
    if not (len(s.coul) == n and len(s.vdw) == n) or n == 0:
        raise ValueError("restraint/coul/vdw lists must be non-empty and equal length")
    for name in ("restraint", "coul", "vdw"):
        v = getattr(s, name)
        if any(x < 0 or x > 1 for x in v):
            raise ValueError(f"{name} lambdas must lie in [0, 1]")
        if any(b < a for a, b in zip(v, v[1:])):
            raise ValueError(f"{name} lambdas must be non-decreasing")


def _expand(restraint, coul, vdw) -> Schedule:
    """Ramp restraint, then Coulomb, then vdW; hold others at end values."""
    states: list[tuple[float, float, float]] = []
    r_end = restraint[-1] if restraint else 0.0
    c_end = coul[-1] if coul else 0.0
    for r in restraint:
        states.append((r, 0.0, 0.0))
    for c in coul:
        states.append((r_end, c, 0.0))
    for v in vdw:
        states.append((r_end, c_end, v))
    out: list[tuple[float, float, float]] = []
    for st in states:
        if not out or st != out[-1]:
            out.append(st)
    s = Schedule([float(x[0]) for x in out], [float(x[1]) for x in out],
                 [float(x[2]) for x in out])
    validate(s)
    return s


def load_schedule(path, leg: Literal["complex", "solvent"]) -> Schedule:
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict) or leg not in data:
        raise ValueError(f"{path}: no '{leg}' section")
    sec = data[leg]
    return _expand(list(sec.get("restraint", [0])), list(sec.get("coul", [])),
                   list(sec.get("vdw", [])))


def default_schedule(leg: Literal["complex", "solvent"]) -> Schedule:
    return load_schedule(SCHEDULE_DIR / "default.yaml", leg)
