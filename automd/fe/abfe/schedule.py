"""Lambda schedules for ABFE legs (restraint, Coulomb, vdW components)."""
from __future__ import annotations

import math
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


def validate(s: Schedule, leg: str | None = None) -> None:
    """Check a schedule; raise ValueError on anything that would give a wrong dG.

    State 0 must be fully coupled, the last state fully decoupled (coul = vdw = 1);
    the complex restraint must end at 1, the solvent restraint must stay 0.
    """
    n = len(s.restraint)
    if not (len(s.coul) == n and len(s.vdw) == n) or n == 0:
        raise ValueError("restraint/coul/vdw lists must be non-empty and equal length")
    for name in ("restraint", "coul", "vdw"):
        v = getattr(s, name)
        if any((not isinstance(x, (int, float))) or not math.isfinite(x) for x in v):
            raise ValueError(f"{name} lambdas must be finite numbers")
        if any(x < 0 or x > 1 for x in v):
            raise ValueError(f"{name} lambdas must lie in [0, 1]")
        if any(b < a for a, b in zip(v, v[1:])):
            raise ValueError(f"{name} lambdas must be non-decreasing")
    if (s.restraint[0], s.coul[0], s.vdw[0]) != (0, 0, 0):
        raise ValueError("state 0 must be fully coupled (restraint = coul = vdw = 0)")
    if s.coul[-1] != 1 or s.vdw[-1] != 1:
        raise ValueError("the last state must be fully decoupled (coul = vdw = 1)")
    if leg == "complex" and s.restraint[-1] != 1:
        raise ValueError("complex leg: the restraint lambdas must end at 1")
    if leg == "solvent" and any(x != 0 for x in s.restraint):
        raise ValueError("solvent leg: the restraint lambdas must all be 0")


def _expand(restraint, coul, vdw, leg=None) -> Schedule:
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
    validate(s, leg)
    return s


def load_schedule(path, leg: Literal["complex", "solvent"]) -> Schedule:
    try:
        data = yaml.safe_load(Path(path).read_text())
    except yaml.YAMLError as e:
        raise ValueError(f"{path}: malformed YAML: {e}") from e
    if not isinstance(data, dict) or leg not in data:
        raise ValueError(f"{path}: no '{leg}' section")
    sec = data[leg]
    if not isinstance(sec, dict) or not sec:
        raise ValueError(f"{path}: '{leg}' section is empty or not a mapping")
    try:
        lists = [sec.get("restraint", [0]), sec.get("coul", []), sec.get("vdw", [])]
        if any(not isinstance(x, list) for x in lists):
            raise ValueError("restraint/coul/vdw must be lists")
        return _expand(*lists, leg=leg)
    except ValueError as e:
        raise ValueError(f"{path} ({leg}): {e}") from e


def default_schedule(leg: Literal["complex", "solvent"]) -> Schedule:
    return load_schedule(SCHEDULE_DIR / "default.yaml", leg)
