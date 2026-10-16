import gzip
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "data" / "runs"
FIGS = ROOT / "figures"
YR = 3.156e7


def load(name):
    path = RUNS / name / "history.data.gz"
    with gzip.open(path, "rt") as f:
        lines = f.read().split("\n")
    names = lines[5].split()
    d = np.genfromtxt(lines[6:])
    mn = d[:, names.index("model_number")]
    keep = np.ones(len(d), bool)
    last = np.inf
    for i in range(len(d) - 1, -1, -1):
        if mn[i] >= last:
            keep[i] = False
        else:
            last = mn[i]
    d = d[keep]
    return {n: d[:, i] for i, n in enumerate(names)}


def summary(c):
    t = c["star_age"]
    v = c["surf_avg_v_rot"]
    shed = np.where(c["mb_mdot_rot"] > 0)[0]
    t_shed = t[shed[-1]] if len(shed) else 0.0
    m = (t > max(t_shed, 1e6)) & (v < 200) & (v > 30)
    tau = -1 / np.polyfit(t[m], np.log(v[m]), 1)[0] if m.sum() > 5 else np.nan
    below = np.where((v < 40) & (t > 3e7))[0]
    t40 = t[below[0]] if len(below) else np.nan
    return dict(dM=c["star_mass"][0] - c["star_mass"][-1], t_shed_end=t_shed, t40=t40, tau=tau,
                v_end=v[-1], M_end=c["star_mass"][-1], Xc_end=c["center_h1"][-1], t_end=t[-1])
