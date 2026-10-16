import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from spindown_model import MSUN, RSUN, YR
from pathlib import Path

FIGS = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parents[1] / "data"
d = np.loadtxt(DATA / "mesa_track_122.dat")
t, R, lam, I = d.T
k = t > 1
B, O = "#2a78d6", "#eb6834"
ta = np.geomspace(1, 1.1e9, 400)
Ra = 1.5 + (6.0 - 1.5) * np.exp(-ta / 1e7)
Ia = 0.1 * 1.3 * MSUN * (Ra * RSUN) ** 2
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))
ax[0].semilogx(ta, Ra, color=O, lw=1.6, ls="--"); ax[0].semilogx(t[k], R[k], color=B, lw=2.2)
ax[0].set_ylabel("radius [Rsun]"); ax[0].set_title("Radius: swells to 5, settles at 1.1 Rsun", loc="left", fontsize=11)
ax[0].text(3, 6.15, "assumed (v0.1-v0.3)", color=O, fontsize=9); ax[0].text(3e3, 5.2, "MESA", color=B, fontsize=9)
ax[1].semilogx(ta, 0.1 + 0 * ta, color=O, lw=1.6, ls="--"); ax[1].semilogx(t[k], lam[k], color=B, lw=2.2)
ax[1].set_ylim(0, 0.11); ax[1].set_ylabel(r"$\lambda = I/(MR^2)$"); ax[1].set_title(r"$\lambda$ stays below half the assumed 0.1", loc="left", fontsize=11)
ax[2].loglog(ta, Ia, color=O, lw=1.6, ls="--"); ax[2].loglog(t[k], I[k], color=B, lw=2.2)
ax[2].set_ylabel(r"moment of inertia $I$ [g cm$^2$]"); ax[2].set_title("Final I is 4x smaller than assumed", loc="left", fontsize=11)
for a in ax:
    a.set_xlabel("time since collision [yr]"); a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25, lw=0.5)
fig.tight_layout(); fig.savefig(FIGS / "fig9_mesa_track.png", dpi=160)
