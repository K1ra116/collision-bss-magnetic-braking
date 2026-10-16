import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from make_collision_product import *
import io, contextlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "parents"
OUT = ROOT / "data" / "product_x000"
OUT.mkdir(exist_ok=True)
with contextlib.redirect_stdout(io.StringIO()):
    p1, p2, s, xa, T = main(D / "profile_0.7.data", D / "parent_0.7.mod",
                            D / "profile_0.6.data", D / "parent_0.6.mod", prefix=f"{OUT}/")
B, O, GR = "#2a78d6", "#eb6834", "#1baf7a"
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
q = lambda p: p.m / p.M
qs = s["m"] / s["M"]
ax[0].semilogy(q(p1), p1.A, color=B, lw=1.6, label="parent 0.7 Msun")
ax[0].semilogy(q(p2), p2.A, color=O, lw=1.6, label="parent 0.6 Msun")
ax[0].semilogy(qs, s["A"], color=GR, lw=2.2, label="product 1.22 Msun")
ax[0].set_ylim(2e13, 3e16); ax[0].set_ylabel(r"entropic variable $A = P/\rho^{5/3}$ [cgs]")
ax[0].set_title("Shock heating raises A everywhere", loc="left", fontsize=11)
ax[1].plot(q(p1), p1.xa[:, 0], color=B, lw=1.6)
ax[1].plot(q(p2), p2.xa[:, 0], color=O, lw=1.6)
ax[1].plot(qs, s["xa"][:, 0], color=GR, lw=0.8, alpha=0.45, label="product, before mixing")
ax[1].plot(qs, xa[:, 0], color=GR, lw=2.2, label="product, after mixing")
ax[1].set_ylabel("hydrogen mass fraction"); ax[1].set_ylim(0.2, 0.8)
ax[1].set_title("The core is a blend of both parent cores", loc="left", fontsize=11)
ax[2].semilogy(p1.r / RSUN, p1.rho, color=B, lw=1.6)
ax[2].semilogy(p2.r / RSUN, p2.rho, color=O, lw=1.6)
ax[2].semilogy(s["r"] / RSUN, s["rho"], color=GR, lw=2.2)
ax[2].set_xlabel("radius [Rsun]"); ax[2].set_ylabel(r"density [g cm$^{-3}$]"); ax[2].set_ylim(1e-6, 1e3)
ax[2].set_title("The product is puffed up to ~3 Rsun", loc="left", fontsize=11)
for a in ax[:2]: a.set_xlabel("enclosed mass fraction m/M")
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25, lw=0.5)
ax[0].legend(frameon=False, fontsize=9, loc="upper left"); ax[1].legend(frameon=False, fontsize=9, loc="lower right")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "fig8_collision_product.png", dpi=160)
