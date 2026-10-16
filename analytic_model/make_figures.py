import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import replace
from spindown_model import *
from pathlib import Path

FIGS = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parents[1] / "data"

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "grid.color": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF})
MSY = MSUN / YR
base = Params()
M, R, lam = base.M0, base.R_ms, base.lam

md = np.geomspace(1e-16, 1e-9, 200)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axvspan(5e-14, 1e-13, color=GRID, lw=0)
ax.text(7e-14, 3e-3, "A/F-star models\n(Vick+2010)", fontsize=8, color=INK2, ha="center", va="bottom")
ax.axvline(1e-10, color=INK2, lw=1, ls=":")
ax.text(1.3e-10, 3e-3, "observed upper limit\n= Ryu+2025 reference", fontsize=8, color=INK2, va="bottom")
ax.axhline(1.0, color=INK2, lw=1, ls="--"); ax.text(1.2e-16, 1.25, "1 Gyr", fontsize=8, color=INK2)
series = [("Matt+2012 (slow rotation)", tsd_matt_slow(M, R, 100.0, md * MSY, lam) / GYR, C1),
          ("Justham+2006", tsd_justham(M, R, 100.0, md * MSY, lam) / GYR, C2),
          ("Ryu+2025 estimate", tsd_ryu(100.0, md * MSY) / GYR, C3)]
for name, y, c in series:
    ax.plot(md, y, color=c, lw=2, label=name)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-16, 1e-9); ax.set_ylim(2e-3, 3e2)
ax.set_xlabel("wind mass-loss rate Ṁ [M☉/yr]"); ax.set_ylabel("spin-down time [Gyr]")
ax.grid(True, which="major", lw=0.6); ax.legend(frameon=False, loc="upper right", fontsize=9)
ax.set_title("Spin-down time for a 100 G dipole (1.3 M☉, 1.5 R☉)", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(FIGS / "fig3_formula_comparison.png", dpi=160)

Bg = np.geomspace(1, 1e4, 60)
runs = {"unc13": [replace(base, B_ms=B, Mdot=1e-13 * MSY) for B in Bg],
        "cpl13": [replace(base, B_ms=B, Mdot=1e-13 * MSY, couple_shed=True) for B in Bg],
        "cpl10": [replace(base, B_ms=B, Mdot=1e-10 * MSY, couple_shed=True) for B in Bg]}
res = {k: evolve(v, t_end=1 * GYR, n_steps=12000) for k, v in runs.items()}
shed = {k: o["M_shed"][-1] / MSUN for k, o in res.items()}
veq = {k: o["veq"][-1] / 1e5 for k, o in res.items()}

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].plot(Bg, shed["unc13"], color=C2, lw=2, label="shed gas not coupled to field")
ax[0].plot(Bg, shed["cpl13"], color=C1, lw=2, label="shed gas coupled to field")
ax[0].set_ylabel("mass shed at breakup [M☉]"); ax[0].set_ylim(0, 0.1)
ax[0].set_title("A. Mass lost during contraction (Ṁ_wind = 1e-13 M☉/yr)", loc="left", fontsize=10)
ax[0].legend(frameon=False, loc="lower left", fontsize=9)
ax[1].plot(Bg, np.maximum(veq["cpl13"], 1), color=C1, lw=2)
ax[1].plot(Bg, np.maximum(veq["cpl10"], 1), color=C3, lw=2)
ax[1].text(1.9e3, 20, "Ṁ_wind = 1e-13 M☉/yr\n(A/F-star models)", fontsize=9, color=INK, ha="right", va="center")
ax[1].text(15, 110, "Ṁ_wind = 1e-10 M☉/yr\n(upper limit)", fontsize=9, color=INK, ha="left", va="center")
ax[1].axhspan(1, 10, color=GRID, lw=0); ax[1].text(60, 2.2, "below v sin i floor (10 km/s)", fontsize=8, color=INK2)
ax[1].axhline(40, color=INK2, lw=1, ls=":"); ax[1].text(1.3, 46, "'fast rotator' 40 km/s", fontsize=8, color=INK2)
ax[1].set_yscale("log"); ax[1].set_ylim(1, 600); ax[1].set_ylabel("equatorial v at 1 Gyr [km/s]")
ax[1].set_title("B. Rotation 1 Gyr after collision (shed gas coupled)", loc="left", fontsize=10)
for a in ax:
    a.set_xscale("log"); a.set_xlabel("surface field after contraction B_ms [G]"); a.grid(True, which="major", lw=0.6)
fig.text(0.01, 0.005, "Matt+2012 torque, dipole fraction 0.3, no field decay, J₀ = 2×10⁵¹ g cm²/s, R: 6 → 1.5 R☉ over 10 Myr.", fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIGS / "fig4_coupled_shedding.png", dpi=160)

nB, nM = 40, 40
Bgrid, Mgrid = np.geomspace(10, 1e5, nB), np.geomspace(1e-16, 1e-9, nM)
plist = [replace(base, B_ms=B, Mdot=m * MSY, couple_shed=True) for m in Mgrid for B in Bgrid]
o = evolve(plist, t_end=1 * GYR, n_steps=12000)
v = (o["veq"][-1] / 1e5).reshape(nM, nB)
fig, ax = plt.subplots(figsize=(7.2, 5))
X, Y = np.meshgrid(Bgrid, Mgrid)
pc = ax.pcolormesh(X, Y, np.log10(np.clip(v, 1, None)), cmap="Blues", shading="auto")
cs = ax.contour(X, Y, v, levels=[10, 40], colors=INK, linewidths=1)
ax.text(2.5e3, 2e-11, "slower than 10 km/s\n(below v sin i floor)", fontsize=9, color=INK, ha="center")
ax.text(60, 4e-15, "still a fast rotator (> 40 km/s)", fontsize=9, color="white")
ax.text(8e4, 7e-10, "black lines: 40 and 10 km/s", fontsize=8, color=INK2, ha="right", va="top")
ax.axhspan(5e-14, 1e-13, facecolor="none", edgecolor=C2, lw=1.5)
ax.text(12, 1.5e-13, "A/F-star models (Vick+2010)", fontsize=8, color=INK)
ax.axvspan(10, 1e4, ymin=0, ymax=0.015, color=C2); ax.text(12, 1.6e-16, "Ryu+2025 field range", fontsize=8, color=INK)
fig.colorbar(pc, ax=ax, label="log₁₀ equatorial v at 1 Gyr [km/s]")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("surface field after contraction B_ms [G]"); ax.set_ylabel("wind mass-loss rate Ṁ [M☉/yr]")
ax.set_title("Rotation 1 Gyr after collision: what field and wind are needed", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(FIGS / "fig5_B_Mdot_map.png", dpi=160)

def first_B(arr, cond):
    idx = np.where(cond(arr))[0]
    return Bg[idx[0]] if idx.size else None
print("fig3 @1e-13:", {n: float(np.interp(1e-13, md, y)) for n, y, _ in series})
print("fig3 @1e-10:", {n: float(np.interp(1e-10, md, y)) for n, y, _ in series})
print("uncoupled shed range:", shed["unc13"].min(), shed["unc13"].max())
print("B for coupled shed < 0.02:", first_B(shed["cpl13"], lambda a: a < 0.02), " < 0.01:", first_B(shed["cpl13"], lambda a: a < 0.01))
print("B for v<40 @1e-13:", first_B(veq["cpl13"], lambda a: a < 40), " v<10:", first_B(veq["cpl13"], lambda a: a < 10))
print("B for v<40 @1e-10:", first_B(veq["cpl10"], lambda a: a < 40), " v<10:", first_B(veq["cpl10"], lambda a: a < 10))
print("veq@1Gyr cpl13 at B=100,1000,10000:", np.interp([100, 1000, 1e4], Bg, veq["cpl13"]))
