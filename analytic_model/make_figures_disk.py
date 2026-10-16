import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import replace
from spindown_model import *
from pathlib import Path

FIGS = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parents[1] / "data"

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "grid.color": GRID, "figure.facecolor": SURF, "axes.facecolor": SURF})
MSY = MSUN / YR
base = Params(Mdot=1e-13 * MSY)

tds = [0.0, 3e6, 1e7, 3e7]
names = ["no disc", "disc for 3 Myr", "disc for 10 Myr", "disc for 30 Myr"]
o = evolve([replace(base, B_ms=1000.0, t_disk=td * YR) for td in tds], t_end=1 * GYR, n_steps=20000)
t = o["t"][1:] / YR
fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
for j in range(4):
    ax[0].plot(t, o["f"][1:, j], color=COLS[j], lw=2, label=names[j])
    ax[1].plot(t, np.maximum(o["veq"][1:, j] / 1e5, 1), color=COLS[j], lw=2)
ax[0].axhline(1, color=INK2, lw=1, ls="--"); ax[0].text(1.5e3, 1.03, "breakup", fontsize=8, color=INK2)
ax[0].set_ylabel("Ω / Ω_crit"); ax[0].set_ylim(0, 1.15); ax[0].legend(frameon=False, loc="lower left", fontsize=9)
ax[1].axhspan(1, 10, color=GRID, lw=0); ax[1].text(1.5e3, 2.5, "below v sin i floor (10 km/s)", fontsize=8, color=INK2)
ax[1].axhline(40, color=INK2, lw=1, ls=":"); ax[1].text(1.5e3, 46, "'fast rotator' 40 km/s", fontsize=8, color=INK2)
ax[1].set_yscale("log"); ax[1].set_ylim(1, 600); ax[1].set_ylabel("equatorial v [km/s]")
ax[1].set_xscale("log"); ax[1].set_xlim(1e3, 1e9); ax[1].set_xlabel("time since collision [yr]")
for a in ax:
    a.grid(True, which="major", lw=0.6); a.axvline(1e7, color=INK2, lw=1, ls=":")
ax[0].text(8.5e6, 0.04, "contraction\ntime (10 Myr)", fontsize=8, color=INK2, ha="right")
fig.suptitle("Disc locking: spin evolution for B_ms = 1000 G, Ṁ_wind = 1e-13 M☉/yr", fontsize=11, x=0.02, ha="left")
fig.tight_layout(); fig.savefig(FIGS / "fig6_disk_time_evolution.png", dpi=160)
print("fig6 veq@1Gyr:", np.round(o["veq"][-1] / 1e5, 1), "shed:", np.round(o["M_shed"][-1] / MSUN, 4))
print("fig6 f max:", np.round(o["f"].max(axis=0), 3))

nB, nT = 44, 40
Bgrid, Tgrid = np.geomspace(10, 1e5, nB), np.geomspace(1e5, 1e9, nT)
plist = [replace(base, B_ms=B, t_disk=td * YR) for td in Tgrid for B in Bgrid]
o2 = evolve(plist, t_end=1 * GYR, n_steps=12000)
shed = (o2["M_shed"][-1] / MSUN).reshape(nT, nB)
veq = (o2["veq"][-1] / 1e5).reshape(nT, nB)
X, Y = np.meshgrid(Bgrid, Tgrid)
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
c0 = ax[0].pcolormesh(X, Y, shed, cmap="Blues", shading="auto")
ax[0].contour(X, Y, shed, levels=[0.02], colors=INK, linewidths=1)
fig.colorbar(c0, ax=ax[0], label="mass shed at breakup [M☉]")
ax[0].set_title("A. Mass lost during contraction", loc="left", fontsize=10)
c1 = ax[1].pcolormesh(X, Y, np.log10(np.clip(veq, 1, None)), cmap="Blues", shading="auto")
ax[1].contour(X, Y, veq, levels=[10, 40], colors=INK, linewidths=1)
fig.colorbar(c1, ax=ax[1], label="log₁₀ equatorial v at 1 Gyr [km/s]")
ax[1].set_title("B. Rotation 1 Gyr after collision", loc="left", fontsize=10)
for a in ax:
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("surface field after contraction B_ms [G]"); a.set_ylabel("disc lifetime [yr]")
    a.axhspan(1e6, 1e7, facecolor="none", edgecolor=COLS[1], lw=1.5)
    a.axvspan(10, 1e4, ymin=0, ymax=0.015, color=COLS[1])
fig.text(0.01, 0.005, "Orange box: disc lifetimes of 1–10 Myr (range used by Sills+2005). Orange tick: Ryu+2025 field range. Black lines: 0.02 M☉ (A); 40 and 10 km/s (B).\n"
         "Ṁ_wind = 1e-13 M☉/yr, shed gas not coupled, dipole fraction 0.3, disc inner edge at corotation, efficiency 1.", fontsize=8, color=INK2)
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(FIGS / "fig7_disk_map.png", dpi=160)

def minB(mask_row):
    idx = np.where(mask_row)[0]
    return round(float(Bgrid[idx[0]])) if idx.size else None
for td in (1e6, 3e6, 1e7, 3e7, 1e8):
    k = np.argmin(abs(Tgrid - td))
    print(f"t_disk={Tgrid[k]:.1e} yr: min B  shed<0.02: {minB(shed[k] < 0.02)}  v<40: {minB(veq[k] < 40)}  v<10: {minB(veq[k] < 10)}")
