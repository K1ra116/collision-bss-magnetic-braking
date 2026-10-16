import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mesa_history import FIGS, load, summary

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.linewidth": 0.8,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
})
TXT, TXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_B = "#2a78d6"
C_0 = "#eb6834"
BAND = ["#f3efe6", "#eef2f6", "#eef5ef"]

b0 = load("field_scan/B0")
b1 = load("field_scan/B1000")
for c in (b0, b1):
    c["R"] = 10 ** c["log_R"]
    c["Teff"] = 10 ** c["log_Teff"]
    c["L"] = 10 ** c["log_L"]
t40 = summary(b1)["t40"]


def smooth(c, key, dex=0.04):
    lt = np.log10(np.maximum(c["star_age"], 1e-9))
    edges = np.arange(-2.0, lt.max() + dex, dex)
    idx = np.digitize(lt, edges)
    ts, ys = [], []
    for k in np.unique(idx):
        m = idx == k
        if m.sum() == 0 or lt[m].mean() < -2:
            continue
        ts.append(10 ** np.median(lt[m])); ys.append(np.median(c[key][m]))
    return np.array(ts), np.array(ys)


def line(ax, c, key, col, lw=1.8, **kw):
    ax.plot(c["star_age"], c[key], color=col, lw=0.6, alpha=0.22)
    ax.plot(*smooth(c, key), color=col, lw=lw, **kw)


T0, T1 = 1e-2, 1.3e9
STAGES = [(T0, 1e2, "1. Bloated product\nreaches critical rotation"),
          (1e2, 3e7, "2. Kelvin-Helmholtz contraction\nsheds mass at critical rotation"),
          (3e7, T1, "3. Main sequence\nmagnetic braking")]


def at(c, key, t):
    return np.interp(np.log10(t), np.log10(np.maximum(c["star_age"], 1e-9)), c[key])


def style(ax, ylabel):
    ax.set_xscale("log")
    ax.set_xlim(T0, T1)
    for (a, b, _), col in zip(STAGES, BAND):
        ax.axvspan(a, b, color=col, zorder=0, lw=0)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_ylabel(ylabel, color=TXT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


fig = plt.figure(figsize=(11, 13.2))
gs = fig.add_gridspec(6, 1, height_ratios=[1.35, 1, 1.15, 0.85, 0.85, 0.9], hspace=0.12)
axs = [fig.add_subplot(gs[i]) for i in range(6)]
for a in axs[1:]:
    a.sharex(axs[0])

a = axs[0]
style(a, "")
a.set_yticks([])
a.spines["left"].set_visible(False)
a.set_ylim(-1, 1)
cmap = plt.get_cmap("RdYlBu")
snap = [1, 1e3, 1e5, 3e6, 3e7, 2.5e8, 9e8]
for t in snap:
    R = at(b1, "R", t)
    Te = at(b1, "Teff", t)
    v = at(b1, "surf_avg_v_rot", t)
    x = np.log10(t)
    rad = 0.05 + 0.042 * R
    col = cmap(np.clip((Te - 4500) / (8500 - 4500), 0, 1))
    ax_w = np.log10(T1) - np.log10(T0)
    fig_aspect = (a.get_position().height * fig.get_figheight()) / (a.get_position().width * fig.get_figwidth())
    wx = rad / 2 * ax_w * fig_aspect
    xs = np.linspace(x - wx, x + wx, 2)
    th = np.linspace(0, 2 * np.pi, 120)
    a.fill(10 ** (x + wx * np.cos(th)), 0.08 + rad * np.sin(th), color=col, ec="#52514e", lw=0.6, zorder=3, clip_on=False)
    a.text(t, -0.55, f"{R:.1f} Rsun\n{Te:,.0f} K\n{v:,.0f} km/s",
           ha="center", va="top", fontsize=7.5, color=TXT2, zorder=4)
for (lo, hi, lab) in STAGES:
    a.text(10 ** ((np.log10(lo) + np.log10(hi)) / 2), 0.98, lab, ha="center", va="top",
           fontsize=9.5, color=TXT, fontweight="bold")
a.text(T0 * 1.3, -1.02, "B = 1000 G star (size shows radius order only, colour shows Teff)",
       fontsize=7.5, color=TXT2, va="bottom")
a.set_ylim(-1.05, 1.0)

a = axs[1]
style(a, "radius [Rsun]")
line(a, b0, "R", C_0, label="no field (B = 0)")
line(a, b1, "R", C_B, label="fossil field B = 1000 G")
a.set_yscale("log")
a.set_yticks([1, 2, 5, 10]); a.set_yticklabels(["1", "2", "5", "10"])
a.annotate("bloated by shock heating\n(~9x the settled radius)", xy=(1, 10.4), xytext=(3e-2, 3.0),
           fontsize=8.5, color=TXT2, arrowprops=dict(arrowstyle="-", color="#8a8984", lw=0.7))
a.annotate("settles on the main sequence (~1.2 Rsun)", xy=(1e8, 1.22), xytext=(3e5, 3.4), fontsize=8.5, color=TXT2,
           arrowprops=dict(arrowstyle="-", color="#8a8984", lw=0.7))
a.legend(loc="upper right", fontsize=8.5, frameon=False, bbox_to_anchor=(1.0, 1.0))

a = axs[2]
style(a, "surface rotation [km/s]")
line(a, b0, "surf_avg_v_rot", C_0)
line(a, b1, "surf_avg_v_rot", C_B)
a.axhline(40, color="#8a8984", lw=0.9, ls=(0, (4, 3)))
a.text(T0 * 1.3, 46, "40 km/s (fast-rotator threshold, Ferraro+2023)", fontsize=8, color=TXT2)
a.set_ylim(0, 380)
a.annotate("B = 0: stays near 270 km/s", xy=(8e8, 268), xytext=(3e2, 330), fontsize=8.5, color=C_0,
           arrowprops=dict(arrowstyle="-", color=C_0, lw=0.7))
a.annotate(f"B = 1000 G: below 40 km/s\nafter {t40 / 1e6:.0f} Myr", xy=(t40, 40), xytext=(1.2e6, 105), fontsize=8.5,
           color=C_B, arrowprops=dict(arrowstyle="-", color=C_B, lw=0.7))
a.plot([t40], [40], "o", ms=6, color=C_B, mec="#fcfcfb", mew=1.5, zorder=5)

a = axs[3]
style(a, "omega / omega_crit")
line(a, b0, "surf_avg_omega_div_omega_crit", C_0)
line(a, b1, "surf_avg_omega_div_omega_crit", C_B)
a.axhline(1, color="#8a8984", lw=0.9, ls=(0, (4, 3)))
a.text(T0 * 1.3, 1.05, "critical (gas leaves the equator)", fontsize=8, color=TXT2)
a.set_ylim(0, 1.35)

a = axs[4]
style(a, "mass [Msun]")
a.plot(b0["star_age"], b0["star_mass"], color=C_0, lw=1.6)
a.plot(b1["star_age"], b1["star_mass"], color=C_B, lw=1.6)
a.set_ylim(1.0, 1.30)
a.text(2e8, 1.035, "1.03 (-0.24)", fontsize=8.5, color=C_0, va="bottom")
a.text(2e8, 1.185, "1.18 (-0.09)", fontsize=8.5, color=C_B, va="bottom")
a.text(T0 * 1.3, 1.285, "initial 1.27 Msun", fontsize=8, color=TXT2, va="center")

a = axs[5]
style(a, "surface field [G]")
line(a, b1, "mb_B_surf", C_B)
a.set_yscale("log"); a.set_ylim(5, 2000)
a.annotate("flux conservation: field grows as the star shrinks\n(B ~ 1/R^2, 10 G -> ~900 G)", xy=(1e5, 115), xytext=(3e-2, 400),
           fontsize=8.5, color=TXT2, arrowprops=dict(arrowstyle="-", color="#8a8984", lw=0.7))
a.set_xlabel("age since collision [yr]", color=TXT)
for ax in axs[:-1]:
    plt.setp(ax.get_xticklabels(), visible=False)

fig.suptitle("Life of the collision product: 0.7 + 0.6 Msun, x = 0.25, wind 1e-13 Msun/yr (MESA)",
             fontsize=13, color=TXT, y=0.995)
fig.text(0.5, 0.008, "Thick lines: median in log-age bins. Faint lines: raw MESA output (jitter from the shedding prescription). Both runs use the same shedding law (w0 = 0.70).",
         ha="center", fontsize=8, color=TXT2)
fig.subplots_adjust(left=0.09, right=0.98, top=0.965, bottom=0.055)
fig.savefig(FIGS / "fig17_life_timeline.png", dpi=150)

fig, ax = plt.subplots(figsize=(8.5, 7))
for c, col, lab in ((b0, C_0, "no field (B = 0)"), (b1, C_B, "fossil field B = 1000 G")):
    m = c["star_age"] > 1e-2
    ax.plot(c["Teff"][m], c["L"][m], color=col, lw=0.6, alpha=0.2)
    _, te = smooth(c, "Teff", 0.05); _, ll = smooth(c, "L", 0.05)
    ax.plot(te, ll, color=col, lw=1.8, label=lab)
marks = [(1, "1 yr"), (1e3, "1 kyr"), (1e5, "100 kyr"), (3e6, "3 Myr"), (3e7, "30 Myr"),
         (6e8, "600 Myr"), (1.2e9, "1.2 Gyr")]
for t, lab in marks:
    for c, col, dx in ((b1, C_B, 1.06), (b0, C_0, 0.94)):
        x, y = at(c, "Teff", t), at(c, "L", t)
        ax.plot(x, y, "o", ms=6, color=col, mec="#fcfcfb", mew=1.3, zorder=5)
    x, y = at(b1, "Teff", t), at(b1, "L", t)
    ax.text(x * 1.035, y * 1.08, lab, fontsize=8, color=TXT2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(9500, 4300)
ax.set_xticks([9000, 8000, 7000, 6000, 5000]); ax.set_xticklabels(["9000", "8000", "7000", "6000", "5000"])
ax.set_yticks([2, 5, 10, 20, 50]); ax.set_yticklabels(["2", "5", "10", "20", "50"])
ax.minorticks_off()
ax.set_xlabel("effective temperature [K]")
ax.set_ylabel("luminosity [Lsun]")
ax.grid(color=GRID, lw=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.text(4900, 60, "1. cool, bloated product\n(convective envelope, critical rotation)", fontsize=9, color=TXT, ha="center")
ax.text(6450, 14, "2. contracts and heats up,\nshedding mass at critical rotation", fontsize=9, color=TXT, ha="center")
ax.text(8300, 9.0, "3. main sequence (blue straggler)\nsettles at a different place\nbecause the remaining mass differs\n(1.18 vs 1.03 Msun)", fontsize=9, color=TXT, ha="center")
ax.legend(loc="lower left", frameon=False, fontsize=9, title="dots: same ages (labels on the blue track)", title_fontsize=8)
ax.set_title("Two tracks of the same collision product", fontsize=12, color=TXT)
fig.tight_layout()
fig.savefig(FIGS / "fig18_hr_track.png", dpi=150)
