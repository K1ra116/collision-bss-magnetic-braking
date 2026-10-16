import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mesa_history import FIGS, load, summary

Bs = [0, 300, 500, 700, 1000, 3000]
cols = dict(zip(Bs, ["#c0392b", "#e08a1e", "#b8a400", "#6aa84f", "#1f5fa8", "#6a3d9a"]))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
out = []
for B in Bs:
    c = load(f"field_scan/B{B}")
    t, v = c["star_age"], c["surf_avg_v_rot"]
    lab = "B = 0" if B == 0 else f"{B} G"
    ax[0].semilogx(t, v, color=cols[B], label=lab, lw=1.2)
    ax[1].semilogx(t, c["star_mass"], color=cols[B], label=lab)
    s = summary(c)
    i1 = min(np.searchsorted(t, 1e9), len(t) - 1)
    out.append((B, s["dM"], s["t40"], s["tau"], v[i1], v[-1]))
    print(f"B {B:4d}  shed {s['dM']:.3f} Msun  t40 {s['t40']:.3e}  tau {s['tau']:.3e}  v(1 Gyr) {v[i1]:6.1f}  v_end {v[-1]:6.1f}")
o = np.array(out[1:])
slope = np.polyfit(np.log10(o[:, 0]), np.log10(o[:, 3]), 1)
print(f"tau ~ B^{slope[0]:.2f}")
for lab, col in (("v(1 Gyr)", 4), ("v_end", 5)):
    lv = np.log10(np.maximum(o[:, col], 1e-3))
    print(f"{lab} = 40 km/s at B ~ {10 ** np.interp(-np.log10(40), -lv, np.log10(o[:, 0])):.0f} G")
ax[2].loglog(o[:, 0], o[:, 3] / 1e9, "o-", color="k", label=f"MESA runs (slope {slope[0]:.2f})")
bb = np.array([200, 4000.])
ax[2].loglog(bb, 10 ** np.polyval(slope, 3) / 1e9 * (bb / 1000) ** -0.87, ":", color="gray",
             label="Matt+2012 torque alone: B^-0.87")
ax[2].loglog(bb, 0.1 * (bb / 100) ** -1.0 * (1e-13 / 1e-10) ** -0.5, "--", color="#c0392b",
             label="Ryu+2025 estimate (B_s = B_ms, 1e-13 Msun/yr)")
ax[2].set_xlabel("B_ms [G]")
ax[2].set_ylabel("spin-down e-folding time [Gyr]")
ax[2].legend(fontsize=8)
ax[0].axhline(40, ls="--", color="gray")
ax[0].text(2e2, 50, "40 km/s", color="gray", fontsize=8)
ax[0].set_ylabel("surface rotation speed [km/s]")
ax[1].set_ylabel("mass [Msun]")
for a in ax[:2]:
    a.set_xlabel("age since collision [yr]")
    a.set_xlim(1e2, 1.3e9)
    a.legend(fontsize=8)
fig.suptitle("Off-axis collision product (0.7+0.6 Msun, x=0.25): fossil-field scan\n"
             "wind 1e-13 Msun/yr, f_dip = 0.3, shedding threshold w0 = 0.70, weak internal coupling, max timestep 3 Myr")
fig.tight_layout()
fig.savefig(FIGS / "fig13_field_scan.png", dpi=140)
