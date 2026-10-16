import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mesa_history import FIGS, load, summary

runs = [(1000, 1e-14, "B1000_mdot1e-14"), (1000, 1e-13, "B1000_mdot1e-13"), (1000, 1e-12, "B1000_mdot1e-12"),
        (1000, 1e-11, "B1000_mdot1e-11"), (200, 1e-10, "B200_mdot1e-10")]
cols = ["#6a3d9a", "#1f5fa8", "#2e8b57", "#e08a1e", "#c0392b"]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
res = []
for (B, md, name), col in zip(runs, cols):
    c = load(f"wind_scan/{name}")
    t, v = c["star_age"], c["surf_avg_v_rot"]
    lab = f"{B} G, {md:.0e} Msun/yr"
    ax[0].semilogx(t, v, color=col, label=lab, lw=1.2)
    ax[1].semilogx(t, c["star_mass"], color=col, label=lab)
    s = summary(c)
    ryu = 0.1e9 * (B / 100.) ** -1 * (md / 1e-10) ** -0.5
    res.append((B, md, s["tau"], s["t40"], ryu))
    print(f"{lab:22s}  shed end {s['t_shed_end']:.2e}  t40 {s['t40']:.3e}  tau {s['tau']:.3e}  Ryu {ryu:.3e}  M_end {s['M_end']:.3f}")
r = np.array([x for x in res if x[0] == 1000])
slope = np.polyfit(np.log10(r[:, 1]), np.log10(r[:, 2]), 1)[0]
print(f"tau ~ Mdot^{slope:.2f} at 1000 G (Matt alone: {-(1 - 2 * 0.2177):.3f}, Ryu: -0.5)")
ax[2].loglog(r[:, 1], r[:, 2] / 1e9, "o-", color="k", label=f"MESA, B_ms = 1000 G (slope {slope:.2f})")
mm = np.array([3e-15, 3e-11])
ax[2].loglog(mm, 0.1 * (1000 / 100.) ** -1 * (mm / 1e-10) ** -0.5, "--", color="#c0392b",
             label="Ryu+2025 estimate, B_s = 1000 G")
s = res[-1]
ax[2].loglog([s[1]], [s[2] / 1e9], "s", color="#c0392b", ms=8, label="MESA, 200 G (Sills+2005 choice)")
ax[2].loglog([s[1]], [s[4] / 1e9], "x", color="#c0392b", ms=9, label="Ryu+2025 estimate, 200 G")
ax[2].set_xlabel("steady wind [Msun/yr]")
ax[2].set_ylabel("spin-down e-folding time [Gyr]")
ax[2].legend(fontsize=8)
ax[0].axhline(40, ls="--", color="gray")
ax[0].set_ylabel("surface rotation speed [km/s]")
ax[1].set_ylabel("mass [Msun]")
for a in ax[:2]:
    a.set_xlabel("age since collision [yr]")
    a.set_xlim(1e2, 1.3e9)
    a.legend(fontsize=8)
fig.suptitle("Off-axis collision product (0.7+0.6 Msun, x=0.25): wind-strength scan\n"
             "f_dip = 0.3, shedding threshold w0 = 0.70 with lag 0.03, weak internal coupling")
fig.tight_layout()
fig.savefig(FIGS / "fig14_wind_scan.png", dpi=140)
