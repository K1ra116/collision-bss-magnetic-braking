import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mesa_history import FIGS, load

runs = [("B = 0", load("field_scan/B0"), "#eb6834"),
        ("B_ms = 1000 G", load("field_scan/B1000"), "#2a78d6")]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for lab, c, col in runs:
    t = c["star_age"]
    ax[0].semilogx(t, c["surf_avg_v_rot"], color=col, label=lab, lw=1.2)
    ax[1].semilogx(t, c["star_mass"], color=col, label=lab)
    ax[2].loglog(t, c["total_angular_momentum"] / c["total_angular_momentum"][0], color=col, label=lab)
ax[0].axhline(40, ls="--", color="gray")
ax[0].text(2e2, 50, "40 km/s", color="gray", fontsize=8)
ax[0].set_ylabel("surface rotation speed [km/s]")
ax[1].set_ylabel("mass [Msun]")
ax[2].set_ylabel("J / J0")
ax[2].set_ylim(1e-3, 1.5)
for a in ax:
    a.set_xlabel("age since collision [yr]")
    a.set_xlim(1e2, 1.3e9)
    a.legend(fontsize=9)
fig.suptitle("Off-axis collision product (0.7+0.6 Msun, x=0.25): with and without a fossil field\n"
             "wind 1e-13 Msun/yr, shedding threshold w0 = 0.70 in both runs")
fig.tight_layout()
fig.savefig(FIGS / "fig12_B0_vs_B1000.png", dpi=140)
