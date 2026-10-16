import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mesa_history import FIGS, load, summary

runs = [("B = 0 (no braking)", "field_scan/B0", "#9a9a96", "-"),
        ("1000 G, field decays in 1e8 yr", "decay/B1000_tB1e8", "#86b6ef", "-"),
        ("1000 G, field decays in 3e8 yr", "decay/B1000_tB3e8", "#2a78d6", "-"),
        ("1000 G, no decay (baseline)", "field_scan/B1000", "#0d366b", "-"),
        ("1000 G, no decay, shed gas not held", "decay/B1000_decoupled", "#e08a1e", "--")]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for lab, p, col, ls in runs:
    c = load(p)
    s = summary(c)
    print(f"{lab:38s} shed {s['dM']:.3f}  t40 {s['t40']:.3e}  v_end {s['v_end']:6.1f}  M_end {s['M_end']:.3f}")
    t = c["star_age"]
    m = t > 1e6
    ax[0].plot(t[m] / 1e9, c["surf_avg_v_rot"][m], color=col, ls=ls, lw=2, label=lab)
    ax[1].plot(t[m] / 1e9, c["star_mass"][m], color=col, ls=ls, lw=2)
s = summary(load("decay/B3000_decoupled"))
print(f"{'3000 G, no decay, shed gas not held':38s} shed {s['dM']:.3f}  t40 {s['t40']:.3e}")
ax[0].axhline(40, color="#888", lw=0.8, ls=":")
ax[0].text(1.05, 45, "40 km/s", fontsize=8, color="#666")
ax[0].set_ylabel("surface rotation speed [km/s]")
ax[1].set_ylabel("stellar mass [Msun]")
for a in ax:
    a.set_xlabel("age [Gyr]")
    a.grid(True, color="#e5e5e5", lw=0.6)
    a.set_axisbelow(True)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
ax[0].legend(fontsize=8, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.98))
ax[0].set_title("Rotation", fontsize=11)
ax[1].set_title("Mass (shedding)", fontsize=11)
ax[0].set_ylim(0, 420)
fig.tight_layout()
fig.savefig(FIGS / "fig16_decay_decouple.png", dpi=150)
