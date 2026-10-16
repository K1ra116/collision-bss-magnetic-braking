import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mesa_history import FIGS, load, summary

fit_runs = [(500, 1e-13, "field_scan/B500"), (700, 1e-13, "field_scan/B700"), (1000, 1e-13, "field_scan/B1000"),
            (3000, 1e-13, "field_scan/B3000"), (1000, 1e-12, "wind_scan/B1000_mdot1e-12"),
            (1000, 1e-11, "wind_scan/B1000_mdot1e-11")]
check_run = (200, 1e-10, "wind_scan/B200_mdot1e-10")
never_runs = [(300, 1e-13, "field_scan/B300"), (1000, 1e-14, "wind_scan/B1000_mdot1e-14")]

runs = [(B, md, summary(load(p))["t40"]) for B, md, p in fit_runs]
check = (check_run[0], check_run[1], summary(load(check_run[2]))["t40"])
for B, md, p in never_runs:
    assert np.isnan(summary(load(p))["t40"])

B, Md, t = map(np.array, zip(*runs))
X = np.column_stack([np.ones(len(B)), np.log10(B / 1000), np.log10(Md / 1e-13)])
c, a, b = np.linalg.lstsq(X, np.log10(t), rcond=None)[0]
pred = 10 ** (X @ [c, a, b])
print(f"t40 = {10 ** c:.3e} yr (B/1000 G)^{a:.2f} (Mdot/1e-13)^{b:.2f}")
print(f"max fit error {100 * np.max(np.abs(pred / t - 1)):.0f}%")
pc = 10 ** (c + a * np.log10(check[0] / 1000) + b * np.log10(check[1] / 1e-13))
print(f"held-out 200 G, 1e-10: predicted {pc:.2e}, MESA {check[2]:.2e}")


def B_needed(T, mdot):
    return 1000 * (T / (10 ** c * (mdot / 1e-13) ** b)) ** (1 / a)


for T in (0.3e9, 0.6e9, 1.0e9):
    print(f"T = {T:.1e} yr: B needed at " + ", ".join(f"{m:.0e} -> {B_needed(T, m):.0f} G"
                                                   for m in (1e-14, 1e-13, 1e-12, 1e-11)))

fig, ax = plt.subplots(figsize=(6.4, 4.6))
md = np.logspace(-14.3, -9.7, 100)
for T, col in ((1.0e9, "#86b6ef"), (0.6e9, "#3987e5"), (0.3e9, "#184f95")):
    ax.loglog(md, B_needed(T, md), color=col, lw=2)
    ax.text(md[-1] * 1.15, B_needed(T, md[-1]), f"{T / 1e9:.1f} Gyr", color="#333", fontsize=9, va="center")
ax.loglog(Md, B, "o", ms=8, color="#222", label="MESA runs: reached 40 km/s")
ax.loglog(check[1], check[0], "D", ms=8, mfc="white", mec="#222", label="held-out check (200 G, 1e-10)")
ax.loglog([n[1] for n in never_runs], [n[0] for n in never_runs], "x", ms=9, mew=2, color="#888",
          label="never below 40 km/s by 1.2 Gyr")
ax.set_xlabel("wind mass-loss rate  Mdot  [Msun/yr]")
ax.set_ylabel("settled surface field  B  [G]  (dipole = 0.3 B)")
ax.set_title("Field and wind needed to spin down below 40 km/s", fontsize=11)
ax.grid(True, which="major", color="#ddd", lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.set_ylim(3, 2e4)
ax.set_xlim(3e-15, 5e-10)
ax.fill_between(md, 3, B_needed(1.0e9, md), color="#f0efec", zorder=0)
ax.text(1e-14, 15, "too weak: still fast\nafter 1 Gyr", fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIGS / "fig15_constraint.png", dpi=150)
