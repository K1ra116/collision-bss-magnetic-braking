import numpy as np
from make_collision_product import *

def collide(p1, p2, x, b1_head):
    r_p = x * (p1.R + p2.R)
    lq = np.log10(p1.M / p2.M)
    b1 = b1_head - C4 * x * lq
    b2 = b1_head + (C5 * x - C6) * lq
    f_L = mass_loss_fraction(p1, p2, r_p)
    s = merge(p1, p2, b1, f_L, b2=b2)
    J_tot = p1.M * p2.M * np.sqrt(2 * G * r_p / (p1.M + p2.M))
    s["J_tot"], s["J_r"] = J_tot, J_tot * (1 - C9 * f_L)
    s["f_L"], s["b1"], s["b2"], s["x"] = f_L, b1, b2, x
    R = s["r_out"][-1]
    s["R"] = R
    s["I"] = np.sum(2 / 3 * s["dm"] * s["r"] ** 2)
    s["lam"] = s["I"] / (s["M"] * R**2)
    s["Omega_c"] = np.sqrt(G * s["M"] / R**3)
    s["J_max_rigid"] = s["I"] * s["Omega_c"]
    return s


def j_profile(s):
    m, dm, M, r = s["m"], s["dm"], s["M"], s["r"]
    cs = np.sqrt(5 / 3 * s["P"] / s["rho"])
    jin = C10 * cs * r * (m / M) ** (1 / 3)
    g = np.sqrt(G * m * r) * (m / M) ** (1 / 3)
    mg = np.linspace(m[1], 0.999 * M, 4000)
    jin_g, g_g = np.interp(mg, m, jin), np.interp(mg, m, g)
    djin, dg = np.gradient(jin_g, mg), np.gradient(g_g, mg)

    def total(i):
        k2 = djin[i] / dg[i]
        k3 = jin_g[i] - k2 * g_g[i]
        j = np.where(m < mg[i], jin, k2 * g + k3)
        return np.sum(j * dm), k2, k3, j

    lo, hi = 5, len(mg) - 5
    if not (total(hi)[0] <= s["J_r"] <= total(lo)[0]):
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if total(mid)[0] > s["J_r"]:
            lo = mid
        else:
            hi = mid
    J, k2, k3, j = total(lo)
    return dict(k1=mg[lo] / M, k2=k2, k3=k3, j=j, jK=np.sqrt(G * m * r))


def report_j(p1, p2, b1_head, J_max_ms=4.1e50):
    print("\n x     k1     k2     surface j/j_Kep   max j/j_Kep   mass holding the excess over J_max(MS)")
    for x in (0.05, 0.1, 0.25, 0.5):
        s = collide(p1, p2, x, b1_head)
        jp = j_profile(s)
        if jp is None:
            print(f" {x:4.2f}  no solution for k1"); continue
        ratio = jp["j"] / jp["jK"]
        Jcum = np.cumsum(jp["j"] * s["dm"])
        Jcap = J_max_ms * (s["M"] / (1.2155 * MSUN)) ** 1.5
        m_star = np.interp(Jcap, Jcum, np.cumsum(s["dm"]))
        k = s["m"] > 0.02 * s["M"]
        print(f" {x:4.2f}  {jp['k1']:.3f}  {jp['k2']:.3f}      {ratio[-1]:.3f}           {ratio[k].max():.3f}        "
              f"{(s['M']-m_star)/MSUN:.3f} Msun ({(1-m_star/s['M'])*100:.1f}% of M), J_r/J_cap = {s['J_r']/Jcap:.1f}")


def write_product(p1, p2, x, b1_head, prefix):
    s = collide(p1, p2, x, b1_head)
    xa, alpha = mix_composition(s)
    T, mu = temperature(s["P"], s["rho"], xa, p1.species)
    T = T * s["tcorr"]
    write_relax_files(s, xa, T, prefix)
    jp = j_profile(s)
    if jp is not None:
        xq = 1.0 - s["m"] / s["M"]
        order = np.argsort(xq)
        with open(prefix + "angular_momentum.dat", "w") as f:
            f.write(f"{len(xq)}\n")
            for i in order:
                f.write(f"{xq[i]:.16e} {jp['j'][i]:.16e}\n")
    return s, xa, T, jp


def table(p1, p2, b1_head):
    print(f"b1(0) = {b1_head:.4f}")
    print(" x=rp/(R1+R2)  f_L    M_r    R[Rsun] r99   r90   lam     I[g cm2]   J_r       J_r/J_max(rigid)  E_r")
    for x in (0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0):
        s = collide(p1, p2, x, b1_head)
        q = np.cumsum(s["dm"]) / s["M"]
        r99, r90 = np.interp(.99, q, s["r_out"]) / RSUN, np.interp(.9, q, s["r_out"]) / RSUN
        print(f"  {x:5.2f}     {s['f_L']:.4f} {s['M']/MSUN:.4f} {s['R']/RSUN:7.2f} {r99:5.2f} {r90:5.2f} {s['lam']:.4f} "
              f"{s['I']:.3e} {s['J_r']:.3e}  {s['J_r']/s['J_max_rigid']:6.2f}   {s['E']:.3e}")


if __name__ == "__main__":
    import sys
    from pathlib import Path

    data = Path(__file__).resolve().parents[1] / "data"
    p1 = Parent(data / "parents/profile_0.7.data", data / "parents/parent_0.7.mod")
    p2 = Parent(data / "parents/profile_0.6.data", data / "parents/parent_0.6.mod")
    fL0 = mass_loss_fraction(p1, p2)
    b1_head = find_b1(p1, p2, fL0, (p1.energy() + p2.energy()) * (1 + C7 * fL0))
    if len(sys.argv) == 1:
        table(p1, p2, b1_head)
        report_j(p1, p2, b1_head)
    else:
        x = float(sys.argv[1])
        out = data / f"product_x{round(100 * x):03d}"
        out.mkdir(exist_ok=True)
        s, xa, T, jp = write_product(p1, p2, x, b1_head, f"{out}/")
        print(f"wrote composition.dat, entropy.dat, angular_momentum.dat to {out}")
        print(f"product mass {s['M']/MSUN:.4f} Msun, R {s['R']/RSUN:.2f} Rsun, J_r {s['J_r']:.3e}, central H {xa[0,0]:.3f}")
