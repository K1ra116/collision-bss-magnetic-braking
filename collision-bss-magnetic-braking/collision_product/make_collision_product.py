import sys
import numpy as np

G = 6.67430e-8
MSUN = 1.9884098706980504e33
RSUN = 6.957e10
KB = 1.380649e-16
MU_AMU = 1.66053906660e-24
A_RAD = 7.565733250e-15

C1, C3, C6, C7, C8 = 0.157, -1.1, 2.5, 2.5, 5.0
C2, C4, C5, C9, C10 = 1.8, 0.5, 5.0, 2.0, 0.6

P_SURF = 1.0e4

T_MIN_ENTROPY_FILE = 2.0e5

SPECIES_AZ = {"h1": (1, 1), "he3": (3, 2), "he4": (4, 2), "c12": (12, 6),
              "n14": (14, 7), "o16": (16, 8), "ne20": (20, 10), "mg24": (24, 12)}


def read_profile(path):
    with open(path) as f:
        lines = f.readlines()
    names = lines[5].split()
    data = np.array([[float(x) for x in l.split()] for l in lines[6:] if l.strip()])
    d = {n: data[::-1, i] for i, n in enumerate(names)}
    return d


def read_mod(path):
    species, rows, mass = None, [], None
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2 and p[0] == "M/Msun":
                mass = float(p[1].replace("D", "E")) * MSUN
            if len(p) > 6 and p[0] == "lnd":
                species = p[6:]
            elif species and len(p) == 7 + len(species) and p[0].isdigit():
                rows.append([float(x.replace("D", "E")) for x in p[1:]])
    arr = np.array(rows)[::-1]
    return species, arr[:, 6:], arr[:, 4], mass


class Parent:
    def __init__(self, profile_path, mod_path):
        pr = read_profile(profile_path)
        self.species, self.xa, dq, self.M = read_mod(mod_path)
        n = len(dq)
        assert len(pr["mass"]) == n, "profile and .mod have different zone counts"
        self.dm = dq * self.M
        m_out = np.cumsum(self.dm)
        self.m = m_out - 0.5 * self.dm
        r_out = 10 ** pr["logR"] * RSUN
        r_in = np.concatenate([[0.0], r_out[:-1]])
        self.r_out = r_out
        self.r = (0.5 * (r_out**3 + r_in**3)) ** (1 / 3)
        self.P = 10 ** pr["logP"]
        self.rho = 10 ** pr["logRho"]
        self.T = 10 ** pr["logT"]
        self.A = np.maximum.accumulate(self.P / self.rho ** (5 / 3))
        self.R = r_out[-1]
        self.tcorr = self.T / temperature(self.P, self.rho, self.xa, self.species)[0]

    def energy(self):
        return np.sum(self.dm * (1.5 * self.P / self.rho - G * self.m / self.r))

    def radius_enclosing(self, frac):
        return np.interp(frac * self.M, np.cumsum(self.dm), self.r_out)


def mass_loss_fraction(p1, p2, r_p=0.0):
    mu = p1.M * p2.M / (p1.M + p2.M)
    return (C1 * mu / (p1.M + p2.M)
            * (p1.radius_enclosing(0.86) + p2.radius_enclosing(0.86))
            / (p1.radius_enclosing(0.5) + p2.radius_enclosing(0.5) + C2 * r_p))


def shocked_A(p, b):
    return p.A + 10.0 ** b * p.P ** C3


def hydrostatic(A, dm, Pc, P_surf):
    n = len(A)
    r_out = np.empty(n)
    Pcell = np.empty(n)
    r3, m, P = 0.0, 0.0, Pc
    four_pi = 4 * np.pi
    for i in range(n):
        dmi, Ai = dm[i], A[i]
        rho = (P / Ai) ** 0.6
        r3_mid = r3 + 0.5 * 3 * dmi / (four_pi * rho)
        dP = G * (m + 0.5 * dmi) * dmi / (four_pi * r3_mid ** (4 / 3))
        Pmid = P - 0.5 * dP
        if Pmid <= P_surf:
            return False, r_out, Pcell
        rho = (Pmid / Ai) ** 0.6
        r3_mid = r3 + 0.5 * 3 * dmi / (four_pi * rho)
        dP = G * (m + 0.5 * dmi) * dmi / (four_pi * r3_mid ** (4 / 3))
        Pmid = P - 0.5 * dP
        if Pmid <= P_surf or P - dP <= P_surf:
            return False, r_out, Pcell
        Pcell[i] = Pmid
        r3 += 3 * dmi / (four_pi * rho)
        r_out[i] = r3 ** (1 / 3)
        m += dmi
        P -= dP
    return True, r_out, Pcell


def solve_structure(A, dm, P_surf=P_SURF):
    M = dm.sum()
    lo, hi = 1e10, 1e22
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        ok, r_out, P = hydrostatic(A, dm, mid, P_surf)
        if ok:
            hi = mid
            best = (r_out.copy(), P.copy())
        else:
            lo = mid
        if hi / lo - 1 < 1e-13:
            break
    r_out, P = best
    rho = (P / A) ** 0.6
    r_in = np.concatenate([[0.0], r_out[:-1]])
    r = (0.5 * (r_out**3 + r_in**3)) ** (1 / 3)
    m = np.cumsum(dm) - 0.5 * dm
    E = np.sum(dm * (1.5 * P / rho - G * m / r))
    return dict(r_out=r_out, r=r, m=m, P=P, rho=rho, E=E, M=M)


def retained_masks(p1, p2, A1, A2, m_lost):
    def exterior_mass(p, Af, Acut):
        Amono = np.maximum.accumulate(Af)
        return p.dm[Amono > Acut].sum()

    lo, hi = min(A1.min(), A2.min()), max(A1.max(), A2.max())
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if exterior_mass(p1, A1, mid) + exterior_mass(p2, A2, mid) > m_lost:
            lo = mid
        else:
            hi = mid
    keep1 = np.maximum.accumulate(A1) <= hi
    keep2 = np.maximum.accumulate(A2) <= hi
    return keep1, keep2


def refine_surface(A, dm, xa, origin, tcorr, n_sub=40):
    w = 0.5 ** np.arange(1, n_sub + 1)
    w[-1] *= 2
    rep = lambda a: np.concatenate([a[:-1], np.repeat(a[-1:], n_sub, axis=0)])
    dm_new = np.concatenate([dm[:-1], dm[-1] * w])
    return rep(A), dm_new, rep(xa), rep(origin), rep(tcorr)


def merge(p1, p2, b1, f_L, heat=True, b2=None):
    if b2 is None:
        b2 = b1 - C6 * np.log10(p1.M / p2.M)
    A1 = shocked_A(p1, b1) if heat else p1.A.copy()
    A2 = shocked_A(p2, b2) if heat else p2.A.copy()
    k1, k2 = retained_masks(p1, p2, A1, A2, f_L * (p1.M + p2.M)) if f_L > 0 else (
        np.ones(len(A1), bool), np.ones(len(A2), bool))
    A = np.concatenate([A1[k1], A2[k2]])
    dm = np.concatenate([p1.dm[k1], p2.dm[k2]])
    xa = np.vstack([p1.xa[k1], p2.xa[k2]])
    tcorr = np.concatenate([p1.tcorr[k1], p2.tcorr[k2]])
    origin = np.concatenate([np.ones(k1.sum(), int), 2 * np.ones(k2.sum(), int)])
    order = np.argsort(A, kind="stable")
    A, dm, xa, origin, tcorr = refine_surface(A[order], dm[order], xa[order], origin[order], tcorr[order])
    s = solve_structure(A, dm)
    s.update(A=A, dm=dm, xa=xa, origin=origin, tcorr=tcorr,
             lost1=p1.dm[~k1].sum(), lost2=p2.dm[~k2].sum())
    return s


def mix_composition(s):
    m, dm, M = s["m"], s["dm"], s["M"]
    alpha = C8 * np.log(s["A"].max() / s["A"].min()) ** 2
    mi, mk = m[:, None], m[None, :]
    g = (np.exp(-alpha / M**2 * (mi - mk) ** 2)
         + np.exp(-alpha / M**2 * (mi + mk) ** 2)
         + np.exp(-alpha / M**2 * (mi + mk - 2 * M) ** 2))
    norm = (g * dm[None, :]).sum(axis=1)
    w = g * (dm / norm)[:, None]
    xa = w.T @ s["xa"]
    return xa / xa.sum(axis=1, keepdims=True), alpha


def temperature(P, rho, xa, species):
    inv_mu = np.zeros(len(P))
    for j, sp in enumerate(species):
        a, z = SPECIES_AZ[sp]
        inv_mu += xa[:, j] * (1 + z) / a
    T = P / (rho * KB * inv_mu / MU_AMU)
    for _ in range(60):
        f = rho * KB * T * inv_mu / MU_AMU + A_RAD * T**4 / 3 - P
        T -= f / (rho * KB * inv_mu / MU_AMU + 4 * A_RAD * T**3 / 3)
    return T, 1 / inv_mu


def find_b1(p1, p2, f_L, target_E):
    f = lambda b: merge(p1, p2, b, f_L)["E"] - target_E
    lo, hi = 20.0, 40.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise RuntimeError(f"b1 not bracketed: f({lo})={flo:.3e}, f({hi})={fhi:.3e}")
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm * flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    return 0.5 * (lo + hi)


def write_relax_files(s, xa, T, prefix=""):
    xq = 1.0 - s["m"] / s["M"]
    order = np.argsort(xq)
    with open(prefix + "composition.dat", "w") as f:
        f.write(f"{len(xq)} {xa.shape[1]}\n")
        for i in order:
            f.write(f"{xq[i]:.16e} " + " ".join(f"{x:.16e}" for x in xa[i]) + "\n")
    keep = [i for i in order if T[i] >= T_MIN_ENTROPY_FILE]
    with open(prefix + "entropy.dat", "w") as f:
        f.write(f"{len(keep)}\n")
        for i in keep:
            f.write(f"{xq[i]:.16e} {s['rho'][i]:.16e} {T[i]:.16e}\n")


def main(profile1, mod1, profile2, mod2, prefix="", mix=True):
    p1, p2 = Parent(profile1, mod1), Parent(profile2, mod2)
    if p2.M > p1.M:
        p1, p2 = p2, p1
    assert p1.species == p2.species
    f_L = mass_loss_fraction(p1, p2)
    E_par = p1.energy() + p2.energy()
    target = E_par * (1 + C7 * f_L)
    b1 = find_b1(p1, p2, f_L, target)
    s = merge(p1, p2, b1, f_L)
    xa, alpha = mix_composition(s)
    xa_unmixed = s["xa"]
    if not mix:
        xa = xa_unmixed
    T, mu = temperature(s["P"], s["rho"], xa, p1.species)
    T = T * s["tcorr"]
    write_relax_files(s, xa, T, prefix)

    lam = np.sum(2 / 3 * s["dm"] * s["r"] ** 2) / (s["M"] * s["r_out"][-1] ** 2)
    ih = p1.species.index("h1")
    print(f"parents           : {p1.M/MSUN:.3f} + {p2.M/MSUN:.3f} Msun")
    print(f"mass loss fraction: {f_L:.4f}  (lost {s['lost1']/MSUN:.4f} + {s['lost2']/MSUN:.4f} Msun)")
    print(f"product mass      : {s['M']/MSUN:.4f} Msun")
    print(f"shock heating b1  : {b1:.4f}   b2: {b1 - C6*np.log10(p1.M/p2.M):.4f}")
    print(f"energy  parents   : {E_par:.4e} erg   remnant: {s['E']:.4e} erg (target {target:.4e})")
    print(f"product radius    : {s['r_out'][-1]/RSUN:.3f} Rsun")
    print(f"radius with 99% M : {np.interp(0.99*s['M'], np.cumsum(s['dm']), s['r_out'])/RSUN:.3f} Rsun")
    print(f"radius with 90% M : {np.interp(0.90*s['M'], np.cumsum(s['dm']), s['r_out'])/RSUN:.3f} Rsun")
    print(f"central P, rho, T : {s['P'][0]:.3e}, {s['rho'][0]:.3e}, {T[0]:.3e}")
    print(f"central H1        : {xa[0, ih]:.4f}   surface H1: {xa[-1, ih]:.4f}"
          f"   (central H1 before mixing: {xa_unmixed[0, ih]:.4f})")
    print(f"lambda = I/(M R^2): {lam:.4f}")
    print(f"mixing alpha      : {alpha:.1f}")
    print(f"wrote {prefix}composition.dat ({len(T)} points) and {prefix}entropy.dat "
          f"({int((T >= T_MIN_ENTROPY_FILE).sum())} points)")
    return p1, p2, s, xa, T


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--no-mix"]
    if len(args) != 4:
        print("usage: python make_collision_product.py profile1 mod1 profile2 mod2 [--no-mix]")
        sys.exit(1)
    main(*args, mix="--no-mix" not in sys.argv)
