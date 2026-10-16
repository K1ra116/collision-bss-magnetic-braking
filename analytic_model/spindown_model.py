import numpy as np
from dataclasses import dataclass

G = 6.674e-8
MSUN = 1.989e33
RSUN = 6.957e10
YR = 3.156e7
GYR = 1e9 * YR

K1, K2, M_EXP = 1.30, 0.0506, 0.2177


@dataclass
class Params:
    M0: float = 1.3 * MSUN
    R_ms: float = 1.5 * RSUN
    R0: float = 6.0 * RSUN
    t_KH: float = 1e7 * YR
    lam: float = 0.1
    J0: float = 2.0e51
    B_ms: float = 100.0
    f_dip: float = 0.3
    t_B: float = np.inf
    Mdot: float = 1e-10 * MSUN / YR
    law: str = "matt"
    couple_shed: bool = False
    t_disk: float = 0.0
    disk_eff: float = 1.0


def radius(t, R0, R_ms, t_KH):
    return R_ms + (R0 - R_ms) * np.exp(-t / t_KH)


def surface_field(t, R, B_ms, R_ms, t_B):
    return B_ms * (R_ms / R) ** 2 * np.exp(-t / t_B)


def omega_crit(M, R):
    return np.sqrt(G * M / R**3)


def lever2_matt(Omega, M, R, B_dip, Mdot):
    v_esc = np.sqrt(2 * G * M / R)
    v_mod = np.sqrt(K2**2 * v_esc**2 + (Omega * R) ** 2)
    return K1**2 * (B_dip**2 * R**2 / (Mdot * v_mod)) ** (2 * M_EXP)


def lever2_justham(Omega, M, R, B_dip, Mdot):
    return B_dip * R**1.25 * Mdot**-0.5 * (G * M) ** -0.25


def lever2(Omega, M, R, B_dip, Mdot, is_justham):
    Mdot = np.maximum(Mdot, 1.0)
    return np.where(is_justham,
                    lever2_justham(Omega, M, R, B_dip, Mdot),
                    lever2_matt(Omega, M, R, B_dip, Mdot))


def wind_torque(Omega, M, R, B_dip, Mdot, is_justham):
    return Mdot * Omega * R**2 * lever2(Omega, M, R, B_dip, Mdot, is_justham)


def disk_torque(Omega, M, R, B_dip, eff=1.0):
    mu = B_dip * R**3
    return eff * mu**2 * Omega**2 / (3 * G * M)


def evolve(params_list, t_end=2 * GYR, t_start=10 * YR, n_steps=30000, track=None):
    P = len(params_list)
    arr = lambda name: np.array([getattr(p, name) for p in params_list], dtype=float)
    M0, R_ms, R0, t_KH = arr("M0"), arr("R_ms"), arr("R0"), arr("t_KH")
    lam, J0, B_ms, f_dip = arr("lam"), arr("J0"), arr("B_ms"), arr("f_dip")
    t_B, Mdot = arr("t_B"), arr("Mdot")
    t_disk, disk_eff = arr("t_disk"), arr("disk_eff")
    is_j = np.array([p.law == "justham" for p in params_list])
    couple = np.array([p.couple_shed for p in params_list])

    ts = np.concatenate([[0.0], np.geomspace(t_start, t_end, n_steps)])
    keys = ("R", "M", "J", "Omega", "f", "B", "veq", "M_shed", "lam")
    out = {k: np.zeros((len(ts), P)) for k in keys}

    if track is not None:
        lt = np.log(np.maximum(track["t"], 1e-30))
        tr_R = lambda t: np.interp(np.log(max(t, 1e-30)), lt, track["R"])
        tr_lam = lambda t: np.interp(np.log(max(t, 1e-30)), lt, track["lam"])
    lam_const = lam

    M, J = M0.copy(), J0.copy()
    M_shed = np.zeros(P)
    for i, t in enumerate(ts):
        if track is not None:
            R = np.full(P, tr_R(t))
            lam = np.full(P, tr_lam(t))
        else:
            R = radius(t, R0, R_ms, t_KH)
        B = surface_field(t, R, B_ms, R_ms, t_B)
        dt = t - ts[i - 1] if i > 0 else t_start
        if i > 0:
            Omega = J / (lam * M * R**2)
            tau = wind_torque(Omega, M, R, f_dip * B, Mdot, is_j)
            J = J * np.exp(-dt * tau / J)
            M = M - Mdot * dt
            I = lam * M * R**2
            a = np.where(t <= t_disk, disk_eff * (f_dip * B * R**3) ** 2 / (3 * G * M * I**2), 0.0)
            J = J / (1.0 + a * J * dt)

        for _ in range(2):
            Oc = omega_crit(M, R)
            excess = np.clip(J - lam * M * R**2 * Oc, 0.0, None)
            dM = excess / ((1.0 - lam) * Oc * R**2)
            if couple.any() and excess.any():
                x = dM.copy()
                for _ in range(14):
                    lev = lever2(Oc, M, R, f_dip * B, x / dt, is_j)
                    lev = np.where(couple, np.maximum(lev, 1.0), 1.0)
                    x = excess / (Oc * R**2 * (lev - lam))
                dM = np.where(excess > 0, x, 0.0)
                carried = np.where(excess > 0, dM * Oc * R**2 * lev, 0.0)
            else:
                carried = dM * Oc * R**2
            M = M - dM
            M_shed = M_shed + dM
            J = J - carried

        Omega = J / (lam * M * R**2)
        for k, v in zip(keys, (R, M, J, Omega, Omega / omega_crit(M, R), B, Omega * R, M_shed, lam)):
            out[k][i] = v

    out["t"] = ts
    return out


def load_mesa_track(path):
    d = np.loadtxt(path)
    return {"t": d[:, 0] * YR, "R": d[:, 1] * RSUN, "lam": d[:, 2]}


def tsd_matt_slow(M, R, B_dip, Mdot, lam):
    return lam * M / (Mdot * lever2_matt(0.0, M, R, B_dip, Mdot))


def tsd_justham(M, R, B_dip, Mdot, lam):
    return lam * M / (Mdot * lever2_justham(0.0, M, R, B_dip, Mdot))


def tsd_ryu(B, Mdot):
    return 0.1 * GYR * (B / 100.0) ** -1 * (Mdot / (1e-10 * MSUN / YR)) ** -0.5


def tsd_disk(M, R, B_dip, Omega0, lam, eff=1.0):
    return 3 * G * lam * M**2 / (eff * B_dip**2 * R**4 * Omega0)
