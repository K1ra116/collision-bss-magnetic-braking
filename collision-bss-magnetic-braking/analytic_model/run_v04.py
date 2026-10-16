import numpy as np
from dataclasses import replace
import spindown_model as sm
from spindown_model import Params, evolve, MSUN, RSUN, YR, GYR, G
from pathlib import Path

FIGS = Path(__file__).resolve().parents[1] / "figures"
DATA = Path(__file__).resolve().parents[1] / "data"

track = sm.load_mesa_track(DATA / "mesa_track_122.dat")
T_END = 1.1 * GYR
base = Params(M0=1.2155 * MSUN, R_ms=1.10 * RSUN, J0=2.0e51, B_ms=100.0, Mdot=1e-13 * MSUN / YR)

def run(ps, **kw):
    return evolve(ps, t_end=T_END, n_steps=20000, track=track, **kw)

def at(o, key, t_yr, j=0):
    return np.interp(t_yr * YR, o["t"], o[key][:, j])

if __name__ == "__main__":
    for ty in (10, 1e3, 1e5, 1e7, 1e8, 1e9):
        R = np.interp(np.log(ty * YR), np.log(track["t"]), track["R"]); lam = np.interp(np.log(ty * YR), np.log(track["t"]), track["lam"])
        Jmax = lam * base.M0 * R**2 * np.sqrt(G * base.M0 / R**3)
        print(f"t={ty:8.0e} yr  R={R/RSUN:5.2f}  lam={lam:.4f}  J_max={Jmax:.2e}  v_crit={np.sqrt(G*base.M0/R)/1e5:.0f} km/s")
    print("\n-- no magnetic coupling of shed gas, wind only (Mdot=1e-13)")
    for J0 in (0.3e51, 0.5e51, 1e51, 2e51, 3e51):
        o = run([replace(base, J0=J0, B_ms=B) for B in (100.0, 1000.0)])
        print(f"J0={J0:.1e}: shed {o['M_shed'][-1,0]/MSUN:.3f} Msun | v(1 Gyr) B=100G: {at(o,'veq',1e9,0)/1e5:5.0f}  B=1kG: {at(o,'veq',1e9,1)/1e5:5.0f} km/s | f(1Gyr)={at(o,'f',1e9,0):.2f}")
    print("\n-- coupled shedding, J0=2e51")
    Bs = [100.0, 300.0, 1000.0, 3000.0]
    o = run([replace(base, B_ms=B, couple_shed=True) for B in Bs])
    for j, B in enumerate(Bs):
        print(f"B_ms={B:6.0f} G: shed {o['M_shed'][-1,j]/MSUN:.4f} Msun, v(1 Gyr)={at(o,'veq',1e9,j)/1e5:5.0f} km/s")
    print("\n-- disc locking + coupled shedding, J0=2e51: v(1 Gyr) [km/s]")
    for td in (0, 1e6, 1e7, 1e8):
        o = run([replace(base, B_ms=B, couple_shed=True, t_disk=td * YR) for B in Bs])
        print(f"t_disk={td:8.0e} yr: " + "  ".join(f"B={B:.0f}: {at(o,'veq',1e9,j)/1e5:5.0f}" for j, B in enumerate(Bs)))
    print("\n-- minimum B_ms for v(1 Gyr) < 40 km/s (coupled shedding)")
    Bgrid = np.geomspace(30, 3e4, 121)
    for td in (0, 1e6, 3e6, 1e7, 3e7, 1e8):
        for eff in (1.0, 0.1):
            o = run([replace(base, B_ms=B, couple_shed=True, t_disk=td * YR, disk_eff=eff) for B in Bgrid])
            v = np.array([at(o, 'veq', 1e9, j) for j in range(len(Bgrid))]) / 1e5
            ok = np.where(v < 40)[0]
            print(f"t_disk={td:8.0e} eff={eff}: B_min = {Bgrid[ok[0]]:.0f} G" if len(ok) else f"t_disk={td:8.0e} eff={eff}: none below 3e4 G")
            if td == 0: break
