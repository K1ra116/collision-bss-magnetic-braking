import numpy as np
from dataclasses import replace
from spindown_model import *

M, R, lam = 1.3*MSUN, 1.5*RSUN, 0.1
Md = 1e-10*MSUN/YR
print("== 1. numerical integration vs analytic (fixed R, B, Mdot; slow rotation; dipole fraction 1)")
for law, fn in (("matt", tsd_matt_slow), ("justham", tsd_justham)):
    p = Params(R0=R, R_ms=R, J0=1e48, B_ms=100.0, f_dip=1.0, law=law)
    o = evolve([p], t_end=5*GYR, n_steps=20000)
    J, t = o["J"][:, 0], o["t"]
    t_num = t[np.argmax(J < J[0]/np.e)]
    t_ana = fn(M, R, 100.0, Md, lam)
    print(f"  {law:8s} numeric {t_num/GYR:.4f} Gyr | analytic {t_ana/GYR:.4f} Gyr | err {abs(t_num-t_ana)/t_ana*100:.2f}%")
print(f"  Ryu+2025 estimate: {tsd_ryu(100.0, Md)/GYR:.4f} Gyr")

print("== 2. exponents (B x2, Mdot x4)")
for name, fn in (("matt", tsd_matt_slow), ("justham", tsd_justham)):
    a = fn(M, R, 100, Md, lam)
    print(f"  {name:8s} B exponent {np.log(fn(M,R,200,Md,lam)/a)/np.log(2):+.3f} | Mdot exponent {np.log(fn(M,R,100,4*Md,lam)/a)/np.log(4):+.3f}")

print("== 3. when the Justham form gives Ryu's 0.1 Gyr")
for Rs, lm in ((1.5, 0.1), (1.0, 0.1), (1.0, 0.2), (1.2, 0.08)):
    print(f"  R={Rs} Rsun lam={lm}: {tsd_justham(M, Rs*RSUN, 100, Md, lm)/GYR:.3f} Gyr")

print("== 4. angular momentum and mass bookkeeping (no braking, no coupling)")
p = Params(B_ms=1e-9)
o = evolve([p], t_end=1*GYR, n_steps=20000)
print(f"  final f = {o['f'][-1,0]:.4f} (must be <= 1), mass shed = {o['M_shed'][-1,0]/MSUN:.4f} Msun")
print(f"  M0 - M_final - wind = {(p.M0 - o['M'][-1,0] - p.Mdot*1*GYR)/MSUN:.4f} Msun (must equal mass shed)")

print("== 5. timestep convergence (coupled, B = 300 G, Mdot = 1e-13)")
p = Params(B_ms=300.0, Mdot=1e-13*MSUN/YR, couple_shed=True)
for n in (3000, 12000, 48000):
    o = evolve([p], t_end=1*GYR, n_steps=n)
    print(f"  n_steps={n:6d}: mass shed {o['M_shed'][-1,0]/MSUN:.5f} Msun, v_eq(1 Gyr) {o['veq'][-1,0]/1e5:.2f} km/s")

print("== 6. disc torque vs analytic Omega0/(1+t/t_d) (fixed R, B; no wind)")
p = Params(R0=R, R_ms=R, J0=0.5*lam*M*R**2*omega_crit(M, R), B_ms=1000.0, f_dip=1.0, Mdot=1e-30, t_disk=np.inf)
o = evolve([p], t_end=1*GYR, n_steps=10000)
Om0 = o["Omega"][0, 0]; td = tsd_disk(M, R, 1000.0, Om0, lam)
print(f"  t_d = {td/YR:.3e} yr")
for tt in (1e4, 1e6, 1e8):
    i = np.argmin(abs(o["t"] - tt*YR))
    num, ana = o["Omega"][i, 0]/Om0, 1/(1 + o["t"][i]/td)
    print(f"  t={tt:.0e} yr: numeric {num:.5f} | analytic {ana:.5f} | err {abs(num-ana)/ana*100:.3f}%")

print("== 7. timestep convergence (5 Myr disc, B = 1000 G, Mdot = 1e-13, coupled)")
p = Params(B_ms=1000.0, Mdot=1e-13*MSUN/YR, couple_shed=True, t_disk=5e6*YR)
for n in (3000, 12000, 48000):
    o = evolve([p], t_end=1*GYR, n_steps=n)
    print(f"  n_steps={n:6d}: mass shed {o['M_shed'][-1,0]/MSUN:.5f} Msun, v_eq(1 Gyr) {o['veq'][-1,0]/1e5:.2f} km/s")
