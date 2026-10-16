from mesa_history import RUNS, load, summary

print(f"{'run':32s} {'dM_tot':>6s} {'shed end':>9s} {'t40':>9s} {'tau':>9s} {'v_end':>6s} {'M_end':>6s} {'Xc_end':>6s} {'t_end':>9s}")
for p in sorted(RUNS.glob("*/*")):
    s = summary(load(p.relative_to(RUNS)))
    print(f"{str(p.relative_to(RUNS)):32s} {s['dM']:6.3f} {s['t_shed_end']:9.2e} {s['t40']:9.2e} {s['tau']:9.2e} "
          f"{s['v_end']:6.1f} {s['M_end']:6.3f} {s['Xc_end']:6.3f} {s['t_end']:9.2e}")
