#!/bin/bash
set -e
[ -f product_t0.mod ] || { echo "product_t0.mod not found: run 2_relax first"; exit 1; }

if [ "$(uname)" == "Darwin" ]; then SED=(sed -i ''); else SED=(sed -i); fi
cp "$MESA_DIR/star/defaults/history_columns.list" history_columns.list
cp "$MESA_DIR/star/defaults/profile_columns.list" profile_columns.list
for c in star_mdot kh_timescale surf_avg_j_rot surf_avg_omega surf_avg_omega_crit \
         surf_avg_omega_div_omega_crit surf_avg_v_rot rotational_mdot_boost \
         center_omega total_angular_momentum i_rot_total; do
  "${SED[@]}" -E "s/^([[:space:]]*)!(${c})([[:space:]]|\$)/\1\2\3/" history_columns.list
done
for c in radius dq omega i_rot j_rot v_rot omega_crit omega_div_omega_crit am_log_nu_omega; do
  "${SED[@]}" -E "s/^([[:space:]]*)!(${c})([[:space:]]|\$)/\1\2\3/" profile_columns.list
done

cp inlist_rot_relax inlist_project
./rn
