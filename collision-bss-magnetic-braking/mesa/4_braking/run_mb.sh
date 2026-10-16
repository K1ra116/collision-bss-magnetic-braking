#!/bin/bash
# Usage: bash run_mb.sh B_ms [w0] [wind] [max_dt_yr] [t_B_yr] [decouple] [f_lag]
# Defaults: w0 = 0.70, wind = 1d-13, max_dt = 3d6, t_B = 0 (no decay), decouple = 0, f_lag = 0.03
set -e
B=${1:?give B_ms in gauss}
W0=${2:-0.70}
WIND=${3:-1d-13}
DT=${4:-3d6}
TB=${5:-0}
DC=${6:-0}
LAG=${7:-0.03}
TAG="B${B}_w${W0/./}_mdot${WIND}_dt${DT}_tB${TB}_dc${DC}_lag${LAG/./}"
[ -f product_rot_t0.mod ] || { echo "product_rot_t0.mod not found: run 3_rotation first"; exit 1; }
sed -e "s/TAG/${TAG}/g" -e "s/BMS_VALUE/${B}d0/" -e "s/W0_VALUE/${W0}d0/" -e "s/WIND_VALUE/${WIND}/" \
    -e "s/DT_VALUE/${DT}/" -e "s/TDECAY_VALUE/${TB}/" -e "s/DC_VALUE/${DC}d0/" -e "s/LAG_VALUE/${LAG}d0/" \
    inlist_mb_template > inlist_project
./rn 2>&1 | tee "log_mb_${TAG}.txt"
