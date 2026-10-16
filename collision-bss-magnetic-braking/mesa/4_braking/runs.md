# Commands for the runs in `data/runs`

Arguments: `bash run_mb.sh B_ms w0 wind max_dt t_B decouple f_lag`

| data folder | command |
|---|---|
| field_scan/B0 | `bash run_mb.sh 0 0.70 1d-13 0 0 0 0` |
| field_scan/B300 | `bash run_mb.sh 300 0.70 1d-13 3d6 0 0 0` |
| field_scan/B500 | `bash run_mb.sh 500 0.70 1d-13 3d6 0 0 0` |
| field_scan/B700 | `bash run_mb.sh 700 0.70 1d-13 3d6 0 0 0` |
| field_scan/B1000 | `bash run_mb.sh 1000 0.70 1d-13 3d6 0 0 0` |
| field_scan/B3000 | `bash run_mb.sh 3000 0.70 1d-13 3d6 0 0 0` |
| wind_scan/B1000_mdot1e-14 | `bash run_mb.sh 1000 0.70 1d-14 3d6 0 0 0.03` |
| wind_scan/B1000_mdot1e-13 | `bash run_mb.sh 1000 0.70 1d-13 3d6 0 0 0.03` |
| wind_scan/B1000_mdot1e-12 | `bash run_mb.sh 1000 0.70 1d-12 3d6 0 0 0.03` |
| wind_scan/B1000_mdot1e-11 | `bash run_mb.sh 1000 0.70 1d-11 1d6 0 0 0.03` |
| wind_scan/B200_mdot1e-10 | `bash run_mb.sh 200 0.70 1d-10 1d6 0 0 0.03` |
| threshold/B0_w085 | `bash run_mb.sh 0 0.85 1d-13 3d6 0 0 0.03` |
| threshold/B1000_w085 | `bash run_mb.sh 1000 0.85 1d-13 3d6 0 0 0.03` |
| decay/B1000_tB1e8 | `bash run_mb.sh 1000 0.70 1d-13 3d6 1d8 0 0.03` |
| decay/B1000_tB3e8 | `bash run_mb.sh 1000 0.70 1d-13 3d6 3d8 0 0.03` |
| decay/B1000_decoupled | `bash run_mb.sh 1000 0.70 1d-13 3d6 0 1 0.03` |
| decay/B3000_decoupled | `bash run_mb.sh 3000 0.70 1d-13 3d6 0 1 0.03` |

`f_lag = 0` is the original shedding prescription. The field scan was run with it; later runs use `f_lag = 0.03`, which reproduces the 1000 G run to within 1%.
The runs were originally launched with four separate scripts; `run_mb.sh` combines them and writes the same inlist settings.

Each history file was copied from `LOGS_mb_*/history.data` and compressed with gzip.
