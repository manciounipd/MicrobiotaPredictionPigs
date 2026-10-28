# prediction_micro_v2 – pipeline

Only the scripts actually used to produce the results in
`/home/enrico/Work/prediction_micro_v2/`. Launch everything from this folder
(`./main.sh`): the `source("functions/...")` paths are relative to it.

Input data: `/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/`

| Step | Script | Uses | Output (in Work dir) |
|---|---|---|---|
| 1 | `01_prepare_data/prepare_data_1.R` | `functions/function.R` | `rdata/info_r1.R` |
| 1 | `01_prepare_data/prepare_data_2.R` | | `rdata/start<date>.RData` |
| 2 | `02_cross_validation/CV_across_system.R` | `functions/function_cv.R` | `CV_multi/` |
| 2 | `02_cross_validation/CV_within.R` | `functions/function_within.R` | `CV_withincage/` |
| 2 | `02_cross_validation/CV_sire.R` | `functions/function_CV_sire.R` | `CV_sire/` |
| 3 | `03_results/results_across_system1.R` | | `res/cv_Between.RData` |
| 3 | `03_results/results_within_system.R` | | `res/cvW.RData`, `rdata/within.RData` |
| 3 | `03_results/results_sire.R` | | `res/cvSire.RData` |
| 3 | `03_results/collect_results_nooMpred.R` | `res/*.RData` | `plot/Between*.png`, `plot/Within*.png`, `plot/Sire*.png` |
| 4 | `04_variance_components/semplice_VCE.R` | `rdata/start2025-11-20.RData` | `VCE_CB/`, `VCE_PB/`, `rdata/all_vceCB.RData`, `rdata/all_vcePB.RData` |
| 4 | `04_variance_components/save_vce.R` | `functions/function_plot.R` | `plot/Donut_Plots*.png`, `table/tab2.csv` |

Notes
- `CV_across_system.R` and `CV_sire.R` load the hard-coded `rdata/start2025-11-18.RData`
  (time point 3); `semplice_VCE.R` uses `start2025-11-20.RData` (time points 1 2 3).
  Change the date if you re-run step 1.
