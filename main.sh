#!/bin/bash
# dove sono script : questa cartella (pipeline/) -> LANCIO DA QUA
# dove sono i files: "/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/"
# dove lavoro : "/home/enrico/Work/prediction_micro_v2/"

cd "$(dirname "$0")"
WORK=/home/enrico/Work/prediction_micro_v2

run_all () {
    TP="$1"   # time point microbiota: "3" oppure "1 2 3"
    LOG="$2"

    echo " script for preparing data "
    Rscript 01_prepare_data/prepare_data_1.R "$TP" | tee $WORK/$LOG
    Rscript 01_prepare_data/prepare_data_2.R

    echo "Now performing Cross Validation analyses.. "
    Rscript 02_cross_validation/CV_across_system.R | tee $WORK/CV_across_system.log
    Rscript 02_cross_validation/CV_within.R        | tee $WORK/CV_within_system.log
    Rscript 02_cross_validation/CV_sire.R          | tee $WORK/CV_Sire_system.log

    echo "Collect results .. "
    Rscript 03_results/results_across_system1.R
    Rscript 03_results/results_sire.R
    Rscript 03_results/results_within_system.R
}

run_all "3"     log_3
## faccio la stessa corsa am con tutti i times point
run_all "1 2 3" log_123

echo "Plot CV results .. "
Rscript 03_results/collect_results_nooMpred.R

echo "Variance components .. "
Rscript 04_variance_components/semplice_VCE.R
Rscript 04_variance_components/save_vce.R
