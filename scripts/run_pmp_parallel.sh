#!/bin/sh
set -a

# To avoid below error
# OpenBLAS blas_thread_init: pthread_create failed for thread XX of 96: Resource temporarily unavailable
export OMP_NUM_THREADS=1

# Working conda env in gates: cdat82_20191107_py27

case_id="v"`date +"%Y%m%d"`
#case_id="v20191115"
case_id="v20191121"

mips='cmip5 cmip6'
#mips='cmip5'
#mips='cmip6'

MCs='ENSO_perf ENSO_tel ENSO_proc'
#MCs='ENSO_perf'
#MCs='ENSO_tel ENSO_proc'
#MCs='ENSO_tel'
#MCs='ENSO_proc'

modnames='all'

for mip in $mips; do
    if [ $mip == 'cmip5' ]; then
        realization='r1i1p1'
        #modnames="BNU-ESM CESM1-FASTCHEM CMCC-CM FGOALS-g2 HadCM3 HadGEM2-CC IPSL-CM5A-LR IPSL-CM5A-MR MIROC4h MIROC5 MPI-ESM-LR MPI-ESM-MR"
        #modnames="HadGEM2-CC HadGEM2-ES INMCM4 IPSL-CM5A-LR IPSL-CM5A-MR IPSL-CM5B-LR MIROC-ESM MIROC-ESM-CHEM MIROC4h MIROC5 MPI-ESM-LR MPI-ESM-MR MPI-ESM-P MRI-CGCM3 MRI-ESM1 NorESM1-M NorESM1-ME"
        modnames="FGOALS-s2"
    elif [ $mip == 'cmip6' ]; then
        realization='r1i1p1f1'
        #modnames="EC-Earth3 EC-Earth3-Veg"
        modnames="CNRM-CM6-1 CNRM-ESM2-1 CanESM5 HadGEM3-GC31-LL MIROC-ES2L NorCPM1 NorESM2-LM UKESM1-0-LL"
    fi
    for MC in $MCs; do
        echo $mip $MC $realization $case_id
        ./parallel_driver.py -p my_Param_ENSO.py --mip $mip --case_id=$case_id --modnames $modnames --metricsCollection $MC --realization $realization >& log_parallel.${mip}.${MC}.all.${case_id}.txt &
        disown
        sleep 1
    done
done
