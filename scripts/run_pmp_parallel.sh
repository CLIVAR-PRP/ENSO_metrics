#!/bin/sh
set -a

# To avoid below error
# OpenBLAS blas_thread_init: pthread_create failed for thread XX of 96: Resource temporarily unavailable
export OMP_NUM_THREADS=1

# Working conda env in gates: cdat82_20191107_py27

case_id="v"`date +"%Y%m%d"`
#case_id="v20191115"
#case_id="v20191121"

#mips='cmip5 cmip6'
#mips='cmip5'
#mips='cmip6'
mips='CLIVAR_LE'

MCs='ENSO_perf ENSO_tel ENSO_proc'
#MCs='ENSO_perf'
#MCs='ENSO_tel ENSO_proc'
#MCs='ENSO_tel'
#MCs='ENSO_proc'
#MCs='test_perf test_proc test_tel'
#MCs='test_proc'

modnames='all'
realization='all'

mkdir -p log/$case_id

for mip in $mips; do
    for MC in $MCs; do
        echo $mip $MC $realization $case_id
        python -u ./parallel_driver.py -p my_Param_ENSO.py --mip $mip --case_id=$case_id --modnames $modnames --metricsCollection $MC --realization $realization >& log/$case_id/log_parallel.${mip}.${MC}.all.${case_id}.txt &
        disown
        sleep 1
    done
done
