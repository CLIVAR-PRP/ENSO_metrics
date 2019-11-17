#!/bin/sh
set -a

# Working conda env in gates: cdat82_20191107_py27

ver=`date +"%Y%m%d-%H%M"`
case_id="v"`date +"%Y%m%d"`
case_id="v20191115"

mips='cmip5 cmip6'
#mips='cmip5'
#mips='cmip6'

#MCs='ENSO_perf ENSO_tel ENSO_proc'
MCs='ENSO_perf'
#MCs='ENSO_tel'
#MCs='ENSO_proc'

modnames='all'

for mip in $mips; do
    if [ $mip == 'cmip5' ]; then
        realization='r1i1p1'
        modnames="BNU-ESM CESM1-FASTCHEM CMCC-CM FGOALS-g2 HadCM3 HadGEM2-CC IPSL-CM5A-LR IPSL-CM5A-MR MIROC4h MIROC5 MPI-ESM-LR MPI-ESM-MR"
    elif [ $mip == 'cmip6' ]; then
        realization='r1i1p1f1'
        modnames="EC-Earth3 EC-Earth3-Veg"
    fi
    for MC in $MCs; do
        echo $mip $MC $realization $case_id
        python parallel_driver.py -p my_Param_ENSO.py --mip $mip --case_id=$case_id --modnames $modnames --metricsCollection $MC --realization $realization >& log_parallel.${mip}.${MC}.all.v${ver}.txt &
        disown
    done
done
