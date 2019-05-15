#!/bin/sh
set -a

# Working conda env in Crunchy: pmp_nightly_20180830

ver=`date +"%Y%m%d"`

mips='cmip5 cmip6'

#MCs='ENSO_perf ENSO_tel ENSO_proc'
MCs='ENSO_tel'

for mip in $mips; do
    for MC in $MCs; do
        echo $mip $MC
        python PMPdriver_EnsoMetrics.py -p my_Param_${mip}_${MC}.py >& log.${mip}.${MC}.all.v${ver}.txt &
        disown
    done
done
