#!/bin/sh
set -ax

# Working conda env in Crunchy: pmp_nightly_20180830

ver=`date +"%Y%m%d"`

mip=cmip5
#mip=cmip6

for MC in 'ENSO_perf' 'ENSO_tel' 'ENSO_proc'; do
    python PMPdriver_EnsoMetrics.py -p my_Param_${mip}_${MC}.py >& log.${mip}.${MC}.all.v${ver}.txt &
    disown
done

