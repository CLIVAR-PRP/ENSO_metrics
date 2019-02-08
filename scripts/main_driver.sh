#!/bin/bash
source activate cdat8
# you must change the PYTHONPATH to match your own environment
export PYTHONPATH='/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages':'/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages/vtk':'/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages/pycf':'/home/yplanton/New_programs/Func':'/home/yplanton/New_programs/Tools':'/home/yplanton/New_programs/Lib':'/home/yplanton/New_programs/lib_cmip_bash':'/home/yplanton/ENSO_metrics/scripts':'/home/yplanton/ENSO_metrics/lib'

#---------------------------------------------------#
# This line can be modified
nbr_years=500 # define the number of model years used to compute the metric
#---------------------------------------------------#

# this loop is to ensure that the metric is computed for all possible slice of $nbr_years, even if a memory error or
# something occurs
for num in $(seq 0 20)
do
    python /home/yplanton/ENSO_metrics/scripts/main_driver.py ${nbr_years}
done
