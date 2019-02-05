#!/bin/bash
source arcivate cdat8
export PYTHONPATH='/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages':'/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages/vtk':'/home/yplanton/miniconda/envs/cdat8/lib/python2.7/site-packages/pycf':'/home/yplanton/New_programs/Func':'/home/yplanton/New_programs/Tools':'/home/yplanton/New_programs/Lib':'/home/yplanton/New_programs/lib_cmip_bash':'/home/yplanton/ENSO_metrics/scripts':'/home/yplanton/ENSO_metrics/lib'
python /home/yplanton/ENSO_metrics/scripts/main_driver.py
