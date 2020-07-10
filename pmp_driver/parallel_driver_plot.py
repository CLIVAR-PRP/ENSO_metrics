#!/usr/bin/env python

"""
Usage example:
1. First realization per model
./parallel_driver.py --mip cmip6 --exp historical --case_id v20200305 --modnames all --realization r1i1p1f1 --metricsCollection ENSO_perf
2. All realizations of individual models
./parallel_driver.py --mip cmip6 --exp historical --case_id v20200305 --modnames all --realization all --metricsCollection ENSO_perf

"""

from __future__ import print_function
#from argparse import RawTextHelpFormatter
#from genutil import StringConstructor
from subprocess import Popen

from PMPdriver_lib import AddParserArgument
from PMPdriver_lib import sort_human

import datetime
#import glob
import json
import os
#import pcmdi_metrics
import sys
import time

from os import makedirs as OS__makedirs
from os.path import exists as OSpath__exists
from os.path import join as OSpath__join

# To avoid below error
# OpenBLAS blas_thread_init: pthread_create failed for thread XX of 96: Resource temporarily unavailable
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# Must be done before any CDAT library is called.
# https://github.com/CDAT/cdat/issues/2213
if 'UVCDAT_ANONYMOUS_LOG' not in os.environ:
    os.environ['UVCDAT_ANONYMOUS_LOG'] = 'no'

# =================================================
# Collect user defined options
# -------------------------------------------------
param = AddParserArgument()

# Metrics Collection
metric_collection = param.metricsCollection

# Pre-defined options
mip = param.mip
exp = param.exp
print('mip:', mip)
print('exp:', exp)

# Check given model option
models = param.modnames
print('models:', models)

# Realizations
realization = param.realization
if ('all' in [r.lower() for r in realization]) or (realization == 'all'):
    realization = 'all'
print('realization: ', realization)

# case id
case_id = param.case_id
print('case_id:', case_id)

path_main = "/p/user_pub/pmp/pmp_results/pmp_v1.1.2"
path_in_json = OSpath__join(path_main, "metrics_results", "enso_metric", mip, exp, case_id, metric_collection)
path_out = OSpath__join(path_main, "graphics", "enso_metric", mip, exp, case_id, metric_collection)

pattern = "_".join([mip, exp, metric_collection, case_id])

# ---------------------------------------------------#
# Adjust model list if "ALL" given
# ---------------------------------------------------#
# read json file
filename_js = OSpath__join(path_in_json, pattern + "_allModels_allRuns.json")
print('filename_js:', filename_js)
with open(filename_js) as ff:
    data_json = json.load(ff)['RESULTS']['model']

# Include all models if conditioned
if ('all' in [m.lower() for m in models]) or (models == 'all'):
    models = sort_human(list(data_json.keys()))

print('models:', models)
print('number of models:', len(models))

# Debug
debug = param.debug
print('debug:', debug)

# =================================================
# Create output directories
# -------------------------------------------------
print("path_out:", path_out)
if not OSpath__exists(path_out):
    try:
        OS__makedirs(path_out)
    except:
        pass

# =================================================
# Generates list of command
# e.g.: python PMPdriver_plot.py --mip cmip5 --exp historical --metricsCollection ENSO_perf --modnames IPSL-CM5A-LR --realization r1i1p1 --case_id v20200305
# -------------------------------------------------
cmds_list = []
for model in models:
    print(' ----- model: ', model, ' ---------------------')
    if realization is "all":
        runs_list = sort_human(list(data_json[model].keys()))
        print('runs_list (all):', runs_list)
    else:
        runs_list = [realization]
    for run in runs_list:
        cmd = ['python', 'PMPdriver_plot.py',
               '--mip', mip, '--exp', exp, '--metricsCollection', metric_collection,
               '--case_id', case_id,
               '--modnames', model,
               '--realization', run]
        cmds_list.append(cmd)

for i, cmd in enumerate(cmds_list):
    print(i+1, ' '.join(cmd))

# =================================================
# Run subprocesses in parallel
# -------------------------------------------------
# log dir
log_dir = os.path.join("log", case_id, metric_collection)

if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# number of tasks to submit at the same time
#num_workers = 8 
num_workers = 10
#num_workers = 30
#num_workers = 25

print("Start : %s" % time.ctime())

# submit tasks and wait for subset of tasks to complete
procs_list = []
for p, cmd in enumerate(cmds_list):
    timenow = time.ctime()
    print(timenow, p, ' '.join(cmd))
    model = cmd[-3]
    run = cmd[-1]
    log_filename = '_'.join(['log_ensoPlot', metric_collection, mip, exp, model, run, case_id])
    log_file = os.path.join(log_dir, log_filename)
    with open(log_file+"_stdout.txt", "wb") as out, open(log_file+"_stderr.txt", "wb") as err:
        procs_list.append(Popen(cmd, stdout=out, stderr=err))
        time.sleep(1)
    if ((p > 0 and p % num_workers == 0) or (p == len(cmds_list)-1)):
        print('wait...')
        for proc in procs_list:
            proc.wait()
        print("Tasks end : %s" % time.ctime())
        procs_list = []

# tasks done
print("End : %s" % time.ctime())
sys.exit('DONE')
