#!/usr/bin/env python

from __future__ import print_function
from genutil import StringConstructor
from pcmdi_metrics.variability_mode.lib import dict_merge

import copy
import glob
import json
import os

# -------------------------------
mip = 'cmip5'
mip = 'cmip6'
exp = 'historical'
#case_id = 'v20191204'
case_id = 'v20200113'
metricsCollection = 'ENSO_perf'
metricsCollection = 'ENSO_proc'
metricsCollection = 'ENSO_tel'
#metricsCollection = 'test_perf'
#metricsCollection = 'test_proc'
#metricsCollection = 'test_tel'
#pmprdir = '/work/lee1043/imsi/result_test'
pmprdir = '/p/user_pub/pmp/pmp_results/pmp_v1.1.2'
# -------------------------------

json_file_dir_template = os.path.join(
    pmprdir,
    '%(output_type)', 'enso_metric',
    '%(mip)', '%(exp)', '%(case_id)', '%(metricsCollection)')
json_file_dir_template = StringConstructor(json_file_dir_template)
json_file_dir = json_file_dir_template(
    output_type='metrics_results', mip=mip, exp=exp, case_id=case_id, metricsCollection=metricsCollection)

json_file_template = '_'.join(['%(mip)_%(exp)_%(metricsCollection)', '%(case_id)', '%(model)', '%(realization)'])
json_file_template = '%(mip)_%(exp)_%(metricsCollection)_%(case_id)_%(model)_%(realization)'
json_file_template = StringConstructor(json_file_template)

# Search for individual JSONs
json_files = sorted(glob.glob(
    os.path.join(
        json_file_dir,
        json_file_template(
            mip=mip, exp=exp, metricsCollection=metricsCollection, case_id=case_id, model='*', realization='*')+'.json')))

# Remove diveDown JSONs and previously generated merged JSONs if included
json_files_revised = copy.copy(json_files)
for j, json_file in enumerate(json_files):
    filename_component = json_file.split('/')[-1].split('.')[0].split('_')
    if 'diveDown' in filename_component:
        json_files_revised.remove(json_file)
    elif 'allModels' in filename_component:
        json_files_revised.remove(json_file)
    elif 'allRuns' in filename_component:
        json_files_revised.remove(json_file)

# Load individual JSON and merge to one big dictionary
for j, json_file in enumerate(json_files_revised):
    print(j, json_file)
    f = open(json_file)
    dict_tmp = json.loads(f.read())
    if j == 0:
        dict_final = dict_tmp.copy()
    else:
        dict_merge(dict_final, dict_tmp)
    f.close()

# Dump final dictionary to JSON
final_json_filename = json_file_template(
    mip=mip, exp=exp, metricsCollection=metricsCollection, case_id=case_id,
    model='allModels', realization='allRuns')+'.json'
final_json_file = os.path.join(json_file_dir, final_json_filename)

with open(final_json_file, 'w') as fp:
    json.dump(dict_final, fp, sort_keys=True, indent=4)
