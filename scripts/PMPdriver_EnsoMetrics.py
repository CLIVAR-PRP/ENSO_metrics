#!/usr/bin/env python
#=================================================
# Dependencies
#-------------------------------------------------
from __future__ import print_function

import collections
import copy
import glob
import json
import os
import pcmdi_metrics
import sys

from collections import defaultdict
from genutil import StringConstructor
#from pcmdi_metrics.pcmdi.pmp_parser import PMPParser
from pcmdi_metrics.driver.pmp_parser import PMPParser
from pmpParser import ReadOptions
from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection

#=================================================
# Collect user defined options
#-------------------------------------------------
param = ReadOptions()

# Pre-defined options
mip = param.mip
exp = param.exp

# Path to model data as string template
#modpath = StringConstructor(param.modpath)
modpath = param.process_templated_argument("modpath")
#modpath_lf = StringConstructor(param.modpath_lf)
modpath_lf = param.process_templated_argument("modpath_lf")

# Check given model option
models = param.modnames

# List up all available models if 'all' given in models
if 'all' in [m.lower() for m in models]:
    models = [p.split('.')[1]
        for p in glob.glob(modpath(
            mip=mip,
            model='*',
            realization='*',
            variable='ts'))]
    # remove duplicates
    models = sorted(list(dict.fromkeys(models)), key=lambda s: s.lower())

print('models:', models)

# Output
outdir = param.process_templated_argument("results_dir")
netcdf_path = outdir(output_type='diagnostic_results')

mc_name = param.metricsCollection 
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

#=================================================
# Prepare loop iteration
#-------------------------------------------------
# Create output directory
for output_type in ['graphics', 'diagnostic_results', 'metrics_results']:
    if not os.path.exists(outdir(output_type=output_type)):
        os.makedirs(outdir(output_type=output_type))
    print(outdir(output_type=output_type))

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print(list_variables)

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
try:
    list_obs.remove('ERSSTv5')
except:
    pass
try:
    list_obs.remove('GPCPv2.3')
except:
    pass
print(list_obs)

#
# finding file and variable name in file for each observations dataset
#
dict_obs = dict()

for obs in list_obs:
    # be sure to add your datasets to EnsoCollectionsLib.ReferenceObservations if needed
    dict_var = ReferenceObservations(obs)['variable_name_in_file']
    dict_obs[obs] = dict()
    for var in list_variables:
        #
        # finding variable name in file
        #
        try: var_in_file = dict_var[var]['var_name']
        except:
            print('\033[95m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m')
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            # finding file for 'obs', 'var'
            file_name = param.reference_data_path[obs].replace('VAR',var0)
            file_areacell = None ## temporary for now
            try:
                file_landmask = param.reference_data_lf_path[obs]
            except:
                file_landmask = None

            try:
                areacell_in_file = dict_var['areacell']['var_name']
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var['landmask']['var_name']
            except:
                landmask_in_file = None
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                list_files = [param.reference_data_path[obs].replace('VAR',var1) for var1 in var_in_file]
                list_areacell = [file_areacell for var1 in var_in_file]
                list_name_area = [areacell_in_file for var1 in var_in_file]
                list_landmask = [param.reference_data_lf_path[obs] for var1 in var_in_file]
                list_name_land = [landmask_in_file for var1 in var_in_file]
            else:
                list_files = file_name
                list_areacell = file_areacell
                list_name_area = areacell_in_file
                list_landmask = file_landmask
                list_name_land = landmask_in_file
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file,
                                  'path + filename_area': list_areacell, 'areaname': list_name_area,
                                  'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}

print('PMPdriver: dict_obs readin end')

#=================================================
# Loop for Models 
#-------------------------------------------------
# Dictionary to save result 
def tree(): return defaultdict(tree)
enso_stat_dic = tree() # Use tree dictionary to avoid declearing everytime

# finding file and variable name in file for each observations dataset
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()['variable_name_in_file']

for mod in models:
    try:
        dict_mod = {mod: {}}
        print('PMPdriver: var loop start for model ', mod)
        for var in list_variables:
            # finding variable name in file
            var_in_file = dict_var[var]['var_name']
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'mod', 'var'
            #
            if mip == 'cmip5': realization = 'r1i1p1'
            elif mip == 'cmip6': realization = 'r1i1p1f1'
            # -- TEMPORARY --
            if mip =='cmip6' and mod in ['CESM2-WACCM']: realization = 'r2i1p1f1'
            # -- TEMPORARY END --
            file_name = modpath(mip=mip, model=mod, realization=realization, variable=var0)
            file_areacell = None ## temporary for now
            file_landmask = modpath_lf(mip=mip, model=mod)
            # -- TEMPORARY --
            if mip == 'cmip6' and mod in ['IPSL-CM6A-LR', 'CNRM-CM6-1']:
                file_landmask = '/work/lee1043/ESGF/CMIP6/CMIP/'+mod+'/sftlf_fx_'+mod+'_historical_r1i1p1f1_gr.nc'
            # -- TEMPORARY END --
            try:
                areacell_in_file = dict_var['areacell']['var_name']
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var['landmask']['var_name']
            except:
                landmask_in_file = None
    
            if isinstance(var_in_file, list):
                list_files = list()
                list_files = [modpath(model=mod, variable=var1) for var1 in var_in_file]
                list_areacell = [file_areacell for var1 in var_in_file]
                list_name_area = [areacell_in_file for var1 in var_in_file]
                list_landmask = [modpath_lf(mip=mip, model=mod) for var1 in var_in_file]
                list_name_land = [landmask_in_file for var1 in var_in_file]
            else:
                list_files = file_name
                list_areacell = file_areacell
                list_name_area = areacell_in_file
                list_landmask = file_landmask
                list_name_land = landmask_in_file
            dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file,
                                  'path + filename_area': list_areacell, 'areaname': list_name_area,
                                  'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}
    
        print('PMPdriver: var loop end')
    
        # dictionary needed by EnsoMetrics.ComputeMetricsLib.ComputeCollection
        dictDatasets = {'model': dict_mod, 'observations': dict_obs}
        # regridding dictionary (only if you want to specify the regridding)
        dict_regrid = {}
        # dict_regrid = {
        #     'regridding': {
        #         'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
        #         'newgrid_name': 'generic 1x1deg'},
        # }
        # Computes the metric collection
        netcdf_name = StringConstructor(param.netcdf_name)(model=mod) 
        netcdf = os.path.join(netcdf_path, netcdf_name)
        dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, netcdf=param.nc_out,
                                             netcdf_name=netcdf, debug=param.debug)
    except Exception as e: 
        print('failed for ', mod)
        print(e)
        if not param.debug:
            pass
  
print('PMPdriver: model loop end')

#=================================================
#  OUTPUT METRICS TO JSON FILE
#-------------------------------------------------
# disclaimer and reference for JSON header
disclaimer = open(
    os.path.join(
        sys.prefix,
        "share",
        "pmp",
        "disclaimer.txt")).read()
if param.metricsCollection == 'MC1':
    reference = "The statistics in this file are based on Bellenger, H et al. Clim Dyn (2014) 42:1999-2018. doi:10.1007/s00382-013-1783-z"
elif param.metricsCollection == 'ENSO_perf':
    reference = "GFDL..."
elif param.metricsCollection == 'ENSO_tel':
    reference = "MC3 for ENSO Teleconnection..."

# First JSON for metrics results
enso_stat_dic['obs'] = dict_obs
enso_stat_dic['model'] = dict_metric
metrics_dictionary = collections.OrderedDict()
metrics_dictionary["DISCLAIMER"] = disclaimer
metrics_dictionary["REFERENCE"] = reference
metrics_dictionary["RESULTS"] = enso_stat_dic

OUT = pcmdi_metrics.io.base.Base(outdir(output_type='metrics_results'), param.json_name+'.json')
OUT.write(
    metrics_dictionary,
    json_structure=["type", "data", "metric", "item", "value or description"],
    indent=4,
    separators=(
        ',',
        ': '),
    sort_keys=True)

# Second JSON for dive down information
diveDown_dictionary = collections.OrderedDict()
diveDown_dictionary["DISCLAIMER"] = disclaimer
diveDown_dictionary["REFERENCE"] = reference
diveDown_dictionary["RESULTS"] = {}
diveDown_dictionary["RESULTS"]["model"] = dict_dive

OUT2 = pcmdi_metrics.io.base.Base(outdir(output_type='metrics_results'), param.json_name+'_diveDown.json')
OUT2.write(
    dict_dive,
    json_structure=["type", "data", "metric", "item", "value or description"],
    indent=4,
    separators=(
        ',',
        ': '),
    sort_keys=True)

sys.exit('done')
