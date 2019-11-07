#!/usr/bin/env python
# =================================================
# Dependencies
# -------------------------------------------------
from __future__ import print_function

import collections
import copy
import glob
import json
import os
import pcmdi_metrics
import pkg_resources
import sys

from collections import defaultdict
from genutil import StringConstructor
from pcmdi_metrics.driver.pmp_parser import PMPParser
from PMPdriver_lib import ReadOptions
from PMPdriver_lib import metrics_to_json
from PMPdriver_lib import sort_human
from PMPdriver_lib import tree
from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection

# =================================================
# Collect user defined options
# -------------------------------------------------
param = ReadOptions()

# Pre-defined options
mip = param.mip
exp = param.exp
print('mip:', mip)
print('exp:', exp)

# Path to model data as string template
modpath = param.process_templated_argument("modpath")
modpath_lf = param.process_templated_argument("modpath_lf")

# Check given model option
models = param.modnames

# Realizations
realization = param.realization
print('realization: ', realization)

# Include all models if conditioned
if ('all' in [m.lower() for m in models]) or (models == 'all'):
    models = [p.split('.')[1]
        for p in glob.glob(modpath(
            mip=mip,
            exp=exp,
            model='*',
            realization=realization,
            variable='ts'))]
    # remove duplicates
    models = sorted(list(dict.fromkeys(models)), key=lambda s: s.lower())

print('models:', models)


# Metrics Collection
mc_name = param.metricsCollection 
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())
print('mc_name:', mc_name)

# Output
outdir_template = param.process_templated_argument("results_dir")
outdir = StringConstructor(str(outdir_template(
    output_type='%(output_type)',
    mip=mip, exp=exp, metricsCollection=mc_name)))
netcdf_path = outdir(output_type='diagnostic_results')
json_name_template = param.process_templated_argument("json_name")
netcdf_name_template = param.process_templated_argument("netcdf_name")

print('outdir:', str(outdir_template(
    output_type='%(output_type)',
    mip=mip, exp=exp, metricsCollection=mc_name))) 
print('netcdf_path:', netcdf_path)

# Switches
debug = param.debug
print('debug:', debug)

# =================================================
# Prepare loop iteration
# -------------------------------------------------
# Environmental setup
try:
    egg_pth = pkg_resources.resource_filename(
        pkg_resources.Requirement.parse("pcmdi_metrics"), "share/pmp")
except Exception:
    # python 2 seems to fail when ran in home directory of source?
    #egg_pth = os.path.join(os.getcwd(), "share", "pmp")
    egg_pth = os.path.join(sys.prefix, "share", "pmp")
print('egg_pth:', egg_pth)

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

            try:
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
                    try:
                        list_landmask = [param.reference_data_lf_path[obs] for var1 in var_in_file]
                    except:
                        list_landmask = None
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
            except:
                print('\033[95m' + 'Observation dataset' + str(obs) + " is not given for variable " + str(var) + '\033[0m')

print('PMPdriver: dict_obs readin end')

# =================================================
# Loop for Models 
# -------------------------------------------------
# finding file and variable name in file for each observations dataset
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()['variable_name_in_file']

print('models:', models)

for mod in models:
    print(' ----- model: ', mod, ' ---------------------')
    print('PMPdriver: var loop start for model ', mod)
    dict_mod = {mod: {}}
    dict_metric[mod], dict_dive[mod] = dict(), dict()

    model_path_list = os.popen(
        'ls '+modpath(mip=mip, exp=exp, model=mod, realization=realization,
        variable='ts')).readlines()

    model_path_list = sort_human(model_path_list)
    if debug:
        print('model_path_list:', model_path_list)

    # Find where run can be gripped from given filename template for modpath
    run_in_modpath = modpath(mip=mip, exp=exp, model=mod, realization=realization,
        variable=var).split('/')[-1].split('.').index(realization)
    # Collect available runs
    runs_list = [model_path.split('/')[-1].split('.')[run_in_modpath] for model_path in model_path_list]
    if debug:
        print('runs_list:', runs_list)

    # =================================================
    # Loop for Realizations
    # -------------------------------------------------
    for run in runs_list:

        print(' --- run: ', run, ' ---')
        dict_mod = {mod: {}}
        #dict_mod[mod][run] = {}

        if debug:
            print('list_variables:', list_variables)
    
        try:
            for var in list_variables:
                print(' --- var: ', var, ' ---')
                # finding variable name in file
                var_in_file = dict_var[var]['var_name']
                if isinstance(var_in_file, list):
                    var0 = var_in_file[0]
                else:
                    var0 = var_in_file
                #
                # finding file for 'mod', 'var'
                #
                file_name = modpath(mip=mip, exp=exp, model=mod, realization=run, variable=var0)
                file_areacell = None ## temporary for now
                file_landmask = modpath_lf(mip=mip, model=mod)
                # -- TEMPORARY --
                if mip == 'cmip6' and mod in ['IPSL-CM6A-LR', 'CNRM-CM6-1']:
                    try:
                        file_landmask = '/work/lee1043/ESGF/CMIP6/CMIP/'+mod+'/sftlf_fx_'+mod+'_historical_r1i1p1f1_gr.nc'
                    except:
                        pass
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
                    list_files = [modpath(mip=mip, exp=exp, model=mod, realization=realization, variable=var1) for var1 in var_in_file]
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

                dict_mod[mod][var] = {
                    'path + filename': list_files, 'varname': var_in_file,
                    'path + filename_area': list_areacell, 'areaname': list_name_area,
                    'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}
        
                print('PMPdriver: var loop end')
            
            # dictionary needed by EnsoMetrics.ComputeMetricsLib.ComputeCollection
            dictDatasets = {'model': dict_mod, 'observations': dict_obs}
            print('dictDatasets:')
            print(json.dumps(dictDatasets, indent=4, sort_keys=True))

            # regridding dictionary (only if you want to specify the regridding)
            dict_regrid = {}
            """
            # Usage of dict_regrid (select option as below):
            dict_regrid = {
                'regridding': {
                    'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                    'newgrid_name': 'generic 1x1deg'},
            }
            """

            # Prepare netcdf file setup
            json_name = json_name_template(mip=mip, exp=exp, metricsCollection=mc_name, model=mod, realization=run)
            netcdf_name = netcdf_name_template(mip=mip, exp=exp, metricsCollection=mc_name, model=mod, realization=run)
            netcdf = os.path.join(netcdf_path, netcdf_name)

            if debug:
                print('file_name:', file_name)
                print('list_files:', list_files)
                print('netcdf_name:', netcdf_name)
                print('json_name:', json_name)

            # Computes the metric collection
            dict_metric[mod][run], dict_dive[mod][run] = ComputeCollection(mc_name, dictDatasets, mod, netcdf=param.nc_out,
                                                     netcdf_name=netcdf, debug=debug)
            if debug:
                print('file_name:', file_name)
                print('list_files:', list_files)
                print('netcdf_name:', netcdf_name)
                print('dict_metric:')
                print(json.dumps(dict_metric, indent=4, sort_keys=True))

            # OUTPUT METRICS TO JSON FILE
            metrics_to_json(mc_name, dict_obs, dict_metric, dict_dive, egg_pth, outdir, json_name, mod=mod, run=run)

        except Exception as e: 
            print('failed for ', mod, run)
            print(e)
            if not debug:
                pass
      
print('PMPdriver: model loop end')

# =================================================
# OUTPUT METRICS TO JSON FILE
# -------------------------------------------------
json_name = json_name_template(mip=mip, exp=exp, metricsCollection=mc_name, model='all', realization='all')
metrics_to_json(mc_name, dict_obs, dict_metric, dict_dive, egg_pth, outdir, json_name)
