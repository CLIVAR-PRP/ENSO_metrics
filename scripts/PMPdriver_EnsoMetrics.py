#!/usr/bin/env python
#=================================================
# Dependencies
#-------------------------------------------------
from cdms2 import open as CDMS2open

import collections
import copy
import json
import os
import pcmdi_metrics
import sys

from collections import defaultdict
from pcmdi_metrics.pcmdi.pmp_parser import PMPParser

from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
#from EnsoMetrics.EnsoMetricsLib import EnsoAlphaLhf, EnsoAlphaLwr, EnsoAlphaSwr, EnsoAlphaThf, EnsoAmpl, EnsoMu, EnsoRMSE, EnsoSeasonality

from pmpParser import ReadOptions

#=================================================
# Collect user defined options
#-------------------------------------------------
param = ReadOptions()

#=================================================
# User input option
#-------------------------------------------------
mc_name = param.metricsCollection 
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

#=================================================
# Prepare loop iteration
#-------------------------------------------------
# setup an output directory 
try:
    os.mkdir(param.outpathjson)
except BaseException:
    pass

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print list_variables

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
print list_obs

##### TEST
#list_obs = ['ERA-Interim', 'HadISST', 'Tropflux']

#
# finding file and variable name in file for each observations dataset
#
dict_obs = dict()
for obs in list_obs:
    # @jiwoo: be sure to add your datasets to EnsoCollectionsLib.ReferenceObservations if needed
    dict_var = ReferenceObservations(obs)['variable_name_in_file']
    dict_obs[obs] = dict()
    for var in list_variables:
        #
        # finding variable name in file
        #
        # @jiwoo: correct / adapt the 'varname' in
        # EnsoCollectionsLib.ReferenceObservations(obs)['variable_name_in_file'][var] if it is not correct or if you
        # changed a name in the xml
        # I usually alias the variable names from observations and models in the xml in order to have the same name
        # for sst (or any other variable) in every xml. This way I don not need to go through this function to know the
        # variable name in file
        try: var_in_file = dict_var[var]['var_name']
        except:
            print var + " in not available for " + str(obs) + " or unscripted"
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'obs', 'var'
            #
            # @jiwoo: pretty easy as I have all variables in one file
            #file_name = find_xml_obs(obs, frequency, var0)
            #file_name = obs+'_'+var0
            file_name = param.obspath[obs].replace('VAR',var0)
            if not os.path.isfile(file_name): print 'File not available: ', file_name
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                for var1 in var_in_file:
                    #list_files.append(file_name)
                    file_name1 = param.obspath[obs].replace('VAR',var1)
                    if not os.path.isfile(file_name1): print 'File not available: ', file_name1
                    list_files.append(file_name1)
            else:
                list_files = file_name
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file}

print dict_obs

#=================================================
# Loop for Observation and Models 
#-------------------------------------------------
# Insert observation at the beginning of the loop ---
list_models = copy.copy(param.modnames)
#list_models.insert(0,'obs')
print list_models

# Dictionary to save result ---
def tree(): return defaultdict(tree)
enso_stat_dic = tree() # Use tree dictionary to avoid declearing everytime

#######################################################
#######################################################
#
# finding file and variable name in file for each observations dataset
#
dict_metric = dict()
dict_var = CmipVariables()['variable_name_in_file']
for mod in list_models:
    dict_mod = {mod: {}}
    # ------------------------------------------------
    # @jiwoo: between these dash the program is a bit ad hoc...
    # it works well for me because I am looking for sst and taux on the ocean grid, and fluxes [lhf, lwr, swr, shf, thf]
    # on the atmosphere grid
    # if you want to use atmosphere only, do not use this or create your own way to find the equivalent between the
    # variable name in the program and the variable name in the file
    for var in list_variables:
        #
        # finding variable name in file
        #
        var_in_file = dict_var[var]['var_name']
        if isinstance(var_in_file, list):
            var0 = var_in_file[0]
        else:
            var0 = var_in_file
        #
        # finding file for 'mod', 'var'
        #
        # @jiwoo: first try in the realm 'O' (for ocean)
        #file_name = find_xml_cmip(mod, project, experiment, ensemble, frequency, realm, var0)
        file_name = param.modpath.replace('MOD',mod).replace('VAR',var0)
        print mod, var0, file_name
        file = CDMS2open(file_name)
        if isinstance(var_in_file, list):
            list_files = list()
            for var1 in var_in_file:
                list_files.append(file_name)
        else:
            list_files = file_name
        # ------------------------------------------------
        dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file}
    # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
    # @jiwoo the ComputeCollection function it still on development and it does not read the observations requirement
    # defined in the metric collection, i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
    # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must include
    # every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its coming from
    # another dataset
    dictDatasets = {'model': dict_mod, 'observations': dict_obs}
    # Computes the metric collection
    dict_metric[mod] = ComputeCollection(mc_name, dictDatasets)
    # Prints the metrics values
    for ii in range (3): print ''
    print str().ljust(5) + str(mod)
    list_metric = dict_metric[mod]['metrics'].keys()
    for metric in list_metric:
        print str().ljust(10) + str(metric)
        metric_dict = dict_metric[mod]['metrics'][metric]['metric_values']
        for ref in metric_dict.keys():
            print str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + str(metric_dict[ref]['value']) + ', error = '\
                  + str(metric_dict[ref]['value_error'])
        raw_dict = dict_metric[mod]['metrics'][metric]['raw_values']['observations']
        for ref in raw_dict.keys():
            print str().ljust(15) + 'raw (diag) obs: ' + str(ref) + ' value = ' + str(raw_dict[ref]['value']) +\
                  ', error = ' + str(raw_dict[ref]['value_error'])
        raw_dict = dict_metric[mod]['metrics'][metric]['raw_values']['model']
        print str().ljust(15) + 'raw (diag) model: ' + str(mod) + ' value = ' + str(raw_dict['value']) +\
                  ', error = ' + str(raw_dict['value_error'])

sys.exit('test')
#######################################################
#######################################################


for mod in models:
    print ' ----- ', mod,' ---------------------'
  
    if mod == 'obs':
        sstName = param.sstNameObs
        tauxName = param.tauxNameObs

        sstFile = param.sstObsPath 
        tauxFile = param.tauxObsPath
    else:
        sstName = param.sstName
        tauxName = param.tauxName

        sstFile = (param.modpath.replace('MOD', mod)).replace('VAR',sstName) ## Will need land mask out at some point...!
        tauxFile = (param.modpath.replace('MOD', mod)).replace('VAR',tauxName)

    print sstFile
    print tauxFile
  
    try:
        for metric in list_metric:

            print metric

            if metric == 'EnsoAmpl':
                nBox = dict_mc['metrics_list'][metric]['regions']['sst']
                tmp_dict = EnsoAmpl(sstFile, sstName, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoSeasonality':
                nBox = dict_mc['metrics_list'][metric]['regions']['sst']
                tmp_dict = EnsoSeasonality(sstFile, sstName, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoMu':
                sstBox = dict_mc['metrics_list'][metric]['regions']['sst']
                tauxBox = dict_mc['metrics_list'][metric]['regions']['taux']
                tmp_dict = EnsoMu(sstFile, tauxFile, sstName, tauxName, sstBox, tauxBox)
                tmp_dict['input_data'] = [sstFile, tauxFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoRMSE' and mod != 'obs':
                nBox = dict_mc['metrics_list'][metric]['regions']['sst']
                tmp_dict = EnsoRMSE(sstFile, sstName, param.sstObsPath, param.sstNameObs, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

    except Exception as e: 
        print 'failed for ', mod
        print(e)
  
#=================================================
#  OUTPUT METRICS TO JSON FILE
#-------------------------------------------------
OUT = pcmdi_metrics.io.base.Base(os.path.abspath(param.outpathjson), param.outnamejson)

disclaimer = open(
    os.path.join(
        sys.prefix,
        "share",
        "pmp",
        "disclaimer.txt")).read()

metrics_dictionary = collections.OrderedDict()
metrics_dictionary["DISCLAIMER"] = disclaimer
metrics_dictionary["REFERENCE"] = "The statistics in this file are based on Bellenger, H et al. Clim Dyn (2014) 42:1999-2018. doi:10.1007/s00382-013-1783-z"
metrics_dictionary["RESULTS"] = enso_stat_dic

OUT.write(
    metrics_dictionary,
    json_structure=["model", "metric", "item", "value or description"],
    indent=4,
    separators=(
        ',',
        ': '),
    sort_keys=True)

sys.exit('done')
