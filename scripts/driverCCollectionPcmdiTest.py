from cdms2 import open as CDMS2open
from os.path import join as join_path
from os import environ
from sys import exit

from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection

# metric collection
mc_name = 'MC2'
##mc_name = 'MC1'
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())
print 'jwlee_debug: mc_name:', mc_name
print 'jwlee_debug: list_metric:', list_metric

# parameters
project = 'CMIP5'
experiment = 'historical'
ensemble = 'r1i1p1'
frequency = 'mon'
realm = 'O'

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print 'list_variables:', list_variables

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
print 'list_obs:', list_obs

################################################
# Below is something should go back to parameter file
obspath = {
    'ERA-Interim': '/work/lee1043/DATA/reanalysis/ERAINT/mon/ERA-Interim_VAR_mo.xml',
    'HadISST': '/clim_obs/obs/ocn/mo/tos/UKMETOFFICE-HadISST-v1-1/130122_HadISST_sst.nc',
    'OISST': '/work/lee1043/DATA/OISST/xmls/OISST_tos_mo.xml',
    'Tropflux': '/work/lee1043/DATA/TropFlux/monthly/xmls/Tropflux_VAR_mo.xml',
    'OAFlux': '/work/lee1043/DATA/OAFlux/xmls/OAFlux_VAR_mo.xml',
}

modpath = '/work/lee1043/ESGF/xmls/cmip5/historical/mo/VAR/cmip5.MOD.historical.r1i1p1.mo.VAR.xml'
################################################

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
            file_name = obspath[obs].replace('VAR',var0)
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                for var1 in var_in_file:
                    file_name = obspath[obs].replace('VAR',var1)
                    list_files.append(file_name)
            else:
                list_files = file_name
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file}

# models
list_models = ['IPSL-CM5B-LR']
print 'jwlee_debug: model loop start for ', list_models
#
# finding file and variable name in file for each observations dataset
#
dict_metric = dict()
dict_var = CmipVariables()['variable_name_in_file']
for mod in list_models:
    print 'jwlee_debug: mod:', mod
    dict_mod = {mod: {}}
    # ------------------------------------------------
    # @jiwoo: between these dash the program is a bit ad hoc...
    # it works well for me because I am looking for sst and taux on the ocean grid, and fluxes [lhf, lwr, swr, shf, thf]
    # on the atmosphere grid
    # if you want to use atmosphere only, do not use this or create your own way to find the equivalent between the
    # variable name in the program and the variable name in the file
    print 'jwlee_debug: var loop start'
    for var in list_variables:
        ###print 'jwlee_debug: var:', var
        #
        # finding variable name in file
        #
        var_in_file = dict_var[var]['var_name']
        ###print 'jwlee_debug: var_in_file:', var_in_file
        if isinstance(var_in_file, list):
            var0 = var_in_file[0]
        else:
            var0 = var_in_file
        #
        # finding file for 'mod', 'var'
        #
        file_name = modpath.replace('MOD',mod).replace('VAR',var0)
        if isinstance(var_in_file, list):
            list_files = list()
            for var1 in var_in_file:
                file_name = modpath.replace('MOD',mod).replace('VAR',var1)
                list_files.append(file_name)
        else:
            list_files = file_name
        ###print 'jwlee_debug: list_files:', list_files
        # ------------------------------------------------
        dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file}
    print 'jwlee_debug: var loop end'
    # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
    # @jiwoo the ComputeCollection function it still on development and it does not read the observations requirement
    # defined in the metric collection, i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
    # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must include
    # every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its coming from
    # another dataset
    dictDatasets = {'model': dict_mod, 'observations': dict_obs}
    # Computes the metric collection
    print 'jwlee_debug: computes metric collection start'
    dict_metric[mod] = ComputeCollection(mc_name, dictDatasets)
    print 'jwlee_debug: computes metric collection end'
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
