from cdms2 import open as CDMS2open
from os.path import join as join_path
from os import environ
from sys import exit

#from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
#from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection

xmldir = environ['XMLDIR']


def find_xml_cmip(model, project, experiment, ensemble, frequency, variable):
    file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                          '_eq_pac_' + str(frequency) + '_regular_1x1_grid.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print('\033[95m' + str().ljust(5) + "CMIP var " + str(variable) + " cannot be found " + '\033[0m')
        print('\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m')
        print('\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m')
        exit(1)
    return file_name


def find_xml_obs(obs,frequency, variable):
    file_name = join_path(xmldir, 'obs_' + str(obs) + '_glob_' + str(frequency) + '_O.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print('\033[95m' + str().ljust(5) + "obs var " + str(variable) + " cannot be found" + '\033[0m')
        print('\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m')
        print('\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m')
        exit(1)
    return file_name


# metric collection
mc_name = 'ENSO_perf'#'MC1'#
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

# parameters
project = 'CMIP5'
experiment = 'piControl'#'historical'
ensemble = 'r1i1p1'
frequency = 'mon'
realm = 'A'

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print('\033[95m' + str(list_variables) + '\033[0m')

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in list(dict_var_obs.keys()):
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
if mc_name == 'MC1':
    list_obs = ['Tropflux']
elif mc_name == 'ENSO_perf':
    list_obs = ['Tropflux','GPCPv2.3']
elif mc_name == 'ENSO_tel':
    list_obs = ['HadISST','GPCPv2.3']
print('\033[95m' + str(list_obs) + '\033[0m')


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
            print('\033[95m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m')
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'obs', 'var'
            #
            # @jiwoo: pretty easy as I have all variables in one file
            file_name = find_xml_obs(obs, frequency, var0)
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                for var1 in var_in_file:
                    list_files.append(file_name)
            else:
                list_files = file_name
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file}

# models
list_models = ['ACCESS1-0', 'ACCESS1-3', 'BNU-ESM', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CESM1-FASTCHEM', 'CESM1-WACCM',
               'CMCC-CESM', 'CMCC-CMS', 'CMCC-CM', 'CNRM-CM5-2', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2', 'FGOALS-s2',
               'FIO-ESM', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H-CC', 'GISS-E2-R-CC',
               'HadGEM2-CC', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC-ESM-CHEM',
               'MIROC-ESM', 'MIROC4h', 'MIROC5', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P', 'MRI-CGCM3', 'NorESM1-ME',
               'NorESM1-M', 'bcc-csm1-1-m', 'bcc-csm1-1', 'inmcm4']
# finding file and variable name in file for each observations dataset
#
dict_metric, dict_dive = dict(), dict()
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
        if var == 'pr':
            var0 = 'pr'
        elif var == 'sst':
            var0 = 'sst'
        elif var == 'taux':
            var0 = 'tauu'
        #
        # finding file for 'mod', 'var'
        #
        # @jiwoo: first try in the realm 'O' (for ocean)
        file_name = find_xml_cmip(mod, project, experiment, ensemble, frequency, var0)
        # ------------------------------------------------
        dict_mod[mod][var] = {'path + filename': file_name, 'varname': var0}
    # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
    # @jiwoo the ComputeCollection function it still on development and it does not read the observations requirement
    # defined in the metric collection, i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
    # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must include
    # every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its coming from
    # another dataset
    dictDatasets = {'model': dict_mod, 'observations': dict_obs}
    # Computes the metric collection
    dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, dive_down=True)
    # Prints the metrics values
    for ii in range(3): print('')
    print('\033[95m' + str().ljust(5) + str(mod) + '\033[0m')
    list_metric = list(dict_metric[mod]['value'].keys())
    for metric in list_metric:
        print('\033[95m' + str().ljust(10) + str(metric) + '\033[0m')
        metric_dict = dict_metric[mod]['value'][metric]['metric']
        for ref in list(metric_dict.keys()):
            print('\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + str(metric_dict[ref]['value']) \
                  + ', error = ' + str(metric_dict[ref]['value_error']) + '\033[0m')

# import plot_frame as PFRAME
# path_plot = '/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Plots'
# list_metric = sorted(dict_metric[list_models[0]]['value'].keys())
# for met in list_metric:
#     list_mod = dict_metric.keys()
#     list_mod.sort(key=lambda v: v.upper())
#     ref = dict_metric[list_models[0]]['value'][met]['metric'].keys()[0]
#     tab = [dict_metric[mod]['value'][met]['metric'][ref]['value'] for mod in list_mod]
#     x_dict = dict((elt,mod) if elt%3==0 else (elt,'') for elt,mod in enumerate(list_mod))
#     print str(max(tab))+' for model '+str(list_mod[tab.index(max(tab))])
#     name_png = path_plot+'/'+met+'_'+str(len(list_mod))+'models'
#     PFRAME.plot_enso_stat(tab, x_axis=[-1,len(list_mod)], x_dico=x_dict, name_in_xlabel=True, draw_all_ylines=True, name=met, colormap='bl_to_darkred', \
#                           path_plus_name_png=name_png, save_ps=False, draw_white_background=True, bg=1)

