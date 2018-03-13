# -*- coding:UTF-8 -*-
from os.path import join as join_path
from os import environ
from sys import exit

# uvcdat based functions:
from cdms2 import open as CDMS2open
from cdutil import area_weights as CDUTILarea_weights
from MV2 import masked_where as MV2masked_where
from MV2 import where as MV2where

# ENSO_metrics package functions:
#from lib.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
#from lib.EnsoComputeMetricsLib import ComputeMetric
from DriverPreprocessingUvcdatToolsLib import CheckTime, CommonPeriod, TimeBounds
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeMetric
import DriverPreprocessing

xmldir = environ['XMLDIR']


def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                          '_glob_' + str(frequency) + '_' + str(realm) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    new_realm = deepcopy(realm)
    if variable not in listvar1:
        if realm == 'O':
            new_realm = 'A'
        elif realm == 'A':
            new_realm = 'O'
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere)
        file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                              '_glob_' + str(frequency) + '_' + str(new_realm) + '.xml')
        xml = CDMS2open(file_name)
        listvar2 = sorted(xml.listvariables())
        if variable not in listvar2:
            print '\033[95m' + "CMIP var " + str(variable) + " cannot be found (realm A and O)" + '\033[0m'
            print '\033[95m' + str().ljust(5) + "file_name = " + str(file_name) + '\033[0m'
            print '\033[95m' + str().ljust(5) + "variables = " + str(listvar1) + '\033[0m'
            print '\033[95m' + str().ljust(5) + "AND" + '\033[0m'
            print '\033[95m' + str().ljust(5) + "variables = " + str(listvar2) + '\033[0m'
            exit(1)
    if new_realm == 'A':
        file_name_area = join_path(xmldir, str(model) + '_areacella_fx_' + str(experiment) + '_r0i0p0.nc')
        var_area = 'areacella'
    elif realm == 'O':
        file_name_area = join_path(xmldir, str(model) + '_areacello_fx_' + str(experiment) + '_r0i0p0.nc')
        var_area = 'areacello'
    return file_name, file_name_area, var_area


def find_xml_obs(obs,frequency, variable):
    file_name = join_path(xmldir, 'obs_' + str(obs) + '_glob_' + str(frequency) + '_O.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print '\033[95m' + "obs var " + str(variable) + " cannot be found" + '\033[0m'
        print '\033[95m' + str().ljust(5) + "file_name = " + str(file_name) + '\033[0m'
        print '\033[95m' + str().ljust(5) + "variables = " + str(listvar1) + '\033[0m'
        exit(1)
    return file_name


# metric collection
mc_name = 'MC1'
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

# parameters
project = 'CMIP5'
experiment = 'historical'
ensemble = 'r1i1p1'
frequency = dict_mc['common_collection_parameters']['frequency']
if frequency == 'daily':
    freq = 'day'
elif frequency == 'monthly':
    freq = 'mon'
else:
    freq = None
realm = 'O'

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print '\033[95m' + str(list_variables) + '\033[0m'

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
print '\033[95m' + str(list_obs) + '\033[0m'


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
            print '\033[95m' + str(var) + " in not available for " + str(obs) + " or unscripted" + '\033[0m'
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'obs', 'var'
            #
            # @jiwoo: pretty easy as I have all variables in one file
            file_name = find_xml_obs(obs, freq, var0)
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                for var1 in var_in_file:
                    list_files.append(file_name)
            else:
                list_files = file_name
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file}


# models
list_models = ['IPSL-CM5B-LR']
#
# finding file and variable name in file for each observations dataset
#
dict_model_metric = dict()
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
        #
        # finding file for 'mod', 'var'
        #
        # @jiwoo: first try in the realm 'O' (for ocean)
        if isinstance(var_in_file, list):
            list_files, list_files_area, list_var_area = list(), list(), list()
            for var1 in var_in_file:
                file_var, file_area, var_area = find_xml_cmip(mod, project, experiment, ensemble, freq, realm, var1)
                list_files.append(file_var)
                list_files_area.append(file_area)
                list_var_area.append(var_area)
        else:
            list_files, list_files_area, list_var_area = find_xml_cmip(mod, project, experiment, ensemble, freq, realm,
                                                                       var_in_file)
        # ------------------------------------------------
        # @jiwoo:
        # I propose to look for the areacell file here and join it to the var file. This way the variable and the
        # areacell can be in different files
        # if you don't have the corresponding areacell or if we do not need them (for example, for atmospheric
        # variables), you can set 'list_files_area' and 'list_var_area' to None
        dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file,
                              'path + filename_area': list_files_area, 'varname_area': list_var_area}
    # ------------------------------------------------
    # @jiwoo:
    # ok, now we are playing
    # I propose to put the loop on metrics here, and then do the preprocessing
    #
    # As a first step you can use the preprocessing tools that I developed in 'EnsoUvcdatToolsLib.py'
    # I moved the useful tools in a new file 'DriverPreprocessingUvcdatToolsLib.py' in the 'scripts' directory
    #
    # Then you can develop your own tools or use other language / software (like cdo) for the preprocessing
    # We (me / PCMDI / ESMValTool) will just have to agree on what keyword means what preprocessing
    #
    # ------------------------------------------------
    # loop over metrics
    #
    # create dictionary to save results for each metric
    dict_metric = dict()
    for metric in list_metric:
        print '\033[95m' + str(mod) + ", metric = " + str(metric) + '\033[0m'
        # read the dictionary of parameters for this metric from 'EnsoCollectionsLib.py'
        dict_param_metric = dict_mc['metrics_list'][metric]
        # list of variables for this metric
        list_var_metric = dict_param_metric['variables']
        # dictionary of preprocessing for this metric
        preprocessing = dict_param_metric['preprocessing']
        # period for model and obs
        period_mod = preprocessing[0]['selection_period_and_region']['period']['model']
        period_obs = preprocessing[0]['selection_period_and_region']['period']['observations']
        # create dictionary to save the array to compute the metric and the preprocessing steps
        dict_data_to_compute_metric = dict()
        dict_preprocessing_steps = dict()
        # ------------------------------------------------
        # loop over variables
        for nvar, var in zip(range(1, len(listvar)+1), listvar):
            print '\033[95m' + str(mod) + ", metric = " + str(metric) + ", variable = " + str(var) + '\033[0m'
            # region for this variable
            region = preprocessing[0]['selection_period_and_region']['regions'][var]
            # wanted units for this variable
            units = preprocessing[0]['selection_period_and_region']['units'][var]
            # file_name (model)
            file_name = dict_mod[mod][var]['path + filename']
            # var_name (model)
            var_name = dict_mod[mod][var]['varname']
            # file_name_area (model)
            file_name_area = dict_mod[mod][var]['path + filename_area']
            # var_name_area (model)
            var_name_area = dict_mod[mod][var]['varname_area']
            # var_name (model)
            # test if a regridding is needed
            regridding = False
            for ii in preprocessing.keys():
                if preprocessing[ii].keys()[0] == 'regridding':
                    grid = preprocessing[ii]['regridding']['newgrid']
                    regridding = True
                    print '\033[95m' + "grid '" + str(grid) + "' is needed" + '\033[0m'
                    # @jiwoo: here begins the hard part...
                    # this 'grid' is for now the name of a dataset
                    # either you find the file for this dataset (using find_xml_obs in this program), read any
                    # variable in the file (using ReadAndSelectRegion defined in 'DriverPreprocessingUvcdatToolsLib.py')
                    # and retrieve the grid (grid = tab.getGrid()) and put it in the preprocessing dictionary instead of
                    # the string in 'newgrid' (preprocessing[ii]['regridding']['newgrid'] = grid)
                    #
                    # the other choice is to do nothing. In this case grid defined in
                    # preprocessing[ii]['regridding']['newgrid_name'] will be created and the model / obs data will be
                    # regridded toward it
                    #
                    # the problem (whatever your choice above) is that model and observations will have a different mask
                    # you will have to create a common mask and apply it later
                    break
            # compute preprocessing (model)
            tab_mod, Nyears_mod, preprocessing_steps = DriverPreprocessing.preprocess(
                file_name, var_name, preprocessing, file_name_area=file_name_area, var_name_area=var_name_area,
                region=region, period=period_mod, frequency=frequency, units=units)
            print '\033[95m' + str(mod) + ", metric = " + str(metric) + ", variable = " + str(var) + ": done"\
                  + '\033[0m'
            dict_preprocessing_steps[var] = {'model': preprocessing_steps}
            # begin the common mask
            if regridding:
                mask = tab_mod.mask
            # ------------------------------------------------
            # loop over observations and do the same as for the model
            # list of observations for this metric and this variable
            list_observations = dict_param_metric['obs_name'][var]
            list_dict_obs = list()
            for nobs, obs in zip(range(len(list_observations)), list_observations):
                print '\033[95m' + str(obs) + ", metric = " + str(metric) + ", variable = " + str(var) + '\033[0m'
                # file_name (obs)
                file_name = dict_obs[obs][var]['path + filename']
                # var_name (obs)
                var_name = dict_obs[obs][var]['varname']
                # compute preprocessing (obs)
                tab_obs, Nyears_obs, preprocessing_steps = DriverPreprocessing.preprocess(
                    file_name, var_name, preprocessing, region=region, period=period_obs, frequency=frequency,
                    units=units)
                # continue the common mask
                if regridding:
                    mask = MV2where(tab_obs.mask, True, mask)
                list_dict_obs.append({'array': tab_obs, 'Nyears': Nyears_obs})
                print '\033[95m' + str(obs) + ", metric = " + str(metric) + ", variable = " + str(var) + ": done"\
                      + '\033[0m'
            dict_preprocessing_steps[var] = {'observations': preprocessing_steps}
            # apply the common mask on this variable
            if regridding:
                tab_mod = MV2masked_where(mask, tab_mod)
                for nobs in range(len(list_observations)):
                    list_dict_obs[nobs]['array'] = MV2masked_where(mask, list_dict_obs[nobs]['array'])
            # put variables in a dictionary (model & obs)
            dict_data_to_compute_metric['modelVar' + str(nvar)] = {'array': tab_mod, 'Nyears': Nyears_mod}
            dict_data_to_compute_metric['obsVar' + str(nvar)] = list_dict_obs
            dict_data_to_compute_metric['obsNameVar' + str(nvar)] = list_observations
            dict_data_to_compute_metric['regionVar' + str(nvar)] = region
        # ------------------------------------------------
        # Check if all variables cover the same period
        if len(listvar) > 0:
            tab = dict_data_to_compute_metric['modelVar1']['array']
            if tab.getTime() is not None:
                #
                # find the common period of the model variables
                #
                time_bounds = TimeBounds(dict_data_to_compute_metric['modelVar1']['array'])
                for nvar in range(2, len(listvar) + 1):
                    time_bounds = CommonPeriod(time_bounds, TimeBounds(
                        dict_data_to_compute_metric['modelVar' + str(nvar)]['array']))
                # get the model variables during this period
                for nvar in range(1, len(listvar) + 1):
                    tab_mod = CheckTime(dict_data_to_compute_metric['modelVar' + str(nvar)]['array'], time_bounds,
                                        frequency=frequency)
                    dict_data_to_compute_metric['modelVar' + str(nvar)]['array'] = tab_mod
                    if frequency == "daily":
                        Nyears = len(tab_mod) / 365
                    elif frequency == "yearly":
                        Nyears = len(tab_mod) / 1
                    else:
                        Nyears = len(tab_mod) / 12
                    dict_data_to_compute_metric['modelVar' + str(nvar)]['Nyears'] = Nyears
                #
                # Same for the observations
                #
                time_bounds = TimeBounds(dict_data_to_compute_metric['obsVar1'][0]['array'])
                for nvar in range(1, len(listvar) + 1):
                    for nobs in range(len(list_observations)):
                        time_bounds = CommonPeriod(time_bounds, TimeBounds(
                            dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['array']))
                for nvar in range(1, len(listvar) + 1):
                    for nobs in range(len(list_observations)):
                        tab_obs = CheckTime(dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['array'],
                                            time_bounds, frequency=frequency)
                        dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['array'] = tab_obs
                        if frequency == "daily":
                            Nyears = len(tab_obs) / 365
                        elif frequency == "yearly":
                            Nyears = len(tab_obs) / 1
                        else:
                            Nyears = len(tab_obs) / 12
                        dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['Nyears'] = Nyears
        # ------------------------------------------------
        # get axes and weights
        for nvar in range(1, len(listvar) + 1):
            tab_mod = dict_data_to_compute_metric['modelVar' + str(nvar)]['array']
            dict_data_to_compute_metric['modelVar' + str(nvar)]['axes'] = tab_mod.getAxisList()
            dict_data_to_compute_metric['modelVar' + str(nvar)]['weights'] = CDUTILarea_weights(tab_mod)
            for nobs in range(len(list_observations)):
                tab_obs = dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['array']
                dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['axes'] = tab_obs.getAxisList()
                dict_data_to_compute_metric['obsVar' + str(nvar)][nobs]['weights'] = CDUTILarea_weights(tab_obs)
        # ------------------------------------------------
        # compute the given metric
        dict_metric[metric] = ComputeMetric(mc_name, metric, mod, **dict_data_to_compute_metric)
        dict_metric[metric]['preprocessing_steps'] = dict_preprocessing_steps
    # ------------------------------------------------
    # save results for each model
    dict_model_metric[mod] = dict_metric
    # prints the metrics values
    for ii in range (3): print ''
    print '\033[95m' + str().ljust(5) + str(mod) + '\033[0m'
    for metric in list_metric:
        print '\033[95m' + str().ljust(10) + str(metric) + '\033[0m'
        dict_tmp1 = dict_model_metric[mod][metric]['metric_values']
        for ref in dict_tmp1.keys():
            print '\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + str(dict_tmp1[ref]['value'])\
                  + ', error = ' + str(dict_tmp1[ref]['value_error']) + '\033[0m'
        dict_tmp1 = dict_model_metric[mod][metric]['raw_values']['observations']
        for ref in dict_tmp1.keys():
            print '\033[95m' + str().ljust(15) + 'raw (diag) obs: ' + str(ref) + ' value = '\
                  + str(dict_tmp1[ref]['value']) + ', error = ' + str(dict_tmp1[ref]['value_error']) + '\033[0m'
        dict_tmp1 = dict_model_metric[mod][metric]['raw_values']['model']
        print '\033[95m' + str().ljust(15) + 'raw (diag) model: ' + str(mod) + ' value = ' + str(dict_tmp1['value'])\
              + ', error = ' + str(dict_tmp1['value_error']) + '\033[0m'
