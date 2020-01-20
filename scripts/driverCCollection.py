from cdms2 import open as CDMS2open
from os.path import join as join_path
from os import environ
from sys import exit

#from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
#from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection

xmldir = environ['XMLDIR']


def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                          '_glob_' + str(frequency) + '_' + str(realm) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
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
            print('\033[95m' + str().ljust(5) + "CMIP var " + str(variable) +\
                  " cannot be found (realm A and O)" + '\033[0m')
            print('\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m')
            print('\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m')
            print('\033[95m' + str().ljust(10) + "AND" + '\033[0m')
            print('\033[95m' + str().ljust(10) + "variables = " + str(listvar2) + '\033[0m')
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
mc_name = 'ENSO_tel'#'MC1'#'ENSO_perf'
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

# parameters
project = 'CMIP5'
experiment = 'historical'
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
list_obs = ['HadISST']#['Tropflux','HadISST']

if mc_name == 'ENSO_tel':
    list_obs = ['GPCPv2.3']
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
list_models = ['CNRM-CM5']#['IPSL-CM5B-LR']
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
        file_name = find_xml_cmip(mod, project, experiment, ensemble, frequency, realm, var0)
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
    # regridding dictionary
    dict_regrid = {
        'regridding': {
            'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
            'newgrid_name': 'generic 1x1deg'},
    }
    # Computes the metric collection
    dict_metric[mod] = ComputeCollection(mc_name, dictDatasets, user_regridding=dict_regrid)
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

