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
            print str().ljust(5) + "CMIP var " + str(variable) + " cannot be found (realm A and O)"
            print str().ljust(10) + "file_name = " + str(file_name)
            print str().ljust(10) + "variables = " + str(listvar1)
            print str().ljust(10) + "AND"
            print str().ljust(10) + "variables = " + str(listvar2)
            exit(1)
    return file_name


def find_xml_obs(obs,frequency, variable):
    file_name = join_path(xmldir, 'obs_' + str(obs) + '_glob_' + str(frequency) + '_O.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print str().ljust(5) + "obs var " + str(variable) + " cannot be found"
        print str().ljust(10) + "file_name = " + str(file_name)
        print str().ljust(10) + "variables = " + str(listvar1)
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

# @jewoo: I am lazy so I am using only one obvervations dataset
#list_obs = ['Tropflux']
list_obs = ['IPSL-CM5B-LR']

#
# finding file and variable name in file for each observations dataset
#
dict_obs = dict()
for obs in list_obs:
    # @jewoo: be sure to add your datasets to EnsoCollectionsLib.ReferenceObservations if needed
#    dict_var = ReferenceObservations(obs)['variable_name_in_file']
    dict_var = CmipVariables()['variable_name_in_file']
    dict_obs[obs] = dict()
    for var in list_variables:
        #
        # finding variable name in file
        #
        # @jewoo: correct / adapt the 'varname' in
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
            # @jewoo: pretty easy as I have all variables in one file
#            file_name = find_xml_obs(obs, frequency, variable)
            file_name = find_xml_cmip(obs, project, experiment, ensemble, frequency, realm, var0)
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = list()
                for var1 in var_in_file:
                    list_files.append(file_name)
            else:
                list_files = file_name
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file}

# models
list_models = ['CNRM-CM5']
#
# finding file and variable name in file for each observations dataset
#
dict_metric = dict()
dict_var = CmipVariables()['variable_name_in_file']
for mod in list_models:
    dict_mod = {mod: {}}
    # ------------------------------------------------
    # @jewoo: between these dash the program is a bit ad hoc...
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
        # @jewoo: first try in the realm 'O' (for ocean)
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
    # @jewoo the ComputeCollection function it still on development and it does not read the observations requirement
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
