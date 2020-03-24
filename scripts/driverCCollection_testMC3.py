# -*- coding:UTF-8 -*-
from __future__ import print_function
from getpass import getuser as GETPASSgetuser
from cdms2 import open as CDMS2open
import json
from os import environ as OSenviron
from os.path import join as OSpath__join
from sys import exit as SYSexit
from sys import path as SYSpath

# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# xmldir
xmldir = OSenviron['XMLDIR']
if user_name == "yplanton":
    SYSpath.insert(0, "/home/" + user_name + "/Obs2Obs/ENSO_metrics/lib")
    path_obs = "/data/" + user_name + "/Obs"
    path_netcdf = "/data/" + user_name + "/ENSO_metrics/v20200311"
else:
    path_obs = "/Users/" + user_name + "/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Netcdf/Obs/Global"
    path_netcdf = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Test"

#from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
#from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection, ComputeCollection_ObsOnly


def find_xml(name, frequency, variable, project='', experiment='', ensemble='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_name, file_area, file_land = find_xml_obs(name, frequency, variable)
    else:
        file_name, file_area, file_land = find_xml_cmip(name, project, experiment, ensemble, frequency, realm, variable)
    return file_name, file_area, file_land


def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    file_name = OSpath__join(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                             '_glob_' + str(frequency) + '_' + str(realm) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        if realm == 'O':
            new_realm = 'A'
        elif realm == 'A':
            new_realm = 'O'
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere)
        file_name = OSpath__join(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                                 '_glob_' + str(frequency) + '_' + str(new_realm) + '.xml')
        xml = CDMS2open(file_name)
        listvar2 = sorted(xml.listvariables())
        if variable not in listvar2:
            print('\033[95m' + str().ljust(5) + "CMIP var " + str(variable) + " cannot be found (realm A and O)"
                  + '\033[0m')
            print('\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m')
            print('\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m')
            print('\033[95m' + str().ljust(10) + "AND" + '\033[0m')
            print('\033[95m' + str().ljust(10) + "variables = " + str(listvar2) + '\033[0m')
            SYSexit("")
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=new_realm)
    else:
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=realm)
    return file_name, file_area, file_land


def find_xml_fx(name, project='', experiment='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        if user_name == "yplanton":
            file_area = None
            if name in ["20CRv2", "NCEP2"]:
                file_land = OSpath__join(path_obs, name + "/land.sfc.gauss.nc")
            elif name == "CMAP":
                file_land = OSpath__join(path_obs, name + "/lsmask_fx_cmap.nc")
            elif name == "GPCPv2.3":
                file_land = OSpath__join(path_obs, name + "/lsmask_fx_gpcpv2.3.nc")
            elif name == "OISSTv2":
                file_land = OSpath__join(path_obs, name + "/lsmask_fx_oisstv2.nc")
            else:
                file_land = None
        else:
            file_area = None
            if name in ["20CRv2", "CMAP", "ERA-Interim", "ERSSTv5", "GPCPv2.3", "NCEP2"]:
                file_land = OSpath__join(path_obs, "lsmask_" + name + ".nc")
            else:
                file_land = None
    else:
        file_area = OSpath__join(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_' +
                                 str(realm) + '_areacell.xml')
        file_land = OSpath__join(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_' +
                                 str(realm) + '_landmask.xml')
    try: xml = CDMS2open(file_area)
    except: file_area = None
    try: xml = CDMS2open(file_land)
    except: file_land = None
    return file_area, file_land


def find_xml_obs(obs, frequency, variable):
    # file_name = OSpath__join(xmldir, 'obs_' + str(obs) + '_glob_' + str(frequency) + '_O.xml')
    file_name = OSpath__join(xmldir, "obs_ENSO_metrics_" + str(obs) + ".xml")
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print('\033[95m' + str().ljust(5) + "obs var " + str(variable) + " cannot be found" + '\033[0m')
        print('\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m')
        print('\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m')
        SYSexit("")
    file_area, file_land = find_xml_fx(obs)
    return file_name, file_area, file_land


def save_json_obs(dict_in, json_name):
    dict_out = {"RESULTS": {"model": dict_in}}
    # save as json file
    with open(json_name + '.json', 'w') as outfile:
        json.dump(dict_out, outfile, sort_keys=True)
    return


# metric collection
mc_name = 'ENSO_tel'  # 'ENSO_perf'  # 'ENSO_proc'  #
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
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
if mc_name == 'MC1':
    list_obs = ['Tropflux']
elif mc_name == 'ENSO_perf':
    list_obs = ['20CRv2', 'AVISO', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'GPCPv2.3', 'HadISST', 'NCEP2',
                'SODA3.4.2', 'Tropflux']
elif mc_name == 'ENSO_tel':
    list_obs = ['20CRv2', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GPCPv2.3', 'HadISST', 'NCEP2', 'SODA3.4.2', 'Tropflux']
elif mc_name == 'ENSO_proc':
    list_obs = ['20CRv2', 'AVISO', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'HadISST', 'NCEP2', 'SODA3.4.2', 'Tropflux']
else:
    list_obs = ['20CRv2', 'AVISO', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'GPCPv2.3', 'HadISST', 'NCEP2',
                'SODA3.4.2', 'Tropflux']
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
        try:
            var_in_file = dict_var[var]['var_name']
        except:
            print('\033[95m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m')
        else:
            try:
                areacell_in_file = dict_var['areacell']['var_name']
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var['landmask']['var_name']
            except:
                landmask_in_file = None
            if isinstance(var_in_file, list):
                list_areacell, list_files, list_landmask, list_name_area, list_name_land = \
                    list(), list(), list(), list(), list()
                for var1 in var_in_file:
                    file_name, file_areacell, file_landmask = find_xml(obs, frequency, var1)
                    list_files.append(file_name)
                    list_areacell.append(file_areacell)
                    list_name_area.append(areacell_in_file)
                    list_landmask.append(file_landmask)
                    list_name_land.append(landmask_in_file)
            else:
                file_name, file_areacell, file_landmask = find_xml(obs, frequency, var_in_file)
                list_files = file_name
                list_areacell = file_areacell
                list_name_area = areacell_in_file
                list_landmask = file_landmask
                list_name_land = landmask_in_file
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file,
                                  'path + filename_area': list_areacell, 'areaname': list_name_area,
                                  'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}

# computing metrics on observations only
dictDatasets = {'observations': dict_obs}
netcdf_name = user_name + "_" + mc_name + "_OBSNAME"
netcdf = OSpath__join(path_netcdf, netcdf_name)
#dict_metric, dict_dive = ComputeCollection_ObsOnly(mc_name, dictDatasets, debug=False, netcdf=True, netcdf_name=netcdf)
#netcdf = netcdf.replace("OBSNAME", "observation")
#save_json_obs(dict_metric, netcdf)


# models
list_models = ['IPSL-CM5B-LR']#['CNRM-CM5']#['IPSL-CM5B-LR']#['CNRM-CM5','IPSL-CM5B-LR']#
ens = 'r1i1p1'
#
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
        try:
            areacell_in_file = dict_var['areacell']['var_name']
        except:
            areacell_in_file = None
        try:
            landmask_in_file = dict_var['landmask']['var_name']
        except:
            landmask_in_file = None
        if isinstance(var_in_file, list):
            list_areacell, list_files, list_landmask, list_name_area, list_name_land = \
                list(), list(), list(), list(), list()
            for var1 in var_in_file:
                file_name, file_areacell, file_landmask = \
                    find_xml(mod, frequency, var1, project=project, experiment=experiment, ensemble=ens,
                             realm=realm)
                list_files.append(file_name)
                list_areacell.append(file_areacell)
                list_name_area.append(areacell_in_file)
                list_landmask.append(file_landmask)
                list_name_land.append(landmask_in_file)
        else:
            file_name, file_areacell, file_landmask = \
                find_xml(mod, frequency, var_in_file, project=project, experiment=experiment, ensemble=ens,
                         realm=realm)
            list_files = file_name
            list_areacell = file_areacell
            list_name_area = areacell_in_file
            list_landmask = file_landmask
            list_name_land = landmask_in_file
        dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file,
                              'path + filename_area': list_areacell, 'areaname': list_name_area,
                              'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}
    # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
    # @jiwoo the ComputeCollection function it still on development and it does not read the observations requirement
    # defined in the metric collection, i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
    # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must include
    # every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its coming from
    # another dataset
    dictDatasets = {'model': dict_mod, 'observations': dict_obs}
    # regridding dictionary (only if you want to specify the regridding)
    dict_regrid = {}
    # dict_regrid = {
    #     'regridding': {
    #         'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
    #         'newgrid_name': 'generic 1x1deg'},
    # }
    # Computes the metric collection
    netcdf_name = user_name + "_" + mc_name + "_" + mod
    netcdf = OSpath__join(path_netcdf, netcdf_name)
    dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, mod, netcdf=True, netcdf_name=netcdf,
                                                         debug=True)
    tmp = sorted(dict_metric[mod]['value'].keys(), key=lambda v: v.upper())
    for kk in tmp:
        print(kk.ljust(13) + ': ' + str(dict_metric[mod]['value'][kk]['metric']))
    # Prints the metrics values
    for ii in range(3): print("")
    print('\033[95m' + str().ljust(5) + str(mod) + '\033[0m')
    list_metric = dict_metric[mod]['value'].keys()
    for metric in list_metric:
        print('\033[95m' + str().ljust(10) + str(metric) + '\033[0m')
        metric_dict = dict_metric[mod]['value'][metric]['metric']
        for ref in metric_dict.keys():
            print('\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + str(metric_dict[ref]['value'])
                  + ', error = ' + str(metric_dict[ref]['value_error']) + '\033[0m')
            if 'value2' in metric_dict[ref].keys():
                print('\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' +
                      str(metric_dict[ref]['value2']) + ', error = ' + str(metric_dict[ref]['value_error2']) +
                      '\033[0m')
            if 'value3' in metric_dict[ref].keys():
                print('\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' +
                      str(metric_dict[ref]['value3']) + ', error = ' + str(metric_dict[ref]['value_error3']) +
                      '\033[0m')

