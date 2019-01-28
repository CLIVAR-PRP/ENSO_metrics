from cdms2 import open as CDMS2open
from copy import deepcopy
import datetime
from inspect import stack as INSPECTstack
import json
from os import environ as OSenviron
from os.path import join as OSpath__join
import sys

sys.path.insert(0, "/home/yplanton/ENSO_metrics/lib")

# ENSO_metrics package
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection

# My (YYP) package
from getfiles_sh_to_py import find_path_and_files
from getfiles_sh_to_py import get_ensembles

#---------------------------------------------------#
# colors for printing
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#---------------------------------------------------#

today = str(datetime.date.today())
today = today.replace("-", "")
xmldir = OSenviron['XMLDIR']
netcdf_path = '/data/yplanton/ENSO_metrics'


def find_xml(name, frequency, variable, project='', experiment='', ensemble='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_name, file_area, file_land = find_xml_obs(name, frequency, variable)
    else:
        file_name, file_area, file_land = find_xml_cmip(name, project, experiment, ensemble, frequency, realm, variable)
    return file_name, file_area, file_land

def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project, rea=realm, var=variable)
    except:
        if realm == 'O':
            new_realm = 'A'
        elif realm == 'A':
            new_realm = 'O'
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere), and conversely
        try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project, rea=new_realm, var=variable)
        except:
            # given variable is neither in realm 'A' nor 'O'
            print bcolors.FAIL + '%%%%%     -----     %%%%%'
            print 'ERROR: function: '+str(INSPECTstack()[0][3])+', line: '+str(INSPECTstack()[0][2])
            print 'given variable cannot be found in either realm A or O: '+str(variable)
            print 'param: '+str(model)+', '+str(project)+', '+str(experiment)+', '+str(ensemble)+', '+str(frequency)+', '+str(realm)
            print '%%%%%     -----     %%%%%' + bcolors.ENDC
            sys.exit('')
        file_area, file_land = find_fx(model, project=project, experiment=experiment, ensemble=ensemble, realm=new_realm)
    else:
        file_area, file_land = find_fx(model, project=project, experiment=experiment, ensemble=ensemble, realm=realm)
    file_name = OSpath__join(pathnc, filenc[0])
    return file_name, file_area, file_land

def find_fx(model, project='', experiment='', ensemble='', realm=''):
    if project in ['CMIP5', 'CMIP6']:
        if project in ['CMIP5']:
            my_ens = 'r0i0p0'
        else:
            my_ens = deepcopy(ensemble)
        if realm == 'A':
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='areacella')
            fland1, fland2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='sftlf')
            file_land = OSpath__join(fland1, fland2[0])
        elif realm == 'O':
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='areacello')
            file_land = None
        file_area = OSpath__join(farea1, farea2[0])
    else:
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=realm)
    try: CDMS2open(file_area)
    except: file_area = None
    try: CDMS2open(file_land)
    except: file_land = None
    return file_area, file_land

def find_xml_fx(name, project='', experiment='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_area = OSpath__join(xmldir, 'obs_' + str(name) + '_areacell.xml')
        file_land = OSpath__join(xmldir, 'obs_' + str(name) + '_landmask.xml')
    else:
        file_area = OSpath__join(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_'
                                   + str(realm) + '_areacell.xml')
        file_land = OSpath__join(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_'
                                   + str(realm) + '_landmask.xml')
    return file_area, file_land

def find_xml_obs(obs, frequency, variable):
    if obs == 'HadISST':
        file_name = OSpath__join(xmldir, 'obs_' + str(obs) + 'v1.1.xml')
    else:
        file_name = OSpath__join(xmldir, 'obs_' + str(obs) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print '\033[95m' + str().ljust(5) + "obs var " + str(variable) + " cannot be found" + '\033[0m'
        print '\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m'
        print '\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m'
        sys.exit('')
    file_area, file_land = find_fx(obs)
    return file_name, file_area, file_land

def save_json(dict_in, json_name, metric=True):
    # reshape dictionary
    liste = sorted(dict_in.keys())
    dict_out = dict()
    for met in list_metric:
        dict1 = dict()
        for ens in liste:
            if metric is True:
                # metadata (nyears)
                dict3 = dict()
                for key4 in dict_in[ens]['metadata']['metrics'][met]['diagnostic'].keys():
                    if key4 not in ['time_frequency', 'units', 'model', 'ref', 'method', 'name']:
                        dict3[key4] = dict_in[ens]['metadata']['metrics'][met]['diagnostic'][key4]['nyears']
                # metrics
                dict4 = dict()
                for key4 in dict_in[ens]['value'][met]['metric'].keys():
                    tmp = dict_in[ens]['value'][met]['metric'][key4]['value']
                    tmp_key = key4.replace("ref_", "")
                    dict4[tmp_key] = {'metric':tmp, 'nyears_obs':dict3[tmp_key]}
                    del tmp, tmp_key
                dict1[ens] = dict4
                del dict3, dict4
            else:
                # dive down diagnostics
                for key4 in dict_in[ens]['value'][met]['diagnostic'].keys():
                    tmp = dict_in[ens]['value'][met]['diagnostic'][key4]['value']
                    if key4 == 'model':
                        dict1[ens] = tmp
                    else:
                        dict1[key4] = tmp
                    del tmp
        dict_out[met] = dict1
        del dict1
    # save as json file
    with open(json_name + '.json', 'w') as outfile:
        json.dump(dict_out, outfile)
    return

# metric collection
mc_name = 'ENSO_perf'#'ENSO_tel'#'MC1'#
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

# parameters
project = 'CMIP6'
experiment = 'hist'
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
print '\033[95m' + str(list_variables) + '\033[0m'

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
    list_obs = ['HadISST','Tropflux','GPCPv2.3']#['Tropflux','GPCPv2.3']#['HadISST']#
elif mc_name == 'ENSO_tel':
    list_obs = ['HadISST','GPCPv2.3']
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
            print '\033[95m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m'
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'obs', 'var'
            #
            # @jiwoo: pretty easy as I have all variables in one file
            file_name, file_areacell, file_landmask = find_xml(obs, frequency, var0)
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
                list_files = [file_name for var1 in var_in_file]
                list_areacell = [file_areacell for var1 in var_in_file]
                list_name_area = [areacell_in_file for var1 in var_in_file]
                list_landmask = [file_landmask for var1 in var_in_file]
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

# models
list_models = ['IPSL-CM6A-LR']#['CNRM-CM5']#['IPSL-CM5B-LR']#['CNRM-CM5','IPSL-CM5B-LR']#
#
# finding file and variable name in file for each observations dataset
#
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()['variable_name_in_file']
for mod in list_models:
    list_ens = get_ensembles(exp='hist', fre=frequency, mod=mod, pro=project, rea=realm) #toto
    #list_ens = list_ens[4:]
    dict_ens, dict_ens_dive = dict(), dict()
    for ens in list_ens:
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
            file_name, file_areacell, file_landmask = find_xml(mod, frequency, var0, project=project,
                                                               experiment=experiment,
                                                               ensemble=ens, realm=realm)
            try:
                areacell_in_file = dict_var['areacell']['var_name']
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var['landmask']['var_name']
            except:
                landmask_in_file = None
            if isinstance(var_in_file, list):
                list_files = [file_name for var1 in var_in_file]
                list_areacell = [file_areacell for var1 in var_in_file]
                list_name_area = [areacell_in_file for var1 in var_in_file]
                list_landmask = [file_landmask for var1 in var_in_file]
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
            del areacell_in_file, file_areacell, file_landmask, file_name, landmask_in_file, list_areacell, list_files,\
                list_landmask, list_name_area, list_name_land, var0, var_in_file
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
        # dict_metric[mod] = ComputeCollection(mc_name, dictDatasets, user_regridding=dict_regrid, debug=False)
        # dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, user_regridding=dict_regrid,
        #                                                     debug=False, dive_down=True)
        netcdf_name = today + '_YANN_PLANTON_' + mc_name + '_' + mod + '_' + experiment + '_' + ens
        netcdf = OSpath__join(netcdf_path, netcdf_name)
        # dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, netcdf=True, netcdf_path=netcdf_path,
        #                                                      netcdf_name=netcdf_name, debug=True, dive_down=True)
        dict_ens[mod + '__' + ens], dict_ens_dive[mod + '__' + ens] = ComputeCollection(mc_name, dictDatasets,
                                                                                        netcdf=True, netcdf_name=netcdf,
                                                                                        debug=True)
        stop
        # save json
        save_json({mod + '__' + ens: dict_ens[mod + '__' + ens]}, netcdf_name, metric=True)
        del dict_mod, dict_regrid, dictDatasets, netcdf, netcdf_name
    dict_metric[mod], dict_dive[mod] = dict_ens, dict_ens_dive
    del dict_ens, dict_ens_dive, list_ens
# ------------------------------------------------
# reshape dictionary
# listm = sorted(dict_metric.keys())
# liste = dict((mod, sorted(dict_metric[mod].keys())) for mod in listm)
# dict_dive, dict_metr = dict(), dict()
# for met in list_metric:
#     dict1, dict2 = dict(), dict()
#     for mod in listm:
#         for ens in liste[mod]:
#             # dive down diagnostics & metadata (nyears)
#             dict3 = dict()
#             for key4 in dict_metric[mod][ens]['value'][met]['diagnostic'].keys():
#                 # dive
#                 tmp = dict_metric[mod][ens]['value'][met]['diagnostic'][key4]['value']
#                 if key4 == 'model':
#                     dict1[mod + '__' + ens] = tmp
#                 else:
#                     dict1[key4] = tmp
#                 del tmp
#                 # meta
#                 tmp = dict_metric[mod][ens]['metadata']['metrics'][met]['diagnostic'][key4]['nyears']
#                 if key4 != 'model':
#                     dict3[key4] = tmp
#                 del tmp
#             # metrics
#             dict4 = dict()
#             for key4 in dict_metric[mod][ens]['value'][met]['metric'].keys():
#                 tmp = dict_metric[mod][ens]['value'][met]['metric'][key4]['value']
#                 tmp_key = key4.replace("ref_", "")
#                 dict4[tmp_key] = {'metric': tmp, 'nyears_obs': dict3[tmp_key]}
#                 del tmp, tmp_key
#             dict2[mod + '__' + ens] = dict4
#             del dict3, dict4
#     dict_dive[met], dict_metr[met] = dict1, dict2
#     del dict1, dict2
# # ------------------------------------------------
# # save as json file
# name = today + '_YANN_PLANTON_' + mc_name + '.json'
# with open(name, 'w') as outfile:
#     json.dump(dict_metr, outfile)
# del name
