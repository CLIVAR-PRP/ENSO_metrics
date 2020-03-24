# -*- coding:UTF-8 -*-
from __future__ import print_function
from copy import deepcopy
from getpass import getuser as GETPASSgetuser
from glob import iglob as GLOBiglob
from cdms2 import open as CDMS2open
from inspect import stack as INSPECTstack
import json
from os import environ as OSenviron
from os.path import join as OSpath__join
from sys import argv as SYSargv
from sys import exit as SYSexit
from sys import path as SYSpath

# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# path
xmldir = OSenviron['XMLDIR']
path_obs = "/data/" + user_name + "/Obs"
path_netcdf = "/data/" + user_name + "/ENSO_metrics/v20200311"
# ENSO_metrics package
SYSpath.insert(0, "/home/" + user_name + "/Test001/ENSO_metrics/lib")
SYSpath.insert(1, "/home/" + user_name + "/Test001/ENSO_metrics/plots")
SYSpath.insert(2, "/home/" + user_name + "/Test001/ENSO_metrics/scripts")
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection

# My (YYP) package
# set new path where to find programs
SYSpath.insert(0, "/home/yplanton/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import find_path_and_files
from getfiles_sh_to_py import get_ensembles


# ---------------------------------------------------#
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
# ---------------------------------------------------#


def find_xml(name, frequency, variable, project='', experiment='', ensemble='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_name, file_area, file_land = find_xml_obs(name, frequency, variable)
    else:
        file_name, file_area, file_land = find_xml_cmip(name, project, experiment, ensemble, frequency, realm, variable)
    return file_name, file_area, file_land


def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                              rea=realm, var=variable)
    except:
        if realm == "O":
            new_realm = "A"
        elif realm == "A":
            new_realm = "O"
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere), and conversely
        try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                                  rea=new_realm, var=variable)
        except:
            # given variable is neither in realm 'A' nor 'O'
            print(bcolors.FAIL + "%%%%%     -----     %%%%%")
            print("ERROR: function: " + str(INSPECTstack()[0][3]) + ", line: " + str(INSPECTstack()[0][2]))
            print("given variable cannot be found in either realm A or O: " + str(variable))
            print("param: " + str(model) + ", " + str(project) + ", " + str(experiment) + ", " + str(ensemble) +
                  ", " + str(frequency) + ", " + str(realm))
            print("%%%%%     -----     %%%%%" + bcolors.ENDC)
            SYSexit("")
        file_area, file_land =\
            find_fx(model, project=project, experiment=experiment, ensemble=ensemble, realm=new_realm)
    else:
        file_area, file_land = find_fx(model, project=project, experiment=experiment, ensemble=ensemble, realm=realm)
    file_name = OSpath__join(pathnc, filenc[0])
    return file_name, file_area, file_land


def find_fx(model, project='', experiment='', ensemble='', realm=''):
    if project in ["CMIP5", "CMIP6"]:
        if project in ['CMIP5']:
            my_ens = "r0i0p0"
        else:
            my_ens = deepcopy(ensemble)
        if realm == "A":
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre="fx", mod=model, pro=project,
                                                 rea=realm, var="areacella")
            fland1, fland2 = find_path_and_files(ens=my_ens, exp=experiment, fre="fx", mod=model, pro=project,
                                                 rea=realm, var="sftlf")
            file_land = OSpath__join(fland1, fland2[0])
        elif realm == "O":
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre="fx", mod=model, pro=project,
                                                 rea=realm, var="areacello")
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
        if user_name == "yplanton":
            file_area = None
            if name in ["20CRv2", "NCEP2"]:
                file_land = OSpath__join(path_obs, name + "/land.sfc.gauss.nc")
            elif name == "ERSSTv5":
                file_land = OSpath__join(path_obs, name + "/lsmask_ERSSTv5.nc")
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
        file_area = OSpath__join(xmldir, str(name) + "_" + str(project) + "_" + str(experiment) + "_r0i0p0_glob_fx_"
                                 + str(realm) + "_areacell.xml")
        file_land = OSpath__join(xmldir, str(name) + "_" + str(project) + "_" + str(experiment) + "_r0i0p0_glob_fx_"
                                 + str(realm) + "_landmask.xml")
    return file_area, file_land


def find_xml_obs(obs, frequency, variable):
    # if obs == 'HadISST':
    #     file_name = OSpath__join(xmldir, 'obs_' + str(obs) + 'v1.1.xml')
    # else:
    #     file_name = OSpath__join(xmldir, 'obs_' + str(obs) + '.xml')
    file_name = OSpath__join(xmldir, "obs_ENSO_metrics_" + str(obs) + ".xml")
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print(bcolors.FAIL + "%%%%%     -----     %%%%%")
        print(str().ljust(5) + "obs var " + str(variable) + " cannot be found")
        print(str().ljust(10) + "file_name = " + str(file_name))
        print(str().ljust(10) + "variables = " + str(listvar1))
        print("%%%%%     -----     %%%%%" + bcolors.ENDC)
        SYSexit("")
    file_area, file_land = find_fx(obs)
    return file_name, file_area, file_land


def save_json(dict_in, json_name, metric_only=True):
    # reshape dictionary
    liste = sorted(dict_in.keys())
    listm = sorted(dict_in[liste[0]]['value'].keys())
    dict_out = dict()
    for met in listm:
        dict1 = dict()
        for ens in liste:
            # metadata (nyears)
            dict_meta = dict()
            for key1 in dict_in[ens]['metadata']['metrics'][met]['diagnostic'].keys():
                if key1 not in ['time_frequency', 'ref', 'method', 'method_nonlinearity', 'name']:
                    if key1 == "units":
                        dict_meta[key1] = dict_in[ens]['metadata']['metrics'][met]['diagnostic'][key1]
                    else:
                        dict_meta[key1] = dict_in[ens]['metadata']['metrics'][met]['diagnostic'][key1]['nyears']
            units = dict_in[ens]['metadata']['metrics'][met]['metric']['units']
            if metric_only is True:
                # metrics
                dict2 = dict()
                for key1 in dict_in[ens]['value'][met]['metric'].keys():
                    tmp = dict_in[ens]['value'][met]['metric'][key1]['value']
                    tmp_key = key1.replace("ref_", "")
                    dict2[tmp_key] = {'metric': tmp, 'nyears_obs': dict_meta[tmp_key], 'units': units}
                    del tmp, tmp_key
            else:
                # metrics
                dict2 = {'metric': {}, 'diagnostic': {}}
                for key1 in dict_in[ens]['value'][met]['metric'].keys():
                    tmp = dict_in[ens]['value'][met]['metric'][key1]['value']
                    tmp_key = key1.replace("ref_", "")
                    dict2['metric'][tmp_key] = {'value': tmp, 'nyears_obs': dict_meta[tmp_key], 'units': units}
                    del tmp, tmp_key
                # dive down diagnostics
                for key1 in dict_in[ens]['value'][met]['diagnostic'].keys():
                    tmp = dict_in[ens]['value'][met]['diagnostic'][key1]['value']
                    if key1 == 'model':
                        dict2['diagnostic'][ens] = \
                            {'value': tmp, 'nyears': dict_meta[key1], 'units': dict_meta['units']}
                    else:
                        dict2['diagnostic'][key1] = \
                            {'value': tmp, 'nyears': dict_meta[key1], 'units': dict_meta['units']}
                    del tmp
            dict1[ens] = dict2
            del dict_meta, dict2
        dict_out[met] = dict1
        del dict1
    # save as json file
    with open(json_name + '.json', 'w') as outfile:
        json.dump(dict_out, outfile, sort_keys=True)
    return


# metric collection
mc_name = "test_tel"
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc["metrics_list"].keys())

# parameters
project = "CMIP5"  # "CMIP6"
experiment = "hist"
frequency = "mon"
realm = "A"

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc["metrics_list"][metric]["variables"]
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print(bcolors.HEADER + str(list_variables) + bcolors.ENDC)

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
if mc_name == 'ENSO_perf':
    list_obs = ['20CRv2', 'AVISO', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'GPCPv2.3', 'HadISST', 'NCEP2',
                'SODA3.4.2', 'Tropflux']
elif mc_name == 'ENSO_tel':
    list_obs = ['20CRv2', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GPCPv2.3', 'HadISST', 'NCEP2', 'SODA3.4.2', 'Tropflux']
elif mc_name == 'ENSO_proc':
    list_obs = ['20CRv2', 'AVISO', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'HadISST', 'NCEP2', 'SODA3.4.2', 'Tropflux']
else:
    list_obs = ['20CRv2', 'AVISO', 'CMAP', 'ERA-Interim', 'ERSSTv5', 'GODAS', 'GPCPv2.3', 'HadISST', 'NCEP2',
                'SODA3.4.2', 'Tropflux']
print(bcolors.HEADER + str(list_obs) + bcolors.ENDC)


#
# finding file and variable name in file for each observations dataset
#
dict_obs = dict()
for obs in list_obs:
    # @jiwoo: be sure to add your datasets to EnsoCollectionsLib.ReferenceObservations if needed
    dict_var = ReferenceObservations(obs)["variable_name_in_file"]
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
        try: var_in_file = dict_var[var]["var_name"]
        except:
            print(bcolors.HEADER + str(var) + " is not available for " + str(obs) + " or unscripted" + bcolors.ENDC)
        else:
            try:
                areacell_in_file = dict_var["areacell"]["var_name"]
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var["landmask"]["var_name"]
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
            dict_obs[obs][var] = {"path + filename": list_files, "varname": var_in_file,
                                  "path + filename_area": list_areacell, "areaname": list_name_area,
                                  "path + filename_landmask": list_landmask, "landmaskname": list_name_land}

# models
list_models = ['IPSL-CM5B-LR']
#
# finding file and variable name in file for each observations dataset
#
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()["variable_name_in_file"]
for mod in list_models:
    # list_ens = get_ensembles(exp=experiment, fre=frequency, mod=mod, pro=project, rea=realm)
    if project == "CMIP5":
        list_ens = ["r1i1p1"]
    else:
        if mod in ["CNRM-CM6-1", "MIROC-ES2L", "UKESM1-0-LL"]:
            list_ens = ["r1i1p1f2"]
        elif mod == "CNRM-ESM2-1":
            list_ens = ["r2i1p1f2"]
        elif mod == "HadGEM3-GC31-LL":
            list_ens = ["r1i1p1f3"]
        else:
            list_ens = ["r1i1p1f1"]
    pattern_out = OSpath__join(path_netcdf, user_name + "_" + mc_name + "_" + mod + "_" + experiment)
    files_in = list(GLOBiglob(pattern_out + "*.json"))
    list_ens = list_ens[len(files_in):]
    dict_ens, dict_ens_dive = dict(), dict()
    for ens in list_ens:
        dict_mod = {mod + '_' + ens: {}}
        # ------------------------------------------------
        # @jiwoo: between these dash the program is a bit ad hoc...
        # it works well for me because I am looking for sst and taux on the ocean grid, and fluxes
        # [lhf, lwr, swr, shf, thf] on the atmosphere grid
        # if you want to use atmosphere only, do not use this or create your own way to find the equivalent between the
        # variable name in the program and the variable name in the file
        for var in list_variables:
            #
            # finding variable name in file
            #
            var_in_file = dict_var[var]["var_name"]
            try:
                areacell_in_file = dict_var["areacell"]["var_name"]
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var["landmask"]["var_name"]
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
            dict_mod[mod + '_' + ens][var] =\
                {"path + filename": list_files, "varname": var_in_file, "path + filename_area": list_areacell,
                 "areaname": list_name_area, "path + filename_landmask": list_landmask, "landmaskname": list_name_land}
            del areacell_in_file, file_areacell, file_landmask, file_name, landmask_in_file, list_areacell, list_files,\
                list_landmask, list_name_area, list_name_land, var_in_file
            # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
            # @jiwoo the ComputeCollection function it still on development and it does not read the observations
            # requirement defined in the metric collection,
            # i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
            # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must
            # include every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its
            # coming from another dataset
        dictDatasets = {"model": dict_mod, "observations": dict_obs}
        # regridding dictionary (only if you want to specify the regridding)
        dict_regrid = {}
        # dict_regrid = {
        #     'regridding': {
        #         'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
        #         'newgrid_name': 'generic 1x1deg'},
        # }
        # Computes the metric collection
        netcdf_name = user_name + "_" + mc_name + "_" + mod + "_" + experiment + "_" + ens
        netcdf = pattern_out + "_" + ens
        dict_ens[mod + "_" + ens], dict_ens_dive[mod + "_" + ens] =\
            ComputeCollection(mc_name, dictDatasets, mod + "_" + ens, netcdf=True, netcdf_name=netcdf, debug=False)
        # save json
        save_json({mod + "_" + ens: dict_ens[mod + "_" + ens]}, netcdf, metric_only=True)
        with open(netcdf + "_raw.json", "w") as outfile:
            json.dump(dict_ens[mod + "_" + ens], outfile, sort_keys=True)
        del dict_mod, dict_regrid, dictDatasets, netcdf, netcdf_name
    dict_metric[mod], dict_dive[mod] = dict_ens, dict_ens_dive
    del dict_ens, dict_ens_dive, files_in, list_ens, pattern_out
print(str(SYSargv[0]) + ": done")
