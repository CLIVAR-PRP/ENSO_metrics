# -*- coding:UTF-8 -*-

from getpass import getuser as GETPASSgetuser
import json
from os.path import join as OSpath__join

# ENSO_metrics package
from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection

# set of functions to find cmip/obs files and save a json file
# to be adapted/changed by users depending on their environments
from driver_tools_lib import find_members, find_xml_cmip, find_xml_obs, save_json


# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()


# ---------------------------------------------------#
#
# param needed to compute
#
# path where to save data
path_netcdf = "/data/" + user_name + "/ENSO_metrics/v20200311"

# metric collection
mc_name = "ENSO_perf"

# project / models
project = "CMIP5"  # "CMIP6"
experiment = "hist"
frequency = "mon"
realm = "A"
list_models = ['IPSL-CM5B-LR']
first_member_only = True

# list of observations
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
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc["metrics_list"].keys())


#
# list of variables needed for the given metric collection
#
list_variables = list()
for metric in list_metric:
    listvar = dict_mc["metrics_list"][metric]["variables"]
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)


#
# finding file and variable name in file for each observational dataset
#
dict_obs = dict()
for obs in list_obs:
    # be sure that the datasets is defined in EnsoCollectionsLib.ReferenceObservations
    dict_var = ReferenceObservations(obs)["variable_name_in_file"]
    dict_obs[obs] = dict()
    for var in list_variables:
        #
        # finding variable name in file
        #
        try: var_in_file = dict_var[var]["var_name"]
        except:
            print(str(var) + " is not available for " + str(obs) + " or unscripted")
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
                    file_name, file_areacell, file_landmask = find_xml_obs(obs, var1)
                    list_files.append(file_name)
                    list_areacell.append(file_areacell)
                    list_name_area.append(areacell_in_file)
                    list_landmask.append(file_landmask)
                    list_name_land.append(landmask_in_file)
            else:
                file_name, file_areacell, file_landmask = find_xml_obs(obs, var_in_file)
                list_files = file_name
                list_areacell = file_areacell
                list_name_area = areacell_in_file
                list_landmask = file_landmask
                list_name_land = landmask_in_file
            dict_obs[obs][var] = {"path + filename": list_files, "varname": var_in_file,
                                  "path + filename_area": list_areacell, "areaname": list_name_area,
                                  "path + filename_landmask": list_landmask, "landmaskname": list_name_land}


#
# finding file and variable name in file for each models
#
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()["variable_name_in_file"]
for mod in list_models:
    list_ens = find_members(experiment, frequency, mod, project, realm, first_only=first_member_only)
    pattern_out = OSpath__join(path_netcdf, user_name + "_" + mc_name + "_" + mod + "_" + experiment)
    dict_ens, dict_ens_dive = dict(), dict()
    for ens in list_ens:
        dict_mod = {mod + '_' + ens: {}}
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
                        find_xml_cmip(experiment, frequency, mod, project, realm, ens, var1)
                    list_files.append(file_name)
                    list_areacell.append(file_areacell)
                    list_name_area.append(areacell_in_file)
                    list_landmask.append(file_landmask)
                    list_name_land.append(landmask_in_file)
            else:
                file_name, file_areacell, file_landmask = \
                    find_xml_cmip(experiment, frequency, mod, project, realm, ens, var_in_file)
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
        dictDatasets = {"model": dict_mod, "observations": dict_obs}
        # Computes the metric collection
        netcdf = pattern_out + "_" + ens
        dict_ens[mod + "_" + ens], dict_ens_dive[mod + "_" + ens] =\
            ComputeCollection(mc_name, dictDatasets, mod + "_" + ens, netcdf=True, netcdf_name=netcdf, debug=False)
        # save json
        save_json({mod + "_" + ens: dict_ens[mod + "_" + ens]}, netcdf, metric_only=True)
        with open(netcdf + "_raw.json", "w") as outfile:
            json.dump(dict_ens[mod + "_" + ens], outfile, sort_keys=True)
        del dict_mod, dictDatasets, netcdf
    dict_metric[mod], dict_dive[mod] = dict_ens, dict_ens_dive
    del dict_ens, dict_ens_dive, list_ens, pattern_out
