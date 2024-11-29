# -*- coding:UTF-8 -*-

from copy import deepcopy
from getpass import getuser as GETPASSgetuser
from cdms2 import open as CDMS2open
from inspect import stack as INSPECTstack
import json
from numpy import array as NUMPYarray
from os import environ as OSenviron
from os.path import join as OSpath__join
from sys import exit as SYSexit
from sys import path as SYSpath

# ENSO_metrics package
from EnsoMetrics.EnsoCollectionsLib import ReferenceObservations
from EnsoPlots.EnsoPlotToolsLib import find_first_member, get_reference, remove_metrics, sort_members

# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# path
xmldir = OSenviron['XMLDIR']
path_obs = "/data/" + user_name + "/Obs"
path_netcdf = "/data/" + user_name + "/ENSO_metrics/v20200311"

# My (YYP) package
# set new path where to find programs
# SYSpath.insert(0, "/home/yplanton/New_programs/lib_cmip_bash")
# from getfiles_sh_to_py import find_path_and_files
# from getfiles_sh_to_py import get_ensembles


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


def find_members(experiment, frequency, model, project, realm, first_only=False):
    """
    Finds member names

    Inputs:
    ------
    :param experiment: string
        experiment name (e.g., "historical", "piControl")
    :param frequency: string
        data frequency: "day" for daily, "mon" for monthly
    :param model: string
        model name (e.g., "CNRM-CM5", "IPSL-CM5A-LR")
    :param project: string
        project name (e.g., "CMIP5", "CMIP6")
    :param realm: string
        data realm: "A" for atmosphere, "O" for ocean
    **Optional arguments:**
    :param first_only: boolean, optional
        True to return only the first member

    Output:
    ------
    :return members: list
        list of member(s) for the given information
    """
    members = get_ensembles(exp=experiment, fre=frequency, mod=model, pro=project, rea=realm)
    if first_only is True:
        members = [find_first_member(members)]
    return members


def find_fx(model, experiment='', project='', realm='', ensemble=''):
    """
    Finds fixed variables, here areacell and sftlf (landmask)

    Inputs:
    ------
    :param model: string
        model name (e.g., "CNRM-CM5", "IPSL-CM5A-LR")
    **Optional arguments:**
    :param experiment: string, optional
        experiment name (e.g., "historical", "piControl")
    :param project: string, optional
        project name (e.g., "CMIP5", "CMIP6")
    :param realm: string, optional
        data realm: "A" for atmosphere, "O" for ocean
    :param ensemble: string, optional
        ensemble name (e.g., "r1i1p1", "r1i1p1f1")

    Outputs:
    -------
    :return file_area: string
        Path and areacell file name corresponding to the given information (e.g., /path/to/file/areacell.xml)
        Set to None if the file cannot be found
    :return file_land: string
        Path and landmask file name corresponding to the given information (e.g., /path/to/file/landmask.xml)
        Set to None if the file cannot be found
    """
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
        else:
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


def find_xml_cmip(experiment, frequency, model, project, realm, ensemble, variable):
    """
    Finds cmip variable file, as well as corresponding areacell and landmask

    Inputs:
    ------
    :param experiment: string
        experiment name (e.g., "historical", "piControl")
    :param frequency: string
        data frequency: "day" for daily, "mon" for monthly
    :param model: string
        model name (e.g., "CNRM-CM5", "IPSL-CM5A-LR")
    :param project: string
        project name (e.g., "CMIP5", "CMIP6")
    :param realm: string
        data realm: "A" for atmosphere, "O" for ocean
    :param ensemble: string
        ensemble name (e.g., "r1i1p1", "r1i1p1f1")
    :param variable: string
        variable name (e.g., "pr", "tos")

    Outputs:
    -------
    :return file_name: string
        Path and file name corresponding to the given information (e.g., /path/to/file/filename.xml)
    :return file_area: string
        Path and areacell file name corresponding to the given information (e.g., /path/to/file/areacell.xml)
        Set to None if the file cannot be found
    :return file_land: string
        Path and landmask file name corresponding to the given information (e.g., /path/to/file/landmask.xml)
        Set to None if the file cannot be found
    """
    try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                              rea=realm, var=variable)
    except:
        if realm == "O":
            new_realm = "A"
        else:
            new_realm = "O"
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere), and conversely
        try: pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                                  rea=new_realm, var=variable)
        except:
            pathnc, filenc = None, [None]
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
    file_name = OSpath__join(pathnc, str(filenc[0]))
    return file_name, file_area, file_land


def find_xml_obs(dataset, variable):
    """
    Finds observational variable file, as well as corresponding areacell and landmask

    Inputs:
    ------
    :param dataset: string
        model name (e.g., "CNRM-CM5", "IPSL-CM5A-LR")
    :param variable: string
        variable name (e.g., "pr", "tos")

    Outputs:
    -------
    :return file_name: string
        Path and file name corresponding to the given information (e.g., /path/to/file/filename.xml)
    :return file_area: string
        Path and areacell file name corresponding to the given information (e.g., /path/to/file/areacell.xml)
        Set to None if the file cannot be found
    :return file_land: string
        Path and landmask file name corresponding to the given information (e.g., /path/to/file/landmask.xml)
        Set to None if the file cannot be found
    """
    file_name = OSpath__join(xmldir, "obs_ENSO_metrics_" + str(dataset) + ".xml")
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print(bcolors.FAIL + "%%%%%     -----     %%%%%")
        print(str().ljust(5) + "obs var " + str(variable) + " cannot be found")
        print(str().ljust(10) + "file_name = " + str(file_name))
        print(str().ljust(10) + "variables = " + str(listvar1))
        print("%%%%%     -----     %%%%%" + bcolors.ENDC)
        SYSexit("")
    file_area, file_land = find_fx(dataset)
    return file_name, file_area, file_land


def find_xml_fx(model, experiment='', project='', realm=''):
    """
    Finds fixed variables, here areacell and sftlf (landmask), mostly used for observational dataset

    Inputs:
    ------
    :param model: string
        model name (e.g., "CNRM-CM5", "IPSL-CM5A-LR")
    **Optional arguments:**
    :param experiment: string, optional
        experiment name (e.g., "historical", "piControl")
    :param project: string, optional
        project name (e.g., "CMIP5", "CMIP6")
    :param realm: string, optional
        data realm: "A" for atmosphere, "O" for ocean

    Outputs:
    -------
    :return file_area: string
        Path and areacell file name corresponding to the given information (e.g., /path/to/file/areacell.xml)
        Set to None if the file cannot be found
    :return file_land: string
        Path and landmask file name corresponding to the given information (e.g., /path/to/file/landmask.xml)
        Set to None if the file cannot be found
    """
    list_obs = list(ReferenceObservations().keys())
    if model in list_obs:
        if user_name == "yplanton":
            file_area = None
            if model in ["20CRv2", "NCEP2"]:
                file_land = OSpath__join(path_obs, model + "/land.sfc.gauss.nc")
            elif model == "ERSSTv5":
                file_land = OSpath__join(path_obs, model + "/lsmask_ERSSTv5.nc")
            elif model == "CMAP":
                file_land = OSpath__join(path_obs, model + "/lsmask_fx_cmap.nc")
            elif model == "GPCPv2.3":
                file_land = OSpath__join(path_obs, model + "/lsmask_fx_gpcpv2.3.nc")
            elif model == "OISSTv2":
                file_land = OSpath__join(path_obs, model + "/lsmask_fx_oisstv2.nc")
            else:
                file_land = None
        else:
            file_area = None
            if model in ["20CRv2", "CMAP", "ERA-Interim", "ERSSTv5", "GPCPv2.3", "NCEP2"]:
                file_land = OSpath__join(path_obs, "lsmask_" + model + ".nc")
            else:
                file_land = None
    else:
        file_area = OSpath__join(xmldir, str(model) + "_" + str(project) + "_" + str(experiment) + "_r0i0p0_glob_fx_"
                                 + str(realm) + "_areacell.xml")
        file_land = OSpath__join(xmldir, str(model) + "_" + str(project) + "_" + str(experiment) + "_r0i0p0_glob_fx_"
                                 + str(realm) + "_landmask.xml")
    return file_area, file_land


def get_metric_values(project, metric_collection, dict_json, dict_mod_mem, reduced_set=True, portraitplot=False):
    """
    Finds fixed variables, here areacell and sftlf (landmask), mostly used for observational dataset

    Inputs:
    ------
    :param project: strings
        project name (e.g., "CMIP5", "CMIP6")
    :param metric_collection: strings
        metric collection (e.g., "ENSO_perf", "ENSO_proc", "ENSO_tel")
    :param dict_json: dictionary
        Dictionary with path and name of json files output of the CLIVAR PRP ENSO recipes package.
    :param dict_mod_mem: dictionary
        Dictionary with every models available and members in the given json files, for the given projects and metric
        collections.
    **Optional arguments:**
    :param reduced_set: boolean, optional
        True to remove extra recipes that are not in the final set chosen by CLIVAR PRP.
        If set to False it removes recipes that are in more than one metric collection.
        Default value is True.
    :param portraitplot: boolean, optional
        True to remove extra recipes that are not in the final set chosen by CLIVAR PRP but keep recipes that are in
        more than one metric collection.
        If set to False it removes recipes that are in more than one metric collection.
        Default value is False.

    Output:
    ------
    :return dict_out: dictionary
        Dictionary with every models available and recipes, member values averaged
    """

    # open and read json file
    data_json = read_json(dict_json[project][metric_collection])
    list_models = list(dict_mod_mem[project].keys())
    dict_out = dict()
    for mod in list_models:
        list_members = sort_members(dict_mod_mem[project][mod])
        dict2 = dict()
        for mem in list_members:
            try:    data_mod = data_json[mod][mem]["value"]
            except: data_mod = None
            list_metrics = list()
            try:    list(data_mod.keys())
            except: pass
            else:   list_metrics += list(data_mod.keys())
            list_metrics = remove_metrics(list_metrics, metric_collection, reduced_set=reduced_set,
                                          portraitplot=portraitplot)
            dict3 = dict()
            if len(list_metrics) > 0:
                for met in list_metrics:
                    dict3[met] = data_mod[met]["metric"][get_reference(metric_collection, met)]["value"]
                dict2[mem] = dict3
            del data_mod, dict3, list_metrics
        # models and recipes
        list_metrics = list()
        for mem in list_members:
            try:    list(dict2[mem].keys())
            except: pass
            else:   list_metrics += list(dict2[mem].keys())
        # average member values if there is more than one available
        dict_met = dict()
        for met in list_metrics:
            tmp = [dict2[mem][met] for mem in list_members
                   if dict2[mem][met] is not None and dict2[mem][met] != 1e20]
            if len(tmp) > 0:
                dict_met[met] = float(tmp[0]) if len(tmp) == 1 else float(NUMPYarray(tmp).mean())
            del tmp
        dict_out[mod] = dict_met
        del dict2, dict_met, list_members, list_metrics
    return dict_out


def get_metric_values_observations(filename_json, obsvation_names, list_met, metric_collection):
    """
    Reads given json file (must have usual jiwoo's structure) and read given obs

    Inputs:
    ------
    :param filename_json: string
        Path and name of a json file output of the CLIVAR 2020 ENSO recipes package.
    :param obsvation_names: list of string
        Names of wanted additional observations for the portrait plot.
    :param list_met: list of string
        List of recipes.
    :param metric_collection: string
        Name of a metric collection.

    Output:
    ------
    :return data: list
        Dictionary output of additional observations metric values.
    """
    data_json = read_json(filename_json)
    dict_out = dict()
    for obs in obsvation_names:
        for met in list_met:
            ref = get_reference(metric_collection, met)
            if obs == "20CRv2":
                if "Ssh" not in met:
                    try:
                        tab = data_json["20CRv2"]["r1i1p1"]["value"][met]["metric"]
                    except:
                        tab = data_json["20CRv2_20CRv2"]["r1i1p1"]["value"][met]["metric"]
            elif obs == "NCEP2":
                if "TauxSsh" in met or "SshSst" in met:
                    tab = data_json["NCEP2_GODAS"]["r1i1p1"]["value"][met]["metric"]
                elif "Ssh" in met:
                    tab = data_json["GODAS"]["r1i1p1"]["value"][met]["metric"]
                else:
                    try:
                        tab = data_json["NCEP2"]["r1i1p1"]["value"][met]["metric"]
                    except:
                        tab = data_json["NCEP2_NCEP2"]["r1i1p1"]["value"][met]["metric"]
            elif obs == "ERA-Interim":
                if "SstMap" in met:
                    tab = {ref: {"value": 0}}
                elif "TauxSsh" in met or "SshSst" in met:
                    tab = data_json["ERA-Interim_SODA3.4.2"]["r1i1p1"]["value"][met]["metric"]
                elif "Ssh" in met:
                    tab = data_json["SODA3.4.2"]["r1i1p1"]["value"][met]["metric"]
                else:
                    try:
                        tab = data_json["ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
                    except:
                        tab = data_json["ERA-Interim_ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
            try:
                val = tab[ref]["value"]
            except:
                val = 1e20
            try:
                dict_out[obs]
            except:
                dict_out[obs] = {met: val}
            else:
                dict_out[obs][met] = val
            try:
                del tab
            except:
                pass
            del ref, val
    return dict_out


def get_mod_mem_json(projects, metric_collections, dict_json, first_only=False):
    """
    Creates a dictionary with every models available in the given json files, for the given projects and metric
    collections.
    Also provides a list of available members. If first_only is True, provides only the first available member.

    dict_out = {'project1': {'model1': ['member1', 'member2', ...],
                             'model2': ['member1', 'member2', ...],
                             ...},
                'project2': {'model1': ['member1', 'member2', ...],
                             'model2': ['member1', 'member2', ...],
                             ...},
                }

    Inputs:
    ------
    :param projects: list of strings
        list project names (e.g., "CMIP5", "CMIP6")
    :param metric_collections: list of strings
        list of metric collections (e.g., "ENSO_perf", "ENSO_proc", "ENSO_tel")
    :param dict_json: dictionary
        Dictionary with path and name of json files output of the CLIVAR PRP ENSO recipes package.
    **Optional arguments:**
    :param first_only: boolean, optional
        True to return only the first member

    Output:
    ------
    :return model_by_proj: dictionary
        Dictionary with every models available and members in the given json files, for the given projects and metric
        collections.
        If first_only is True, provides only the first available member.
    """
    # get members by model by project from json file
    # only recipes from models/members chosen here will be used
    # all recipes from models/members chosen here will be used (ensures that if a model/member is not available for one
    # or several metric collections, the corresponding line will still be created in the portraitplot)
    model_by_proj = dict()
    for proj in projects:
        list_models = list()
        dict_members = dict()
        for mc in metric_collections:
            # read json files
            tmp = read_json(dict_json[proj][mc])
            # list models
            list_models += list(tmp.keys())
            # members
            for mod in list(tmp.keys()):
                try:
                    dict_members[mod]
                except:
                    dict_members[mod] = list(tmp[mod].keys())
                else:
                    dict_members[mod] += list(tmp[mod].keys())
            del tmp
        list_models = sorted(list(set(list_models)), key=lambda v: v.upper())
        list_to_remove = ["EC-EARTH", "FIO-ESM", "GFDL-CM2p1", "HadGEM2-AO", "CIESM", "E3SM-1-1-ECA", "FGOALS-g3",
                          "MCM-UA-1-0"]
        # EC-EARTH: incorrect time coordinate
        # FIO-ESM, HadCM3: grid issues
        # GFDL-CM2p1: hfls not published
        # HadGEM2-AO: rlus and rsus not published
        # E3SM-1-1-ECA: Experimental stage
        # CIESM, FGOALS-g3: ???
        # MCM-UA-1-0: unit issue with pr
        for mod in list_to_remove:
            while mod in list_models:
                list_models.remove(mod)
        for mod in list_models:
            list_members = sorted(list(set(dict_members[mod])), key=lambda v: v.upper())
            if first_only is True:
                list_members = [find_first_member(list_members)]
            try:
                model_by_proj[proj]
            except:
                model_by_proj[proj] = {mod: list_members}
            else:
                try:
                    model_by_proj[proj][mod]
                except:
                    model_by_proj[proj][mod] = list_members
                else:
                    print("this model should not be here")
            del list_members
        del dict_members, list_models, list_to_remove
    return model_by_proj


def read_json(filename_json):
    """
    Reads given json file (must have usual jiwoo's structure)

    Input:
    -----
    :param filename_json: string
        Path and name of a json file output of the CLIVAR PRP ENSO recipes package.

    Output:
    ------
    :return data: dictionary
        Dictionary output of the CLIVAR PRP ENSO recipes package, first level is models, second is members.
    """
    with open(filename_json) as ff:
        data = json.load(ff)
    ff.close()
    data = data["RESULTS"]["model"]
    return data


def save_json(dict_in, json_name, metric_only=True):
    """
    Saves given dictionary under given name in a json file

    Inputs:
    ------
    :param dict_in: dictionary
        data to save in a json file
    :param json_name: string
        Path and file name where to save the given data (e.g., /path/to/file/jsonname.json)
    **Optional arguments:**
    :param metric_only: boolean, optional
        True to save only the metric values

    Output:
    ------
    :return:
    """
    # reshape dictionary
    liste = sorted(dict_in.keys())
    listm = sorted(dict_in[liste[0]]['value'].keys())
    dict_out = dict()
    for met in listm:
        dict1 = dict()
        for ens in liste:
            # metadata (nyears)
            dict_meta = dict()
            for key1 in list(dict_in[ens]['metadata']['recipes'][met]['diagnostic'].keys()):
                if key1 not in ['time_frequency', 'ref', 'method', 'method_nonlinearity', 'name']:
                    if key1 == "units":
                        dict_meta[key1] = dict_in[ens]['metadata']['recipes'][met]['diagnostic'][key1]
                    else:
                        dict_meta[key1] = dict_in[ens]['metadata']['recipes'][met]['diagnostic'][key1]['nyears']
            units = dict_in[ens]['metadata']['recipes'][met]['metric']['units']
            if metric_only is True:
                # recipes
                dict2 = dict()
                for key1 in list(dict_in[ens]['value'][met]['metric'].keys()):
                    tmp = dict_in[ens]['value'][met]['metric'][key1]['value']
                    tmp_key = key1.replace("ref_", "")
                    dict2[tmp_key] = {'metric': tmp, 'nyears_obs': dict_meta[tmp_key], 'units': units}
                    del tmp, tmp_key
            else:
                # recipes
                dict2 = {'metric': {}, 'diagnostic': {}}
                for key1 in list(dict_in[ens]['value'][met]['metric'].keys()):
                    tmp = dict_in[ens]['value'][met]['metric'][key1]['value']
                    tmp_key = key1.replace("ref_", "")
                    dict2['metric'][tmp_key] = {'value': tmp, 'nyears_obs': dict_meta[tmp_key], 'units': units}
                    del tmp, tmp_key
                # dive down diagnostics
                for key1 in list(dict_in[ens]['value'][met]['diagnostic'].keys()):
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
    if ".json" not in json_name:
        json_name += ".json"
    with open(json_name, "w") as outfile:
        json.dump(dict_out, outfile, sort_keys=True)
    return
