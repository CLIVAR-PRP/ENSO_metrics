# -*- coding:UTF-8 -*-
from __future__ import print_function
from copy import deepcopy
from getpass import getuser as GETPASSgetuser
from cdms2 import open as CDMS2open
from inspect import stack as INSPECTstack
import json
from os import environ as OSenviron
from os.path import join as OSpath__join
from sys import exit as SYSexit
from sys import path as SYSpath
# ENSO_metrics package
from EnsoMetrics.EnsoCollectionsLib import ReferenceObservations

# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# path
xmldir = OSenviron['XMLDIR']
path_obs = "/data/" + user_name + "/Obs"
path_netcdf = "/data/" + user_name + "/ENSO_metrics/v20200311"

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
        if "r1i1p1" in members:
            members = "r1i1p1"
        elif "r1i1p1f1" in members:
            members = "r1i1p1f1"
        elif "r1i1p1f2" in members:
            members = "r1i1p1f2"
        else:
            tmp = deepcopy(members)
            members = list()
            for mem in tmp:
                for ii in range(1, 10):
                    if "r"+str(ii)+"i" in mem:
                        members.append(mem.replace("r"+str(ii)+"i", "r"+str(ii).zfill(2)+"i"))
                    else:
                        members.append(mem)
            members = sorted(list(set(members)), key=lambda v: v.upper())[0].replace("r0", "r")
    return members


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
    list_obs = ReferenceObservations().keys()
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
    if ".json" not in json_name:
        json_name += ".json"
    with open(json_name, "w") as outfile:
        json.dump(dict_out, outfile, sort_keys=True)
    return
