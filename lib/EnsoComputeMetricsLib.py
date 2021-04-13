# -*- coding:UTF-8 -*-
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json

# ENSO_metrics package functions:
from .EnsoCollectionsLib import defCollection, ReferenceObservations
from . import EnsoErrorsWarnings
from .EnsoMetricsLib import BiasPrLatRmse, BiasPrLonRmse, BiasPrMapRmse, BiasLhfLatRmse, BiasLhfLonRmse,\
    BiasLhfMapRmse, BiasLwrLatRmse, BiasLwrLonRmse, BiasLwrMapRmse, BiasShfLatRmse, BiasShfLonRmse, BiasShfMapRmse,\
    BiasSwrLatRmse, BiasSwrLonRmse, BiasSwrMapRmse, BiasSshLatRmse, BiasSshLonRmse, BiasSshMapRmse, BiasSstLatRmse,\
    BiasSstLonRmse, BiasSstMapRmse, BiasSstSkLonRmse, BiasTauxLatRmse, BiasTauxLonRmse, BiasTauxMapRmse,\
    BiasTauyLatRmse, BiasTauyLonRmse, BiasTauyMapRmse, BiasThfLatRmse, BiasThfLonRmse, BiasThfMapRmse, EnsoAmpl,\
    EnsoDiversity, EnsodSstOce, EnsoDuration, EnsoFbSshSst, EnsoFbSstLhf, EnsoFbSstLwr, EnsoFbSstShf, EnsoFbSstSwr,\
    EnsoFbSstTaux, EnsoFbSstThf, EnsoFbTauxSsh, EnsoPrMap, EnsoPrMapDjf, EnsoPrMapJja, EnsoPrDjfTel, EnsoPrJjaTel,\
    EnsoSeasonality, EnsoSlpMap, EnsoSlpMapDjf, EnsoSlpMapJja, EnsoSstDiversity, EnsoPrLonRmse, EnsoSshLonRmse,\
    EnsoSstLonRmse, EnsoTauxLonRmse, EnsoTauyLonRmse, EnsoSstMap, EnsoSstMapDjf, EnsoSstMapJja, EnsoSstSkew,\
    EnsoPrTsRmse, EnsoSshTsRmse, EnsoSstTsRmse, EnsoTauxTsRmse, EnsoTauyTsRmse, grad_lat_pr, grad_lat_sst, grad_lon_pr,\
    grad_lon_sst, NinaPrMap, NinaSlpMap, NinaSstDiv, NinaSstDivRmse, NinaSstDur, NinaSstLonRmse, NinaSstMap,\
    NinaSstTsRmse, NinoPrMap, NinoSlpMap, NinoSstDiv, NinoSstDiversity, NinoSstDivRmse, NinoSstDur, NinoSstLonRmse,\
    NinoSstMap, NinoSstTsRmse, SeasonalPrLatRmse, SeasonalPrLonRmse, SeasonalSshLatRmse, SeasonalSshLonRmse,\
    SeasonalSstLatRmse, SeasonalSstLonRmse, SeasonalTauxLatRmse, SeasonalTauxLonRmse, SeasonalTauyLatRmse,\
    SeasonalTauyLonRmse
from .EnsoToolsLib import math_metric_computation
from .KeyArgLib import default_arg_values


# sst only datasets (not good for surface temperature teleconnection)
sst_only = [
    "C-GLORSv5", "CFSR", "COBE", "COBE1", "COBE-1", "COBEv1", "COBE2", "COBE-2", "COBEv2", "ERSSTv3b", "ERSSTv4",
    "ERSSTv5", "GODAS", "HadISST", "HadISST1", "HadISST-1", "HadISSTv1", "HadISST1.1", "HadISST1-1", "HadISST-1.1",
    "HadISST-1-1", "HadISSTv1.1", "OAFlux", "ORAS4", "ORAS5", "SODA3.3.2" "SODA3.4.2", "SODA3.11.2", "SODA3.12.2",
    "Tropflux", "TropFlux", "Tropflux1", "TropFlux1", "Tropflux-1", "TropFlux-1", "Tropfluxv1", "TropFluxv1",
    "Tropflux1.0", "TropFlux1.0", "Tropflux1-0", "TropFlux1-0", "Tropflux-1.0", "TropFlux-1.0", "Tropflux-1-0",
    "TropFlux-1-0", "Tropfluxv1.0", "TropFluxv1.0"]


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric collection
#
def ComputeCollection(metricCollection, dictDatasets, modelName, user_regridding={}, debug=False, dive_down=False,
                      netcdf=False, netcdf_name="", observed_fyear=None, observed_lyear=None, modeled_fyear=None,
                      modeled_lyear=None, obs_interpreter=None):
    """
    The ComputeCollection() function computes all the diagnostics / metrics associated with the given Metric Collection

    Inputs:
    ------
    :param metricCollection: string
        name of a Metric Collection, must be defined in EnsoCollectionsLib.defCollection()
    :param dictDatasets: dict
        dictionary containing all information needed to compute the Metric Collection for one model and observations
        it must be like:
        dictDatasets = {
            'model': {
                'modelName': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                },
            },
            'observations': {
                # obsName1 can be 'HadISST' if all variables come from HadISST, or it could be 'ERA-Interim + HadISST'
                # if the variables come from ERA-Interim AND HadISST,...
                'obsName1': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                    'variable n': {
                        'path + filename': ['path_to_file1/filename1', 'path_to_file2/filename2', ...],
                        'varname': ['variable_name_in_file1', 'variable_name_in_file2', ...],
                        # this last one is not compulsory if 'obsName1' is defined in
                        # EnsoCollectionsLib.ReferenceObservations()
                        'algebric_calculation': [+1,-1,..., '/'],

                    },
                },
                'obsName2': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                },
            },
        }
    :param user_regridding: dict, optional
        regridding parameters selected by the user and not the "climate experts" (i.e. the program)
        to see parameters:
        help(EnsoUvcdatToolsLib.TwoVarRegrid)
        help(EnsoUvcdatToolsLib.Regrid)
        e.g.:
        user_regridding = {
            'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                'newgrid_name': 'generic 1x1deg'},
        }
    :param debug: boolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param dive_down: boolean, optional
        default value = False dive_down are not saved in a dictionary
        If you want to save the dive down diagnostics set it to True
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' root name of the saved NetCDFs
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='USER_DATE_METRICCOLLECTION_MODEL'
    :param observed_fyear: integer, optional
        first year to use for observational datasets, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param observed_lyear: integer, optional
        last year to use for observational datasets, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param modeled_fyear: integer, optional
        first year to use for CMIP simulations, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param modeled_lyear: integer, optional
        last year to use for CMIP simulations, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'modeled_period' defined in EnsoCollectionsLib.defCollection is used
    :param obs_interpreter: string, optional
        special variable interpreter for all observational datasets
        the only possibility is 'CMIP' to interpret all observational's variables as CMIP (datasets have been CMORized)
        default value = None, observational datasets are considered not CMORized and will be interpreted as defined in
        EnsoCollectionsLib.ReferenceObservations

    :return: MCvalues: dict
        name of the Metric Collection, Metrics, value, value_error, units, ...
        MCvalues = {
            'MetricCollection': {
                'information about the MetricCollection': 'descriptions',
                'metrics': {
                    'metric1': {
                        'metric_values': {
                            'obsName1': {
                                'value': 'value of the metric',
                                'value_error': 'estimation of the error on the metric',
                            },
                            'obsName2': {'value': val, 'value_error': err},
                            ...
                        },
                        'information about the metric': 'description of how this metric if computed from model and
                                                        observations values',
                        'raw_values': {
                            'model': {
                                'value': 'model value of the diagnostic',
                                'value_error': 'estimation of the error on the diagnostic',
                                'nyears': 'number of years used for the computation of the diagnostic',
                                'time_period': 'period used for the computation of the diagnostic',
                            },
                            'observations': {
                                'obsName1': {'value': val, 'value_error': err, 'nyears': ny, 'time_period': bnds},
                                'obsName2': {'value': val, 'value_error': err, 'nyears': ny, 'time_period': bnds},
                                ...
                            },
                            'units': 'units of the diagnostic', 'method': 'method used to compute the diagnostic',
                            'time_frequency': 'data frequency used to compute the diagnostic',
                            'ref': 'reference paper', ...
                        },
                    },
                    'metric2': {
                        ...
                    },
                    ...
                },
            },
        }
    """
    dict_mc = defCollection(metricCollection)
    dict_col_meta = {
        "name": dict_mc["long_name"], "description_of_the_collection": dict_mc["description"], "metrics": {},
    }
    dict_col_dd_meta = {
        "name": dict_mc["long_name"], "description_of_the_collection": dict_mc["description"], "metrics": {},
    }
    dict_col_valu = dict()
    dict_col_dd_valu = dict()
    dict_m = dict_mc["metrics_list"]
    list_metrics = sorted(list(dict_m.keys()), key=lambda v: v.upper())
    for metric in list_metrics:
        try:  # try per metric
            print("\033[94m" + str().ljust(5) + "ComputeCollection: metric = " + str(metric) + "\033[0m")
            # sets arguments for this metric
            list_variables = dict_m[metric]["variables"]
            dict_regions = dict_m[metric]["regions"]
            # model name, file, variable name in file
            try:
                modelFile1 = dictDatasets["model"][modelName][list_variables[0]]["path + filename"]
            except:
                modelFile1 = ""
            try:
                modelVarName1 = dictDatasets["model"][modelName][list_variables[0]]["varname"]
            except:
                modelVarName1 = ""
            try:
                modelFileArea1 = dictDatasets["model"][modelName][list_variables[0]]["path + filename_area"]
            except:
                modelFileArea1, modelAreaName1 = None, None
            else:
                modelAreaName1 = dictDatasets["model"][modelName][list_variables[0]]["areaname"]
            try:
                modelFileLandmask1 = dictDatasets["model"][modelName][list_variables[0]]["path + filename_landmask"]
            except:
                modelFileLandmask1, modelLandmaskName1 = None, None
            else:
                modelLandmaskName1 = dictDatasets["model"][modelName][list_variables[0]]["landmaskname"]
            # observations name(s), file(s), variable(s) name in file(s)
            obsNameVar1, obsFile1, obsVarName1, obsFileArea1, obsAreaName1 = list(), list(), list(), list(), list()
            obsFileLandmask1, obsLandmaskName1, obsInterpreter1 = list(), list(), list()
            for obs in sorted(list(dictDatasets["observations"].keys()), key=lambda v: v.upper()):
                try:
                    dictDatasets["observations"][obs][list_variables[0]]
                except:
                    pass
                else:
                    obsNameVar1.append(obs)
                    obsFile1.append(dictDatasets["observations"][obs][list_variables[0]]["path + filename"])
                    obsVarName1.append(dictDatasets["observations"][obs][list_variables[0]]["varname"])
                    try:
                        obsFileArea1.append(
                            dictDatasets["observations"][obs][list_variables[0]]["path + filename_area"])
                    except:
                        obsFileArea1.append(None)
                        obsAreaName1.append(None)
                    else:
                        obsAreaName1.append(dictDatasets["observations"][obs][list_variables[0]]["areaname"])
                    try:
                        obsFileLandmask1.append(
                            dictDatasets["observations"][obs][list_variables[0]]["path + filename_landmask"])
                    except:
                        obsFileLandmask1.append(None)
                        obsLandmaskName1.append(None)
                    else:
                        obsLandmaskName1.append(dictDatasets["observations"][obs][list_variables[0]]["landmaskname"])
                    try:
                        obsInterpreter1.append(
                            dictDatasets["observations"][obs][list_variables[0]]["obs_interpreter"])
                    except:
                        obsInterpreter1.append(obs)
            # same if a second variable is needed
            # this time in the form of a keyarg dictionary
            arg_var2 = {
                "modelFileArea1": modelFileArea1, "modelAreaName1": modelAreaName1,
                "modelFileLandmask1": modelFileLandmask1, "modelLandmaskName1": modelLandmaskName1,
                "obsFileArea1": obsFileArea1, "obsAreaName1": obsAreaName1, "obsFileLandmask1": obsFileLandmask1,
                "obsLandmaskName1": obsLandmaskName1, "observed_fyear": observed_fyear,
                "observed_lyear": observed_lyear, "modeled_fyear": modeled_fyear, "modeled_lyear": modeled_lyear,
                "obsInterpreter1": obsInterpreter1}
            if len(list_variables) > 1:
                try:
                    arg_var2["modelFile2"] = dictDatasets["model"][modelName][list_variables[1]]["path + filename"]
                except:
                    arg_var2["modelFile2"] = ""
                try:
                    arg_var2["modelVarName2"] = dictDatasets["model"][modelName][list_variables[1]]["varname"]
                except:
                    arg_var2["modelVarName2"] = ""
                arg_var2["regionVar2"] = dict_regions[list_variables[1]]
                try:
                    arg_var2["modelFileArea2"] = \
                        dictDatasets["model"][modelName][list_variables[1]]["path + filename_area"]
                except:
                    arg_var2["modelFileArea2"], arg_var2["modelAreaName2"] = None, None
                else:
                    arg_var2["modelAreaName2"] = dictDatasets["model"][modelName][list_variables[1]]["areaname"]
                try:
                    arg_var2["modelFileLandmask2"] = \
                        dictDatasets["model"][modelName][list_variables[1]]["path + filename_landmask"]
                except:
                    arg_var2["modelFileLandmask2"], arg_var2["modelLandmaskName2"] = None, None
                else:
                    arg_var2["modelLandmaskName2"] = dictDatasets["model"][modelName][list_variables[1]]["landmaskname"]
                obsNameVar2, obsFile2, obsVarName2, obsFileArea2, obsAreaName2 = list(), list(), list(), list(), list()
                obsFileLandmask2, obsLandmaskName2, obsInterpreter2 = list(), list(), list()
                for obs in sorted(list(dictDatasets["observations"].keys()), key=lambda v: v.upper()):
                    try:
                        dictDatasets["observations"][obs][list_variables[1]]
                    except:
                        pass
                    else:
                        obsNameVar2.append(obs)
                        obsFile2.append(dictDatasets["observations"][obs][list_variables[1]]["path + filename"])
                        obsVarName2.append(dictDatasets["observations"][obs][list_variables[1]]["varname"])
                        try:
                            obsFileArea2.append(
                                dictDatasets["observations"][obs][list_variables[1]]["path + filename_area"])
                        except:
                            obsFileArea2.append(None)
                            obsAreaName2.append(None)
                        else:
                            obsAreaName2.append(dictDatasets["observations"][obs][list_variables[1]]["areaname"])
                        try:
                            obsFileLandmask2.append(
                                dictDatasets["observations"][obs][list_variables[1]]["path + filename_landmask"])
                        except:
                            obsFileLandmask2.append(None)
                            obsLandmaskName2.append(None)
                        else:
                            obsLandmaskName2.append(
                                dictDatasets["observations"][obs][list_variables[1]]["landmaskname"])
                        try:
                            obsInterpreter2.append(
                                dictDatasets["observations"][obs][list_variables[0]]["obs_interpreter"])
                        except:
                            obsInterpreter2.append(obs)
                arg_var2["obsNameVar2"] = obsNameVar2
                arg_var2["obsFile2"] = obsFile2
                arg_var2["obsVarName2"] = obsVarName2
                arg_var2["obsFileArea2"] = obsFileArea2
                arg_var2["obsAreaName2"] = obsAreaName2
                arg_var2["obsFileLandmask2"] = obsFileLandmask2
                arg_var2["obsLandmaskName2"] = obsLandmaskName2
                arg_var2['obsInterpreter2'] = obsInterpreter2
            # computes the metric
            if modelFile1 is None or len(modelFile1) == 0 or (isinstance(modelFile1, list) and None in modelFile1) or \
                    (len(list_variables) > 1 and
                     (arg_var2["modelFile2"] is None or len(arg_var2["modelFile2"]) == 0 or
                      (isinstance(arg_var2["modelFile2"], list) and None in arg_var2["modelFile2"]))):
                print("\033[94m" + str().ljust(5) + "ComputeCollection: " + str(metricCollection) + ", metric "
                      + str(metric) + " not computed" + "\033[0m")
                print("\033[94m" + str().ljust(10) + "reason(s):" + "\033[0m")
                if modelFile1 is None or len(modelFile1) == 0:
                    print("\033[94m" + str().ljust(11) + "no modeled " + list_variables[0] + " given" + "\033[0m")
                if isinstance(modelFile1, list) and None in modelFile1:
                    for ff, vv in zip(modelFile1, modelVarName1):
                        if ff is None or vv is None:
                            print("\033[94m" + str().ljust(11) + "no modeled " + str(vv) + " given" + "\033[0m")
                if (len(list_variables) > 1 and arg_var2["modelFile2"] is None) or \
                        (len(list_variables) > 1 and len(arg_var2["modelFile2"]) == 0):
                    print("\033[94m" + str().ljust(11) + "no modeled " + list_variables[1] + " given" + "\033[0m")
                if isinstance(arg_var2["modelFile2"], list) and None in arg_var2["modelFile2"]:
                    for ff, vv in zip(arg_var2["modelFile2"], arg_var2["modelVarName2"]):
                        if ff is None or vv is None:
                            print("\033[94m" + str().ljust(11) + "no modeled " + str(vv) + " given" + "\033[0m")
            elif obsFile1 is None or len(obsFile1) == 0 or (isinstance(obsFile1, list) and None in obsFile1) or \
                    (len(list_variables) > 1 and
                     (arg_var2["obsFile2"] is None or len(arg_var2["obsFile2"]) == 0 or
                      (isinstance(arg_var2["obsFile2"], list) and None in arg_var2["obsFile2"]))):
                print("\033[94m" + str().ljust(5) + "ComputeCollection: " + str(metricCollection) + ", metric "
                      + str(metric) + " not computed" + "\033[0m")
                print("\033[94m" + str().ljust(10) + "reason(s):" + "\033[0m")
                if obsFile1 is None or len(obsFile1) == 0:
                    print("\033[94m" + str().ljust(11) + "no observed " + list_variables[0] + " given" + "\033[0m")
                if isinstance(obsFile1, list) and None in obsFile1:
                    for ff, vv in zip(obsFile1, obsVarName1):
                        if ff is None or vv is None:
                            print("\033[94m" + str().ljust(11) + "no observed " + str(vv) + " given" + "\033[0m")
                if (len(list_variables) > 1 and arg_var2["obsFile2"] is None) or \
                        (len(list_variables) > 1 and len(arg_var2["obsFile2"]) == 0):
                    print("\033[94m" + str().ljust(11) + "no observed " + list_variables[1] + " given" + "\033[0m")
                if isinstance(arg_var2["obsFile2"], list) and None in arg_var2["obsFile2"]:
                    for ff, vv in zip(arg_var2["obsFile2"], arg_var2["obsVarName2"]):
                        if ff is None or vv is None:
                            print("\033[94m" + str().ljust(11) + "no observed " + str(vv) + " given" + "\033[0m")
            else:
                valu, vame, dive, dime = ComputeMetric(
                    metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                    dict_regions[list_variables[0]], user_regridding=user_regridding, debug=debug, netcdf=netcdf,
                    netcdf_name=netcdf_name, obs_interpreter=obs_interpreter, **arg_var2)
                keys1 = list(valu.keys())
                keys2 = list(set([kk.replace("value", "").replace("__", "").replace("_error", "")
                                  for ll in list(valu[keys1[0]].keys()) for kk in list(valu[keys1[0]][ll].keys())]))
                if len(keys2) > 1:
                    for kk in keys2:
                        mm, dd = dict(), dict()
                        keys3 = list(valu["metric"].keys())
                        for ll in keys3:
                            mm[ll] = {"value": valu["metric"][ll][kk + "__value"],
                                      "value_error": valu["metric"][ll][kk + "__value_error"]}
                        keys3 = list(valu["diagnostic"].keys())
                        for ll in keys3:
                            dd[ll] = {"value": valu["diagnostic"][ll][kk + "__value"],
                                      "value_error": valu["diagnostic"][ll][kk + "__value_error"]}
                        dict_col_valu[metric + kk] = {"metric": mm, "diagnostic": dd}
                        mm = dict((ii, vame["metric"][ii]) for ii in list(vame["metric"].keys()) if "units" not in ii)
                        mm["units"] = vame["metric"][kk + "__units"]
                        dict_col_meta["metrics"][metric + kk] = {"metric": mm, "diagnostic": vame["diagnostic"]}
                        dict_col_dd_valu[metric + kk], dict_col_dd_meta["metrics"][metric + kk] = dive, dime
                        del mm, dd
                else:
                    dict_col_valu[metric], dict_col_meta["metrics"][metric] = valu, vame
                    dict_col_dd_valu[metric], dict_col_dd_meta["metrics"][metric] = dive, dime
        except Exception as e:
            print(e)
            pass
    if dive_down is True:
        return {"value": dict_col_valu, "metadata": dict_col_meta}, \
               {"value": dict_col_dd_valu, "metadata": dict_col_dd_meta}
    else:
        return {"value": dict_col_valu, "metadata": dict_col_meta}, {}


def group_json_obs(pattern, json_name_out, metric_name):
    list_files = sorted(list(GLOBiglob(pattern)), key=lambda v: v.upper())
    for file1 in list_files:
        with open(file1) as ff:
            data = json.load(ff)
        ff.close()
        data = data["RESULTS"]["model"]
        for dataset in sorted(list(data.keys()), key=lambda v: v.upper()):
            try:    dict_out
            except: dict_out = deepcopy(data)
            else:
                if dataset in list(dict_out.keys()):
                    if metric_name == "all":
                        for met in sorted(list(data[dataset]["r1i1p1"]["value"].keys()), key=lambda v: v.upper()):
                            if met not in list(dict_out[dataset]["r1i1p1"]["value"].keys()):
                                dict_out[dataset]["r1i1p1"]["value"][met] = data[dataset]["r1i1p1"]["value"][met]
                            if met not in list(dict_out[dataset]["r1i1p1"]["metadata"]["metrics"].keys()):
                                dict_out[dataset]["r1i1p1"]["metadata"]["metrics"][met] = \
                                    data[dataset]["r1i1p1"]["metadata"]["metrics"][met]
                    else:
                        dict_out[dataset]["r1i1p1"]["value"][metric_name] = data[dataset]["r1i1p1"]["value"][
                            metric_name]
                else:
                    dict_out[dataset] = data[dataset]
        # OSremove(file1)
        del data, dataset, ff
    save_json_obs(dict_out, json_name_out)
    return


def save_json_obs(dict_in, json_name):
    dict_out = {"RESULTS": {"model": dict_in}}
    # save as json file
    with open(json_name + ".json", "w") as outfile:
        json.dump(dict_out, outfile, sort_keys=True)
    return


def ComputeCollection_ObsOnly(metricCollection, dictDatasets, user_regridding={}, debug=False, dive_down=False,
                              netcdf=False, netcdf_name="", observed_fyear=None, observed_lyear=None,
                              modeled_fyear=None, modeled_lyear=None, obs_interpreter=None):
    dict_mc = defCollection(metricCollection)
    dict_col_meta = {
        "name": dict_mc["long_name"], "description_of_the_collection": dict_mc["description"], "metrics": {},
    }
    dict_col_dd_meta = {
        "name": dict_mc["long_name"], "description_of_the_collection": dict_mc["description"], "metrics": {},
    }
    dict_col_valu = dict()
    dict_col_dd_valu = dict()
    dict_m = dict_mc["metrics_list"]
    list_metrics = sorted(list(dict_m.keys()), key=lambda v: v.upper())
    for metric in list_metrics:
        print("\033[94m" + str().ljust(5) + "ComputeCollection: metric = " + str(metric) + "\033[0m")
        # sets arguments for this metric
        list_variables = dict_m[metric]["variables"]
        dict_regions = dict_m[metric]["regions"]
        # observations name(s), file(s), variable(s) name in file(s)
        obsNameVar1, obsFile1, obsVarName1, obsFileArea1, obsAreaName1 = list(), list(), list(), list(), list()
        obsFileLandmask1, obsLandmaskName1, obsInterpreter1 = list(), list(), list()
        for obs in sorted(list(dictDatasets["observations"].keys()), key=lambda v: v.upper()):
            try:
                dictDatasets["observations"][obs][list_variables[0]]
            except:
                pass
            else:
                obsNameVar1.append(obs)
                obsFile1.append(dictDatasets["observations"][obs][list_variables[0]]["path + filename"])
                obsVarName1.append(dictDatasets["observations"][obs][list_variables[0]]["varname"])
                try:
                    obsFileArea1.append(dictDatasets["observations"][obs][list_variables[0]]["path + filename_area"])
                except:
                    obsFileArea1.append(None)
                    obsAreaName1.append(None)
                else:
                    obsAreaName1.append(dictDatasets["observations"][obs][list_variables[0]]["areaname"])
                try:
                    obsFileLandmask1.append(
                        dictDatasets["observations"][obs][list_variables[0]]["path + filename_landmask"])
                except:
                    obsFileLandmask1.append(None)
                    obsLandmaskName1.append(None)
                else:
                    obsLandmaskName1.append(dictDatasets["observations"][obs][list_variables[0]]["landmaskname"])
                try:
                    obsInterpreter1.append(dictDatasets["observations"][obs][list_variables[0]]["obs_interpreter"])
                except:
                    obsInterpreter1.append(obs)
        # same if a second variable is needed
        obsNameVar2, obsFile2, obsVarName2, obsFileArea2, obsAreaName2 = list(), list(), list(), list(), list()
        obsFileLandmask2, obsLandmaskName2, obsInterpreter2 = list(), list(), list()
        if len(list_variables) > 1:
            for obs in sorted(list(dictDatasets["observations"].keys()), key=lambda v: v.upper()):
                try:
                    dictDatasets["observations"][obs][list_variables[1]]
                except:
                    pass
                else:
                    obsNameVar2.append(obs)
                    obsFile2.append(dictDatasets["observations"][obs][list_variables[1]]["path + filename"])
                    obsVarName2.append(dictDatasets["observations"][obs][list_variables[1]]["varname"])
                    try:
                        obsFileArea2.append(
                            dictDatasets["observations"][obs][list_variables[1]]["path + filename_area"])
                    except:
                        obsFileArea2.append(None)
                        obsAreaName2.append(None)
                    else:
                        obsAreaName2.append(dictDatasets["observations"][obs][list_variables[1]]["areaname"])
                    try:
                        obsFileLandmask2.append(
                            dictDatasets["observations"][obs][list_variables[1]]["path + filename_landmask"])
                    except:
                        obsFileLandmask2.append(None)
                        obsLandmaskName2.append(None)
                    else:
                        obsLandmaskName2.append(
                            dictDatasets["observations"][obs][list_variables[1]]["landmaskname"])
                    try:
                        obsInterpreter2.append(
                            dictDatasets["observations"][obs][list_variables[0]]["obs_interpreter"])
                    except:
                        obsInterpreter2.append(obs)
        # observations as model
        print(obsNameVar1)
        print(obsNameVar2)
        for ii in range(len(obsFileArea1)):
            modelName = obsNameVar1[ii]
            modelFile1 = obsFile1[ii]
            modelVarName1 = obsVarName1[ii]
            modelFileArea1 = obsFileArea1[ii]
            modelAreaName1 = obsAreaName1[ii]
            modelFileLandmask1 = obsFileLandmask1[ii]
            modelLandmaskName1 = obsLandmaskName1[ii]
            modelInterpreter1 = obsInterpreter1[ii]
            arg_var2 = {
                "modelFileArea1": modelFileArea1, "modelAreaName1": modelAreaName1,
                "modelFileLandmask1": modelFileLandmask1, "modelLandmaskName1": modelLandmaskName1,
                "modelInterpreter1": modelInterpreter1, "obsFileArea1": obsFileArea1, "obsAreaName1": obsAreaName1,
                "obsFileLandmask1": obsFileLandmask1, "obsLandmaskName1": obsLandmaskName1,
                "obsInterpreter1": obsInterpreter1, "observed_fyear": observed_fyear, "observed_lyear": observed_lyear,
                "modeled_fyear": modeled_fyear, "modeled_lyear": modeled_lyear}
            if len(list_variables) == 1:
                nbr = 1
            else:
                nbr = len(obsNameVar2)
            for jj in range(nbr):
                try:
                    obsNameVar2[jj]
                except:
                    modelName2 = deepcopy(modelName)
                else:
                    modelName2 = modelName + "_" + obsNameVar2[jj]
                    arg_var2["modelFile2"] = obsFile2[jj]
                    arg_var2["modelVarName2"] = obsVarName2[jj]
                    arg_var2["modelFileArea2"] = obsFileArea2[jj]
                    arg_var2["modelAreaName2"] = obsAreaName2[jj]
                    arg_var2["modelFileLandmask2"] = obsFileLandmask2[jj]
                    arg_var2["modelLandmaskName2"] = obsLandmaskName2[jj]
                    arg_var2["modelInterpreter2"] = obsInterpreter2[jj]
                    arg_var2["regionVar2"] = dict_regions[list_variables[1]]
                    arg_var2["obsNameVar2"] = obsNameVar2
                    arg_var2["obsFile2"] = obsFile2
                    arg_var2["obsVarName2"] = obsVarName2
                    arg_var2["obsFileArea2"] = obsFileArea2
                    arg_var2["obsAreaName2"] = obsAreaName2
                    arg_var2["obsFileLandmask2"] = obsFileLandmask2
                    arg_var2["obsLandmaskName2"] = obsLandmaskName2
                    arg_var2["obsInterpreter2"] = obsInterpreter2
                if netcdf is True:
                    netcdf_name_out = netcdf_name.replace("OBSNAME", modelName2)
                else:
                    netcdf_name_out = ""
                json_name_out = netcdf_name.replace("OBSNAME", "tmp1_" + modelName2 + "_" + metric)
                if ("EnsoSstMap" in metric and modelName2 in sst_only) or \
                        len(list(GLOBiglob(json_name_out + ".json"))) == 1:
                    pass
                else:
                    print(modelName2 + "_as_model")
                    valu, vame, dive, dime = ComputeMetric(
                        metricCollection, metric, modelName2, modelFile1, modelVarName1, obsNameVar1,
                        obsFile1, obsVarName1, dict_regions[list_variables[0]], user_regridding=user_regridding,
                        debug=debug, netcdf=netcdf, netcdf_name=netcdf_name_out, obs_interpreter=obs_interpreter,
                        **arg_var2)
                    keys1 = list(valu.keys())
                    keys2 = list(set([kk.replace("value", "").replace("__", "").replace("_error", "")
                                      for ll in list(valu[keys1[0]].keys()) for kk in list(valu[keys1[0]][ll].keys())]))
                    if len(keys2) > 1:
                        for kk in keys2:
                            mm1, dd1 = dict(), dict()
                            keys3 = list(valu["metric"].keys())
                            for ll in keys3:
                                mm1[ll] = {"value": valu["metric"][ll][kk + "__value"],
                                           "value_error": valu["metric"][ll][kk + "__value_error"]}
                            keys3 = list(valu["diagnostic"].keys())
                            for ll in keys3:
                                dd1[ll] = {"value": valu["diagnostic"][ll][kk + "__value"],
                                           "value_error": valu["diagnostic"][ll][kk + "__value_error"]}
                            mm2 = dict((ll,
                                        vame["metric"][ll]) for ll in list(vame["metric"].keys()) if "units" not in ll)
                            mm2["units"] = vame["metric"][kk + "__units"]
                            dict1 = {"metric": mm1, "diagnostic": dd1}
                            dict2 = {"metric": mm2, "diagnostic": vame["diagnostic"]}
                            try:
                                dict_col_valu[modelName2]
                            except:
                                dict_col_valu[modelName2] = {metric + kk: dict1}
                                dict_col_meta[modelName2] = {"metrics": {metric + kk: dict2}}
                                dict_col_dd_valu[modelName2] = {metric + kk: dive}
                                dict_col_dd_meta[modelName2] = {"metrics": {metric + kk: dime}}
                            else:
                                dict_col_valu[modelName2][metric + kk] = dict1
                                dict_col_meta[modelName2]["metrics"][metric + kk] = dict2
                                dict_col_dd_valu[modelName2][metric + kk] = dive
                                dict_col_dd_meta[modelName2]["metrics"][metric + kk] = dime
                            del dd1, dict1, dict2, keys3, mm1, mm2
                    else:
                        try:
                            dict_col_valu[modelName2]
                        except:
                            dict_col_valu[modelName2] = {metric: valu}
                            dict_col_meta[modelName2] = {"metrics": {metric: vame}}
                            dict_col_dd_valu[modelName2] = {metric: dive}
                            dict_col_dd_meta[modelName2] = {"metrics": {metric: dime}}
                        else:
                            dict_col_valu[modelName2][metric] = valu
                            dict_col_meta[modelName2]["metrics"][metric] = vame
                            dict_col_dd_valu[modelName2][metric] = dive
                            dict_col_dd_meta[modelName2]["metrics"][metric] = dime
                    # save json
                    dict_out = {modelName2: {
                        "r1i1p1": {"value": dict_col_valu[modelName2], "metadata": dict_col_meta[modelName2]}}}
                    save_json_obs(dict_out, json_name_out)
                    del dict_out, dime, dive, keys1, keys2, valu
                del json_name_out, modelName2, netcdf_name_out
            del arg_var2, modelName, modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, \
                modelLandmaskName1, modelInterpreter1, nbr
        # read all jsons and group them
        pattern = netcdf_name.replace("OBSNAME", "tmp1_*_" + metric + ".json")
        json_name_out = netcdf_name.replace("OBSNAME", "tmp2_" + metric)
        group_json_obs(pattern, json_name_out, metric)
        del dict_regions, json_name_out, list_variables, obsAreaName1, obsAreaName2, obsFile1, obsFile2, obsFileArea1, \
            obsFileArea2, obsFileLandmask1, obsFileLandmask2, obsInterpreter1, obsInterpreter2, obsLandmaskName1,\
            obsLandmaskName2, obsNameVar1, obsNameVar2, obsVarName1, obsVarName2, pattern
    # read all jsons and group them
    pattern = netcdf_name.replace("OBSNAME", "tmp2_*.json")
    json_name_out = netcdf_name.replace("OBSNAME", "observation")
    group_json_obs(pattern, json_name_out, "all")
    return
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric
#
dict_oneVar_modelAndObs = {
    "BiasPrLatRmse": BiasPrLatRmse, "BiasPrLonRmse": BiasPrLonRmse, "BiasPrRmse": BiasPrMapRmse,
    "BiasLhfLatRmse": BiasLhfLatRmse, "BiasLhfLonRmse": BiasLhfLonRmse, "BiasLhfRmse": BiasLhfMapRmse,
    "BiasLwrLatRmse": BiasLwrLatRmse, "BiasLwrLonRmse": BiasLwrLonRmse, "BiasLwrRmse": BiasLwrMapRmse,
    "BiasShfLatRmse": BiasShfLatRmse, "BiasShfLonRmse": BiasShfLonRmse, "BiasShfRmse": BiasShfMapRmse,
    "BiasSshLatRmse": BiasSshLatRmse, "BiasSshLonRmse": BiasSshLonRmse, "BiasSshRmse": BiasSshMapRmse,
    "BiasSstLatRmse": BiasSstLatRmse, "BiasSstLonRmse": BiasSstLonRmse, "BiasSstRmse": BiasSstMapRmse,
    "BiasSwrLatRmse": BiasSwrLatRmse, "BiasSwrLonRmse": BiasSwrLonRmse, "BiasSwrRmse": BiasSwrMapRmse,
    "BiasTauxLatRmse": BiasTauxLatRmse, "BiasTauxLonRmse": BiasTauxLonRmse, "BiasTauxRmse": BiasTauxMapRmse,
    "BiasTauyLatRmse": BiasTauyLatRmse, "BiasTauyLonRmse": BiasTauyLonRmse, "BiasTauyRmse": BiasTauyMapRmse,
    "BiasThfLatRmse": BiasThfLatRmse, "BiasThfLonRmse": BiasThfLonRmse, "BiasThfRmse": BiasThfMapRmse,
    "BiasSstSkLonRmse": BiasSstSkLonRmse,
    "EnsoSstLonRmse": EnsoSstLonRmse, "NinaSstLonRmse": NinaSstLonRmse, "NinoSstTsRmse": NinoSstTsRmse,
    "EnsoSstMap": EnsoSstMap, "EnsoSstMapDjf": EnsoSstMapDjf, "EnsoSstMapJja": EnsoSstMapJja, "NinaSstMap": NinaSstMap,
    "NinoSstMap": NinoSstMap,
    "EnsoSstTsRmse": EnsoSstTsRmse, "NinaSstTsRmse": NinaSstTsRmse, "NinoSstLonRmse": NinoSstLonRmse,
    "NinaSstDivRmse": NinaSstDivRmse, "NinoSstDivRmse": NinoSstDivRmse,
    "SeasonalPrLatRmse": SeasonalPrLatRmse, "SeasonalPrLonRmse": SeasonalPrLonRmse,
    "SeasonalSshLatRmse": SeasonalSshLatRmse, "SeasonalSshLonRmse": SeasonalSshLonRmse,
    "SeasonalSstLatRmse": SeasonalSstLatRmse, "SeasonalSstLonRmse": SeasonalSstLonRmse,
    "SeasonalTauxLatRmse": SeasonalTauxLatRmse, "SeasonalTauxLonRmse": SeasonalTauxLonRmse,
    "SeasonalTauyLatRmse": SeasonalTauyLatRmse, "SeasonalTauyLonRmse": SeasonalTauyLonRmse,
}

dict_twoVar_modelAndObs = {
    "EnsoPrMap": EnsoPrMap, "EnsoPrMapDjf": EnsoPrMapDjf, "EnsoPrMapJja": EnsoPrMapJja,
    "EnsoPrDjfTel": EnsoPrDjfTel, "EnsoPrJjaTel": EnsoPrJjaTel,
    "EnsoSlpMap": EnsoSlpMap, "EnsoSlpMapDjf": EnsoSlpMapDjf, "EnsoSlpMapJja": EnsoSlpMapJja,
    "NinaPrMap": NinaPrMap, "NinaSlpMap": NinaSlpMap, "NinoPrMap": NinoPrMap, "NinoSlpMap": NinoSlpMap,
    "EnsoPrLonRmse": EnsoPrLonRmse, "EnsoPrTsRmse": EnsoPrTsRmse,
    "EnsoSshLonRmse": EnsoSshLonRmse, "EnsoSshTsRmse": EnsoSshTsRmse,
    "EnsoTauxLonRmse": EnsoTauxLonRmse, "EnsoTauxTsRmse": EnsoTauxTsRmse,
    "EnsoTauyLonRmse": EnsoTauyLonRmse, "EnsoTauyTsRmse": EnsoTauyTsRmse,
}

dict_oneVar = {
    "EnsoAmpl": EnsoAmpl, "EnsoDuration": EnsoDuration, "EnsoDiversity": EnsoDiversity,
    "EnsoSeasonality": EnsoSeasonality, "EnsoSstDiversity": EnsoSstDiversity, "EnsoSstSkew": EnsoSstSkew,
    "grad_lat_pr": grad_lat_pr, "grad_lat_sst": grad_lat_sst, "grad_lon_pr": grad_lon_pr, "grad_lon_sst": grad_lon_sst,
    "NinaSstDiv": NinaSstDiv, "NinaSstDur": NinaSstDur, "NinoSstDiv": NinoSstDiv, "NinoSstDiversity": NinoSstDiversity,
    "NinoSstDur": NinoSstDur,
}

dict_twoVar = {
    "EnsoFbSshSst": EnsoFbSshSst, "EnsoFbSstLhf": EnsoFbSstLhf, "EnsoFbSstLwr": EnsoFbSstLwr,
    "EnsoFbSstShf": EnsoFbSstShf, "EnsoFbSstSwr": EnsoFbSstSwr, "EnsoFbSstTaux": EnsoFbSstTaux,
    "EnsoFbSstThf": EnsoFbSstThf, "EnsoFbTauxSsh": EnsoFbTauxSsh, "EnsodSstOce": EnsodSstOce,
}


def ComputeMetric(metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                  regionVar1, modelFileArea1="", modelAreaName1="", modelFileLandmask1="", modelLandmaskName1="",
                  modelInterpreter1=None, obsFileArea1="", obsAreaName1="", obsFileLandmask1="", obsLandmaskName1="",
                  obsInterpreter1=None, modelFile2="", modelVarName2="", modelFileArea2="", modelAreaName2="",
                  modelFileLandmask2="", modelLandmaskName2="", modelInterpreter2=None, obsNameVar2="", obsFile2="",
                  obsVarName2="", obsFileArea2="", obsAreaName2="", obsFileLandmask2="", obsLandmaskName2="",
                  regionVar2="", obsInterpreter2=None, user_regridding={}, debug=False, netcdf=False, netcdf_name="",
                  observed_fyear=None, observed_lyear=None, modeled_fyear=None, modeled_lyear=None,
                  obs_interpreter=None):
    """
    :param metricCollection: string
        name of a Metric Collection, must be defined in EnsoCollectionsLib.defCollection()
    :param metric: string
        name of a Metric, must be defined in EnsoCollectionsLib.defCollection()
    :param modelName: string
        path_to/filename of a model
    :param modelFile1: string or list of strings
        model file number 1
        e.g.: 'model_file.nc' or ['model_file1.nc','model_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param modelVarName1: string or list of strings
        model variable number 1, must be defined in EnsoCollectionsLib.CmipVariables()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param obsNameVar1: string
        name of the observations for variable number 1, must be defined in EnsoCollectionsLib.ReferenceObservations()
    :param obsFile1: string or list of strings
        observations file number 1
        e.g.: 'observations_file.nc' or ['observations_file1.nc','observations_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param obsVarName1: string or list of strings
        observations variable number 1, must be defined in EnsoCollectionsLib.ReferenceObservations()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param regionVar1: string
        name of box ('tropical_pacific') for variable number 1, must be defined in EnsoCollectionsLib.ReferenceRegions()
    :param modelFileArea1: string, optional
        path_to/filename of the model areacell for variable number 1
    :param modelAreaName1: string, optional
        name the model areacell (e.g. 'areacella' or 'areacello') for variable number 1
    :param modelFileLandmask1: string, optional
        path_to/filename of the model landmask for variable number 1
    :param modelLandmaskName1: string, optional
        name the model landmask (e.g. 'landmask' or 'sftlf') for variable number 1
    :param modelInterpreter1: string, optional
        special variable interpreter for model variable number 1
        this is useful if the given model is not a CMIP model but an observational dataset
        the only possibility is 'CMIP' to interpret the variable number 1 as CMIP (dataset has been CMORized)
        default value = None, dataset will be interpreted as defined in EnsoCollectionsLib.ReferenceObservations if
        recognized as observations, else EnsoCollectionsLib.defCollection
    :param obsFileArea1: string, optional
        path_to/filename of the observations areacell for variable number 1
    :param obsAreaName1: string, optional
        name the observations areacell (e.g. 'areacella' or 'areacello') for variable number 1
    :param obsFileLandmask1: string, optional
        path_to/filename of the observations landmask for variable number 1
    :param obsLandmaskName1: string, optional
        name the observations landmask (e.g. 'landmask' or 'sftlf') for variable number 1
    :param obsInterpreter1: string, optional
        special variable interpreter for observational variable number 1
        the only possibility is 'CMIP' to interpret all observational's variable number 1 as CMIP (dataset has been
        CMORized)
        default value = None, observational dataset is considered not CMORized and will be interpreted as defined in
        EnsoCollectionsLib.ReferenceObservations
    :param modelFile2: string or list of strings, optional
        model file number 2
        e.g.: 'model_file.nc' or ['model_file1.nc','model_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param modelVarName2: string or list of strings, optional
        model variable number 2, must be defined in EnsoCollectionsLib.CmipVariables()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param modelFileArea2: string, optional
        path_to/filename of the model areacell for variable number 2
    :param modelAreaName2: string, optional
        name the model areacell (e.g. 'areacella' or 'areacello') for variable number 2
    :param modelFileLandmask2: string, optional
        path_to/filename of the model landmask for variable number 2
    :param modelLandmaskName2: string, optional
        name the model landmask (e.g. 'landmask' or 'sftlf') for variable number 2
    :param modelInterpreter2: string, optional
        special variable interpreter for model variable number 2
        this is useful if the given model is not a CMIP model but an observational dataset
        the only possibility is 'CMIP' to interpret the variable number 2 as CMIP (dataset has been CMORized)
        default value = None, dataset will be interpreted as defined in EnsoCollectionsLib.ReferenceObservations if
        recognized as observations, else EnsoCollectionsLib.defCollection
    :param obsNameVar2: string, optional
        name of the observations for variable number 2, must be defined in EnsoCollectionsLib.ReferenceObservations()
    :param obsFile2: string or list of strings, optional
        observations file number 2
        e.g.: 'observations_file.nc' or ['observations_file1.nc','observations_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param obsVarName2: string or list of strings, optional
        observations variable number 2, must be defined in EnsoCollectionsLib.ReferenceObservations()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param obsFileArea2: string, optional
        path_to/filename of the observations areacell for variable number 2
    :param obsAreaName2: string, optional
        name the observations areacell (e.g. 'areacella' or 'areacello') for variable number 2
    :param obsFileLandmask2: string, optional
        path_to/filename of the observations landmask for variable number 2
    :param obsLandmaskName2: string, optional
        name the observations landmask (e.g. 'landmask' or 'sftlf') for variable number 2
    :param regionVar2: string
        name of box ('tropical_pacific') for variable number 2, must be defined in EnsoCollectionsLib.ReferenceRegions()
    :param obsInterpreter2: string, optional
        special variable interpreter for observational variable number 2
        the only possibility is 'CMIP' to interpret all observational's variable number 2 as CMIP (dataset has been
        CMORized)
        default value = None, observational dataset is considered not CMORized and will be interpreted as defined in
        EnsoCollectionsLib.ReferenceObservations
    :param user_regridding: dict, optional
        regridding parameters selected by the user and not the "climate experts" (i.e. the program)
        to see parameters:
        help(EnsoUvcdatToolsLib.TwoVarRegrid)
        help(EnsoUvcdatToolsLib.Regrid)
        e.g.:
        user_regridding = {
            'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                'newgrid_name': 'generic 1x1deg'},
        }
    :param debug: boolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' root name of the saved NetCDFs
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='USER_DATE_METRICCOLLECTION_MODEL'
    :param observed_fyear: integer, optional
        first year to use for observational datasets, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param observed_lyear: integer, optional
        last year to use for observational datasets, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param modeled_fyear: integer, optional
        first year to use for CMIP simulations, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'observed_period' defined in EnsoCollectionsLib.defCollection is used
    :param modeled_lyear: integer, optional
        last year to use for CMIP simulations, given to overrule value defined in EnsoCollectionsLib.py
        default value = None, 'modeled_period' defined in EnsoCollectionsLib.defCollection is used
    :param obs_interpreter: string, optional
        special variable interpreter for all observational datasets
        the only possibility is 'CMIP' to interpret all observational's variables as CMIP (datasets have been CMORized)
        default value = None, observational datasets are considered not CMORized and will be interpreted as defined in
        EnsoCollectionsLib.ReferenceObservations

    :return:
    """
    tmp_metric = deepcopy(metric)
    metric = metric.replace("_1", "").replace("_2", "").replace("_3", "").replace("_4", "").replace("_5", "")
    # retrieving keyargs from EnsoCollectionsLib.defCollection
    dict_mc = defCollection(metricCollection)
    list_obs = sorted(list(ReferenceObservations().keys()), key=lambda v: v.upper())

    # common_collection_parameters
    keyarg = dict()
    for arg in list(dict_mc["common_collection_parameters"].keys()):
        keyarg[arg] = dict_mc["common_collection_parameters"][arg]
    for arg in list(dict_mc["metrics_list"][tmp_metric].keys()):
        keyarg[arg] = dict_mc["metrics_list"][tmp_metric][arg]
    # if "metric_computation" is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its
    # default value
    try:
        keyarg["metric_computation"]
    except:
        keyarg["metric_computation"] = default_arg_values("metric_computation")
    # if "modeled_period" is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg["time_bounds_mod"] = keyarg["modeled_period"]
    except:
        keyarg["time_bounds_mod"] = default_arg_values("time_bounds_mod")
    # YYP !!! experimental period defined bt user !!!
    if isinstance(modeled_fyear, int) is True and isinstance(modeled_lyear, int) is True:
        keyarg["time_bounds_mod"] = (str(modeled_fyear) + "-01-01 00:00:00", str(modeled_lyear) + "-12-31 23:59:60.0")
    # if "modeled_period" is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg["time_bounds_obs"] = keyarg["observed_period"]
    except:
        keyarg["time_bounds_obs"] = default_arg_values("time_bounds_obs")
    # YYP !!! experimental period defined bt user !!!
    if isinstance(observed_fyear, int) is True and isinstance(observed_lyear, int) is True:
        keyarg["time_bounds_obs"] = (str(observed_fyear) + "-01-01 00:00:00", str(observed_lyear) + "-12-31 23:59:60.0")
    # if the user gave a specific regridding Tool / method, use it
    if tmp_metric in list(user_regridding.keys()):
        keyarg["regridding"] = user_regridding[tmp_metric]
    elif "regridding" in list(user_regridding.keys()):
        keyarg["regridding"] = user_regridding["regridding"]

    # if model is an observation
    if modelName.split("_")[0] in list_obs and obs_interpreter != "CMIP":
        keyarg["project_interpreter_mod_var1"] = \
            modelName.split("_")[0] if modelInterpreter1 is None else deepcopy(modelInterpreter1)
        if modelFile2 != "":
            keyarg["project_interpreter_mod_var2"] = \
                modelName.split("_")[1] if modelInterpreter2 is None else deepcopy(modelInterpreter2)
    else:
        keyarg["project_interpreter_mod_var1"] = keyarg["project_interpreter"]
        if modelFile2 != "":
            keyarg["project_interpreter_mod_var2"] = keyarg["project_interpreter"]
        keyarg["time_bounds_mod"] = deepcopy(keyarg["time_bounds_obs"])

    # obsName could be a list if the user wants to compare the model with a set of observations
    # if obsName is just a name (string) it is put in a list
    if isinstance(obsNameVar1, str) is True or isinstance(obsNameVar1, str) is True:
        obsNameVar1 = [obsNameVar1]
    if isinstance(obsFile1, str) is True or isinstance(obsFile1, str) is True:
        obsFile1 = [obsFile1]
    if isinstance(obsVarName1, str) is True or isinstance(obsVarName1, str) is True:
        obsVarName1 = [obsVarName1]
    # Var2, do the same as for Var1
    if isinstance(obsNameVar2, str) is True or isinstance(obsNameVar2, str) is True:
        obsNameVar2 = [obsNameVar2]
    if isinstance(obsFile2, str) is True or isinstance(obsFile2, str) is True:
        obsFile2 = [obsFile2]
    if isinstance(obsVarName2, str) is True or isinstance(obsVarName2, str) is True:
        obsVarName2 = [obsVarName2]

    dict_metric_val = dict()
    dict_diagnostic = dict()
    dict_diagnostic_metadata = dict()
    dict_dive_down = dict()
    dict_dive_down_metadata = dict()

    multimetric = False

    # test files
    if isinstance(modelFile1, list):
        noerror = all(isinstance(file1, str) is True for file1 in modelFile1)
    else:
        noerror = True if isinstance(modelFile1, str) is True else False
    if noerror is False:
        if isinstance(modelFile1, list):
            tmp = str([modelVarName1[kk] for kk, file1 in enumerate(modelFile1) if isinstance(file1, str) is False])
        else:
            tmp = str(modelVarName1)
        tmperr = "model var1 (" + tmp + ") not given"
    if metric in list(dict_twoVar.keys()) + list(dict_twoVar_modelAndObs.keys()):
        if isinstance(modelFile2, list):
            noerror2 = all(isinstance(file1, str) is True for file1 in modelFile2)
        else:
            noerror2 = True if isinstance(modelFile2, str) is True else False
        if noerror2 is False:
            if isinstance(modelFile2, list):
                tmp = str([modelVarName2[kk] for kk, file1 in enumerate(modelFile2) if isinstance(file1, str) is False])
            else:
                tmp = str(modelVarName2)
            try:
                tmperr
            except:
                tmperr = "model var2(" + tmp + ") not given"
            else:
                tmperr = tmperr + " ; model var2(" + tmp + ") not given"
        noerror = False if noerror is False or noerror2 is False else True
    if noerror is False:
        dict_metrics = {
            "metric": {"value": None, "value_error": None}, "diagnostic": {"value": None, "value_error": None}}
        dict_metadata = {
            "metric": {"name": metric, "method": None, "datasets": modelName, "units": None},
            "diagnostic": {{modelName: {"method": None, "name": None, "ref": None, "time_frequency": None,
                                        "units": None, "keyerror": tmperr}}},
        }
        dict_dive_down, dict_dive_down_metadata = {}, {}
    else:
        if metric in list(dict_oneVar_modelAndObs.keys()) or metric in list(dict_twoVar_modelAndObs.keys()):
            #
            # this part regroups all diagnostics comparing model and obs (rmse)
            # so the diagnostic is the metric
            #
            description_metric = "The metric is the statistical value between the model and the observations"
            diagnostic1 = dict()
            for ii in range(len(obsNameVar1)):
                # computes the diagnostic/metric
                if metric in list(dict_oneVar_modelAndObs.keys()):
                    output_name = obsNameVar1[ii]
                    if output_name != modelName:
                        if "EnsoSstMap" in metric and output_name in sst_only:
                            pass
                        else:
                            print("\033[94m" + str().ljust(5) + "ComputeMetric: oneVarRMSmetric, " + metric + " = "
                                  + modelName + " and " + output_name + "\033[0m")
                            keyarg["project_interpreter_mod"] = deepcopy(keyarg["project_interpreter_mod_var1"])
                            keyarg["project_interpreter_obs"] = deepcopy(obsInterpreter1[ii])
                            diagnostic1[output_name] = dict_oneVar_modelAndObs[metric](
                                modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1,
                                modelLandmaskName1, obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii],
                                obsFileLandmask1[ii], obsLandmaskName1[ii], regionVar1, dataset1=modelName,
                                dataset2=output_name, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                                metname=tmp_metric, **keyarg)
                            del keyarg["project_interpreter_mod"]
                            del keyarg["project_interpreter_obs"]
                elif metric in list(dict_twoVar_modelAndObs.keys()):
                    for jj in range(len(obsNameVar2)):
                        output_name = obsNameVar1[ii] + "_" + obsNameVar2[jj]
                        if output_name != modelName:
                            print("\033[94m" + str().ljust(5) + "ComputeMetric: twoVarRMSmetric, " + metric + " = " +
                                  modelName + " and " + output_name + "\033[0m")
                            diagnostic1[output_name] = dict_twoVar_modelAndObs[metric](
                                modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1,
                                modelLandmaskName1, modelFile2, modelVarName2, modelFileArea2, modelAreaName2,
                                modelFileLandmask2, modelLandmaskName2, obsFile1[ii], obsVarName1[ii], obsFileArea1[ii],
                                obsAreaName1[ii], obsFileLandmask1[ii], obsLandmaskName1[ii], obsFile2[jj],
                                obsVarName2[jj], obsFileArea2[jj], obsAreaName2[jj], obsFileLandmask2[jj],
                                obsLandmaskName2[jj], regionVar1, regionVar2, dataset1=modelName, dataset2=output_name,
                                debug=debug, netcdf=netcdf, netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
            for obs in list(diagnostic1.keys()):
                # puts metric values in its proper dictionary
                if "value" in list(diagnostic1[obs].keys()):
                    dict_metric_val[obs] = {
                        "value": diagnostic1[obs]["value"], "value_error": diagnostic1[obs]["value_error"]}
                    dict_diagnostic[modelName] = {"value": None, "value_error": None}
                    dict_diagnostic[obs] = {"value": None, "value_error": None}
                else:
                    multimetric = True
                    lkeys = list(set([key.split("__")[0] for key in list(diagnostic1[obs].keys()) if "value" in key]))
                    for key in lkeys:
                        try:
                            dict_metric_val[obs]
                        except:
                            dict_metric_val[obs] = {key + "__value": diagnostic1[obs][key + "__value"],
                                                    key + "__value_error": diagnostic1[obs][key + "__value_error"]}
                            dict_diagnostic[modelName] = {key + "__value": None, key + "__value_error": None}
                            dict_diagnostic[obs] = {key + "__value": None, key + "__value_error": None}
                        else:
                            dict_metric_val[obs][key + "__value"] = diagnostic1[obs][key + "__value"]
                            dict_metric_val[obs][key + "__value_error"] = diagnostic1[obs][key + "__value_error"]
                            dict_diagnostic[modelName][key + "__value"] = None
                            dict_diagnostic[modelName][key + "__value_error"] = None
                            dict_diagnostic[obs][key + "__value"] = None
                            dict_diagnostic[obs][key + "__value_error"] = None
                if "dive_down_diag" in list(diagnostic1[obs].keys()):
                    dict_dive_down[modelName] = diagnostic1[obs]["dive_down_diag"]["model"]
                    dict_dive_down[obs] = diagnostic1[obs]["dive_down_diag"]["observations"]
                    dict1 = {}
                    for elt in list(diagnostic1[obs]["dive_down_diag"].keys()):
                        if elt not in ["model", "observations"]:
                            dict1[elt] = diagnostic1[obs]["dive_down_diag"][elt]
                    dict_dive_down_metadata[obs] = dict1
                    del dict1
                # puts diagnostic metadata in its proper dictionary
                dict_diagnostic_metadata[modelName] = {"name": modelName, "nyears": diagnostic1[obs]["nyears_model"],
                                                       "time_period": diagnostic1[obs]["time_period_model"]}
                dict_diagnostic_metadata[obs] = {"name": obs, "nyears": diagnostic1[obs]["nyears_observations"],
                                                 "time_period": diagnostic1[obs]["time_period_observations"]}
                if "events_model" in list(diagnostic1[obs].keys()):
                    dict_diagnostic_metadata[modelName]["events"] = diagnostic1[obs]["events_model"]
                    dict_diagnostic_metadata[obs]["events"] = diagnostic1[obs]["events_observations"]
                if "keyerror" in list(diagnostic1[obs].keys()):
                    dict_diagnostic_metadata[modelName]["keyerror"] = diagnostic1[obs]["keyerror"]
            units = diagnostic1[obs]["units"]
            diagnostic1 = diagnostic1[obs]
        else:
            #
            # model diagnostic
            #
            keyarg["time_bounds"] = keyarg["time_bounds_mod"]
            keyarg["project_interpreter_var1"] = keyarg["project_interpreter_mod_var1"]
            if metric in list(dict_oneVar.keys()):
                # computes diagnostic that needs only one variable
                print("\033[94m" + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(modelName) + "\033[0m")
                diagnostic1 = dict_oneVar[metric](
                    modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                    regionVar1, dataset=modelName, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                    metname=tmp_metric, **keyarg)
            elif metric in list(dict_twoVar.keys()):
                # computes diagnostic that needs two variables
                print("\033[94m" + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(modelName) + "\033[0m")
                keyarg["project_interpreter_var2"] = keyarg["project_interpreter_mod_var2"]
                diagnostic1 = dict_twoVar[metric](
                    modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                    regionVar1, modelFile2, modelVarName2, modelFileArea2, modelAreaName2, modelFileLandmask2,
                    modelLandmaskName2, regionVar2, dataset=modelName, debug=debug, netcdf=netcdf,
                    netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
            else:
                diagnostic1 = None
                list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": metric",
                                str().ljust(5) + "unknown metric name: " + str(metric)]
                EnsoErrorsWarnings.my_error(list_strings)
            # puts metric / diagnostic values in its proper dictionary
            dict_diagnostic[modelName] = {"value": diagnostic1["value"], "value_error": diagnostic1["value_error"]}
            if "nonlinearity" in list(diagnostic1.keys()):
                dict_diagnostic[modelName]["nonlinearity"] = diagnostic1["nonlinearity"]
                dict_diagnostic[modelName]["nonlinearity_error"] = diagnostic1["nonlinearity_error"]
            if "dive_down_diag" in list(diagnostic1.keys()):
                dict_dive_down[modelName] = diagnostic1["dive_down_diag"]["value"]
                for elt in list(diagnostic1["dive_down_diag"].keys()):
                    if elt not in ["value"]:
                        try:
                            dict_dive_down_metadata[modelName]
                        except:
                            dict_dive_down_metadata[modelName] = {elt: diagnostic1["dive_down_diag"][elt]}
                        else:
                            dict_dive_down_metadata[modelName][elt] = diagnostic1["dive_down_diag"][elt]
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata[modelName] = {
                "name": modelName, "nyears": diagnostic1["nyears"], "time_period": diagnostic1["time_period"],
            }
            if "events_model" in list(diagnostic1.keys()):
                dict_diagnostic_metadata[modelName]["events"] = diagnostic1["events"]
            if "keyerror" in list(diagnostic1.keys()):
                dict_diagnostic_metadata[modelName]["keyerror"] = diagnostic1["keyerror"]
            #
            # observations diag
            #
            diag_obs = dict()
            keyarg["time_bounds"] = keyarg["time_bounds_obs"]
            for ii in range(len(obsNameVar1)):
                # sets observations
                obs1 = obsNameVar1[ii]
                keyarg["project_interpreter_var1"] = deepcopy(obsInterpreter1[ii])
                if metric in list(dict_oneVar.keys()):
                    output_name = obs1
                    if output_name != modelName:
                        print("\033[94m" + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(output_name) +
                              "\033[0m")
                        diag_obs[output_name] = dict_oneVar[metric](
                            obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                            obsLandmaskName1[ii], regionVar1, dataset=output_name, debug=debug, netcdf=netcdf,
                            netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
                elif metric in list(dict_twoVar.keys()):
                    for jj in range(len(obsNameVar2)):
                        obs2 = obsNameVar2[jj]
                        output_name = obs1 + "_" + obs2
                        keyarg["project_interpreter_var2"] = deepcopy(obsInterpreter2[jj])
                        if output_name != modelName:
                            print("\033[94m" + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(output_name) +
                                  "\033[0m")
                            diag_obs[output_name] = dict_twoVar[metric](
                                obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                                obsLandmaskName1[ii], regionVar1, obsFile2[jj], obsVarName2[jj], obsFileArea2[jj],
                                obsAreaName2[jj], obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar2,
                                dataset=output_name, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                                metname=tmp_metric, **keyarg)
            for obs in list(diag_obs.keys()):
                # computes the metric
                metric_val, metric_err, description_metric = math_metric_computation(
                    diagnostic1["value"], diagnostic1["value_error"], obs=diag_obs[obs]["value"],
                    obs_err=diag_obs[obs]["value_error"], keyword=keyarg["metric_computation"])
                dict_metric_val[obs] = {"value": metric_val, "value_error": metric_err}
                # puts metric / diagnostic values in its proper dictionary
                dict_diagnostic[obs] = {"value": diag_obs[obs]["value"], "value_error": diag_obs[obs]["value_error"]}
                if "nonlinearity" in list(diag_obs[obs].keys()):
                    dict_diagnostic[obs]["nonlinearity"] = diag_obs[obs]["nonlinearity"]
                    dict_diagnostic[obs]["nonlinearity_error"] = diag_obs[obs]["nonlinearity_error"]
                if "dive_down_diag" in list(diag_obs[obs].keys()):
                    dict_dive_down[obs] = diag_obs[obs]["dive_down_diag"]["value"]
                    for elt in list(diag_obs[obs]["dive_down_diag"].keys()):
                        if elt not in ["value"]:
                            try:
                                dict_dive_down_metadata[obs]
                            except:
                                dict_dive_down_metadata[obs] = {elt: diag_obs[obs]["dive_down_diag"][elt]}
                            else:
                                dict_dive_down_metadata[obs][elt] = diag_obs[obs]["dive_down_diag"][elt]
                # puts diagnostic metadata in its proper dictionary
                dict_diagnostic_metadata[obs] = {
                    "name": obs, "nyears": diag_obs[obs]["nyears"], "time_period": diag_obs[obs]["time_period"],
                }
                if "events_model" in list(diag_obs[obs].keys()):
                    dict_diagnostic_metadata[obs]["events"] = diag_obs[obs]["events"]
                if "keyerror" in list(diag_obs[obs].keys()):
                    dict_diagnostic_metadata[obs]["keyerror"] = diag_obs[obs]["keyerror"]
            if keyarg["metric_computation"] in ["ratio", "relative_difference"]:
                units = diagnostic1["units"] + " / " + diagnostic1["units"]
            elif keyarg["metric_computation"] in ["abs_relative_difference"]:
                units = "%"
            else:
                units = diagnostic1["units"]
            try: del keyarg["time_bounds"]
            except: pass
        # finishes to fill the diagnostic dictionary
        list_keys = ["method", "name", "ref", "time_frequency", "units"]
        for key in list_keys:
            if key == "method":
                dict_diagnostic_metadata[key] = diagnostic1["method"]
                dict_dive_down_metadata[key] = diagnostic1["method"]
                try:
                    diagnostic1["nonlinearity"]
                except:
                    pass
                else:
                    dict_diagnostic_metadata["method_nonlinearity"] = diagnostic1["method_nonlinearity"]
            else:
                dict_diagnostic_metadata[key] = diagnostic1[key]
                dict_dive_down_metadata[key] = diagnostic1[key]
        # creates the output dictionaries
        dict_metrics = {"metric": dict_metric_val, "diagnostic": dict_diagnostic}
        datasets = modelName + "; "
        for ii in range(len(obsNameVar1)):
            if isinstance(obsVarName1[ii], list):
                datasets = datasets + obsNameVar1[ii] + "'s "
                for jj in range(len(obsVarName1[ii])):
                    datasets += obsVarName1[ii][jj]
                    if jj != len(obsVarName1[ii]) - 1:
                        datasets += " & "
            else:
                datasets = datasets + obsNameVar1[ii] + "'s " + obsVarName1[ii]
            if ii != len(obsNameVar1) - 1:
                datasets = datasets + "; "
        if len(obsNameVar2) > 0:
            datasets = datasets + "; "
            for ii in range(len(obsNameVar2)):
                if isinstance(obsVarName2[ii], list):
                    datasets = datasets + obsNameVar2[ii] + "'s "
                    for jj in range(len(obsVarName2[ii])):
                        datasets += obsVarName2[ii][jj]
                        if jj != len(obsVarName2[ii]) - 1:
                            datasets += " & "
                else:
                    datasets = datasets + obsNameVar2[ii] + "'s " + obsVarName2[ii]
                if ii != len(obsNameVar2) - 1:
                    datasets = datasets + "; "
        if multimetric is True:
            tmp = {"name": metric, "method": description_metric, "datasets": datasets}
            for key in lkeys:
                tmp[key + "__units"] = diagnostic1[key + "__units"]
            dict_metadata = {"metric": tmp, "diagnostic": dict_diagnostic_metadata}
        else:
            dict_metadata = {
                "metric": {"name": metric, "method": description_metric, "datasets": datasets, "units": units},
                "diagnostic": dict_diagnostic_metadata}
    return dict_metrics, dict_metadata, dict_dive_down, dict_dive_down_metadata
# ---------------------------------------------------------------------------------------------------------------------#
