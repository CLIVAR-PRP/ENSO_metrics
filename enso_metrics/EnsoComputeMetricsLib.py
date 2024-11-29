# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions to compute collections or recipes
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
from inspect import getmembers, isfunction
from typing import Union

# local functions
from enso_metrics import recipes
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
all_recipe = dict((k[0], k[1]) for k in getmembers(recipes, isfunction))


def compute_metric(
        metric: str,
        metric_param: dict,
        project: str,
        model: str,
        experiment: str,
        member: str,
        dict_mod: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, None]]]]],
        dict_obs: dict[
            str, dict[
                str, dict[
                    str, Union[int, float, str, list[str], None, dict[
                        str, Union[int, float, None]]]]]],
        **kwargs):
    #
    # -- compute model and observations one after the others
    #
    # put model and obs in a single dictionary
    dict_a = {model: dict_mod}
    dict_a.update(dict_obs)
    # loop on datasets
    for dat, dic in dict_a.items():
        if metric not in list(all_recipe.keys()):
            continue
        # update metric_param
        metric_p = deepcopy(metric_param)
        metric_p["time_bounds"] = metric_param["time_bounds_mod"] if dat == model else metric_param["time_bounds_obs"]
        diagnostic = all_recipe[metric](
                    dic, dataset=modelName, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                    metname=tmp_metric, **keyarg)
            else:
                diagnostic1 = None
                list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": metric",
                                str().ljust(5) + "unknown metric name: " + str(metric)]
                EnsoErrorsWarnings.my_error(list_strings)
            #
            # observations diag
            #
            diag_obs = dict()
            keyarg["time_bounds"] = keyarg["time_bounds_obs"]
            for ii in range(len(obsNameVar1)):
                # sets observations
                keyarg["project_interpreter_var1"] = \
                    "CMIP" if obs_interpreter == "CMIP" else deepcopy(obsInterpreter1[ii])
                if metric in list(dict_oneVar.keys()):
                    output_name = deepcopy(obsNameVar1[ii])
                    if output_name != modelName:
                        print("\033[94m" + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(output_name) +
                              "\033[0m")
                        diag_obs[output_name] = dict_oneVar[metric](
                            obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                            obsLandmaskName1[ii], regionVar1, dataset=output_name, debug=debug, netcdf=netcdf,
                            netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
                    del output_name
                elif metric in list(dict_twoVar.keys()):
                    for jj in range(len(obsNameVar2)):
                        output_name = deepcopy(obsNameVar1[ii]) + "_" + deepcopy(obsNameVar2[jj])
                        keyarg["project_interpreter_var2"] = \
                            "CMIP" if obs_interpreter == "CMIP" else deepcopy(obsInterpreter2[jj])
                        # if obsNameVar1[ii] == "Tropflux":
                        # if ("AVISO" in obsNameVar2[jj] and obsNameVar1[ii] == "Tropflux") or \
                        #         obsNameVar2[jj] == obsNameVar1[ii]:
                        if output_name != modelName:
                            print("\033[94m" + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(output_name) +
                                  "\033[0m")
                            diag_obs[output_name] = dict_twoVar[metric](
                                obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                                obsLandmaskName1[ii], regionVar1, obsFile2[jj], obsVarName2[jj], obsFileArea2[jj],
                                obsAreaName2[jj], obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar2,
                                dataset=output_name, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                                metname=tmp_metric, **keyarg)
                        del output_name
                        del keyarg["project_interpreter_var2"]
                del keyarg["project_interpreter_var1"]
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
                "diagnostic": dict_diagnostic_metadata,
            }
    return dict_metrics, dict_metadata, dict_dive_down, dict_dive_down_metadata
# ---------------------------------------------------------------------------------------------------------------------#
