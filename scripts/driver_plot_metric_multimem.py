# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create divedown plots for ENSO_metrics
# netCDF files are not provided in the package to do these plots
# Updated json files (needed to create this plot) can be downloaded from the page "Summary statistics in Interactive
# Portrait Plots" at https://cmec.llnl.gov/results/enso/
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoPlots.EnsoMetricPlot import plotter_multimem_1obs

# ---------------------------------------------------#


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "telecon_mix"
project = "CMIP6"
model = "MIROC6"  # "CNRM-CM6-1"  # "NorCPM1"  #
experiment = "historical"
# computation version, 'v20200427' for models and 'v20201231' for obs are provided with the package
version_mod = "v20200427"
version_obs = "v20201231"
# json files
dict_json = {
    "CMIP5": {
        "ENSO_perf": "share/EnsoMetrics/cmip5_historical_ENSO_perf_" + version_mod + "_allModels_allRuns.json",
        "ENSO_proc": "share/EnsoMetrics/cmip5_historical_ENSO_proc_" + version_mod + "_allModels_allRuns.json",
        "ENSO_tel": "share/EnsoMetrics/cmip5_historical_ENSO_tel_" + version_mod + "_allModels_allRuns.json"},
    "CMIP6": {
        "ENSO_perf": "share/EnsoMetrics/cmip6_historical_ENSO_perf_" + version_mod + "_allModels_allRuns.json",
        "ENSO_proc": "share/EnsoMetrics/cmip6_historical_ENSO_proc_" + version_mod + "_allModels_allRuns.json",
        "ENSO_tel": "share/EnsoMetrics/cmip6_historical_ENSO_tel_" + version_mod + "_allModels_allRuns.json"},
    "obs2obs": {
        "ENSO_perf": "share/EnsoMetrics/obs2obs_historical_ENSO_perf_" + version_obs + "_allObservations.json",
        "ENSO_proc": "share/EnsoMetrics/obs2obs_historical_ENSO_proc_" + version_obs + "_allObservations.json",
        "ENSO_tel": "share/EnsoMetrics/obs2obs_historical_ENSO_tel_" + version_obs + "_allObservations.json"}}
# path
# path = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/v20210105_metrics_macbook"
path_j = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Data"
path_n = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Data/" + \
    str(project).lower() + "/historical/" + str(metric_collection)
path_p = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Plots"
# figure pattern
figure_pattern = str(project.lower()) + "_" + str(experiment) + "_" + str(metric_collection) + "_METRIC_" + \
    str(model) + "_NBR_MEM"
# json pattern
json_pattern = str(metric_collection) + "_v20210611_allModels_allRuns.json"
# netcdf pattern
netcdf_pattern = str(metric_collection) + "_" + str(project) + "_" + str(model) + "_" + str(experiment) + "_MEM" + \
    "_METRIC.nc"
# list metrics
list_metrics = ["telecon_pr_ano_djf", "telecon_pr_amp_djf"]
# list observations
list_observations = ["OISSTv2_CMAP", "OISSTv2_GPCPv2.3"]
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
for met in list_metrics:
    for obs in list_observations:
        # observation name without dots
        obs2 = obs.replace(".", "")
        # get json file
        json_file = list(GLOBiglob(OSpath__join(path_j, json_pattern)))[0]
        with open(json_file) as ff:
            data_json = json.load(ff)["RESULTS"]["model"][model]
        ff.close()
        # list members
        list_members = sorted(list(data_json.keys()), key=lambda v: v.upper())
        # get NetCDF files
        netcdf_file = OSpath__join(path_n, netcdf_pattern.replace("METRIC", met))
        netcdf_file = dict(
            (str(model) + "_" + str(k1), list(GLOBiglob(netcdf_file.replace("MEM", k1)))[0]) for k1 in list_members)
        # get diagnostic values for the given model and observations
        diagnostic_values = dict(
            (k2, data_json[k1]["value"][met]["diagnostic"][k2]["value"])
            for k1 in list_members for k2 in list(data_json[k1]["value"][met]["diagnostic"].keys())
            if k2 == model or k2 == str(model) + "_" + str(k1))
        diagnostic_values[obs] = data_json[list_members[0]]["value"][met]["diagnostic"][obs]["value"]
        diagnostic_units = data_json[list_members[0]]["metadata"][met]["diagnostic"]["units"]
        # get metric values computed with the given model and observations
        metric_values = dict(
            (str(model) + "_" + str(k1), data_json[k1]["value"][met]["metric"][obs]["value"]) for k1 in list_members)
        metric_units = data_json[list_members[0]]["metadata"][met]["metric"]["units"]
        # figure name
        name_png = OSpath__join(path_p, figure_pattern.replace("METRIC", met))
        name_png = name_png.replace("NBR_MEM", str(len(list_members)) + "mem") + "_" + str(obs2)
        # plot
        plotter_multimem_1obs(
            metric_collection, met, project, model, experiment, netcdf_file, diagnostic_values, diagnostic_units,
            metric_values, metric_units, name_png, reference=obs, member=list_members)
        # delete
        del diagnostic_units, diagnostic_values, ff, json_file, list_members, metric_values, metric_units, name_png, \
            netcdf_file, obs2
