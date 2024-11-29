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
from EnsoPlots.EnsoMetricPlot import plotter_1mod_1obs
# ---------------------------------------------------#


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "telecon_mix"
project = "CMIP6"
model = "NorCPM1"
experiment = "historical"
member = "r1i1p1f1"
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
        "ENSO_tel":  "share/EnsoMetrics/obs2obs_historical_ENSO_tel_" + version_obs + "_allObservations.json"}}
# path
# path = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/v20210105_metrics_macbook"
path_j = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Data"
path_n = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Data/" + \
         str(project).lower() + "/historical/" + str(metric_collection)
path_p = "/Users/yannplanton/Documents/Yann/Fac/2019_2021_postdoc_NOAA/2021_06_11_teleconnections/Plots"
# figure pattern
figure_pattern = str(project.lower()) + "_" + str(experiment) + "_" + str(metric_collection) + "_METRIC_" + \
                 str(model) + "_" + str(member)
# json pattern
# json_pattern = "yplanton_" + str(metric_collection) + "_" + str(model) + "_vs_OBSERVATION_METRIC.json"
json_pattern = str(metric_collection) + "_v20210611_allModels_allRuns.json"
# netcdf pattern
# netcdf_pattern = "yplanton_" + str(metric_collection) + "_" + str(model) + "_vs_OBSERVATION_METRIC.nc"
netcdf_pattern = str(metric_collection) + "_" + str(project) + "_" + str(model) + "_" + str(experiment) + "_" + \
                 str(member) + "_METRIC.nc"
# list recipes
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
        json_file = json_pattern.replace("OBSERVATION", obs2).replace("METRIC", met)
        json_file = list(GLOBiglob(OSpath__join(path_j, json_file)))[0]
        with open(json_file) as ff:
            data_json = json.load(ff)["RESULTS"]["model"][model][member]
        ff.close()
        # get NetCDF file
        netcdf_file = netcdf_pattern.replace("OBSERVATION", obs2).replace("METRIC", met)
        netcdf_file = list(GLOBiglob(OSpath__join(path_n, netcdf_file)))[0]
        # get diagnostic values for the given model and observations
        dict_dia = data_json["value"][met]["diagnostic"]
        diagnostic_values = dict((key1, dict_dia[key1]["value"]) for key1 in list(dict_dia.keys()))
        diagnostic_units = data_json["metadata"][met]["diagnostic"]["units"]
        # get metric values computed with the given model and observations
        dict_met = data_json["value"][met]["metric"]
        metric_values = dict((key1, dict_met[key1]["value"]) for key1 in list(dict_met.keys()))
        metric_units = data_json["metadata"][met]["metric"]["units"]
        # figure name
        name_png = OSpath__join(path_p, figure_pattern.replace("METRIC", met) + "_" + str(obs2))
        # plot
        plotter_1mod_1obs(
            metric_collection, met, project, model, experiment, netcdf_file, diagnostic_values, diagnostic_units,
            metric_values, metric_units, name_png, reference=obs, member=member)
        # delete
        del diagnostic_units, diagnostic_values, dict_dia, dict_met, ff, json_file, metric_values, metric_units, \
            name_png, netcdf_file, obs2