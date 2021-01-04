# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create divedown plots for ENSO_metrics
# netCDF files are not provided in the package to do these plots
# Updated json files (needed to create this plot) can be downloaded from the page "Summary statistics in Interactive
# Portrait Plots" at https://cmec.llnl.gov/results/enso/
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#

from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from os.path import join as OSpath__join

# ENSO_metrics functions
from EnsoPlots.EnsoMetricPlot import main_plotter
from EnsoPlots.EnsoPlotToolsLib import remove_metrics


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "ENSO_tel"
project = "CMIP5"  # "obs2obs"  #
model = "CNRM-CM5"  # "ERA-Interim_SODA3.4.2"  # "ERA-Interim_ERA-Interim"  # "ERA-Interim"  #
experiment = "historical"
member = "r1i1p1"
dataname = deepcopy(model) if project == "obs2obs" else model + "_" + member
plot_ref = True if project == "obs2obs" else False
# True to use the set of metric in the BAMS paper
# More metric have been computed and tested but not kept
reduced_set = True  # False  #
# computation version, 'v20200427' is provided with the package
version = "v20200427"
# json files
dict_json = {
    "CMIP5": {
        "ENSO_perf": "share/EnsoMetrics/cmip5_historical_ENSO_perf_" + version + "_allModels_allRuns.json",
        "ENSO_proc": "share/EnsoMetrics/cmip5_historical_ENSO_proc_" + version + "_allModels_allRuns.json",
        "ENSO_tel": "share/EnsoMetrics/cmip5_historical_ENSO_tel_" + version + "_allModels_allRuns.json"},
    "CMIP6": {
        "ENSO_perf": "share/EnsoMetrics/cmip6_historical_ENSO_perf_" + version + "_allModels_allRuns.json",
        "ENSO_proc": "share/EnsoMetrics/cmip6_historical_ENSO_proc_" + version + "_allModels_allRuns.json",
        "ENSO_tel": "share/EnsoMetrics/cmip6_historical_ENSO_tel_" + version + "_allModels_allRuns.json"},
    "obs2obs": {
        "ENSO_perf": "share/EnsoMetrics/obs2obs_ENSO_perf_" + version + ".json",
        "ENSO_proc": "share/EnsoMetrics/obs2obs_ENSO_proc_" + version + ".json",
        "ENSO_tel":  "share/EnsoMetrics/obs2obs_ENSO_tel_" + version + ".json"}}

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2020_05_report"
path_nc = OSpath__join(path_main, "Data/" + project.lower() + "/" + experiment + "/" + metric_collection)
# figure name
path_out = ""
dataname2 = dataname.replace("GPCPv2.3", "GPCPv23").replace("SODA3.4.2", "SODA342")
figure_name = project.lower() + "_" + experiment + "_" + metric_collection + "_" + dataname2
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
with open(dict_json[project][metric_collection]) as ff:
    data_json = json.load(ff)['RESULTS']['model'][model][member]
ff.close()
del ff
# get metric names
list_metrics = sorted(list(data_json["value"].keys()), key=lambda v: v.upper())
if reduced_set is True:
    metrics = remove_metrics(list_metrics, metric_collection, reduced_set=reduced_set)
# pattern of netCDF files
pattern = project.lower() + "_" + experiment + "_" + metric_collection + "_" + version + "_" + dataname
#
# Loop on metrics
#
for met in ["EnsoPrMapDjfRmse"]:#list_metrics:
    print(met)
    # get NetCDF file name
    met2 = met.replace("Rmse", "") if metric_collection in ["ENSO_tel"] and "Map" in met else deepcopy(met)
    filename_nc = pattern + "_" + met2 + ".nc"
    filename_nc = list(GLOBiglob(OSpath__join(path_nc, filename_nc)))[0]
    # get diagnostic values for the given model and observations
    dict_dia = data_json["value"][met]["diagnostic"]
    diagnostic_values = dict((key1, dict_dia[key1]["value"]) for key1 in list(dict_dia.keys()))
    diagnostic_units = data_json["metadata"]["metrics"][met]["diagnostic"]["units"]
    # get metric values computed with the given model and observations
    if metric_collection in ["ENSO_tel"] and "Map" in met:
        list1, list2 = [met.replace("Rmse", "Corr"), met], ["metric", "metric"]
        dict_met = data_json["value"]
        metric_values = dict((key1, {model: [1-dict_met[su][ty][key1]["value"] if "Corr" in su else
                                             dict_met[su][ty][key1]["value"]for su, ty in zip(list1, list2)]})
                             for key1 in list(dict_met[list1[0]]["metric"].keys()))
        metric_units = [data_json["metadata"]["metrics"][su]["metric"]["units"] for su in list1]
        del list1, list2
    else:
        dict_met = data_json["value"][met]["metric"]
        metric_values = dict((key1, {model: dict_met[key1]["value"]}) for key1 in list(dict_met.keys()))
        metric_units = data_json["metadata"]["metrics"][met]["metric"]["units"]
    # figure name
    name_png = figure_name + "_" + met
    # this function needs:
    #      - the name of the metric collection: metric_collection
    #      - the name of the metric: metric
    #      - the name of the model: model
    #      - name of the experiment: experiment
    #      - name of the netCDF file name and path: filename_nc
    #      - a dictionary containing the diagnostic values: diagnostic_values (e.g., {"ERA-Interim": 1, "Tropflux": 1.1,
    #                                                                                 model: 1.5})
    #      - the diagnostic units: diagnostic_units
    #      - a dictionary containing the metric values: metric_values (e.g., {"ERA-Interim": {model: 1.5},
    #                                                                         "Tropflux": {model: 1.36}})
    #      - the metric units: metric_units
    #      - (optional) the member name, not needed if project aims to compare observational datasets (obs2obs): member
    #      - (optional) the path where to save the plots: path_png
    #      - (optional) the name of the plots: name_png
    #      - (optional) if the project aims to compare observational datasets (obs2obs): plot_ref
    main_plotter(metric_collection, met2, model, experiment, filename_nc, diagnostic_values, diagnostic_units,
                 metric_values, metric_units, member=member, path_png=path_out, name_png=figure_name, plot_ref=plot_ref)
    del diagnostic_values, diagnostic_units, dict_dia, dict_met, filename_nc, met2, metric_values, metric_units, \
        name_png
