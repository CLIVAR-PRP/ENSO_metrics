# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function
from glob import iglob as GLOBiglob
import json
from os.path import join as OSpath__join
# ENSO_metrics functions
#from EnsoCollectionsLib import defCollection
from EnsoMetrics.EnsoCollectionsLib import defCollection
from EnsoMetricPlot import main_plotter


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "ENSO_tel"
metric_collection = "ENSO_perf"
metric_collection = "ENSO_proc"
project = "cmip5"
model = "CNRM-CM5"
#model = "IPSL-CM5A-LR"
experiment = "historical"
member = "r1i1p1"
modname = model + "_" + member

#path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_12_report"
path_main = "/work/lee1043/cdat/pmp/ENSO_metrics/plot_test_json_20191218"
path_in = OSpath__join(path_main, "Data_lee")
path_out = OSpath__join(path_main, "Plots_wiki")

expe = "hist" if experiment == "historical" else "pi"
pattern = project + "_" + experiment + "_" + metric_collection + "_v2019????"


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
filename_js = list(GLOBiglob(OSpath__join(path_in, pattern + "_allModels_allRuns_modified.json")))[0]
print('filename_js:', filename_js)
with open(filename_js) as ff:
    data_json = json.load(ff)['RESULTS']['model'][model][member]
ff.close()
del ff, filename_js
# loop on metrics
metrics = sorted(defCollection(metric_collection)['metrics_list'].keys(), key=lambda v: v.upper())
for met in metrics:
    print('met:', met)
    # get NetCDF file name
    path_nc = OSpath__join(path_in, project + "/" + experiment + "/" + metric_collection)
    filename_nc = list(GLOBiglob(OSpath__join(path_nc, pattern + "_" + model + "_" + member + "_" + met + ".nc")))
    if len(filename_nc) != 1:
        print('    Pass for metric', met, ': no NC file detected.')
        pass
    else:
        filename_nc = filename_nc[0]
        print('filename_nc:', filename_nc)
        # get diagnostic values for the given model and observations
        if metric_collection == "ENSO_tel" and "Map" in met:
            dict_dia = data_json["value"][met+"Corr"]["diagnostic"]
            diagnostic_values = dict((key1, None) for key1 in dict_dia.keys())
            diagnostic_units = ""
        else:
            dict_dia = data_json["value"][met]["diagnostic"]
            diagnostic_values = dict((key1, dict_dia[key1]["value"]) for key1 in dict_dia.keys())
            diagnostic_units = data_json["metadata"]["metrics"][met]["diagnostic"]["units"]
        # get metric values computed with the given model and observations
        if metric_collection == "ENSO_tel" and "Map" in met:
            list1, list2 = [met+"Corr", met+"Rmse"], ["diagnostic", "metric"]
            dict_met = data_json["value"]
            metric_values = dict((key1, {model: [dict_met[su][ty][key1]["value"] for su, ty in zip(list1, list2)]})
                                 for key1 in dict_met[list1[0]]["metric"].keys())
            metric_units = [data_json["metadata"]["metrics"][su]["metric"]["units"] for su in list1]
        else:
            dict_met = data_json["value"][met]["metric"]
            metric_values = dict((key1, {model: dict_met[key1]["value"]}) for key1 in dict_met.keys())
            metric_units = data_json["metadata"]["metrics"][met]["metric"]["units"]
        # figure name
        figure_name = project + "_" + experiment + "_" + metric_collection + "_" + model + "_" + member + "_" + met
        # this function needs:
        #      - the name of the metric collection: metric_collection
        #      - the name of the metric: metric
        #      - the name of the model: modname (!!!!! this must be the name given when computed because it is the name used
        #                                              for in the netCDF files and in the json file !!!!!)
        #      - name of the experiment: experiment
        #      - name of the netCDF file name and path: filename_nc
        #      - a dictionary containing the diagnostic values: diagnostic_values (e.g., {"ERA-Interim": 1, "Tropflux": 1.1,
        #                                                                                 modname: 1.5})
        #      - the diagnostic units: diagnostic_units
        #      - a dictionary containing the metric values: metric_values (e.g., {"ERA-Interim": {modname: 1.5},
        #                                                                         "Tropflux": {modname: 1.36}})
        #      - the metric units: metric_units
        #      - (optional) the path where to save the plots: path_out
        #      - (optional) the name of the plots: name_png
        main_plotter(metric_collection, met, model, experiment, filename_nc, diagnostic_values,
                     diagnostic_units, metric_values, metric_units, path_png=path_out, name_png=figure_name)
