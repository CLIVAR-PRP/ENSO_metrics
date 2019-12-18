# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function

import json
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoMetrics.EnsoCollectionsLib import defCollection
from EnsoMetricPlot import main_plotter


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "ENSO_perf"
mip = "cmip5"
model = "CNRM-CM5"
#experiment = "piControl" # "historical"
experiment = "historical"
member = "r1i1p1"
modname = model + "_" + member

# path_in = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_08_report/Data"
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_08_report/Plots"
#path_in = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/Eval_CNRM/v20190901/Data"
#path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/Eval_CNRM/v20190901/Plots"
path_in_js = "/work/lee1043/imsi/result_test/metrics_results/enso_metric/cmip5/historical/v20191204/"+metric_collection
path_in_nc = "/work/lee1043/imsi/result_test/diagnostic_results/enso_metric/cmip5/historical/v20191204/"+metric_collection
path_out = "/work/lee1043/imsi/result_test/graphics/enso_metric/cmip5/historical/v20191204/"+metric_collection

expe = "hist" if experiment == "historical" else "pi"
# pattern = "yplanton_" + metric_collection + "_" + model + "_" + expe + "_" + member
#pattern = "yplanton_" + metric_collection + "_" + model + "_" + experiment + "_" + member
pattern = "_".join([mip, experiment, metric_collection, "v20191204", model, member])


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
# filename_js = OSpath__join(path_in, "yplanton_" + metric_collection + "_all_dataset_" + expe + "_raw.json")
#filename_js = OSpath__join(path_in, "yplanton_" + metric_collection + "_all_dataset_" + experiment + "_raw.json")
filename_js = OSpath__join(path_in_js, pattern + ".json")
with open(filename_js) as ff:
    print('filename_js:', filename_js)
    #data_json = json.load(ff)
    data_json_raw = json.load(ff)
    data_json = data_json_raw["RESULTS"]["model"][model][member]
ff.close()
del ff, filename_js
# loop on metrics
metrics = sorted(defCollection(metric_collection)['metrics_list'].keys(), key=lambda v: v.upper())
for met in metrics:
    print('met:', met)
    filename_nc = OSpath__join(path_in_nc, pattern + "_" + met + ".nc")
    print('filename_nc:', filename_nc)
    # diagnostic_values = dict((key, data_json['value'][met]["diagnostic"][key]["value"])
    #                          for key in data_json['value'][met]["diagnostic"].keys())
    # metric_values = dict((key, data_json['value'][met]["metric"][key]["value"])
    #                      for key in data_json['value'][met]["metric"].keys())
    # diagnostic_units = data_json['metadata']['metrics'][met]['diagnostic']['units']
    # metric_units = data_json['metadata']['metrics'][met]['metric']['units']
    diagnostic_values = dict((key1, data_json['value'][met]["diagnostic"][key1]["value"])
                             for key1 in list(data_json['value'][met]["diagnostic"].keys()))
    """
    metric_values = dict((key1, dict((key2, data_json['value'][met]["metric"][key1][key2]["value"])
                                     for key2 in list(data_json['value'][met]["metric"][key1].keys())))
                         for key1 in list(data_json['value'][met]["metric"].keys()))
    """
    """
    for key1 in list(data_json['value'][met]["metric"].keys()):
        print('key1:', key1)
        for key2 in list(data_json['value'][met]["metric"][key1].keys()):
            print('key2:', key2)
    """
    metric_values = dict((key2, {model: data_json['value'][met]["metric"][key2]["value"]})
                          for key2 in data_json['value'][met]["metric"].keys())

    diagnostic_units = data_json['metadata']['metrics'][met]['diagnostic']['units']
    metric_units = data_json['metadata']['metrics'][met]['metric']['units']
    # this function needs:
    #      - the name of the metric collection: metric_collection
    #      - the name of the metric: metric
    #      - the name of the model: modname (!!!!! this must be the name given when computed because it is the name used
    #                                              for in the netCDF files and in the json file !!!!!)
    #      - name of the experiment: experiment
    #      - name of the netCDF file name and path: filename_nc
    #      - a dictionary containing the diagnostic values: diagnostic_values (e.g., {'ERA-Interim': 1, 'Tropflux': 1.1,
    #                                                                                 modname: 1.5})
    #      - the diagnostic units: diagnostic_units
    #      - a dictionary containing the metric values: metric_values (e.g., {'ERA-Interim': 1.5, 'Tropflux': 1.36})
    #      - the metric units: metric_units
    #      - the path where to save the plots: path_out
    #main_plotter(metric_collection, met, modname, experiment, filename_nc, diagnostic_values,
    main_plotter(metric_collection, met, model, experiment, filename_nc, diagnostic_values,
                 diagnostic_units, metric_values, metric_units, path_png=path_out)



