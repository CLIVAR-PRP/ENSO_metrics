# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
import json
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
from EnsoMetricPlot import main_plotter


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "ENSO_perf"
models = ["IPSL-CM6A-LR"] * 32 + ["IPSL-CM5A-LR"] * 6 #["CNRM-CM5", "CNRM-CM6-1", "CNRM-CM6-1"]# "CNRM-ESM2-1", "IPSL-CM6A-LR"]
experiment = "historical" #"piControl" #
members = ["r"+str(ii+1)+"i1p1f1" for ii in range(32)] + ["r"+str(ii+1)+"i1p1" for ii in range(6)] #["r1i1p1", "r1i1p1f2", "r1i1p2f2"]#, "r1i1p1f2", "r1i1p1f1"]


path_in = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/eval_IPSL-CM6/2019_06_24_eval/Data"
#"/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/Eval_CNRM/v20190901/Data"
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/eval_IPSL-CM6/2019_06_24_eval/Plots_v2"
#"/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/Eval_CNRM/v20190901/Plots_v2"

expe = "hist" if experiment == "historical" else "pi"
pattern = "yplanton_" + metric_collection + "_modelname_" + experiment + "_membername"


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
filename_js = OSpath__join(path_in, "yplanton_" + metric_collection + "_all_dataset_" + experiment + "_raw.json")
with open(filename_js) as ff:
    data_json = json.load(ff)
ff.close()
del ff, filename_js
# loop on metrics
metrics = sorted(list(defCollection(metric_collection)['metrics_list'].keys()), key=lambda v: v.upper())
if metric_collection == "ENSO_perf":
    # metrics = [
    #     "BiasPrLatRmse", "BiasPrLonRmse", "BiasSstLatRmse", "BiasSstLonRmse", "BiasTauxLatRmse", "BiasTauxLonRmse",
    #     "EnsoAmpl", "EnsoDuration", "EnsoPrTsRmse", "EnsoSeasonality", "EnsoSstLonRmse", "EnsoSstSkew", "EnsoSstTsRmse",
    #     "EnsoTauxTsRmse", "NinoSstDiversity_1", "NinoSstDiversity_2", "SeasonalPrLatRmse", "SeasonalPrLonRmse",
    #     "SeasonalSstLatRmse", "SeasonalSstLonRmse", "SeasonalTauxLatRmse", "SeasonalTauxLonRmse"]
    metrics = ["SeasonalPrLonRmse", "SeasonalSstLonRmse", "SeasonalTauxLonRmse"]
elif metric_collection == "ENSO_proc":
    # metrics = ["EnsoAmpl", "EnsoFbSshSst", "EnsoFbSstLhf", "EnsoFbSstLwr", "EnsoFbSstShf", "EnsoFbSstSwr",
    #            "EnsoFbSstTaux", "EnsoFbSstThf", "EnsoFbTauxSsh", "EnsodSstOce_1", "EnsodSstOce_2"]
    metrics = ["EnsoFbSshSst", "EnsoFbSstLhf", "EnsoFbSstLwr", "EnsoFbSstShf", "EnsoFbSstSwr",
               "EnsoFbSstTaux", "EnsoFbSstThf", "EnsoFbTauxSsh"]
for met in metrics:
    print(met)
    filename_nc, modname = list(), list()
    for ii in range(len(models)):
        tmp = pattern.replace("modelname", models[ii]).replace("membername", members[ii]) + "_" + met + ".nc"
        filename_nc.append(OSpath__join(path_in, tmp))
        modname.append(models[ii] + "_" + members[ii])
    del tmp
    diagnostic_values = dict((key1, data_json['value'][met]["diagnostic"][key1]["value"])
                             for key1 in list(data_json['value'][met]["diagnostic"].keys()))
    metric_values = dict((key1, dict((key2, data_json['value'][met]["metric"][key1][key2]["value"])
                                     for key2 in list(data_json['value'][met]["metric"][key1].keys())))
                         for key1 in list(data_json['value'][met]["metric"].keys()))
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
    main_plotter(metric_collection, met, modname, experiment, filename_nc, diagnostic_values, diagnostic_units,
                 metric_values, metric_units, path_png=path_out)

