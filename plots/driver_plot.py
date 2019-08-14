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
metric_collection = "ENSO_test"
model = "CNRM-CM5"
experiment = "historical"
member = "r1i1p1"
modname = model + "_" + member

path_in = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_08_report/Data"
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_08_report/Plots"

expe = "hist" if experiment == "historical" else "pi"
pattern = "yplanton_" + metric_collection + "_" + model + "_" + expe + "_" + member

# this is used because I failed to save the metrics units... But it should be in the json file
dict_units = {
    "BiasPrLatRmse": "mm/day", "BiasPrLonRmse": "mm/day", "BiasPrRmse": "mm/day",
    "BiasSstLatRmse": "C", "BiasSstLonRmse": "C", "BiasSstRmse": "C",
    "BiasTauxLatRmse": "1e-3 N/m2", "BiasTauxLonRmse": "1e-3 N/m2", "BiasTauxRmse": "1e-3 N/m2",
    "EnsoAmpl": "C", "EnsoDuration": "months", "EnsodSstOce": "C/C", "EnsoSeasonality": "C/C", "EnsoSstSkew": "C",
    "EnsoFbSshSst": "C/cm",
    "EnsoFbSstLhf": "W/m2/C", "EnsoFbSstLwr": "W/m2/C", "EnsoFbSstShf": "W/m2/C", "EnsoFbSstSwr": "W/m2/C",
    "EnsoFbSstThf": "W/m2/C", "EnsoFbSstTaux": "1e-3 N/m2/C", "EnsoFbTauxSsh": "1e3 cm/N/m2",
    "EnsoPrMap": "mm/day/C", "EnsoSlpMap": "hPa/C", "EnsoSstMap": "C/C", "EnsoSstLonRmse": "C/C",
    "EnsoPrTsRmse": "mm/day/C", "EnsoSstTsRmse": "C/C", "EnsoTauxTsRmse": "1e-3 N/m2/C",
    "NinaPrMap": "mm/day", "NinaSlpMap": "hPa", "NinaSstMap": "C", "NinaSstLonRmse": "C", "NinaSstDur": "months",
    "NinaPrTsRmse": "mm/day", "NinaSstTsRmse": "C", "NinaTauxTsRmse": "1e-3 N/m2",
    "NinoPrMap": "mm/day", "NinoSlpMap": "hPa", "NinoSstMap": "C", "NinoSstLonRmse": "C", "NinoSstDur": "months",
    "NinoPrTsRmse": "mm/day", "NinoSstTsRmse": "C", "NinoTauxTsRmse": "1e-3 N/m2",
    "NinoSstDiversity": "longitude",
    "SeasonalPrLatRmse": "mm/day", "SeasonalPrLonRmse": "mm/day",
    "SeasonalSstLatRmse": "C", "SeasonalSstLonRmse": "C",
    "SeasonalTauxLatRmse": "1e-3 N/m2", "SeasonalTauxLonRmse": "1e-3 N/m2",
}

# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
filename_js = OSpath__join(path_in, "yplanton_" + metric_collection + "_all_dataset_" + expe + ".json")
with open(filename_js) as ff:
    data_json = json.load(ff)
ff.close()
del ff, filename_js
# loop on metrics
metrics = sorted(defCollection(metric_collection)['metrics_list'].keys(), key=lambda v: v.upper())
for met in metrics:
    print met
    filename_nc = OSpath__join(path_in, pattern + "_" + met + ".nc")
    diagnostic_values = dict((key, data_json[met][model + "_" + member]["diagnostic"][key]["value"])
                             for key in data_json[met][model + "_" + member]["diagnostic"].keys())
    metric_values = dict((key, data_json[met][modname]["metric"][key]["value"])
                         for key in data_json[met][model + "_" + member]["metric"].keys())
    key = data_json[met][model + "_" + member]["diagnostic"].keys()[0]
    diagnostic_units = data_json[met][model + "_" + member]["diagnostic"][key]["units"]
    #key = data_json[met][model + "_" + member]["metric"].keys()[0]
    #metric_units = data_json[met][model + "_" + member]["metric"][key]["units"]
    metric_units = dict_units[met]  # trick... units should come from the json file
    main_plotter(metric_collection, met, model, experiment, member, filename_nc, diagnostic_values,
                 diagnostic_units, metric_values, metric_units, path_png=path_out)



