# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function
<<<<<<< HEAD
from glob import iglob as GLOBiglob
import json
=======

# Run matplotlib background to prevent 
# display localhost error after console disconnected
# and to speed up
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# Import other libs
from glob import iglob as GLOBiglob
import json
from os import makedirs as OS__makedirs
from os.path import exists as OSpath__exists
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
from os.path import join as OSpath__join
# ENSO_metrics functions
#from EnsoCollectionsLib import defCollection
from EnsoMetrics.EnsoCollectionsLib import defCollection
from EnsoMetricPlot import main_plotter
<<<<<<< HEAD

=======
import sys

from PMPdriver_lib import AddParserArgument
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946

# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
<<<<<<< HEAD
=======
param = AddParserArgument()

# Metrics Collection
metric_collection = param.metricsCollection

# Pre-defined options
mip = param.mip
exp = param.exp

# model
if param.modnames is None:
    model = "IPSL-CM5A-LR"
else:
    model = param.modnames[0]

# Realizations
run = param.realization

# case id
case_id = param.case_id

# Switches
debug = param.debug

"""
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
metric_collection = "ENSO_perf"
#metric_collection = "ENSO_tel"
#metric_collection = "ENSO_proc"

mip = "cmip5"
<<<<<<< HEAD
model = "CNRM-CM5"
#model = "IPSL-CM5A-LR"
experiment = "historical"
member = "r1i1p1"
modname = model + "_" + member

case_id = "v20200305"
debug = True

# ---------------------------------------------------#
=======
exp = "historical"
model = "IPSL-CM5A-LR"
run = "r1i1p1"

case_id = "v20200305"
debug = True
"""

# ---------------------------------------------------#
# Check Arguments
# ---------------------------------------------------#
print("metric_collection:", metric_collection)
print("mip:", mip)
print("exp:", exp)
print("model:", model)
print("run:", run)
print("case_id:", case_id)
print("debug:", debug)
# ---------------------------------------------------#
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946

path_main = "/p/user_pub/pmp/pmp_results/pmp_v1.1.2"
path_in_json = OSpath__join(path_main, "metrics_results", "enso_metric", mip, exp, case_id, metric_collection)
path_in_nc = OSpath__join(path_main, "diagnostic_results", "enso_metric", mip, exp, case_id, metric_collection)

if debug:
    path_main = "/work/lee1043/imsi/result_test"
path_out = OSpath__join(path_main, "graphics", "enso_metric", mip, exp, case_id, metric_collection)

<<<<<<< HEAD
if not os.path.exists(path_out):
    os.makedirs(path_out)
    print("path_out:", path_out)
=======
if not OSpath__exists(path_out):
    try:
        OS__makedirs(path_out)
        print("path_out:", path_out)
    except:
        pass
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946

pattern = "_".join([mip, exp, metric_collection, case_id])

# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
filename_js = OSpath__join(path_in_json, pattern + "_allModels_allRuns.json")
print('filename_js:', filename_js)
with open(filename_js) as ff:
<<<<<<< HEAD
    data_json = json.load(ff)['RESULTS']['model'][model][member]
=======
    data_json = json.load(ff)['RESULTS']['model'][model][run]
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
ff.close()
del ff, filename_js
# loop on metrics
metrics = sorted(defCollection(metric_collection)['metrics_list'].keys(), key=lambda v: v.upper())
for met in metrics:
    try:
        print('met:', met)
        # get NetCDF file name
<<<<<<< HEAD
        filename_nc = OSpath__join(path_in_nc, pattern + "_" + model + "_" + member + "_" + met + ".nc")
=======
        filename_nc = OSpath__join(path_in_nc, pattern + "_" + model + "_" + run + "_" + met + ".nc")
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
        print("filename_nc:", filename_nc)
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
<<<<<<< HEAD
        figure_name = "_".join([mip, exp, metric_collection, model, member, met])
=======
        figure_name = "_".join([mip, exp, metric_collection, model, run, met])
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
        # this function needs:
        #      - the name of the metric collection: metric_collection
        #      - the name of the metric: metric
        #      - the name of the model: modname (!!!!! this must be the name given when computed because it is the name used
        #                                              for in the netCDF files and in the json file !!!!!)
        #      - name of the exp: exp
        #      - name of the netCDF file name and path: filename_nc
        #      - a dictionary containing the diagnostic values: diagnostic_values (e.g., {"ERA-Interim": 1, "Tropflux": 1.1,
        #                                                                                 modname: 1.5})
        #      - the diagnostic units: diagnostic_units
        #      - a dictionary containing the metric values: metric_values (e.g., {"ERA-Interim": {modname: 1.5},
        #                                                                         "Tropflux": {modname: 1.36}})
        #      - the metric units: metric_units
        #      - (optional) the path where to save the plots: path_out
        #      - (optional) the name of the plots: name_png
        main_plotter(metric_collection, met, model, exp, filename_nc, diagnostic_values,
<<<<<<< HEAD
                     diagnostic_units, metric_values, metric_units, member=member, path_png=path_out, 
=======
                     diagnostic_units, metric_values, metric_units, member=run, path_png=path_out, 
>>>>>>> 7492f16b3aee130baff54a1c4dc6adf27c1b5946
                     name_png=figure_name)
    except Exception as e:
        print("## ERROR:", e)
        pass
