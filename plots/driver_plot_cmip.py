# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function
from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from math import isnan as MATHisnan
from numpy import array as NUMPYarray
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
from EnsoMetricPlot import cmip_plotter


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collections = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
projects = ["cmip5", "cmip6"]
my_project = "cmip5"
experiment = "historical"
multi_member_ave = False

mod_to_remove = ["CIESM", "E3SM-1-1-ECA", "EC-EARTH", "FIO-ESM", "FGOALS-g3", "GFDL-CM2p1", "HadGEM2-AO", "MCM-UA-1-0"]
mem_to_remove = ["r101i1p1f1", "r102i1p1f1"]

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2020_05_report"
path_out = OSpath__join(path_main, "Plots_wiki")

# json file name for CMIP
version = "v20200427"
dict_json = {
    "cmip5": {
        "ENSO_perf": OSpath__join(path_main, "Data/cmip5/" + experiment + "/cmip5_" + experiment + "_ENSO_perf_" +
                                  version + "_allModels_allRuns.json"),
        "ENSO_proc": OSpath__join(path_main, "Data/cmip5/" + experiment + "/cmip5_" + experiment + "_ENSO_proc_" +
                                  version + "_allModels_allRuns.json"),
        "ENSO_tel":  OSpath__join(path_main, "Data/cmip5/" + experiment + "/cmip5_" + experiment + "_ENSO_tel_" +
                                  version + "_allModels_allRuns.json")},
    "cmip6": {
        "ENSO_perf": OSpath__join(path_main, "Data/cmip6/" + experiment + "/cmip6_" + experiment + "_ENSO_perf_" +
                                  version + "_allModels_allRuns.json"),
        "ENSO_proc": OSpath__join(path_main, "Data/cmip6/" + experiment + "/cmip6_" + experiment + "_ENSO_proc_" +
                                  version + "_allModels_allRuns.json"),
        "ENSO_tel":  OSpath__join(path_main, "Data/cmip6/" + experiment + "/cmip6_" + experiment + "_ENSO_tel_" +
                                  version + "_allModels_allRuns.json")}}


# ---------------------------------------------------#
# Functions
# ---------------------------------------------------#
def find_first_member(members):
    """
    Finds the first member

    Inputs:
    ------
    :param members: list of string
        List of members.

    Output:
    ------
    :return mem: string
        First member of the given list.
    """
    if "r1i1p1" in members:
        mem = "r1i1p1"
    elif "r1i1p1f1" in members:
        mem = "r1i1p1f1"
    elif "r1i1p1f2" in members:
        mem = "r1i1p1f2"
    else:
        tmp = deepcopy(members)
        members = list()
        for mem in tmp:
            for ii in range(1, 10):
                if "r"+str(ii)+"i" in mem:
                    members.append(mem.replace("r"+str(ii)+"i", "r"+str(ii).zfill(2)+"i"))
                else:
                    members.append(mem)
        del tmp
        mem = sorted(list(set(members)), key=lambda v: v.upper())[0].replace("r0", "r")
    return mem


def find_observed_diag(lobs, lmet, dict_mod, dict_val):
    val = None
    found = False
    for mod in dict_mod.keys():
        for mem in dict_mod[mod]:
            if mod in dict_val.keys() and mem in dict_val[mod].keys() and lmet in dict_val[mod][mem]["value"].keys():
                if lobs in dict_val[mod][mem]["value"][lmet]["diagnostic"].keys():
                    val = dict_val[mod][mem]["value"][lmet]["diagnostic"][lobs]["value"]
                    found = True
                    break
        if found is True:
            break
    return val


def find_units(lmet, dict_mod, dict_val, val_type):
    uni = None
    found = False
    for mod in dict_mod.keys():
        for mem in dict_mod[mod]:
            if mod in dict_val.keys() and mem in dict_val[mod].keys() and\
                    lmet in dict_val[mod][mem]["metadata"]["metrics"].keys():
                uni = dict_val[mod][mem]["metadata"]["metrics"][lmet][val_type]["units"]
                found = True
                break
        if found is True:
            break
    return uni
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# find models and members (ensures that same models and members are used for each metric)
models = dict()
for pro in projects:
    dict1 = dict()
    for mc in metric_collections:
        with open(dict_json[pro][mc]) as ff:
            data_json = json.load(ff)["RESULTS"]["model"]
        ff.close()
        for mod in data_json.keys():
            try:    dict1[mod]
            except: dict1[mod] = data_json[mod].keys()
            else:   dict1[mod] += data_json[mod].keys()
        del data_json, ff
    dict2 = dict()
    for mod in dict1.keys():
        # remove given models
        if mod not in mod_to_remove:
            # remove duplicated and given members
            dict2[mod] = [mem for mem in sorted(list(set(dict1[mod])), key=lambda v: v.upper())
                          if mem not in mem_to_remove]
    models[pro] = dict2
    del dict1, dict2

# loop on metric collections
for mc in metric_collections:
    print(mc)
    # loop on metrics
    metrics = sorted(defCollection(mc)['metrics_list'].keys(), key=lambda v: v.upper())
    metrics = [met + "Rmse" if "Map" in met else met for met in metrics]
    for met in metrics:
        print(str().ljust(5)+met)
        # find every possible observations for the given metric
        observations = list()
        for pro in projects:
            with open(dict_json[pro][mc]) as ff:
                data_json = json.load(ff)["RESULTS"]["model"]
            ff.close()
            for mod in models[pro].keys():
                for mem in models[pro][mod]:
                    if mod in data_json.keys() and mem in data_json[mod].keys() and\
                            met in data_json[mod][mem]["value"].keys():
                        observations += data_json[mod][mem]["value"][met]["metric"].keys()
            del data_json, ff
        observations = sorted(list(set(observations)), key=lambda v: v.upper())
        # get diagnostic and metric values
        diagnostic_values, metric_values = dict(), dict()
        diagnostic_units, metric_units = "", ""
        for pro in projects:
            # open json
            with open(dict_json[pro][mc]) as ff:
                data_json = json.load(ff)["RESULTS"]["model"]
            ff.close()
            # loop on models
            list_dia, dict_met = list(), dict()
            for mod in models[pro].keys():
                members = models[pro][mod] if multi_member_ave is True else [find_first_member(models[pro][mod])]
                # loop on members
                l1, l2 = list(), dict()
                for mem in models[pro][mod]:
                    if mod in data_json.keys() and mem in data_json[mod].keys() and\
                            met in data_json[mod][mem]["value"].keys():
                        dict1 = data_json[mod][mem]["value"][met]
                        # keep diagnostic value only if it is not None
                        tmp = dict1["diagnostic"][mod + "_" + mem]["value"]
                        if tmp is not None and MATHisnan(tmp) is False and tmp < 1e10:
                            l1.append(dict1["diagnostic"][mod + "_" + mem]["value"])
                        for obs in observations:
                            if obs in dict1["metric"].keys():
                                tmp = dict1["metric"][obs]["value"]
                                if tmp is not None and MATHisnan(tmp) is False and tmp < 1e10:
                                    try:    l2[obs]
                                    except: l2[obs] = [tmp]
                                    else:   l2[obs] += [tmp]
                        del dict1, tmp
                # set model value to
                #      None if there is no data available
                #      first value if only one value is available
                #      multi-member mean if more than one value is available
                tmp1 = None if len(l1) == 0 else (l1[0] if len(l1) == 1 else NUMPYarray(l1).mean())
                tmp2 = dict((obs, l2[obs][0] if len(l2[obs]) == 1 else NUMPYarray(l2[obs]).mean()) for obs in l2.keys())
                # keep model value only if it is not None
                if tmp1 is not None:
                    list_dia.append(tmp1)
                for obs in observations:
                    if obs in tmp2.keys():
                        try:    dict_met[obs]
                        except: dict_met[obs] = [tmp2[obs]]
                        else:   dict_met[obs] += [tmp2[obs]]
                del l1, l2, members, tmp1, tmp2
            diagnostic_values[pro] = list_dia
            if pro == projects[-1]:
                for obs in observations:
                    tmp = find_observed_diag(obs, met, models[pro], data_json)
                    if tmp is not None:
                        try:    diagnostic_values["ref"]
                        except: diagnostic_values["ref"] = {obs: tmp}
                        else:   diagnostic_values["ref"][obs] = tmp
                    del tmp
                diagnostic_units =find_units(met, models[pro], data_json, "diagnostic")
                metric_units = find_units(met, models[pro], data_json, "metric")
            metric_values[pro] = dict_met
            del dict_met, ff, list_dia
        del observations
        # figure name
        figure_name = OSpath__join(path_out, my_project + "_" + experiment + "_" + mc + "_" + my_project + "_" + met)
        if multi_member_ave is True:
            figure_name += "_multi_member_mean"
        else:
            figure_name += "_first_member"
        # this function needs:
        #      - the name of the metric collection: metric_collection (e.g., "ENSO_perf")
        #      - the name of the metric: metric (e.g., "EnsoAmpl")
        #      - a dictionary containing the diagnostic values: diagnostic_values
        #        e.g., {"ref": {"ERA-Interim": 1, "Tropflux": 1.1}, "cmip5": [0.9, ...], "cmip6": [0.9, ...]}
        #      - the diagnostic units: diagnostic_units (e.g., "C")
        #      - a dictionary containing the metric values: metric_values
        #        e.g., {"cmip5": {"ERA-Interim": [10, ...], "Tropflux": [20, ...]},
        #               "cmip6": {"ERA-Interim": [10, ...], "Tropflux": [20, ...]}}
        #      - the metric units: metric_units (e.g., "%")
        #      - multi-member average information: True if multi-member mean is used, False otherwise
        #      - the path and name of the plots: name_png (e.g., "path/to/directory/file_name")
        cmip_plotter(mc, met, experiment, diagnostic_values, diagnostic_units, metric_values, metric_units,
                     multi_member_ave, figure_name)
