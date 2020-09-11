# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots to compare groups of models
#      FIG. 3 in Planton et al. 2020: Evaluating climate models with the CLIVAR 2020 ENSO metrics package. BAMS
# It uses the first available member of each model or all members of each model and averages them
# Updated json files (needed to create this plot) can be downloaded from the page "Summary statistics in Interactive
# Portrait Plots" at https://cmec.llnl.gov/results/enso/
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#

from copy import deepcopy
from numpy import array as NUMPYarray
from numpy import mean as NUMPYmean
from numpy import moveaxis as NUMPYmoveaxis
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYma__masked_where
from os.path import join as OSpath__join

# set of functions to find cmip/obs files and save a json file
# to be adapted/changed by users depending on their environments
from driver_tools_lib import get_metric_values, get_mod_mem_json

# ENSO_metrics functions
from EnsoPlots.EnsoPlotTemplate import plot_projects_comparison
from EnsoPlots.EnsoPlotToolsLib import bootstrap, sort_metrics


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
# metric collections to plot
list_metric_collections = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
# CMIP experiment
experiment = "historical"
# project to use, here both CMIP5 and CMIP6 models will be used
list_projects = ["CMIP6", "CMIP5"]
# True to use the set of metric in the BAMS paper
# More metric have been computed and tested but not kept
reduced_set = True  # False  #
# False to projects defined in 'list_project'
# If set to True, all projects defined in 'list_project' will be used as one and will be compared to a given selection
# of models (see 'my_project' and 'my_selection')
big_ensemble = False  # True  #
# marker colors
if big_ensemble is False:
    colors = ["r", "dodgerblue"]
else:
    colors = ["orange", "forestgreen"]
# True to use the first available member only
# If set to False, all members will be used and the metric values computed for all members of each model will be
# averaged
first_member = True  # False  #
# If 'big_ensemble' is set to True, 'BAMS_teleconnection' will be compared to 'CMIP', all projects defined in
# 'list_project' used as one
my_project = ["BAMS_teleconnection", "CMIP"]
# Definition of selection to use if 'big_ensemble' is set to True
my_selection = {
    "BAMS_teleconnection": [
        'CESM2', 'CESM2-FV2', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'CMCC-CM', 'CNRM-CM5', 'CNRM-CM5-2', 'EC-Earth3',
        'EC-Earth3-Veg', 'FGOALS-f3-L', 'FGOALS-s2', 'GFDL-CM4', 'GFDL-ESM4', 'MIROC-ES2L', 'MIROC6', 'NESM3',
        'NorESM2-MM']}
# List of additional observations
# the reading part is very 'ad hoc', do not change the obs!
list_obs = ["20CRv2", "NCEP2", "ERA-Interim"]
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
# figure name
path_out = ""
figure_name = "metrics_intercomparison_" + str(len(list_metric_collections)) + "metric_collections_" + version
if len(list_projects) == 1:
    figure_name += "_" + str(list_projects[0])
else:
    figure_name += "_" + str(len(list_projects)) + "cmip"
if big_ensemble is False:
    figure_name += "_" + list_projects[1] + "_vs_" + list_projects[0]
else:
    figure_name += "_" + my_project[1] + "_vs_" + my_project[0]
if first_member is True:
    figure_name += "_first_member"
else:
    figure_name += "_members_averaged"
if reduced_set is False:
    figure_name += "_all_metrics"
figure_name = OSpath__join(path_out, figure_name)
# ---------------------------------------------------#


# ---------------------------------------------------#
# Functions
# ---------------------------------------------------#
def common_save(dict_in, dict_out={}):
    for mod in list(dict_in.keys()):
        try:    dict_out[mod]
        except: dict_out[mod] = dict_in[mod]
        else:
            for met in list(dict_in[mod].keys()):
                dict_out[mod][met] = dict_in[mod][met]
    return dict_out
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# get members by model by project from json file
# only metrics from models/members chosen here will be used
# all metrics from models/members chosen here will be used (ensures that if a model/member is not available for one or
# several metric collections, the corresponding line will still be created in the portraitplot)
model_by_proj = get_mod_mem_json(list_projects, list_metric_collections, dict_json, first_only=first_member)
# read json file
dict_met = dict()
if big_ensemble is False:
    for proj in list_projects:
        dict_mc = dict()
        for mc in list_metric_collections:
            dict1 = get_metric_values(proj, mc, dict_json, model_by_proj, reduced_set=reduced_set)
            # save in common dictionary
            dict_mc = common_save(dict1, dict_out=dict_mc)
        dict_met[proj] = dict_mc
        del dict_mc
else:
    dict_mc = dict()
    for proj in list_projects:
        for mc in list_metric_collections:
            dict1 = get_metric_values(proj, mc, dict_json, model_by_proj, reduced_set=reduced_set)
            # save in common dictionary
            dict_mc = common_save(dict1, dict_out=dict_mc)
    dict_met["CMIP"] = dict_mc
    # put the selected models in a separate key
    for mod in my_selection[my_project[0]]:
        try:    dict_met[my_project[0]]
        except: dict_met[my_project[0]] = {mod: dict_met["CMIP"][mod]}
        else:   dict_met[my_project[0]][mod] = dict_met["CMIP"][mod]
    del dict_mc


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
if ' ':
    list_metrics = list()
    for k1 in list(dict_met.keys()):
        for k2 in list(dict_met[k1].keys()):
            list_metrics += list(dict_met[k1][k2].keys())
    list_metrics = sort_metrics(list(set(list_metrics)))
    opposed_groups = deepcopy(list_projects) if big_ensemble is False else deepcopy(my_project)
    # mean metric evaluation
    tab_bst, tab_val = list(), list()
    for met in list_metrics:
        tab_tmp = list()
        for grp in opposed_groups:
            tab = list()
            for mod in list(dict_met[grp].keys()):
                if met in list(dict_met[grp][mod].keys()):
                    if dict_met[grp][mod][met] is not None and dict_met[grp][mod][met] != 1e20:
                        tab.append(dict_met[grp][mod][met])
            tab = NUMPYarray(tab)
            tab_tmp.append(NUMPYma__masked_invalid(tab).compressed())
            del tab
        tab1, tab2 = list(), list()
        for ii in range(len(tab_tmp)):
            tab1.append(float(NUMPYmean(tab_tmp[ii])))
            nbr = nbr = len(tab_tmp[1]) if ii==0 else len(tab_tmp[0])
            bst = bootstrap(tab_tmp[ii], nech=nbr)
            tab2.append(bst)
            del bst, nbr
        tab_bst.append(tab2)
        tab_val.append(tab1)
    tab_bst = NUMPYmoveaxis(NUMPYarray(tab_bst), 0, 1)
    tab_bst = NUMPYma__masked_where(tab_bst == 1e20, tab_bst)
    tab_val = NUMPYmoveaxis(NUMPYarray(tab_val), 0, -1)
    tmp = NUMPYmoveaxis(NUMPYarray([tab_val[1], tab_val[1]]), 0, 1)
    tab_bst = tab_bst / tmp
    tab_val = tab_val / tab_val[1]
    # plot project comparison
    plot_projects_comparison(tab_val, figure_name, xticklabel=list_metrics, yticklabel=opposed_groups[1].upper(),
                             colors=colors, tab_bst=tab_bst, legend=opposed_groups, chigh=True, cfram=True)
    del list_metrics, opposed_groups, tab_bst, tab_val, tmp
