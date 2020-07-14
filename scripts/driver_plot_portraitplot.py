# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create portraitplot
#      FIG. 2 in Planton et al. 2020: Evaluating climate models with the CLIVAR 2020 ENSO metrics package. BAMS
# It uses the first available member of each model or all members of each model and averages them
# Updated json files (needed to create this plot) can be downloaded from the page "Summary statistics in Interactive
# Portrait Plots" at https://cmec.llnl.gov/results/enso/
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function
from numpy import mean as NUMPYmean
from numpy import std as NUMPYstd
from numpy.ma import array as NUMPYma__array
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros
from os.path import join as OSpath__join
import string

# set of functions to find cmip/obs files and save a json file
# to be adapted/changed by users depending on their environments
from driver_tools_lib import get_metric_values, get_metric_values_observations, get_mod_mem_json

# ENSO_metrics functions
from EnsoPlots.EnsoPlotTemplate import plot_portraitplot
from EnsoPlots.EnsoPlotToolsLib import sort_metrics, sort_models


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
# metric collections to plot
list_metric_collections = ["ENSO_perf", "ENSO_tel", "ENSO_proc"]
# CMIP experiment
experiment = "historical"
# project to use, here both CMIP5 and CMIP6 models will be used
list_projects = ["CMIP6", "CMIP5"]
# list of additional observations
# the reading part is very 'ad hoc', do not change the obs!
list_observations = ["20CRv2", "NCEP2", "ERA-Interim"]
# True to use the set of metric in the BAMS paper
# More metric have been computed and tested but not kept
reduced_set = True  # False  #
# True to use the first available member only
# If set to False, all members will be used and the metric values computed for all members of each model will be
# averaged
first_member = True  # False  #
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
figure_name = "portraitplot_" + str(len(list_metric_collections)) + "metric_collections_" + version
if len(list_projects) == 1:
    figure_name += "_" + str(list_projects[0])
else:
    figure_name += "_" + str(len(list_projects)) + "cmip"
if first_member is True:
    figure_name += "_first_member"
else:
    figure_name += "_members_averaged"
if reduced_set is False:
    figure_name += "_all_metrics"
figure_name = OSpath__join(path_out, figure_name)
# Metric collection names on the figure
metric_collection_names_for_plot = {"ENSO_perf": "Performance", "ENSO_proc": "Processes", "ENSO_tel": "Telecon."}
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
tab_all, tab_all_act, x_names = list(), list(), list()
for mc in list_metric_collections:
    dict_met = dict()
    for proj in list_projects:
        dict1 = get_metric_values(proj, mc, dict_json, model_by_proj, reduced_set=reduced_set, portraitplot=True)
        # save in common dictionary
        for mod in dict1.keys():
            dict_met[mod] = dict1[mod]
        del dict1
    # models and metrics
    tmp_models = sorted([str(mod) for mod in list(dict_met.keys())], key=lambda v: v.upper())
    my_metrics = list()
    for mod in tmp_models:
        try: list(dict_met[mod].keys())
        except: pass
        else: my_metrics += list(dict_met[mod].keys())
    my_metrics = sort_metrics(sorted(list(set(my_metrics)), key=lambda v: v.upper()))
    my_models = sort_models(tmp_models)
    # read other observational datasets compared to the reference
    dict_ref_met = get_metric_values_observations(dict_json["obs2obs"][mc], list_observations, my_metrics, mc)
    # number of line to add to the array (CMIP mean, reference, other observational datasets,...)
    plus = 3 + len(list_observations)
    # fill array
    tab = NUMPYma__zeros((len(my_models) + plus, len(my_metrics)))
    for ii, mod in enumerate(my_models):
        for jj, met in enumerate(my_metrics):
            if met not in dict_met[mod].keys() or dict_met[mod][met] is None:
                tab[ii + plus, jj] = 1e20
            else:
                tab[ii + plus, jj] = dict_met[mod][met]
    tab = NUMPYma__masked_invalid(tab)
    tab = NUMPYmasked_where(tab == 1e20, tab)
    # add values to the array (CMIP mean, reference, other observational datasets,...)
    for jj, met in enumerate(my_metrics):
        tmp = tab[plus:, jj].compressed()
        mea = float(NUMPYmean(tmp))
        std = float(NUMPYstd(tmp))
        del tmp
        for ii, dd in enumerate(list_observations + ["reference"] + list_projects):
            if dd in list_observations:
                val = dict_ref_met[dd][met]
            elif dd in list_projects:
                tmp = [tab[kk + plus, jj] for kk, mod in enumerate(my_models) if mod in list(model_by_proj[dd].keys())]
                tmp = NUMPYma__masked_invalid(NUMPYma__array(tmp))
                tmp = NUMPYmasked_where(tmp == 1e20, tmp).compressed()
                val = float(NUMPYmean(tmp))
                del tmp
            else:
                val = 0
            tab[ii, jj] = val
            del val
        # normalize
        tab[:, jj] = (tab[:, jj] - mea) / std
        del mea, std
    tab = NUMPYma__masked_invalid(tab)
    tab = NUMPYmasked_where(tab > 1e3, tab)
    tab_all.append(tab)
    x_names.append(my_metrics)
    if mc == list_metric_collections[0]:
        y_names = ["("+dd+")" for dd in list_observations] + ["(reference)"] + list_projects +\
                  ["* " + mod if mod in list(model_by_proj["CMIP6"].keys()) else mod for mod in my_models]
    del dict_met, dict_ref_met, my_metrics, my_models, plus, tab, tmp_models


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#if ' ':
    numbering = [ii+") " for ii in list(string.ascii_lowercase)]
    # plot
    title = [numbering[ii] + metric_collection_names_for_plot[mc]
             for ii, mc in enumerate(list_metric_collections)]
    text = "* = CMIP6\nmodel"
    levels = list(range(-2, 3))
    plot_portraitplot(tab_all, figure_name, xticklabel=x_names, yticklabel=y_names, title=title, my_text=text,
                      levels=levels, cfram=True, chigh=True)
    del levels, numbering, text, title
