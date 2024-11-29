# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots of the correlation inter recipes or inter models
#      FIG. 6 in Planton et al. 2020: Evaluating climate models with the CLIVAR 2020 ENSO recipes package. BAMS
# It uses the first available member of each model or all members of each model and averages them
# Updated json files (needed to create this plot) can be downloaded from the page "Summary statistics in Interactive
# Portrait Plots" at https://cmec.llnl.gov/results/enso/
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#

from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros
from os.path import join as OSpath__join
from scipy.stats import linregress as SCIPYstats__linregress

# set of functions to find cmip/obs files and save a json file
# to be adapted/changed by users depending on their environments
from driver_tools_lib import get_metric_values, get_mod_mem_json

# ENSO_metrics functions
from EnsoPlots.EnsoPlotTemplate import plot_metrics_correlations
from EnsoPlots.EnsoPlotToolsLib import sort_metrics


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
# True to use the first available member only
# If set to False, all members will be used and the metric values computed for all members of each model will be
# averaged
first_member = True  # False  #
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
# figure name
path_out = ""
figure_name = "metrics_correlations_" + str(len(list_metric_collections)) + "metric_collections_" + version_mod
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
# ---------------------------------------------------#


# ---------------------------------------------------#
# Functions
# ---------------------------------------------------#
def compute_correlation(tab_in):
    """
    Computes correlations

    Input:
    -----
    :param tab_in: `cdms2` variable
        A `cdms2` variable containing the data to be analysed.

    Outputs:
    -------
    :return rval: `cdms2` variable
        A `cdms2` variable containing the correlation coefficients between the values along the first axis.
    :return pval: `cdms2` variable
        A `cdms2` variable containing the two-sided p-value for a hypothesis test whose null hypothesis is that the
        slope is zero, using Wald Test with t-distribution of the test statistic. I.e., if the absolute value of the
        correlation is smaller than the p-value, it means that the correlation is not significant.
    """
    pval = NUMPYma__zeros((len(tab_in), len(tab_in)))
    rval = NUMPYma__zeros((len(tab_in), len(tab_in)))
    for ii in range(len(tab_in)):
        tmp1 = tab_in[ii]
        for jj in range(len(tab_in)):
            tmp2 = tab[jj]
            tmpf1 = NUMPYmasked_where(tmp2.mask, tmp1)
            tmpf2 = NUMPYmasked_where(tmpf1.mask, tmp2)
            tmpf1 = tmpf1.flatten().compressed()
            tmpf2 = tmpf2.flatten().compressed()
            slope, intercept, r_value, p_value, std_err = SCIPYstats__linregress(tmpf1, tmpf2)
            pval[ii, jj] = float(p_value)
            rval[ii, jj] = float(r_value)
            del intercept, p_value, r_value, slope, std_err, tmp2, tmpf1, tmpf2
        del tmp1
    return rval, pval
# ---------------------------------------------------#


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# get members by model by project from json file
# only recipes from models/members chosen here will be used
# all recipes from models/members chosen here will be used (ensures that if a model/member is not available for one or
# several metric collections, the corresponding line will still be created in the portraitplot)
model_by_proj = get_mod_mem_json(list_projects, list_metric_collections, dict_json, first_only=first_member)
# read json file
dict_met = dict()
for proj in list_projects:
    for mc in list_metric_collections:
        dict1 = get_metric_values(proj, mc, dict_json, model_by_proj, reduced_set=reduced_set)
        # save in common dictionary
        for mod in list(dict1.keys()):
            try:    dict_met[mod]
            except: dict_met[mod] = dict1[mod]
            else:
                for met in list(dict1[mod].keys()):
                    dict_met[mod][met] = dict1[mod][met]
        del dict1


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
if ' ':
    list_metrics = list()
    for k1 in list(dict_met.keys()):
        list_metrics += list(dict_met[k1].keys())
    list_metrics = sort_metrics(list(set(list_metrics)))
    list_models = list(dict_met.keys())
    # fill 2D-array with metric values
    tab = NUMPYma__zeros((len(list_metrics), len(list_models)))
    for ii, met in enumerate(list_metrics):
        for jj, mod in enumerate(list_models):
            if met not in list(dict_met[mod].keys()) or dict_met[mod][met] is None:
                tab[ii, jj] = 1e20
            else:
                tab[ii, jj] = dict_met[mod][met]
    tab = NUMPYma__masked_invalid(tab)
    tab = NUMPYmasked_where(tab == 1e20, tab)
    # compute inter model correlations
    rval, pval = compute_correlation(tab)
    # plot recipes correlations
    plot_metrics_correlations(rval, figure_name, list_metrics, tab_pval=pval, cfram=True, chigh=True)
