# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots of the correlation inter metrics or inter models
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
import cmocean
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from numpy import arange as NUMPYarange
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros
from os.path import join as OSpath__join
from scipy.stats import linregress as SCIPYstats__linregress
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoPlotLib import plot_param


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
experiment = "historical"  # "piControl" #
list_project = ["cmip5", "cmip6"]
my_project = ["CMIP"]
big_ensemble = True
reduced_set = True  # False  #

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_10_report"
path_in = OSpath__join(path_main, "Data_grouped")
path_out = OSpath__join(path_main, "Plots_v5")

expe = "hist" if experiment == "historical" else "pi"


# ---------------------------------------------------#
# Functions
# ---------------------------------------------------#
def add_suffix(list_met):
    list_metrics1 = deepcopy(list_met)
    list_metrics2 = deepcopy(list_met)
    for met in list_metrics2:
        if "Map" in met:
            while met in list_metrics1:
                list_metrics1.remove(met)
            for suffix in ["Corr", "Rmse"]:
                list_metrics1.append(met + suffix)
    return sorted(list_metrics1, key=lambda v: v.upper())


def common_save(dict_in, dict_out={}):
    for mod in dict_in.keys():
        try:
            dict_out[mod]
        except:
            dict_out[mod] = dict_in[mod]
        else:
            for met in dict_in[mod].keys():
                try:
                    dict_out[mod][met]
                except:
                    dict_out[mod][met] = dict_in[mod][met]
                else:
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                        ": metric already in output",
                        str().ljust(5) + "metric '" + str(met) + "' for model '" + str(mod) +
                        "' is already in output dictionary",
                        str().ljust(10) + "in  input = " + str(dict_in[mod][met]),
                        str().ljust(10) + "in output = " + str(dict_out[mod][met]),
                    ]
                    EnsoErrorsWarnings.MyError(list_strings)
    return dict_out


def compute_correlation(tab):
    """
    Performs hierarchical/agglomerative clustering on the given array.

    Input:
    -----
    :param tab: `cdms2` variable
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
    pval = NUMPYma__zeros((len(tab), len(tab)))
    rval = NUMPYma__zeros((len(tab), len(tab)))
    for ii in range(len(tab)):
        tmp1 = tab[ii]
        for jj in range(len(tab)):
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


def get_reference(metric_collection, metric):
    if metric_collection == "ENSO_tel" and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    return plot_param(metric_collection, my_met)['metric_reference']



def get_ref(metric):
    for mc in metric_collection:
        list_met = sorted(defCollection(mc)['metrics_list'].keys(), key=lambda v: v.upper())
        if mc == "ENSO_tel" and "Map" in metric:
            my_met = metric.replace("Corr", "").replace("Rmse", "")
        else:
            my_met = deepcopy(metric)
        if my_met in list_met:
            break
    return plot_param(mc, my_met)['metric_reference']


def plot_correlation(tab_rval, name_plot, xy_names, tab_pval=None, write_corr=False, title='correlations'):
    """
    Plots the correlations matrix.

    Inputs:
    ------
    :param tab_rval: `cdms2` variable
        A `cdms2` variable containing the correlation coefficients.
    :param name_plot: string
        Path to and name of the output plot.
    :param xy_names: list of string
        List of the names to put as x and y ticks.
    **Optional arguments:**
    :param tab_pval: `cdms2` variable
        A `cdms2` variable containing the p-value of the correlation coefficients. It is used to mask the correlations
        that are not significant.
        Default is None (all correlation coefficients are plotted).
    :param title: boolean
        True to write the correlation value in each box
    :param title: string
        Title of the plot.

    Output:
    ------
    """
    if tab_pval is not None:
        mask1 = NUMPYmasked_where(tab_rval < tab_pval, tab_rval).mask.astype('f')
        mask2 = NUMPYmasked_where(tab_rval > -tab_pval, tab_rval).mask.astype('f')
        mask = mask1 + mask2
        tab_plot = NUMPYmasked_where(mask == 2, tab_rval)
    else:
        tab_plot = deepcopy(tab_rval)
    fig, ax = plt.subplots(figsize=(0.5 * len(tab_rval), 0.5 * len(tab_rval)))
    # shading & colorbar
    # cs = ax.contourf(xx, xx, tab_plot, levels=levels, extend="both", cmap="cmo.balance")
    # cbar = plt.colorbar(cs, orientation="horizontal", ticks=labels, pad=0.35, extend="both")
    levels = MaxNLocator(nbins=20).tick_values(-1, 1)
    cmap = plt.get_cmap("cmo.balance")
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cs = plt.pcolormesh(tab_plot, cmap=cmap, norm=norm)
    # title
    plt.title(title, fontsize=30, y=1.01, loc='center')
    # x axis
    xticks = [ii + 0.5 for ii in range(len(tab_plot))]
    xlabel = [met for met in xy_names]
    plt.xticks(xticks, xlabel, fontsize=15)
    plt.xticks(rotation=90)
    # y axis
    plt.yticks(xticks, xlabel, fontsize=15)
    # text
    if write_corr is True:
        for ii in range(len(tab_plot)):
            for jj in range(len(tab_plot)):
                plt.text(ii + 0.5, jj + 0.5, str(round(tab_rval[ii, jj], 1)), fontsize=10,
                         horizontalalignment='center', verticalalignment='center')
    if tab_pval is not None:
        for ii in range(len(tab_plot)):
            nbr1 = str(sum([1 for jj in range(len(tab_plot[ii])) if tab_plot[ii][jj] < 0]))
            nbr2 = str(int(len(tab_plot) - 1 - sum(tab_plot[ii].mask.astype('f')))).zfill(2)
            plt.text(len(tab_plot) + 0.5, ii + 0.5, nbr1, fontsize=18, horizontalalignment='center',
                     verticalalignment='center')
            plt.text(len(tab_plot) + 0.5, -0.5, "nbr < 0", fontsize=15, horizontalalignment='center',
                     verticalalignment='top', rotation=90)
            plt.text(len(tab_plot) + 1 + 0.5, ii + 0.5, nbr2, fontsize=18, horizontalalignment='center',
                     verticalalignment='center')
            plt.text(len(tab_plot) + 1 + 0.5, -0.5, "nbr significant", fontsize=15, horizontalalignment='center',
                     verticalalignment='top', rotation=90)
    # color bar
    levels = [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)]
    x2 = ax.get_position().x1
    y1 = ax.get_position().y0
    y2 = ax.get_position().y1
    if tab_pval is not None:
        cax = plt.axes([x2 + 0.09, y1, 0.02, y2 - y1])
    else:
        cax = plt.axes([x2 + 0.02, y1, 0.02, y2 - y1])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
    cbar.ax.tick_params(labelsize=18)
    # save fig
    plt.savefig(name_plot, bbox_inches='tight')
    plt.close()
    return


def read_data(project, metric_collection):
    #lpath = OSpath__join(path_in, project + "/" + experiment + "/" + metric_collection)
    lname = project + "_" + experiment + "_" + metric_collection + "_v2019????_modified.json"
    filename_js = list(GLOBiglob(OSpath__join(path_in, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    return data["RESULTS"]["model"]


def remove_metrics(list_met, metric_collection):
    list_met1 = deepcopy(list_met)
    if reduced_set is True:
        if metric_collection == "ENSO_perf":
            # to_remove = ['BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse', 'NinaSstDur_1',
            #              'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
            #              'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
            #              'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
            #              'NinoSstTsRmse_2']
            to_remove = ['BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse', 'NinaSstDur_1',
                         'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                         'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                         'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                         'NinoSstTsRmse_2', "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
        elif metric_collection == "ENSO_proc":
            # to_remove = ['EnsoAmpl', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf',
            #              'EnsoFbTauxSsh']
            to_remove = ['BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality', 'EnsoSstLonRmse',
                         'EnsoSstSkew', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf',
                         'EnsoFbSstSwr']
        else:
            to_remove = ['EnsoAmpl', 'EnsoPrMapStd', 'EnsoSlpMapCorr', 'EnsoSlpMapRmse', 'EnsoSlpMapStd',
                         'EnsoSstMapStd', 'EnsoSstLonRmse',
                         'NinaPrMap_1Corr', 'NinaPrMap_1Rmse', 'NinaPrMap_1Std',
                         'NinaPrMap_2Corr', 'NinaPrMap_2Rmse', 'NinaPrMap_2Std',
                         'NinaSlpMap_1Corr', 'NinaSlpMap_1Rmse', 'NinaSlpMap_1Std',
                         'NinaSlpMap_2Corr', 'NinaSlpMap_2Rmse', 'NinaSlpMap_2Std',
                         'NinaSstLonRmse_1', 'NinaSstLonRmse_2',
                         'NinaSstMap_1Corr', 'NinaSstMap_1Rmse', 'NinaSstMap_1Std',
                         'NinaSstMap_2Corr', 'NinaSstMap_2Rmse', 'NinaSstMap_2Std',
                         'NinoPrMap_1Corr', 'NinoPrMap_1Rmse', 'NinoPrMap_1Std',
                         'NinoPrMap_2Corr', 'NinoPrMap_2Rmse', 'NinoPrMap_2Std',
                         'NinoSlpMap_1Corr', 'NinoSlpMap_1Rmse', 'NinoSlpMap_1Std',
                         'NinoSlpMap_2Corr', 'NinoSlpMap_2Rmse', 'NinoSlpMap_2Std',
                         'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                         'NinoSstMap_1Corr', 'NinoSstMap_1Rmse', 'NinoSstMap_1Std',
                         'NinoSstMap_2Corr', 'NinoSstMap_2Rmse', 'NinoSstMap_2Std']
    else:
        if metric_collection == "ENSO_perf":
            to_remove = []
        elif metric_collection == "ENSO_proc":
            to_remove = ['BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality', 'EnsoSstLonRmse',
                         'EnsoSstSkew']
        else:
            to_remove = ['EnsoAmpl', 'EnsoSstLonRmse', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinoSstLonRmse_1',
                         'NinoSstLonRmse_2']
    for met in to_remove:
        while met in list_met1:
            list_met1.remove(met)
    # !!!!! temporary: start !!!!!
    # # ssh metrics are not computed yet (ask jiwoo)
    # list_met2 = deepcopy(list_met1)
    # for met in list_met2:
    #     if "Ssh" in met:
    #         while met in list_met1:
    #             list_met1.remove(met)
    # del list_met2
    # # slp metrics are wrong (error in observation?)
    # list_met2 = deepcopy(list_met1)
    # for met in list_met2:
    #     if "Slp" in met:
    #         while met in list_met1:
    #             list_met1.remove(met)
    # del list_met2
    # !!!!! temporary: end !!!!!
    return list_met1


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
dict_met = dict()
for proj in list_project:
    if big_ensemble is not True or (big_ensemble is True and proj == list_project[0]):
        dict_mc = dict()
    for mc in metric_collection:
        # read json
        data_json = read_data(proj, mc)
        list_models = sorted(data_json.keys(), key=lambda v: v.upper())
        # read metrics
        list_models = sorted(data_json.keys(), key=lambda v: v.upper())
        dict1 = dict()
        for mod in list_models:
            data_mod = data_json[mod]["value"]
            list_metrics = sorted(data_mod.keys(), key=lambda v: v.upper())
            list_metrics = remove_metrics(list_metrics, mc)
            dict2 = dict()
            for met in list_metrics:
                try:
                    ref = get_reference(mc, met)
                except:
                    ref = "Tropflux"
                data_met = data_mod[met]["metric"]
                list_ref = sorted(data_met.keys(), key=lambda v: v.upper())
                if ref in list_ref:
                    my_ref = deepcopy(ref)
                else:
                    if "ERA-Interim" in list_ref:
                        my_ref = "ERA-Interim"
                    elif "ERA-Interim_ERA-Interim" in list_ref:
                        my_ref = "ERA-Interim_ERA-Interim"
                    else:
                        list_strings = [
                            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                            ": cannot fing a proper reference",
                            str().ljust(5) + "project '" + proj + "', MC '" + mc + "', model '" + mod + "', metric '" +
                            met + "'",
                            str().ljust(10) + "input references = " + str(list_ref)]
                        EnsoErrorsWarnings.MyError(list_strings)
                    # list_strings = [
                    #     "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                    #     ": reference changed",
                    #     str().ljust(5) + "project '" + proj + "', MC '" + mc + "', model '" + mod + "', metric '" +
                    #     met + "'",
                    #     str().ljust(10) + "reference set to '" + my_ref +"' instead of '" + ref + "'"]
                    # EnsoErrorsWarnings.MyWarning(list_strings)
                if data_met[my_ref]["value"] is None:
                    dict2[met] = 1e20
                else:
                    dict2[met] = data_met[my_ref]["value"]
                del data_met, my_ref, list_ref, ref
            dict1[mod] = dict2
            del data_mod, dict2
        # save in common dictionary
        dict_mc = common_save(dict1, dict_out=dict_mc)
        del data_json, dict1, list_metrics, list_models
    if big_ensemble is not True:
        dict_met[proj] = dict_mc
        del dict_mc
if big_ensemble is True:
    dict_met = deepcopy(dict_mc)
    del dict_mc


# # show dictionary levels
# lev1 = sorted(dict_met.keys(), key=lambda v: v.upper())
# print "level1 (" + str(len(lev1)) + ") = " + str(lev1)
# print ""
# # check metrics
# if big_ensemble is True:
#     list1 = list(set([len(dict_met[key1].keys()) for key1 in lev1]))
#     for key1 in lev1:
#         if len(dict_met[key1].keys()) == max(list1):
#             list_metrics = sorted(dict_met[key1].keys(), key=lambda v: v.upper())
#             pass
#     # !!!!! temporary: start !!!!!
#     # some models are not used in all metric collection (ask jiwoo)
#     for key2 in list_metrics:
#         for key1 in lev1:
#             if key2 not in dict_met[key1].keys():
#                 dict_met[key1][key2] = dict((ref, 1e20) for ref in dict_met["ACCESS1-0"][key2].keys())
#     list1 = list(set([len(dict_met[key1].keys()) for key1 in lev1]))
#     # !!!!! temporary: end !!!!!
# else:
#     list1 = list(set([len(dict_met[key1][key2].keys()) for key1 in lev1 for key2 in dict_met[key1].keys()]))
#     for key1 in lev1:
#         for key2 in dict_met[key1].keys():
#             if len(dict_met[key1][key2].keys()) == max(list1):
#                 list_metrics = sorted(dict_met[key1][key2].keys(), key=lambda v: v.upper())
#                 pass
# if len(list1) != 1:
#     for key1 in lev1:
#         lev2 = sorted(dict_met[key1].keys(), key=lambda v: v.upper())
#         if big_ensemble is True:
#             if len(lev2) != len(list_metrics):
#                 print key1.rjust(15) + " (" + str(len(lev2)).zfill(2) + "), missing = " +\
#                       str(list(set(list_metrics) - set(lev2)))
#         else:
#             for key2 in lev2:
#                 lev3 = sorted(dict_met[key1][key2].keys(), key=lambda v: v.upper())
#                 if len(lev3) != len(list_metrics):
#                     print key2.rjust(15) + " (" + str(len(lev3)).zfill(2) + "), missing = " +\
#                           str(list(set(list_metrics) - set(lev3)))
#                 del lev3
#         del lev2
#     list_strings = [
#         "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
#         ": not the same number of metrics for all models",
#         str().ljust(5) + "metric number = " + str(list1),
#     ]
#     EnsoErrorsWarnings.MyError(list_strings)
# else:
#     print str().ljust(5) + "metrics (" + str(len(list_metrics)) + ") = " + str(list_metrics)
# del list1
# for ii in range(3): print ""


# # set reference observation
# dict_out = dict()
# for met in list_metrics:
#     ref = get_ref(met)
#     for key1 in lev1:
#         if big_ensemble is True:
#             try: tmp = dict_met[key1][met][ref]
#             except:
#                 ref2 = sorted(dict_met[key1][met].keys(), key=lambda v: v.upper())[0]
#                 list_strings = [
#                     "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
#                     ": reference (" + ref + ") not available",
#                     str().ljust(5) + key1.rjust(15) + ", " + str(met),
#                     str().ljust(5) + "another reference is used (" + ref2 + ")",
#                 ]
#                 # EnsoErrorsWarnings.MyWarning(list_strings)
#                 tmp = dict_met[key1][met][ref2]
#             try: dict_out[key1]
#             except: dict_out[key1] = {met: tmp}
#             else: dict_out[key1][met] = tmp
#             del tmp
#         else:
#             lev2 = sorted(dict_met[key1].keys(), key=lambda v: v.upper())
#             for key2 in lev2:
#                 try: tmp = dict_met[key1][key2][met][ref]
#                 except:
#                     ref2 = sorted(dict_met[key1][key2].keys(), key=lambda v: v.upper())[0]
#                     list_strings = [
#                         "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
#                         ": reference (" + ref + ") not available",
#                         str().ljust(5) + key1.rjust(15) + ", " + str(met),
#                         str().ljust(5) + "another reference is used (" + ref2 + ")",
#                     ]
#                     # EnsoErrorsWarnings.MyWarning(list_strings)
#                     tmp = dict_met[key1][key2][ref2]
#                 try: dict_out[key1]
#                 except: dict_out[key1] = {key2: {met: tmp}}
#                 else:
#                     try: dict_out[key1][key2]
#                     except: dict_out[key1][key2] = {met: tmp}
#                     else: dict_out[key1][key2][met] = tmp
#                 del tmp
#             del lev2
#     del ref

# if ' ':
#     met = "EnsodSstOce_2"
#     tmp = [dict_out[mod][met] for jj, mod in enumerate(lev1) if dict_out[mod][met] < 1e10]
#     print met + " ranges from "+str(round(min(tmp), 1))+" to "+str(round(max(tmp), 1))
#     myset = [mod for jj, mod in enumerate(lev1) if dict_out[mod][met] < 1.]
#     print "model weakly driven by ocean: "+str(myset)
#     met = "EnsoFbSstThf"
#     for ii in NUMPYarange(0, -15.01, -0.1):
#         tmp = [mod for jj, mod in enumerate(lev1) if (dict_out[mod][met] > ii and dict_out[mod][met] < 1e10)]
#         allmod = True
#         for mod in myset:
#             if mod not in tmp:
#                 allmod = False
#         if allmod is True:
#             break
#     print "model weakly driven by fluxes (fb>"+str(round(ii, 1))+"): " + str(tmp)
#     stop

# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
dict_out = deepcopy(dict_met)
lev1 = sorted(dict_out.keys(), key=lambda v: v.upper())
# correlation
if ' ':
    if big_ensemble is True:
        dict_out = deepcopy(dict_met)
        lev1 = sorted(dict_out.keys(), key=lambda v: v.upper())
        list_metrics = sorted(dict_out[lev1[0]].keys(), key=lambda v: v.upper())
        tab = NUMPYma__zeros((len(list_metrics), len(lev1)))
        for ii, met in enumerate(list_metrics):
            for jj, mod in enumerate(lev1):
                try: dict_out[mod][met]
                except: tab[ii, jj] = 1e20
                else: tab[ii, jj] = dict_out[mod][met]
                # tab[ii, jj] = 1 - dict_out[mod][met] if "Corr" in met else dict_out[mod][met]
        tab = NUMPYmasked_where(tab == 1e20, tab)
        # compute inter model correlations
        rval, pval = compute_correlation(tab)
        name_plot = OSpath__join(path_out, "correlations_inter_metrics_" + str(len(list_metrics)).zfill(2) +
                                 "metrics_" + str(len(lev1)).zfill(2) + "models")
        title = "inter metric correlations"
        plot_correlation(rval, name_plot, list_metrics, tab_pval=pval, write_corr=True, title=title)
        # # compute inter model correlations
        # tab = tab.reorder('10')
        # rval, pval = compute_correlation(tab)
        # name_plot = OSpath__join(path_out, "correlations_inter_models_" + str(len(list_metrics)).zfill(2) +
        #                          "metrics_" + str(len(lev1)).zfill(2) + "models")
        # title = "inter model correlations"
        # plot_correlation(rval, name_plot, list_metrics, tab_pval=pval, write_corr=True, title=title)

# if ' ':
#     lev1.remove("BNU-ESM")
#     lev1.remove("GFDL-CM2p1")
#     lev1.remove("GFDL-CM4")
#     lev1.remove("HadCM3")
#     lev1.remove("HadGEM2-AO")
#     tmp1 = [dict_out[mod]["EnsoSeasonality"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.25 / Q2 ~ 0.34
#     tmp2 = [dict_out[mod]["BiasSstLonRmse"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.64 / Q2 ~ 1.00
#     tmp3 = [dict_out[mod]["BiasPrLonRmse"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.73 / Q2 ~ 1.00
#     tmp4 = [dict_out[mod]["SeasonalSstLonRmse"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.18 / Q2 ~ 0.21
#     tmp5 = [dict_out[mod]["SeasonalPrLonRmse"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.44 / Q2 ~ 0.63
#     tmp6 = [dict_out[mod]["EnsoFbSstTaux"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.39 / Q2 ~ 0.47
#     tmp7 = [dict_out[mod]["EnsoFbSstThf"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.23 / Q2 ~ 0.48
#     tmp8 = [dict_out[mod]["EnsoFbSstSwr"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.54 / Q2 ~ 1.18
#     tmp9 = [dict_out[mod]["EnsodSstOce_2"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.17 / Q2 ~ 0.33
#     tmp10 = [dict_out[mod]["EnsoSstSkew"] for jj, mod in enumerate(lev1)]  # Q1 ~ 0.52 / Q2 ~ 0.82
#     selections = [
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp2[ii] < 1.0],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp3[ii] < 1.0],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp4[ii] < 0.25],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp5[ii] < SCIPYstats__scoreatpercentile(tmp5, 50)],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp6[ii] < 0.5],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp7[ii] < 0.5],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp1[ii] < 0.25 and tmp8[ii] < 1.0],
#         [lev1[ii] for ii in range(len(tmp1)) if tmp6[ii] < SCIPYstats__scoreatpercentile(tmp6, 25)
#          and tmp7[ii] < SCIPYstats__scoreatpercentile(tmp7, 25)],
#         [lev1[ii] for ii in range(len(tmp1))
#          if tmp6[ii] < SCIPYstats__scoreatpercentile(tmp6, 50) and tmp9[ii] < SCIPYstats__scoreatpercentile(tmp9, 25)],
#         [lev1[ii] for ii in range(len(tmp1))
#          if tmp7[ii] < SCIPYstats__scoreatpercentile(tmp7, 50) and tmp9[ii] < SCIPYstats__scoreatpercentile(tmp9, 25)],
#         [lev1[ii] for ii in range(len(tmp1)) if
#          tmp1[ii] < SCIPYstats__scoreatpercentile(tmp1, 50) and tmp10[ii] < SCIPYstats__scoreatpercentile(tmp10, 50)],
#     ]
#     for ii, tmp in enumerate(selections):
#         print "selection" + str(ii+1).zfill(2) + " (" + str(len(tmp)).zfill(2) + "): "+str(tmp)

