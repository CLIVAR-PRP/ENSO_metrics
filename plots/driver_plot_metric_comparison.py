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
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from numpy import array as NUMPYarray
from numpy import mean as NUMPYmean
from numpy import moveaxis as NUMPYmoveaxis
from numpy import sort as NUMPYsort
from numpy.ma import masked_where as NUMPYma__masked_where
from numpy.random import randint as NUMPYrandom__randint
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoPlotLib import plot_param
from EnsoPlotToolsLib import minmax_plot


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
experiment = "historical"  # "piControl" #
member = "r1i1p1"
list_project = ["cmip6", "cmip5"]
my_project = ["select8", "CMIP"]
big_ensemble = False  # True
reduced_set = True  # False  #
dict_selection = {
    #
    # EnsoSeasonality within 25% of obs &
    #
    # BiasSstLonRmse less than 1°C
    "select1": ["ACCESS1-0", "BCC-ESM1", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "GISS-E2-R", "GISS-E2-R-CC",
                "MIROC4h", "NorESM1-M"],
    # BiasPrLonRmse less than 1 mm/day
    "select2": ["BCC-CSM1-1", "BCC-ESM1", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2",
                "NorESM1-M", "NorESM1-ME"],
    # SeasonalSstLonRmse less than 0.25°C
    "select3": ["BCC-CSM1-1", "BCC-ESM1", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2",
                "NorESM1-M", "NorESM1-ME"],
    # SeasonalPrLonRmse less than Q2 (median)
    "select4": ["CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R", "GISS-E2-H-CC",
                "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstTaux within 50% of obs
    "select5": ["CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R", "GISS-E2-R-CC",
                "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstThf within 50% of obs
    "select6": ["ACCESS1-0", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R",
                "GISS-E2-R-CC", "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstSwr within 100% of obs
    "select7": ["ACCESS1-0", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R",
                "GISS-E2-R-CC", "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    #
    # EnsoFbSstTaux & EnsoFbSstThf within Q1 (approximately)
    #
    "select8": ["CCSM4", "CESM1-FASTCHEM", "CNRM-CM5", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
                "GISS-E2-R", "GISS-E2-R-CC", "MIROC6"],
    #
    # EnsoFbSstTaux within Q2 & EnsodSstOce_2 within Q1
    #
    "select9": ["CESM1-BGC", "CESM2", "CNRM-CM5", "EC-Earth3-Veg", "FGOALS-g2", "GISS-E2-R-CC", "MIROC4h", "NorESM1-ME",
                "SAM0-UNICON"],
    #
    # EnsoFbSstThf within Q2 & EnsodSstOce_2 within Q1
    #
    "select10": ["ACCESS1-0", "CESM1-BGC", "CESM2", "CNRM-CM5", "EC-Earth3-Veg", "FGOALS-g2", "GFDL-CM3", "GISS-E2-H",
                 "GISS-E2-R-CC", "MIROC4h", "NorESM1-ME"],
    #
    # EnsoSeasonality within Q2 & EnsoSstSkew within Q2
    #
    "select11": ["BCC-CSM1-1-M", "CCSM4", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2-WACCM", "CNRM-CM5",
                 "GFDL-ESM2G", "GISS-E2-R", "MIROC4h", "MIROC6", "NorESM1-ME"],
}

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


def bootstrap(tab, num_samples=100000, alpha=0.05, nech=None, statistic=NUMPYmean):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    n = len(tab)
    if nech is None:
        nech = deepcopy(n)
    idx = NUMPYrandom__randint(0, n, (num_samples, nech))
    samples = tab[idx]
    stat = NUMPYsort(statistic(samples, 1))
    return [stat[int((alpha/2.0)*num_samples)], stat[int((1-alpha/2.0)*num_samples)]]


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


def plot_metrics(tab_val, name_plot, title="", x_names=None, y_name="", colors=None, tab_bst=None, legend=None):
    fig, ax = plt.subplots(figsize=(0.5 * len(tab_val[0]), 4))
    # title
    plt.title(title, fontsize=15, y=1.01, loc='center')
    # x axis
    axis = list(range(len(tab_val[0])))
    ax.set_xticks(axis)
    if isinstance(x_names, list):
        ax.set_xticklabels(x_names)
    else:
        ax.set_xticklabels(axis)
    ax.set_xlim([min(axis) - 0.5, max(axis) + 0.5])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    # y axis
    tmp = [tab_val.min(), tab_val.max()]
    if tab_bst is not None:
        tmp += [tab_bst.min(), tab_bst.max()]
    ytick = minmax_plot(tmp)
    ax.set_yticks(ytick)
    ax.set_yticklabels(ytick)
    # ax.set_ylim([min(ytick), max(ytick)])
    ax.set_ylim([min(ytick), 2.5])
    ax.set_ylabel(y_name, fontsize=15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    # plot marker
    for ii in range(len(tab_val)):
        if isinstance(colors, list):
            col = colors[len(colors) - 1 - ii]
        else:
            col = "k"
        ind = len(tab_val) - 1 - ii
        ax.scatter(axis, list(tab_val[ind]), s=60, c=col, marker="D")
        if tab_bst is not None:
            for jj in range(len(tab_bst[ind])):
                if tab_bst[ind][jj][0] > 0 and tab_bst[ind][jj][1] > 0:
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tab_bst[ind][jj][0], tab_bst[ind][jj][0]], c=col, lw=1))
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tab_bst[ind][jj][1], tab_bst[ind][jj][1]], c=col, lw=1))
                    ax.add_line(Line2D([jj, jj], [tab_bst[ind][jj][0], tab_bst[ind][jj][1]], c=col, lw=1))
        del col
    # text
    if legend is not None:
        x1, x2 = ax.get_xlim()
        dx = (x2 - x1) / 100.
        y1, y2 = ax.get_ylim()
        dy = (y2 - y1) / 100.
        for ii in range(len(legend)):

            if isinstance(colors, list):
                col = colors[len(colors) - 1 - ii]
            else:
                col = "k"
            font = {'color': col, 'weight': 'normal', 'size': 15}
            ax.text(x2 - 2 * dx, y2 - (ii + 1) * 7 * dy, legend[len(legend) - 1 - ii], horizontalalignment="right",
                    verticalalignment="center", fontdict=font)
            del col, font
        ax.add_line(Line2D([25 * dx - 0.3, 25 * dx + 0.3], [2.3, 2.3], c="k", lw=1))
        ax.add_line(Line2D([25 * dx - 0.3, 25 * dx + 0.3], [1.7, 1.7], c="k", lw=1))
        ax.add_line(Line2D([25 * dx, 25 * dx], [1.7, 2.3], c="k", lw=1))
        font = {'weight': 'normal', 'size': 15}
        ax.text(26 * dx, 2., "95% confidence interval\n(Monte Carlo Method)", horizontalalignment="left",
                verticalalignment="center", fontdict=font, transform=ax.transData)
    # save fig
    plt.savefig(name_plot, bbox_inches='tight')
    plt.close()
    return


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
dict_met = dict()
for proj in list_project:
    if big_ensemble is not True or (big_ensemble is True and proj == list_project[0]):
        dict_mc = dict()
    for mc in metric_collection:
        # get metrics list
        list_metrics = sorted(defCollection(mc)['metrics_list'].keys(), key=lambda v: v.upper())
        if reduced_set is True:
            if mc == "ENSO_perf":
                # to_remove = ['BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse', 'NinaSstDur_1',
                #              'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                #              'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                #              'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                #              'NinoSstTsRmse_2']
                # to_remove = ['BiasSstLatRmse', 'BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse',
                #              'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                #              'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                #              'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                #              'NinoSstTsRmse_2', "SeasonalSstLatRmse", "SeasonalTauxLatRmse", "SeasonalTauxLonRmse"]
                to_remove = ['BiasSstLatRmse', 'BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse',
                             'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                             'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                             'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                             'NinoSstTsRmse_2', "SeasonalSstLatRmse", "SeasonalTauxLatRmse", "SeasonalTauxLonRmse"]
            elif mc == "ENSO_proc":
                to_remove = ['EnsoAmpl', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstSwr',
                             'EnsoFbSstShf']
            else:
                to_remove = ['EnsoAmpl', 'EnsoSlpMap', 'EnsoSstLonRmse', 'NinaPrMap_1', 'NinaPrMap_2', 'NinaSlpMap_1',
                             'NinaSlpMap_2',
                             'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstMap_1', 'NinaSstMap_2', 'NinoPrMap_1',
                             'NinoPrMap_2', 'NinoSlpMap_1', 'NinoSlpMap_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                             'NinoSstMap_1', 'NinoSstMap_2']
        else:
            if mc == "ENSO_perf":
                # to_remove = []
                to_remove = ['BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoTauxTsRmse', 'SeasonalTauxLatRmse',
                             'SeasonalTauxLonRmse']
            elif mc == "ENSO_proc":
                to_remove = ['EnsoAmpl']
            else:
                to_remove = ['EnsoAmpl', 'EnsoSstLonRmse', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2',
                             'NinoSstLonRmse_1', 'NinoSstLonRmse_2']
        for met in to_remove:
            while met in list_metrics:
                list_metrics.remove(met)
        if mc == "ENSO_tel":
            list_metrics = add_suffix(list_metrics)
        # # !!!!! temporary: start !!!!!
        # # ssh metrics are not computed yet (ask jiwoo)
        # list_metrics2 = deepcopy(list_metrics)
        # for met in list_metrics2:
        #     if "Ssh" in met:
        #         while met in list_metrics:
        #             list_metrics.remove(met)
        # del list_metrics2
        # # slp metrics are wrong (error in observation?)
        # list_metrics2 = deepcopy(list_metrics)
        # for met in list_metrics2:
        #     if "Slp" in met:
        #         while met in list_metrics:
        #             list_metrics.remove(met)
        # del list_metrics2
        # # !!!!! temporary: end !!!!!
        # read json
        # lpath = OSpath__join(path_in, proj + "/" + experiment + "/" + mc)
        lpath = deepcopy(path_in)
        # lname = proj + "_" + experiment + "_" + mc + "_v2019????.json"
        lname = proj + "_" + experiment + "_" + mc + "_v2019*.json"
        filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
        with open(filename_js) as ff:
            data = json.load(ff)
        ff.close()
        list_models = sorted(data["RESULTS"]["model"].keys(), key=lambda v: v.upper())
        # read metrics
        dict1 = dict()
        for mod in list_models:
            dict2 = dict()
            for met in list_metrics:
                try: data["RESULTS"]["model"][mod]["value"][met]["metric"]
                except: tmp = 1e20
                else:
                    key1 = data["RESULTS"]["model"][mod]["value"][met]["metric"].keys()[0]
                    tmp = 1e20 if data["RESULTS"]["model"][mod]["value"][met]["metric"][key1]["value"] is None \
                        else data["RESULTS"]["model"][mod]["value"][met]["metric"][key1]["value"]
                dict2[met] = tmp
                del tmp
                # try:
                #     defCollection(mc)["metrics_list"][met]["metric_computation"]
                # except:
                #     tmp = data["RESULTS"]["model"][mod]["value"][met]["metric"]
                #     dict2[met] = dict(
                #         (key, 1e20 if (("Taux" in met and mod == "BCC-ESM1")
                #                        or tmp[key]["value"] is None) else tmp[key]["value"]) for key in tmp.keys())
                #     del tmp
                # else:
                #     tmp = data["RESULTS"]["model"][mod]["value"][met]["diagnostic"]
                #     list1 = tmp.keys()
                #     list1.remove(mod)
                #     mod_val = tmp[mod]["value"]
                #     dict3 = dict()
                #     for key in list1:
                #         obs_val = tmp[key]["value"]
                #         if ("Taux" in met and mod == "BCC-ESM1") or obs_val is None or mod_val is None:
                #             dict3[key] = 1e20
                #         else:
                #             dict3[key] = abs(float(mod_val - obs_val) / obs_val)
                #         del obs_val
                #     dict2[met] = dict3
                #     del dict3, list1, mod_val, tmp
            dict1[mod] = dict2
            del dict2
        # save in common dictionary
        dict_mc = common_save(dict1, dict_out=dict_mc)
        del data, dict1, ff, filename_js, list_metrics, list_models, lname, lpath
    if big_ensemble is not True:
        dict_met[proj] = dict_mc
        del dict_mc
if big_ensemble is True:
    dict_met = deepcopy(dict_mc)
    del dict_mc


# show dictionary levels
lev1 = sorted(dict_met.keys(), key=lambda v: v.upper())
print "level1 (" + str(len(lev1)) + ") = " + str(lev1)
print ""
# check metrics
if big_ensemble is True:
    list1 = list(set([len(dict_met[key1].keys()) for key1 in lev1]))
    for key1 in lev1:
        if len(dict_met[key1].keys()) == max(list1):
            list_metrics = sorted(dict_met[key1].keys(), key=lambda v: v.upper())
            pass
    # !!!!! temporary: start !!!!!
    # some models are not used in all metric collection (ask jiwoo)
    for key2 in list_metrics:
        for key1 in lev1:
            if key2 not in dict_met[key1].keys():
                dict_met[key1][key2] = dict((ref, 1e20) for ref in dict_met["ACCESS1-0"][key2].keys())
    list1 = list(set([len(dict_met[key1].keys()) for key1 in lev1]))
    # !!!!! temporary: end !!!!!
else:
    list1 = list(set([len(dict_met[key1][key2].keys()) for key1 in lev1 for key2 in dict_met[key1].keys()]))
    for key1 in lev1:
        for key2 in dict_met[key1].keys():
            if len(dict_met[key1][key2].keys()) == max(list1):
                list_metrics = sorted(dict_met[key1][key2].keys(), key=lambda v: v.upper())
                pass
    # !!!!! temporary: start !!!!!
    # some models are not used in all metric collection (ask jiwoo)
    for key3 in list_metrics:
        for key1 in lev1:
            for key2 in dict_met[key1].keys():
                if key3 not in dict_met[key1][key2].keys():
                    if isinstance(dict_met["cmip5"]["ACCESS1-0"][key3], dict):
                        dict_met[key1][key2][key3] = dict((ref, 1e20)
                                                          for ref in dict_met["cmip5"]["ACCESS1-0"][key3].keys())
                    else:
                        dict_met[key1][key2][key3] = 1e20
    list1 = list(set([len(dict_met[key1][key2].keys()) for key1 in lev1 for key2 in dict_met[key1].keys()]))
    # !!!!! temporary: end !!!!!
if len(list1) != 1:
    for key1 in lev1:
        lev2 = sorted(dict_met[key1].keys(), key=lambda v: v.upper())
        if big_ensemble is True:
            if len(lev2) != len(list_metrics):
                print key1.rjust(15) + " (" + str(len(lev2)).zfill(2) + "), missing = " +\
                      str(list(set(list_metrics) - set(lev2)))
        else:
            for key2 in lev2:
                lev3 = sorted(dict_met[key1][key2].keys(), key=lambda v: v.upper())
                if len(lev3) != len(list_metrics):
                    print key2.rjust(15) + " (" + str(len(lev3)).zfill(2) + "), missing = " +\
                          str(list(set(list_metrics) - set(lev3)))
                del lev3
        del lev2
    list_strings = [
        "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
        ": not the same number of metrics for all models",
        str().ljust(5) + "metric number = " + str(list1),
    ]
    EnsoErrorsWarnings.MyError(list_strings)
else:
    print str().ljust(5) + "metrics (" + str(len(list_metrics)) + ") = " + str(list_metrics)
del list1
for ii in range(3): print ""


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
#                 EnsoErrorsWarnings.MyWarning(list_strings)
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
#                     EnsoErrorsWarnings.MyWarning(list_strings)
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
dict_out = deepcopy(dict_met)


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
list_metrics = [
    "BiasPrLatRmse", "BiasPrLonRmse", "BiasSstLonRmse", "SeasonalPrLatRmse", "SeasonalPrLonRmse", "SeasonalSstLonRmse",
    "EnsoAmpl", "EnsoSeasonality", "EnsoSstSkew", "EnsoDuration", "NinoSstDiversity_2", "EnsodSstOce_2", "EnsoFbSshSst",
    "EnsoFbSstTaux", "EnsoFbSstThf", "EnsoFbTauxSsh", "EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoSstMapCorr",
    "EnsoSstMapRmse", "EnsoSstLonRmse", "EnsoSstTsRmse",
]
# mean metric evaluation
if ' ':
    if big_ensemble is True:
        tab_bst, tab_val = list(), list()
        for met in list_metrics:
            tab1, tab2 = list(), list()
            for grp in my_project:
                if grp in dict_selection.keys():
                    tab = list()
                    for mod in dict_selection[grp]:
                        if dict_out[mod][met] != 1e20:
                            if "Corr" in met:
                                tab.append(1 - dict_out[mod][met])
                            else:
                                tab.append(dict_out[mod][met])
                    nbr = len(tab)
                    tab = NUMPYmean(tab)
                    tab1.append(tab)
                    tab2.append([1e20, 1e20])
                    del tab
                else:
                    tab = NUMPYarray([dict_out[mod][met] for mod in lev1 if dict_out[mod][met] != 1e20])
                    bst = bootstrap(tab, nech=nbr)
                    tab = NUMPYmean(tab)
                    tab1.append(tab)
                    tab2.append(bst)
                    del bst, nbr, tab
            tab_bst.append(tab2)
            tab_val.append(tab1)
        tab_bst = NUMPYmoveaxis(NUMPYarray(tab_bst), 0, 1)
        tab_bst = NUMPYma__masked_where(tab_bst == 1e20, tab_bst)
        tab_val = NUMPYmoveaxis(NUMPYarray(tab_val), 0, -1)
        # figure_name = OSpath__join(path_out, "metrics_comparison_" + str(len(list_metrics)).zfill(2) +
        #                            "metrics_" + str(len(my_project)).zfill(2) + "selections_"+my_project[0] + "_v2")
        figure_name = OSpath__join(path_out, "metrics_comparison_" + str(len(list_metrics)).zfill(2) +
                                   "metrics_" + str(len(my_project)).zfill(2) + "selections_" + my_project[0])
        title = "metrics comparison"
        colors = ["r", "dodgerblue"]
        plot_metrics(tab_val, figure_name, title=title, x_names=list_metrics, y_name="", colors=colors, tab_bst=tab_bst,
                     legend=my_project)
        del colors, figure_name, tab_bst, tab_val, title
    else:
        tab_bst, tab_val = list(), list()
        for met in list_metrics:
            tab1, tab2 = list(), list()
            for grp in list_project:
                tmp = dict_out[grp]
                tab = NUMPYarray([tmp[mod][met] for mod in tmp.keys() if tmp[mod][met] != 1e20])
                tab1.append(float(NUMPYmean(tab)))
                if grp == "cmip6":
                    nbr = len(tab)
                    tab2.append([1e20, 1e20])
                else:
                    bst = bootstrap(tab, nech=nbr)
                    tab2.append(bst)
                    del bst, nbr
                del tab, tmp
            tab_bst.append(tab2)
            tab_val.append(tab1)
        tab_bst = NUMPYmoveaxis(NUMPYarray(tab_bst), 0, 1)
        tab_bst = NUMPYma__masked_where(tab_bst == 1e20, tab_bst)
        tab_val = NUMPYmoveaxis(NUMPYarray(tab_val), 0, -1)
        # figure_name = OSpath__join(path_out, "metrics_comparison_" + str(len(list_metrics)).zfill(2) +
        #                            "metrics_" + str(len(my_project)).zfill(2) + "selections_"+my_project[0] + "_v2")
        figure_name = OSpath__join(path_out, "metrics_comparison_" + str(len(list_metrics)).zfill(2) +
                                   "metrics_cmip5_vs_cmip6")
        title = "metrics comparison"
        colors = ["r", "dodgerblue"]
        legend = [proj.upper() for proj in list_project]
        plot_metrics(tab_val, figure_name, title=title, x_names=list_metrics, y_name="", colors=colors, tab_bst=tab_bst,
                     legend=legend)
        nbr_bet, nbr_wor = 0, 0
        for ii, met in enumerate(list_metrics):
            if tab_val[0][ii] < min(tab_bst[1][ii]):
                print list_project[0] + " significantly better: " + met
                nbr_bet += 1
            elif tab_val[0][ii] > max(tab_bst[1][ii]):
                print list_project[0] + " significantly  worse: " + met
                nbr_wor += 1
        print list_project[0] + " significantly better (" + str(nbr_bet).zfill(2) + ") and worse (" +\
              str(nbr_wor).zfill(2) + ")"
        del colors, figure_name, legend, tab_bst, tab_val, title

models = [
    "AWI-CM-1-1-MR", "BCC-CSM2-MR", "BCC-ESM1", "CAMS-CSM1-0", "FGOALS-f3-L", "FGOALS-g3", "CanESM5", "IITM-ESM",
    "CNRM-CM6-1", "CNRM-ESM2-1", "E3SM-1-0", "EC-Earth3", "EC-Earth3-Veg", "FIO-ESM-2-0", "IPSL-CM6A-LR", "MIROC6",
    "MIROC-ES2L", "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "UKESM1-0-LL", "MRI-ESM2-0", "GISS-E2-1-G", "GISS-E2-1-H",
    "CESM2", "CESM2-WACCM", "NorCPM1", "NorESM2-LM", "GFDL-AM4", "GFDL-CM4", "GFDL-ESM4", "NESM3", "SAM0-UNICON",
    "MCM-UA-1-0"
]

