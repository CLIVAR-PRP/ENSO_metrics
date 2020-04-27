# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots to compare groups of models
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from __future__ import print_function
import cmocean
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from numpy import array as NUMPYarray
from numpy import mean as NUMPYmean
from numpy import moveaxis as NUMPYmoveaxis
from numpy import sort as NUMPYsort
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
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
my_project = ["12 models", "CMIP"]
big_ensemble = False  # True  #
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
    # "select8": ["CCSM4", "CESM1-FASTCHEM", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
    #             "GISS-E2-R", "GISS-E2-R-CC", "MIROC-ES2L", "MIROC6"],
    # "select8": ["CCSM4", "CESM1-FASTCHEM", "CNRM-CM5", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
    #             "GISS-E2-R", "GISS-E2-R-CC", "MIROC6"],
    "12 models": ['CCSM4', 'CESM1-BGC', 'CESM1-FASTCHEM', 'CNRM-CM5', 'GFDL-ESM2M', 'GISS-E2-1-G', 'GISS-E2-1-G-CC',
                  'GISS-E2-1-H', 'GISS-E2-R', 'GISS-E2-R-CC', 'MIROC-ES2L', 'MIROC6'],
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

# path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_12_report"
# path_in = OSpath__join(path_main, "Data_grouped")
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/Review/r01"
path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2020_05_report"
path_in = OSpath__join(path_main, "Data")
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/Review/r01"

expe = "hist" if experiment == "historical" else "pi"

met_o1 = ["BiasPrLatRmse", "BiasPrLonRmse", "BiasSshLatRmse", "BiasSshLonRmse", "BiasSstLatRmse", "BiasSstLonRmse",
          "BiasTauxLatRmse", "BiasTauxLonRmse", "SeasonalPrLatRmse", "SeasonalPrLonRmse", "SeasonalSshLatRmse",
          "SeasonalSshLonRmse", "SeasonalSstLatRmse", "SeasonalSstLonRmse", "SeasonalTauxLatRmse",
          "SeasonalTauxLonRmse"]
met_o2 = ["EnsoSstLonRmse", "EnsoPrTsRmse", "EnsoSstTsRmse", "EnsoTauxTsRmse", "EnsoAmpl", "EnsoSeasonality",
          "EnsoSstSkew", "EnsoDuration", "EnsoSstDiversity_1", "EnsoSstDiversity_2", "NinoSstDiversity_1",
          "NinoSstDiversity_2"]
met_o3 = ["EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoPrMapStd", "EnsoPrMapDjfCorr", "EnsoPrMapDjfRmse", "EnsoPrMapDjfStd",
          "EnsoPrMapJjaCorr", "EnsoPrMapJjaRmse", "EnsoPrMapJjaStd", "EnsoSlpMapCorr", "EnsoSlpMapRmse",
          "EnsoSlpMapStd", "EnsoSlpMapDjfCorr", "EnsoSlpMapDjfRmse", "EnsoSlpMapDjfStd", "EnsoSlpMapJjaCorr",
          "EnsoSlpMapJjaRmse", "EnsoSlpMapJjaStd", "EnsoSstMapCorr", "EnsoSstMapRmse", "EnsoSstMapStd",
          "EnsoSstMapDjfCorr", "EnsoSstMapDjfRmse", "EnsoSstMapDjfStd", "EnsoSstMapJjaCorr", "EnsoSstMapJjaRmse",
          "EnsoSstMapJjaStd"]
met_o4 = ["EnsodSstOce_1", "EnsodSstOce_2", "EnsoFbSstThf", "EnsoFbSstSwr", "EnsoFbSstLhf", "EnsoFbSstLwr",
          "EnsoFbSstShf", "EnsoFbSstTaux", "EnsoFbTauxSsh", "EnsoFbSshSst"]
met_order = met_o1 + met_o2 + met_o3 + met_o4


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


def bootstrap(tab, num_samples=1000000, alpha=0.05, nech=None, statistic=NUMPYmean):
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
                        "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                        ": metric already in output",
                        str().ljust(5) + "metric '" + str(met) + "' for model '" + str(mod) +
                        "' is already in output dictionary",
                        str().ljust(10) + "in  input = " + str(dict_in[mod][met]),
                        str().ljust(10) + "in output = " + str(dict_out[mod][met]),
                    ]
                    EnsoErrorsWarnings.my_error(list_strings)
    return dict_out


def get_reference(metric_collection, metric):
    if metric_collection in ["ENSO_tel", "test_tel"] and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    return plot_param(metric_collection, my_met)['metric_reference']


def plot_metrics(tab_val, name_plot, title="", x_names=None, y_name="", colors=None, tab_bst=None, legend=None,
                 xticklabel="", cname=False, chigh=False, cfram=False):
    fig, ax = plt.subplots(figsize=(0.5 * len(tab_val[0]), 4))
    mylab = ["EnsoFbSstTaux", "EnsoFbSstThf"]
    # title
    plt.title(title, fontsize=20, y=1.01, loc='left')
    # x axis
    axis = list(range(len(tab_val[0])))
    ax.set_xticks(axis)
    if isinstance(x_names, list):
        # ax.set_xticklabels(x_names)
        ax.set_xticklabels([""] * len(x_names))
        for ll, txt in enumerate(x_names):
            if txt in met_o1 or txt + "_1" in met_o1 or txt + "_2" in met_o1:
                cc = "yellowgreen"
            elif txt in met_o2 or txt + "_1" in met_o2 or txt + "_2" in met_o2:
                cc = "plum"
            elif txt in met_o3 or txt + "_1" in met_o3 or txt + "_2" in met_o3:
                cc = "gold"
            else:
                cc = "turquoise"
            boxdict = dict(lw=0, facecolor=cc, pad=3, alpha=1)
            if txt in mylab:
                if cname is True:
                    ax.text(ll, 0, txt, fontsize=14, ha='right', va='top', rotation=45, color=cc, weight="bold")
                elif chigh is True:
                    ax.text(ll, -0.05, txt, fontsize=14, ha='right', va='top', rotation=45, color="k", weight="bold",
                            bbox=boxdict)
                else:
                    ax.text(ll, 0, txt, fontsize=14, ha='right', va='top', rotation=45, color="k", weight="bold")
            else:
                if cname is True:
                    ax.text(ll, 0, txt, fontsize=12, ha='right', va='top', rotation=45, color=cc)
                elif chigh is True:
                    ax.text(ll, -0.05, txt, fontsize=12, ha='right', va='top', rotation=45, color="k", bbox=boxdict)
                else:
                    ax.text(ll, 0, txt, fontsize=12, ha='right', va='top', rotation=45, color="k")
    else:
        ax.set_xticklabels(axis)
    if cfram is True:
        nn = 0
        lic, lix = list(), list()
        for cc, tmp1 in zip(["yellowgreen", "plum", "gold", "turquoise"], [met_o1, met_o2, met_o3, met_o4]):
            tmp2 = [txt for ll, txt in enumerate(x_names) if txt in tmp1 or txt + "_1" in tmp1 or txt + "_2" in tmp1]
            if len(tmp2) > 0:
                lic += [cc, cc]
                if nn == 0:
                    lix += [[-0.4, len(tmp2) + 0.5], [-0.4, len(tmp2) - 0.5]]
                    nn += len(tmp2) - 0.5
                elif nn + len(tmp2) > len(x_names) - 1:
                    lix += [[nn, nn + len(tmp2) - 0.1], [nn, nn + len(tmp2) - 0.1]]
                else:
                    lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                    nn += len(tmp2)
        lis = ["-"] * len(lic)
        liw = [5] * len(lic)
        liy = [[2, 2], [0, 0]] * int(round(float(len(lic)) / 2))
        for lc, ls, lw, lx, ly in zip(lic, lis, liw, lix, liy):
            line = Line2D(lx, ly, c=lc, lw=lw, ls=ls, zorder=10)
            line.set_clip_on(False)
            ax.add_line(line)
    ax.set_xlim([min(axis) - 0.5, max(axis) + 0.5])
    ax.tick_params(axis="x", labelsize=12, labelrotation=90)
    # y axis
    ax.set_yticks([0.5, 1.5], minor=True)
    ax.set_yticks([0, 1, 2], minor=False)
    ax.set_yticklabels(["reference", xticklabel, "2 * " + xticklabel],
                       fontdict={"fontsize": 12, "fontweight": "normal"})
    ax.set_ylim([0, 2])
    ax.set_ylabel(y_name, fontsize=15)
    # plot marker
    for ii in range(len(tab_val)):
        if isinstance(colors, list):
            col = colors[len(colors) - 1 - ii]
        else:
            col = "k"
        ind = len(tab_val) - 1 - ii
        if tab_bst is not None:
            for jj in range(len(tab_bst[ind])):
                tmp1, tmp2 = tab_val[ind][jj], tab_bst[ind][jj]
                if ind == 0:
                    tmp3, tmp4 = tab_val[1][jj], tab_bst[1][jj]
                else:
                    tmp3, tmp4 = tab_val[0][jj], tab_bst[0][jj]
                if jj in [6, 7]:
                    print(tmp1, tmp2, tmp3, tmp4)
                if (min(tmp4) <= tmp1 <= max(tmp4)) or (min(tmp2) <= tmp3 <= max(tmp2)):
                    ax.plot([jj], [tmp1], markersize=13, color="none", marker="D", fillstyle="none",
                            markeredgecolor=col, markeredgewidth=3, zorder=2)
                else:
                    ax.scatter(jj, tmp1, s=200, c=col, marker="D", zorder=2)
                if tmp2[0] > 0 and tmp2[1] > 0:
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tmp2[0], tmp2[0]], c=col, lw=2, zorder=3))
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tmp2[1], tmp2[1]], c=col, lw=2, zorder=3))
                    if ii == 0:
                        tmpl = [jj - 0.05, jj - 0.05]
                    else:
                        tmpl = [jj + 0.05, jj + 0.05]
                    ax.add_line(Line2D(tmpl, [tmp2[0], tmp2[1]], c=col, lw=2, zorder=3))
                    del tmpl
                del tmp1, tmp2, tmp3, tmp4
        else:
            ax.scatter(axis, list(tab_val[ind]), s=200, c=col, marker="D", zorder=2)
        del col
    # grid
    ax.grid(linestyle="--", linewidth=1, axis="y", which="both", zorder=1)
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
            font = {'color': col, 'weight': 'bold', 'size': 15}
            ax.text(x2 - 2 * dx, y2 - (ii + 1) * 8 * dy, legend[len(legend) - 1 - ii], horizontalalignment="right",
                    verticalalignment="center", fontdict=font)
            del col, font
        xxx, ddx, yyy, ddy = x1 + (2 * dx), deepcopy(dx), 1.75, 0.2
        ax.add_line(Line2D([xxx - ddx, xxx + ddx], [yyy + ddy, yyy + ddy], c=colors[1], lw=2))
        ax.add_line(Line2D([xxx - ddx, xxx + ddx], [yyy - ddy, yyy - ddy], c=colors[1], lw=2))
        ax.add_line(Line2D([xxx, xxx], [yyy - ddy, yyy + ddy], c=colors[1], lw=2))
        dicttext = {"horizontalalignment": "left", "verticalalignment": "center",
                    "fontdict": {'weight': 'normal', 'size': 12}, "transform": ax.transData}
        ax.text(xxx + dx, yyy, "95% confidence interval of MMM\n(Monte Carlo sampling method)", **dicttext)
        arrowdict = dict(facecolor="k", width=2, headwidth=10, headlength=10, shrink=0.0)
        ax.annotate("", xy=(-0.05, 0.05), xycoords='axes fraction', xytext=(-0.05, 0.42), fontsize=13,
                    rotation="vertical", ha="center", va='bottom', arrowprops=arrowdict)
        ax.text(-0.07, 0.42, "improved", fontsize=13, rotation="vertical", ha="center", va='top',
                transform=ax.transAxes)
        ax.annotate("", xy=(-0.05, 0.95), xycoords='axes fraction', xytext=(-0.05, 0.58), fontsize=13,
                    rotation="vertical", ha="center", va='top', arrowprops=arrowdict)
        ax.text(-0.07, 0.58, "worsened", fontsize=13, rotation="vertical", ha="center", va='bottom',
                transform=ax.transAxes)
    # save fig
    plt.savefig(name_plot, bbox_inches='tight')
    plt.close()
    return


def read_data(project, metric_collection):
    # lname = project + "_" + experiment + "_" + metric_collection + "_v2019????_modified.json"
    # filename_js = list(GLOBiglob(OSpath__join(path_in, lname)))[0]
    lpath = OSpath__join(path_in, project + "/" + experiment)
    lname = project + "_" + experiment + "_" + metric_collection + "_v20200430.json"
    filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
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
            to_remove = ['BiasSshLatRmse', 'BiasSshLonRmse', 'BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse',
                         'EnsoSstDiversity_1', 'EnsoTauxTsRmse', 'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1',
                         'NinaSstLonRmse_2', 'NinaSstTsRmse_1', 'NinaSstTsRmse_2', 'NinoSstDiversity_1',
                         'NinoSstDiversity_2', 'NinoSstDur_1', 'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                         'NinoSstTsRmse_1', 'NinoSstTsRmse_2', "SeasonalSshLatRmse", "SeasonalSshLonRmse",
                         "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
        elif metric_collection == "ENSO_proc":
            # to_remove = ['EnsoAmpl', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf',
            #              'EnsoFbTauxSsh']
            to_remove = ['BiasSshLonRmse', 'BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality',
                         'EnsoSstLonRmse', 'EnsoSstSkew', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr',
                         'EnsoFbSstShf', 'EnsoFbSstSwr']
        else:
            to_remove = ['EnsoAmpl', 'EnsoPrMapCorr', 'EnsoPrMapRmse', 'EnsoPrMapStd', 'EnsoPrMapDjfStd',
                         'EnsoPrMapJjaStd', 'EnsoSlpMapCorr', 'EnsoSlpMapRmse', 'EnsoSlpMapStd', 'EnsoSlpMapDjfCorr',
                         'EnsoSlpMapDjfRmse', 'EnsoSlpMapDjfStd', 'EnsoSlpMapJjaCorr', 'EnsoSlpMapJjaRmse',
                         'EnsoSlpMapJjaStd', 'EnsoSstMapCorr', 'EnsoSstMapRmse', 'EnsoSstMapStd', 'EnsoSstMapDjfStd',
                         'EnsoSstMapJjaStd', 'EnsoSstLonRmse',
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
            to_remove = ['BiasSshLonRmse', 'BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality',
                         'EnsoSstLonRmse', 'EnsoSstSkew']
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
        dict1 = dict()
        for mod in list_models:
            data_mod = data_json[mod][data_json[mod].keys()[0]]["value"]
            list_metrics = sorted(data_mod.keys(), key=lambda v: v.upper())
            list_metrics = remove_metrics(list_metrics, mc)
            dict2 = dict()
            for met in list_metrics:
                if mc == "ENSO_tel":
                    try:
                        ref = get_reference(mc, met)
                    except:
                        ref = get_reference(mc.replace("ENSO", "test"), met)
                else:
                    ref = get_reference(mc, met)
                data_met = data_mod[met]["metric"]
                list_ref = sorted(data_met.keys(), key=lambda v: v.upper())
                my_ref = deepcopy(ref)
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


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
dict_out = deepcopy(dict_met)
lev1 = sorted(dict_out.keys(), key=lambda v: v.upper())
# mean metric evaluation
if ' ':
    if big_ensemble is True:
        list_metrics = sorted(dict_out[lev1[0]].keys(), key=lambda v: v.upper())
        list_metrics = [met for met in met_order if met in list_metrics]
        tab_bst, tab_val = list(), list()
        for met in list_metrics:
            tab1, tab2 = list(), list()
            for grp in my_project:
                if grp in dict_selection.keys():
                    tab = list()
                    for mod in dict_selection[grp]:
                        if dict_out[mod][met] != 1e20:
                            # if "Corr" in met:
                            #     tab.append(1 - dict_out[mod][met])
                            # else:
                            #     tab.append(dict_out[mod][met])
                            tab.append(dict_out[mod][met])
                    nbr = len(tab)
                    tab = NUMPYarray(tab)
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
        tmp = NUMPYmoveaxis(NUMPYarray([tab_val[1], tab_val[1]]), 0, 1)
        tab_bst = tab_bst / tmp
        tab_val = tab_val / tab_val[1]
        figure_name = OSpath__join(path_out, "Figure_07a_cmip_vs_" + my_project[0].replace(" ", "") + "_20200430")
        title = "a) Mean metric values of a subset of models relative to CMIP"
        colors = ["orange", "forestgreen"]
        if reduced_set is True:
            list_metrics2 = [met.replace("_1", "").replace("_2", "") for met in list_metrics]
        else:
            figure_name += "_all_metrics"
            list_metrics2 = deepcopy(list_metrics)
        list_names2 = [str(len(dict_selection[my_project[0]])) + " models", my_project[1]]
        plot_metrics(tab_val, figure_name, title=title, x_names=list_metrics2, y_name="", colors=colors,
                     tab_bst=tab_bst, legend=list_names2, xticklabel="CMIP", chigh=True, cfram=True)
        del colors, figure_name, tab_bst, tab_val, title
    else:
        list_metrics = sorted(dict_out[lev1[0]][dict_out[lev1[0]].keys()[0]].keys(), key=lambda v: v.upper())
        list_metrics = [met for met in met_order if met in list_metrics]
        tab_bst, tab_val = list(), list()
        for met in list_metrics:
            tab_tmp = list()
            for grp in list_project:
                tmp = dict_out[grp]
                tab = NUMPYarray([tmp[mod][met] for mod in tmp.keys() if tmp[mod][met] != 1e20])
                tab_tmp.append(NUMPYma__masked_invalid(tab).compressed())
                del tab, tmp
            tab1, tab2 = list(), list()
            for ii in range(len(tab_tmp)):
                tab1.append(float(NUMPYmean(tab_tmp[ii])))
                if ii == 0:
                    nbr = len(tab_tmp[1])
                else:
                    nbr = len(tab_tmp[0])
                bst = bootstrap(tab_tmp[ii], nech=nbr)
                tab2.append(bst)
                del bst, nbr
            tab_bst.append(tab2)
            tab_val.append(tab1)
        tab_bst = NUMPYmoveaxis(NUMPYarray(tab_bst), 0, 1)
        tab_val = NUMPYmoveaxis(NUMPYarray(tab_val), 0, -1)
        tmp = NUMPYmoveaxis(NUMPYarray([tab_val[1], tab_val[1]]), 0, 1)
        tab_bst = tab_bst / tmp
        tab_val = tab_val / tab_val[1]
        figure_name = OSpath__join(path_out, "Figure_03_cmip5_vs_cmip6_20200430")
        title = "Mean metric values relative to CMIP5"
        colors = ["r", "dodgerblue"]
        legend = [proj.upper() for proj in list_project]
        if reduced_set is True:
            list_metrics2 = [met.replace("_1", "").replace("_2", "") for met in list_metrics]
        else:
            figure_name += "_all_metrics"
            list_metrics2 = deepcopy(list_metrics)
        plot_metrics(tab_val, figure_name, title=title, x_names=list_metrics2, y_name="", colors=colors,
                     tab_bst=tab_bst, legend=legend, xticklabel="CMIP5", chigh=True, cfram=True)
        stop
        nbr_bet, nbr_wor = 0, 0
        for ii, met in enumerate(list_metrics):
            if tab_val[0][ii] < tab_val[1][ii]:
                print(list_project[0] + " better: " + met)
                nbr_bet += 1
            elif tab_val[0][ii] > tab_val[1][ii]:
                print(list_project[0] + "  worse: " + met)
                nbr_wor += 1
        print(list_project[0] + " better (" + str(nbr_bet).zfill(2) + ") and worse (" + str(nbr_wor).zfill(2) + ")")
        print("")
        nbr_bet, nbr_wor = 0, 0
        for ii, met in enumerate(list_metrics):
            if tab_val[0][ii] < min(tab_bst[1][ii]):
                print(list_project[0] + " significantly better: " + met)
                nbr_bet += 1
            elif tab_val[0][ii] > max(tab_bst[1][ii]):
                print(list_project[0] + " significantly  worse: " + met)
                nbr_wor += 1
        print(list_project[0] + " significantly better (" + str(nbr_bet).zfill(2) + ") and worse (" +\
              str(nbr_wor).zfill(2) + ")")
        del colors, figure_name, legend, tab_bst, tab_val, title

models = [
    "AWI-CM-1-1-MR", "BCC-CSM2-MR", "BCC-ESM1", "CAMS-CSM1-0", "FGOALS-f3-L", "FGOALS-g3", "CanESM5", "IITM-ESM",
    "CNRM-CM6-1", "CNRM-ESM2-1", "E3SM-1-0", "EC-Earth3", "EC-Earth3-Veg", "FIO-ESM-2-0", "IPSL-CM6A-LR", "MIROC6",
    "MIROC-ES2L", "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "UKESM1-0-LL", "MRI-ESM2-0", "GISS-E2-1-G", "GISS-E2-1-H",
    "CESM2", "CESM2-WACCM", "NorCPM1", "NorESM2-LM", "GFDL-AM4", "GFDL-CM4", "GFDL-ESM4", "NESM3", "SAM0-UNICON",
    "MCM-UA-1-0"
]

