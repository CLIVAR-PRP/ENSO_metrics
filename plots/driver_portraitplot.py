# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      plot portraitplots based on modified the metric values
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
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgba
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from numpy import arange as NUMPYarange
from numpy import linspace as NUMPYlinspace
from numpy import mean as NUMPYmean
from numpy import moveaxis as NUMPYmoveaxis
from numpy import std as NUMPYstd
from numpy.ma import array as NUMPYma__array
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros
from os.path import join as OSpath__join
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile
# ENSO_metrics functions
import EnsoErrorsWarnings
from EnsoPlotLib import plot_param


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collections = ["ENSO_perf", "ENSO_tel", "ENSO_proc"]
# metric_collections = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
experiment = "historical"  # "piControl" #
list_project = ["cmip5", "cmip6"]
selection = "select8"

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_12_report"
path_in = OSpath__join(path_main, "Data_grouped")
# path_out = OSpath__join(path_main, "Plots_v5")
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/v02"
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_12_09_AGU/Poster"

expe = "hist" if experiment == "historical" else "pi"

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
dict_mc = {"ENSO_perf": "ENSO performance", "ENSO_proc": "ENSO processes", "ENSO_tel": "ENSO teleconnections"}


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
def remove_metrics(list_met, metric_collection):
    list_met1 = deepcopy(list_met)
    if mc == "ENSO_perf":
        # to_remove = ['BiasTauxLatRmse', 'BiasTauxLonRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse',
        #              'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2',
        #              'NinaSstTsRmse_1', 'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1', 'NinoSstDur_2',
        #              'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1', 'NinoSstTsRmse_2']
        to_remove = ['BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse', 'NinaSstDur_1',
                     'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                     'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                     'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                     'NinoSstTsRmse_2', "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
    elif mc == "ENSO_proc":
        to_remove = ['EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf', 'EnsoFbSstSwr']
    else:
        # to_remove = ['EnsoPrMapStd', 'EnsoSlpMapCorr', 'EnsoSlpMapRmse', 'EnsoSlpMapStd',
        #              'EnsoSstMapStd', 'EnsoSstLonRmse',
        #              'NinaPrMap_1Corr', 'NinaPrMap_1Rmse', 'NinaPrMap_1Std',
        #              'NinaPrMap_2Corr', 'NinaPrMap_2Rmse', 'NinaPrMap_2Std',
        #              'NinaSlpMap_1Corr', 'NinaSlpMap_1Rmse', 'NinaSlpMap_1Std',
        #              'NinaSlpMap_2Corr', 'NinaSlpMap_2Rmse', 'NinaSlpMap_2Std',
        #              'NinaSstMap_1Corr', 'NinaSstMap_1Rmse', 'NinaSstMap_1Std',
        #              'NinaSstMap_2Corr', 'NinaSstMap_2Rmse', 'NinaSstMap_2Std',
        #              'NinaSstLonRmse_1', 'NinaSstLonRmse_2',
        #              'NinoPrMap_1Corr', 'NinoPrMap_1Rmse', 'NinoPrMap_1Std',
        #              'NinoPrMap_2Corr', 'NinoPrMap_2Rmse', 'NinoPrMap_2Std',
        #              'NinoSlpMap_1Corr', 'NinoSlpMap_1Rmse', 'NinoSlpMap_1Std',
        #              'NinoSlpMap_2Corr', 'NinoSlpMap_2Rmse', 'NinoSlpMap_2Std',
        #              'NinoSstMap_1Corr', 'NinoSstMap_1Rmse', 'NinoSstMap_1Std',
        #              'NinoSstMap_2Corr', 'NinoSstMap_2Rmse', 'NinoSstMap_2Std',
        #              'NinoSstLonRmse_1', 'NinoSstLonRmse_2']
        to_remove = ['EnsoPrMapStd', 'EnsoSlpMapStd', 'EnsoSstMapStd',
                     'NinaPrMap_1Corr', 'NinaPrMap_1Rmse', 'NinaPrMap_1Std',
                     'NinaPrMap_2Corr', 'NinaPrMap_2Rmse', 'NinaPrMap_2Std',
                     'NinaSlpMap_1Corr', 'NinaSlpMap_1Rmse', 'NinaSlpMap_1Std',
                     'NinaSlpMap_2Corr', 'NinaSlpMap_2Rmse', 'NinaSlpMap_2Std',
                     'NinaSstMap_1Corr', 'NinaSstMap_1Rmse', 'NinaSstMap_1Std',
                     'NinaSstMap_2Corr', 'NinaSstMap_2Rmse', 'NinaSstMap_2Std',
                     'NinaSstLonRmse_1', 'NinaSstLonRmse_2',
                     'NinoPrMap_1Corr', 'NinoPrMap_1Rmse', 'NinoPrMap_1Std',
                     'NinoPrMap_2Corr', 'NinoPrMap_2Rmse', 'NinoPrMap_2Std',
                     'NinoSlpMap_1Corr', 'NinoSlpMap_1Rmse', 'NinoSlpMap_1Std',
                     'NinoSlpMap_2Corr', 'NinoSlpMap_2Rmse', 'NinoSlpMap_2Std',
                     'NinoSstMap_1Corr', 'NinoSstMap_1Rmse', 'NinoSstMap_1Std',
                     'NinoSstMap_2Corr', 'NinoSstMap_2Rmse', 'NinoSstMap_2Std',
                     'NinoSstLonRmse_1', 'NinoSstLonRmse_2', ]
    # remove given metrics
    list_met2 = deepcopy(list_met1)
    for met in list_met2:
        if met in to_remove:
            while met in list_met1:
                list_met1.remove(met)
    del list_met2
    # !!!!! temporary: start !!!!!
    # # ssh metrics are not computed yet (ask jiwoo)
    # list_met2 = deepcopy(list_met1)
    # for met in list_met2:
    #     if "Ssh" in met:
    #         while met in list_met1:
    #             list_met1.remove(met)
    # del list_met2
    # slp metrics are wrong (error in observation?)
    list_met2 = deepcopy(list_met1)
    for met in list_met2:
        if "Slp" in met:
            while met in list_met1:
                list_met1.remove(met)
    del list_met2
    # I do not use Nina metrics
    list_met2 = deepcopy(list_met1)
    for met in list_met2:
        if "Nina" in met:
            while met in list_met1:
                list_met1.remove(met)
    del list_met2
    # I do not use std metrics (ratio of std in maps)
    list_met2 = deepcopy(list_met1)
    for met in list_met2:
        if "Std" in met:
            while met in list_met1:
                list_met1.remove(met)
    del list_met2
    # !!!!! temporary: end !!!!!
    return list_met1


def get_reference(metric_collection, metric):
    if metric_collection == "ENSO_tel" and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    return plot_param(metric_collection, my_met)['metric_reference']


def my_colorbar(mini=-1., maxi=1., nbins=20):
    levels = MaxNLocator(nbins=nbins).tick_values(mini, maxi)
    cmap = plt.get_cmap("cmo.balance")  # plt.get_cmap("coolwarm")
    newcmp1 = cmap(NUMPYlinspace(0.15, 0.85, 256))
    newcmp2 = cmap(NUMPYlinspace(0.0, 1.0, 256))
    #newcmp1[120:136, :] = to_rgba('snow')
    newcmp1 = ListedColormap(newcmp1)
    newcmp1.set_bad(color='grey')
    newcmp1.set_over(newcmp2[-30])  # color='darkred')
    newcmp1.set_under(newcmp2[29])
    norm = BoundaryNorm(levels, ncolors=newcmp1.N)
    return newcmp1, norm


def portraitplot(tab, name_plot, x_names, y_names, title='portraitplot', write_metrics=False, levels=None):
    if levels is None:
        levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    fig, ax = plt.subplots(figsize=(0.5 * len(tab[0]), 0.5 * len(tab)))
    # shading
    cmap, norm = my_colorbar(mini=min(levels), maxi=max(levels))
    cs = plt.pcolormesh(tab, cmap=cmap, norm=norm)
    # title
    yy1, yy2 = ax.get_ylim()
    dy = 1 + (0.5 / (yy2 - yy1))
    ax.set_title(title, fontsize=40, y=dy, loc="center")
    # x axis
    ticks = [ii + 0.5 for ii in range(len(x_names))]
    ax.set_xticks(ticks)
    ax.set_xticklabels(x_names)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    # y axis
    ticks = [ii + 0.5 for ii in range(len(y_names))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(y_names)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    # lines
    for ii in range(1, len(tab)):
        ax.axhline(ii, color='k', linestyle='-', linewidth=1)
    for ii in range(1, len(tab[0])):
        ax.axvline(ii, color='k', linestyle='-', linewidth=1)
    # text
    if write_metrics is True:
        for jj in range(len(tab[0])):
            for ii in range(len(tab)):
                if tab.mask[ii, jj] == False:
                    plt.text(jj + 0.5, ii + 0.5, str(round(tab[ii, jj], 1)), fontsize=10, horizontalalignment='center',
                             verticalalignment='center')
    # color bar
    x2 = ax.get_position().x1
    y1 = ax.get_position().y0
    y2 = ax.get_position().y1
    cax = plt.axes([x2 + 0.02, y1, 0.02, y2 - y1])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
    cbar.ax.tick_params(labelsize=20)
    # save fig
    plt.savefig(name_plot, bbox_inches='tight')
    plt.close()
    return


def multiportraitplot(tab, name_plot, xlabel, ylabel, x_names, y_names, title='portraitplot', write_metrics=False,
                      my_text="", levels=None):
    if levels is None:
        levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    nbrl = sum([len(tab[ii]) for ii in range(len(tab))]) + (len(tab) - 1) * 3
    fig = plt.figure(0, figsize=(0.5 * len(tab[0][0]), 0.5 * nbrl))
    gs = GridSpec(nbrl, 1)
    # colorbar
    cmap, norm = my_colorbar(mini=min(levels), maxi=max(levels))
    count = 0
    for kk, tmp in enumerate(tab):
        ax = plt.subplot(gs[count: count + len(tmp), 0])
        # shading
        cs = ax.pcolormesh(tmp, cmap=cmap, norm=norm)
        # title
        yy1, yy2 = ax.get_ylim()
        dy = 0.5 / (yy2 - yy1)
        # try: ax.set_title(title[kk], fontsize=40, y=1+dy, loc="center")
        # except: ax.set_title(title, fontsize=40, y=1+dy, loc="center")
        try: ax.set_title(title[kk], fontsize=40, y=1+dy, loc="left")
        except: ax.set_title(title, fontsize=40, y=1+dy, loc="left")
        # x axis
        ticks = [ii + 0.5 for ii in range(len(x_names))]
        ax.set_xticks(ticks)
        if kk != len(tab) - 1:
            ax.set_xticklabels([""] * len(ticks))
        else:
            ax.set_xticklabels(x_names)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(20)
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
        try: ax.set_xlabel(xlabel[kk], y=dy, fontsize=40)
        except: ax.set_xlabel(xlabel, y=dy, fontsize=40)
        # y axis
        ticks = [ii + 0.5 for ii in range(len(y_names[kk]))]
        ax.set_yticks(ticks)
        ax.set_yticklabels(y_names[kk])
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(20)
        xx1, xx2 = ax.get_xlim()
        dx = 0.5 / (xx2 - xx1)
        try: ax.set_ylabel(ylabel[kk], fontsize=40, verticalalignment='top')
        except: ax.set_ylabel(ylabel, fontsize=40, verticalalignment='top')
        ax.yaxis.set_label_coords(-20 * dx, 0.5)
        # lines
        for ii in range(1, len(tmp)):
            ax.axhline(ii, color='k', linestyle='-', linewidth=1)
        for ii in range(1, len(tmp[0])):
            ax.axvline(ii, color='k', linestyle='-', linewidth=1)
        # text
        if write_metrics is True:
            for jj in range(len(tmp[0])):
                for ii in range(len(tmp)):
                    if tmp.mask[ii, jj] == False:
                        plt.text(jj + 0.5, ii + 0.5, str(round(tmp[ii, jj], 1)), fontsize=10,
                                 horizontalalignment='center', verticalalignment='center')
        if kk == 0:
            x2 = ax.get_position().x1
            y2 = ax.get_position().y1
        if kk == len(tab) - 1:
            y1 = ax.get_position().y0
        count += len(tmp) + 3
    # plt.text(0., -14 * dy, my_text, fontsize=30, horizontalalignment='right', verticalalignment='center',
    #          transform=ax.transAxes)
    plt.text(0., -14 * dy, my_text, fontsize=25, horizontalalignment='right', verticalalignment='center',
             transform=ax.transAxes)
    # color bar
    cax = plt.axes([x2 + 0.01, y1, 0.015, y2 - y1])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
    # cbar.ax.set_yticklabels(levels, fontsize=30, weight='bold')
    cbar.ax.set_yticks(levels)
    fontdict = {"fontsize": 30, "fontweight": "bold"}
    cbar.ax.set_yticklabels(["-2 $\sigma$", "-1", "MMM", "1", "2 $\sigma$"], fontdict=fontdict)
    cax.annotate('better', xy=(2.5, 0.05), xycoords='axes fraction', xytext=(2.5, 0.35), fontsize=30, weight="bold",
                 rotation="vertical", horizontalalignment="center", verticalalignment='bottom',
                 arrowprops=dict(facecolor="k", width=8, headwidth=30, headlength=30, shrink=0.05))
    cax.annotate('worse', xy=(2.5, 0.95), xycoords='axes fraction', xytext=(2.5, 0.65), fontsize=30, weight="bold",
                 rotation="vertical", horizontalalignment="center", verticalalignment='top',
                 arrowprops=dict(facecolor="k", width=8, headwidth=30, headlength=30, shrink=0.05))
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
    return filename_js, data["RESULTS"]["model"]


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
tab_all, y_names_all = list(), list()
for mc in metric_collections:
    dict1 = dict()
    dict_mod = dict()
    for proj in list_project:
        # open and read json file
        filename_js, data_json = read_data(proj, mc)
        # read metrics
        list_models = sorted(data_json.keys(), key=lambda v: v.upper())
        dict_mod[proj] = list_models
        for mod in list_models:
            data_mod = data_json[mod]["value"]
            list_metrics = sorted(data_mod.keys(), key=lambda v: v.upper())
            list_metrics = remove_metrics(list_metrics, mc)
            dict2 = dict()
            for met in list_metrics:
                try: ref = get_reference(mc, met)
                except: ref = "Tropflux"
                data_met = data_mod[met]["metric"]
                list_ref = sorted(data_met.keys(), key=lambda v: v.upper())
                if ref in list_ref:
                    my_ref = deepcopy(ref)
                else:
                    if "ERA-Interim" in list_ref:
                        my_ref = "ERA-Interim"
                    elif "ERA-Interim_ERA-Interim" in list_ref:
                        my_ref = "ERA-Interim_ERA-Interim"
                    # else:
                    #     list_strings = [
                    #         "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    #         ": cannot fing a proper reference",
                    #         str().ljust(5) + "project '" + proj + "', MC '" + mc + "', model '" + mod + "', metric '" +
                    #         met + "'",
                    #         str().ljust(10) + "input references = " + str(list_ref)]
                    #     EnsoErrorsWarnings.my_error(list_strings)
                    # list_strings = [
                    #     "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    #     ": reference changed",
                    #     str().ljust(5) + "project '" + proj + "', MC '" + mc + "', model '" + mod + "', metric '" +
                    #     met + "'",
                    #     str().ljust(10) + "reference set to '" + my_ref +"' instead of '" + ref + "'"]
                    # EnsoErrorsWarnings.my_warning(list_strings)
                if data_met[my_ref]["value"] is None:
                    dict2[met] = 1e20
                else:
                    dict2[met] = data_met[my_ref]["value"]
                del data_met, my_ref, list_ref, ref
            dict1[mod] = dict2
            del data_mod, dict2, list_metrics
        del data_json, list_models
    # shape data
    # my_models = ["ACCESS1-0", "ACCESS1-3", "BCC-CSM1-1", "BCC-CSM1-1-M", "BCC-ESM1", "BNU-ESM", "CCSM4", "CESM1-BGC",
    #              "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2", "CESM2-WACCM", "CMCC-CESM", "CMCC-CM",
    #              "CMCC-CMS", "CNRM-CM5", "CNRM-CM5-2", "CSIRO-Mk3-6-0", "CanCM4", "CanESM2", "EC-Earth3-Veg",
    #              "FGOALS-g2", "GFDL-CM2p1", "GFDL-CM3", "GFDL-CM4", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H",
    #              "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H", "HadCM3",
    #              "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "INMCM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
    #              "IPSL-CM6A-LR", "MIROC4h", "MIROC5", "MIROC6", "MIROC-ESM", "MIROC-ESM-CHEM", "MPI-ESM-LR",
    #              "MPI-ESM-MR", "MPI-ESM-P", "MRI-ESM2-0", "NorESM1-M", "NorESM1-ME", "SAM0-UNICON"]
    # my_models = ["ACCESS1-0", "ACCESS1-3", "BCC-CSM1-1", "BCC-CSM1-1-M", "BCC-ESM1", "CCSM4", "CESM1-BGC",
    #              "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2", "CESM2-WACCM", "CMCC-CESM", "CMCC-CM",
    #              "CMCC-CMS", "CNRM-CM5", "CNRM-CM5-2", "CSIRO-Mk3-6-0", "CanCM4", "CanESM2", "EC-Earth3-Veg",
    #              "FGOALS-g2", "GFDL-CM2p1", "GFDL-CM3", "GFDL-CM4", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H",
    #              "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
    #              "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "INMCM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
    #              "IPSL-CM6A-LR", "MIROC4h", "MIROC5", "MIROC6", "MIROC-ESM", "MIROC-ESM-CHEM", "MPI-ESM-LR",
    #              "MPI-ESM-MR", "MPI-ESM-P", "MRI-ESM2-0", "NorESM1-M", "NorESM1-ME", "SAM0-UNICON"]
    # my_models = ["ACCESS1-0", "ACCESS1-3", "BCC-CSM1-1", "BCC-CSM1-1-M", "BCC-ESM1", "CCSM4", "CESM1-BGC",
    #              "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2", "CESM2-WACCM", "CMCC-CESM", "CMCC-CM",
    #              "CMCC-CMS", "CNRM-CM5", "CNRM-CM5-2", "CNRM-CM6-1", "CNRM-ESM2-1", "CSIRO-Mk3-6-0", "CanCM4",
    #              "CanESM2", "CanESM5", "EC-Earth3-Veg", "FGOALS-g2", "FGOALS-s2", "GFDL-CM2p1", "GFDL-CM3", "GFDL-CM4",
    #              "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-H", "GISS-E2-H-CC",
    #              "GISS-E2-1-H", "GISS-E2-R", "GISS-E2-R-CC", "HadGEM2-AO", "HadGEM2-CC",
    #              "HadGEM2-ES", "HadGEM3-GC31-LL", "INMCM4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",
    #              "IPSL-CM6A-LR", "MIROC4h", "MIROC5", "MIROC6", "MIROC-ES2L", "MIROC-ESM", "MIROC-ESM-CHEM",
    #              "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P", "MRI-CGCM3", "MRI-ESM1", "MRI-ESM2-0", "NorESM1-M",
    #              "NorESM1-ME", "NorCPM1", "NorESM2-LM", "SAM0-UNICON", "UKESM1-0-LL"]
    my_models = ['ACCESS1-0', 'ACCESS1-3', 'BCC-CSM1-1', 'BCC-CSM1-1-M', 'BCC-CSM2-MR', 'BCC-ESM1', 'CAMS-CSM1-0',
                 'CanCM4', 'CanESM2',
                 'CanESM5', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CESM1-FASTCHEM', 'CESM1-WACCM', 'CESM2', 'CESM2-WACCM',
                 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CNRM-CM5-2', 'CNRM-CM6-1', 'CNRM-CM6-1-HR',
                 'CNRM-ESM2-1', 'CSIRO-Mk3-6-0', 'EC-Earth3-Veg', 'FGOALS-g2', 'FGOALS-s2', 'GFDL-CM3', 'GFDL-CM4',
                 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-1-G', 'GISS-E2-1-G-CC', 'GISS-E2-H', 'GISS-E2-H-CC',
                 'GISS-E2-1-H', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES',
                 'HadGEM3-GC31-LL', 'INMCM4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'IPSL-CM6A-LR',
                 'MIROC4h', 'MIROC5', 'MIROC6', 'MIROC-ES2L', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR',
                 'MPI-ESM-MR', 'MPI-ESM-P', 'MRI-ESM2-0', 'NESM3', 'NorESM1-M', 'NorESM1-ME', 'NorCPM1', 'NorESM2-LM',
                 'SAM0-UNICON', 'UKESM1-0-LL']
    for mod in my_models:
        try: dict1[mod]
        except: pass
        else:
            try: my_metrics
            except: my_metrics = dict1[mod].keys()
            else: my_metrics += dict1[mod].keys()
    my_metrics = sorted(list(set(my_metrics)), key=lambda v: v.upper())
    if mc == "ENSO_perf":
        my_metrics = [
            "BiasPrLatRmse", "BiasPrLonRmse", "BiasSstLonRmse", "BiasTauxLonRmse", "SeasonalPrLatRmse",
            "SeasonalPrLonRmse", "SeasonalSstLonRmse", "SeasonalTauxLonRmse", "EnsoSstLonRmse", "EnsoSstTsRmse",
            "EnsoAmpl", "EnsoSeasonality", "EnsoSstSkew", "EnsoDuration", "NinoSstDiversity_2"]
    elif mc == "ENSO_proc":
        my_metrics = [
            "BiasSstLonRmse", "BiasTauxLonRmse", "EnsoSstLonRmse", "EnsoAmpl", "EnsoSeasonality", "EnsoSstSkew",
            "EnsodSstOce_2", "EnsoFbSstThf", "EnsoFbSstTaux", "EnsoFbTauxSsh", "EnsoFbSshSst"]
    else:
        my_metrics = [
            "EnsoAmpl", "EnsoSstLonRmse", "EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoSstMapCorr", "EnsoSstMapRmse"]
    my_metrics = list(reversed(my_metrics))
    # tab = NUMPYma__zeros((len(my_models) + 5, len(my_metrics)))
    tab = NUMPYma__zeros((len(my_models) + 3, len(my_metrics)))
    for ii, mod in enumerate(my_models):
        for jj, met in enumerate(my_metrics):
            try: dict1[mod][met]
            except: tab[ii, jj] = 1e20
            else: tab[ii, jj] = dict1[mod][met]
    tab = NUMPYmasked_where(tab == 1e20, tab)
    for jj, met in enumerate(my_metrics):
        tmp = tab[:-3, jj].compressed()
        med = float(SCIPYstats__scoreatpercentile(list(tmp), 50))
        mea = float(NUMPYmean(tmp))
        std = float(NUMPYstd(tmp))
        tmp = [tab[ii, jj] for ii, mod in enumerate(my_models) if mod in dict_mod["cmip5"]]
        tmp = NUMPYma__array(tmp).compressed()
        sel1 = float(NUMPYmean(tmp))
        tmp = [tab[ii, jj] for ii, mod in enumerate(my_models) if mod in dict_mod["cmip6"]]
        tmp = NUMPYma__array(tmp).compressed()
        sel2 = float(NUMPYmean(tmp))
        tmp = [tab[ii, jj] for ii, mod in enumerate(my_models) if mod in dict_selection[selection]]
        tmp = NUMPYma__array(tmp).compressed()
        sel3 = float(NUMPYmean(tmp))
        # tab[len(my_models), jj] = med
        # tab[len(my_models) + 1, jj] = mea
        # tab[len(my_models) + 2, jj] = sel1
        # tab[len(my_models) + 3, jj] = sel2
        # tab[len(my_models) + 4, jj] = sel3
        # tab[len(my_models), jj] = mea
        # tab[len(my_models) + 1, jj] = sel1
        # tab[len(my_models) + 2, jj] = sel2
        tab[len(my_models), jj] = sel1
        tab[len(my_models) + 1, jj] = sel2
        tab[len(my_models) + 2, jj] = 0
        # tab[:, jj] = (tab[:, jj] - med) / med
        tab[:, jj] = (tab[:, jj] - mea) / std
        print(met.rjust(20) + " = " + str(round(med, 1)))
        del mea, med, sel1, sel2, tmp
    tab = NUMPYmoveaxis(tab, 0, 1)
    # plot
    # levels = [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)]
    levels = list(range(-2, 3))  # [round(ii, 1) for ii in NUMPYarange(-1.5, 1.6, 0.5)]
    figure_name = OSpath__join(path_out, "portraitplot_" + str(mc) + "_modified_v2")
    title = str(dict_mc[mc]) + ": normalized by median (each row)"
    # x_names = ["* " + mod if mod in dict_mod["cmip6"] else mod for mod in my_models] +\
    #           ["(median)", "(CMIP)", "(CMIP5)", "(CMIP6)", "(selection)"]
    x_names = ["* " + mod if mod in dict_mod["cmip6"] else mod for mod in my_models] + \
              ["(CMIP5)", "(CMIP6)", "(reference)"]
              # ["(CMIP MMM)", "(CMIP5)", "(CMIP6)"]
    #portraitplot(tab, figure_name, x_names, my_metrics, title=title, levels=levels)
    tab_all.append(tab)
    y_names_all.append([met.replace("_2", "") for met in my_metrics])
    if mc == metric_collections[0]:
        x_names_all = deepcopy(x_names)
    del dict1, figure_name, my_metrics, my_models, tab, title, x_names


# all portraitplots in one plot
if ' ':
    # plot
    figure_name = OSpath__join(path_out, "portraitplot_" + str(len(metric_collections)) +
                               "metric_collections_modified_v4")
    tmp_num = ["a) ", "b) ", "c) "]
    title = [tmp_num[ii] + dict_mc[mc] + " metric collection" for ii, mc in enumerate(metric_collections)]
    # title = ""
    xlabel = ""
    ylabel = ""
    # ylabel = [dict_mc[mc].replace(" ", "\n") if mc == "ENSO_tel" else dict_mc[mc] for mc in metric_collections]
    text = "* = CMIP6 model"
    multiportraitplot(tab_all, figure_name, xlabel, ylabel, x_names_all, y_names_all, title=title, my_text=text,
                      levels=levels)
