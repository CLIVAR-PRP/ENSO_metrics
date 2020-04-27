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
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros
from os.path import join as OSpath__join
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile
import string
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoPlotLib import plot_param


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collections = ["ENSO_perf", "ENSO_tel", "ENSO_proc"]
experiment = "historical"
list_project = ["cmip5", "cmip6"]
reduced_set = True  # False  #

# path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_12_report"
# path_in = OSpath__join(path_main, "Data_grouped")
path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2020_05_report"
path_in = OSpath__join(path_main, "Data")
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/Review/r01"

expe = "hist" if experiment == "historical" else "pi"

dict_mc = {"ENSO_perf": "Performance", "ENSO_proc": "Processes", "ENSO_tel": "Telecon."}

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

mod_order = [
    'ACCESS1-0', 'ACCESS1-3', 'ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM1-1', 'BCC-CSM1-1-M', 'BCC-CSM2-MR', 'BCC-ESM1',
    'BNU-ESM', 'CAMS-CSM1-0', 'CanCM4', 'CanESM2', 'CanESM5', 'CanESM5-CanOE', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5',
    'CESM2', 'CESM2-FV2', 'CESM1-FASTCHEM', 'CESM1-WACCM', 'CESM2-WACCM', 'CESM2-WACCM-FV2', 'CMCC-CESM', 'CMCC-CM',
    'CMCC-CMS', 'CNRM-CM5', 'CNRM-CM5-2', 'CNRM-CM6-1', 'CNRM-CM6-1-HR', 'CNRM-ESM2-1', 'CSIRO-Mk3-6-0',
    'CSIRO-Mk3L-1-2', 'E3SM-1-0', 'E3SM-1-1', 'EC-EARTH', 'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-g2', 'FGOALS-s2',
    'FIO-ESM', 'GFDL-CM2p1', 'GFDL-CM3', 'GFDL-CM4', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GFDL-ESM4', 'GISS-E2-1-G',
    'GISS-E2-1-G-CC', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-1-H', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadCM3', 'HadGEM2-AO',
    'HadGEM2-CC', 'HadGEM2-ES', 'HadGEM3-GC31-LL', 'INMCM4', 'INM-CM4-8', 'INM-CM5-0', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR',
    'IPSL-CM5B-LR', 'IPSL-CM6A-LR', 'KACE-1-0-G', 'MIROC4h', 'MIROC5', 'MIROC6', 'MIROC-ESM', 'MIROC-ESM-CHEM',
    'MIROC-ES2L', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P', 'MPI-ESM-1-2-HAM', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR',
    'MRI-CGCM3', 'MRI-ESM1', 'MRI-ESM2-0', 'NESM3', 'NorESM1-M', 'NorESM1-ME', 'NorCPM1', 'NorESM2-LM', 'NorESM2-MM',
    'SAM0-UNICON', 'TaiESM1', 'UKESM1-0-LL']


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
def remove_metrics(list_met, metric_collection):
    list_met1 = deepcopy(list_met)
    if mc == "ENSO_perf":
        to_remove = ['BiasSshLatRmse', 'BiasSshLonRmse', 'BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse',
                     'EnsoSstDiversity_1', 'EnsoTauxTsRmse', 'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1',
                     'NinaSstLonRmse_2', 'NinaSstTsRmse_1', 'NinaSstTsRmse_2', 'NinoSstDiversity_1',
                     'NinoSstDiversity_2', 'NinoSstDur_1', 'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                     'NinoSstTsRmse_1', 'NinoSstTsRmse_2', "SeasonalSshLatRmse", "SeasonalSshLonRmse",
                     "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
    elif mc == "ENSO_proc":
        to_remove = ['BiasSshLonRmse', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf', 'EnsoFbSstSwr']
    else:
        to_remove = ['EnsoPrMapCorr', 'EnsoPrMapRmse', 'EnsoPrMapStd', 'EnsoPrMapDjfStd', 'EnsoPrMapJjaStd',
                     'EnsoSlpMapCorr', 'EnsoSlpMapRmse', 'EnsoSlpMapStd', 'EnsoSlpMapDjfCorr', 'EnsoSlpMapDjfRmse',
                     'EnsoSlpMapDjfStd', 'EnsoSlpMapJjaCorr', 'EnsoSlpMapJjaRmse', 'EnsoSlpMapJjaStd', 'EnsoSstMapCorr',
                     'EnsoSstMapRmse', 'EnsoSstMapStd', 'EnsoSstMapDjfStd', 'EnsoSstMapJjaStd',
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
    if metric_collection in ["ENSO_tel", "test_tel"] and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    return plot_param(metric_collection, my_met)['metric_reference']


def my_colorbar(mini=-1., maxi=1., nbins=20, pptype=0):
    levels = MaxNLocator(nbins=nbins).tick_values(mini, maxi)
    if pptype in [0, 3, 4]:
        cmap = plt.get_cmap("cmo.balance")  # plt.get_cmap("coolwarm")
        newcmp1 = cmap(NUMPYlinspace(0.15, 0.85, 256))
        newcmp2 = cmap(NUMPYlinspace(0.0, 1.0, 256))
        newcmp1 = ListedColormap(newcmp1)
        newcmp1.set_over(newcmp2[-30])  # color='darkred')
        newcmp1.set_under(newcmp2[29])
    else:
        cmap = plt.get_cmap("cmo.amp")  # plt.get_cmap("coolwarm")
        newcmp1 = cmap(NUMPYlinspace(0.0, 1.0, 256))
        newcmp2 = cmap(NUMPYlinspace(0.0, 1.0, 256))
        newcmp1 = ListedColormap(newcmp1)
        newcmp1.set_over(newcmp2[-1])
        newcmp1.set_under(color="w")
    newcmp1.set_bad(color="k")
    norm = BoundaryNorm(levels, ncolors=newcmp1.N)
    return newcmp1, norm


def multiportraitplot(tab, name_plot, xlabel, ylabel, x_names, y_names, title='portraitplot', write_metrics=False,
                      my_text="", levels=None, cname=False, chigh=False, cfram=False, pptype=0):
    if levels is None:
        levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    fontdict = {"fontsize": 40, "fontweight": "bold"}
    nbr_space = 2
    nbrc = sum([len(tab[ii][0]) for ii in range(len(tab))]) + (len(tab) - 1) * nbr_space
    fig = plt.figure(0, figsize=(0.5 * nbrc, 0.5 * len(tab[0])))
    gs = GridSpec(1, nbrc)
    # colorbar
    cmap, norm = my_colorbar(mini=min(levels), maxi=max(levels), pptype=pptype)
    count = 0
    for kk, tmp in enumerate(tab):
        ax = plt.subplot(gs[0, count: count + len(tmp[0])])
        # shading
        cs = ax.pcolormesh(tmp, cmap=cmap, norm=norm)
        # title
        xx1, xx2 = ax.get_xlim()
        dx = 0.5 / (xx2 - xx1)
        yy1, yy2 = ax.get_ylim()
        dy = 0.5 / (yy2 - yy1)
        try: ax.set_title(title[kk], fontdict=fontdict, y=1+dy, loc="center")
        except: ax.set_title(title, fontdict=fontdict, y=1+dy, loc="center")
        # x axis
        ticks = [ii + 0.5 for ii in range(len(x_names[kk]))]
        ax.set_xticks(ticks)
        ax.set_xticklabels([] * len(ticks))
        for ll, txt in enumerate(x_names[kk]):
            if txt in met_o1 or txt + "_1" in met_o1 or txt + "_2" in met_o1:
                cc = "yellowgreen"
            elif txt in met_o2 or txt + "_1" in met_o2 or txt + "_2" in met_o2:
                cc = "plum"
            elif txt in met_o3 or txt + "_1" in met_o3 or txt + "_2" in met_o3:
                cc = "gold"
            else:
                cc = "turquoise"
            if cname is True:
                ax.text(ll + 0.5, 0, txt, fontsize=15, ha='right', va='top', rotation=45, color=cc)
            elif chigh is True:
                boxdict = dict(lw=0, facecolor=cc, pad=3, alpha=1)
                ax.text(ll + 0.5, -0.2, txt, fontsize=15, ha='right', va='top', rotation=45, color="k", bbox=boxdict)
            else:
                ax.text(ll + 0.5, 0, txt, fontsize=20, ha='right', va='top', rotation=45, color="k")
        if cfram is True:
            tmp1 = [met_o1, met_o2, met_o3, met_o4]
            nn = 0
            lix = [[0, 0]]
            for tt in tmp1:
                tmp2 = [txt for ll, txt in enumerate(x_names[kk]) if
                        txt in tt or txt + "_1" in tt or txt + "_2" in tt]
                nn += len(tmp2)
                if len(tmp2) > 0:
                    lix += [[nn, nn]]
                del tmp2
            liy = [[0, len(tab[0])]] * len(lix)
            lic, lis = ["k"] * len(lix), ["-"] * len(lix)
            for lc, ls, lx, ly in zip(lic, lis, lix, liy):
                line = Line2D(lx, ly, c=lc, lw=7, ls=ls, zorder=10)
                line.set_clip_on(False)
                ax.add_line(line)
            nn = 0
            lic, lix = list(), list()
            for uu, tt in enumerate(tmp1):
                tmp2 = [txt for ll, txt in enumerate(x_names[kk]) if
                        txt in tt or txt + "_1" in tt or txt + "_2" in tt]
                if len(tmp2) > 0:
                    if uu == 0:
                        cc = "yellowgreen"
                    elif uu == 1:
                        cc = "plum"
                    elif uu == 2:
                        cc = "gold"
                    else:
                        cc = "turquoise"
                    lic += [cc, cc]
                    if nn > 0:
                        lix += [[nn + 0.2, nn + len(tmp2)], [nn + 0.2, nn + len(tmp2)]]
                    else:
                        lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                    nn += len(tmp2)
                    del cc
                del tmp2
            liy = [[len(tab[0]), len(tab[0])], [0, 0]] * int(float(len(lix)) / 2)
            lis = ["-"] * len(lix)
            for mm, (lc, ls, lx, ly) in enumerate(zip(lic, lis, lix, liy)):
                if mm < 2:
                    line = Line2D([lx[0] + 0.05, lx[1]], ly, c=lc, lw=10, ls=ls, zorder=10)
                elif mm > len(lis) - 3:
                    line = Line2D([lx[0], lx[1] - 0.05], ly, c=lc, lw=10, ls=ls, zorder=10)
                else:
                    line = Line2D(lx, ly, c=lc, lw=10, ls=ls, zorder=10)
                line.set_clip_on(False)
                ax.add_line(line)
        try: ax.set_xlabel(xlabel[kk], y=dy, fontsize=40)
        except: ax.set_xlabel(xlabel, y=dy, fontsize=40)
        # y axis
        ticks = [ii + 0.5 for ii in range(len(y_names))]
        ax.set_yticks(ticks)
        if kk != 0:
            ax.set_yticklabels([""] * len(ticks))
        else:
            ax.text(-5 * dx, -1 * dy, my_text, fontsize=25, ha='right', va='top', transform=ax.transAxes)
            ax.tick_params(axis="y", labelsize=20)
            ax.set_yticklabels(y_names)
        try: ax.set_ylabel(ylabel[kk], fontsize=40, va='top')
        except: ax.set_ylabel(ylabel, fontsize=40, va='top')
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
                        plt.text(jj + 0.5, ii + 0.5, str(round(tmp[ii, jj], 1)), fontsize=10, ha='center', va='center')
        if kk == len(tab) - 1:
            x2 = ax.get_position().x1
            y1 = ax.get_position().y0
            y2 = ax.get_position().y1
        count += len(tmp[0]) + nbr_space
    # color bar
    cax = plt.axes([x2 + 0.03, y1, 0.02, y2 - y1])
    if pptype == 0:
        cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
        cbar.ax.set_yticklabels(["-2 $\sigma$", "-1", "MMM", "1", "2 $\sigma$"], fontdict=fontdict)
        cax.annotate("", xy=(3.7, 0.06), xycoords='axes fraction', xytext=(3.7, 0.45), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='bottom',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(5.2, 0.45, "closer to reference", fontsize=40, rotation="vertical", ha="center", va='top',
                 weight="bold")
        cax.annotate("", xy=(3.7, 0.94), xycoords='axes fraction', xytext=(3.7, 0.55), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='top',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(5.2, 0.55, "further from reference", fontsize=40, rotation="vertical", ha="center", va='bottom',
                 weight="bold")
    elif pptype == 1:
        cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='max', aspect=40)
        cbar.ax.set_yticklabels(levels, fontdict=fontdict)
        cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.00), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='top',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.5, "further from reference", fontsize=40, rotation="vertical", ha="center", va='center',
                 weight="bold")
    elif pptype == 2:
        cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='max', aspect=40)
        cbar.ax.set_yticklabels(["0", "1 $\sigma$", "2 $\sigma$", "3 $\sigma$", "4 $\sigma$"], fontdict=fontdict)
        cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.00), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='top',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.5, "further from reference", fontsize=40, rotation="vertical", ha="center", va='center',
                 weight="bold")
    elif pptype == 3:
        cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
        cbar.ax.set_yticklabels(["-2", "-1", "MMM", "1", "2"], fontdict=fontdict)
        cax.annotate("", xy=(5, 0.06), xycoords='axes fraction', xytext=(5, 0.45), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='bottom',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.45, "closer to reference", fontsize=40, rotation="vertical", ha="center", va='top',
                 weight="bold")
        cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.55), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='top',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.55, "further from reference", fontsize=40, rotation="vertical", ha="center", va='bottom',
                 weight="bold")
    elif pptype == 4:
        cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend='both', aspect=40)
        cbar.ax.set_yticklabels(["-2 $\sigma$", "-1", "MMM", "1", "2 $\sigma$"], fontdict=fontdict)
        cax.annotate("", xy=(5, 0.06), xycoords='axes fraction', xytext=(5, 0.45), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='bottom',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.45, "closer to reference", fontsize=40, rotation="vertical", ha="center", va='top',
                 weight="bold")
        cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.55), fontsize=40,
                     weight="bold", rotation="vertical", ha="center", va='top',
                     arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
        cax.text(6.5, 0.55, "further from reference", fontsize=40, rotation="vertical", ha="center", va='bottom',
                 weight="bold")
    plt.savefig(name_plot, bbox_inches='tight')
    plt.close()
    return


def my_pp(ax, index, tab, title, xnames, ynames, cname=False, chigh=False, cfram=False, cfram_bot=False,
          cfram_top=False, my_text="", plot_title=False, plot_xname=False, write_metrics=False):
    fontdict = {"fontsize": 40, "fontweight": "bold"}
    xx1, xx2 = ax.get_xlim()
    dx = 0.5 / (xx2 - xx1)
    yy1, yy2 = ax.get_ylim()
    dy = 0.5 / (yy2 - yy1)
    # title
    if plot_title is True:
        ax.set_title(title, fontdict=fontdict, y=1+dy, loc="center")
    # x axis
    ticks = [ii + 0.5 for ii in range(len(xnames))]
    ax.set_xticks(ticks)
    ax.set_xticklabels([] * len(ticks))
    if plot_xname is True:
        for ll, txt in enumerate(xnames):
            if txt in met_o1 or txt + "_1" in met_o1 or txt + "_2" in met_o1:
                cc = "yellowgreen"
            elif txt in met_o2 or txt + "_1" in met_o2 or txt + "_2" in met_o2:
                cc = "plum"
            elif txt in met_o3 or txt + "_1" in met_o3 or txt + "_2" in met_o3:
                cc = "gold"
            else:
                cc = "turquoise"
            if cname is True:
                ax.text(ll + 0.5, 0, txt, fontsize=15, ha='right', va='top', rotation=45, color=cc)
            elif chigh is True:
                boxdict = dict(lw=0, facecolor=cc, pad=3, alpha=1)
                ax.text(ll + 0.5, -0.2, txt, fontsize=15, ha='right', va='top', rotation=45, color="k", bbox=boxdict)
            else:
                ax.text(ll + 0.5, 0, txt, fontsize=20, ha='right', va='top', rotation=45, color="k")
    if cfram is True:
        tmp1 = [met_o1, met_o2, met_o3, met_o4]
        nn = 0
        lix = [[0, 0]]
        for tt in tmp1:
            tmp2 = [txt for ll, txt in enumerate(xnames) if
                    txt in tt or txt + "_1" in tt or txt + "_2" in tt]
            nn += len(tmp2)
            if len(tmp2) > 0:
                lix += [[nn, nn]]
            del tmp2
        liy = [[0, len(tab)]] * len(lix)
        lic, lis = ["k"] * len(lix), ["-"] * len(lix)
        for lc, ls, lx, ly in zip(lic, lis, lix, liy):
            line = Line2D(lx, ly, c=lc, lw=7, ls=ls, zorder=10)
            line.set_clip_on(False)
            ax.add_line(line)
        nn = 0
        lic, lix = list(), list()
        for uu, tt in enumerate(tmp1):
            tmp2 = [txt for ll, txt in enumerate(xnames) if
                    txt in tt or txt + "_1" in tt or txt + "_2" in tt]
            if len(tmp2) > 0:
                if uu == 0:
                    cc = "yellowgreen"
                elif uu == 1:
                    cc = "plum"
                elif uu == 2:
                    cc = "gold"
                else:
                    cc = "turquoise"
                lic += [cc, cc]
                if nn > 0:
                    lix += [[nn + 0.2, nn + len(tmp2)], [nn + 0.2, nn + len(tmp2)]]
                else:
                    lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                nn += len(tmp2)
                del cc
            del tmp2
        liy = [[len(tab), len(tab)], [0, 0]] * int(float(len(lix)) / 2)
        lis = ["-"] * len(lix)
        for mm, (lc, ls, lx, ly) in enumerate(zip(lic, lis, lix, liy)):
            if mm < 2:
                line = Line2D([lx[0] + 0.05, lx[1]], ly, c=lc, lw=10, ls=ls, zorder=10)
            elif mm > len(lis) - 3:
                line = Line2D([lx[0], lx[1] - 0.05], ly, c=lc, lw=10, ls=ls, zorder=10)
            else:
                line = Line2D(lx, ly, c=lc, lw=10, ls=ls, zorder=10)
            line.set_clip_on(False)
            ax.add_line(line)
    # y axis
    ticks = [ii + 0.5 for ii in range(len(ynames))]
    ax.set_yticks(ticks)
    if index != 0:
        ax.set_yticklabels([""] * len(ticks))
    else:
        ax.text(-5 * dx, -1 * dy, my_text, fontsize=25, ha='right', va='top', transform=ax.transAxes)
        ax.tick_params(axis="y", labelsize=20)
        ax.set_yticklabels(ynames)
    ax.yaxis.set_label_coords(-20 * dx, 0.5)
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
                    plt.text(jj + 0.5, ii + 0.5, str(round(tab[ii, jj], 1)), fontsize=10, ha='center', va='center')
    return


def multiportraitplot_grouped(tab, name_plot, xnames, ynames, title='portraitplot', write_metrics=False, my_text="",
                              levels=None, cname=False, chigh=False, cfram=False):
    if levels is None:
        levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    fontdict = {"fontsize": 40, "fontweight": "bold"}
    nbr_space = 2
    nbrc = sum([len(tab[0][ii][0]) for ii in range(len(tab[0]))]) + (len(tab[0]) - 1) * nbr_space
    nbrl = len(tab[0][0]) + 1
    fig = plt.figure(0, figsize=(0.5 * nbrc, 0.5 * nbrl))
    gs = GridSpec(nbrl, nbrc)
    countl = 0
    for pp in range(len(tab)):
        # colorbar
        cmap, norm = my_colorbar(mini=min(levels[pp]), maxi=max(levels[pp]), pptype=pp+1)
        countc = 0
        for kk, tmp in enumerate(tab[pp]):
            ax = plt.subplot(gs[countl: countl + 5, countc: countc + len(tmp[0])])
            tmp1 = tmp[-5:]
            cs = ax.pcolormesh(tmp1, cmap=cmap, norm=norm)
            my_pp(ax, kk, tmp1, title[pp][kk], xnames[kk], ynames[-5:], cname=cname, chigh=chigh, cfram=cfram,
                  cfram_top=True, plot_title=True, write_metrics=write_metrics)
            y2_cs = ax.get_position().y1
            if ' ':
                # dots to show the missing models
                ax = plt.subplot(gs[countl + 5: countl + 8, countc: countc + len(tmp[0])])
                ax.axis("off")
                x1, x2 = ax.get_xlim()
                y1, y2 = ax.get_ylim()
                l1 = [(x2-x1)/2.] * 3
                l2 = [tt*(y2-y1)/100. for tt in range(40, 61, 10)]#[tt*(y2-y1)/3. for tt in range(1,3)] + [(y2-y1)/2.]
                ax.scatter(l1, l2, marker="s", s=80, facecolors="k", edgecolors="k", clip_on=False)
                del l1, l2, x1, x2, y1, y2
            ax = plt.subplot(gs[countl + 8: countl + 13, countc: countc + len(tmp[0])])
            tmp1 = tmp[3:8]
            cs = ax.pcolormesh(tmp1, cmap=cmap, norm=norm)
            plot_xname = True if pp+1 == len(tab) else False
            my_tt = deepcopy(my_text) if pp+1 == len(tab) else ""
            my_pp(ax, kk, tmp1, title[pp][kk], xnames[kk], ynames[3:8], cname=cname, chigh=chigh, cfram=cfram,
                  cfram_bot=True, my_text=my_tt, plot_xname=plot_xname, write_metrics=write_metrics)
            x2_cs = ax.get_position().x1
            y1_cs = ax.get_position().y0
            countc += len(tmp[0]) + nbr_space
            del ax, my_tt, plot_xname, tmp1
        countl += 20
        # color bar
        cax = plt.axes([x2_cs + 0.03, y1_cs, 0.02, y2_cs - y1_cs])
        if pp+1 == 1:
            cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels[pp], pad=0.05, extend='max',
                                aspect=40)
            cbar.ax.set_yticklabels(levels[pp], fontdict=fontdict)
            cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.00), fontsize=40, weight="bold",
                         rotation="vertical", ha="center", va='top',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.5, "further from ref.", fontsize=40, rotation="vertical", ha="left", va='center',
                     weight="bold")
        elif pp+1 == 2:
            cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels[pp], pad=0.05, extend='max',
                                aspect=40)
            tmp = [str(tt)+" $\sigma$" if tt!=0 else str(tt) for tt in levels[pp]]
            cbar.ax.set_yticklabels(tmp, fontdict=fontdict)
            cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.00), fontsize=40, weight="bold",
                         rotation="vertical", ha="center", va='top',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.5, "further from ref.", fontsize=40, rotation="vertical", ha="left", va='center',
                     weight="bold")
            del tmp
        elif pp+1 == 3:
            cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels[pp], pad=0.05, extend='both',
                                aspect=40)
            tmp = [str(tt) if tt!=0 else "MMM" for tt in levels[pp]]
            cbar.ax.set_yticklabels(tmp, fontdict=fontdict)
            cax.annotate("", xy=(5, 0.06), xycoords='axes fraction', xytext=(5, 0.43), fontsize=40, weight="bold",
                         rotation="vertical", ha="center", va='bottom',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.45, "closer", fontsize=40, rotation="vertical", ha="left", va='top', weight="bold")
            cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.57), fontsize=40,
                         weight="bold", rotation="vertical", ha="center", va='top',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.55, "further", fontsize=40, rotation="vertical", ha="left", va='bottom', weight="bold")
            del tmp
        elif pp+1 == 4:
            cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels[pp], pad=0.05, extend='both',
                                aspect=40)
            tmp = [str(tt)+" $\sigma$" if tt!=0 else "MMM" for tt in levels[pp]]
            cbar.ax.set_yticklabels(tmp, fontdict=fontdict)
            cax.annotate("", xy=(5, 0.06), xycoords='axes fraction', xytext=(5, 0.43), fontsize=40, weight="bold",
                         rotation="vertical", ha="center", va='bottom',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.45, "closer", fontsize=40, rotation="vertical", ha="left", va='top', weight="bold")
            cax.annotate("", xy=(5, 0.94), xycoords='axes fraction', xytext=(5, 0.57), fontsize=40, weight="bold",
                         rotation="vertical", ha="center", va='top',
                         arrowprops=dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0))
            cax.text(6.5, 0.55, "further", fontsize=40, rotation="vertical", ha="left", va='bottom', weight="bold")
            del tmp
        del cax, cmap, countc, cs, norm, x2_cs, y1_cs, y2_cs
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
    return filename_js, data["RESULTS"]["model"]


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
tab_all, tab_ed1, tab_ed2, tab_ed3, y_names_all = list(), list(), list(), list(), list()
for mc in metric_collections:
    dict1 = dict()
    dict_mod = dict()
    for proj in list_project:
        # open and read json file
        filename_js, data_json = read_data(proj, mc)
        # read metrics
        list_models = sorted(data_json.keys(), key=lambda v: v.upper())
        print(mc, proj, len(list_models))
        dict_mod[proj] = list_models
        for mod in list_models:
            data_mod = data_json[mod][data_json[mod].keys()[0]]["value"]
            list_metrics = sorted(data_mod.keys(), key=lambda v: v.upper())
            if reduced_set is True:
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
                if data_met[my_ref]["value"] is None:
                    dict2[met] = 1e20
                else:
                    dict2[met] = data_met[my_ref]["value"]
                del data_met, my_ref, list_ref, ref
            dict1[mod] = dict2
            del data_mod, dict2
        del data_json, list_models
    # other references
    path_ref = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Test"
    lname = "yplanton_" + mc + "_observation.json"
    filename_js = list(GLOBiglob(OSpath__join(path_ref, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    data_json = data["RESULTS"]["model"]
    if mc == "ENSO_tel":
        lname = "yplanton_" + mc.replace("ENSO", "test") + "_observation.json"
        filename_js = list(GLOBiglob(OSpath__join(path_ref, lname)))[0]
        with open(filename_js) as ff:
            data = json.load(ff)
        ff.close()
        data_json2 = data["RESULTS"]["model"]
    else:
        data_json2 = dict()
    del data, ff, filename_js, lname, path_ref
    list_ref2, list_ref3, list_ref4 = dict(), dict(), dict()
    for met in list_metrics:
        if "MapDjf" in met or "MapJja" in met:
            ref = get_reference(mc.replace("ENSO", "test"), met)
        else:
            ref = get_reference(mc, met)
        # era
        if "TauxSsh" in met or "SshSst" in met:
            tab = data_json["ERA-Interim_SODA3.4.2"]["r1i1p1"]["value"][met]["metric"]
        elif "Ssh" in met:
            tab = data_json["SODA3.4.2"]["r1i1p1"]["value"][met]["metric"]
        elif "MapDjf" in met or "MapJja" in met:
            try:
                tab = data_json2["ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
            except:
                tab = data_json2["ERA-Interim_ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
        else:
            try:    tab = data_json["ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
            except: tab = data_json["ERA-Interim_ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
        if "SstMap" in met:
            list_ref2[met] = 0
        else:
            try:    list_ref2[met] = tab[ref]["value"]
            except: list_ref2[met] = 1e20
        del tab
        # ncep2
        if "TauxSsh" in met or "SshSst" in met:
            tab = data_json["NCEP2_GODAS"]["r1i1p1"]["value"][met]["metric"]
        elif "Ssh" in met:
            tab = data_json["GODAS"]["r1i1p1"]["value"][met]["metric"]
        elif "MapDjf" in met or "MapJja" in met:
            try:
                tab = data_json2["NCEP2"]["r1i1p1"]["value"][met]["metric"]
            except:
                tab = data_json2["NCEP2_NCEP2"]["r1i1p1"]["value"][met]["metric"]
        else:
            try:    tab = data_json["NCEP2"]["r1i1p1"]["value"][met]["metric"]
            except: tab = data_json["NCEP2_NCEP2"]["r1i1p1"]["value"][met]["metric"]
        try:    list_ref3[met] = tab[ref]["value"]
        except: list_ref3[met] = 1e20
        del tab
        # 20cr
        if "MapDjf" in met or "MapJja" in met:
            try:    tab = data_json2["20CRv2"]["r1i1p1"]["value"][met]["metric"]
            except: tab = data_json2["20CRv2_20CRv2"]["r1i1p1"]["value"][met]["metric"]
        elif "Ssh" not in met:
            try:    tab = data_json["20CRv2"]["r1i1p1"]["value"][met]["metric"]
            except: tab = data_json["20CRv2_20CRv2"]["r1i1p1"]["value"][met]["metric"]
        try:    list_ref4[met] = tab[ref]["value"]
        except: list_ref4[met] = 1e20
        try: del tab
        except: pass
        del ref
    # shape data
    my_models = sorted(dict1.keys(), key=lambda v: v.upper())
    my_models = [met for met in mod_order if met in my_models]
    tmp = list(set(dict1.keys()) - set(my_models))
    if len(tmp) > 0:
        my_models += tmp
    for mod in my_models:
        try: dict1[mod]
        except: pass
        else:
            try:    my_metrics
            except: my_metrics = dict1[mod].keys()
            else:   my_metrics += dict1[mod].keys()
    my_metrics = sorted(list(set(my_metrics)), key=lambda v: v.upper())
    my_metrics = [met for met in met_order if met in my_metrics]
    my_models = list(reversed(my_models))
    plus = 6
    tab = NUMPYma__zeros((len(my_models) + plus, len(my_metrics)))
    for ii, mod in enumerate(my_models):
        for jj, met in enumerate(my_metrics):
            try: dict1[mod][met]
            except: tab[ii + plus, jj] = 1e20
            else: tab[ii + plus, jj] = dict1[mod][met]
    tab = NUMPYma__masked_invalid(tab)
    tab = NUMPYmasked_where(tab == 1e20, tab)
    tab1, tab2, tab3 = deepcopy(tab), deepcopy(tab), deepcopy(tab)
    for jj, met in enumerate(my_metrics):
        tmp = tab[plus:, jj].compressed()
        med = float(SCIPYstats__scoreatpercentile(list(tmp), 50))
        mea = float(NUMPYmean(tmp))
        mini = float(min(tmp))
        std = float(NUMPYstd(tmp))
        tmp = [tab[ii + plus, jj] for ii, mod in enumerate(my_models) if mod in dict_mod["cmip5"]]
        tmp = NUMPYma__array(tmp).compressed()
        sel1 = float(NUMPYmean(tmp))
        tmp = [tab[ii + plus, jj] for ii, mod in enumerate(my_models) if mod in dict_mod["cmip6"]]
        tmp = NUMPYma__array(tmp).compressed()
        sel2 = float(NUMPYmean(tmp))
        #list_v = [0, sel2, sel1]
        list_v = [list_ref4[met], list_ref3[met], list_ref2[met], 0, sel2, sel1]
        for ii, val in enumerate(list_v):
            tab[ii, jj] = val
        tab[:, jj] = (tab[:, jj] - mea) / std
        # educational v1
        for ii, val in enumerate(list_v):
            tab1[ii, jj] = val
        # educational v2
        for ii, val in enumerate(list_v):
            tab2[ii, jj] = val
        tab2[:, jj] = tab2[:, jj] / std
        # educational v3
        for ii, val in enumerate(list_v):
            tab3[ii, jj] = val
        tab3[:, jj] = tab3[:, jj] - mea
        del mea, med, sel1, sel2, tmp
    tab = NUMPYma__masked_invalid(tab)
    tab = NUMPYmasked_where(tab > 1e3, tab)
    x_names = ["(20CRv2)", "(NCEP2)", "(ERA-Interim)", "(reference)", "(CMIP6)", "(CMIP5)"] +\
              ["* " + mod if mod in dict_mod["cmip6"] else mod for mod in my_models]
    tab_all.append(tab)
    tab_ed1.append(tab1)
    tab_ed2.append(tab2)
    tab_ed3.append(tab3)
    if reduced_set is True:
        y_names_all.append([met.replace("_2", "") for met in my_metrics])
    else:
        y_names_all.append(my_metrics)
    if mc == metric_collections[0]:
        x_names_all = deepcopy(x_names)
    del dict1, my_metrics, my_models, tab, x_names

# all portraitplots in one plot
if ' ':
    numbering = [ii+") " for ii in list(string.ascii_lowercase)]
    # plot
    figure_name =\
        OSpath__join(path_out, "Figure_02_portraitplot_" + str(len(metric_collections)) + "metric_collections_20200430")
    if reduced_set is False:
        figure_name += "_all_metrics"
    title = [numbering[ii] + dict_mc[mc] for ii, mc in enumerate(metric_collections)]
    xlabel = ""
    ylabel = ""
    text = "* = CMIP6\nmodel"
    levels = list(range(-2, 3))
    multiportraitplot(tab_all, figure_name, xlabel, ylabel, y_names_all, x_names_all, title=title, my_text=text,
                      levels=levels, cname=False, chigh=True, cfram=True, pptype=0)
    if reduced_set is True:
        # plot education v1
        figure_name = OSpath__join(path_out, "Figure_02_education_v1_20200430")
        levels = list(range(0, 5))
        multiportraitplot(tab_ed1, figure_name, xlabel, ylabel, y_names_all, x_names_all, title=title, my_text=text,
                          levels=levels, cname=False, chigh=True, cfram=True, pptype=1)
        # plot education v2
        figure_name = OSpath__join(path_out, "Figure_02_education_v2_20200430")
        levels = list(range(0, 5))
        multiportraitplot(tab_ed2, figure_name, xlabel, ylabel, y_names_all, x_names_all, title=title, my_text=text,
                          levels=levels, cname=False, chigh=True, cfram=True, pptype=2)
        # plot education v3
        figure_name = OSpath__join(path_out, "Figure_02_education_v3_20200430")
        levels = list(range(-2, 3))
        multiportraitplot(tab_ed3, figure_name, xlabel, ylabel, y_names_all, x_names_all, title=title, my_text=text,
                          levels=levels, cname=False, chigh=True, cfram=True, pptype=3)
        # plot education v4
        figure_name = OSpath__join(path_out, "Figure_02_education_v4_20200430")
        levels = list(range(-2, 3))
        multiportraitplot(tab_all, figure_name, xlabel, ylabel, y_names_all, x_names_all, title=title, my_text=text,
                          levels=levels, cname=False, chigh=True, cfram=True, pptype=4)
        # plot education grouped
        figure_name = OSpath__join(path_out, "Figure_A4_education_20200430")
        levels = [list(range(0, 5))] * 2 + [list(range(-2, 3))] * 2
        title = [[numbering[(jj * 3) + ii] + dict_mc[mc] for ii, mc in enumerate(metric_collections)] for jj in
                 range(4)]
        multiportraitplot_grouped([tab_ed1, tab_ed2, tab_ed3, tab_all], figure_name, y_names_all, x_names_all,
                                  title=title,
                                  my_text=text, levels=levels, cname=False, chigh=True, cfram=True)
    del figure_name, levels, text, title, xlabel, ylabel
