# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots of the correlation inter metrics or inter models
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from numpy import mean as NUMPYmean
from numpy import std as NUMPYstd
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoPlotLib import plot_param


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
experiment = "historical"  # "piControl" #
member = "r1i1p1"
list_project = ["cmip6", "cmip5"]
my_project = ["select8", "CMIP"]
big_ensemble = True  # False  #
reduced_set = True  # False  #

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_10_report"
path_in = OSpath__join(path_main, "Data_grouped")
# path_out = OSpath__join(path_main, "Plots_v5")
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_12_09_AGU/Poster"
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/v02"

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


def parallelplot(tab, obs_val, ind=0, labels="", title="", uppertitle="", colors=None, legend=None, plot_legend=True,
                 ylim=None):
    ax = plt.subplot(gs[0, ind])
    if colors is None:
        colors = ["dodgerblue", "k"]
    elif isinstance(colors, dict):
        colors = [colors[leg] for leg in legend]
    col = colors[0]
    boxproperties = {
        "boxprops": dict(linestyle="-", linewidth=3, color=col),
        "capprops": dict(linestyle="-", linewidth=3, color=col),
        "flierprops": dict(marker="o", markersize=6.0, markeredgecolor=col, markerfacecolor=col, markeredgewidth=0),
        "meanprops":  dict(marker="D", markersize=15.0, markeredgecolor=col, markerfacecolor=col, markeredgewidth=0),
        "medianprops": dict(linestyle="-", linewidth=0, color=col),
        "whiskerprops": dict(linestyle="-", linewidth=3, color=col),
    }
    ax.set_title(title, fontsize=20, y=1.05, loc='center')
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlim([-1, 1])
    ax.tick_params(axis="x", labelsize=15, rotation=90)
    ax.set_ylim(ylim)
    ax.set_yticks(ylim)
    ax.set_yticklabels(ylim)
    ax.tick_params(axis="y", labelsize=15)
    ax.boxplot(tab, positions=[0], whis=[5, 95], widths=0.7, labels=[labels], showmeans=True, showfliers=True,
               zorder=4, **boxproperties)
    col = colors[1]
    if obs_val is not None and obs_val != 1e20:
        ax.plot([0.2], obs_val, ls='None', marker="<", mec=col, mew=1, mfc=col, markersize=18, zorder=5, clip_on=False)
    if plot_legend is True:
        ax = plt.subplot(gs[0, ind + 1: ind + 4])
        ax.axis('off')
        coltmp = list(reversed(colors))
        legtmp = list(reversed(legend))
        martmp = ["<", "D"]
        msitmp = [20, 15]
        lines = [Line2D([0], [0], marker=martmp[kk], c="w", markerfacecolor=coltmp[kk], markersize=msitmp[kk])
                 for kk, leg in enumerate(legtmp)]
        ax.legend(lines, legtmp, fontsize=18, bbox_to_anchor=(0.5, 0.5), loc="center", ncol=1)
    ax.text(0.5, 1.3, uppertitle, fontsize=30, transform=ax.transAxes, horizontalalignment="center",
            verticalalignment="bottom")
    return


def transition(text, yy=0.6, ind=0, nint=2, uppertitle=""):
    ax = plt.subplot(gs[0, ind: ind + nint])
    ax.axis("off")
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    dictarrow = dict(facecolor="k", width=5, headwidth=20, headlength=20)
    ax.annotate("", xy=(0.9, 0.45), xycoords='data', xytext=(0.1, 0.45), arrowprops=dictarrow)
    fondict = dict(horizontalalignment="center", verticalalignment="bottom")
    ax.text(0.5, yy, text, fontsize=19, **fondict)
    ax.text(0.5, 1.3, uppertitle, fontsize=30, transform=ax.transAxes, **fondict)
    return


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
dict_dia, dict_met, dict_uni = dict(), dict(), dict()
for proj in list_project:
    if big_ensemble is not True or (big_ensemble is True and proj == list_project[0]):
        dict_di, dict_me = dict(), dict()
    for mc in metric_collection:
        # get metrics list
        list_metrics = sorted(defCollection(mc)['metrics_list'].keys(), key=lambda v: v.upper())
        if reduced_set is True:
            if mc == "ENSO_perf":
                to_remove = ['BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse', 'EnsoTauxTsRmse',
                             'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
                             'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDur_1',
                             'NinoSstDur_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1',
                             'NinoSstTsRmse_2', "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
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
        # !!!!! temporary: start !!!!!
        # # ssh metrics are not computed yet (ask jiwoo)
        # list_metrics2 = deepcopy(list_metrics)
        # for met in list_metrics2:
        #     if "Ssh" in met:
        #         while met in list_metrics:
        #             list_metrics.remove(met)
        # del list_metrics2
        # slp metrics are wrong (error in observation?)
        list_metrics2 = deepcopy(list_metrics)
        for met in list_metrics2:
            if "Slp" in met:
                while met in list_metrics:
                    list_metrics.remove(met)
        del list_metrics2
        # !!!!! temporary: end !!!!!
        # read json
        lpath = OSpath__join(path_in, proj + "/" + experiment + "/" + mc)
        # lpath = deepcopy(path_in)
        lname = proj + "_" + experiment + "_" + mc + "_v2019????.json"
        # lname = proj + "_" + experiment + "_" + mc + "_v2019*.json"
        filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
        with open(filename_js) as ff:
            data = json.load(ff)
        ff.close()
        # units
        if proj == "cmip5":
            dicttmp = data["RESULTS"]["model"]["CNRM-CM5"]["metadata"]["metrics"]
            keystmp = dicttmp.keys()
            for met in keystmp:
                uni = dicttmp[met]["metric"]["units"] if ("Rmse" in met or "Corr" in met)\
                    else dicttmp[met]["diagnostic"]["units"]
                dict_uni[met] = uni.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
            del dicttmp, keystmp
        list_models = sorted(data["RESULTS"]["model"].keys(), key=lambda v: v.upper())
        # read metrics
        dict11, dict12 = dict(), dict()
        for mod in list_models:
            dict21, dict22 = dict(), dict()
            for met in list_metrics:
                if proj == "cmip5":
                    dicttmp1 = data["RESULTS"]["model"]["CNRM-CM5"]["value"][met]["metric"]
                else:
                    dicttmp1 = data["RESULTS"]["model"]["CNRM-CM6-1"]["value"][met]["metric"]
                keystmp = dicttmp1.keys()
                # metric
                try: data["RESULTS"]["model"][mod]["value"][met]["metric"]
                except: tmp = dict((key1, 1e20) for key1 in keystmp)
                else:
                    dicttmp2 = data["RESULTS"]["model"][mod]["value"][met]["metric"]
                    tmp = dict((key1, 1e20 if (key1 not in dicttmp2.keys() or dicttmp2[key1]["value"] is None) else
                                dicttmp2[key1]["value"]) for key1 in keystmp)
                dict21[met] = tmp
                # diagnostic
                if proj == "cmip5":
                    dicttmp1 = data["RESULTS"]["model"]["CNRM-CM5"]["value"][met]["diagnostic"]
                else:
                    dicttmp1 = data["RESULTS"]["model"]["CNRM-CM6-1"]["value"][met]["diagnostic"]
                try: data["RESULTS"]["model"][mod]["value"][met]["diagnostic"]
                except: tmp = dict((key1, 1e20 if key1 == mod else dicttmp1[key1]["value"]) for key1 in keystmp + [mod])
                else:
                    dicttmp2 = data["RESULTS"]["model"][mod]["value"][met]["diagnostic"]
                    tmp = dict((key1, dicttmp1[key1]["value"] if (key1 not in dicttmp2.keys() and key1 != mod) else
                                (1e20 if dicttmp2[key1]["value"] is None else dicttmp2[key1]["value"]))
                               for key1 in keystmp + [mod])
                dict22[met] = tmp
                del dicttmp1, dicttmp2, keystmp, tmp
            dict11[mod] = dict21
            dict12[mod] = dict22
            del dict21, dict22
        # save in common dictionary
        dict_di = common_save(dict12, dict_out=dict_di)
        dict_me = common_save(dict11, dict_out=dict_me)
        del data, dict11, dict12, ff, filename_js, list_metrics, list_models, lname, lpath
    if big_ensemble is not True:
        dict_dia[proj] = dict_di
        dict_met[proj] = dict_me
        del dict_di, dict_me
if big_ensemble is True:
    dict_dia = deepcopy(dict_di)
    dict_met = deepcopy(dict_me)
    del dict_di, dict_me


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
    # !!!!! temporary: start !!!!!
    # some models are not used in all metric collection (ask jiwoo)
    for key2 in list_metrics:
        for key1 in lev1:
            if key2 not in dict_dia[key1].keys():
                dict_dia[key1][key2] = dict((ref, 1e20) for ref in dict_dia["ACCESS1-0"][key2].keys())
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


# ---------------------------------------------------#
# Plot
# ---------------------------------------------------#
if ' ':
    figure_name = OSpath__join(path_out, "creating_portraitplot")
    mimadict = {
        "BiasPrLatRmse": [0, 5], "BiasPrLonRmse": [0, 3], "BiasSstLonRmse": [0, 3], "EnsoSstLonRmse": [0.0, 0.6],
        "EnsoAmpl": [0.4, 1.6], "EnsodSstOce": [0, 3], "EnsoSeasonality": [0.5, 2.5], "EnsoSstSkew": [-0.8, 0.8],
        "NinoSstDiversity": [0, 70], "EnsoFbSshSst": [0.1, 0.3], "EnsoFbSstTaux": [0, 15], "EnsoFbSstThf": [-25, 5],
        "EnsoFbTauxSsh": [0.1, 0.4], "EnsoPrMapRmse": [0.00, 0.35], "EnsoPrMapCorr": [0, 1],
    }
    # my_col = ["dodgerblue", "k"]
    my_col = ["forestgreen", "k"]
    my_leg = ["CMIP", "obs"]
    dict_out = deepcopy(dict_met)
    # list1 = ["BiasSstLonRmse", "EnsoFbSstTaux"]
    # list2 = [dict_met, dict_dia]
    list1 = ["EnsoPrMapRmse", "EnsoFbSstTaux", "EnsoPrMapCorr"]
    list2 = [dict_met, dict_dia, dict_met]
    numtmp = ["a) ", "b) ", "c) "]
    nbr = sum([4 if "Rmse" in met else 7 for ii, met in enumerate(list1)]) + 3 * len(list1)
    fig = plt.figure(0, figsize=(1 * nbr, 4))
    gs = GridSpec(1, nbr)
    nbr = 0
    for ii, met in enumerate(list1):
        ref = get_ref(met)
        # diag
        if "Rmse" not in met:
            if "Corr" in met:
                tname = "correlation"
                xlabel = deepcopy(met)
                my_text = "1 - corr"
                val = [list2[ii][mod][met][ref] for mod in list2[ii].keys() if list2[ii][mod][met][ref] != 1e20]
                obs = 1
            else:
                tname = "statistic"
                xlabel = met + "\n(" + dict_uni[met] + ")"
                my_text = r'abs$\left(\frac{mod - obs}{obs}\right)$'
                val = [list2[ii][mod][met][mod] for mod in list2[ii].keys() if list2[ii][mod][met][ref] != 1e20]
                obs = list2[ii][list2[ii].keys()[0]][met][ref]
            my_ylim = mimadict[met]
            parallelplot(val, obs, ind=nbr, labels=xlabel, title=tname, colors=my_col, legend=my_leg,
                         plot_legend=False, ylim=my_ylim)
            nbr += 1
            transition(my_text, yy=0.55, ind=nbr, nint=2)
            nbr += 2
            del my_text, my_ylim, tname, xlabel
        # metric
        if "Rmse" in met:
            xlabel = met + "\n(" + dict_uni[met] + ")"
            my_ylim = mimadict[met]
            val = [list2[ii][mod][met][ref] for mod in list2[ii].keys() if list2[ii][mod][met][ref] != 1e20]
            uptitle = ""
        elif "Corr" in met:
            xlabel = "error"
            val = [1 - jj for jj in val]
            my_ylim = [0.0, 0.7]
            uptitle = numtmp[ii] + "correlation based metrics"
        else:
            xlabel = "% of error"
            val = [abs(100. * (jj - obs) / obs) for jj in val]
            my_ylim = [0, 80]
            uptitle = numtmp[ii] + "statistic based metrics"
        tname = "metric"
        obs = 0
        my_text =r'$\frac{metric - MMM}{\sigma}$'
        parallelplot(val, obs, ind=nbr, labels=xlabel, title=tname, uppertitle=uptitle, colors=my_col, legend=my_leg,
                     plot_legend=False, ylim=my_ylim)
        nbr += 1
        if "Rmse" in met:
            uptitle = numtmp[ii] + "RMSE based metrics"
        else:
            uptitle = ""
        transition(my_text, yy=0.55, ind=nbr, nint=2, uppertitle=uptitle)
        nbr += 2
        del my_text, my_ylim, tname, uptitle, xlabel
        # normalized
        tname = "portrait plot"
        xlabel = "displayed\nmetric"
        mmm = float(NUMPYmean(val))
        std = float(NUMPYstd(val))
        val = [float(jj - mmm) / std for jj in val]
        my_ylim = [-2.5, 2.5]
        plotleg = True if ii == 0 else False
        obs = 1e20  # float(obs - mmm) / std
        parallelplot(val, obs, ind=nbr, labels=xlabel, title=tname, colors=my_col, legend=my_leg,
                     plot_legend=plotleg, ylim=my_ylim)
        nbr += 4
        del mmm, my_ylim, obs, std, tname, xlabel
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()



