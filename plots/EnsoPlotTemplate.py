# -*- coding:UTF-8 -*-
import cmocean
from copy import deepcopy
from math import ceil as MATHceil
from math import floor as MATHfloor
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from numpy import arange as NUMPYarange
from numpy import array as NUMPYarray
from numpy import linspace as NUMPYlinspace
from numpy import meshgrid as NUMPYmeshgrid
from numpy.ma import masked_where as NUMPYmasked_where

# ENSO_metrics functions
from EnsoMetrics.EnsoCollectionsLib import ReferenceRegions
from .EnsoPlotToolsLib import create_labels, create_levels, format_metric, minimaxi, minmax_plot, my_average,\
    my_bootstrap, my_legend, my_mask, my_mask_map, read_diag, read_var, return_metrics_type, shading_levels

colors_sup = ["r", "lime", "peru", "gold", "forestgreen", "sienna", "gold"]
dict_col = {"REF": "k", "CMIP": "forestgreen", "CMIP3": "orange", "CMIP5": "dodgerblue", "CMIP6": "r"}

met_names = {
    "BiasPrLatRmse": "double_ITCZ_bias", "BiasPrLonRmse": "eq_PR_bias",
    "BiasSshLatRmse": "lat_SSH_bias", "BiasSshLonRmse": "eq_SSH_bias",
    "BiasSstLatRmse": "lat_SST_bias", "BiasSstLonRmse": "eq_SST_bias",
    "BiasTauxLatRmse": "lat_Taux_bias", "BiasTauxLonRmse": "eq_Taux_bias",
    "SeasonalPrLatRmse": "double_ITCZ_sea_cycle", "SeasonalPrLonRmse": "eq_PR_sea_cycle",
    "SeasonalSshLatRmse": "lat_SSH_sea_cycle", "SeasonalSshLonRmse": "eq_SSH_sea_cycle",
    "SeasonalSstLatRmse": "lat_SST_sea_cycle", "SeasonalSstLonRmse": "eq_SST_sea_cycle",
    "SeasonalTauxLatRmse": "lat_Taux_sea_cycle", "SeasonalTauxLonRmse": "eq_Taux_sea_cycle",
    "EnsoPrLonRmse": "ENSO_pattern_PR", "EnsoSshLonRmse": "ENSO_pattern_SSH", "EnsoSstLonRmse": "ENSO_pattern",
    "EnsoTauxLonRmse": "ENSO_pattern_Taux", "EnsoPrTsRmse": "ENSO_lifecycle_PR", "EnsoSshTsRmse": "ENSO_lifecycle_SSH",
    "EnsoSstTsRmse": "ENSO_lifecycle", "EnsoTauxTsRmse": "ENSO_lifecycle_Taux",
    "EnsoAmpl": "ENSO_amplitude", "EnsoSeasonality": "ENSO_seasonality", "EnsoSstSkew": "ENSO_asymmetry",
    "EnsoDuration": "ENSO_duration", "EnsoSstDiversity": "ENSO_diversity", "EnsoSstDiversity_1": "ENSO_diversity",
    "EnsoSstDiversity_2": "ENSO_diversity", "EnsoPrMapCorr": "Dec_PR_teleconnection_CORR",
    "EnsoPrMapRmse": "Dec_PR_teleconnection", "EnsoPrMapStd": "Dec_PR_teleconnection_STD",
    "EnsoPrMapDjfCorr": "DJF_PR_teleconnection_CORR", "EnsoPrMapDjfRmse": "DJF_PR_teleconnection",
    "EnsoPrMapDjfStd": "DJF_PR_teleconnection_STD", "EnsoPrMapJjaCorr": "JJA_PR_teleconnection_CORR",
    "EnsoPrMapJjaRmse": "JJA_PR_teleconnection", "EnsoPrMapJjaStd": "JJA_PR_teleconnection_STD",
    "EnsoSlpMapRmse": "Dec_SLP_teleconnection", "EnsoSlpMapStd": "Dec_SLP_teleconnection_STD",
    "EnsoSlpMapDjfCorr": "DJF_SLP_teleconnection_CORR", "EnsoSlpMapDjfRmse": "DJF_SLP_teleconnection",
    "EnsoSlpMapDjfStd": "DJF_SLP_teleconnection_STD", "EnsoSlpMapJjaCorr": "JJA_SLP_teleconnection_CORR",
    "EnsoSlpMapJjaRmse": "JJA_SLP_teleconnection", "EnsoSlpMapJjaStd": "JJA_SLP_teleconnection_STD",
    "EnsoSstMapRmse": "Dec_TS_teleconnection", "EnsoSstMapStd": "Dec_TS_teleconnection_STD",
    "EnsoSstMapDjfCorr": "DJF_TS_teleconnection_CORR", "EnsoSstMapDjfRmse": "DJF_TS_teleconnection",
    "EnsoSstMapDjfStd": "DJF_TS_teleconnection_STD", "EnsoSstMapJjaCorr": "JJA_TS_teleconnection_CORR",
    "EnsoSstMapJjaRmse": "JJA_TS_teleconnection", "EnsoSstMapJjaStd": "JJA_TS_teleconnection_STD",
    "EnsoFbSstTaux": "SST-Taux_feedback", "EnsoFbTauxSsh": "Taux-SSH_feedback", "EnsoFbSshSst": "SSH-SST_feedback",
    "EnsoFbSstThf": "SST-NHF_feedback", "EnsodSstOce": "ocean_driven_SST", "EnsodSstOce_1": "ocean_driven_SST",
    "EnsodSstOce_2": "ocean_driven_SST"}

mod_nicknames = ["CMIP5", "CMIP6"]

article_fig = False  # True
plot_for_wiki = False  # True


def cmip_boxplot(dict_param, dict_values, units, reference, val_type, my_text, figure_name):
    # concatenate every mips
    keys = sorted([kk for kk in list(dict_values.keys()) if kk != "ref"], key=lambda v: v.upper())
    vall = list()
    for kk in keys:
        vall += dict_values[kk][reference] if val_type == "metric" else dict_values[kk]
    if len(vall) > 0:
        # add every mips in the list
        legend = ["CMIP"] + [kk.upper() for kk in keys]
        colors = [dict_col[kk] for kk in legend]
        vall = [vall] + [dict_values[kk][reference] if val_type == "metric" else dict_values[kk] for kk in keys]
        # number of 'columns' (boxplot, markers, etc)
        nbrc = len(vall)
        if val_type == "metric":
            nbrc += int(round(len(keys)*(len(keys) - 1) / 2.))
        # plot param
        title = "Metric values" if val_type == "metric" else "Scalar diagnostic values"
        yname = dict_param["title"][0] if isinstance(dict_param["title"], list) is True else dict_param["title"]
        if units != "":
            uni = units.replace("1e2", "10$^2$").replace("1e3", "10$^3$")
            uni = uni.replace("1e-2", "10$^{-2}$").replace("1e-3", "10$^{-3}$")
            uni = uni.replace("m2", "m$^2$")
            yname += " (" + uni + ")"
            del uni
        tmp = vall[0]
        if val_type != "metric":
            tmp += [dict_values["ref"][reference]]
        tick_labels = minmax_plot(tmp)
        mini, maxi = min(tick_labels), max(tick_labels)
        fig, ax = plt.subplots(1, 1, figsize=(max(4, nbrc), 4), sharex="col", sharey="row")
        # title
        ax.set_title(title, fontsize=15, y=1.01, loc="left")
        # # y axis
        ax.set_yticks(tick_labels)
        ax.set_ylim(ymin=mini, ymax=maxi)
        ax.set_ylabel(yname, fontsize=15)
        for axis_tick_label in ax.get_yticklabels():
            axis_tick_label.set_fontsize(12)
        # boxplots
        for ii, (cc, tab) in enumerate(zip(colors, vall)):
            boxproperties = {
                "boxprops": dict(linestyle="-", linewidth=2, color=cc),
                "capprops": dict(linestyle="-", linewidth=2, color=cc),
                "flierprops": dict(marker="o", markersize=5.0, markeredgecolor=cc, markerfacecolor=cc,
                                   markeredgewidth=0),
                "meanprops": dict(marker="D", markersize=10.0, markeredgecolor=cc, markerfacecolor=cc,
                                  markeredgewidth=0),
                "medianprops": dict(linestyle="-", linewidth=2, color=cc),
                "whiskerprops": dict(linestyle="-", linewidth=2, color=cc)}
            tmp = [[1e20, 1e20]] * ii + [tab] + [[1e20, 1e20]] * (nbrc-1-ii)
            ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=True, **boxproperties)
            del boxproperties, tmp
        # bootstrap
        x1, x2 = ax.get_xlim()
        dx = (x2 - x1) / 100.
        y1, y2 = ax.get_ylim()
        dy = (y2 - y1) / 100.
        if val_type == "metric":
            for ii, tab1 in enumerate(vall[1:]):
                for jj, tab2 in enumerate(vall[2+ii:]):
                    cc1, cc2 = colors[ii + 1], colors[ii + jj + 2]
                    bst1, bst2, mea1, mea2 = my_bootstrap(tab1, tab2)
                    fillstyle = "none" if (min(bst2) <= mea1 <= max(bst2)) or (min(bst1) <= mea2 <= max(bst1)) else\
                        "full"
                    xxx = len(vall) + ii + 1
                    ax.plot([xxx], [mea1], markersize=10, color=cc1, marker="D", fillstyle=fillstyle,
                            markeredgecolor=cc1, markeredgewidth=3, zorder=2)
                    ax.plot([xxx], [mea2], markersize=10, color=cc2, marker="D", fillstyle=fillstyle,
                            markeredgecolor=cc2, markeredgewidth=3, zorder=2)
                    ax.add_line(Line2D([xxx - 0.3, xxx + 0.3], [min(bst1), min(bst1)], c=cc1, lw=2, zorder=3))
                    ax.add_line(Line2D([xxx - 0.3, xxx + 0.3], [max(bst1), max(bst1)], c=cc1, lw=2, zorder=3))
                    ax.add_line(Line2D([xxx - 0.05, xxx - 0.05], [min(bst1), max(bst1)], c=cc1, lw=2, zorder=3))
                    ax.add_line(Line2D([xxx - 0.3, xxx + 0.3], [min(bst2), min(bst2)], c=cc2, lw=2, zorder=3))
                    ax.add_line(Line2D([xxx - 0.3, xxx + 0.3], [max(bst2), max(bst2)], c=cc2, lw=2, zorder=3))
                    ax.add_line(Line2D([xxx + 0.05, xxx + 0.05], [min(bst2), max(bst2)], c=cc2, lw=2, zorder=3))
                    del bst1, bst2, cc1, cc2, fillstyle, mea1, mea2, xxx
            ax.plot([x2 + 5 * dx], [y1 + 50 * dy], markersize=8, color="k", marker="D", fillstyle="full",
                    markeredgecolor="k", markeredgewidth=2, clip_on=False, zorder=2)
            ax.text(x2 + 10 * dx, y1 + 50 * dy, r'significantly $\neq$', fontsize=12, color="k", ha="left", va="center")
            ax.plot([x2 + 5 * dx], [y1 + 43 * dy], markersize=8, color="k", marker="D", fillstyle="none",
                    markeredgecolor="k", markeredgewidth=2, clip_on=False, zorder=2)
            ax.text(x2 + 10 * dx, y1 + 43 * dy, r'not $\neq$', fontsize=12, color="k", ha="left", va="center")
            ax.text(x2 + 2 * dx, y1 + 34 * dy, "at the 95% confidence level\n(Monte Carlo sampling method)",
                    fontsize=10, color="k", ha="left", va="center")
        # reference
        if val_type != "metric":
            ax.axhline(dict_values["ref"][reference], c="k", ls="-", lw=4, zorder=1)
        # legend
        if val_type == "metric":
            legend = ["ref: " + reference + " = 0"] + legend
        else:
            legend = ["ref: " + reference] + legend
        lines = [Line2D([0], [0], marker="D", color="w", mfc=cc, markersize=10) for cc in colors]
        lines = [Line2D([0], [0], color=dict_col["REF"], ls="-", lw=4)] + lines
        ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        ax.text(x2 + 2 * dx, y1, my_text, fontsize=12, color="k", ha="left", va="bottom")
        # grid
        ax.grid(ls="--", lw=1, which="major", axis="y")
        # save
        plt.savefig(figure_name, bbox_inches="tight")
        plt.close()
        del ax, colors, fig, legend, lines, maxi, mini, nbrc, tick_labels, title, yname
    else:
        # plot param
        title = "Metric values" if val_type == "metric" else "Scalar diagnostic values"
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), sharex="col", sharey="row")
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        # title
        ax.set_title(title, fontsize=15, y=1.01, loc="left")
        ax.set_xlim(xmin=0, xmax=1)
        ax.set_ylim(ymin=0, ymax=1)
        txt = "No metric\nvalues?" if val_type == "metric" else "No scalar diagnostic\nfor this metric"
        ax.text(0.5, 0.5, txt, fontsize=20, color="k", ha="center", va="center")
        # save
        plt.savefig(figure_name, bbox_inches="tight")
        plt.close()
        del ax, fig, title, txt
    return


def my_boxplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
               metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None,
               regions=None, shading=False, plot_ref=False):
    # get data
    variables = dict_param["varpattern"]
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii + 1), regions[vv].replace("nino", "N"))
    if isinstance(variables, str) is True or isinstance(variables, str) is True:
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values, models2=models2,
                 member=member, shading=shading)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    title = dict_param["title"]
    if isinstance(title, str) is True or isinstance(title, str) is True:
        title = [title] * nbr_panel
    yname = dict_param["yname"]
    if isinstance(yname, str) is True or isinstance(yname, str) is True:
        yname = [yname] * nbr_panel
        one_yaxis = True
    else:
        one_yaxis = False
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    elif isinstance(filename_nc, dict) is True and shading is True:
        units = tab_mod[0][0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    if "legend" in list(dict_param.keys()):
        legend = dict_param["legend"]
    else:
        legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                           shading=shading)
        if plot_for_wiki is True:
            if "CNRM-CM5" in legend[0]:
                legend[0] = "model"
            elif "CNRM-CM5" in legend[1]:
                legend[1] = "model"
        if plot_ref is True:
            legend = [legend[0]]
    if "custom_label" in list(dict_param.keys()):
        custom_label = dict_param["custom_label"]
    else:
        custom_label = None
    nbrl = nbr_panel // 2
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    legco = ["k", "dodgerblue"]
    if isinstance(model, list) is True:
        for ii in range(len(model)-1):
            legco.append(colors_sup[ii])
    lines = [Line2D([0], [0], marker="o", c="w", markerfacecolor=cc, markersize=12) for cc in legco]
    if shading is True:
        tmp1 = list()
        for ii in range(len(tab_mod)):
            tmp2 = list()
            for jj in range(nbr_val):
                tmp3 = list()
                for kk in range(len(tab_mod[ii])):
                    # tmp3.append(tab_mod[ii][kk][jj])
                    tmp3 += list(tab_mod[ii][kk][jj])
                tmp2.append(tmp3)
            tmp1.append(tmp2)
        tab_mod = deepcopy(tmp1)
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        if plot_ref is True:
            tmp = tab_obs
        else:
            tmp = tab_mod + tab_obs
    else:
        tmp = deepcopy(tab_obs)
        if plot_ref is False:
            for ii in range(len(tab_mod)):
                tmp += tab_mod[ii]
    if custom_label is None:
        tick_labels = minmax_plot(tmp, metric=plot_metric)
        mini, maxi = min(tick_labels), max(tick_labels)
    else:
        mini, maxi = minimaxi(tmp)
        tick_labels = list(range(int(MATHfloor(mini)), int(MATHceil(maxi)) + 1))
        tick_labels, label = create_labels(custom_label, tick_labels)
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii // 2, ii % 2]
        # title
        ax.set_title(title[ii], fontsize=15, y=1.01, loc="left")
        # x axis
        for axis_tick_label in ax.get_xticklabels():
            axis_tick_label.set_fontsize(12)
        # y axis
        ax.set_yticks(tick_labels)
        if custom_label is not None:
            ax.set_yticklabels(label)
        ax.set_ylim(ymin=mini, ymax=maxi)
        if (one_yaxis is True and ii % 2 == 0) or one_yaxis is False:
            ylabel = yname[ii]
            if units != "":
                ylabel = ylabel + " (" + units + ")"
            ax.set_ylabel(ylabel, fontsize=15)
        for axis_tick_label in ax.get_yticklabels():
            axis_tick_label.set_fontsize(12)
        # boxplots
        boxproperties = {
            "boxprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "capprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=legco[0], markerfacecolor=legco[0],
                               markeredgewidth=0),
            "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=legco[0], markerfacecolor=legco[0],
                              markeredgewidth=0),
            "medianprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "whiskerprops": dict(linestyle="-", linewidth=2, color=legco[0])}
        if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
            ax.boxplot([tab_obs[ii], [1e20, 1e20]], whis=[5, 95], labels=["", ""], showmeans=True, showfliers=False,
                       **boxproperties)
            boxproperties = {
                "boxprops": dict(linestyle="-", linewidth=2, color=legco[1]),
                "capprops": dict(linestyle="-", linewidth=2, color=legco[1]),
                "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=legco[1], markerfacecolor=legco[1],
                                   markeredgewidth=0),
                "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=legco[1], markerfacecolor=legco[1],
                                  markeredgewidth=0),
                "medianprops": dict(linestyle="-", linewidth=2, color=legco[1]),
                "whiskerprops": dict(linestyle="-", linewidth=2, color=legco[1]),
            }
            if plot_ref is False:
                ax.boxplot([[1e20, 1e20], tab_mod[ii]], whis=[5, 95], labels=["", ""], showmeans=True, showfliers=False,
                           **boxproperties)
            # my text
            if plot_metric is True:
                # relative space
                x1, x2 = ax.get_xlim()
                dx = (x2 - x1) / 100.
                y1, y2 = ax.get_ylim()
                dy = (y2 - y1) / 100.
                txt = format_metric(metric_type, metval, metric_units)
                ax.text(x2 - (2 * dx), y2 - (6 * dy), txt, fontsize=12, color="k", horizontalalignment="right",
                        verticalalignment="center")
            # legend
            if (nbr_panel == 1 and ii == 0) or (nbr_panel != 1 and ii == 1):
                ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        else:
            tmp = [tab_obs[ii]] + [1e20, 1e20] * len(tab_mod)
            ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=False, **boxproperties)
            if plot_ref is False:
                for kk in range(len(tab_mod)):
                    boxproperties = {
                        "boxprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                        "capprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                        "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=legco[kk+1],
                                           markerfacecolor=legco[kk+1], markeredgewidth=0),
                        "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=legco[kk+1],
                                          markerfacecolor=legco[kk+1], markeredgewidth=0),
                        "medianprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                        "whiskerprops": dict(linestyle="-", linewidth=2, color=legco[kk+1])}
                    tmp = [[1e20, 1e20]] * (kk + 1) + [tab_mod[kk][ii]] + [[1e20, 1e20]] * (len(tab_mod) - 1 - kk)
                    ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=False,
                               **boxproperties)
            # legend
            if (nbr_panel == 1 and ii == 0) or (nbr_panel != 1 and ii == 1):
                if plot_metric is True:
                    for jj in range(1, len(legend)):
                        legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj-1]) + " " + metric_units + ")"
                ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        # grid
        ax.grid(linestyle="--", linewidth=1, which="major", axis="y")
        if ii == nbr_panel - 1 and plot_ref is True and (ii+1)%2 == 0 and plot_ref is True:
            x1, x2 = ax.get_xlim()
            dx = (x2 - x1) / 100.
            y1, y2 = ax.get_ylim()
            dy = (y2 - y1) / 100.
            ax.text(x2 + 2 * dx, y2 - 20 * dy, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
    if ii + 1 < (nbrc*nbrl):
        if nbrl == 1 and nbrc != 1:
            ax = axes[-1]
        else:
            ax = axes[-1, ii-1]
        ax.axis("off")
        if plot_ref is True:
            x1, x2 = ax.get_xlim()
            dx = (x2 - x1) / 100.
            y1, y2 = ax.get_ylim()
            ax.text(x1 - 18 * dx, y2, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_curve(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
             metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None,
             regions=None, shading=False, plot_ref=False):
    # get data
    variables = dict_param["varpattern"]
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii+1), regions[vv].replace("nino", "N"))
    if isinstance(variables, str) is True or isinstance(variables, str) is True:
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname =\
        read_var(deepcopy(variables), filename_nc, model, reference, metric_variables, metric_values, models2=models2,
                 member=member, shading=shading)
    if metric_type is not None and (isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True):
        plot_metric = True
    else:
        plot_metric = False
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        axis = list(NUMPYarray(tab_mod[0].coords[list(tab_mod[0].coords.keys())[0]]))
    elif isinstance(filename_nc, dict) is True and shading is True:
        axis = list(NUMPYarray(tab_mod[0][0][0].coords[list(tab_mod[0][0][0].coords.keys())[0]]))
    else:
        axis = list(NUMPYarray(tab_mod[0][0].coords[list(tab_mod[0][0].coords.keys())[0]]))
    # figure initialization
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    elif isinstance(filename_nc, dict) is True and shading is True:
        units = tab_mod[0][0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    if units != "":
        yname = yname + " (" + units + ")"
    if "colors" in list(dict_param.keys()):
        linecolors = dict_param["colors"]
    else:
        linecolors = {"model": ["dodgerblue"], "reference": ["k"]}
    if "linestyles" in list(dict_param.keys()):
        linestyles = dict_param["linestyles"]
    else:
        linestyles = {"model": ["-"], "reference": ["-"]}
    if "legend" in list(dict_param.keys()):
        legend = dict_param["legend"]
    else:
        legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                           shading=shading)
        if plot_for_wiki is True:
            if "CNRM-CM5" in legend[0]:
                legend[0] = "model"
            elif "CNRM-CM5" in legend[1]:
                legend[1] = "model"
    nbr = len(tab_mod[0]) if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True else\
        (len(tab_mod[0][0][0]) if isinstance(filename_nc, dict) is True and shading is True else len(tab_mod[0][0]))
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        tmp = tab_mod + tab_obs
    else:
        tmp = deepcopy(tab_obs)
        for ii in range(len(tab_mod)):
            tmp += tab_mod[ii]
    ytick_labels = minmax_plot(tmp, metric=plot_metric)
    # plot
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        if xname == "months" and nbr == 72:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig, ax = plt.subplots(figsize=(4, 4))
        plot_curve(tab_mod, tab_obs, ax, title, axis, xname, yname, ytick_labels, linecolors, linestyles, metric_type,
                   metval, metric_units, model=model, member=member, obsname=obsname, legend=legend, multimodel=False,
                   plot_metric=plot_metric, shading=shading, plot_ref=plot_ref, method=method)
    else:
        if metric_type is not None:
            plot_metric = True
        if xname == "months" and nbr == 72:
            nbrl = nbr_val
            nbrc = 1
            fig, axes = plt.subplots(nbrl, nbrc, figsize=(8, 4 * nbrl), sharex="col", sharey="row")
        else:
            nbrl = nbr_val // 2
            nbrc = 1 if nbr_val == 1 else 2
            fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
        old_leg = deepcopy(legend)
        for kk in range(len(tab_obs)):
            if isinstance(filename_nc, dict) is True and shading is True:
                tab_tmp = [[tab_mod[jj][ll][kk] for ll in range(len(tab_mod[jj]))] for jj in range(len(tab_mod))]
            else:
                tab_tmp = [tab_mod[jj][kk] for jj in range(len(tab_mod))]
            if nbr_val == 1:
                ax = axes
            elif (nbrl == 1 and nbrc != 1) or (nbrl != 1 and nbrc == 1):
                ax = axes[kk % 2]
            else:
                ax = axes[kk // 2, kk % 2]
            if "legend" in list(dict_param.keys()):
                title_tmp = title + ": " + old_leg[kk]
            else:
                title_tmp = title
            lcol = {"model": ["dodgerblue"], "reference": ["k"]}
            lsty = {"model": ["-"], "reference": ["-"]}
            for jj in range(len(filename_nc) - 1):
                lcol["model"].append(colors_sup[jj])
                lsty["model"].append("-")
            if (nbrc == 2 and kk == 1) or (nbrc == 1 and kk == 0) or nbr_val == 1:
                plot_legend = True
            else:
                plot_legend = False
            plot_curve(tab_tmp, [tab_obs[kk]], ax, title_tmp, axis, xname, yname, ytick_labels, lcol, lsty, metric_type,
                       metval, metric_units, model=model, member=member, obsname=obsname, legend=legend,
                       multimodel=True, plot_metric=plot_metric, plot_legend=plot_legend, shading=shading,
                       plot_ref=plot_ref, method=method)
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_dotplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
               metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None,
               regions=None, shading=False, plot_ref=False):
    # get data
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii + 1), regions[vv].replace("nino", "N"))
    diag_mod, diag_obs, metval, obsname =\
        read_diag(diagnostic_values, metric_values, model, reference, metric_variables, member=member)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization if isinstance(filename_nc, str) is True or isinstance(filename_nc, unicode) is True
    title = dict_param["title"]
    yname = dict_param["yname"]
    if diagnostic_units != "":
        yname = yname + " (" + diagnostic_units + ")"
    if len(yname) < 20:
        for kk in list(regions.keys()):
            if kk in yname.lower():
                yname = regions[kk] + " " + yname
    if "colors" in list(dict_param.keys()):
        mcolors = dict_param["colors"]
    else:
        mcolors = ["k", "dodgerblue"]
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                mcolors.append(colors_sup[ii])
    if "markers" in list(dict_param.keys()):
        markers = dict_param["markers"]
    else:
        markers = ["D", "o"]
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                markers.append("o")
    if "legend" in list(dict_param.keys()):
        legend = dict_param["legend"]
    else:
        legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                           shading=shading)
        if plot_for_wiki is True:
            if "CNRM-CM5" in legend[0]:
                legend[0] = "model"
            elif "CNRM-CM5" in legend[1]:
                legend[1] = "model"
        if plot_ref is True:
            legend = [legend[0]]
    fig, ax = plt.subplots(figsize=(4, 4))
    if plot_ref is True:
        tab = [diag_obs]
    else:
        if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
            tab = [diag_obs, diag_mod]
        else:
            tab = [diag_obs] + diag_mod
    # title
    ax.set_title(title, fontsize=15, y=1.01, loc="left")
    # x axis
    label_ticks = [-0.5] + list(range(len(tab))) + [len(tab) - 0.5]
    label = [""] * len(label_ticks)
    plt.xticks(label_ticks, label)
    plt.xlim(min(label_ticks), max(label_ticks))
    for axis_tick_label in ax.get_xticklabels():
        axis_tick_label.set_fontsize(12)
    # y axis
    tick_labels = minmax_plot(tab, metric=plot_metric)
    plt.yticks(tick_labels, tick_labels)
    plt.ylim(min(tick_labels), max(tick_labels))
    ax.set_ylabel(yname, fontsize=15)
    for axis_tick_label in ax.get_yticklabels():
        axis_tick_label.set_fontsize(12)
    if min(tick_labels) < 0 and max(tick_labels) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # dots
    for ii in range(len(tab)):
        ax.scatter([ii], tab[ii], s=80, c=mcolors[ii], marker=markers[ii], clip_on=False)
    x1, x2 = ax.get_xlim()
    dx = (x2 - x1) / 100.
    y1, y2 = ax.get_ylim()
    dy = (y2 - y1) / 100.
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        # my text
        if plot_metric is True:
            txt = format_metric(metric_type, metval, metric_units)
            ax.text(x2 - (2 * dx), y2 - (6 * dy), txt, fontsize=12, color='k', horizontalalignment='right',
                     verticalalignment='center')
        # legend
        lines = [Line2D([0], [0], marker=markers[kk], c="w", markerfacecolor=mcolors[kk], markersize=12)
                 for kk in range(len(mcolors))]
        ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    else:
        # legend
        lines = [Line2D([0], [0], marker=markers[kk], c="w", markerfacecolor=mcolors[kk], markersize=12)
                 for kk in range(len(mcolors))]
        for jj in range(1, len(legend)):
            legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj - 1]) + " " + metric_units + ")"
        ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    if plot_ref is True:
        ax.text(x2 + 2 * dx, y2 - 20 * dy, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
    plt.grid(linestyle='--', linewidth=1, which='major')
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()
    return


def my_dot_to_box(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
                  metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None,
                  diagnostic_units=None, regions=None, shading=False, plot_ref=False):
    # get data
    diag_mod, diag_obs, metval, obsname = \
        read_diag(diagnostic_values, metric_values, model, reference, metric_variables, shading=shading, member=member)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    title = dict_param["title"]
    if isinstance(title, str) is True or isinstance(title, str) is True:
        title = [title] * nbr_panel
    yname = dict_param["yname"]
    if diagnostic_units != "":
        yname = yname + " (" + diagnostic_units + ")"
    if len(yname) < 20:
        for kk in list(regions.keys()):
            if kk in yname.lower():
                yname = regions[kk] + " " + yname
    if "colors" in list(dict_param.keys()):
        mcolors = dict_param["colors"]
    else:
        mcolors = ["k", "dodgerblue"]
        if isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True:
            for ii in range(len(filename_nc) - 1):
                mcolors.append(colors_sup[ii])
    if "legend" in list(dict_param.keys()):
        legend = dict_param["legend"]
    else:
        legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                           shading=shading)
        if plot_for_wiki is True:
            legend[0] = "model"
    fig, ax = plt.subplots(figsize=(4, 4))
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        tab = [diag_mod]
    else:
        tab = diag_mod
    lines = [Line2D([0], [0], marker="o", c="w", markerfacecolor=cc, markersize=12) for cc in mcolors]
    # x axis
    for axis_tick_label in ax.get_xticklabels():
        axis_tick_label.set_fontsize(12)
    # y axis
    tmp = [diag_obs] + [min(my_mask(tt, remove_masked=True)) for tt in tab] +\
          [max(my_mask(tt, remove_masked=True)) for tt in tab]
    tick_labels = minmax_plot(tmp, metric=plot_metric)
    ax.set_yticks(tick_labels)
    ax.set_ylim(ymin=min(tick_labels), ymax=max(tick_labels))
    ax.set_ylabel(yname, fontsize=15)
    for axis_tick_label in ax.get_yticklabels():
        axis_tick_label.set_fontsize(12)
    # plot
    ax.axhline(diag_obs, color=mcolors[0], linestyle='-', linewidth=2)
    for ii in range(len(tab)):
        # boxplots
        boxproperties = {
            "boxprops": dict(linestyle="-", linewidth=2, color=mcolors[ii+1]),
            "capprops": dict(linestyle="-", linewidth=2, color=mcolors[ii+1]),
            "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=mcolors[ii+1], markerfacecolor=mcolors[ii+1],
                               markeredgewidth=0),
            "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=mcolors[ii+1], markerfacecolor=mcolors[ii+1],
                              markeredgewidth=0),
            "medianprops": dict(linestyle="-", linewidth=2, color=mcolors[ii+1]),
            "whiskerprops": dict(linestyle="-", linewidth=2, color=mcolors[ii+1]),
        }
        tmp = [[1e20, 1e20]] * ii + [my_mask(tab[ii], remove_masked=True)] + [[1e20, 1e20]] * (len(tab) - 1 - ii)
        ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=True, **boxproperties)
    # legend
    if plot_metric is True:
        for jj in range(1, len(legend)):
            legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj-1]) + " " + metric_units + ")"
    ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    # grid
    ax.grid(linestyle="--", linewidth=1, which="major", axis="y")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_hovmoeller(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
                  metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None,
                  diagnostic_units=None, regions=None, shading=False, plot_ref=False):
    # get data
    variables = dict_param["varpattern"]
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii + 1), regions[vv].replace("nino", "N"))
    if isinstance(variables, str) is True or isinstance(variables, str) is True:
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values, models2=models2,
                 member=member,
                 shading=shading)
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        tim = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
        lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
    elif isinstance(filename_nc, dict) is True and shading is True:
        tim = list(NUMPYarray(tab_mod[0][0][0].coords[tab_mod[0][0][0].dims[0]]))
        lon = list(NUMPYarray(tab_mod[0][0][0].coords[tab_mod[0][0][0].dims[1]]))
    else:
        tim = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[0]]))
        lon = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[1]]))
    nbr_years = len(tim) / 12.
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if isinstance(filename_nc, list) or (isinstance(filename_nc, dict) is True and shading is True):
        nbr_panel = nbr_panel + ((len(tab_mod) - 1) * nbr_val)
    if plot_ref is True:
        nbr_panel = nbr_panel // 2
    title = dict_param["title"]
    if isinstance(filename_nc, list):
        if plot_ref is True:
            title = [obsname]
        else:
            if isinstance(member, list) is True and len(member) == len(model):
                title = [obsname] + [mod + mem for mod, mem in zip(model, member)]
            else:
                title = [obsname] + model
        if nbr_val > 1:
            title = title * nbr_val
    elif isinstance(filename_nc, dict) is True:
        if isinstance(member, list) is True and len(member) == len(model):
            title = [obsname] + [mod.upper() + mem + " (" + str(len(models2[mod])) + ")"
                                 for mod, mem in zip(model, member)]
        else:
            title = [obsname] + [mod.upper() + "(" + str(len(models2[mod])) + ")" for mod in model]
        if nbr_val > 1:
            title = title * nbr_val
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    elif isinstance(filename_nc, dict) is True and shading is True:
        units = tab_mod[0][0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    if units != "":
        zname = zname + " (" + units + ")"
    colorbar = "cmo." + dict_param["colorbar"]
    labelbar = dict_param["label"]
    if shading is True and len(model) + 1 == 3:
        nbrl = nbr_panel // 3
        nbrc = 1 if nbr_panel == 1 else 3
    else:
        nbrl = nbr_panel // 2
        nbrc = 1 if nbr_panel == 1 else 2
    if plot_ref is True:
        nbrl = deepcopy(nbr_panel)
        nbrc = 1
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbr_years * nbrl), sharex="col", sharey="row")
    hspa1 = 0.3 / nbr_years
    hspa2 = 0.01 / nbr_years
    plt.subplots_adjust(hspace=hspa1, wspace=0.1)
    xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    xlabel_ticks, xlabel = create_labels(xname, xlabel_ticks)
    ylabel_ticks = list(range(int(MATHfloor(min(tim))), int(MATHceil(max(tim))) + 1))
    ylabel_ticks, ylabel = create_labels(yname, ylabel_ticks)
    tab = list()
    legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, shading=shading)
    if plot_for_wiki is True:
        legend[1] = "model"
    if plot_ref is True:
        legend = [legend[0]]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            if plot_ref is False:
                tab.append(tab_mod[ii])
    elif isinstance(filename_nc, dict) is True and shading is True:
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            for kk in range(len(tab_mod)):
                tmp = [tab_mod[kk][jj][ii] for jj in range(len(tab_mod[kk]))]
                tab.append(my_average(tmp, axis=0))
    else:
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            if plot_ref is False:
                for kk in range(len(tab_mod)):
                    tab.append(tab_mod[kk][ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif (nbrl == 1 and nbrc != 1) or (nbrl != 1 and nbrc == 1):
            if nbrl == 1 and nbrc != 1:
                ax = axes[ii % nbrc]
            else:
                ax = axes[ii % nbrl]
        else:
            ax = axes[ii // nbrc, ii % nbrc]
        # title
        if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
            tt = title[ii] if plot_ref is True else title[ii // 2]
            ax.set_title(tt, fontsize=15, y=1. + hspa2, loc="left")
            if ii in [0, 1] and len(legend) - 1 >= ii % 2:
                ax.text(0.5, 1. + hspa1, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes)
        else:
            if nbrc == 2:
                if ii % nbrc == 0:
                    ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="left")
                else:
                    ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="right")
            else:
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="center")
            if nbr_val > 1:
                if nbrc == 2:
                    if ii % (len(filename_nc) + 1) == 0:
                        ttx2 = ax.get_position().x1
                    elif (ii - 1) % (len(filename_nc) + 1) == 0:
                        ttx1 = ax.get_position().x0
                        tty2 = ax.get_position().y1
                        ax.text(ttx2 + (ttx1 - ttx2) / 2., tty2 + hspa2 / 2.,
                                dict_param["title"][(ii - 1) / (len(filename_nc) + 1)], fontsize=15, weight="bold",
                                horizontalalignment="center", verticalalignment="center", transform=fig.transFigure)
                else:
                    if ii % nbrc == 1:
                        ttx2 = ax.get_position().x1
                        ttx1 = ax.get_position().x0
                        tty2 = ax.get_position().y1
                        ax.text(ttx2 + (ttx1 - ttx2) / 2., tty2 + (hspa2 * 4),
                                dict_param["title"][(ii - 1) / (len(filename_nc) + 1)], fontsize=15, weight="bold",
                                horizontalalignment="center", verticalalignment="center", transform=fig.transFigure)
        # x axis
        ax.set_xlim(xmin=min(lon), xmax=max(lon))
        if ii >= nbr_panel - nbrc:
            ax.set_xticks(xlabel_ticks)
            ax.set_xticklabels(xlabel)
            ax.set_xlabel(xname, fontsize=15)
        for axis_tick_label in ax.get_xticklabels():
            axis_tick_label.set_fontsize(12)
        # y axis
        ax.set_ylim(ymin=min(tim), ymax=max(tim))
        if ii % nbrc == 0:
            ax.set_yticks(ylabel_ticks)
            ax.set_yticklabels(ylabel)
            ax.set_ylabel(yname, fontsize=15)
        for axis_tick_label in ax.get_yticklabels():
            axis_tick_label.set_fontsize(12)
        # hovmoeller
        levels = create_levels(labelbar)
        xx, yy = NUMPYmeshgrid(lon, tim)
        cs = ax.contourf(xx, yy, tab[ii], levels=levels, extend="both", cmap=colorbar)
        if ii == 0 and plot_ref is True:
            tx1, tx2 = ax.get_xlim()
            dx = (tx2 - tx1) / 100.
            ty1, ty2 = ax.get_ylim()
            ax.text(tx2 + 2 * dx, ty2, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
        if ii == nbr_panel - nbrc:
            x1 = ax.get_position().x0
        if ii == nbr_panel - 1:
            x2 = ax.get_position().x1
            y1 = ax.get_position().y0
    # add colorbar
    if nbr_years == 1 and nbrl == 1:
        cax = plt.axes([x1, -0.1, x2 - x1, 0.04])
    else:
        if nbr_panel in [5, 6]:
            if nbr_years == 6:
                cax = plt.axes([x1, y1, x2 - x1, 0.005])
            elif nbr_years == 1:
                cax = plt.axes([x1, y1 - 0.08, x2 - x1, 0.02])
            else:
                cax = plt.axes([x1, y1 - 0.12, x2 - x1, 0.02])
        elif nbrl == 2 and nbrc <= 2 and nbr_years == 6:
            cax = plt.axes([x1, 0.09, x2 - x1, 0.006])
        elif nbr_panel in [7, 8] and nbr_years == 6:
            cax = plt.axes([x1, 0.1, x2 - x1, 0.002])
        elif nbr_panel in [11, 12] and nbr_years == 6:
            cax = plt.axes([x1, 0.1, x2 - x1, 0.002])
        elif nbr_panel in [17, 18]:
            if nbr_years < 1:
                cax = plt.axes([x1, 0.07, x2 - x1, 0.01])
            else:
                cax = plt.axes([x1, 0.08, x2 - x1, 0.005])
        else:
            cax = plt.axes([x1, 0.05, x2-x1, 0.015])
    cbar = plt.colorbar(cs, cax=cax, orientation="horizontal", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(zname, fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_map(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
           metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None,
           regions=None, shading=False, plot_ref=False):
    # get data
    variables = dict_param["varpattern"]
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii + 1), regions[vv].replace("nino", "N"))
    if isinstance(variables, str) is True or isinstance(variables, str) is True:
        nbr_val = 1
    else:
        nbr_val = len(variables)
    if "africaSE" in variables or "america" in variables or "asiaS" in variables or "oceania" in variables:
        met_in_file = True
        my_reg = "africaSE" if "africaSE" in variables else (
            "americaN" if "americaN" in variables else (
                "americaS" if "americaS" in variables else ("asiaS" if "asiaS" in variables else "oceania")))
    elif isinstance(variables, list) is True and (
            "africaSE" in variables[0] or "america" in variables[0] or "asiaS" in variables[0] or
            "oceania" in variables[0]):
        met_in_file = True
        my_reg = "africaSE" if "africaSE" in variables[0] else (
            "americaN" if "americaN" in variables[0] else (
                "americaS" if "americaS" in variables[0] else ("asiaS" if "asiaS" in variables[0] else "oceania")))
    elif isinstance(variables, list) is True and "_nina_" in variables[0] and "_nino_" in variables[1]:
        met_in_file = True
        my_reg = ""
    else:
        met_in_file = False
        my_reg = None
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values, models2=models2,
                 member=member, shading=shading, met_in_file=met_in_file, met_type=metric_type, met_pattern=my_reg)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
        lat = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
    elif isinstance(filename_nc, dict) is True and shading is True:
        lon = list(NUMPYarray(tab_mod[0][0][0].coords[tab_mod[0][0][0].dims[1]]))
        lat = list(NUMPYarray(tab_mod[0][0][0].coords[tab_mod[0][0][0].dims[0]]))
    else:
        lon = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[1]]))
        lat = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[0]]))
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True:
        nbr_panel = nbr_panel + ((len(tab_mod) - 1) * nbr_val)
    if plot_ref is True:
        nbr_panel = nbr_panel // 2
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        nbr_mod = 2
    else:
        nbr_mod = len(model) + 1
    if isinstance(filename_nc, list) is True:
        title = [obsname]
        if plot_ref is False:
            if isinstance(member, list) is True and len(member) == len(model):
                title += [mod + mem for mod, mem in zip(model, member)]
            else:
                if plot_for_wiki is True:
                    title += ["model"]
                else:
                    title += model
        if nbr_val > 1:
            title = title * nbr_val
    elif isinstance(filename_nc, dict) is True:
        # title = [obsname] + [mod.upper() + " (" + str(len(models2[mod])) + ")" for mod in model]
        tmp_let = ["b) ", "c) ", "d) "]
        title = ["a) ref: " + obsname]
        if isinstance(member, list) is True and len(member) == len(model):
            title += [tmp_let[ii] + mod.upper() + mem + " (" + str(len(models2[mod])) + ")"
                      for ii, (mod, mem) in enumerate(zip(model, member))]
        else:
            title += [tmp_let[ii] + mod.upper() + "(" + str(len(models2[mod])) + ")" for ii, mod in enumerate(model)]
        if nbr_val > 1:
            title = title * nbr_val
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    elif isinstance(filename_nc, dict) is True and shading is True:
        units = tab_mod[0][0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
    if units != "":
        zname = zname + " (" + units + ")"
    colorbar = "cmo." + dict_param["colorbar"]
    labelbar = dict_param["label"]
    if "maskland" in list(dict_param.keys()):
        maskland = dict_param["maskland"]
    else:
        maskland = False
    if "maskocean" in list(dict_param.keys()):
        maskocean = dict_param["maskocean"]
    else:
        maskocean = False
    if shading is True and len(model) + 1 == 3:
        nbrl = nbr_panel // 3
        nbrc = 1 if nbr_panel == 1 else 3
    else:
        nbrl = nbr_panel // 2
        nbrc = 1 if nbr_panel == 1 else 2
    if plot_ref is True:
        nbrl = deepcopy(nbr_panel)
        nbrc = 1
    if article_fig is True:
        if "EnsoPrMap" not in figure_name:
            nbrc = 1
            nbrl = 3
    if (isinstance(variables, str) is True and (
            "reg_pr_over_sst_map" in variables or "reg_slp_over_sst_map" in variables or
            "reg_ts_over_sst_map" in variables or "djf_map__" in variables or "jja_map__" in variables)) or\
        (isinstance(variables, list) is True and ("djf_map__" in variables[0] or "jja_map__" in variables[0])):
        fig, axes = plt.subplots(nbrl, nbrc, figsize=(6 * nbrc, 6 * nbrl), sharex="col", sharey="row", 
                                 subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
    else:
        fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row", 
                                 subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
    hspa1 = 0.1
    hspa2 = 0.01
    if ((nbrc == 2 and nbrl == 2) or (nbrc == 1 and plot_ref is True)) and isinstance(variables, list) is True and\
            my_reg in ["africaSE", "americaN", "americaS", "asiaS", "oceania"]:
        if my_reg == "africaSE":
            hspace = 0.3
        elif my_reg == "americaN":
            hspace = 0.1
        elif my_reg == "americaS":
            hspace = 0.4
        elif my_reg == "asiaS":
            hspace = 0.1
        else:
            hspace = 0.0
        plt.subplots_adjust(hspace=hspace, wspace=0.2)
    elif ((nbrc == 2 and nbrl == 2 and plot_metric is True) or (nbrc == 1 and plot_ref is True)) and\
            isinstance(variables, list) is True and my_reg == "":
        plt.subplots_adjust(hspace=-0.58, wspace=0.2)
    elif nbr_panel / float(nbrc) <= 2:
        if nbrc == 3 and nbr_val > 1:
            plt.subplots_adjust(hspace=-0.70, wspace=0.2)
        else:
            plt.subplots_adjust(hspace=-0.75, wspace=0.2)
    elif nbr_panel / float(nbrc) <= 3:
        plt.subplots_adjust(hspace=-0.8, wspace=0.2)
    elif nbr_panel / float(nbrc) <= 4:
        plt.subplots_adjust(hspace=-0.85, wspace=0.2)
    elif nbr_panel / float(nbrc) <= 6:
        plt.subplots_adjust(hspace=-0.9, wspace=0.2)
    else:
        plt.subplots_adjust(hspace=hspa1, wspace=0.2)
    if article_fig is True and nbrl == 3:
        plt.subplots_adjust(hspace=-0.7, wspace=0.2)
    if "BiasPrLatRmse" in figure_name and article_fig is True and nbrc == 1 and nbrl == 3:
        plt.subplots_adjust(hspace=-0.80, wspace=0.2)
    xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    xlabel_ticks, xlabel = create_labels(xname, xlabel_ticks)
    ylabel_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
    ylabel_ticks, ylabel = create_labels(yname, ylabel_ticks)
    tab = list()
    legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                       shading=shading)
    if plot_for_wiki is True:
        if "CNRM-CM5" in legend[0]:
            legend[0] = "model"
        elif "CNRM-CM5" in legend[1]:
            legend[1] = "model"
    if plot_ref is True:
        legend = [legend[0]]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        for ii in range(len(tab_mod)):
            tab.append(tab_obs[ii])
            if plot_ref is False:
                tab.append(tab_mod[ii])
    elif isinstance(filename_nc, dict) is True and shading is True:
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            for kk in range(len(tab_mod)):
                if "africaSE" in variables or "america" in variables or "asiaS" in variables or "oceania" in variables:
                    tmp = list()
                    for jj in range(len(tab_mod[kk])):
                        tmp.append(my_mask_map(tab_obs[ii], tab_mod[kk][jj][ii]))
                else:
                    tmp = [tab_mod[kk][jj][ii] for jj in range(len(tab_mod[kk]))]
                tab.append(my_average(tmp, axis=0))
    else:
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            if plot_ref is False:
                for kk in range(len(tab_mod)):
                    tab.append(tab_mod[kk][ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif (nbrl == 1 and nbrc != 1) or (nbrl != 1 and nbrc == 1):
            if nbrl == 1 and nbrc != 1:
                ax = axes[ii % nbrc]
            else:
                ax = axes[ii % nbrl]
        else:
            ax = axes[ii // nbrc, ii % nbrc]
        if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
            tt = title[ii] if plot_ref is True else title[ii // 2]
            ax.set_title(tt, fontsize=15, y=1., loc="left")
            if ii in [0, 1] and len(legend) -1 >= ii % 2:
                ax.text(0.5, 1. + (0.15 * (len(lon) + 10)) / len(lat), legend[ii % 2], fontsize=15, weight="bold",
                        horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
        else:
            if nbrc == 2 or plot_metric is True:
                if ii == 0 or (ii % nbr_mod == 0):
                    location = "left"
                else:
                    location = "left"#"right"
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc=location)
                del location
            else:
                if "BiasPrLatRmse" in figure_name and article_fig is True and nbrc == 1 and nbrl == 3:
                    ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="left")
                else:
                    ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="center")
            if nbr_val > 1:
                if nbrc == 2:
                    if ii % nbrc == 0:
                        ax.text(1.1, 1. + 20 * hspa2, dict_param["title"][ii // nbrc],
                                fontsize=15, weight="bold", horizontalalignment="center", verticalalignment="center",
                                transform=ax.transAxes)
                    # if ii % (len(filename_nc) + 1) == 0:
                    #     ttx2 = ax.get_position().x1
                    # elif (ii - 1) % (len(filename_nc) + 1) == 0:
                    #     ttx1 = ax.get_position().x0
                    #     tty2 = ax.get_position().y1
                    #     ax.text(ttx2 + (ttx1 - ttx2) / 2., tty2 + hspa2 / 2.,
                    #             dict_param["title"][(ii - 1) / (len(filename_nc) + 1)], fontsize=15, weight="bold",
                    #             horizontalalignment="center", verticalalignment="center", transform=fig.transFigure)
                else:
                    if ii % nbrc == 1:
                        ax.text(0.5, 1. + 60 * hspa2, dict_param["title"][(ii - 1) / (len(filename_nc) + 1)],
                                fontsize=15, weight="bold", horizontalalignment="center", verticalalignment="center",
                                transform=ax.transAxes)
        # map
        xx, yy = NUMPYmeshgrid(lon, lat)
        # set extent
        if lat[-1] - lat[0] < 40:
            ax.set_extent([lon[0], lon[-1], lat[0] - 5, lat[-1] + 5], crs=ccrs.PlateCarree())
        else:
            ax.set_extent([lon[0], lon[-1], lat[0], lat[-1]], crs=ccrs.PlateCarree())  
        # draw coastlines
        ax.coastlines()
        # fill continents
        if maskland:
            ax.add_feature(cfeature.LAND, color="gainsboro")
        if maskocean:
            ax.add_feature(cfeature.OCEAN, color="white")
        # adjust the yticks to convert longitude over 180 to negative to properly add gridlines
        xlabel_ticks_adjusted = [i if i < 180 else i - 360 for i in xlabel_ticks]
        # draw parallels and meridians by adding grid lines only at specified ticks
        gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', color='k')
        gl.xlocator = mticker.FixedLocator(xlabel_ticks_adjusted)
        gl.ylocator = mticker.FixedLocator(ylabel_ticks)
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()        
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}
        # contour plot
        levels = create_levels(labelbar)
        cs = ax.contourf(xx, yy, tab[ii], levels=levels, extend="both", cmap=colorbar, transform=ccrs.PlateCarree())  
        # my text
        if (ii > 0 and plot_metric is True and isinstance(variables, list) is False) or\
                (isinstance(variables, list) is True and "nina" in variables[0] and "nino" in variables[1] and
                 (ii+1) % 2 == 0):
            if isinstance(metric_type, str) is True or isinstance(metric_type, str) is True:
                txt = format_metric(metric_type, my_average(metval[ii - 1], remove_masked=True), metric_units)
                ax.text(0.5, 0.5, txt, fontsize=12, color="k", horizontalalignment="center",
                        verticalalignment="center")
            elif isinstance(metric_type, list) is True:
                for jj in range(len(metric_type)):
                    if shading is True:
                        tmp = [metval[ii - 1][kk][jj] for kk in range(len(metval[ii - 1]))]
                        tmp = my_average(tmp, remove_masked=True)
                    else:
                        if isinstance(metval[0], list) is True:
                            tmp = metval[ii // 2][jj]
                        else:
                            tmp = metval[jj]
                    txt = format_metric(metric_type[jj], tmp, metric_units[jj])
                    if my_reg in ["africaSE", "americaN", "americaS", "asiaS", "oceania"]:
                        if my_reg in ["africaSE"]:
                            xxx, yyy = 0.00, -0.05 - jj * 0.07
                        elif my_reg in ["americaN"]:
                            xxx, yyy = 0.00, -0.14 - jj * 0.10
                        elif my_reg in ["americaS"]:
                            xxx, yyy = 0.00, -0.12 - jj * 0.06
                        elif my_reg in ["asiaS"]:
                            xxx, yyy = 0.00, -0.13 - jj * 0.09
                        else:
                            xxx, yyy = 0.00, -0.15 - jj * 0.10
                    elif "reg_pr_over_sst_map" in variables or "reg_slp_over_sst_map" in variables or\
                            "reg_ts_over_sst_map" in variables or "reg_pr_over_sst_djf_map" in variables or\
                            "reg_slp_over_sst_djf_map" in variables or "reg_ts_over_sst_djf_map" in variables or\
                            "reg_pr_over_sst_jja_map" in variables or "reg_slp_over_sst_jja_map" in variables or\
                            "reg_ts_over_sst_jja_map" in variables or (isinstance(variables, list) is True and (
                            "djf_map__" in variables[0] or "jja_map__" in variables[0])):
                        xxx, yyy = 0.00, -0.30 - jj * 0.18
                    else:
                        xxx, yyy = -0.12, 1.26 - jj * 0.16
                    ax.text(xxx, yyy, txt, fontsize=11, color="k", horizontalalignment="left",
                            verticalalignment="center", transform=ax.transAxes)
        if ii == 0 and plot_ref is True:
            tx1, tx2 = ax.get_xlim()
            dx = (tx2 - tx1) / 100.
            ty1, ty2 = ax.get_ylim()
            ax.text(tx2 + 2 * dx, ty2, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
        if plot_ref is True:
            if regions is not None:
                if (isinstance(variables, list) is True  and
                        ("djf_map__" in variables[0] or "jja_map__" in variables[0])) or (
                        "djf_map__" in variables or "jja_map__" in variables):
                    lreg = ReferenceRegions("nino3.4")
                else:
                    lreg = ReferenceRegions(regions[metric_variables[0]])
                lons = [lreg["longitude"][0]] * 2 + [lreg["longitude"][1]] * 2
                lats = list(lreg["latitude"]) + list(reversed(list(lreg["latitude"])))
                x, y = lons, lats
                ax.add_patch(Polygon(list(zip(x, y)), edgecolor="k", linewidth=3, linestyle="-", facecolor="none"))
        # if ii == 0 and plot_metric is True:
        #     x1 = ax.get_position().x1
        #     y2 = ax.get_position().y1
        #     if nbrl == 1:
        #         ax2 = axes[(ii + 1) % 2]
        #     else:
        #         ax2 = axes[int(round((ii + 1) / 2)), (ii + 1) % 2]
        #     x2 = ax2.get_position().x0
        #     txt = format_metric(metric_type, metval, metric_units)
        #     ax.text(x1 + ((x2 - x1) / 2.), y2 + 0.1, txt, fontsize=12, color="k", horizontalalignment="center",
        #             verticalalignment="center", transform=fig.transFigure)
        if ii == 0:
            x1 = ax.get_position().x0
        if ii == nbr_panel - 1:
            x2 = ax.get_position().x1
            y1 = ax.get_position().y0
    # add colorbar
    if my_reg in ["africaSE", "americaN", "americaS", "asiaS", "oceania"]:
        if my_reg in ["africaSE"]:
            if isinstance(variables, list) is True:
                cax = plt.axes([x1, y1 - 0.08, x2 - x1, 0.02])
            else:
                cax = plt.axes([x1, y1 - 0.15, x2 - x1, 0.04])
        elif my_reg in ["americaN"]:
            if isinstance(variables, list) is True:
                cax = plt.axes([x1, y1 - 0.10, x2 - x1, 0.02])
            else:
                cax = plt.axes([x1, y1 - 0.22, x2 - x1, 0.05])
        elif my_reg in ["americaS"]:
            if isinstance(variables, list) is True:
                cax = plt.axes([x1, y1 - 0.10, x2 - x1, 0.025])
            else:
                cax = plt.axes([x1, y1 - 0.2, x2 - x1, 0.035])
        elif my_reg in ["asiaS"]:
            if isinstance(variables, list) is True:
                cax = plt.axes([x1, y1 - 0.10, x2 - x1, 0.02])
            else:
                cax = plt.axes([x1, y1 - 0.2, x2 - x1, 0.05])
        else:
            if isinstance(variables, list) is True:
                cax = plt.axes([x1, y1 - 0.10, x2 - x1, 0.02])
            else:
                cax = plt.axes([x1, y1 - 0.2, x2 - x1, 0.05])
    elif nbrl == 2:
        if isinstance(variables, list) is True and ("djf_map__" in variables[0] or "jja_map__" in variables[0]):
            cax = plt.axes([x1, y1 - 0.09, x2 - x1, 0.018])
        else:
            cax = plt.axes([x1, y1 - 0.05, x2 - x1, 0.015])
    elif nbrl == 3:
        if article_fig is True and nbrl == 3 and (
                "reg_pr_over_sst_map" in variables or "reg_pr_over_sst_djf_map" in variables or
                "reg_pr_over_sst_jja_map" in variables):
            cax = plt.axes([x1, y1 + 0.15, x2 - x1, 0.015])
        else:
            cax = plt.axes([x1, y1 + 0.2, x2 - x1, 0.015])
    elif nbrl == 4:
        cax = plt.axes([x1, y1 + 0.21, x2 - x1, 0.01])
    elif nbrl == 6:
        cax = plt.axes([x1, y1 + 0.22, x2 - x1, 0.005])
    else:
        if "reg_pr_over_sst_map" in variables or "reg_slp_over_sst_map" in variables or\
                "reg_ts_over_sst_map" in variables or "reg_pr_over_sst_djf_map" in variables or\
                "reg_slp_over_sst_djf_map" in variables or "reg_ts_over_sst_djf_map" in variables or\
                "reg_pr_over_sst_jja_map" in variables or "reg_slp_over_sst_jja_map" in variables or\
                "reg_ts_over_sst_jja_map" in variables:
            cax = plt.axes([x1, y1 - 0.18, x2 - x1, 0.035])
        else:
            cax = plt.axes([x1, y1 - 0.1, x2 - x1, 0.03])
    cbar = plt.colorbar(cs, cax=cax, orientation="horizontal", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(zname, fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_scatterplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, models2=None, member=None,
                   metric_type=None, metric_values=None, metric_units=None, diagnostic_values=None,
                   diagnostic_units=None, regions=None, shading=False, plot_ref=False):
    # get data
    variables = dict_param["varpattern"]
    method = dict_param["method"]
    if isinstance(metric_variables, list) is True and regions is not None:
        for ii, vv in enumerate(metric_variables):
            method = method.replace("REGION" + str(ii + 1), regions[vv].replace("nino", "N"))
    if isinstance(variables, str) is True or isinstance(variables, str) is True:
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values, models2=models2,
                 member=member,
                 shading=shading)
    if metric_type is not None and (isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True):
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if (isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True) and nbr_panel > 1:
        nbr_panel = nbr_panel + len(filename_nc) - 1
    if plot_ref is True:
        nbr_panel = nbr_panel // 2
    title = dict_param["title"]
    if (isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True) and nbr_panel > 1:
        if isinstance(filename_nc, dict):
            if "EnsoFbSstSwr" in figure_name and article_fig is True:
                tmp_let = ["b) ", "c) ", "d) ", "e) "]
                title = [tmp_let[0] + "ref: Tropflux"]
                if isinstance(member, list) is True and len(member) == len(model):
                    title += [tmp_let[ii+1] + mod.upper() + mem + " (" + str(len(models2[mod])) + ")"
                              for ii, (mod, mem) in enumerate(zip(model, member))]
                else:
                    title += [tmp_let[ii+1] + mod.upper() + "(" + str(len(models2[mod])) + ")"
                              for ii, mod in enumerate(model)]
            else:
                title = [obsname]
                if plot_ref is False:
                    if isinstance(member, list) is True and len(member) == len(model):
                        title += [mod.upper() + mem + " (" + str(len(models2[mod])) + ")"
                                  for mod, mem in zip(model, member)]
                    else:
                        title += [mod.upper() + "(" + str(len(models2[mod])) + ")" for mod in model]
        else:
            title = [obsname]
            if plot_ref is False:
                if isinstance(member, list) is True and len(member) == len(model):
                    title += [mod + mem for mod, mem in zip(model, member)]
                else:
                    title += model
    if isinstance(title, str) is True or isinstance(title, str) is True:
        title = [title] * nbr_panel
    xname = dict_param["xname"]
    if isinstance(xname, str) is True or isinstance(xname, str) is True:
        xname = [xname] * nbr_panel
        one_xaxis = True
    else:
        one_xaxis = False
    yname = dict_param["yname"]
    if isinstance(yname, str) is True or isinstance(yname, str) is True:
        yname = [yname] * nbr_panel
        one_yaxis = True
    else:
        one_yaxis = False
    if "colors" in list(dict_param.keys()):
        mcolors = dict_param["colors"]
    else:
        mcolors = ["k", "dodgerblue"]
        # if dict_param["nbr_panel"] == len(nbr_val):
        #     mcolors = list(reversed(mcolors))
        if isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True:
            for ii in range(len(filename_nc) - 1):
                mcolors.append(colors_sup[ii])
        if plot_ref is True:
            mcolors = ["k"]
    if "markers" in list(dict_param.keys()):
        markers = dict_param["markers"]
    else:
        markers = ["D", "."]
        # if dict_param["nbr_panel"] == len(nbr_val):
        #     markers = list(reversed(markers))
        if isinstance(filename_nc, list) is True or isinstance(filename_nc, dict) is True:
            for ii in range(len(filename_nc) - 1):
                markers.append(".")
        if plot_ref is True:
            markers = ["D"]
    if "legend" in list(dict_param.keys()):
        legend = dict_param["legend"]
    else:
        legend = my_legend(model, obsname, filename_nc, models2=models2, member=member, plot_metric=plot_metric,
                           shading=shading)
        if plot_for_wiki is True:
            legend[1] = "model"
        if plot_ref is True:
            legend = [legend[0]]
    keys1 = ["", "_neg", "_pos"]
    keys2 = ["all", "x<0", "x>0"]
    keys3 = [[None, None], [None, 0], [0, None]]
    keys4 = ["black", "dodgerblue", "red"]
    lines = [Line2D([0], [0], marker=markers[kk], c="w", markerfacecolor=mcolors[kk], markersize=12)
             for kk in range(len(mcolors))]
    if shading is True and nbr_panel == 3:
        nbrl = nbr_panel // 3
        nbrc = 1 if nbr_panel == 1 else 3
    else:
        nbrl = 1 if nbr_panel == 1 else nbr_panel // 2
        nbrc = 1 if nbr_panel == 1 else 2
    if plot_ref is True:
        nbrl = deepcopy(nbr_panel)
        nbrc = 1
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    plt.subplots_adjust(hspace=0.3, wspace=0.1)
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        tab1 = tab_obs[::2]
        tab2 = tab_obs[1::2]
        if plot_ref is False:
            tab1 += tab_mod[::2]
            tab2 += tab_mod[1::2]
    elif isinstance(filename_nc, dict) is True and shading is True:
        tab1 = tab_obs[::2]
        tab2 = tab_obs[1::2]
        for kk in range(len(tab_mod)):
            tmp1, tmp2 = list(), list()
            for jj in range(len(tab_mod[kk])):
                tmp1 += list(tab_mod[kk][jj][::2][0])
                tmp2 += list(tab_mod[kk][jj][1::2][0])
            tab1.append(tmp1)
            tab2.append(tmp2)
            del tmp1, tmp2
    else:
        tab1 = tab_obs[::2]
        tab2 = tab_obs[1::2]
        if plot_ref is False:
            for ii in range(len(filename_nc)):
                tab1 += tab_mod[ii][::2]
                tab2 += tab_mod[ii][1::2]
    xtick_labels = minmax_plot(tab1, metric=False)
    ytick_labels = minmax_plot(tab2, metric=plot_metric)
    if nbr_panel == nbr_val / 2.:
        markers = list(reversed(markers))
        mcolors = list(reversed(mcolors))
        tab1 = list(reversed(tab1))
        tab2 = list(reversed(tab2))
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif (nbrl == 1 and nbrc != 1) or (nbrc == 1 and nbrl != 1):
            if nbrl == 1 and nbrc != 1:
                ax = axes[ii % nbrc]
            else:
                ax = axes[ii % nbrl]
        else:
            ax = axes[ii // nbrc, ii % nbrc]
        # title
        if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
            ax.set_title(title[ii // 2], fontsize=15, y=1.01, loc="left")
            if nbr_panel == len(tab_mod) and ii in [0, 1]:
                ax.text(0.5, 1.15, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes)
        else:
            ax.set_title(title[ii], fontsize=15, y=1.01, loc="left")
        # x axis
        if "EnsoFbSshSst" in figure_name and article_fig is True:
            ax.set_xticks([-40, -20, 0, 20, 40], minor=False)
            ax.set_xticklabels([-40, -20, 0, 20, 40])
            ax.set_xticks([-30, -10, 10, 30], minor=True)
            ax.set_xlim([-40, 40])
        elif ("EnsoFbSstTaux" in figure_name or "EnsoFbSstSwr" in figure_name) and article_fig is True:
            ax.set_xticks([-6, -3, 0, 3, 6], minor=False)
            ax.set_xticklabels([-6, -3, 0, 3, 6])
            ax.set_xticks([-4.5, -1.5, 1.5, 4.5], minor=True)
            ax.set_xlim([-6, 6])
        elif "EnsoFbTauxSsh" in figure_name and article_fig is True:
            ax.set_xticks([-100, -50, 0, 50, 100], minor=False)
            ax.set_xticklabels([-100, -50, 0, 50, 100])
            ax.set_xticks([-75, -25, 25, 75], minor=True)
            ax.set_xlim([-125, 125])
        else:
            ax.set_xticks(xtick_labels)
            ax.set_xlim(xmin=min(xtick_labels), xmax=max(xtick_labels))
        if (one_xaxis is True and (ii >= (nbrc * nbrl) - nbrc)) or one_xaxis is False:
            xlabel = xname[ii]
            for kk in list(regions.keys()):
                if kk in xlabel.lower():
                    xlabel = regions[kk] + " " + xlabel
            units = tab_obs[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
            if units != "":
                xlabel = xlabel + " (" + units + ")"
            ax.set_xlabel(xlabel, fontsize=15)
        for axis_tick_label in ax.get_xticklabels():
            axis_tick_label.set_fontsize(12)
        # y axis
        if "EnsoFbSshSst" in figure_name and article_fig is True:
            ax.set_yticks([-6, -3, 0, 3, 6], minor=False)
            ax.set_yticklabels([-6, -3, 0, 3, 6])
            ax.set_yticks([-4.5, -1.5, 1.5, 4.5], minor=True)
            ax.set_ylim([-6, 6])
        elif "EnsoFbSstTaux" in figure_name and article_fig is True:
            ax.set_yticks([-100, -50, 0, 50, 100], minor=False)
            ax.set_yticklabels([-100, -50, 0, 50, 100])
            ax.set_yticks([-75, -25, 25, 75], minor=True)
            ax.set_ylim([-125, 125])
        elif "EnsoFbTauxSsh" in figure_name and article_fig is True:
            ax.set_yticks([-30, -15, 0, 15, 30], minor=False)
            ax.set_yticklabels([-30, -15, 0, 15, 30])
            ax.set_yticks([-22.5, -7.5, 0, 7.5, 22.5], minor=True)
            ax.set_ylim([-37.5, 37.5])
        else:
            ax.set_yticks(ytick_labels)
            ax.set_ylim(ymin=min(ytick_labels), ymax=max(ytick_labels))
        # ax.set_yticks([-150, -75, 0, 75, 150], minor=False)
        # ax.set_yticklabels([-150, -75, 0, 75, 150])
        # ax.set_yticks([-112.5, -37.5, 37.5, 112.5], minor=True)
        # ax.set_ylim([-150, 150])
        # ax.set_yticks([-100, -50, 0, 50, 100], minor=False)
        # ax.set_yticklabels([-100, -50, 0, 50, 100])
        # ax.set_yticks([-75, -25, 25, 75], minor=True)
        # ax.set_ylim([-100, 100])
        if (one_yaxis is True and ii % nbrc == 0) or one_yaxis is False:
            ylabel = yname[ii]
            for kk in list(regions.keys()):
                if kk in ylabel.lower():
                    ylabel = regions[kk] + " " + ylabel
            units = tab_obs[1].units.replace("C", "$^\circ$C").replace("long", "$^\circ$lon")
            if units != "":
                ylabel = ylabel + " (" + units + ")"
            ax.set_ylabel(ylabel, fontsize=15)
        for axis_tick_label in ax.get_yticklabels():
            axis_tick_label.set_fontsize(12)
        # scatterplots and slopes
        x1, x2 = ax.get_xlim()
        dx = (x2 - x1) / 100.
        y1, y2 = ax.get_ylim()
        dy = (y2 - y1) / 100.
        if (nbr_panel > 1 and nbr_panel == nbr_val / 2.) or (nbr_panel == len(legend) and title[0] == "nonlinarity"):
            # multiple panel
            ax.scatter(tab1[ii], tab2[ii], s=10, c="k", marker=markers[ii])
            for jj in range(len(keys1)):
                if isinstance(filename_nc, dict):
                    intercept, slope = list(), list()
                    if ii == 0:
                        if "slope" + keys1[jj] in list(tab_obs[0].attrs.keys()) and\
                                "intercept" + keys1[jj] in list(tab_obs[0].attrs.keys()):
                            intercept.append(tab_obs[0].attrs["intercept" + keys1[jj]])
                            slope.append(tab_obs[0].attrs["slope" + keys1[jj]])
                    else:
                        for kk in range(len(tab_mod[ii-1])):
                            if "slope" + keys1[jj] in list(tab_mod[ii-1][kk][0].attrs.keys()) and\
                                    "intercept" + keys1[jj] in list(tab_mod[ii-1][kk][0].attrs.keys()):
                                intercept.append(tab_mod[ii-1][kk][0].attrs["intercept" + keys1[jj]])
                                slope.append(tab_mod[ii-1][kk][0].attrs["slope" + keys1[jj]])
                    if len(intercept) == 1:
                        intercept = intercept[0]
                        slope = slope[0]
                    else:
                        intercept = float(my_average(intercept, remove_masked=True))
                        slope = float(my_average(slope, remove_masked=True))
                    col = keys4[jj]
                    xx = keys3[jj]
                    if xx[0] is None:
                        xx[0] = x1
                    if xx[1] is None:
                        xx[1] = x2
                    yy = [kk * slope + intercept for kk in xx]
                    ax.plot(xx, yy, lw=2, c=col)
                    txt = "slope(" + keys2[jj] + ") = " + "{0:+.2f}".format(round(slope, 2))
                    # ax.text(dx + x1, ((97 - 5 * jj) * dy) + y1, txt, fontsize=12, color=col,
                    #         horizontalalignment="left", verticalalignment="center")
                    ax.text(2 * dx + x1, ((93 - 6 * jj) * dy) + y1, txt, fontsize=12, color=col,
                            horizontalalignment="left", verticalalignment="center")
                else:
                    if "slope" + keys1[jj] in list(tab1[ii].attrs.keys()) and\
                            "intercept" + keys1[jj] in list(tab1[ii].attrs.keys()):
                        col = keys4[jj]
                        slope = tab1[ii].attrs["slope" + keys1[jj]]
                        intercept = tab1[ii].attrs["intercept" + keys1[jj]]
                        xx = keys3[jj]
                        if xx[0] is None:
                            xx[0] = x1
                        if xx[1] is None:
                            xx[1] = x2
                        yy = [kk * slope + intercept for kk in xx]
                        ax.plot(xx, yy, lw=2, c=col)
                        txt = "slope(" + keys2[jj] + ") = " + "{0:+.2f}".format(round(slope, 2))
                        # ax.text(dx + x1, ((93 - 6 * jj) * dy) + y1, txt, fontsize=12, color=col,
                        #         horizontalalignment="left", verticalalignment="center")
                        ax.text(2 * dx + x1, ((93 - 6 * jj) * dy) + y1, txt, fontsize=12, color=col,
                                horizontalalignment="left", verticalalignment="center")
        else:
            tmp = list()
            for jj in range(len(tab1)):
                col = mcolors[jj]
                ax.scatter(tab1[jj], tab2[jj], s=10, c=col, marker=markers[jj])
                if isinstance(filename_nc, dict):
                    intercept, slope = list(), list()
                    if jj == len(tab1) - 1:
                        if "slope" + keys1[0] in list(tab_obs[0].attrs.keys()) and\
                                "intercept" + keys1[0] in list(tab_obs[0].attrs.keys()):
                            intercept.append(tab_obs[0].attrs["intercept" + keys1[0]])
                            slope.append(tab_obs[0].attrs["slope" + keys1[0]])
                    else:
                        for kk in range(len(tab_mod[len(tab_mod) - 1 - jj])):
                            if "slope" + keys1[0] in list(tab_mod[len(tab_mod) - 1 - jj][kk][0].attrs.keys()) and\
                                    "intercept" + keys1[0] in list(tab_mod[len(tab_mod) - 1 - jj][kk][0].attrs.keys()):
                                intercept.append(tab_mod[len(tab_mod) - 1 - jj][kk][0].attrs["intercept" + keys1[0]])
                                slope.append(tab_mod[len(tab_mod) - 1 - jj][kk][0].attrs["slope" + keys1[0]])
                    if len(intercept) == 1:
                        intercept = intercept[0]
                        slope = slope[0]
                    else:
                        intercept = float(my_average(intercept, remove_masked=True))
                        slope = float(my_average(slope, remove_masked=True))
                    xx = [x1, x2]
                    yy = [kk * slope + intercept for kk in xx]
                    if jj == 0:
                        ax.plot(xx, yy, lw=4, c=col)
                    else:
                        ax.plot(xx, yy, lw=2, c=col)
                    tmp.append("slope = " + "{0:+.2f}".format(round(slope, 2)))
                else:
                    if "slope" in list(tab1[jj].attrs.keys()) and "intercept" in list(tab1[jj].attrs.keys()):
                        slope = tab1[jj].attrs["slope"]
                        intercept = tab1[jj].attrs["intercept"]
                        xx = [x1, x2]
                        yy = [kk * slope + intercept for kk in xx]
                        ax.plot(xx, yy, lw=2, c=col)
                        tmp.append("slope = " + "{0:+.2f}".format(round(slope, 2)))
            for jj in range(len(tmp)):
                col = mcolors[len(tmp) - 1 - jj]
                ax.text(2 * dx + x1, ((93 - 6 * jj) * dy) + y1, tmp[len(tmp) - 1 - jj], fontsize=12,
                        color=col, horizontalalignment="left", verticalalignment="center")
                # ax.text(2 * dx + x1, ((13 - 6 * jj) * dy) + y1, tmp[len(tmp) - 1 - jj], fontsize=12,
                #         color=col, horizontalalignment="left", verticalalignment="center")
            if nbr_panel == 1 or ii == nbrc - 1:
                if metric_type is not None:
                    for jj in range(1, len(legend)):
                        if isinstance(model, str) is True or isinstance(model, str) is True:
                            tmp = deepcopy(metval)
                        elif shading is True:
                            tmp = my_average(metval[jj - 1], remove_masked=True)
                        else:
                            tmp = metval[jj - 1]
                        # metric_units = "N/m2/$^\circ$C / N/m2/$^\circ$C"
                        # metric_units = "cm/N/m2 / cm/N/m2"
                        # legend[jj] = legend[jj] + " (" + "{0:.2f}".format(tmp) + " " + metric_units + ")"
                        if ("EnsoFbSshSst" in figure_name or "EnsoFbSstTaux" in figure_name or\
                                "EnsoFbTauxSsh" in figure_name) and article_fig is True:
                            metric_units = "%"
                            legend[jj] = legend[jj] + " (" + "{0:}".format(int(round(tmp))) + " " + metric_units + ")"
                if ("EnsoFbSshSst" in figure_name or "EnsoFbSstTaux" in figure_name or "EnsoFbTauxSsh" in figure_name) \
                        and article_fig is True:
                    ax.legend(lines, legend, bbox_to_anchor=(1, 0), loc="lower right", ncol=1)
                else:
                    ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        # my text
        if plot_metric is True:
            txt = format_metric(metric_type, metval, metric_units)
            ax.text(x2 - (2 * dx), y1 + (6 * dy), txt, fontsize=12, color="k", horizontalalignment="right",
                    verticalalignment="center")
        if ii == 0 and plot_ref is True:
            ax.text(x2 + 2 * dx, y2 - 20 * dy, str(method, "utf-8"), fontsize=12, color="k", ha="left", va="top")
        # grid
        ax.grid(linestyle="--", linewidth=1, which="major")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def plot_curve(tab_mod, tab_obs, ax, title, axis, xname, yname, ytick_labels, linecolors, linestyles, metric_type,
               metval, metric_units, model='', member=None, obsname='', legend=[], multimodel=False, plot_metric=False,
               plot_legend=False, shading=False, plot_ref=False, method=""):
    # title
    ax.set_title(title, fontsize=15, y=1.01, loc="left")
    # x axis
    label_ticks = list(range(int(MATHfloor(min(axis))), int(MATHceil(max(axis))) + 1))
    label_ticks, label = create_labels(xname, label_ticks)
    ax.set_xticks(label_ticks)
    ax.set_xticklabels(label)
    ax.set_xlim([min(axis), max(axis)])
    # ax.set_xlim([-13, 13])
    ax.set_xlabel(xname, fontsize=15)
    for axis_tick_label in ax.get_xticklabels():
        axis_tick_label.set_fontsize(12)
    # y axis
    ax.set_yticks(ytick_labels)
    ax.set_yticklabels(ytick_labels)
    ax.set_ylim([min(ytick_labels), max(ytick_labels)])
    # ax.set_yticks([0, 3, 6, 9], minor=False)
    # ax.set_yticklabels([0, 3, 6, 9])
    # ax.set_yticks([1.5, 4.5, 7.5], minor=True)
    # ax.set_ylim([0, 9.5])
    # ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0], minor=False)
    # ax.set_yticklabels([-1.0, -0.5, 0.0, 0.5, 1.0])
    # ax.set_yticks([-0.75, -0.25, 0.25, 0.75], minor=True)
    # ax.set_ylim([-1.1, 1.1])
    # ax.add_line(Line2D([29, 39], [0.25, 0.25], c="orange", lw=2))
    # ax.scatter([29], 0.25, s=80, c="orange", marker="<", zorder=10)
    # ax.scatter([39], 0.25, s=80, c="orange", marker=">", zorder=10)
    # ax.text(34, 0.1, "duration", fontsize=18, color="orange", horizontalalignment='center',
    #          verticalalignment='center')
    ax.set_ylabel(yname, fontsize=15)
    for axis_tick_label in ax.get_yticklabels():
        axis_tick_label.set_fontsize(12)
    if min(ytick_labels) < 0 and max(ytick_labels) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # plot curves
    if len(tab_mod) + len(tab_obs) > 2:
        lw = 4  # 2  #
    else:
        lw = 4
    if plot_ref is False:
        for ii, tab in enumerate(tab_mod):
            if shading is True:
                tab_sh = shading_levels(tab, axis=0)
                # # !!!!! temporary: start !!!!!
                # if ii != len(tab_mod) - 1:
                #     ax.fill_between(axis, list(tab_sh[0]), list(tab_sh[3]), facecolor=linecolors["model"][ii],
                #                     alpha=0.3)
                #     ax.fill_between(axis, list(tab_sh[1]), list(tab_sh[2]), facecolor=linecolors["model"][ii],
                #                     alpha=0.4)
                # # !!!!! temporary: end !!!!!
                # ax.fill_between(axis, list(tab_sh[0]), list(tab_sh[3]), facecolor=linecolors["model"][ii], alpha=0.3)
                # ax.fill_between(axis, list(tab_sh[1]), list(tab_sh[2]), facecolor=linecolors["model"][ii], alpha=0.4)
                ax.plot(axis, list(tab_sh[4]), lw=lw, color=linecolors["model"][ii], ls=linestyles["model"][ii])
            else:
                ax.plot(axis, list(tab), c=linecolors["model"][ii], lw=lw, ls=linestyles["model"][ii])
    for ii, tab in enumerate(tab_obs):
        ax.plot(axis, list(tab), c=linecolors["reference"][ii], lw=lw, ls=linestyles["reference"][ii])
    # relative space
    x1, x2 = ax.get_xlim()
    dx = (x2 - x1) / 100.
    y1, y2 = ax.get_ylim()
    dy = (y2 - y1) / 100.
    if multimodel is False:
        # legend
        if len(linecolors["model"]) == 1 and len(linecolors["reference"]) == 1:
            if plot_ref is False:
                legen = deepcopy(legend)
                legco = [linecolors["reference"][0], linecolors["model"][0]]
            else:
                legen = [legend[0]]
                legco = [linecolors["reference"][0]]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in legco]
            ax.legend(lines, legen, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
        else:
            legtxt = deepcopy(legend)
            if plot_ref is False:
                if isinstance(member, str) is True or isinstance(member, str) is True:
                    legtxt += [model + " " + member, obsname]
                else:
                    legtxt += [model, obsname]
                legls = [linestyles["model"][0], linestyles["reference"][0]]
            else:
                legtxt += [obsname]
                legls = [linestyles["reference"][0]]
            if plot_for_wiki is True:
                legtxt = ["model" if "CNRM-CM5" in tt else tt for tt in legtxt]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in
                     linecolors["model"]]
            lines = lines + [Line2D([0], [0], c='k', lw=2, ls=ls) for ls in legls]
            ax.legend(lines, legtxt, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
        if plot_ref is True:
            ax.text(x2 + 2 * dx, y2 - 50 * dy, str(method, "utf-8"), fontsize=12, color="k", ha="left",
                     va="top")
        # my text
        if plot_metric is True:
            txt = format_metric(metric_type, metval, metric_units)
            ax.text(x2 - (2 * dx), y2 - (6 * dy), txt, fontsize=12, color='k', ha='right', va='center')
    else:
        # legend
        if plot_legend is True:
            legco = linecolors["reference"] + linecolors["model"]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in legco]
            if plot_metric is True:
                for jj in range(1, len(legend)):
                    if shading is True:
                        tmp = my_average(metval[jj - 1], remove_masked=True)
                    else:
                        tmp = metval[jj - 1]
                    legend[jj] = legend[jj] + " (" + "{0:.2f}".format(tmp) + " " + metric_units + ")"
            # ax.legend(lines, legend, bbox_to_anchor=(0, 1), loc="upper left", ncol=1)
            legend = ["model" if "CNRM-CM5" in tt else tt for tt in legend]
            ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(linestyle='--', linewidth=1, which='major')
    return


def plot_metrics_correlations(tab_rval, figure_name, xy_names, tab_pval=None, write_corr=False, title="", cfram=False,
                              chigh=False, bold_names=[], save_eps=False):
    """
    Plots the correlations matrix

    Inputs:
    ------
    :param tab_rval: `cdms2` variable
        A `cdms2` variable containing the correlation coefficients.
    :param figure_name: string
        Path to and name of the output plot.
    :param xy_names: list of string
        List of the names to put as x and y ticks.
    **Optional arguments:**
    :param tab_pval: `cdms2` variable, optional
        A `cdms2` variable containing the p-value of the correlation coefficients. It is used to mask the correlations
        that are not significant.
        Default is None (all correlation coefficients are plotted).
    :param write_corr: boolean, optional
        True to write the correlation value in each box.
        Default is False (correlation not written).
    :param title: string, optional
        Title of the plot.
        Default is "" (no title).
    :param chigh: boolean, optional
        True to highlight labels in colors for each group of metrics.
        Default is False (names written in black).
    :param cfram: boolean, optional
        True to color the frame for each group of metrics.
        Default is False (frame is black).
    :param bold_names: list, optional
        List of names to write in bold (must be in xy_names).
        Default is False (names written in normal font).
    :param save_eps: boolean, optional
        True to save the plot in eps format instead of png.
        Default is False (plot saved in png format).

    Output:
    ------
    """
    met_o1, met_o2, met_o3, met_o4 = return_metrics_type()
    if tab_pval is not None:
        mask1 = NUMPYmasked_where(tab_rval < tab_pval, tab_rval).mask.astype("f")
        mask2 = NUMPYmasked_where(tab_rval > -tab_pval, tab_rval).mask.astype("f")
        mask = mask1 + mask2
        tab_plot = NUMPYmasked_where(mask == 2, tab_rval)
    else:
        tab_plot = deepcopy(tab_rval)
    fig, ax = plt.subplots(figsize=(0.5 * len(tab_rval), 0.5 * len(tab_rval)))
    # shading & colorbar
    levels = MaxNLocator(nbins=20).tick_values(-1, 1)
    cmap = plt.get_cmap("cmo.balance")
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cs = plt.pcolormesh(tab_plot, cmap=cmap, norm=norm)
    # title
    plt.title(title, fontsize=30, y=1.01, loc="center")
    # x axis
    xticks = [ii + 0.5 for ii in range(len(tab_plot))]
    xlabel = [met for met in xy_names]
    ax.set_xticks(xticks)
    ax.set_xticklabels([""] * len(xticks))
    ax.tick_params(axis="x", labelsize=15, labelrotation=90)
    for ll, txt in enumerate(xlabel):
        cc = "yellowgreen" if txt in met_o1 else ("plum" if txt in met_o2 else
                                                  ("gold" if txt in met_o3 else "turquoise"))
        weight = "bold" if txt in bold_names else "normal"
        if chigh is True:
            boxdict = dict(lw=0, facecolor=cc, pad=1, alpha=1)
            ax.text(ll + 0.5, -0.2, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k",
                    weight=weight, bbox=boxdict)
            ax.text(-0.4, ll + 0.5, met_names[txt], fontsize=18, ha="right", va="center", color="k", weight=weight,
                    bbox=boxdict)
        else:
            ax.text(ll + 0.5, -0.2, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k",
                    weight=weight)
            ax.text(-0.4, ll + 0.5, met_names[txt], fontsize=18, ha="right", va="center", color="k", weight=weight)
    if cfram is True:
        tmp1 = [txt for ll, txt in enumerate(xlabel) if txt in met_o1]
        n1 = len(tmp1)
        tmp2 = [txt for ll, txt in enumerate(xlabel) if txt in met_o2]
        n2 = n1 + len(tmp2)
        tmp3 = [txt for ll, txt in enumerate(xlabel) if txt in met_o3]
        n3 = n2 + len(tmp3)
        tmp4 = [txt for ll, txt in enumerate(xlabel) if txt in met_o4]
        n4 = n3 + len(tmp4)
        lic, nc = list(), 0
        if len(tmp1) > 0:
            lic += ["yellowgreen"] * 4
            nc += 2
        if len(tmp2) > 0:
            lic += ["plum"] * 4
            nc += 2
        if len(tmp3) > 0:
            lic += ["gold"] * 4
            nc += 2
        if len(tmp4) > 0:
            lic += ["turquoise"] * 4
            nc += 2
        lic = ["k"] * nc + lic
        lis = ["-"] * len(lic)
        liw = [4] * nc + [10] * (len(lic) - nc)
        # horizontal and vertical black lines
        lix, liy = list(), list()
        if len(tmp1) > 0:
            lix += [[n1, n1], [0, n4]]
            liy += [[0, n4], [n1, n1]]
        if len(tmp2) > 0:
            lix += [[n2, n2], [0, n4]]
            liy += [[0, n4], [n2, n2]]
        if len(tmp3) > 0:
            lix += [[n3, n3], [0, n4]]
            liy += [[0, n4], [n3, n3]]
        if len(tmp4) > 0:
            lix += [[n1, n1], [0, n4]]
            liy += [[0, n4], [n1, n1]]
        # horizontal and vertical colored frame
        if len(tmp1) > 0:
            lix += [[0, 0], [n4, n4], [0, n1], [0, n1]]
            liy += [[0, n1], [0, n1], [0, 0], [n4, n4]]
        if len(tmp2) > 0:
            add = 0.18 if len(lix) > 0 else 0
            lix += [[0, 0], [n4, n4], [n1 + add, n2], [n1 + add, n2]]
            liy += [[n1 + add, n2], [n1 + add, n2], [0, 0], [n4, n4]]
        if len(tmp3) > 0:
            add = 0.18 if len(lix) > 0 else 0
            lix += [[0, 0], [n4, n4], [n2 + add, n3], [n2 + add, n3]]
            liy += [[n2 + add, n3], [n2 + add, n3], [0, 0], [n4, n4]]
        if len(tmp4) > 0:
            add = 0.18 if len(lix) > 0 else 0
            lix += [[0, 0], [n4, n4], [n3 + add, n4], [n3 + add, n4]]
            liy += [[n3 + add, n4], [n3 + add, n4], [0, 0], [n4, n4]]
        for lc, ls, lw, lx, ly in zip(lic, lis, liw, lix, liy):
            line = Line2D(lx, ly, c=lc, lw=lw, ls=ls, zorder=10)
            line.set_clip_on(False)
            ax.add_line(line)
    # y axis
    ax.set_yticks(xticks)
    ax.set_yticklabels([""] * len(xticks))
    ax.tick_params(axis="y", labelsize=15)
    # text
    if write_corr is True:
        for ii in range(len(tab_plot)):
            for jj in range(len(tab_plot)):
                plt.text(ii + 0.5, jj + 0.5, str(round(tab_rval[ii, jj], 1)), fontsize=10, ha="center", va="center")
    if tab_pval is not None:
        ax.text(len(tab_plot) + 1, 0, "nbr corr < 0", fontsize=15, ha="right", va="top", rotation=45)
        ax.text(len(tab_plot) + 2, 0, "nbr corr > 0", fontsize=15, ha="right", va="top", rotation=45)
        for ii in range(len(tab_plot)):
            nbr1 = str(sum([1 for jj in range(len(tab_plot[ii])) if tab_plot[ii][jj] < 0]))
            nbr2 = str(sum([1 for jj in range(len(tab_plot[ii])) if tab_plot[ii][jj] > 0])).zfill(2)
            ax.text(len(tab_plot) + 0.5, ii + 0.5, nbr1, fontsize=15, ha="center", va="center")
            ax.text(len(tab_plot) + 1 + 0.5, ii + 0.5, nbr2, fontsize=15, ha="center", va="center")
    # color bar
    levels = [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)]
    x2 = ax.get_position().x1
    y1 = ax.get_position().y0
    y2 = ax.get_position().y1
    if tab_pval is not None:
        cax = plt.axes([x2 + 0.09, y1, 0.02, y2 - y1])
    else:
        cax = plt.axes([x2 + 0.02, y1, 0.02, y2 - y1])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend="both", aspect=40)
    cbar.ax.tick_params(labelsize=18)
    # save fig
    if save_eps is True:
        plt.savefig(figure_name + ".eps", bbox_inches="tight", format="eps")
    else:
        plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_colorbar(cmap, mini=-1., maxi=1., nbins=20):
    """
    Modifies given colobar instance (removes the extremes on each side of the colorbar)

    Inputs:
    ------
    :param cmap: colormap instance
        Colormap instance as defined in matplotlib (see matplotlib.cm.get_cmap)
    **Optional arguments:**
    :param mini: float
        Minimum value of the colorbar.
    :param maxi: float
        Maximum value of the colorbar.
    :param nbins: integer
        Number of interval in the colorbar.

    Outputs:
    -------
    :return newcmp1: object
        Colormap, baseclass for all scalar to RGBA mappings
    :return norm: object
        Normalize, a class which can normalize data into the [0.0, 1.0] interval.
    """
    levels = MaxNLocator(nbins=nbins).tick_values(mini, maxi)
    newcmp1 = cmap(NUMPYlinspace(0.15, 0.85, 256))
    newcmp2 = cmap(NUMPYlinspace(0.0, 1.0, 256))
    newcmp1 = ListedColormap(newcmp1)
    newcmp1.set_over(newcmp2[-30])
    newcmp1.set_under(newcmp2[29])
    newcmp1.set_bad(color="k")  # missing values in black
    norm = BoundaryNorm(levels, ncolors=newcmp1.N)
    return newcmp1, norm


def plot_portraitplot(tab, figure_name, xticklabel=[], yticklabel=[], title=[], write_metrics=False, my_text="",
                      levels=None, cfram=False, chigh=False, save_eps=False, nbr_space=2):
    """
    Plot the portraitplot (as in BAMS paper)

    Inputs:
    ------
    :param tab: list of masked_array
        List of masked_array containing metric collections values.
    :param figure_name: string
        Name of the output figure.
    **Optional arguments:**
    :param xticklabel: list of string, optional
        List of the names to put as x ticks (metric names).
        Default is [] nothing will be written).
    :param yticklabel: list of string, optional
        List of the names to put as y ticks (model names).
        Default is [] nothing will be written).
    :param title: list of string, optional
        List of metric collection's title.
        Default is [], no title will be written).
    :param write_metrics: boolean, optional
        True to write the metric value in each box.
        Default is False (value not written).
    :param my_text: string, optional
        Text to add at the bottom right of the plot (I use it to indicate how CMIP6 models are marked in the plot).
        Default is "" (no text will be written).
    :param levels: list of floats, optional
        Levels of the colorbar, if None is given, colobar ranges from -1 to 1.
        Default is None ([-1.0, -0.5, 0.0, 0.5, 1.0] will be used).
    :param chigh: boolean, optional
        True to highlight labels in colors for each group of metrics.
        Default is False (names written in black).
    :param cfram: boolean, optional
        True to color the frame for each group of metrics.
        Default is False (frame is black).
    :param save_eps: boolean, optional
        True to save the plot in eps format instead of png.
        Default is False (plot saved in png format).
    :param nbr_space: integer, optional
        Number of blank space between two metric collections.
        Default is 2.

    Output:
    ------
    """
    met_o1, met_o2, met_o3, met_o4 = return_metrics_type()
    if levels is None:
        levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    fontdict = {"fontsize": 40, "fontweight": "bold"}
    # nbr of columns of the portraitplot
    nbrc = sum([len(tab[ii][0]) for ii in range(len(tab))]) + (len(tab) - 1) * nbr_space
    # figure definition
    fig = plt.figure(0, figsize=(0.5 * nbrc, 0.5 * len(tab[0])))
    gs = GridSpec(1, nbrc)
    # adapt the colorbar
    cmap, norm = my_colorbar(plt.get_cmap("cmo.balance"), mini=min(levels), maxi=max(levels))
    # loop on metric collections
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
        except: pass
        # x axis
        ticks = [ii + 0.5 for ii in range(len(tmp[0]))]
        ax.set_xticks(ticks)
        ax.set_xticklabels([] * len(ticks))
        if len(xticklabel[kk]) == len(tmp[0]):
            for ll, txt in enumerate(xticklabel[kk]):
                cc = "yellowgreen" if txt in met_o1 else ("plum" if txt in met_o2 else
                                                          ("gold" if txt in met_o3 else "turquoise"))
                if chigh is True:
                    boxdict = dict(lw=0, facecolor=cc, pad=1, alpha=1)
                    ax.text(ll + 0.5, -0.2, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k",
                            bbox=boxdict)
                else:
                    ax.text(ll + 0.5, -0.2, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k")
        if cfram is True:
            tmp1 = [met_o1, met_o2, met_o3, met_o4]
            # draw vertical black lines to separate metric types
            nn = 0
            lix = [[0, 0]]
            for tt in tmp1:
                tmp2 = [txt for ll, txt in enumerate(xticklabel[kk]) if txt in tt]
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
            # draw horizontal colored lines to indicate metric types
            nn = 0
            lic, lix = list(), list()
            for uu, tt in enumerate(tmp1):
                tmp2 = [txt for ll, txt in enumerate(xticklabel[kk]) if
                        txt in tt or txt + "_1" in tt or txt + "_2" in tt]
                if len(tmp2) > 0:
                    cc = "yellowgreen" if uu == 0 else ("plum" if uu == 1 else ("gold" if uu == 2 else "turquoise"))
                    lic += [cc, cc]
                    if nn > 0:
                        lix += [[nn + 0.2, nn + len(tmp2)], [nn + 0.2, nn + len(tmp2)]]
                    else:
                        lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                    nn += len(tmp2)
                    del cc
                del tmp2
            liy = [[len(tab[0]), len(tab[0])], [0, 0]] * (len(lix) // 2)
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
        ticks = [ii + 0.5 for ii in range(len(tmp))]
        ax.set_yticks(ticks)
        if kk != 0 or len(yticklabel) != len(tmp):
            ax.set_yticklabels([""] * len(ticks))
        else:
            ax.text(-5 * dx, -1 * dy, my_text, fontsize=25, ha="right", va="top", transform=ax.transAxes)
            ax.tick_params(axis="y", labelsize=20)
            ax.set_yticklabels(yticklabel)
            ax.yaxis.set_label_coords(-20 * dx, 0.5)
        # grid (create squares around metric values)
        for ii in range(1, len(tmp)):
            ax.axhline(ii, color="k", linestyle="-", linewidth=1)
        for ii in range(1, len(tmp[0])):
            ax.axvline(ii, color="k", linestyle="-", linewidth=1)
        # write metric value in each square (standardized value!)
        if write_metrics is True:
            for jj in range(len(tmp[0])):
                for ii in range(len(tmp)):
                    if tmp.mask[ii, jj] == False:
                        plt.text(jj + 0.5, ii + 0.5, str(round(tmp[ii, jj], 1)), fontsize=10, ha="center", va="center")
        if kk == len(tab) - 1:
            x2 = ax.get_position().x1
            y1 = ax.get_position().y0
            y2 = ax.get_position().y1
        count += len(tmp[0]) + nbr_space
    # color bar
    cax = plt.axes([x2 + 0.03, y1, 0.02, y2 - y1])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=levels, pad=0.05, extend="both", aspect=40)
    cbar.ax.set_yticklabels(["-2 $\sigma$", "-1", "MMV", "1", "2 $\sigma$"], fontdict=fontdict)
    dict_arrow = dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0)
    dict_txt = dict(fontsize=40, rotation="vertical", ha="center", weight="bold")
    cax.annotate("", xy=(3.7, 0.06), xycoords="axes fraction", xytext=(3.7, 0.45), arrowprops=dict_arrow)
    cax.text(5.2, 0.45, "closer to reference", va="top", **dict_txt)
    cax.annotate("", xy=(3.7, 0.94), xycoords="axes fraction", xytext=(3.7, 0.55), arrowprops=dict_arrow)
    cax.text(5.2, 0.55, "further from reference", va="bottom", **dict_txt)
    # save fig
    if save_eps is True:
        plt.savefig(figure_name + ".eps", bbox_inches="tight", format="eps")
    else:
        plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def plot_projects_comparison(tab_val, figure_name, title="", x_name="", y_name="", xticklabel=[], yticklabel="",
                             colors=None, tab_bst=None, legend=None, chigh=False, cfram=False,
                             bold_names=[], save_eps=False):
    """
    Plots the projects comparison, markers for each metric, solid markers mean that the difference if significant at the
    95% confidence level

    Inputs:
    ------
    :param tab_val: list
        2D array with metric values averaged by projects, projects in the first axis, metrics in the second.
    :param figure_name: string
        Path to and name of the output plot.
    **Optional arguments:**
    :param title: string, optional
        Title of the plot.
        Default is "" (no title).
    :param x_name: string, optional
        Names of x axis.
        Default is "" (no name written).
    :param y_name: string, optional
        Names of y axis.
        Default is "" (no name written).
    :param xticklabel: list of string, optional
        List of the names to put as x ticks.
        Default is [] (numbers from 0 to len(tab_val[0]) will be written).
    :param yticklabel: string, optional
        Name of the group used to normalized plot.
        Default is "" (no name written).
    :param colors: list of string, optional
        List of colors (e.g., "k", "r") must be the size of the first dimension of tab_val.
        Default is None (every group will be plotted in black).
    :param tab_bst: list, optional
        3D array with confidence interval on the metric values averaged by projects, projects in the first axis, metrics
        in the second, confidence interval in the third.
        Default is None (no confidence interval plotted).
    :param legend: list of string, optional
        Name of the groups to plot the legend.
        Default is None (no legend plotted).
    :param chigh: boolean, optional
        True to highlight labels in colors for each group of metrics.
        Default is False (names written in black).
    :param cfram: boolean, optional
        True to color the frame for each group of metrics.
        Default is False (frame is black).
    :param bold_names: list, optional
        List of names to write in bold (must be in xy_names).
        Default is False (names written in normal font).
    :param save_eps: boolean, optional
        True to save the plot in eps format instead of png.
        Default is False (plot saved in png format).

    Output:
    ------
    :return:
    """
    met_o1, met_o2, met_o3, met_o4 = return_metrics_type()
    fig, ax = plt.subplots(figsize=(0.5 * len(tab_val[0]), 4))
    # title
    plt.title(title, fontsize=20, y=1.01, loc='left')
    # x axis
    axis = list(range(len(tab_val[0])))
    ax.set_xticks(axis)
    if isinstance(xticklabel, list):
        ax.set_xticklabels([""] * len(xticklabel))
        for ll, txt in enumerate(xticklabel):
            cc = "yellowgreen" if txt in met_o1 else ("plum" if txt in met_o2 else
                                                      ("gold" if txt in met_o3 else "turquoise"))
            if chigh is True:
                boxdict = dict(lw=0, facecolor=cc, pad=1, alpha=1)
                ax.text(ll + 0.5, -0.03, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k",
                        bbox=boxdict)
            else:
                ax.text(ll + 0.5, -0.03, met_names[txt], fontsize=18, ha="right", va="top", rotation=45, color="k")
    else:
        ax.set_xticklabels(axis)
    if cfram is True:
        nn = 0
        lic, lix = list(), list()
        for cc, tmp1 in zip(["yellowgreen", "plum", "gold", "turquoise"], [met_o1, met_o2, met_o3, met_o4]):
            tmp2 = [txt for ll, txt in enumerate(xticklabel) if txt in tmp1]
            if len(tmp2) > 0:
                lic += [cc, cc]
                if nn == 0:
                    lix += [[-0.4, len(tmp2) + 0.5], [-0.4, len(tmp2) - 0.5]]
                    nn += len(tmp2) - 0.5
                elif nn + len(tmp2) > len(xticklabel) - 1:
                    lix += [[nn, nn + len(tmp2) - 0.1], [nn, nn + len(tmp2) - 0.1]]
                else:
                    lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                    nn += len(tmp2)
        lis = ["-"] * len(lic)
        liw = [5] * len(lic)
        liy = [[2, 2], [0, 0]] * (len(lic) // 2)
        for lc, ls, lw, lx, ly in zip(lic, lis, liw, lix, liy):
            line = Line2D(lx, ly, c=lc, lw=lw, ls=ls, zorder=10)
            line.set_clip_on(False)
            ax.add_line(line)
    ax.set_xlim([min(axis) - 0.5, max(axis) + 0.5])
    ax.tick_params(axis="x", labelsize=12, labelrotation=90)
    ax.set_xlabel(x_name, fontsize=15)
    # y axis
    ax.set_yticks([0.5, 1.5], minor=True)
    ax.set_yticks([0, 1, 2], minor=False)
    ax.set_yticklabels(["reference", yticklabel, "2 * " + yticklabel],
                       fontdict={"fontsize": 12, "fontweight": "normal"})
    ax.set_ylim([0, 2])
    ax.set_ylabel(y_name, fontsize=15)
    # plot marker
    for ii in range(len(tab_val)):
        col = colors[len(colors) - 1 - ii] if isinstance(colors, list) else "k"
        ind = len(tab_val) - 1 - ii
        if tab_bst is not None:
            for jj in range(len(tab_bst[ind])):
                tmp1, tmp2 = tab_val[ind][jj], tab_bst[ind][jj]
                if ind == 0:
                    tmp3, tmp4 = tab_val[1][jj], tab_bst[1][jj]
                else:
                    tmp3, tmp4 = tab_val[0][jj], tab_bst[0][jj]
                if (min(tmp4) <= tmp1 <= max(tmp4)) or (min(tmp2) <= tmp3 <= max(tmp2)):
                    ax.plot([jj], [tmp1], markersize=13, color="none", marker="D", fillstyle="none",
                            markeredgecolor=col, markeredgewidth=3, zorder=2)
                else:
                    ax.scatter(jj, tmp1, s=200, c=col, marker="D", zorder=2)
                if tmp2[0] > 0 and tmp2[1] > 0:
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tmp2[0], tmp2[0]], c=col, lw=2, zorder=3))
                    ax.add_line(Line2D([jj - 0.3, jj + 0.3], [tmp2[1], tmp2[1]], c=col, lw=2, zorder=3))
                    tmpl = [jj - 0.05, jj - 0.05] if ii == 0 else [jj + 0.05, jj + 0.05]
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
            col = colors[len(colors) - 1 - ii] if isinstance(colors, list) else "k"
            font = {"color": col, "weight": "bold", "size": 15}
            ax.text(x2 - 2 * dx, y2 - (ii + 1) * 8 * dy, legend[len(legend) - 1 - ii], ha="right", va="center",
                    fontdict=font)
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
    if save_eps is True:
        plt.savefig(figure_name + ".eps", bbox_inches="tight", format="eps")
    else:
        plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return
