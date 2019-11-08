# -*- coding:UTF-8 -*-
import cmocean
from copy import deepcopy
from math import ceil as MATHceil
from math import floor as MATHfloor
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import array as NUMPYarray
from numpy import meshgrid as NUMPYmeshgrid
# ENSO
from EnsoPlotToolsLib import create_labels, create_levels, format_metric, minmax_plot, read_diag, read_var

colors_sup = ["r", "lime", "peru", "gold", "forestgreen", "sienna", "gold"]

mod_nicknames = ["CM5", "CM6", "CM6_atm5"]

def my_boxplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
               metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    title = dict_param["title"]
    if isinstance(title, str):
        title = [title] * nbr_panel
    yname = dict_param["yname"]
    if isinstance(yname, str):
        yname = [yname] * nbr_panel
        one_yaxis = True
    else:
        one_yaxis = False
    if isinstance(filename_nc, str):
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        if isinstance(filename_nc, str):
            legend = ["ref: " + obsname, model]
        else:
            legend = ["ref: " + obsname] + mod_nicknames # model
    nbrl = int(round(nbr_panel / 2.))
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    legco = ["k", "dodgerblue"]
    if isinstance(filename_nc, list):
        for ii in range(len(filename_nc)-1):
            legco.append(colors_sup[ii])
    lines = [Line2D([0], [0], marker="o", c="w", markerfacecolor=cc, markersize=12) for cc in legco]
    if isinstance(filename_nc, str):
        tmp = tab_mod + tab_obs
    else:
        tmp = deepcopy(tab_obs)
        for ii in range(len(tab_mod)):
            tmp += tab_mod[ii]
    tick_labels = minmax_plot(tmp, metric=plot_metric)
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        # title
        ax.set_title(title[ii], fontsize=15, y=1.01, loc="left")
        # x axis
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # y axis
        ax.set_yticks(tick_labels)
        ax.set_ylim(ymin=min(tick_labels), ymax=max(tick_labels))
        if (one_yaxis is True and ii % 2 == 0) or one_yaxis is False:
            ylabel = yname[ii]
            if units != "":
                ylabel = ylabel + " (" + units + ")"
            ax.set_ylabel(ylabel, fontsize=15)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # boxplots
        boxproperties = {
            "boxprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "capprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=legco[0], markerfacecolor=legco[0],
                               markeredgewidth=0),
            "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=legco[0], markerfacecolor=legco[0],
                              markeredgewidth=0),
            "medianprops": dict(linestyle="-", linewidth=2, color=legco[0]),
            "whiskerprops": dict(linestyle="-", linewidth=2, color=legco[0]),
        }
        if isinstance(filename_nc, basestring):
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
            ax.boxplot([[1e20, 1e20], tab_mod[ii]], whis=[5, 95], labels=["", ""], showmeans=True, showfliers=False,
                       **boxproperties)
            # my text
            if plot_metric is True:
                # relative space
                x1, x2 = ax.get_xlim()
                dx = (x2 - x1) / 100.
                y1, y2 = ax.get_xlim()
                dy = (y2 - y1) / 100.
                txt = format_metric(metric_type, metval, metric_units)
                ax.text(x2 - (2 * dx), y2 - (5 * dy), txt, fontsize=12, color="k", horizontalalignment="right",
                        verticalalignment="center")
            # legend
            if (nbr_panel == 1 and ii == 0) or (nbr_panel != 1 and ii == 1):
                ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        else:
            tmp = [tab_obs[ii]] + [1e20, 1e20] * len(tab_mod)
            ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=False, **boxproperties)
            for kk in range(len(tab_mod)):
                boxproperties = {
                    "boxprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                    "capprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                    "flierprops": dict(marker="o", markersize=2.0, markeredgecolor=legco[kk+1],
                                       markerfacecolor=legco[kk+1], markeredgewidth=0),
                    "meanprops": dict(marker="D", markersize=8.0, markeredgecolor=legco[kk+1],
                                      markerfacecolor=legco[kk+1], markeredgewidth=0),
                    "medianprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                    "whiskerprops": dict(linestyle="-", linewidth=2, color=legco[kk+1]),
                }
                tmp = [[1e20, 1e20]] * (kk + 1) + [tab_mod[kk][ii]] + [[1e20, 1e20]] * (len(tab_mod) - 1 - kk)
                ax.boxplot(tmp, whis=[5, 95], labels=[""] * len(tmp), showmeans=True, showfliers=False, **boxproperties)
            # legend
            if (nbr_panel == 1 and ii == 0) or (nbr_panel != 1 and ii == 1):
                if plot_metric is True:
                    for jj in range(1, len(legend)):
                        legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj-1]) + " " + metric_units + ")"
                ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        # grid
        ax.grid(linestyle="--", linewidth=1, which="major", axis="y")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_curve(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
             metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    if isinstance(variables, str):
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname =\
        read_var(deepcopy(variables), filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None and isinstance(filename_nc, str) is True:
        plot_metric = True
    else:
        plot_metric = False
    if isinstance(filename_nc, str):
        axis = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].coords.keys()[0]]))
    else:
        axis = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].coords.keys()[0]]))
    # figure initialization
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    if isinstance(filename_nc, str):
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if units != "":
        yname = yname + " (" + units + ")"
    if "colors" in dict_param.keys():
        linecolors = dict_param["colors"]
    else:
        linecolors = {"model": ["r"], "reference": ["k"]}
    if "linestyles" in dict_param.keys():
        linestyles = dict_param["linestyles"]
    else:
        linestyles = {"model": ["-"], "reference": ["-"]}
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        if isinstance(filename_nc, str):
            legend = ["ref: " + obsname, model]
        else:
            legend = ["ref: " + obsname] + mod_nicknames # model
    nbr = len(tab_mod[0]) if isinstance(filename_nc, str) else len(tab_mod[0][0])
    if isinstance(filename_nc, str):
        tmp = tab_mod + tab_obs
    else:
        tmp = deepcopy(tab_obs)
        for ii in range(len(tab_mod)):
            tmp += tab_mod[ii]
    ytick_labels = minmax_plot(tmp, metric=plot_metric)
    # plot
    if isinstance(filename_nc, str):
        if xname == "months" and nbr == 72:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig, ax = plt.subplots(figsize=(4, 4))
        plot_curve(tab_mod, tab_obs, ax, title, axis, xname, yname, ytick_labels, linecolors, linestyles, metric_type,
                   metval, metric_units, model=model, obsname=obsname, legend=legend, multimodel=False,
                   plot_metric=plot_metric)
    else:
        if metric_type is not None:
            plot_metric = True
        if xname == "months" and nbr == 72:
            nbrl = nbr_val
            nbrc = 1
            fig, axes = plt.subplots(nbrl, nbrc, figsize=(8, 4 * nbrl), sharex="col", sharey="row")
        else:
            nbrl = int(round(nbr_val / 2.))
            nbrc = 1 if nbr_val == 1 else 2
            fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
        old_leg = deepcopy(legend)
        legend = ["ref: " + obsname] + mod_nicknames # model
        for kk in range(len(tab_obs)):
            tab_tmp = [tab_mod[jj][kk] for jj in range(len(tab_mod))]
            if nbr_val == 1:
                ax = axes
            elif (nbrl == 1 and nbrc != 1) or (nbrl != 1 and nbrc == 1):
                ax = axes[kk % 2]
            else:
                ax = axes[kk / 2, kk % 2]
            if "legend" in dict_param.keys():
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
                       metval, metric_units, model=model, obsname=obsname, legend=legend, multimodel=True,
                       plot_metric=plot_metric, plot_legend=plot_legend)
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()
    return


def my_dotplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
               metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    diag_mod, diag_obs, metval, obsname =\
        read_diag(diagnostic_values, metric_values, model, reference, metric_variables)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization if isinstance(filename_nc, str):
    title = dict_param["title"]
    yname = dict_param["yname"]
    if diagnostic_units != "":
        yname = yname + " (" + diagnostic_units + ")"
    if len(yname) < 20:
        for kk in regions.keys():
            if kk in yname.lower():
                yname = regions[kk] + " " + yname
    if "colors" in dict_param.keys():
        mcolors = dict_param["colors"]
    else:
        mcolors = ["k", "dodgerblue"]
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                mcolors.append(colors_sup[ii])
    if "markers" in dict_param.keys():
        markers = dict_param["markers"]
    else:
        markers = ["D", "o"]
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                markers.append("o")
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        if isinstance(filename_nc, str):
            legend = ["ref: " + obsname, model]
        else:
            legend = ["ref: " + obsname] + mod_nicknames # model
    fig, ax = plt.subplots(figsize=(4, 4))
    if isinstance(filename_nc, str):
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
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    # y axis
    tick_labels = minmax_plot(tab, metric=plot_metric)
    plt.yticks(tick_labels, tick_labels)
    plt.ylim(min(tick_labels), max(tick_labels))
    ax.set_ylabel(yname, fontsize=15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if min(tick_labels) < 0 and max(tick_labels) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # dots
    for ii in range(len(tab)):
        ax.scatter([ii], tab[ii], s=80, c=mcolors[ii], marker=markers[ii], clip_on=False)
    if isinstance(filename_nc, str):
        # my text
        if plot_metric is True:
            x1, x2 = ax.get_xlim()
            dx = (x2 - x1) / 100.
            y1, y2 = ax.get_ylim()
            dy = (y2 - y1) / 100.
            txt = format_metric(metric_type, metval, metric_units)
            plt.text(x2 - (2 * dx), y2 - (5 * dy), txt, fontsize=12, color='k', horizontalalignment='right',
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
    plt.grid(linestyle='--', linewidth=1, which='major')
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()
    return


def my_hovmoeller(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
                  metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    if isinstance(variables, str):
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if isinstance(filename_nc, str):
        tim = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
        lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
    else:
        tim = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[0]]))
        lon = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[1]]))
    nbr_years = len(tim) / 12.
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if isinstance(filename_nc, list):
        nbr_panel = nbr_panel + ((len(tab_mod) - 1) * nbr_val)
    title = dict_param["title"]
    if isinstance(filename_nc, list):
        if nbr_val == 1:
            title = [obsname] + mod_nicknames # model[ii]
        else:
            title = [obsname] + mod_nicknames # list()
            title = title * nbr_val
            # for tt in dict_param["title"]:
            #     title.append(tt + ": " + obsname)
            #     for mod in mod_nicknames:#model:
            #         title.append(tt + ": " + mod)
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    if isinstance(filename_nc, str):
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if units != "":
        zname = zname + " (" + units + ")"
    colorbar = "cmo." + dict_param["colorbar"]
    labelbar = dict_param["label"]
    nbrl = int(round(nbr_panel / 2.))
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbr_years * nbrl), sharex="col", sharey="row")
    hspa1 = 0.3 / nbr_years
    hspa2 = 0.01 / nbr_years
    plt.subplots_adjust(hspace=hspa1, wspace=0.1)
    xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    xlabel_ticks, xlabel = create_labels(xname, xlabel_ticks)
    ylabel_ticks = list(range(int(MATHfloor(min(tim))), int(MATHceil(max(tim))) + 1))
    ylabel_ticks, ylabel = create_labels(yname, ylabel_ticks)
    tab = list()
    if isinstance(filename_nc, str):
        legend = ["ref: " + obsname, model]
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            tab.append(tab_mod[ii])
    else:
        legend = ["ref: " + obsname] + model
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            for kk in range(len(tab_mod)):
                tab.append(tab_mod[kk][ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        # title
        if isinstance(filename_nc, str):
            ax.set_title(title[ii / 2], fontsize=15, y=1. + hspa2, loc="left")
        else:
            if ii % 2 == 0:
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="left")
            else:
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="right")
            if nbr_val > 1:
                if ii % (len(filename_nc) + 1) == 0:
                    ttx2 = ax.get_position().x1
                elif (ii - 1) % (len(filename_nc) + 1) == 0:
                    ttx1 = ax.get_position().x0
                    tty2 = ax.get_position().y1
                    ax.text(ttx2 + (ttx1 - ttx2) / 2., tty2 + hspa2 / 2,
                            dict_param["title"][(ii - 1) / (len(filename_nc) + 1)], fontsize=15, weight="bold",
                            horizontalalignment="center", verticalalignment="center", transform=fig.transFigure)
        if isinstance(filename_nc, str):
            if ii in [0, 1]:
                ax.text(0.5, 1. + hspa1, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes)
        # x axis
        ax.set_xlim(xmin=min(lon), xmax=max(lon))
        if ii >= nbr_panel - 2:
            ax.set_xticks(xlabel_ticks)
            ax.set_xticklabels(xlabel)
            ax.set_xlabel(xname, fontsize=15)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # y axis
        ax.set_ylim(ymin=min(tim), ymax=max(tim))
        if ii % 2 == 0:
            ax.set_yticks(ylabel_ticks)
            ax.set_yticklabels(ylabel)
            ax.set_ylabel(yname, fontsize=15)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # hovmoeller
        levels = create_levels(labelbar)
        xx, yy = NUMPYmeshgrid(lon, tim)
        cs = ax.contourf(xx, yy, tab[ii], levels=levels, extend="both", cmap=colorbar)
        if ii == nbr_panel - 2:
            x1 = ax.get_position().x0
        elif ii == nbr_panel - 1:
            x2 = ax.get_position().x1
    # add colorbar
    if nbr_years == 1 and nbrl == 1:
        cax = plt.axes([x1, -0.1, x2 - x1, 0.04])
    else:
        if nbr_panel in [5, 6]:
            if nbr_years == 6:
                cax = plt.axes([x1, 0.09, x2 - x1, 0.005])
            else:
                cax = plt.axes([x1, 0.03, x2 - x1, 0.02])
        elif nbr_panel in [3, 4] and nbr_years == 6:
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
            cax = plt.axes([x1, 0.0, x2-x1, 0.02])
    cbar = plt.colorbar(cs, cax=cax, orientation="horizontal", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(zname, fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_map(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
           metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    if isinstance(variables, str):
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    if metric_units is not None:
        metric_units = metric_units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if isinstance(filename_nc, str):
        lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
        lat = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
    else:
        lon = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[1]]))
        lat = list(NUMPYarray(tab_mod[0][0].coords[tab_mod[0][0].dims[0]]))
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if isinstance(filename_nc, list):
        nbr_panel = nbr_panel + ((len(tab_mod) - 1) * nbr_val)
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    if isinstance(filename_nc, list):
        if nbr_val == 1:
            title = [obsname] + mod_nicknames # model
        else:
            title = [obsname] + mod_nicknames #list()
            title = title * nbr_val
            # for tt in dict_param["title"]:
            #     title.append(tt + ": " + obsname)
            #     for mod in mod_nicknames:#model:
            #         title.append(tt + ": " + mod)
    if isinstance(filename_nc, str):
        units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    else:
        units = tab_mod[0][0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if units != "":
        zname = zname + " (" + units + ")"
    colorbar = "cmo." + dict_param["colorbar"]
    labelbar = dict_param["label"]
    if "maskland" in dict_param.keys():
        maskland = dict_param["maskland"]
    else:
        maskland = False
    nbrl = int(round(nbr_panel / 2.))
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    hspa1 = 0.1
    hspa2 = 0.01
    if nbr_panel in [3, 4]:
        plt.subplots_adjust(hspace=-0.75, wspace=0.2)
    elif nbr_panel in [5, 6]:
        plt.subplots_adjust(hspace=-0.8, wspace=0.2)
    elif nbr_panel in [7, 8]:
        plt.subplots_adjust(hspace=-0.85, wspace=0.2)
    elif nbr_panel in [11, 12]:
        plt.subplots_adjust(hspace=-0.9, wspace=0.2)
    else:
        plt.subplots_adjust(hspace=hspa1, wspace=0.2)
    xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    xlabel_ticks, xlabel = create_labels(xname, xlabel_ticks)
    ylabel_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
    ylabel_ticks, ylabel = create_labels(yname, ylabel_ticks)
    tab = list()
    if isinstance(filename_nc, str):
        legend = ["ref: " + obsname, model]
        for ii in range(len(tab_mod)):
            tab.append(tab_obs[ii])
            tab.append(tab_mod[ii])
    else:
        legend = ["ref: " + obsname] + model
        for ii in range(len(tab_obs)):
            tab.append(tab_obs[ii])
            for kk in range(len(tab_mod)):
                tab.append(tab_mod[kk][ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        if isinstance(filename_nc, list):
            if ii % 2 == 0:
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="left")
            else:
                ax.set_title(title[ii], fontsize=15, y=1. + hspa2, loc="right")
            if nbr_val > 1:
                if ii % (len(filename_nc) + 1) == 0:
                    ax.text(1.1, 1.02, dict_param["title"][ii / (len(filename_nc) + 1)], fontsize=15,
                            weight="bold", horizontalalignment="center", verticalalignment="bottom",
                            transform=ax.transAxes)
        else:
            ax.set_title(title[ii / 2], fontsize=15, y=1. + hspa2, loc="left")
        if isinstance(filename_nc, str):
            if ii in [0, 1]:
                ax.text(0.5, 1. + hspa1 * 6, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes)
        # map
        xx, yy = NUMPYmeshgrid(lon, lat)
        if lat[-1] - lat[0] < 40:
            locmap = Basemap(projection="cyl", llcrnrlat=lat[0] - 5, urcrnrlat=lat[-1] + 5, llcrnrlon=lon[0],
                             urcrnrlon=lon[-1], ax=ax)
        else:
            locmap = Basemap(projection="cyl", llcrnrlat=lat[0], urcrnrlat=lat[-1], llcrnrlon=lon[0], urcrnrlon=lon[-1],
                             ax=ax)
        # draw coastlines
        locmap.drawcoastlines()
        # fill continents
        if maskland is True:
            locmap.fillcontinents(color="gainsboro")
        # draw parallels
        locmap.drawparallels(ylabel_ticks, labels=[1, 0, 0, 0], fontsize=12, dashes=[3, 1], linewidth=1)
        # draw meridians
        locmap.drawmeridians(xlabel_ticks, labels=[0, 0, 0, 1], fontsize=12, dashes=[3, 1], linewidth=1)
        #cs = locmap.pcolormesh(xx, yy, tab[ii], vmin=min(labelbar), vmax=max(labelbar), cmap=colorbar)
        levels = create_levels(labelbar)
        cs = locmap.contourf(xx, yy, tab[ii], levels=levels, extend="both", cmap=colorbar)
        # my text
        if ii == 0 and plot_metric is True:
            x1 = ax.get_position().x1
            y2 = ax.get_position().y1
            if nbrl == 1:
                ax2 = axes[(ii + 1) % 2]
            else:
                ax2 = axes[(ii + 1) / 2, (ii + 1) % 2]
            x2 = ax2.get_position().x0
            txt = format_metric(metric_type, metval, metric_units)
            ax.text(x1 + ((x2 - x1) / 2.), y2 + 0.1, txt, fontsize=12, color="k", horizontalalignment="center",
                    verticalalignment="center", transform=fig.transFigure)
        if ii == nbr_panel - 2:
            x1 = ax.get_position().x0
        elif ii == nbr_panel - 1:
            x2 = ax.get_position().x1
            y1 = ax.get_position().y0
    # add colorbar
    if nbr_panel in [3, 4]:
        cax = plt.axes([x1, y1 + 0.21, x2 - x1, 0.015])
    elif nbr_panel in [5, 6]:
        cax = plt.axes([x1, y1 + 0.2, x2 - x1, 0.015])
    elif nbr_panel in [7, 8]:
        cax = plt.axes([x1, y1 + 0.21, x2 - x1, 0.01])
    elif nbr_panel in [11, 12]:
        cax = plt.axes([x1, y1 + 0.22, x2 - x1, 0.005])
    else:
        cax = plt.axes([x1, y1 + 0.2, x2 - x1, 0.03])
    cbar = plt.colorbar(cs, cax=cax, orientation="horizontal", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(zname, fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_scatterplot(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
                   metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    if isinstance(variables, str):
        nbr_val = 1
    else:
        nbr_val = len(variables)
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None and isinstance(filename_nc, str) is True:
        plot_metric = True
    else:
        plot_metric = False
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    if isinstance(filename_nc, list) and nbr_panel > 1:
        nbr_panel = nbr_panel + len(filename_nc) - 1
    title = dict_param["title"]
    if isinstance(filename_nc, list) and nbr_panel > 1:
        title = [obsname] + mod_nicknames#model
    if isinstance(title, str):
        title = [title] * nbr_panel
    xname = dict_param["xname"]
    if isinstance(xname, str):
        xname = [xname] * nbr_panel
        one_xaxis = True
    else:
        one_xaxis = False
    yname = dict_param["yname"]
    if isinstance(yname, str):
        yname = [yname] * nbr_panel
        one_yaxis = True
    else:
        one_yaxis = False
    if "colors" in dict_param.keys():
        mcolors = dict_param["colors"]
    else:
        mcolors = ["k", "dodgerblue"]
        # if dict_param["nbr_panel"] == len(nbr_val):
        #     mcolors = list(reversed(mcolors))
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                mcolors.append(colors_sup[ii])
    if "markers" in dict_param.keys():
        markers = dict_param["markers"]
    else:
        markers = ["D", "."]
        # if dict_param["nbr_panel"] == len(nbr_val):
        #     markers = list(reversed(markers))
        if isinstance(filename_nc, list):
            for ii in range(len(filename_nc) - 1):
                markers.append(".")
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        if isinstance(filename_nc, basestring):
            legend = ["ref: " + obsname, model]
        else:
            legend = ["ref: " + obsname] + mod_nicknames# model
    keys1 = ["", "_neg", "_pos"]
    keys2 = ["all", "x<0", "x>0"]
    keys3 = [[None, None], [None, 0], [0, None]]
    keys4 = ["black", "dodgerblue", "red"]
    lines = [Line2D([0], [0], marker=markers[kk], c="w", markerfacecolor=mcolors[kk], markersize=12)
             for kk in range(len(mcolors))]
    nbrl = int(round(nbr_panel / 2.))
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    plt.subplots_adjust(hspace=0.3, wspace=0.1)
    if isinstance(filename_nc, basestring):
        tab1 = tab_obs[::2] + tab_mod[::2]
        tab2 = tab_obs[1::2] + tab_mod[1::2]
    else:
        tab1 = tab_obs[::2]
        tab2 = tab_obs[1::2]
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
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        # title
        if isinstance(filename_nc, basestring):
            ax.set_title(title[ii / 2], fontsize=15, y=1.01, loc="left")
            if nbr_panel == len(tab_mod) and ii in [0, 1]:
                ax.text(0.5, 1.15, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes)
        else:
            ax.set_title(title[ii], fontsize=15, y=1.01, loc="left")
        # x axis
        ax.set_xticks(xtick_labels)
        ax.set_xlim(xmin=min(xtick_labels), xmax=max(xtick_labels))
        if (one_xaxis is True and (ii >= (nbrc * nbrl) - 2)) or one_xaxis is False:
            xlabel = xname[ii]
            for kk in regions.keys():
                if kk in xlabel.lower():
                    xlabel = regions[kk] + " " + xlabel
            units = tab1[ii].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
            if units != "":
                xlabel = xlabel + " (" + units + ")"
            ax.set_xlabel(xlabel, fontsize=15)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # y axis
        ax.set_yticks(ytick_labels)
        ax.set_ylim(ymin=min(ytick_labels), ymax=max(ytick_labels))
        if (one_yaxis is True and ii % 2 == 0) or one_yaxis is False:
            ylabel = yname[ii]
            for kk in regions.keys():
                if kk in ylabel.lower():
                    ylabel = regions[kk] + " " + ylabel
            units = tab2[ii].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
            if units != "":
                ylabel = ylabel + " (" + units + ")"
            ax.set_ylabel(ylabel, fontsize=15)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        # scatterplots and slopes
        x1, x2 = ax.get_xlim()
        dx = (x2 - x1) / 100.
        y1, y2 = ax.get_ylim()
        dy = (y2 - y1) / 100.
        if (nbr_panel > 1 and nbr_panel == nbr_val / 2.) or (nbr_panel == len(legend)):
            ax.scatter(tab1[ii], tab2[ii], s=10, c="k", marker=markers[ii])
            for jj in range(len(keys1)):
                if "slope" + keys1[jj] in tab1[ii].attrs.keys() and "intercept" + keys1[jj] in tab1[ii].attrs.keys():
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
                    ax.text(dx + x1, ((97 - 5 * jj) * dy) + y1, txt, fontsize=12, color=col, horizontalalignment="left",
                            verticalalignment="center")
        else:
            tmp = list()
            for jj in range(len(tab1)):
                col = mcolors[jj]
                ax.scatter(tab1[jj], tab2[jj], s=10, c=col, marker=markers[jj])
                if "slope" in tab1[jj].attrs.keys() and "intercept" in tab1[jj].attrs.keys():
                    slope = tab1[jj].attrs["slope"]
                    intercept = tab1[jj].attrs["intercept"]
                    xx = [x1, x2]
                    yy = [kk * slope + intercept for kk in xx]
                    ax.plot(xx, yy, lw=2, c=col)
                    tmp.append("slope = " + "{0:+.2f}".format(round(slope, 2)))
            for jj in range(len(tmp)):
                col = mcolors[len(tmp) - 1 - jj]
                ax.text(dx + x1, ((97 - 5 * jj) * dy) + y1, tmp[len(tmp) - 1 - jj], fontsize=12,
                        color=col, horizontalalignment="left", verticalalignment="center")
            if nbr_panel == 1 or ii == 1:
                if metric_type is not None:
                    for jj in range(1, len(legend)):
                        legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj - 1]) + " " + metric_units + ")"
                ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
        # my text
        if plot_metric is True:
            txt = format_metric(metric_type, metval, metric_units)
            ax.text(x2 - (2 * dx), y2 - (5 * dy), txt, fontsize=12, color="k", horizontalalignment="right",
                    verticalalignment="center")
        # grid
        ax.grid(linestyle="--", linewidth=1, which="major")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def plot_curve(tab_mod, tab_obs, ax, title, axis, xname, yname, ytick_labels, linecolors, linestyles, metric_type,
               metval, metric_units, model='', obsname='', legend=[], multimodel=False, plot_metric=False,
               plot_legend=False):
    # title
    ax.set_title(title, fontsize=15, y=1.01, loc="left")
    # x axis
    label_ticks = list(range(int(MATHfloor(min(axis))), int(MATHceil(max(axis))) + 1))
    label_ticks, label = create_labels(xname, label_ticks)
    ax.set_xticks(label_ticks)
    ax.set_xticklabels(label)
    ax.set_xlim([min(axis), max(axis)])
    ax.set_xlabel(xname, fontsize=15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    # y axis
    ax.set_yticks(ytick_labels)
    ax.set_yticklabels(ytick_labels)
    ax.set_ylim([min(ytick_labels), max(ytick_labels)])
    ax.set_ylabel(yname, fontsize=15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if min(ytick_labels) < 0 and max(ytick_labels) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # plot curves
    if len(tab_mod) + len(tab_obs) > 2:
        lw = 2
    else:
        lw = 4
    for ii, tab in enumerate(tab_mod):
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
            legco = [linecolors["reference"][0], linecolors["model"][0]]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in legco]
            ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
        else:
            legtxt = legend + [model, obsname]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in
                     linecolors["model"]]
            legls = [linestyles["model"][0], linestyles["reference"][0]]
            lines = lines + [Line2D([0], [0], c='k', lw=2, ls=ls) for ls in legls]
            ax.legend(lines, legtxt, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
        # my text
        if plot_metric is True:
            txt = format_metric(metric_type, metval, metric_units)
            plt.text(x2 - (2 * dx), y2 - (5 * dy), txt, fontsize=12, color='k', horizontalalignment='right',
                     verticalalignment='center')
    else:
        # legend
        if plot_legend is True:
            legco = linecolors["reference"] + linecolors["model"]
            lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in legco]
            if plot_metric is True:
                for jj in range(1, len(legend)):
                    legend[jj] = legend[jj] + " (" + "{0:.2f}".format(metval[jj - 1]) + " " + metric_units + ")"
            ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(linestyle='--', linewidth=1, which='major')
    return