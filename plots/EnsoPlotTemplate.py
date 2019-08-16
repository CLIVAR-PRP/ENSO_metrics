# -*- coding:UTF-8 -*-
import cmocean
from math import ceil as MATHceil
from math import floor as MATHfloor
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import array as NUMPYarray
from numpy import meshgrid as NUMPYmeshgrid
# ENSO
from EnsoPlotToolsLib import create_labels, format_metric, minmax_plot, read_diag, read_var


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
    units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        legend = ["ref: " + obsname, model]
    nbrl = int(round(nbr_panel / 2.))
    nbrc = 1 if nbr_panel == 1 else 2
    fig, axes = plt.subplots(nbrl, nbrc, figsize=(4 * nbrc, 4 * nbrl), sharex="col", sharey="row")
    # boxprops = dict(linestyle='-', linewidth=2, color='k')
    # capprops = dict(linestyle='-', linewidth=2, color='k')
    # marprops = dict(marker='o', markersize=2.0, markeredgecolor='k', markerfacecolor='k', markeredgewidth=0)
    # meaprops = dict(marker='D', markersize=8.0, markeredgecolor='k', markerfacecolor='k', markeredgewidth=0)
    # medprops = dict(linestyle='-', linewidth=2, color='k')
    # wisprops = dict(linestyle='-', linewidth=2, color='k')
    boxproperties = {
        "model": {
            "boxprops": dict(linestyle="-", linewidth=2, color="k"),
            "capprops": dict(linestyle="-", linewidth=2, color="k"),
            "flierprops": dict(marker="o", markersize=2.0, markeredgecolor="k", markerfacecolor="k", markeredgewidth=0),
            "meanprops": dict(marker="D", markersize=8.0, markeredgecolor="k", markerfacecolor="k", markeredgewidth=0),
            "medianprops": dict(linestyle="-", linewidth=2, color="k"),
            "whiskerprops": dict(linestyle="-", linewidth=2, color="k"),
        },
        "reference": {
            "boxprops": dict(linestyle="-", linewidth=2, color="r"),
            "capprops": dict(linestyle="-", linewidth=2, color="r"),
            "flierprops": dict(marker="o", markersize=2.0, markeredgecolor="r", markerfacecolor="r", markeredgewidth=0),
            "meanprops": dict(marker="D", markersize=8.0, markeredgecolor="r", markerfacecolor="r", markeredgewidth=0),
            "medianprops": dict(linestyle="-", linewidth=2, color="r"),
            "whiskerprops": dict(linestyle="-", linewidth=2, color="r"),
        },
    }
    legco = ["k", "r"]
    lines = [Line2D([0], [0], marker="o", c="w", markerfacecolor=cc, markersize=12) for cc in legco]
    tick_labels = minmax_plot(tab_mod + tab_obs, metric=plot_metric)
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
        ax.boxplot([tab_obs[ii], [1e20, 1e20]], whis=[5, 95], labels=["", ""], showmeans=True, showfliers=False,
                   **boxproperties["reference"])
        ax.boxplot([[1e20, 1e20], tab_mod[ii]], whis=[5, 95], labels=["", ""], showmeans=True, showfliers=False,
                   **boxproperties["model"])
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
        # grid
        ax.grid(linestyle="--", linewidth=1, which="major", axis="y")
    plt.savefig(figure_name, bbox_inches="tight")
    plt.close()
    return


def my_curve(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
             metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    tab_mod, tab_obs, metval, obsname =\
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    axis = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].coords.keys()[0]]))
    # figure initialization
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if units != "":
        yname = yname + " (" + units + ")"
    if "colors" in dict_param.keys():
        linecolors = dict_param["colors"]
    else:
        linecolors = {"model": ["k"], "reference": ["r"]}
    if "linestyles" in dict_param.keys():
        linestyles = dict_param["linestyles"]
    else:
        linestyles = {"model": ["-"], "reference": ["-"]}
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        legend = ["ref: " + obsname, model]
    if xname == "months" and len(tab_mod[0]) == 72:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
    # title
    ax.set_title(title, fontsize=15, y=1.01, loc="left")
    # x axis
    label_ticks = list(range(int(MATHfloor(min(axis))), int(MATHceil(max(axis))) + 1))
    label_ticks, label = create_labels(xname, label_ticks)
    plt.xticks(label_ticks, label)
    plt.xlim(min(axis), max(axis))
    ax.set_xlabel(xname, fontsize=15)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    # y axis
    tick_labels = minmax_plot(tab_mod+tab_obs, metric=plot_metric)
    plt.yticks(tick_labels, tick_labels)
    plt.ylim(min(tick_labels), max(tick_labels))
    ax.set_ylabel(yname, fontsize=15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if min(tick_labels) < 0 and max(tick_labels) > 0:
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
    # legend
    if len(linecolors["model"]) == 1 and len(linecolors["reference"]) == 1:
        # legcol = [linecolors["model"][0], linecolors["reference"][0]]
        # for ii, txt in enumerate(legend):
        #     plt.text((2 * dx) + x1, y2 - ((ii + 1) * 5 * dy), txt, fontsize=12, color=legcol[ii],
        #              horizontalalignment='left', verticalalignment='center')
        legco = [linecolors["reference"][0], linecolors["model"][0]]
        lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in legco]
        ax.legend(lines, legend, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    else:
        legtxt = legend + [model, obsname]
        lines = [Line2D([0], [0], marker='o', c='w', markerfacecolor=cc, markersize=12) for cc in linecolors["model"]]
        legls = [linestyles["model"][0], linestyles["reference"][0]]
        lines = lines + [Line2D([0], [0], c='k', lw=2, ls=ls) for ls in legls]
        ax.legend(lines, legtxt, bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    # my text
    if plot_metric is True:
        txt = format_metric(metric_type, metval, metric_units)
        plt.text(x2 - (2 * dx), y2 - (5 * dy), txt, fontsize=12, color='k', horizontalalignment='right',
                 verticalalignment='center')
    plt.grid(linestyle='--', linewidth=1, which='major')
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
    # figure initialization
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
        mcolors = ["r", "k"]
    if "markers" in dict_param.keys():
        markers = dict_param["markers"]
    else:
        markers = ["D", "o"]
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        legend = ["ref: " + obsname, model]
    fig, ax = plt.subplots(figsize=(4, 4))
    # title
    ax.set_title(title, fontsize=15, y=1.01, loc="left")
    # x axis
    label_ticks = [-0.5, 0, 1, 1.5]
    label = [""] * len(label_ticks)
    plt.xticks(label_ticks, label)
    plt.xlim(min(label_ticks), max(label_ticks))
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    # y axis
    tab = [diag_obs, diag_mod]
    tick_labels = minmax_plot(tab, metric=plot_metric)
    plt.yticks(tick_labels, tick_labels)
    plt.ylim(min(tick_labels), max(tick_labels))
    ax.set_ylabel(yname, fontsize=15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    if min(tick_labels) < 0 and max(tick_labels) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # dots
    for ii in range(2):
        ax.scatter([ii], tab[ii], s=80, c=mcolors[ii], marker=markers[ii])
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
    plt.grid(linestyle='--', linewidth=1, which='major')
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()
    return


def my_hovmoeller(model, filename_nc, dict_param, reference, metric_variables, figure_name, metric_type=None,
                  metric_values=None, metric_units=None, diagnostic_values=None, diagnostic_units=None, regions=None):
    # get data
    variables = dict_param["varpattern"]
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    tim = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
    lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
    nbr_years = len(tim) / 12.
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if units != "":
        zname = zname + " (" + units + ")"
    colorbar = "cmo." + dict_param["colorbar"]
    labelbar = dict_param["label"]
    legend = ["ref: " + obsname, model]
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
    for ii in range(len(tab_mod)):
        tab.append(tab_obs[ii])
        tab.append(tab_mod[ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        # title
        ax.set_title(title[ii / 2], fontsize=15, y=1. + hspa2, loc="left")
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
        diff = float(labelbar[1] - labelbar[0])
        if diff in [0.5, 5]:
            mult = 5
        elif diff in [0.2, 0.4, 0.6, 0.8, 1.0, 2, 4, 6, 8, 10]:
            mult = 4
        else:
            mult = 6
        delta = diff / mult
        levels = [round(kk + jj * delta, 2) for kk in labelbar[:-1] for jj in range(mult)] + [labelbar[-1]]
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
    tab_mod, tab_obs, metval, obsname = \
        read_var(variables, filename_nc, model, reference, metric_variables, metric_values)
    if metric_type is not None:
        plot_metric = True
    else:
        plot_metric = False
    lon = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[1]]))
    lat = list(NUMPYarray(tab_mod[0].coords[tab_mod[0].dims[0]]))
    # figure initialization
    nbr_panel = dict_param["nbr_panel"]
    title = dict_param["title"]
    xname = dict_param["xname"]
    yname = dict_param["yname"]
    zname = dict_param["zname"]
    units = tab_mod[0].units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
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
    else:
        plt.subplots_adjust(hspace=hspa1, wspace=0.2)
    xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    xlabel_ticks, xlabel = create_labels(xname, xlabel_ticks)
    ylabel_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
    ylabel_ticks, ylabel = create_labels(yname, ylabel_ticks)
    legend = ["ref: " + obsname, model]
    tab = list()
    for ii in range(len(tab_mod)):
        tab.append(tab_obs[ii])
        tab.append(tab_mod[ii])
    for ii in range(nbr_panel):
        if nbr_panel == 1:
            ax = axes
        elif nbrl == 1 and nbrc != 1:
            ax = axes[ii % 2]
        else:
            ax = axes[ii / 2, ii % 2]
        ax.set_title(title[ii / 2], fontsize=15, y=1. + hspa2, loc="left")
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
        cs = locmap.pcolormesh(xx, yy, tab[ii], vmin=min(labelbar), vmax=max(labelbar), cmap=colorbar)
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
        mcolors = ["k", "r"]
        if nbr_panel == len(tab_mod):
            mcolors = list(reversed(mcolors))
    if "markers" in dict_param.keys():
        markers = dict_param["markers"]
    else:
        markers = [".", "D"]
        if nbr_panel == len(tab_mod):
            markers = list(reversed(markers))
    if "legend" in dict_param.keys():
        legend = dict_param["legend"]
    else:
        legend = ["ref: " + obsname, model]
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
    xtick_labels = minmax_plot(tab_mod[::2] + tab_obs[::2], metric=False)
    ytick_labels = minmax_plot(tab_mod[1::2] + tab_obs[1::2], metric=plot_metric)
    tab1 = tab_mod[::2] + tab_obs[::2]
    tab2 = tab_mod[1::2] + tab_obs[1::2]
    if nbr_panel == len(tab_mod):
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
        ax.set_title(title[ii / 2], fontsize=15, y=1.01, loc="left")
        if nbr_panel == len(tab_mod) and ii in [0, 1]:
            ax.text(0.5, 1.15, legend[ii % 2], fontsize=15, weight="bold", horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes)
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
        if nbr_panel == len(tab_mod):
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
            for jj in range(len(tab1)):
                col = mcolors[jj]
                ax.scatter(tab1[jj], tab2[jj], s=10, c=col, marker=markers[jj])
                if "slope" in tab1[jj].attrs.keys() and "intercept" in tab1[jj].attrs.keys():
                    slope = tab1[jj].attrs["slope"]
                    intercept = tab1[jj].attrs["intercept"]
                    xx = [x1, x2]
                    yy = [kk * slope + intercept for kk in xx]
                    ax.plot(xx, yy, lw=2, c=col)
                    txt = "slope = " + "{0:+.2f}".format(round(slope, 2))
                    ax.text(dx + x1, ((97 - 5 * jj) * dy) + y1, txt, fontsize=12, color=col, horizontalalignment="left",
                            verticalalignment="center")
            if nbr_panel == 1 or ii == 1:
                ax.legend(lines, list(reversed(legend)), bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
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
