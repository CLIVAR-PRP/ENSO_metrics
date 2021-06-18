# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Provide templates for the dive down plots of each type of metric
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
import cmocean
from inspect import stack as INSPECTstack
from math import ceil as MATHceil
from math import floor as MATHfloor
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import arange as NUMPYarange
from numpy import array as NUMPYarray
from numpy import isnan as NUMPYisnan
from numpy import mean as NUMPYmean
from numpy import meshgrid as NUMPYmeshgrid
from numpy import percentile as NUMPYpercentile
import string
# ENSO_metrics functions
from EnsoMetrics.EnsoCollectionsLib import ReferenceRegions
from EnsoMetrics import EnsoErrorsWarnings
from .EnsoPlotToolsLib import create_labels, create_levels, create_lines, create_round_string
from .EnsoPlotCdatToolsLib import member_average, member_range, minimaxi, my_mask, read_data

# ---------------------------------------------------#


# ---------------------------------------------------#
# Parameters
# ---------------------------------------------------#
# colors for projects
dict_col = {"REF": "k", "CMIP": "forestgreen", "CMIP3": "orange", "CMIP5": "dodgerblue", "CMIP6": "r"}
# supplementary colors
colors_sup = ["lime", "peru", "gold", "forestgreen", "sienna", "gold"]
# line style for events
dict_sty = {"nina": "dotted", "nino": "-"}
# figure numbering
numbering = [str(ii) for ii in string.ascii_lowercase]
# metric names
met_names = {
    "BiasLhfLatRmse": "lat_LHF_bias", "BiasLhfLonRmse": "eq_LHF_bias",
    "BiasLwrLatRmse": "lat_LWR_bias", "BiasLwrLonRmse": "eq_LWR_bias",
    "BiasPrLatRmse": "double_ITCZ_bias", "BiasPrLonRmse": "eq_PR_bias",
    "BiasShfLatRmse": "lat_SHF_bias", "BiasShfLonRmse": "eq_SHF_bias",
    "BiasSshLatRmse": "lat_SSH_bias", "BiasSshLonRmse": "eq_SSH_bias",
    "BiasSstLatRmse": "lat_SST_bias", "BiasSstLonRmse": "eq_SST_bias",
    "BiasSwrLatRmse": "lat_SWR_bias", "BiasSwrLonRmse": "eq_SWR_bias",
    "BiasTauxLatRmse": "lat_Taux_bias", "BiasTauxLonRmse": "eq_Taux_bias",
    "BiasTauyLatRmse": "lat_Tauy_bias", "BiasTauyLonRmse": "eq_Tauy_bias",
    "BiasThfLatRmse": "lat_NHF_bias", "BiasThfLonRmse": "eq_NHF_bias",
    "SeasonalLhfLatRmse": "lat_LHF_sea_cycle", "SeasonalLhfLonRmse": "eq_LHF_sea_cycle",
    "SeasonalLwrLatRmse": "lat_LWR_sea_cycle", "SeasonalLwrLonRmse": "eq_LWR_sea_cycle",
    "SeasonalPrLatRmse": "double_ITCZ_sea_cycle", "SeasonalPrLonRmse": "eq_PR_sea_cycle",
    "SeasonalShfLatRmse": "lat_SHF_sea_cycle", "SeasonalShfLonRmse": "eq_SHF_sea_cycle",
    "SeasonalSshLatRmse": "lat_SSH_sea_cycle", "SeasonalSshLonRmse": "eq_SSH_sea_cycle",
    "SeasonalSstLatRmse": "lat_SST_sea_cycle", "SeasonalSstLonRmse": "eq_SST_sea_cycle",
    "SeasonalSwrLatRmse": "lat_SWR_sea_cycle", "SeasonalSwrLonRmse": "eq_SWR_sea_cycle",
    "SeasonalTauxLatRmse": "lat_Taux_sea_cycle", "SeasonalTauxLonRmse": "eq_Taux_sea_cycle",
    "SeasonalTauyLatRmse": "lat_Tauy_sea_cycle", "SeasonalTauyLonRmse": "eq_Tauy_sea_cycle",
    "SeasonalThfLatRmse": "lat_NHF_sea_cycle", "SeasonalThfLonRmse": "eq_NHF_sea_cycle",
    "EnsoLhfLonRmse": "ENSO_pattern_LHF", "EnsoLwrLonRmse": "ENSO_pattern_LWR", "EnsoPrLonRmse": "ENSO_pattern_PR",
    "EnsoShfLonRmse": "ENSO_pattern_SHF", "EnsoSshLonRmse": "ENSO_pattern_SSH", "EnsoSstLonRmse": "ENSO_pattern",
    "EnsoSwrLonRmse": "ENSO_pattern_SWR", "EnsoTauxLonRmse": "ENSO_pattern_Taux",
    "EnsoTauyLonRmse": "ENSO_pattern_Tauy", "EnsoThfLonRmse": "ENSO_pattern_NHF", "EnsoLhfTsRmse": "ENSO_lifecycle_LHF",
    "EnsoLwrTsRmse": "ENSO_lifecycle_LWR", "EnsoPrTsRmse": "ENSO_lifecycle_PR", "EnsoShfTsRmse": "ENSO_lifecycle_SHF",
    "EnsoSshTsRmse": "ENSO_lifecycle_SSH", "EnsoSstTsRmse": "ENSO_lifecycle", "EnsoSwrTsRmse": "ENSO_lifecycle_SWR",
    "EnsoTauxTsRmse": "ENSO_lifecycle_Taux", "EnsoTauyTsRmse": "ENSO_lifecycle_Tauy",
    "EnsoThfTsRmse": "ENSO_lifecycle_NHF", "EnsoAmpl": "ENSO_SST_amplitude", "EnsoSeasonality": "ENSO_SST_seasonality",
    "EnsoSstSkew": "ENSO_SST_asymmetry", "EnsoDuration": "ENSO_duration", "EnsoSstDiversity": "ENSO_SST_diversity",
    "EnsoSstDiversity_1": "ENSO_SST_diversity", "EnsoSstDiversity_2": "ENSO_SST_diversity",
    "EnsoPrMapCorr": "Dec_PR_teleconnection_CORR", "EnsoPrMapRmse": "Dec_PR_teleconnection",
    "EnsoPrMapStd": "Dec_PR_teleconnection_STD", "EnsoPrMapDjfCorr": "DJF_PR_teleconnection_CORR",
    "EnsoPrMapDjfRmse": "DJF_PR_teleconnection", "EnsoPrMapDjfStd": "DJF_PR_teleconnection_STD",
    "EnsoPrMapJjaCorr": "JJA_PR_teleconnection_CORR", "EnsoPrMapJjaRmse": "JJA_PR_teleconnection",
    "EnsoPrMapJjaStd": "JJA_PR_teleconnection_STD", "EnsoSlpMapRmse": "Dec_SLP_teleconnection",
    "EnsoSlpMapStd": "Dec_SLP_teleconnection_STD", "EnsoSlpMapDjfCorr": "DJF_SLP_teleconnection_CORR",
    "EnsoSlpMapDjfRmse": "DJF_SLP_teleconnection", "EnsoSlpMapDjfStd": "DJF_SLP_teleconnection_STD",
    "EnsoSlpMapJjaCorr": "JJA_SLP_teleconnection_CORR", "EnsoSlpMapJjaRmse": "JJA_SLP_teleconnection",
    "EnsoSlpMapJjaStd": "JJA_SLP_teleconnection_STD", "EnsoSstMapRmse": "Dec_TS_teleconnection",
    "EnsoSstMapStd": "Dec_TS_teleconnection_STD", "EnsoSstMapDjfCorr": "DJF_TS_teleconnection_CORR",
    "EnsoSstMapDjfRmse": "DJF_TS_teleconnection", "EnsoSstMapDjfStd": "DJF_TS_teleconnection_STD",
    "EnsoSstMapJjaCorr": "JJA_TS_teleconnection_CORR", "EnsoSstMapJjaRmse": "JJA_TS_teleconnection",
    "EnsoSstMapJjaStd": "JJA_TS_teleconnection_STD", "EnsoFbSstTaux": "SST-Taux_feedback",
    "EnsoFbTauxSsh": "Taux-SSH_feedback", "EnsoFbSshSst": "SSH-SST_feedback", "EnsoFbSstThf": "SST-NHF_feedback",
    "EnsodSstOce": "ocean_driven_SST", "EnsodSstOce_1": "ocean_driven_SST", "EnsodSstOce_2": "ocean_driven_SST",
    "telecon_pr_djf": "DJF_PR_region_telecon",
}


# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def plot_telecon(model, project, experiment, nc_file, dict_param, reference, figure_name, dict_diagnostic_values,
                 diagnostic_units, dict_metric_values, metric_units, metric_variables=None, metric_regions=None,
                 member=None):
    # get data
    nc_var = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk:
            nc_var += dict_param[kk]["varpattern"]
    nc_var_extra = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk:
            if "varpattern_extra" in list(dict_param[kk].keys()):
                nc_var_extra += dict_param[kk]["varpattern_extra"]
    dict_mod, dict_obs, dia_mod, dia_obs, met_val, att_glo, dict_mod_extra, dict_obs_extra = read_data(
        nc_file, nc_var, model, reference, dict_diagnostic_values, dict_metric_values, member=member,
        netcdf_var_extra=nc_var_extra)
    met_val = create_round_string(met_val)
    # get regions
    list_regions = metric_regions[metric_variables[-1]]
    # initialization of the plot
    sizex, sizey = 6, 8
    fig = plt.figure(0, figsize=(sizex, sizey))
    gs = GridSpec(sizey * 4, sizex * 4, figure=fig)
    fontsize = 10
    # ---------------------------------------------------#
    # maps
    # ---------------------------------------------------#
    pkey = "01_plot"
    colorbar = dict_param[pkey]["colorbar"]
    labelbar = dict_param[pkey]["label"]
    sizex, sizey = 10, 5
    deltx, delty = 1, 1
    counter = 0
    y_pos = 0
    for ii, vv in enumerate(dict_param[pkey]["varpattern"]):
        for jj, (d1, d2) in enumerate(zip([dict_mod, dict_obs], [model, reference])):
            # array and attributes
            tab = d1[vv]["array"]
            tax = d1[vv]["attributes"]["time_period"]
            if "arraySTD" in list(d1[vv]["attributes"].keys()):
                compare = "$\sigma$ = " + create_round_string(float(d1[vv]["attributes"]["arraySTD"]))
                if d2 == model:
                    if str(vv) + str(reference) + "_RMSE" in list(att_glo.keys()) or \
                            str(vv) + str(reference) + "_" + str(reference) + "_RMSE" in list(att_glo.keys()):
                        if str(vv) + str(reference) + "_RMSE" in list(att_glo.keys()):
                            att_na = str(vv) + str(reference) + "_RMSE"
                        else:
                            att_na = str(vv) + str(reference) + "_" + str(reference) + "_RMSE"
                        compare += "; rmse = " + create_round_string(float(att_glo[att_na]))
                        del att_na
                    if str(vv) + str(reference) + "_CORR" in list(att_glo.keys()) or \
                            str(vv) + str(reference) + "_" + str(reference) + "_CORR" in list(att_glo.keys()):
                        if str(vv) + str(reference) + "_CORR" in list(att_glo.keys()):
                            att_na = str(vv) + str(reference) + "_CORR"
                        else:
                            att_na = str(vv) + str(reference) + "_" + str(reference) + "_CORR"
                        compare += "; R = " + create_round_string(float(att_glo[att_na]))
                        del att_na
            else:
                compare = None
            if "nina" in vv and "nina_years" in list(d1[vv]["attributes"].keys()):
                nbr_ev = str(len(d1[vv]["attributes"]["nina_years"].split(", "))).zfill(2)
            elif "nino" in vv and "nino_years" in list(d1[vv]["attributes"].keys()):
                nbr_ev = str(len(d1[vv]["attributes"]["nino_years"].split(", "))).zfill(2)
            else:
                nbr_ev = None
            lat = list(tab.getLatitude()[:])
            lon = list(tab.getLongitude()[:])
            # ticks
            xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
            xlabel_ticks, xlabel = create_labels(dict_param[pkey]["xname"], xlabel_ticks)
            ylabel_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
            ylabel_ticks, ylabel = create_labels(dict_param[pkey]["yname"], ylabel_ticks)
            # ax
            ax = plt.subplot(gs[y_pos: y_pos + sizey, jj * (sizex + deltx):  jj * (sizex + deltx) + sizex])
            locmap = Basemap(projection="cyl", llcrnrlat=lat[0], urcrnrlat=lat[-1], llcrnrlon=lon[0], urcrnrlon=lon[-1],
                             ax=ax)
            # draw coastlines
            locmap.drawcoastlines(linewidth=0.5)
            # x-y axes
            ax.set_xticks(xlabel_ticks[1: -1], minor=False)
            ax.set_xticks([kk - (xlabel_ticks[1] - xlabel_ticks[0]) / 2. for kk in xlabel_ticks], minor=True)
            if ii == len(dict_param[pkey]["varpattern"]) - 1:
                ax.set_xticklabels(xlabel[1: -1])
            else:
                ax.set_xticklabels(["" for kk in xlabel[1: -1]])
            ax.set_xlim(min(xlabel_ticks), max(xlabel_ticks))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            ax.set_yticks(ylabel_ticks, minor=False)
            ax.set_yticks([kk - (ylabel_ticks[1] - ylabel_ticks[0]) / 2. for kk in ylabel_ticks], minor=True)
            if jj == 0:
                ax.set_yticklabels(ylabel)
            else:
                ax.set_yticklabels(["" for kk in ylabel])
            ax.set_ylim(-90, 90)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            ax.tick_params(axis="both", direction="in", which="both", bottom=True, top=True, left=True, right=True)
            # title
            x1, x2 = ax.get_xlim()
            y1, y2 = ax.get_ylim()
            dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
            ax.scatter(x1 + 4. * dx, y2 - 10. * dy, s=150, marker="s", c="darkgrey", zorder=2)
            ax.text(x1 + 4. * dx, y2 - 10. * dy, numbering[counter], fontsize=fontsize, weight="bold", zorder=6,
                    ha="center", va="center")
            if ii == 0:
                color = dict_col[project.upper()] if jj == 0 else dict_col["REF"]
                t1 = str(tax.split(", ")[0].split("-")[0].split("'")[1]).zfill(4)
                t2 = str(tax.split(", ")[1].split("-")[0].split("'")[1]).zfill(4)
                if d2 == model and isinstance(member, str) is True:
                    txt = str(d2) + " " + str(member) + "\n" + str(t1) + "-" + str(t2)
                else:
                    txt = str(d2) + "\n" + str(t1) + "-" + str(t2)
                ax.text(0.5, 1.1, txt, fontsize=fontsize, color=color, ha="center", va="bottom", weight="bold",
                        transform=ax.transAxes)
                del color, t1, t2, txt
            if jj == 0:
                txt = "La Nina" if "nina" in vv else "El Nino"
                ax.text(-0.1 * 22. / sizex, 0.5, txt, fontsize=fontsize, color="k", ha="right", va="center",
                        rotation=90, weight="bold", transform=ax.transAxes)
                del txt
            # comparison numbers
            if isinstance(compare, str) is True:
                ax.text(x1 + 1.1 * dx, y1 + 2 * dy, compare, fontsize=6, zorder=6, ha="left", va="bottom",
                        bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
            # number of events
            if isinstance(nbr_ev, str) is True:
                ax.text(x2 - 1.1 * dx, y1 + 2 * dy, nbr_ev, fontsize=6, zorder=6, ha="right", va="bottom",
                        bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
            # shading
            xx, yy = NUMPYmeshgrid(lon, lat)
            levels = create_levels(labelbar)
            cs = locmap.contourf(xx, yy, tab, levels=levels, extend="both", cmap=colorbar)
            # plot regions
            for reg in list_regions:
                my_reg = ReferenceRegions(reg)
                if "polygon" in list(my_reg.keys()) and my_reg["polygon"] is True:
                    lats, lons = my_reg["latitude"], my_reg["longitude"]
                else:
                    lats = list(my_reg["latitude"]) + list(reversed(list(my_reg["latitude"])))
                    lons = [list(my_reg["longitude"])[0]] * 2 + [list(my_reg["longitude"])[1]] * 2
                xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
                ax.add_patch(Polygon(xy, linewidth=1, edgecolor="grey", facecolor="none", zorder=4))
                if max(lons) > 360:
                    lons = list(NUMPYarray(lons) - 360.)
                    xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
                    ax.add_patch(Polygon(xy, linewidth=1, edgecolor="grey", facecolor="none", zorder=4))
                del lats, lons, my_reg, xy
            counter += 1
            if ii == 0 and jj == 1:
                x1_fig = ax.get_position().x0
                x2_fig = ax.get_position().x1
                y2_fig = ax.get_position().y1
            if ii == len(dict_param[pkey]["varpattern"]) - 1 and jj == 0:
                x1_fig3 = ax.get_position().x0
            if ii == len(dict_param[pkey]["varpattern"]) - 1 and jj == 1:
                y1_fig = ax.get_position().y0
            # delete
            del dx, dy, lat, levels, locmap, lon, tab, tax, x1, x2, xlabel_ticks, xlabel, xx, y1, y2, ylabel_ticks, \
                ylabel, yy
        y_pos += sizey + delty
    cax = plt.axes([x2_fig + (x2_fig - x1_fig) / 25., y1_fig, (x2_fig - x1_fig) / 25., y2_fig - y1_fig])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(dict_param[pkey]["zname"] + " (" + str(d1[vv]["attributes"]["units"]) + ")", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    # method and note
    lines = list()
    txt = dict_param[pkey]["method"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    lines += ["Numbers at the bottom right of the maps: number of events"]
    txt = dict_param[pkey]["note"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x1_fig3 - (x2_fig - x1_fig) / 3.5, y1_fig - (y2_fig - y1_fig) / 8., txt, fontsize=fontsize, color="k",
            ha="left", va="top", transform=fig.transFigure)
    y_pos += 2
    # delete
    del ax, cax, cbar, lines, txt
    # ---------------------------------------------------#
    # dotplot
    # ---------------------------------------------------#
    pkey = "02_plot"
    sizex, sizey = 22, 4
    deltx, delty = 1, 1
    mini, maxi = minimaxi([NUMPYpercentile(d1[vv]["array"], [25, 75], axis=0) for vv in dict_param[pkey]["varpattern"]
                           for d1 in [dict_mod, dict_obs]])
    maxi = max([abs(mini), maxi])
    tmp = MATHceil(maxi * 2)
    if tmp % 2 == 0 and tmp > 3:
        dy = int(round(0.25 * tmp, 0))
        y_ticks = list(range(-dy, dy + 1, dy))
        y_label = [str(ii) for ii in y_ticks]
    else:
        dy = 0.25 * tmp
        y_ticks = [round(ii, 2) for ii in NUMPYarange(-dy, dy + 0.1, dy)]
        y_label = list()
        for ii in y_ticks:
            if len(str(dy).split(".")[1]) == 2:
                y_label.append("{0:.2f}".format(round(ii, 2)))
            else:
                y_label.append("{0:.1f}".format(round(ii, 1)))
    for ii, vv in enumerate(dict_param[pkey]["varpattern"]):
        # ax
        ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
        # x-y axes
        ax.set_xticks(list(range(len(list_regions))), minor=False)
        ax.set_xticklabels(["" for kk in list_regions])
        ax.set_xlim(-0.5, len(list_regions) - 0.5)
        ax.set_yticks(y_ticks, minor=False)
        ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
        ax.set_yticklabels(y_label)
        ax.set_ylim(-maxi, maxi)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        ax.tick_params(axis="both", direction="in", which="both", bottom=True, top=True, left=True, right=True)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
        i2, i3 = [y1 for kk in i1], [y2 for kk in i1]
        ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
        ax.hlines([0], [x1], [x2], color="grey", linestyle="-", lw=1, zorder=1)
        # title
        dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
        ax.scatter(x1 + 2. * dx, y2 - 12. * dy, s=150, marker="s", c="darkgrey", zorder=2)
        ax.text(x1 + 2. * dx, y2 - 12. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
                va="center")
        txt = "La Nina" if "nina" in vv else "El Nino"
        ax.text(-0.1 * 22. / sizex, 0.5, txt, fontsize=fontsize, color="k", ha="right", va="center", weight="bold",
                rotation=90, transform=ax.transAxes)
        for jj, d1 in enumerate([dict_mod, dict_obs]):
            # array and lat-lon
            tab = [list(d1[vv]["array"][:, kk]) for kk in range(len(d1[vv]["array"][0]))]
            # markers
            inds = [round(kk, 2) for kk in NUMPYarange(-0.25 + jj * 0.5, len(tab) - 1 + jj * 0.5, 1)]
            color = dict_col[project.upper()] if jj == 0 else dict_col["REF"]
            quartile1, quartile3 = NUMPYpercentile(tab, [25, 75], axis=1)
            mean = NUMPYmean(tab, axis=1)
            if isinstance(dict_mod_extra, dict) is True and isinstance(dict_obs_extra, dict) is True and \
                    "varpattern_extra" in list(dict_param[pkey].keys()) and \
                    dict_param[pkey]["varpattern_extra"][ii] in dict_mod_extra.keys() and \
                    dict_param[pkey]["varpattern_extra"][ii] in dict_obs_extra.keys():
                if jj == 0:
                    tab_sig = dict_mod_extra[dict_param[pkey]["varpattern_extra"][ii]]["array"]
                else:
                    tab_sig = dict_obs_extra[dict_param[pkey]["varpattern_extra"][ii]]["array"]
            else:
                tab_sig = [1] * len(mean)
            for i1, i2, i3 in zip(inds, mean, tab_sig):
                if i3 == 1:
                    ax.plot([i1], [i2], markersize=4, color=color, marker="D", fillstyle="full", markeredgecolor=color,
                            markeredgewidth=1, zorder=3)
                else:
                    ax.plot([i1], [i2], markersize=4, color="none", marker="D", fillstyle="none", markeredgecolor=color,
                            markeredgewidth=1, zorder=3)
            ax.vlines(inds, quartile1, quartile3, color=color, linestyle="-", lw=1, zorder=2)
        counter += 1
        if ii == 0:
            x1_fig = ax.get_position().x0
            x2_fig = ax.get_position().x1
            y2_fig = ax.get_position().y1
            ax.plot([x1], [y1 - 15. * dy], markersize=5, color="k", marker="D", fillstyle="full",
                    markeredgecolor="k", markeredgewidth=1, zorder=3, clip_on=False)
            ax.text(x1 + 2. * dx, y1 - 15. * dy, "significant anomalies", fontsize=fontsize, ha="left", va="center")
            ax.plot([x1 + 50. * dx], [y1 - 15. * dy], markersize=5, color="none", marker="D", fillstyle="none",
                    markeredgecolor="k", markeredgewidth=1, zorder=3, clip_on=False)
            ax.text(x1 + 52. * dx, y1 - 15. * dy, "not significant", fontsize=fontsize, ha="left",
                    va="center")
        if ii == len(dict_param[pkey]["varpattern"]) - 1:
            y1_fig = ax.get_position().y0
            for kk, reg in enumerate(list_regions):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, ha="right", va="top", rotation=45)
        y_pos += sizey + delty
    # zname
    txt = dict_param[pkey]["zname"] + "\n(" + str(d1[vv]["attributes"]["units"]) + ")"
    ax.text(x2_fig + (x2_fig - x1_fig) / 15., y1_fig + (y2_fig - y1_fig) / 2., txt, fontsize=fontsize, color="k",
            ha="center", va="center", rotation=90, transform=fig.transFigure)
    # method and note
    lines = list()
    txt = dict_param[pkey]["method"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    lines += []
    txt = dict_param[pkey]["note"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x1_fig - (x2_fig - x1_fig) / 3.5, y1_fig - (y2_fig - y1_fig) / 8., txt, fontsize=fontsize, color="k",
            ha="left", va="top", transform=fig.transFigure)
    y_pos += 2
    # delete
    del ax, lines, txt
    # ---------------------------------------------------#
    # dotplot
    # ---------------------------------------------------#
    pkey = "03_plot"
    sizex, sizey = 22, 4
    deltx, delty = 1, 1
    mini, maxi = minimaxi([d1[vv]["array"] for vv in dict_param[pkey]["varpattern"] for d1 in [dict_mod, dict_obs]])
    maxi = max([abs(mini), maxi])
    tmp = MATHceil(maxi * 2)
    if tmp % 2 == 0 and tmp > 3:
        dy = int(round(0.25 * tmp, 0))
        y_ticks = list(range(-dy, dy + 1, dy))
        y_label = [str(ii) for ii in y_ticks]
    else:
        dy = 0.25 * tmp
        y_ticks = [round(ii, 2) for ii in NUMPYarange(-dy, dy + 0.1, dy)]
        y_label = list()
        for ii in y_ticks:
            if len(str(dy).split(".")[1]) == 2:
                y_label.append("{0:.2f}".format(round(ii, 2)))
            else:
                y_label.append("{0:.1f}".format(round(ii, 1)))
    # ax
    ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
    # x-y axes
    ax.set_xticks(list(range(len(list_regions))), minor=False)
    ax.set_xticklabels(["" for kk in list_regions])
    ax.set_xlim(-0.5, len(list_regions) - 0.5)
    ax.set_yticks(y_ticks, minor=False)
    ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
    ax.set_yticklabels(y_label)
    ax.set_ylim(-maxi, maxi)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    ax.tick_params(axis="both", direction="in", which="both", bottom=True, top=True, left=True, right=True)
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
    i2, i3 = [y1 for kk in i1], [y2 for kk in i1]
    ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
    ax.hlines([0], [x1], [x2], color="grey", linestyle="-", lw=1, zorder=1)
    # title
    dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
    ax.scatter(x1 + 2. * dx, y2 - 12. * dy, s=150, marker="s", c="darkgrey", zorder=2)
    ax.text(x1 + 2. * dx, y2 - 12. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
            va="center")
    ax.text(-0.1 * 22. / sizex, 0.5, "both", fontsize=fontsize, color="k", ha="right", va="center", weight="bold",
            rotation=90, transform=ax.transAxes)
    for ii, vv in enumerate(dict_param[pkey]["varpattern"]):
        # arrays
        tab_mod = list(dict_mod[vv]["array"])
        tab_obs = list(dict_obs[vv]["array"])
        if isinstance(dict_mod_extra, dict) is True and isinstance(dict_obs_extra, dict) is True and \
                "varpattern_extra" in list(dict_param[pkey].keys()) and \
                dict_param[pkey]["varpattern_extra"][ii] in dict_mod_extra.keys() and \
                dict_param[pkey]["varpattern_extra"][ii] in dict_obs_extra.keys():
            tab_sig = list(dict_obs_extra[dict_param[pkey]["varpattern_extra"][ii]]["array"])
        else:
            tab_sig = [1] * len(tab_mod)
        inds = [round(kk, 2) for kk in NUMPYarange(-0.25 + 0.5 * ii, len(tab_mod) - 1 + 0.5 * ii, 1)]
        for i1, i2, i3, i4 in zip(inds, tab_mod, tab_obs, tab_sig):
            if i4 == 1:
                ax.plot([i1], [i2], markersize=4, color=dict_col[project.upper()], marker="D", fillstyle="full",
                        markeredgecolor=dict_col[project.upper()], markeredgewidth=1, zorder=3)
                ax.plot([i1], [i3], markersize=4, color=dict_col["REF"], marker="D", fillstyle="full",
                        markeredgecolor=dict_col["REF"], markeredgewidth=1, zorder=3)
                if (i2 <= 0 and i3 <= 0) or (i2 >= 0 and i3 >= 0):
                    ax.vlines(i1, y1, y2, color="springgreen", linestyle="-", lw=5, zorder=1)
                else:
                    ax.vlines(i1, y1, y2, color="lightcoral", linestyle="-", lw=5, zorder=1)
            if i1 == inds[0]:
                txt = "LN" if "nina" in vv else "EN"
                ax.text(i1, y2 + 13 * dy, txt, fontsize=fontsize / 1.3, color="k", ha="center", va="center",
                        rotation=90)
        counter += 1
        if ii == 0:
            x1_fig = ax.get_position().x0
            x2_fig = ax.get_position().x1
            y2_fig = ax.get_position().y1
        if ii == len(dict_param[pkey]["varpattern"]) - 1:
            y1_fig = ax.get_position().y0
            for kk, reg in enumerate(list_regions):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, ha="right", va="top", rotation=45)
        y_pos += sizey + delty
    # zname
    txt = dict_param[pkey]["zname"] + "\n(" + str(d1[vv]["attributes"]["units"]) + ")"
    ax.text(x2_fig + (x2_fig - x1_fig) / 15., y1_fig + (y2_fig - y1_fig) / 2., txt, fontsize=fontsize, color="k",
            ha="center", va="center", rotation=90, transform=fig.transFigure)
    # method and note
    lines = list()
    txt = dict_param[pkey]["method"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=70)
    lines += []
    txt = dict_param[pkey]["note"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=70)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x1_fig - (x2_fig - x1_fig) / 10., y1_fig - (y2_fig - y1_fig) / 2., txt, fontsize=fontsize, color="k",
            ha="left", va="top", transform=fig.transFigure)
    y_pos += 2
    plt.savefig(figure_name, bbox_inches="tight")
    plt.savefig(str(figure_name) + ".eps", bbox_inches="tight", format="eps")
    plt.close()
    return


def plot_telecon_1mem(model, project, experiment, nc_file, dict_param, reference, figure_name, dict_diagnostic_values,
                      diagnostic_units, dict_metric_values, metric_units, metric_variables=None, metric_regions=None,
                      member=None):
    # get data
    nc_var = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk:
            nc_var += dict_param[kk]["varpattern"]
    nc_var_extra = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk:
            if "varpattern_extra" in list(dict_param[kk].keys()):
                nc_var_extra += dict_param[kk]["varpattern_extra"]
    dict_mod, dict_obs, dia_mod, dia_obs, met_val, att_glo, dict_mod_extra, dict_obs_extra = read_data(
        nc_file, nc_var, model, reference, dict_diagnostic_values, dict_metric_values, member=member,
        netcdf_var_extra=nc_var_extra)
    met_val = create_round_string(met_val)
    # get regions
    list_regions = metric_regions[metric_variables[-1]]
    # initialization of the plot
    nbrx, nbry = 32, 32
    fig = plt.figure(0, figsize=(nbrx / 4, nbry / 4))
    gs = GridSpec(nbry, nbrx, figure=fig)
    fontsize = 10
    # ---------------------------------------------------#
    # maps
    # ---------------------------------------------------#
    pkey = "01_plot"
    colorbar = dict_param[pkey]["colorbar"]
    labelbar = dict_param[pkey]["label"]
    sizex, sizey = 14, 7
    deltx, delty = 1, 1
    counter = 0
    y_pos = 0
    for ii, vv in enumerate(dict_param[pkey]["varpattern"]):
        for jj, (d1, d2) in enumerate(zip([dict_mod, dict_obs], [model, reference])):
            # array and attributes
            tab = d1[vv]["array"]
            tax = d1[vv]["attributes"]["time_period"]
            if "arraySTD" in list(d1[vv]["attributes"].keys()):
                compare = "$\sigma$ = " + create_round_string(float(d1[vv]["attributes"]["arraySTD"]))
                if d2 == model:
                    if str(vv) + str(reference) + "_RMSE" in list(att_glo.keys()) or \
                            str(vv) + str(reference) + "_" + str(reference) + "_RMSE" in list(att_glo.keys()):
                        if str(vv) + str(reference) + "_RMSE" in list(att_glo.keys()):
                            att_na = str(vv) + str(reference) + "_RMSE"
                        else:
                            att_na = str(vv) + str(reference) + "_" + str(reference) + "_RMSE"
                        compare += "; rmse = " + create_round_string(float(att_glo[att_na]))
                        del att_na
                    if str(vv) + str(reference) + "_CORR" in list(att_glo.keys()) or \
                            str(vv) + str(reference) + "_" + str(reference) + "_CORR" in list(att_glo.keys()):
                        if str(vv) + str(reference) + "_CORR" in list(att_glo.keys()):
                            att_na = str(vv) + str(reference) + "_CORR"
                        else:
                            att_na = str(vv) + str(reference) + "_" + str(reference) + "_CORR"
                        compare += "; R = " + create_round_string(float(att_glo[att_na]))
                        del att_na
            else:
                compare = None
            if "nina" in vv and "nina_years" in list(d1[vv]["attributes"].keys()):
                nbr_ev = str(len(d1[vv]["attributes"]["nina_years"].split(", "))).zfill(2)
            elif "nino" in vv and "nino_years" in list(d1[vv]["attributes"].keys()):
                nbr_ev = str(len(d1[vv]["attributes"]["nino_years"].split(", "))).zfill(2)
            else:
                nbr_ev = None
            lat = list(tab.getLatitude()[:])
            lon = list(tab.getLongitude()[:])
            # ticks
            xlabel_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
            xlabel_ticks, xlabel = create_labels(dict_param[pkey]["xname"], xlabel_ticks)
            ylabel_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
            ylabel_ticks, ylabel = create_labels(dict_param[pkey]["yname"], ylabel_ticks)
            # ax
            ax = plt.subplot(gs[y_pos: y_pos + sizey, jj * (sizex + deltx):  jj * (sizex + deltx) + sizex])
            locmap = Basemap(projection="cyl", llcrnrlat=lat[0], urcrnrlat=lat[-1], llcrnrlon=lon[0], urcrnrlon=lon[-1],
                             ax=ax)
            # draw coastlines
            locmap.drawcoastlines(linewidth=0.5)
            # x-y axes
            ax.set_xticks(xlabel_ticks[1: -1], minor=False)
            ax.set_xticks([kk - (xlabel_ticks[1] - xlabel_ticks[0]) / 2. for kk in xlabel_ticks], minor=True)
            if ii == len(dict_param[pkey]["varpattern"]) - 1:
                ax.set_xticklabels(xlabel[1: -1])
            else:
                ax.set_xticklabels(["" for kk in xlabel[1: -1]])
            ax.set_xlim(min(xlabel_ticks), max(xlabel_ticks))
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            ax.set_yticks(ylabel_ticks, minor=False)
            ax.set_yticks([kk - (ylabel_ticks[1] - ylabel_ticks[0]) / 2. for kk in ylabel_ticks], minor=True)
            if jj == 0:
                ax.set_yticklabels(ylabel)
            else:
                ax.set_yticklabels(["" for kk in ylabel])
            ax.set_ylim(-90, 90)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            ax.tick_params(axis="both", direction="in", which="both", bottom=True, top=True, left=True, right=True)
            # title
            x1, x2 = ax.get_xlim()
            y1, y2 = ax.get_ylim()
            dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
            ax.scatter(x1 + 3. * dx, y2 - 7. * dy, s=150, marker="s", c="darkgrey", zorder=2)
            ax.text(x1 + 3. * dx, y2 - 7. * dy, numbering[counter], fontsize=fontsize, weight="bold", zorder=6,
                    ha="center", va="center")
            if ii == 0:
                color = dict_col[project.upper()] if jj == 0 else dict_col["REF"]
                t1 = str(tax.split(", ")[0].split("-")[0].split("'")[1]).zfill(4)
                t2 = str(tax.split(", ")[1].split("-")[0].split("'")[1]).zfill(4)
                if d2 == model and isinstance(member, str) is True:
                    txt = str(d2) + " " + str(member) + "\n" + str(t1) + "-" + str(t2)
                else:
                    txt = str(d2) + "\n" + str(t1) + "-" + str(t2)
                ax.text(0.5, 1.1, txt, fontsize=fontsize, color=color, ha="center", va="bottom", weight="bold",
                        transform=ax.transAxes)
                del color, t1, t2, txt
            # comparison numbers
            if isinstance(compare, str) is True:
                ax.text(x1 + 1.1 * dx, y1 + 2 * dy, compare, fontsize=6, zorder=6, ha="left", va="bottom",
                        bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
            # number of events
            if isinstance(nbr_ev, str) is True:
                ax.text(x2 - 1.1 * dx, y1 + 2 * dy, nbr_ev, fontsize=6, zorder=6, ha="right", va="bottom",
                        bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
            # shading
            xx, yy = NUMPYmeshgrid(lon, lat)
            levels = create_levels(labelbar)
            cs = locmap.contourf(xx, yy, tab, levels=levels, extend="both", cmap=colorbar)
            # plot regions
            for reg in list_regions:
                my_reg = ReferenceRegions(reg)
                if "polygon" in list(my_reg.keys()) and my_reg["polygon"] is True:
                    lats, lons = my_reg["latitude"], my_reg["longitude"]
                else:
                    lats = list(my_reg["latitude"]) + list(reversed(list(my_reg["latitude"])))
                    lons = [list(my_reg["longitude"])[0]] * 2 + [list(my_reg["longitude"])[1]] * 2
                xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
                ax.add_patch(Polygon(xy, linewidth=1, edgecolor="grey", facecolor="none", zorder=1))
                if min(lons) < 0 or max(lons) > 360:
                    if min(lons) < 0:
                        lons = list(NUMPYarray(lons) + 360.)
                    else:
                        lons = list(NUMPYarray(lons) - 360.)
                    xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
                    ax.add_patch(Polygon(xy, linewidth=1, edgecolor="grey", facecolor="none", zorder=1))
                del lats, lons, my_reg, xy
            counter += 1
            if jj == 0:
                x11 = ax.get_position().x0
                y11 = ax.get_position().y0
                y22 = ax.get_position().y1
            if jj == 1:
                x22 = ax.get_position().x1
            if ii == 0 and jj == 1:
                x1_fig = ax.get_position().x0
                x2_fig = ax.get_position().x1
                y2_fig = ax.get_position().y1
            if ii == len(dict_param[pkey]["varpattern"]) - 1 and jj == 1:
                y1_fig = ax.get_position().y0
            # delete
            del dx, dy, lat, levels, locmap, lon, tab, tax, x1, x2, xlabel_ticks, xlabel, xx, y1, y2, ylabel_ticks, \
                ylabel, yy
        txt = "La Nina" if "nina" in vv else "El Nino"
        ax.text(x11 - 2.1 * (x22 - x11) / (2 * sizex + deltx), y11 + (y22 - y11) / 2., txt, fontsize=fontsize,
                color="k", ha="right", va="center", rotation=90, weight="bold", transform=fig.transFigure)
        del txt
        y_pos += sizey + delty
    cax = plt.axes([x2_fig + (x2_fig - x1_fig) / 25., y1_fig, (x2_fig - x1_fig) / 25., y2_fig - y1_fig])
    cbar = plt.colorbar(cs, cax=cax, orientation="vertical", ticks=labelbar, pad=0.35, extend="both")
    cbar.set_label(dict_param[pkey]["zname"] + " (" + str(d1[vv]["attributes"]["units"]) + ")", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    # method and note
    lines = list()
    txt = dict_param[pkey]["method"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    lines += ["Numbers at the bottom right of the maps: number of events"]
    txt = dict_param[pkey]["note"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=80)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x11 - 2.1 * (x22 - x11) / (2 * sizex + deltx), y11 - 1.3 * (y22 - y11) / sizey, txt, fontsize=fontsize,
            color="k", ha="left", va="top", transform=fig.transFigure)
    y_pos += 3
    # delete
    del ax, cax, cbar, lines, txt
    # ---------------------------------------------------#
    # dotplot
    # ---------------------------------------------------#
    pkey = "02_plot"
    dict_nbr_ev = dict()
    for v1 in dict_param[pkey]["varpattern"]:
        for jj, (d1, d2) in enumerate(zip([dict_mod, dict_obs], [model, reference])):
            if "nina" in v1 and "nina_years" in list(d1[v1]["attributes"].keys()):
                nbr_ev = len(d1[v1]["attributes"]["nina_years"].split(", "))
            elif "nino" in v1 and "nino_years" in list(d1[v1]["attributes"].keys()):
                nbr_ev = len(d1[v1]["attributes"]["nino_years"].split(", "))
            else:
                nbr_ev = None
            if v1 in list(dict_nbr_ev.keys()):
                dict_nbr_ev[v1][d2] = nbr_ev
            else:
                dict_nbr_ev[v1] = {d2: nbr_ev}
            del nbr_ev
    if isinstance(dict_mod_extra, dict) is True and isinstance(dict_obs_extra, dict) is True and \
            "varpattern_extra" in list(dict_param[pkey].keys()):
        tmp1 = [
            minimaxi(my_mask(dict_obs[v1]["array"], dict_obs_extra[v2]["array"]))
            for v1, v2 in zip(dict_param[pkey]["varpattern"], dict_param[pkey]["varpattern_extra"][-2:])]
        tmp2 = [minimaxi([my_mask(dict_mod_extra[v1]["array"][0], dict_obs_extra[v2]["array"]),
                          my_mask(dict_mod_extra[v1]["array"][1], dict_obs_extra[v2]["array"])])
                for v1, v2 in zip(dict_param[pkey]["varpattern_extra"][:2], dict_param[pkey]["varpattern_extra"][-2:])]
        list_reg1 = [reg for ii, reg in enumerate(list_regions)
                     if dict_obs_extra[dict_param[pkey]["varpattern_extra"][-2]]["array"][ii] == 1] + \
                    [reg for ii, reg in enumerate(list_regions)
                     if dict_obs_extra[dict_param[pkey]["varpattern_extra"][-1]]["array"][ii] == 1]
        list_reg1 = [reg for reg in list_regions if reg in list_reg1]
    else:
        tmp1 = [minimaxi(dict_obs[v1]["array"]) for v1 in dict_param[pkey]["varpattern"]]
        tmp2 = [minimaxi(dict_mod[v1]["array"]) for v1 in dict_param[pkey]["varpattern_extra"][:2]]
        list_reg1 = [reg for reg in list_regions]
    sizex, sizey = len(list_reg1) if len(list_reg1) < nbrx else len(list_reg1) / 2, 5
    deltx, delty = 1, 1
    mini, maxi = minimaxi([tmp1, tmp2])
    maxi = max([abs(mini), maxi])
    tmp = MATHceil(maxi * 2)
    if (tmp % 2 == 0 and tmp > 3) or tmp > 20:
        if 20 < tmp < 400:
            tmp = MATHceil(tmp / 10.)
            tmp += tmp % 4
            tmp = tmp * 10
        elif tmp > 400:
            tmp = MATHceil(tmp / 100.)
            tmp += tmp % 4
            tmp = tmp * 100
        dy = int(round(0.25 * tmp, 0))
        y_ticks = list(range(-dy, dy + 1, dy))
        y_label = [str(ii) for ii in y_ticks]
    else:
        dy = 0.25 * tmp
        y_ticks = [round(ii, 2) for ii in NUMPYarange(-dy, dy + 0.1, dy)]
        y_label = list()
        for ii in y_ticks:
            if len(str(dy).split(".")[1]) == 2:
                y_label.append("{0:.2f}".format(round(ii, 2)))
            else:
                y_label.append("{0:.1f}".format(round(ii, 1)))
    for ii, (v1, v2, v3, v4) in enumerate(
            zip(dict_param[pkey]["varpattern"], dict_param[pkey]["varpattern_extra"][:2],
                dict_param[pkey]["varpattern_extra"][2: 4], dict_param[pkey]["varpattern_extra"][4:])):
        # ax
        ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
        # x-y axes
        ax.set_xticks(list(range(len(list_reg1))), minor=False)
        ax.set_xticklabels(["" for kk in list_reg1])
        ax.set_xlim(-0.5, len(list_reg1) - 0.5)
        ax.set_yticks(y_ticks, minor=False)
        ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
        ax.set_yticklabels(y_label)
        ax.set_ylim(-1.1 * maxi, 1.1 * maxi)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        ax.tick_params(axis="both", direction="in", which="both", bottom=False, top=False, left=True, right=True)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
        i2, i3 = [y1 for kk in i1], [y2 for kk in i1]
        ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
        ax.hlines([0], [x1], [x2], color="grey", linestyle="-", lw=1, zorder=1)
        # title
        dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
        ax.scatter(x1 + 1.8 * dx, y2 - 10. * dy, s=150, marker="s", c="darkgrey", zorder=2)
        ax.text(x1 + 1.8 * dx, y2 - 10. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
                va="center")
        for jj, (d1, d2) in enumerate(zip([dict_mod, dict_obs], [model, reference])):
            # array and lat-lon
            tab1 = d1[v1]["array"]
            tab4 = dict_obs_extra[v4]["array"]
            # markers and range
            inds = [reg for kk, reg in enumerate(list_regions) if tab4[kk] == 1]
            color = dict_col[project.upper()] if jj == 0 else dict_col["REF"]
            nbr = 0
            if d2 == model:
                tab2 = dict_mod_extra[v2]["array"].reorder("10")
                tab3 = dict_obs_extra[v3]["array"]
                for kk, (i1, i4, i2, i3) in enumerate(zip(tab1, tab4, tab2, tab3)):
                    if i4 == 1:
                        i5 = list_reg1.index(inds[nbr])
                        if i3 == 0:
                            ax.plot([i5], [i1], markersize=6, color=color, marker="D", fillstyle="full",
                                    markeredgecolor=color, markeredgewidth=2, zorder=3)
                        else:
                            ax.plot([i5], [i1], markersize=6, color="white", marker="D", fillstyle="full",
                                    markeredgecolor=color, markeredgewidth=2, zorder=3)
                        ax.plot([i5, i5], [min(i2), max(i2)], color=color, lw=1, ls="-", zorder=2)
                        ax.plot([i5 - 0.3, i5 + 0.3], [min(i2), min(i2)], color=color, lw=1, ls="-", zorder=5)
                        ax.plot([i5 - 0.3, i5 + 0.3], [max(i2), max(i2)], color=color, lw=1, ls="-", zorder=5)
                        nbr += 1
                del tab2, tab3
            else:
                for kk, (i1, i4) in enumerate(zip(tab1, tab4)):
                    if i4 == 1:
                        i5 = list_reg1.index(inds[nbr])
                        ax.plot([i5 - 0.4, i5 + 0.4], [i1, i1], color=color, lw=2, ls="-", zorder=4)
                        nbr += 1
            del color, inds, nbr, tab1, tab4
        x11 = ax.get_position().x0
        x22 = ax.get_position().x1
        y11 = ax.get_position().y0
        y22 = ax.get_position().y1
        txt = "La Nina" if "nina" in v1 else "El Nino"
        ax.text(x11 - 2.1 * (x22 - x11) / sizex, y11 + (y22 - y11) / 2., txt, fontsize=fontsize, color="k",
                ha="right", va="center", rotation=90, weight="bold", transform=fig.transFigure)
        del txt
        counter += 1
        if ii == len(dict_param[pkey]["varpattern"]) - 1:
            y1_fig = ax.get_position().y0
            for kk, reg in enumerate(list_reg1):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, ha="right", va="top", rotation=45)
        if ii == 0:
            y2_fig = ax.get_position().y1
        else:
            x1_fig = ax.get_position().x0
            x2_fig = ax.get_position().x1
            y1_fig = ax.get_position().y0
        y_pos += sizey + delty
    # zname
    txt = dict_param[pkey]["zname"] + "\n(" + str(d1[vv]["attributes"]["units"]) + ")"
    ax.text(x2_fig + (x2_fig - x1_fig) / 15., y1_fig + (y2_fig - y1_fig) / 2., txt, fontsize=fontsize, color="k",
            ha="center", va="center", rotation=90, transform=fig.transFigure)
    # method and note
    lines = list()
    txt = dict_param[pkey]["method"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=90)
    lines += [""]
    txt = dict_param[pkey]["note"].replace("MET_MET", str(att_glo["metric_method"]))
    txt = txt.replace("MET_VAL", str(met_val)).replace("MET_UNI", str(metric_units))
    lines = create_lines(txt, line_o=lines, threshold=90)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x11 - 2.1 * (x22 - x11) / sizex, y11 - 2.5 * (y22 - y11) / sizey, txt, fontsize=fontsize,
            color="k", ha="left", va="top", transform=fig.transFigure)
    y_pos += 2
    # delete
    del ax, lines, txt
    plt.savefig(figure_name, bbox_inches="tight")
    plt.savefig(str(figure_name) + ".eps", bbox_inches="tight", format="eps")
    plt.close()
    return


def plot_telecon_multimem(model, project, experiment, nc_file, dict_param, reference, figure_name,
                          dict_diagnostic_values, diagnostic_units, dict_metric_values, metric_units,
                          metric_variables=None, metric_regions=None, member=None):
    # get data
    nc_var = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk and kk != "01_plot":  # do not read the variables linked to the first plot (map)
            nc_var += dict_param[kk]["varpattern"]
    nc_var_extra = list()
    for kk in list(dict_param.keys()):
        if "_plot" in kk and kk != "01_plot":
            if "varpattern_extra" in list(dict_param[kk].keys()):
                nc_var_extra += dict_param[kk]["varpattern_extra"]
    dict_mod, dict_obs, att_glo, dict_mod_extra, dict_obs_extra = dict(), dict(), dict(), dict(), dict()
    for mem in member:
        v1, dict_obs, _, _, _, v6, v7, v8 = read_data(
            nc_file[str(model) + "_" + str(mem)], nc_var, model, reference, dict_diagnostic_values, dict_metric_values,
            member=mem, netcdf_var_extra=nc_var_extra)
        dict_mod[mem], att_glo[mem], dict_mod_extra[mem], dict_obs_extra[mem] = v1, v6, v7, v8
    # get regions
    list_regions = metric_regions[metric_variables[-1]]
    # find regions with significant anomalies
    pkey = "02_plot"
    dict_reg_per_ev = dict()
    for v1 in dict_param[pkey]["varpattern_extra"][-2:]:
        sum_reg = dict_obs_extra[member[0]][v1]["array"]
        for mem in member[1:]:
            sum_reg += dict_obs_extra[mem][v1]["array"]
        dict_reg_per_ev[v1] = sum_reg
        del sum_reg
    significant_regions = [reg for ii, reg in enumerate(list_regions)
                           if dict_reg_per_ev[dict_param[pkey]["varpattern_extra"][-2]][ii] +
                           dict_reg_per_ev[dict_param[pkey]["varpattern_extra"][-1]][ii] >= 1]
    # initialization of the plot
    nbrx, nbry = 32, 80
    fig = plt.figure(0, figsize=(nbrx / 4, nbry / 4))
    gs = GridSpec(nbry, nbrx, figure=fig)
    fontsize = 10
    # ---------------------------------------------------#
    # maps
    # ---------------------------------------------------#
    y_pos = 0
    counter = 0
    sizex, sizey = 30, 15
    # ticks
    lat = [-90, 90]
    lon = [0, 360]
    x_ticks = list(range(int(MATHfloor(min(lon))), int(MATHceil(max(lon))) + 1))
    x_ticks, x_label = create_labels("longitude", x_ticks)
    y_ticks = list(range(int(MATHfloor(min(lat))), int(MATHceil(max(lat))) + 1))
    y_ticks, y_label = create_labels("latitude", y_ticks)
    # ax
    ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
    locmap = Basemap(projection="cyl", llcrnrlat=lat[0], urcrnrlat=lat[-1], llcrnrlon=lon[0], urcrnrlon=lon[-1], ax=ax)
    # draw coastlines
    locmap.drawcoastlines(linewidth=0.5)
    # x-y axes
    ax.set_xticks(x_ticks[1: -1], minor=False)
    ax.set_xticks([kk - (x_ticks[1] - x_ticks[0]) / 2. for kk in x_ticks], minor=True)
    ax.set_xticklabels(x_label[1: -1])
    ax.set_xlim(min(x_ticks), max(x_ticks))
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    ax.set_yticks(y_ticks, minor=False)
    ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
    ax.set_yticklabels(y_label)
    ax.set_ylim(-90, 90)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    ax.tick_params(axis="both", direction="in", which="both", bottom=True, top=True, left=True, right=True)
    # title
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
    ax.scatter(x1 + 1.5 * dx, y2 - 3. * dy, s=150, marker="s", c="darkgrey", zorder=2)
    ax.text(x1 + 1.5 * dx, y2 - 3. * dy, numbering[counter], fontsize=fontsize, weight="bold", zorder=6,
            ha="center", va="center")
    # datasets
    for ii, (d1, d2) in enumerate(zip([dict_mod[member[0]], dict_obs], [model, reference])):
        color = dict_col[project.upper()] if d2 == model else dict_col["REF"]
        tax = d1[dict_param[pkey]["varpattern"][0]]["attributes"]["time_period"]
        t1 = str(tax.split(", ")[0].split("-")[0].split("'")[1]).zfill(4)
        t2 = str(tax.split(", ")[1].split("-")[0].split("'")[1]).zfill(4)
        if d2 == model:
            txt = str(d2) + " " + str(experiment) + " " + str(len(member)) + " members\n" + str(t1) + "-" + str(t2)
        else:
            txt = str(d2) + "\n" + str(t1) + "-" + str(t2)
        ax.text(0.25 + ii * 0.5, 1.05, txt, fontsize=fontsize, color=color, ha="center", va="bottom", weight="bold",
                transform=ax.transAxes)
        del color, t1, t2, tax, txt
    # draw regions
    for reg in significant_regions:
        reg2 = reg.replace("_L", "").replace("_O", "")
        rot2 = 90 if reg in ["WSAF"] else (45 if reg in ["ESAF"] else 0)
        my_reg = ReferenceRegions(reg)
        if "polygon" in list(my_reg.keys()) and my_reg["polygon"] is True:
            lats, lons = my_reg["latitude"], my_reg["longitude"]
        else:
            lats = list(my_reg["latitude"]) + list(reversed(list(my_reg["latitude"])))
            lons = [list(my_reg["longitude"])[0]] * 2 + [list(my_reg["longitude"])[1]] * 2
        xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
        ax.add_patch(Polygon(xy, linewidth=1, edgecolor="darkorange", facecolor="none", zorder=1))
        if min(lons) < 0 or max(lons) > 360:
            if min(lons) < 0:
                lons = list(NUMPYarray(lons) + 360.)
            else:
                lons = list(NUMPYarray(lons) - 360.)
            xy = NUMPYarray([[l1, l2] for l1, l2 in zip(lons, lats)])
            ax.add_patch(Polygon(xy, linewidth=1, edgecolor="darkorange", facecolor="none", zorder=1))
            i1, i2 = NUMPYmean(lats), NUMPYmean(lons)
        else:
            i1, i2 = NUMPYmean(lats), NUMPYmean(lons)
        ax.text(i2, i1, reg2, fontsize=fontsize, weight="bold", color="sienna", ha="center", va="center", rotation=rot2,
                zorder=2)
        del i1, i2, lats, lons, my_reg, xy
    # method
    x1, x2, y1, y2 = ax.get_position().x0, ax.get_position().x1, ax.get_position().y0, ax.get_position().y1
    lines = list()
    txt = att_glo[member[0]]["metric_method"]
    txt = str(txt).split("the metric is")[0] + "regions are selected if they" + \
        str(txt).split("Monte Carlo resampling) (regions must")[1].replace(")", "")
    lines = create_lines(txt, line_o=lines, threshold=85)
    txt = ""
    for ii, elt in enumerate(lines):
        txt += elt
        if ii != len(lines) - 1:
            txt += "\n"
    ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 - 1.5 * (y2 - y1) / sizey, txt, fontsize=fontsize,
            color="k", ha="left", va="top", transform=fig.transFigure)
    counter += 1
    y_pos += sizey + 6
    del ax, dx, dy, lat, lines, locmap, lon, sizex, sizey, txt, x1, x2, x_label, x_ticks, y1, y2, y_label, y_ticks
    # ---------------------------------------------------#
    # proposed metric number 1
    # ---------------------------------------------------#
    sizex, sizey = len(significant_regions) if len(significant_regions) < nbrx else len(significant_regions) / 2, 5
    delty = 1
    y2_fig = 1
    # range
    tmp, lv1, lv2 = list(), dict_param[pkey]["varpattern"][:2], dict_param[pkey]["varpattern_extra"][-2:]
    for d1, d2 in zip([dict_mod, dict_obs], [model, reference]):
        if d2 == model:
            for mem in member:
                if isinstance(dict_obs_extra, dict) is True and "varpattern_extra" in list(dict_param[pkey].keys()):
                    tmp.append([minimaxi(my_mask(d1[mem][v1]["array"], dict_reg_per_ev[v2]))
                                for v1, v2 in zip(lv1, lv2)])
                else:
                    tmp.append([minimaxi(d1[mem][v1]["array"]) for v1 in lv1])
        else:
            if isinstance(dict_obs_extra, dict) is True and "varpattern_extra" in list(dict_param[pkey].keys()):
                tmp.append([minimaxi(my_mask(d1[v1]["array"], dict_reg_per_ev[v2]))
                            for v1, v2 in zip(lv1, lv2)])
            else:
                tmp.append([minimaxi(d1[v1]["array"]) for v1 in lv1])
    mini, maxi = minimaxi(tmp)
    maxi = max([abs(mini), maxi])
    tmp = MATHceil(maxi * 2)
    if (tmp % 2 == 0 and tmp > 3) or tmp > 20:
        if 20 < tmp < 400:
            tmp = MATHceil(tmp / 10.)
            tmp += tmp % 4
            tmp = tmp * 10
        elif tmp > 400:
            tmp = MATHceil(tmp / 100.)
            tmp += tmp % 4
            tmp = tmp * 100
        dy = int(round(0.25 * tmp, 0))
        y_ticks = list(range(-dy, dy + 1, dy))
        y_label = [str(ii) for ii in y_ticks]
    else:
        dy = 0.25 * tmp
        y_ticks = [round(ii, 2) for ii in NUMPYarange(-dy, dy + 0.1, dy)]
        y_label = list()
        for ii in y_ticks:
            if len(str(dy).split(".")[1]) == 2:
                y_label.append("{0:.2f}".format(round(ii, 2)))
            else:
                y_label.append("{0:.1f}".format(round(ii, 1)))
    # template
    nbr_wrong, nbr_total = 0, 0
    for ii, (v1, v2) in enumerate(zip(lv1, lv2)):
        # ax
        ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
        # x-y axes
        ax.set_xticks(list(range(len(significant_regions))), minor=False)
        ax.set_xticklabels([""] * len(significant_regions))
        ax.set_xlim(-0.5, len(significant_regions) - 0.5)
        ax.set_yticks(y_ticks, minor=False)
        ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
        ax.set_yticklabels(y_label)
        ax.set_ylim(-1.1 * maxi, 1.1 * maxi)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        ax.tick_params(axis="both", direction="in", which="both", bottom=False, top=False, left=True, right=True)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
        i2, i3 = [y1] * len(i1), [y2] * len(i1)
        ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
        ax.hlines([0], [x1], [x2], color="grey", linestyle="-", lw=1, zorder=1)
        # title
        dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
        ax.scatter(x1 + 1.8 * dx, y2 - 10. * dy, s=150, marker="s", c="darkgrey", zorder=2)
        ax.text(x1 + 1.8 * dx, y2 - 10. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
                va="center")
        for jj, (d1, d2) in enumerate(zip([dict_mod, dict_obs], [model, reference])):
            # array and lat-lon
            tab2 = dict_reg_per_ev[v2]
            inds = [reg for kk, reg in enumerate(list_regions) if tab2[kk] >= 1]
            color = dict_col[project.upper()] if d2 == model else dict_col["REF"]
            nbr = 0
            if d2 == model:
                tab1 = [[d1[mem][v1]["array"][kk] for mem in member] for kk, reg in enumerate(list_regions)]
                tab3 = dict_obs[v1]["array"]
                # marker and IQR
                for kk, (i1, i2, i3) in enumerate(zip(tab1, tab2, tab3)):
                    if i2 >= 1:
                        i4 = significant_regions.index(inds[nbr])
                        mean = NUMPYmean(i1)
                        m1, m2 = minimaxi(i1)
                        if NUMPYisnan(mean):
                            ax.text(i4, y1 + (y2 - y1) / 2., "no data", fontsize=fontsize, ha="center", va="center",
                                    color=color, rotation=90, bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
                            nbr_wrong += 1
                        else:
                            if m1 <= i3 <= m2:
                                ax.plot([i4], [mean], markersize=6, color=color, marker="D", fillstyle="full",
                                        markeredgecolor=color, markeredgewidth=2, zorder=3)
                            else:
                                ax.plot([i4], [mean], markersize=6, color="white", marker="D", fillstyle="full",
                                        markeredgecolor=color, markeredgewidth=2, zorder=3)
                                nbr_wrong += 1
                            ax.plot([i4, i4], [m1, m2], color=color, lw=1, ls="-", zorder=2)
                            ax.plot([i4 - 0.3, i4 + 0.3], [m1, m1], color=color, lw=1, ls="-", zorder=5)
                            ax.plot([i4 - 0.3, i4 + 0.3], [m2, m2], color=color, lw=1, ls="-", zorder=5)
                        nbr += 1
                        nbr_total += 1
                        del i4, m1, m2, mean
                del tab3
            else:
                tab1 = d1[v1]["array"]
                # line
                for kk, (i1, i2) in enumerate(zip(tab1, tab2)):
                    if i2 >= 1:
                        i4 = significant_regions.index(inds[nbr])
                        ax.plot([i4 - 0.4, i4 + 0.4], [i1, i1], color=color, lw=2, ls="-", zorder=4)
                        nbr += 1
                        del i4
            del color, inds, nbr, tab1, tab2
        if ii == len(lv1) - 1:
            for kk, reg in enumerate(significant_regions):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, color="k", ha="right", va="top", rotation=45)
        x1 = ax.get_position().x0
        x2 = ax.get_position().x1
        y1 = ax.get_position().y0
        y2 = ax.get_position().y1
        txt = "La Nina" if "nina" in v1 else "El Nino"
        ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 + (y2 - y1) / 2., txt, fontsize=fontsize, color="k",
                ha="right", va="center", rotation=90, weight="bold", transform=fig.transFigure)
        if ii == 0:
            ax.text(0., 1.03, "Q1: Is the obs composite within the range of ensemble's composites?", fontsize=fontsize,
                    color="k", ha="left", va="bottom", weight="bold", transform=ax.transAxes)
            y2_fig = ax.get_position().y1
        else:
            # y name for both panels
            txt = dict_param[pkey]["zname"].replace("composite & bst", "ensemble's\nmean & range") + \
                " (" + str(dict_obs[lv1[0]]["attributes"]["units"]) + ")"
            ax.text(x2 + (x2 - x1) / 15., y1 + (y2_fig - y1) / 2., txt, fontsize=fontsize,
                    color="k",  ha="center", va="center", rotation=90, transform=fig.transFigure)
            # method and note
            txt = "The metric is the percentage of regions in which the observed composite does not fall within " + \
                "ensemble's range of composite value (metric value = " + \
                str(create_round_string(nbr_wrong * 100. / nbr_total)) + "%)"
            lines = create_lines(txt, line_o=[], threshold=85)
            txt = ""
            for kk, elt in enumerate(lines):
                txt += elt
                if kk != len(lines) - 1:
                    txt += "\n"
            ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 - 2.5 * (y2 - y1) / sizey, txt, fontsize=fontsize,
                    color="k", ha="left", va="top", transform=fig.transFigure)
            del lines
        counter += 1
        y_pos += sizey + delty
        del ax, dx, dy, i1, i2, i3, txt, x1, x2, y1, y2
    y_pos += sizey + delty
    del delty, lv1, lv2, maxi, mini, nbr_total, nbr_wrong, sizex, sizey, tmp, y_label, y_ticks, y2_fig
    # ---------------------------------------------------#
    # proposed metric number 2
    # ---------------------------------------------------#
    sizex, sizey = len(significant_regions) if len(significant_regions) < nbrx else len(significant_regions) / 2, 5
    delty = 1
    y2_fig = 1
    lv1, lv2, lv3 = dict_param[pkey]["varpattern"][: 2], dict_param[pkey]["varpattern_extra"][-4: -2], \
        dict_param[pkey]["varpattern_extra"][-2:]
    tr1 = 3. / 4
    tr2 = MATHfloor(len(member) * tr1)
    # range
    maxi = len(member)
    maxi += maxi % 4
    dy = int(round(0.25 * maxi, 0))
    y_ticks = list(range(0, maxi + 1, dy))
    y_label = [str(ii) for ii in y_ticks]
    # template
    nbr_wrong, nbr_total = 0, 0
    for ii, (v1, v2, v3) in enumerate(zip(lv1, lv2, lv3)):
        # ax
        ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
        # x-y axes
        ax.set_xticks(list(range(len(significant_regions))), minor=False)
        ax.set_xticklabels([""] * len(significant_regions))
        ax.set_xlim(-0.5, len(significant_regions) - 0.5)
        ax.set_yticks(y_ticks, minor=False)
        ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
        ax.set_yticklabels(y_label)
        ax.set_ylim(0, maxi)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        ax.tick_params(axis="both", direction="in", which="both", bottom=False, top=False, left=True, right=True)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
        i2, i3 = [y1] * len(i1), [y2] * len(i1)
        ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
        ax.hlines([tr2], [x1], [x2], color="green", linestyle="-", lw=1, zorder=1)
        # title
        dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
        ax.scatter(x1 + 1.8 * dx, y2 - 10. * dy, s=150, marker="s", c="darkgrey", zorder=2)
        ax.text(x1 + 1.8 * dx, y2 - 10. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
                va="center")
        # array
        tab1 = dict_mod[member[0]][v1]["array"]
        tab2 = [[dict_obs_extra[mem][v2]["array"][kk] for mem in member] for kk, reg in enumerate(list_regions)]
        tab3 = dict_reg_per_ev[v3]
        inds = [reg for kk, reg in enumerate(list_regions) if tab3[kk] >= 1]
        # histogram
        color = dict_col[project.upper()]
        nbr = 0
        for kk, (i1, i2, i3) in enumerate(zip(tab1, tab2, tab3)):
            if i3 >= 1:
                i4 = significant_regions.index(inds[nbr])
                sum_mem = sum(i2)
                xy = NUMPYarray([[i4 - 0.3, y1], [i4 - 0.3, sum_mem], [i4 + 0.3, sum_mem], [i4 + 0.3, y1]])
                if NUMPYisnan(i1) or isinstance(i1, float) is False:
                    # no model data
                    ax.text(i4, y1 + (y2 - y1) / 2., "no data", fontsize=fontsize, ha="center", va="center",
                            color=color, rotation=90, bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
                    nbr_wrong += 1
                else:
                    if sum_mem < tr2:
                        ax.add_patch(Polygon(xy, linewidth=2, edgecolor=color, facecolor="none", zorder=4))
                        nbr_wrong += 1
                    else:
                        ax.add_patch(Polygon(xy, linewidth=2, edgecolor=color, facecolor=color, zorder=4))
                nbr += 1
                nbr_total += 1
                del i4, sum_mem, xy
        if ii == len(lv1) - 1:
            for kk, reg in enumerate(significant_regions):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, color="k", ha="right", va="top", rotation=45)
        x1 = ax.get_position().x0
        x2 = ax.get_position().x1
        y1 = ax.get_position().y0
        y2 = ax.get_position().y1
        txt = "La Nina" if "nina" in v1 else "El Nino"
        ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 + (y2 - y1) / 2., txt, fontsize=fontsize, color="k",
                ha="right", va="center", rotation=90, weight="bold", transform=fig.transFigure)
        if ii == 0:
            ax.text(0., 1.03, "Q2: Is the obs composite within the modeled range\n(Monte Carlo resampling) in " +
                    str(int(round(tr1 * 100, 0))) + "% of the members?", fontsize=fontsize, color="k", ha="left",
                    va="bottom", weight="bold", transform=ax.transAxes)
            y2_fig = ax.get_position().y1
        else:
            # y name for both panels
            txt = "number of members"
            ax.text(x2 + (x2 - x1) / 15., y1 + (y2_fig - y1) / 2., txt, fontsize=fontsize,
                    color="k", ha="center", va="center", rotation=90, transform=fig.transFigure)
            # method and note
            txt = att_glo[member[0]]["metric_method"].split("within modeled range ")[1].split(" (regions must have")[0]
            if "(90 of the Monte Carlo resampling)" in txt:
                txt = txt.replace("(90 of the Monte Carlo resampling)", "(90% of the Monte Carlo resampling)")
            txt = "The metric is the percentage of regions in which the observed composite does not fall within the" + \
                " modeled range " + str(txt) + " in at least " + str(int(round(tr1 * 100, 0))) + "% of the members " + \
                "(metric value = " + str(create_round_string(nbr_wrong * 100. / nbr_total)) + "%)"
            lines = create_lines(txt, line_o=[], threshold=85)
            txt = ""
            for kk, elt in enumerate(lines):
                txt += elt
                if kk != len(lines) - 1:
                    txt += "\n"
            ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 - 2.5 * (y2 - y1) / sizey, txt, fontsize=fontsize,
                    color="k", ha="left", va="top", transform=fig.transFigure)
            del lines
        counter += 1
        y_pos += sizey + delty
        del ax, dx, dy, i1, i2, txt, x1, x2, y1, y2
    y_pos += sizey + delty
    del delty, lv1, lv2, maxi, nbr_total, nbr_wrong, sizex, sizey, tr1, tr2, y_label, y_ticks, y2_fig
    # ---------------------------------------------------#
    # proposed metric number 3
    # ---------------------------------------------------#
    sizex, sizey = len(significant_regions) if len(significant_regions) < nbrx else len(significant_regions) / 2, 5
    delty = 1
    y2_fig = 1
    lv1, lv2, lv3 = dict_param[pkey]["varpattern"][:2], dict_param[pkey]["varpattern"][2: 4],\
        dict_param[pkey]["varpattern_extra"][-2:]
    tr1 = 90
    # range
    dict_mem_ave, dict_mem_ran = dict(), dict()
    for v1, v2, v3 in zip(lv1, lv2, lv3):
        if "nina" in v1 and "nina_years" in list(dict_obs[v1]["attributes"].keys()):
            if isinstance(dict_obs[v1]["attributes"]["nina_years"], str) is True:
                nbr_ev = len(dict_obs[v1]["attributes"]["nina_years"].split(", "))
            else:
                nbr_ev = len(dict_obs[v2]["attributes"]["nina_years"])
        elif "nino" in v1 and "nino_years" in list(dict_obs[v1]["attributes"].keys()):
            if isinstance(dict_obs[v1]["attributes"]["nino_years"], str) is True:
                nbr_ev = len(dict_obs[v1]["attributes"]["nino_years"].split(", "))
            else:
                nbr_ev = len(dict_obs[v1]["attributes"]["nino_years"])
        else:
            nbr_ev = None
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": ENSO years not defined",
                str().ljust(5) + "nina_years year or nino_years not defined as attribute in " + str(v1),
                str().ljust(10) + "attributes: " +
                str(sorted(list(dict_obs[v1]["attributes"].keys()), key=lambda v: v.upper()))]
            EnsoErrorsWarnings.my_error(list_strings)
        dict_mem_ave[v1] = my_mask(member_average(dict_mod, v1), dict_reg_per_ev[v3])
        dict_mem_ran[v1] = member_range(dict_mod, v2, nbr_ev, lsig_level=tr1, larray_mask=dict_reg_per_ev[v3])
        del nbr_ev
    tmp = list()
    for d2 in [model, reference]:
        for v1 in lv1:
            if d2 == model:
                tmp.append(minimaxi(dict_mem_ave[v1]))
                tmp.append(minimaxi(dict_mem_ran[v1]))
            else:
                tmp.append(minimaxi(dict_obs[v1]["array"]))
    mini, maxi = minimaxi(tmp)
    maxi = max([abs(mini), maxi])
    tmp = MATHceil(maxi * 2)
    if (tmp % 2 == 0 and tmp > 3) or tmp > 20:
        if 20 < tmp < 400:
            tmp = MATHceil(tmp / 10.)
            tmp += tmp % 4
            tmp = tmp * 10
        elif tmp > 400:
            tmp = MATHceil(tmp / 100.)
            tmp += tmp % 4
            tmp = tmp * 100
        dy = int(round(0.25 * tmp, 0))
        y_ticks = list(range(-dy, dy + 1, dy))
        y_label = [str(ii) for ii in y_ticks]
    else:
        dy = 0.25 * tmp
        y_ticks = [round(ii, 2) for ii in NUMPYarange(-dy, dy + 0.1, dy)]
        y_label = list()
        for ii in y_ticks:
            if len(str(dy).split(".")[1]) == 2:
                y_label.append("{0:.2f}".format(round(ii, 2)))
            else:
                y_label.append("{0:.1f}".format(round(ii, 1)))
    # template
    nbr_wrong, nbr_total = 0, 0
    for ii, (v1, v2, v3) in enumerate(zip(lv1, lv2, lv3)):
        # ax
        ax = plt.subplot(gs[y_pos: y_pos + sizey, 0: sizex])
        # x-y axes
        ax.set_xticks(list(range(len(significant_regions))), minor=False)
        ax.set_xticklabels([""] * len(significant_regions))
        ax.set_xlim(-0.5, len(significant_regions) - 0.5)
        ax.set_yticks(y_ticks, minor=False)
        ax.set_yticks([kk - (y_ticks[1] - y_ticks[0]) / 2. for kk in y_ticks], minor=True)
        ax.set_yticklabels(y_label)
        ax.set_ylim(-1.1 * maxi, 1.1 * maxi)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        ax.tick_params(axis="both", direction="in", which="both", bottom=False, top=False, left=True, right=True)
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        i1 = [round(kk, 2) for kk in NUMPYarange(x1 + 1, x2, 1)]
        i2, i3 = [y1] * len(i1), [y2] * len(i1)
        ax.vlines(i1, i2, i3, color="grey", linestyle="-", lw=1, zorder=1)
        ax.hlines([0], [x1], [x2], color="grey", linestyle="-", lw=1, zorder=1)
        # title
        dx, dy = (x2 - x1) / 100., (y2 - y1) / 100.
        ax.scatter(x1 + 1.8 * dx, y2 - 10. * dy, s=150, marker="s", c="darkgrey", zorder=2)
        ax.text(x1 + 1.8 * dx, y2 - 10. * dy, numbering[counter], fontsize=10, weight="bold", zorder=6, ha="center",
                va="center")
        # array
        tab1 = dict_mem_ave[v1]
        tab2 = dict_mem_ran[v1].reorder("10")
        tab3 = dict_obs[v1]["array"]
        tab4 = dict_reg_per_ev[v3]
        inds = [reg for kk, reg in enumerate(list_regions) if tab4[kk] >= 1]
        co1, co2 = dict_col[project.upper()], dict_col["REF"]
        # marker and range
        nbr = 0
        for kk, (i1, i2, i3, i4) in enumerate(zip(tab1, tab2, tab3, tab4)):
            if i4 >= 1:
                i5 = significant_regions.index(inds[nbr])
                # model
                if NUMPYisnan(i1) or isinstance(i1, float) is False:
                    # no model data
                    ax.text(i5, y1 + (y2 - y1) / 2., "no data", fontsize=fontsize, ha="center", va="center", color=co1,
                            rotation=90, bbox=dict(lw=0, facecolor="white", pad=1, alpha=1))
                    nbr_wrong += 1
                else:
                    if min(i2) <= i3 <= max(i2):
                        # obs within modeled range
                        ax.plot([i5], [i1], markersize=6, color=co1, marker="D", fillstyle="full",
                                markeredgecolor=co1, markeredgewidth=2, zorder=3)
                    else:
                        # obs not within modeled range
                        ax.plot([i5], [i1], markersize=6, color="white", marker="D", fillstyle="full",
                                markeredgecolor=co1, markeredgewidth=2, zorder=3)
                        nbr_wrong += 1
                    # modeled range
                    ax.plot([i5, i5], [min(i2), max(i2)], color=co1, lw=1, ls="-", zorder=2)
                    ax.plot([i5 - 0.3, i5 + 0.3], [min(i2), min(i2)], color=co1, lw=1, ls="-", zorder=5)
                    ax.plot([i5 - 0.3, i5 + 0.3], [max(i2), max(i2)], color=co1, lw=1, ls="-", zorder=5)
                # obs
                ax.plot([i5 - 0.4, i5 + 0.4], [i3, i3], color=co2, lw=2, ls="-", zorder=4)
                nbr += 1
                nbr_total += 1
                del i5
        # write region names on x axis
        if ii == len(lv1) - 1:
            for kk, reg in enumerate(significant_regions):
                ax.text(kk + 2 * dx, y1 - 2 * dy, reg, fontsize=fontsize, color="k", ha="right", va="top", rotation=45)
        # write event type on the left
        x1 = ax.get_position().x0
        x2 = ax.get_position().x1
        y1 = ax.get_position().y0
        y2 = ax.get_position().y1
        txt = "La Nina" if "nina" in v1 else "El Nino"
        ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 + (y2 - y1) / 2., txt, fontsize=fontsize, color="k",
                ha="right", va="center", rotation=90, weight="bold", transform=fig.transFigure)
        if ii == 0:
            # write the question asked
            ax.text(0., 1.03, "Q3: Is the obs composite within the modeled range\n(Monte Carlo resampling) computed " +
                    "with the ensemble?", fontsize=fontsize, color="k", ha="left", va="bottom", weight="bold",
                    transform=ax.transAxes)
            y2_fig = ax.get_position().y1
        else:
            # y name for both panels
            txt = "number of members"
            ax.text(x2 + (x2 - x1) / 15., y1 + (y2_fig - y1) / 2., txt, fontsize=fontsize,
                    color="k", ha="center", va="center", rotation=90, transform=fig.transFigure)
            # method and note
            txt = att_glo[member[0]]["metric_method"].split("within modeled range ")[1].split(" (regions must have")[0]
            if "(90 of the Monte Carlo resampling)" in txt:
                txt = txt.replace("(90 of the Monte Carlo resampling)", "(90% of the Monte Carlo resampling)")
            txt = "The metric is the percentage of regions in which the observed composite does not fall within the" + \
                " modeled range computed with the ensemble (" + str(tr1) + "% of the Monte Carlo resampling) " + \
                "(metric value = " + str(create_round_string(nbr_wrong * 100. / nbr_total)) + "%)"
            lines = create_lines(txt, line_o=[], threshold=85)
            txt = ""
            for kk, elt in enumerate(lines):
                txt += elt
                if kk != len(lines) - 1:
                    txt += "\n"
            ax.text(x1 - 2.1 * (x2 - x1) / sizex, y1 - 2.5 * (y2 - y1) / sizey, txt, fontsize=fontsize,
                    color="k", ha="left", va="top", transform=fig.transFigure)
            del lines
        counter += 1
        y_pos += sizey + delty
        del ax, co1, co2, dx, dy, i1, i2, i3, inds, nbr, tab1, tab2, tab3, tab4, txt, x1, x2, y1, y2
    y_pos += sizey + 1
    del delty, dict_mem_ave, dict_mem_ran, lv1, lv2, lv3, maxi, mini, nbr_total, nbr_wrong, sizex, sizey, tmp, tr1, \
        y_label, y_ticks, y2_fig
    plt.savefig(figure_name, bbox_inches="tight")
    plt.savefig(str(figure_name) + ".eps", bbox_inches="tight", format="eps")
    plt.close()
    return
