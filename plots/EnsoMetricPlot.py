# -*- coding:UTF-8 -*-
from __future__ import print_function
from copy import deepcopy
from datetime import datetime
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoMetrics.EnsoPlotLib import plot_param
from EnsoPlotTemplate import cmip_boxplot, my_boxplot, my_curve, my_dotplot, my_dot_to_box, my_hovmoeller, my_map,\
    my_scatterplot


dict_plot = {"boxplot": my_boxplot, "curve": my_curve, "dot": my_dotplot, "dot_to_box": my_dot_to_box,
             "hovmoeller": my_hovmoeller, "map": my_map, "scatterplot": my_scatterplot}


def cmip_plotter(metric_collection, metric, experiment, diagnostic_values, diagnostic_units, metric_values,
                 metric_units, multi_member_ave, figure_name):
    """
    Organizes plots for CMIP ensembles

    Inputs:
    ------
    :param metric_collection: string
        name of the metric collection (e.g., "ENSO_perf")
    :param metric: string
        name of the metric (e.g., "EnsoAmpl")
    :param diagnostic_values: dictionary
        dictionary containing the diagnostic values
        e.g., {"ref": {"ERA-Interim": 1, "Tropflux": 1.1}, "cmip5": [0.9, ...], "cmip6": [0.9, ...]}
    :param diagnostic_units: string
        diagnostic units (e.g., "C")
    :param metric_values: dictionary
        dictionary containing the metric values
        e.g., {"cmip5": {"ERA-Interim": [10, ...], "Tropflux": [20, ...]},
               "cmip6": {"ERA-Interim": [10, ...], "Tropflux": [20, ...]}}
    :param metric_units: string
        metric units (e.g., "%")
    :param multi_member_ave: boolean
        True if multi-member mean is used, False otherwise
    :param figure_name: string
        path and name of the plots (e.g., "path/to/directory/file_name")

    Output:
    ------
    :return:
    """
    lmet = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "") if "Map" in metric else deepcopy(metric)
    dict_param = plot_param(metric_collection, lmet)
    my_param = dict_param["diagnostic"]
    reference = dict_param["metric_reference"]
    diagnostic_units = diagnostic_units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    metric_units = metric_units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if multi_member_ave is True:
        info = "ensemble-mean\n" + experiment + " simulations"
    else:
        info = "only first\n" + experiment + " simulation"
    # metric
    cmip_boxplot(my_param, metric_values, metric_units, reference, "metric", info, figure_name + "_divedown01")
    # diagnostic
    cmip_boxplot(my_param, diagnostic_values, diagnostic_units, reference, "diagnostic", info,
                 figure_name + "_divedown02")
    return


def main_plotter(metric_collection, metric, model, experiment, filename_nc, diagnostic_values, diagnostic_units,
                 metric_values, metric_units, member=None, path_png=None, name_png=None, models2=None, shading=False,
                 plot_ref=False):
    dict_param = plot_param(metric_collection, metric)
    list_var = dict_param['metric_variables']
    met_type = dict_param['metric_computation']
    reference = dict_param['metric_reference']
    dict_reg = dict_param['metric_regions']
    if isinstance(diagnostic_units, basestring) is True:
        diagnostic_units = diagnostic_units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    elif isinstance(diagnostic_units, list) is True:
        for ii, uni in enumerate(diagnostic_units):
            diagnostic_units[ii] = uni.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if isinstance(metric_units, basestring) is True:
        metric_units = metric_units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    elif isinstance(metric_units, list) is True:
        for ii, uni in enumerate(metric_units):
            metric_units[ii] = uni.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if isinstance(name_png, basestring) is False:
        name_png = metric_collection + "_" + metric + "_" + experiment + "_"
        if isinstance(model, basestring):
            name_png = name_png + model
        elif isinstance(model, list) is True and shading is True:
            name_png = name_png + str(len(model)) + "projects"
        else:
            name_png = name_png + str(len(model)) + "models"
        if member is not None:
            if (isinstance(model, basestring) is True and member not in model) or isinstance(model, list) is True:
                name_png = name_png + "_" + member
    # diagnostic
    dict_diag = dict_param['diagnostic']
    if isinstance(path_png, basestring):
        fig_name = OSpath__join(path_png, name_png + "_diagnostic")
    else:
        fig_name = deepcopy(name_png)
    plt_typ = dict_diag['plot_type']
    if plt_typ == "dot" and shading is True:
        plt_typ = "dot_to_box"
    t1 = datetime.now()
    print(str().ljust(20) + plt_typ + " " + str(t1.hour).zfill(2) + ":" + str(t1.minute).zfill(2))
    dict_plot[plt_typ](
        model, filename_nc, dict_diag, reference, list_var, fig_name + "_divedown01", models2=models2,
        member=member, metric_type=met_type, metric_values=metric_values, metric_units=metric_units,
        diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=shading)
    if plot_ref is True:
        model_new = model.replace("GPCPv2.3", "GPCPv23").replace("SODA3.4.2", "SODA342")
        fig_name_ref = fig_name.replace(model_new, "reference")
        dict_plot[plt_typ](
            model, filename_nc, dict_diag, reference, list_var, fig_name_ref + "_divedown01", models2=None,
            member=None, metric_type=None, metric_values=metric_values, metric_units=metric_units,
            diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=False,
            plot_ref=plot_ref)
    dt = datetime.now() - t1
    dt = str(int(round(dt.seconds / 60.)))
    print(str().ljust(30) + "took " + dt + " minute(s)")
    # dive downs
    list_dd = sorted([key for key in dict_param.keys() if "dive_down" in key], key=lambda v: v.upper())
    for ii, dd in enumerate(list_dd):
        dict_diag = dict_param[dd]
        plt_typ = dict_diag['plot_type']
        t1 = datetime.now()
        print(str().ljust(20) + plt_typ + " " + str(t1.hour).zfill(2) + ":" + str(t1.minute).zfill(2))
        if metric_collection in ["ENSO_tel", "test_tel"] and "Map" in metric:
            metype = deepcopy(met_type)
        else:
            metype = None
        dict_plot[plt_typ](
            model, filename_nc, dict_diag, reference, list_var, fig_name + "_divedown" + str(ii+2).zfill(2),
            models2=models2, member=member, metric_type=metype, metric_values=metric_values, metric_units=metric_units,
            diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=shading)
        if plot_ref is True:
            dict_plot[plt_typ](
                model, filename_nc, dict_diag, reference, list_var, fig_name_ref + "_divedown" + str(ii + 2).zfill(2),
                models2=None, member=None, metric_type=None, metric_values=metric_values, metric_units=metric_units,
                diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=False,
                plot_ref=plot_ref)
        dt = datetime.now() - t1
        dt = str(int(round(dt.seconds / 60.)))
        print(str().ljust(30) + "took " + dt + " minute(s)")

