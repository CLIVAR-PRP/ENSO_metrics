# -*- coding:UTF-8 -*-
from copy import deepcopy
from datetime import datetime
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoPlotLib import plot_param
from EnsoPlotTemplate import my_boxplot, my_curve, my_dotplot, my_dot_to_box, my_hovmoeller, my_map, my_scatterplot


dict_plot = {"boxplot": my_boxplot, "curve": my_curve, "dot": my_dotplot, "dot_to_box": my_dot_to_box,
             "hovmoeller": my_hovmoeller, "map": my_map, "scatterplot": my_scatterplot}


def main_plotter(metric_collection, metric, model, experiment, filename_nc, diagnostic_values,
                 diagnostic_units, metric_values, metric_units, path_png=None, name_png=None, models2=None, shading=False):
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
    print str().ljust(20) + plt_typ + " " + str(t1.hour).zfill(2) + ":" + str(t1.minute).zfill(2)
    dict_plot[plt_typ](
        model, filename_nc, dict_diag, reference, list_var, fig_name + "_divedown01", models2=models2,
        metric_type=met_type, metric_values=metric_values, metric_units=metric_units,
        diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=shading)
    dt = datetime.now() - t1
    dt = str(int(round(dt.seconds / 60.)))
    print str().ljust(30) + "took " + dt + " minute(s)"
    # dive downs
    list_dd = sorted([key for key in dict_param.keys() if "dive_down" in key], key=lambda v: v.upper())
    for ii, dd in enumerate(list_dd):
        dict_diag = dict_param[dd]
        plt_typ = dict_diag['plot_type']
        t1 = datetime.now()
        print str().ljust(20) + plt_typ + " " + str(t1.hour).zfill(2) + ":" + str(t1.minute).zfill(2)
        if metric_collection == "ENSO_tel" and "Map" in metric:
            metype = deepcopy(met_type)
        else:
            metype = None
        dict_plot[plt_typ](
            model, filename_nc, dict_diag, reference, list_var, fig_name + "_divedown" + str(ii+2).zfill(2),
            models2=models2, metric_type=metype, metric_values=metric_values, metric_units=metric_units,
            diagnostic_values=diagnostic_values, diagnostic_units=diagnostic_units, regions=dict_reg, shading=shading)
        dt = datetime.now() - t1
        dt = str(int(round(dt.seconds / 60.)))
        print str().ljust(30) + "took " + dt + " minute(s)"

