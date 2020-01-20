# -*- coding:UTF-8 -*-
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoPlotLib import plot_param
from EnsoPlotTemplate import my_boxplot, my_curve, my_dotplot, my_hovmoeller, my_map, my_scatterplot


dict_plot = {"boxplot": my_boxplot, "curve": my_curve, "dot": my_dotplot, "hovmoeller": my_hovmoeller, "map": my_map,
             "scatterplot": my_scatterplot}


def main_plotter(metric_collection, metric, model, experiment, filename_nc, diagnostic_values,
                 diagnostic_units, metric_values, metric_units, path_png=''):
    dict_param = plot_param(metric_collection, metric)
    list_var = dict_param['metric_variables']
    met_type = dict_param['metric_computation']
    reference = dict_param['metric_reference']
    dict_reg = dict_param['metric_regions']
    if isinstance(model, str):
        pattern = metric_collection + "_" + metric + "_" + experiment + "_" + model
    else:
        pattern = metric_collection + "_" + metric + "_" + experiment + "_" + str(len(model)) + "models"
    # diagnostic
    dict_diag = dict_param['diagnostic']
    fig_name = OSpath__join(path_png, pattern + "_diagnostic")
    print(str().ljust(20) + dict_diag['plot_type'])
    dict_plot[dict_diag['plot_type']](
        model, filename_nc, dict_diag, reference, list_var, fig_name, metric_type=met_type,
        metric_values=metric_values, metric_units=metric_units, diagnostic_values=diagnostic_values,
        diagnostic_units=diagnostic_units, regions=dict_reg)
    # dive downs
    list_dd = sorted([key for key in list(dict_param.keys()) if "dive_down" in key], key=lambda v: v.upper())
    for ii, dd in enumerate(list_dd):
        dict_diag = dict_param[dd]
        fig_name = OSpath__join(path_png, pattern + "_divedown" + str(ii+1).zfill(2))
        print(str().ljust(20) + dict_diag['plot_type'])
        dict_plot[dict_diag['plot_type']](
            model, filename_nc, dict_diag, reference, list_var, fig_name, metric_type=None,
            metric_values=metric_values, metric_units=metric_units, diagnostic_values=diagnostic_values,
            diagnostic_units=diagnostic_units, regions=dict_reg)

