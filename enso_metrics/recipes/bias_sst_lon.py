# -*- coding:UTF-8 -*-
"""
Metric designed to evaluate the temperature bias along the equator.

Old name: BiasSstLonRmse (https://github.com/CLIVAR-PRP/ENSO_metrics/wiki/BiasSstLonRmse)

Diagnostic:
Equatorial Pacific climatological (time and meridional 5°S-5°N average) sea surface temperature (SST)
Metric:
Root mean square error (RMSE) comparing given dataset to given reference
Commanded reference(s):
     - COBE2: https://psl.noaa.gov/data/gridded/data.cobe2.html
     - ERSSTv5: https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
     - HadISST: https://www.metoffice.gov.uk/hadobs/hadisst/

Supplementary diagnostic(s):
     - Tropical Pacific climatological (time average) sea surface temperature (SST)
     - zonal gradient of equatorial Pacific sea surface temperature (SST)
"""
# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
from inspect import stack as inspect__stack
from json import dumps as json__dumps
from os.path import basename as os__path__basename
from typing import Union

# local functions
from enso_metrics.tools.default import default_arg_values, input_dictionary_formater, processors_dictionary_formater, \
    set_default_epoch, set_default_str, set_instance
from enso_metrics.wrapper import basics
from enso_metrics.wrapper import processors as pr
# ---------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Diagnostic
    #
    "epoch": ("1980-01-01", "2014-12-31"),
    "recipe": None,
    "supplementary": False,
    # base region / variable for the metric
    "region1_main": "equatorial_pacific",
    "variable1": "ts",
    # supplementary
    "region1_supp1": "tropical_pacific",
    "region1_supp2": "grad_west",
    "region1_supp3": "grad_east",
    # output see processors.saver for options
    "kwargs_saver": {"extension": "nc"},
    #
    # -- Metric
    #
    "metric_method": "rmse",
}
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Description of the main processing to do
# ---------------------------------------------------------------------------------------------------------------------#
processors_main = {
    "d1_variable1_lon": {
        "region": "region1_main",
        "variable": "variable1",
        "to_do": {
            "1__masker": {"tolerance": 0, "kwargs_where": {}},
            "2__selector": {
                "cf_dims": ["T"],
                "bounds_t": "epoch",
                "bounds_z": None,
                "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}},
                "kwargs_select_t": {"kwargs_sel": {}},
                "kwargs_select_z": {"kwargs_sel": {}}},
            "3__averager": {"cf_dims": ["T"]},
            "4__regridder": {
                "cf_dims": ["XY"],
                "grid_xy": "uniform_1x1",
                "grid_z": None,
                "kwargs_regrid_xy": {"method": "conservative", "tool": "regrid2"},
                "kwargs_regrid_z": {"tool": "xgcm", "kwargs_regridder_vertical": {}}},
            "5__selector": {
                "cf_dims": ["XY"],
                "bounds_t": None,
                "bounds_z": None,
                "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}},
                "kwargs_select_t": {"kwargs_sel": {}},
                "kwargs_select_z": {"kwargs_sel": {}}},
            "6__averager": {"cf_dims": ["Y"]},
        },
    },
}
processors_supp = {
    "d2_variable1_map": {
        "region": "region1_supp1",
        "variable": "variable1",
        "to_do": {
            "1__masker": {"tolerance": 0, "kwargs_where": {}},
            "2__selector": {
                "cf_dims": ["T"],
                "bounds_t": "epoch",
                "bounds_z": None,
                "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}},
                "kwargs_select_t": {"kwargs_sel": {}},
                "kwargs_select_z": {"kwargs_sel": {}}},
            "3__averager": {"cf_dims": ["T"]},
            "4__regridder": {
                "cf_dims": ["XY"],
                "grid_xy": "uniform_1x1",
                "grid_z": None,
                "kwargs_regrid_xy": {"method": "conservative", "tool": "regrid2"},
                # "kwargs_regrid_xy": {"method": "bilinear", "tool": "xesmf"},
                "kwargs_regrid_z": {"tool": "xgcm", "kwargs_regridder_vertical": {}}},
            "5__selector": {
                "cf_dims": ["XY"],
                "bounds_t": None,
                "bounds_z": None,
                "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}},
                "kwargs_select_t": {"kwargs_sel": {}},
                "kwargs_select_z": {"kwargs_sel": {}}},
        },
    },
    "d3_variable1_grad": {
        "region": "region1_supp2",
        "variable": "variable1",
        "to_do": {
            "1__masker": {"tolerance": 0, "kwargs_where": {}},
            "2__remover": {
                "kwargs_selector": {
                    "cf_dims": ["XY"],
                    "region": "region1_supp3",
                    "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}}},
                "kwargs_averager": {"cf_dims": ["XY"]},
            },
            "2__selector": {
                "cf_dims": ["T", "XY"],
                "bounds_t": "epoch",
                "bounds_z": None,
                "kwargs_select_xy": {"mask_only": False, "kwargs_sel": {}, "kwargs_where": {}},
                "kwargs_select_t": {"kwargs_sel": {}},
                "kwargs_select_z": {"kwargs_sel": {}}},
            "3__averager": {"cf_dims": ["T", "XY"]},
        },
    },
}
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
def diagnostic(
        input_param: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, str, list[str], None]]]]],
        dataset: str = "unknown",
        epoch: Union[tuple[str, str], slice] = None,
        experiment: str = "unknown",
        member: str = "unknown",
        project: str = "unknown",
        recipe_add: str = None,
        region1_main: str = None,
        region1_supp1: str = None,
        region1_supp2: str = None,
        region1_supp3: str = None,
        supplementary: bool = False,
        variable1: str = None,
        kwargs_saver: str = None,
        **kwargs):
    basics.log_info(inspect__stack(), "# " + "-" * 30 + " > enter metric")
    # get default parameters if they are not given
    needed_kwarg = ["regions_param", "variables_param"]
    for k in needed_kwarg:
        if k not in list(kwargs.keys()):
            kwargs[k] = default_arg_values(k)
    # check region, statistic and variable values as they are required
    epoch = set_default_epoch(epoch, default["epoch"])
    region1_main = set_default_str(region1_main, list(kwargs["regions_param"].keys()), default["region1_main"])
    region1_supp1 = set_default_str(region1_supp1, list(kwargs["regions_param"].keys()), default["region1_supp1"])
    region1_supp2 = set_default_str(region1_supp2, list(kwargs["regions_param"].keys()), default["region1_supp2"])
    region1_supp3 = set_default_str(region1_supp3, list(kwargs["regions_param"].keys()), default["region1_supp3"])
    variable1 = set_default_str(variable1, list(kwargs["variables_param"].keys()), default["variable1"])
    kwargs_saver = set_instance(kwargs_saver, dict, False, {})
    for k, v in default["kwargs_saver"].items():
        if k not in list(kwargs_saver.keys()):
            kwargs_saver[k] = v
    # check / format dict_input
    input_param = input_dictionary_formater(input_param)
    # update processor
    to_replace = {"epoch": epoch, "region1_main": region1_main, "region1_supp1": region1_supp1,
                  "region1_supp2": region1_supp2, "region1_supp3": region1_supp3, "variable1": variable1}
    processors_m = processors_dictionary_formater(processors_main, to_replace)
    processors_s = processors_dictionary_formater(processors_supp, to_replace)
    # fake loop to be able to break out if an error occurs
    for _ in range(1):
        # ------------------------------------------------
        # 1. Read files
        # ------------------------------------------------
        # 1.1 Create a dictionary with variables and files
        # print("input_param")
        # print(json__dumps(input_param, indent=4))
        input_dataset = pr.reader(input_param, variable1, **kwargs)
        # print("processors.reader")
        # print(type(input_dataset))
        # print(list(input_dataset.keys()))
        # for k1, d1 in input_dataset.items():
        #     print(str(k1).rjust(15), type(d1))
        #     for k2, d2 in d1.items():
        #         print(str(k2).rjust(20), type(d2))
        #         if isinstance(d2, dict):
        #             print(json__dumps(d2, indent=4))
        # process
        # print("processors_m")
        # print(json__dumps(processors_m, indent=4))
        processed_main = pr.loop(processors_m, input_dataset, input_param, **kwargs)
        processed_supp = {}
        if isinstance(supplementary, bool) and supplementary:
            processed_supp = pr.loop(processors_s, input_dataset, input_param, **kwargs)
        # save as netCDF
        processed_main = {**processed_main, **processed_supp}
        l1 = ["dataset", "experiment", "member", "project", "recipe_add"]
        l2 = [dataset, experiment, member, project, recipe_add]
        tmp_kwargs = {k: v for k, v in zip(l1, l2) if isinstance(v, str)}
        pr.saver(processed_main, recipe=os__path__basename(__file__).split(".py")[0], **kwargs_saver, **tmp_kwargs)


def metric(
        input_dataset: str,
        input_reference: Union[str, dict[str, str]],
        metric_method: str = None,
        metric_dictionary: dict = None,
        metric_dictionary_keys: tuple[str] = None,
        **kwargs) -> dict:
    basics.log_info(inspect__stack(), "# " + "-" * 30 + " > enter metric")
    # get default parameters if they are not given
    needed_kwarg = ["metric_param"]
    for k in needed_kwarg:
        if k not in list(kwargs.keys()):
            kwargs[k] = default_arg_values(k)
    # check region, statistic and variable values as they are required
    metric_method = set_default_str(metric_method, kwargs["metric_param"], default["metric_method"])
    # compute metric
    metric_dictionary = pr.comparator(input_dataset, input_reference, metric_method,
                                      metric_dictionary=metric_dictionary,
                                      metric_dictionary_keys=metric_dictionary_keys)
    return metric_dictionary
# ---------------------------------------------------------------------------------------------------------------------#
