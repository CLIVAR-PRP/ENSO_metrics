# -*- coding:UTF-8 -*-
from copy import deepcopy
# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from inspect import stack as inspect__stack
from typing import Any, Union

# local functions
from enso_metrics.tools import prints
from enso_metrics.definitions.regions import regions_param
from enso_metrics.definitions.variables import variables_cmip, variables_observation, variables_param
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def default_arg_values(arg):
    default_values = {
        "detrending": False,
        "enso_definition": {
            "duration_min": 5,
            "interannual_anomalies": True,
            "normalization": True,
            "region_ev": "nino3.4",
            "season_ev": "NDJ",
            "smoothing": False,
            "threshold": 0.5,
        },
        "frequency": None,
        "metric_computation": "difference",
        "min_time_steps": None,
        "normalization": False,
        "project_interpreter": "CMIP",
        "regions_param": regions_param(),
        "regridding": False,
        "smoothing": False,
        "threshold_ep_ev": -140,
        "time_bounds": None,
        "time_bounds_mod": None,
        "time_bounds_obs": None,
        "variables_cmip": variables_cmip(),
        "variables_observation": variables_observation(),
        "variables_param": variables_param(),
        "wait_definition": {
            "detect": "en_to_en",
            "method": "peaks",
            "smoothing": {
                "method": "triangle",
                "window": 5,
            },
        },
    }
    if arg not in list(default_values.keys()):
        prints.unknown_key_arg(arg, inspect__stack())
    return default_values[arg]


def input_dictionary_formater(
        dict_input: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, None]]]]],
        **kwargs) -> dict[str, dict[str, Union[str, list[str], None, dict[str, Union[int, float]]]]]:
    dict_o = {}
    # format output dictionary
    for k1 in list(dict_input.keys()):
        # get file and variable name(s) from input dictionary
        files, names = dict_input[k1]["file_name"], dict_input[k1]["variable"]
        # if multiple netCDF variables (names) are required to compute given internal variable (k1), inputs are lists
        # if not, inputs are str. To facilitate the process, lists are create anyway
        if isinstance(names, list) is False:
            files, names = [files], [names]
        # skip if at least one file or variable names is None
        if None in files or None in names:
            continue
        # fill dictionary
        dict_o[k1] = {"file_name": files, "variable": names}
        dict_o[k1]["area"] = dict_input[k1]["area"] if "area" in list(dict_input[k1].keys()) else None
        dict_o[k1]["mask"] = dict_input[k1]["mask"] if "mask" in list(dict_input[k1].keys()) else "estimate"
        for k2 in ["variable_offset", "variable_scaling"]:
            dict_o[k1][k2] = {}
            for k3 in names:
                if k2 in list(dict_input[k1].keys()) and isinstance(dict_input[k1][k2], (float, int)) is True:
                    dict_o[k1][k2][k3] = deepcopy(dict_input[k1][k2])
                elif k2 in list(dict_input[k1].keys()) and isinstance(dict_input[k1][k2], dict) is True and \
                        k3 in list(dict_input[k1][k2].keys()):
                    dict_o[k1][k2][k3] = deepcopy(dict_input[k1][k2][k3])
                else:
                    dict_o[k1][k2][k3] = 0 if k2 == "variable_offset" else 1
        if "variable_computation" in list(dict_input[k1].keys()):
            dict_o[k1]["variable_computation"] = deepcopy(dict_input[k1]["variable_computation"])
        else:
            # format add offset
            ao = sum([dict_o[k1]["variable_scaling"][k2] for k2 in names])
            a3 = ""
            if ao != 0:
                a3 = " - " if ao < 0 else " + "
                a3 += str(abs(ao))
            # format scale factors and variables names
            if len(list(set([dict_o[k1]["variable_scaling"][k2] for k2 in names]))) == 1:
                # -- there is only one scale factor, factorize it
                # format scale factor
                a1 = ""
                if dict_o[k1]["variable_scaling"][names[0]] != 1:
                    a1 = str(dict_o[k1]["variable_scaling"][names[0]]) + " * "
                # format variables names
                a2 = " + ".join(names)
                if len(names) > 1 and a1 != "":
                    a2 = " (" + str(a2) + " )"
            else:
                # -- multiple scale factors, write them in from of each variable name
                a2 = ""
                for k2 in names:
                    # scale factor
                    sf = dict_o[k1]["variable_scaling"][k2]
                    a1 = "" if k2 == names[0] else " "
                    if (k2 == names[0] and sf != 1) or k2 != names[0]:
                        a1 += "- " if sf < 0 else "+ "
                        if sf != 1:
                            a1 += str(abs(sf))
                    # scale factor * variable name
                    a2 += str(a1) + " * " + str(k2)
                a1 = ""
            # computation description
            dict_o[k1]["variable_computation"] = str(a1) + str(a2) + str(a3)
    return dict_o


def set_default_str(input_value: str, defined_values: list[str], optional_default: str, **kwargs) -> str:
    n = 0
    while input_value not in defined_values:
        input_value = deepcopy(optional_default)
        if n > 0:
            input_value = defined_values[0]
        n += 1
    return input_value


def set_instance(input_value: Any, test_type: type, test_bool: bool, default_value: Any) -> Any:
    if isinstance(input_value, test_type) is test_bool:
        input_value = deepcopy(default_value)
    return input_value
# ---------------------------------------------------------------------------------------------------------------------#
