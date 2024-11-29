# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import array as NUMPYarray
from numpy import mean as NUMPYmean
from numpy import square as NUMPYsquare
from numpy import std as NUMPYstd
from numpy import unravel_index as NUMPYunravel_index
from numpy import var as NUMPYvar
from scipy.stats import skew as SCIPYstats__skew
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile
# ENSO_metrics package functions:
from . import EnsoErrorsWarnings


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of functions without CDAT
#
def add_up_errors(list_keyerror):
    """
    #################################################################################
    Description:
    Adds in one string the given keyerrors from the list of keyerrors
    #################################################################################

    :param list_keyerror: list of string
        list of keyerrors (strings or None) listing encountered errors

    :return keyerror: string
        string of all encountered errors concatenated
    """
    keyerror = ""
    for key in list_keyerror:
        if len(keyerror) > 0 and key is not None:
            keyerror += " ; "
        if key is not None:
            keyerror += str(key)
    if keyerror == "":
        keyerror = None
    return keyerror


def find_xy_min_max(tab, return_val="both"):
    """
    #################################################################################
    Description:
    Finds in tab the position (t,x,y,z) of the minimum (return_val='mini') or the maximum (return_val='maxi') or both
    values (if return_val is neither 'mini' nor 'maxi')
    Returned position(s) are not the position in tab but in the (t,x,y,z) space defined by tab axes
    #################################################################################

    :param tab: masked_array
        array for which you would like to know the position (t,x,y,z) of the minimum and/or the maximum values
    :param return_val: string, optional
        'mini' to return the position of the minimum value
        'maxi' to return the position of the maximum value
        to return both minimum and maximum values, pass anything else
        default value = 'both', returns both minimum and maximum values

    :return: minimum/maximum position or both minimum and maximum positions, int, float or list
        position(s) in the (t,x,y,z) space defined by tab axes of the minimum and/or maximum values of tab
    """
    mini = NUMPYunravel_index(tab.argmin(), tab.shape)
    maxi = NUMPYunravel_index(tab.argmax(), tab.shape)
    list_ax = NUMPYarray([tab.getAxis(ii)[:] for ii in range(len(mini))])
    axis_min = [list_ax[ii][mini[ii]] for ii in range(len(mini))]
    axis_max = [list_ax[ii][maxi[ii]] for ii in range(len(mini))]
    if return_val == "mini":
        if len(axis_min) == 1:
            tab_out = axis_min[0]
        else:
            tab_out = axis_min
    elif return_val == "maxi":
        if len(axis_max) == 1:
            tab_out = axis_max[0]
        else:
            tab_out = axis_max
    else:
        if len(axis_min) == 1 and len(axis_max) == 1:
            tab_out = [axis_min[0], axis_max[0]]
        else:
            tab_out = [axis_min, axis_max]
    return tab_out


def math_metric_computation(model, model_err, obs=None, obs_err=None, keyword="difference"):
    """
    #################################################################################
    Description:
    Computes the metric value, i.e., distance between a model and an observational dataset
    #################################################################################

    :param model: float or None
        scalar value computed with a model
    :param model_err: float or None
        error value on the scalar value computed with a model
    :param obs: float or None, optional
        scalar value computed with an observational dataset
        default value is None
    :param obs_err: float or None, optional
        error value on the scalar value computed with an observational dataset
        default value is None
    :param keyword: string, optional
        name of a mathematical method to apply to compute the metric value (distance between a model and an
        observational dataset): 'difference', 'ratio', 'relative_difference', 'abs_relative_difference'
        default value is 'difference'
    :return metric: float or None
        distance between a model and an observational dataset or None if 'model' and/or 'obs' is/are None
    :return metric_err: float or None
        error on the distance between a model and an observational dataset or None if 'model' and/or 'obs' and/or
        'metric_err' and/or 'obs_err' is/are None
    :return description_metric: string
        description of the mathematical method used to compute the metric value
    """
    if keyword not in ["difference", "ratio", "relative_difference", "abs_relative_difference"]:
        metric, metric_err, description_metric = \
            None, None, "unknown keyword for the mathematical computation of the metric: " + str(keyword)
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": keyword",
                        str().ljust(5) + description_metric]
        EnsoErrorsWarnings.my_warning(list_strings)
    else:
        if model is not None and obs is not None:
            if keyword == "difference":
                description_metric = \
                    "The metric is the difference between model and observations values (M = model - obs)"
                metric = model - obs
            elif keyword == "ratio":
                description_metric = "The metric is the ratio between model and observations values (M = model / obs)"
                metric = model / float(obs)
            elif keyword == "relative_difference":
                description_metric = "The metric is the relative difference between model and observations values " + \
                                     "(M = [model-obs] / obs)"
                metric = (model - obs) / float(obs)
            else:
                description_metric = "The metric is the absolute value of the relative difference between model " + \
                                     "and observations values (M = 100 * abs[[model-obs] / obs])"
                metric = 100. * abs((model - obs) / float(obs))
        else:
            metric, description_metric = None, ""
        if model_err is not None or obs_err is not None:
            if keyword == "difference":
                # mathematical definition of the error on addition / subtraction
                metric_err = model_err + obs_err
            elif keyword == "ratio":
                # mathematical definition of the error on division
                if model is not None and obs is not None:
                    metric_err = float((obs * model_err + model * obs_err) / NUMPYsquare(obs))
                else:
                    metric_err = None
            else:
                # mathematical definition of the error on division
                if model is not None and obs is not None:
                    metric_err = float((obs * (model_err + obs_err) + (model - obs) * obs_err) / NUMPYsquare(obs))
                    if keyword == "abs_relative_difference":
                        metric_err = 100. * metric_err
                else:
                    metric_err = None
        else:
            metric_err = None
    return metric, metric_err, description_metric


def percentage_val_eastward(val_longitude, metric_name, region, threshold=-140):
    """
    #################################################################################
    Description:
    Computes the percentage of given values (longitude in val_longitude) eastward of 'threshold'
    #################################################################################

    :param val_longitude: list
        list of longitude
    :param metric_name: string, optional
        name of the metric calling the function
    :param region: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param threshold: float, optional
        threshold to define the westward boundary of the region

    :return ep: float
        percentage of given values eastward of 'threshold'
    """
    keyerror = None
    pos_lon = "pos" if all(ii >= 0 for ii in val_longitude) is True else \
        ("neg" if all(ii <= 0 for ii in val_longitude) is True else False)
    if pos_lon is False:
        keyerror = "longitude in lon_axis are neither all positive nor all negative"
        list_strings = ["ERROR " + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": unknown longitude",
                        str().ljust(5) + metric_name + ": " + str(region) + ": " + keyerror + ": " + str(val_longitude)]
        EnsoErrorsWarnings.my_warning(list_strings)
        ep_event = None
    else:
        if pos_lon == "pos":
            ep_event = [1 for val in val_longitude if val > 360 + threshold]
        else:
            ep_event = [1 for val in val_longitude if val > threshold]
        ep_event = sum(ep_event) * 100. / len(val_longitude)
    return ep_event, keyerror


def simple_stats(arr, axis=0):
    ave = NUMPYmean(arr, axis=axis)
    ske = SCIPYstats__skew(arr, axis=axis)
    std = NUMPYstd(arr, axis=axis)
    var = NUMPYvar(arr, axis=axis)
    return ave, ske, std, var


def statistical_dispersion(tab, method="IQR"):
    """
    #################################################################################
    Description:
    Computes the statistical dispersion of the distribution
    #################################################################################

    :param tab: list or `cdms2` variable
        A list or a `cdms2` variable containing the data to be analysed
    :param method: string, optional
        method to compute the statistical dispersion
        'IQR': interquartile range, IQR = Q3 - Q1
        'MAD': median absolute deviation, MAD = median([Xi - median(tab)])
        Default is 'IQR'

    :return stat_disp: float
        statistical_dispersion
    """
    known_methods = sorted(["IQR", "MAD"])
    if method not in known_methods:
        # if 'method' is not defined -> raise error
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": unknown method",
            str().ljust(5) + "method " + str(method) + " is not defined",
            str().ljust(10) + "known methods: " + str(known_methods)]
        EnsoErrorsWarnings.my_error(list_strings)
    if method == "IQR":
        stat_disp = abs(float(SCIPYstats__scoreatpercentile(tab, 75) - SCIPYstats__scoreatpercentile(tab, 25)))
    else:
        med = float(SCIPYstats__scoreatpercentile(tab, 50))
        stat_disp = float(SCIPYstats__scoreatpercentile([abs(ii - med) for ii in tab], 50))
    return stat_disp


def string_in_dict(string_or_list, dictionary, inspect_stack):
    """
    #################################################################################
    Description:
    Tests if 'string_or_list' is in the given 'dictionary'
    #################################################################################

    :param string_or_list: string or list
        key or list of keys to look for in 'dictionary'
    :param dictionary: dict
        dictionary in which the given keys are looked for
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    # test input parameters
    if not isinstance(dictionary, dict):
        EnsoErrorsWarnings.object_type_error("dictionary", "dictionary", type(dictionary), INSPECTstack())
    if isinstance(string_or_list, str):
        # 'string_or_list' is a string
        if string_or_list not in list(dictionary.keys()):
            # key 'string_or_list' is not in 'dictionary' -> raise error
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect_stack) + ": item not included",
                str().ljust(5) + " key: " + str(string_or_list) + " is not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    elif isinstance(string_or_list, list):
        # 'string_or_list' is a list
        key_not_included = list()
        for key in string_or_list:
            if key not in list(dictionary.keys()):
                # lists keys that are not in 'dictionary'
                key_not_included.append(key)
        if key_not_included:
            # if 'key_not_included' is not empty -> raise error
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect_stack) + ": item not included",
                str().ljust(5) + " key(s): " + str(key_not_included) + " are (is) not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))]
            EnsoErrorsWarnings.my_error(list_strings)
    else:
        # 'string_or_list' is neither a string nor a list -> raise error
        EnsoErrorsWarnings.object_type_error("string_or_list", "[string, list]", type(string_or_list), INSPECTstack())
    return
