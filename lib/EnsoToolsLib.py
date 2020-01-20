# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import array as NUMPYarray
from numpy import unravel_index as NUMPYunravel_index
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
    keyerror = ''
    for key in list_keyerror:
        if len(keyerror) > 0 and key is not None:
            keyerror += " ; "
        if key is not None:
            keyerror += str(key)
    return keyerror


def FindXYMinMax(tab, return_val='both'):
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
    if return_val == 'mini':
        if len(axis_min) == 1:
            tab_out = axis_min[0]
        else:
            tab_out = axis_min
    elif return_val == 'maxi':
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


def percentage_val_eastward(val_longitude, metric, region, threshold=-140):
    """
    #################################################################################
    Description:
    Computes the percentage of given values (longitude in val_longitude) eastward of 'threshold'
    #################################################################################

    :param val_longitude: list
        list of longitude
    :param threshold: float, optional
        threshold to define the westward boundary of the region
    :return ep: float
        percentage of given values eastward of 'threshold'
    """
    keyerror = None
    pos_lon = 'pos' if all(ii >= 0 for ii in val_longitude) is True else \
        ('neg' if all(ii <= 0 for ii in val_longitude) is True else False)
    if pos_lon is False:
        keyerror = "longitude in lon_axis are neither all positive nor all negative"
        list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": unknown longitude",
                        str().ljust(5) + metric + ": " + str(region) + ": " + keyerror + ": " + str(val_longitude)]
        EnsoErrorsWarnings.MyWarning(list_strings)
        ep_event = None
    else:
        if pos_lon == 'pos':
            ep_event = [1 for val in val_longitude if val > 360 + threshold]
        else:
            ep_event = [1 for val in val_longitude if val > threshold]
        ep_event = sum(ep_event) * 100. / len(val_longitude)
    return ep_event, keyerror


def statistical_dispersion(tab, method='IQR'):
    """
    #################################################################################
    Description:
    Computes the statistical dispersion of the distribution
    #################################################################################

    :param tab: list or `cdms2` variable
        A list or a `cdms2` variable containing the data to be analysed.
    :param method: string, optional
        method to compute the statistical dispersion
        'IQR': interquartile range, IQR = Q3 - Q1
        'MAD': median absolute deviation, MAD = median([Xi - median(tab)])
        Default is 'IQR'.
    :return stat_disp: float
        statistical_dispersion
    """
    known_methods = sorted(['IQR', 'MAD'])
    if method not in known_methods:
        # if 'method' is not defined -> raise error
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": unknown method",
            str().ljust(5) + "method " + str(method) + " is not defined",
            str().ljust(10) + "known methods: " + str(known_methods)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if method == 'IQR':
        stat_disp = abs(float(SCIPYstats__scoreatpercentile(tab, 75) - SCIPYstats__scoreatpercentile(tab, 25)))
    else:
        med = float(SCIPYstats__scoreatpercentile(tab, 50))
        stat_disp = float(SCIPYstats__scoreatpercentile([abs(ii - med) for ii in tab], 50))
    return stat_disp


def StringInDict(string_or_list, dictionary, inspect_stack):
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
        EnsoErrorsWarnings.ObjectTypeError('dictionary', 'dictionary', type(dictionary), INSPECTstack())

    if isinstance(string_or_list, str):
        # 'string_or_list' is a string
        if string_or_list not in list(dictionary.keys()):
            # key 'string_or_list' is not in 'dictionary' -> raise error
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(inspect_stack) + ": item not included",
                str().ljust(5) + " key: " + str(string_or_list) + " is not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))
            ]
            EnsoErrorsWarnings.MyError(list_strings)
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
                "ERROR" + EnsoErrorsWarnings.MessageFormating(inspect_stack) + ": item not included",
                str().ljust(5) + " key(s): " + str(key_not_included) + " are (is) not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    else:
        # 'string_or_list' is neither a string nor a list -> raise error
        EnsoErrorsWarnings.ObjectTypeError('string_or_list', '[string, list]', type(string_or_list), INSPECTstack())
    return

