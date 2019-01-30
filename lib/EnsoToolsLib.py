# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
# ENSO_metrics package functions:
import EnsoErrorsWarnings


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of functions without UVCDAT
#
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

    if isinstance(string_or_list, basestring):
        # 'string_or_list' is a string
        if string_or_list not in dictionary.keys():
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
            if key not in dictionary.keys():
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
            ep_event = [1 for tt in range(len(val_longitude)) if val_longitude[tt] > 360 + threshold]
        else:
            ep_event = [1 for tt in range(len(val_longitude)) if val_longitude[tt] > threshold]
        ep_event = sum(ep_event) * 100. / len(val_longitude)
    return ep_event, keyerror

