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

