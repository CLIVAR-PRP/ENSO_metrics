# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from typing import Any

# local functions
from enso_metrics.tools.default import set_instance
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def combine_dict_levels(
        dict_i: dict,
        dict_o: dict = None,
        tuple_k: tuple = None,
        tuple_k_last: tuple = None,
        **kwargs) -> (dict, tuple, tuple):
    # set None input to it's default value
    dict_o = set_instance(dict_o, dict, False, {})
    tuple_k = set_instance(tuple_k, tuple, False, ())
    tuple_k_last = set_instance(tuple_k_last, tuple, False, ())
    # loop through nested levels
    if isinstance(dict_i, dict) is True:
        list_keys = sorted(list(dict_i.keys()), key=lambda s: s.lower())
        for k in list_keys:
            dict_o, tuple_k, tuple_k_last = combine_dict_levels(
                dict_i[k], dict_o=dict_o, tuple_k=tuple_k + (k,), tuple_k_last=tuple_k_last + (list_keys[-1],))
    else:
        # save values
        kk = "__".join(tuple_k)
        put_in_dict(dict_o, dict_i, kk)
        # remove relevant keys from the tuples of keys
        tuple_k, tuple_k_last = tuple_for_dict(tuple_k, tuple_k_last)
    return dict_o, tuple_k, tuple_k_last


def put_in_dict(dict_i: dict, value: Any, *args):
    """
    Put value in the dictionary

    Input:
    ------
    :param dict_i: dict
        Dictionary in which the value must be added
    :param value: Any
        Value to add in the dictionary
        If it is a list, it will be appended to the list already inside the dictionary
    :param args: str
        Non keyword arguments used as keys in the dictionary
    """
    # put value in the dictionary
    _dict = dict_i
    for k in args:
        if k == args[-1]:
            if k in list(_dict.keys()) and isinstance(_dict[k], list) is True and isinstance(value, list) is True:
                _dict[k] += value
            else:
                _dict[k] = value
        else:
            if k not in list(_dict.keys()):
                _dict[k] = {}
            _dict = _dict[k]


def sort_dict(
        dict_i: dict,
        dict_o: dict = None,
        tuple_k: tuple = None,
        tuple_k_last: tuple = None,
        **kwargs) -> (dict, tuple, tuple):
    """
    Sort dictionary

    Input:
    ------
    :param dict_i: dict
        Dictionary to sort
    :param dict_o: dict or None, optional
        Dictionary in which output values will be stored
    :param tuple_k: tuple or None, optional
        List of keys, in order, for the output nested dictionary
    :param tuple_k_last: tuple or None, optional
        List of the last key of each nested level, in order, to keep track of the position within the input nested
        dictionary

    Outputs:
    --------
    :return dict_o: dict
        Sorted dictionary
    :return tuple_k: tuple
        List of keys, in order, for the output nested dictionary
    :return tuple_k_last: tuple
        List of the last key of each nested level, in order, to keep track of the position within the input nested
        dictionary
    """
    # set None input to it's default value
    dict_o = set_instance(dict_o, dict, False,{})
    tuple_k = set_instance(tuple_k, tuple, False, ())
    tuple_k_last = set_instance(tuple_k_last, tuple, False, ())
    # loop through nested levels
    if isinstance(dict_i, dict) is True:
        list_keys = sorted(list(dict_i.keys()), key=lambda s: s.lower())
        for k in list_keys:
            dict_o, tuple_k, tuple_k_last = sort_dict(
                dict_i[k], dict_o=dict_o, tuple_k=tuple_k + (k,), tuple_k_last=tuple_k_last + (list_keys[-1],))
    else:
        # save values
        put_in_dict(dict_o, dict_i, *tuple_k)
        # remove relevant keys from the tuples of keys
        tuple_k, tuple_k_last = tuple_for_dict(tuple_k, tuple_k_last)
    return dict_o, tuple_k, tuple_k_last


def tuple_for_dict(tuple_of_keys: tuple, tuple_of_last_key: tuple) -> (tuple, tuple):
    """
    Remove keys that reached the end of the list

    Input:
    ------
    :param tuple_of_keys: tuple
        Keys for nested dictionary
    :param tuple_of_last_key: tuple
        Keys of the last key of each nested level

    Output:
    -------
    :return tuple_of_keys: tuple
        Keys for nested dictionary, with last key(s) removed
    :param tuple_of_last_key: tuple
        Keys of the last key of each nested level, with last key(s) removed
    """
    # reverse the order of the tuple_of_last_key
    list_r = list(reversed(tuple_of_last_key))
    # check if the last key of a level has been reached, if yes, remove this level from the tuples
    for k in list_r:
        if tuple_of_keys[-1] == k:
            tuple_of_keys = tuple_of_keys[:-1]
            tuple_of_last_key = tuple_of_last_key[:-1]
    # remove last item in both tuples as it has been saved in the output dictionary
    tuple_of_keys = tuple_of_keys[:-1]
    tuple_of_last_key = tuple_of_last_key[:-1]
    return tuple_of_keys, tuple_of_last_key
# ---------------------------------------------------------------------------------------------------------------------#
