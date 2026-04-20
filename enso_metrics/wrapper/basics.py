# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools used by wrapper
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
import logging
from os.path import isdir as os__path__isdir
from os.path import join as os__path__join
from re import split as re__split
from string import ascii_lowercase as string__ascii_lowercase
from typing import Any, Hashable, Union
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
log = logging.getLogger(__name__)
logging.disable(logging.ERROR)
letters_order = string__ascii_lowercase + string__ascii_lowercase.upper()


def description_writer(description: str, text: str, restart_count: bool = False) -> str:
    """
    Add input ‘text’ to ‘description’ using my template

    Input:
    ------
    :param description: str
        Current description of the computation steps; e.g., description = ';; a) SST ;; b) masked over land'
    :param text: str
        New computation step to add; e.g., computation = 'selected in region'
    :param restart_count: bool, optional
        True to restart the count of computation steps; e.g., restart_count = False
        Default is False

    Output:
    -------
    :return: str
        Updated description of the computation steps
    """
    # find next letter number
    counter = 0
    if restart_count is False and isinstance(description, str) is True and len(description) > 0:
        letter = str(description.split(";; ")[-1].split(") ")[0])
        if letter.isalpha() and letter.islower():
            counter = ord(letter) - 96
        elif letter.isalpha() and letter.isupper():
            counter = str(ord(letter) - 38)
    # get next letter
    letter = letters_order[counter]
    # update description
    description += ";; " + str(letter) + ") " + str(text)
    return description


def flatten_dict(
        dict_i: dict[str, Any],
        parent_key: str = "",
        separator: str ="__",
        **kwargs) -> dict[str, str]:
    """
    Flatten nested dictionary by joining keys with a separator.
    This function recursively traverses a dictionary and flattens it so that nested keys are concatenated into a single
    key using the specified separator.

    Input:
    ------
    :param dict_i: dict[str, Any]
        The input dictionary to flatten. Values can be nested dictionaries or any other type.
    :param parent_key: str, optional
        The base key to prepend to all keys (used internally during recursion).
        Default is an empty string.
    :param separator: str, optional
        The separator used to join nested keys.
        Default is "__".
    **kwargs - Discarded

    Output:
    -------
    :return: dict[str, str]
        A flattened dictionary where nested keys are merged into a single level.
    """
    dict_o = {}
    for k, v in dict_i.items():
        new_key = f"{parent_key}{separator}{k}" if parent_key else k
        if isinstance(v, dict):
            dict_o.update(flatten_dict(v, new_key, separator=separator))
        else:
            dict_o[new_key] = str(v)
    return dict_o


def generate_filename(
        filename_netcdf: str = None,
        extension: str = "nc",
        path: str = None,
        project: str = None,
        dataset: str = None,
        experiment: str = None,
        member: str = None,
        recipe: str = None,
        recipe_add: str = None,
        **kwargs) -> str:
    """
    Generate a file name based on inputs.
    If filename_netcdf is provided:
        - path may be added (if path is given and filename_netcdf does not contain one)
        - extension may be added (if extension is given and filename_netcdf does not contain extension)
        - recipe may be added before extension (if recipe is given and filename_netcdf does not contain recipe)
    If filename_netcdf is NOT provided:
        - join all available keywords
        - path may be added (if path is given)
        - extension may be added (if extension is given)

    Input:
    ------
    :param filename_netcdf: str, optional
        User defined output filename_netcdf.
        If not provided other keywords must be provided
        Default is None
    :param extension: str, optional
        File extension (usually 'nc' or 'json')
        Default is 'nc'
    :param path: str, optional
        User defined output path.
        If path not included in filename_netcdf or filename_netcdf not provided, path must be provided
        Default is None
    :param project: str, optional
        Project name.
        Default is None
    :param dataset: str, optional
        Dataset name.
        Default is None
    :param experiment: str, optional
        Experiment name.
        Default is None
    :param member: str, optional
        Member name.
        Default is None
    :param recipe: str, optional
        Recipe name.
        Default is None
    :param recipe_add: str, optional
        Something to add to the recipe name.
        Default is None
    :param kwargs: dict[str, str], optional

    Output:
    -------
    :return: str
        Output filename for the output file
    """
    recipe_t = copy__deepcopy(recipe)
    if isinstance(recipe_t,  str) and isinstance(recipe_add, str):
        recipe_t += recipe_add
    # check if filename_netcdf was given
    if isinstance(filename_netcdf, str) and len(filename_netcdf) > 0:
        fo = filename_netcdf
    else:
        # regroup all values
        lt = [project, dataset, experiment, member, recipe_t] + list(kwargs.values())
        # join str
        fo = "_".join([k for k in lt if isinstance(k, str) and len(k) > 0])
    # check if path is given
    if len(fo.split("/")) < 2 and isinstance(path, str) and len(path) > 0 and os__path__isdir(path):
        fo = os__path__join(path, fo)
    # check if extension is on the filename
    if isinstance(extension, str) and len(extension) > 0 and fo.split(".")[-1] != extension:
        fo += "." + str(extension)
    # check if recipe name is in the filename
    fo_ext = fo.split(".")[-1]
    if isinstance(recipe_t, str) and len(recipe_t) > 0 and recipe_t not in fo and 1 < len(fo_ext) < 4:
        fo = fo.replace("." + str(fo_ext), "_" + str(recipe_t) + "." + str(fo_ext))
    return fo


def is_dim(dim: Any) -> bool:
    """
    Do a test on dimension type to decide if it is probably one

    Input:
    ------
    :param dim: Any
        Object to test, usually a string

    Output:
    -------
    :return: bool
    """
    return dim is not None and isinstance(dim, (Hashable, str)) is True


def log_debug(stack, message: str):
    log.debug(" function " + str(stack[0][3]) + " ; line " + str(stack[0][2]) + "\n" + str(message))


def log_details(list_names: list[str], list_params: list[Any], dict_o: dict = None) -> dict[str, str]:
    if dict_o is None:
        dict_o = {}
    for k1, k2 in zip(list_names, list_params):
        dict_o[str(k1) + ".type"] = str(type(k2))
        if isinstance(k2, dict):
            dict_o[str(k1) + ".keys"] = ", ".join(sorted(list(k2.keys()), key=lambda v: v.lower()))
        elif isinstance(k2, str) is True or k2 is None:
            dict_o[k1] = str(k2)
    return dict_o


def log_error(stack, message: str):
    log.error(" function " + str(stack[0][3]) + " ; line " + str(stack[0][2]) + "\n" + str(message))


def log_info(stack, message: str):
    message_o = " function " + str(stack[0][3]) + " ; line " + str(stack[0][2])
    if isinstance(message, str) is True and len(message) > 0:
        message_o += "\n" + str(message)
    log.info(message_o)


def merge_metadata(meta_dict: dict[str, dict[str, str]], sep: str = "__") -> dict[str, str]:
    """
    Merge metadata dictionaries for multiple variables.

    Input:
    ------
    :param meta_dict: dict[str, dict[str, str]]
        Dictionary mapping variable names to their metadata dictionaries.
        Example: meta_dict = {
            "t": {"name": "Temperature", "model": "IPSL", "unit": "K"},
            "p": {"name": "Pressure", "model": "IPSL"},
        }
    :param sep: str, optional
        Separator used when expanding keys with variable names.
        Example: sep = "__"
        output = {"model": "IPSL", "name__t": "Temperature", "name__p": "Pressure", "unit__t": "K"}
        Default is "__"

    Output:
    -------
    :return: dict[str, str]
        A merged dictionary following these rules:
        1. If a key exists in ALL variables AND all values are identical: keep a single entry: key: value
        2. Otherwise: create one entry per variable (key{sep}variable_name: value)

    Notes
    -----
        - Keys missing from some variables are treated as "not common", and therefore expanded per variable.
        - Output dictionary is sorted by key.
    """
    # --- Step 1: collect all unique keys across all variables ---
    all_keys = set()
    for metadata in meta_dict.values():
        all_keys.update(metadata.keys())
    # --- Step 2: process each key independently ---
    merged = {}
    for key in all_keys:
        # Gather values for this key across variables
        values_by_var = {var: metadata[key] for var, metadata in meta_dict.items() if key in metadata}
        # --- Step 3: check if key is common and identical everywhere ---
        is_common = len(values_by_var) == len(meta_dict)
        is_identical = len(set(values_by_var.values())) == 1
        if is_common and is_identical:
            # Same key/value in all variables → keep once
            merged[key] = next(iter(values_by_var.values()))
        else:
            # Otherwise: expand per variable
            for var, value in values_by_var.items():
                new_key = f"{key}{sep}{var}"
                merged[new_key] = value
    # --- Step 4: return a sorted dictionary ---
    return dict(sorted(merged.items()))


def set_nested_value(
        dict_i: dict[str, Any],
        value: Any,
        *args) -> dict[str, Any]:
    """
    Set a value in a nested dictionary using a list of keys.
    This function creates intermediate dictionaries as needed and assigns the given value at the deepest level defined
    by `keys`.

    Input:
    ------
    :param dict_i: dict[str, Any]
        The dictionary to modify (it will be updated in place).
    :param value: Any
        The value to assign at the nested location.
    *args: str
        A list of keys representing the path where the value should be set.

    Output:
    -------
    :return: dict[str, Any]
        The updated dictionary (same object as input).
    """
    current = dict_i
    # Traverse all keys except the last one
    for key in args[:-1]:
        # If the key doesn't exist or is not a dict, replace it with a dict
        if key not in current or not isinstance(current[key], dict):
            current[key] = {}
        # Move one level deeper
        current = current[key]
    # Set the value at the final key
    current[args[-1]] = value
    return dict_i


def split_time_bound(time_bound: str) -> list[str]:
    """
    Split input ‘time_bound’ using re.split
    https://docs.python.org/3/library/re.html#re.split

    Input:
    ------
    :param time_bound: str
        Time bound; e.g., time_bound = '1980-01-01 12:00:00'

    Output:
    -------
    :return: list[str]
        Input ‘time_bound’ split using ' ', ':', '-'
    """
    return re__split("[ :-]", time_bound)


def write_coordinates(
        lat1: Union[float, int],
        lat2: Union[float, int],
        lon1: Union[float, int],
        lon2: Union[float, int]) -> str:
    """
    Format input coordinates using my template

    Input:
    ------
    :param lat1: float, int
        Minimum latitude of the region
    :param lat2: float, int
        Maximum latitude of the region
    :param lon1: float, int
        Minimum longitude of the region
    :param lon2: float, int
        Maximum longitude of the region

    Output:
    -------
    :return:
    """
    # convert to integer (if applicable)
    lat1 = int(lat1) if lat1 == int(lat1) else lat1
    lat2 = int(lat2) if lat2 == int(lat2) else lat2
    lon1 = int(lon1) if lon1 == int(lon1) else lon1
    lon2 = int(lon2) if lon2 == int(lon2) else lon2
    # write latitudes
    if lat1 >= 0 and lat2 >= 0:
        lat_txt = str(lat1) + "-" + str(lat2) + "N"
    elif lat1 <= 0 and lat2 <= 0:
        lat_txt = str(abs(lat2)) + "-" + str(abs(lat1)) + "S"
    else:
        lat_txt = str(abs(lat1)) + "S-" + str(lat2) + "N"
    # write longitudes
    if (0 <= lon1 <= 180 and 0 <= lon2 <= 180) or (lon1 == 0 and lon2 == 360):
        lon_txt = str(lon1) + "-" + str(lon2) + "E"
    elif 180 <= lon1 <= 360 and 180 <= lon2 <= 360:
        lon_txt = str(360 - lon2) + "-" + str(360 - lon1) + "W"
    elif 0 <= lon1 <= 180 <= lon2 <= 360:
        lon_txt = str(lon1) + "E-" + str(360 - lon2) + "W"
    elif -180 < lon1 < 0 <= lon2 <= 360:
        lon_txt = str(360 + lon1) + "W-" + str(lon2) + "E"
    else:
        lon_txt = str(lon1) + "-" + str(lon2)
    return "[" + str(lat_txt) + " ; " + str(lon_txt) + "]"
# ---------------------------------------------------------------------------------------------------------------------#
