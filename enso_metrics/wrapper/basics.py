# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools used by wrapper
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
import logging
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
