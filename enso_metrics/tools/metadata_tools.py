# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools for metadata
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from string import ascii_lowercase as string__ascii_lowercase

# local functions
from enso_metrics.tools.default import set_instance
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def method_writer(method: str, text: str, variable="") -> str:
    # initialize method and counter
    lc = 1
    if isinstance(method, str) is True and len(method) > 0:
        lc = int(method.split(";; ")[-1].split(") ")[0])
    else:
        method = str(variable)
    # get new letter
    letter = string__ascii_lowercase[lc]
    # update method
    method += ";; " + str(letter) + ") " + str(text)
    return method
# ---------------------------------------------------------------------------------------------------------------------#
