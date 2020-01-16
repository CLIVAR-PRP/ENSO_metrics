# -*- coding:UTF-8 -*-
# ENSO_metrics package functions:
from EnsoCollectionsLib import ReferenceRegions


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


def closest_grid(region, nlat, nlon):
    """
    #################################################################################
    Description:
    Finds the closest regular (e.g., 0.5x0.5deg, 1x1deg) generic grid that would fit the given data
    #################################################################################

    :param region: string
        name of the region (must be defined in EnsoCollectionsLib.ReferenceRegions)
    :param nlat: integer
        number of cells in latitude
    :param nlon: integer
        number of cells in longitude

    :return grid: string
        generic grid name (will be used by EnsoUvcdatToolsLib.Regrid to generate the grid)
    """
    res = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75]
    region_ref = ReferenceRegions(region)
    lats = region_ref['latitude']
    dy = float(abs(max(lats) - min(lats))) / nlat
    lyy = [abs(dy - ii) for ii in res]
    lyy = res[lyy.index(min(lyy))]
    lons = region_ref['longitude']
    dx = float(abs(max(lons) - min(lons))) / nlon
    lxx = [abs(dx - ii) for ii in res]
    lxx = res[lxx.index(min(lxx))]
    if lxx == lyy:
        grid = "generic_" + str(lxx) + "x" + str(lxx) + "deg"
    else:
        dx = abs(lxx + lyy) / 2.
        lxx = [abs(dx - ii) for ii in res]
        lxx = res[lxx.index(min(lxx))]
        grid = "generic_" + str(lxx) + "x" + str(lxx) + "deg"
    return grid
