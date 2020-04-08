# -*- coding:UTF-8 -*-
from __future__ import print_function
from calendar import monthrange
from copy import deepcopy
from datetime import date
from inspect import stack as INSPECTstack
import ntpath
from numpy import array as NParray
from numpy import exp as NPexp
from numpy import histogram as NPhistogram
from numpy import isnan as NPisnan
from numpy import nan as NPnan
from numpy import nonzero as NPnonzero
from numpy import ones as NPones
from numpy import product as NPproduct
from numpy import where as NPwhere
from numpy.ma.core import MaskedArray as NPma__core__MaskedArray
from os.path import isdir as OSpath_isdir
from os.path import isfile as OSpath__isfile
from os.path import join as OSpath__join
from os.path import split as OSpath__split
from scipy.signal import detrend as SCIPYsignal_detrend
from scipy.stats import skew as SCIPYstats__skew
from sys import prefix as SYS_prefix

# ENSO_metrics package functions:
from EnsoCollectionsLib import CmipVariables
from EnsoCollectionsLib import ReferenceObservations
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings
from EnsoToolsLib import add_up_errors, find_xy_min_max, string_in_dict

# uvcdat based functions:
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import createRectGrid as CDMS2createRectGrid
from cdms2 import createUniformLatitudeAxis as CDMS2createUniformLatitudeAxis
from cdms2 import createUniformLongitudeAxis as CDMS2createUniformLongitudeAxis
from cdms2 import createVariable as CDMS2createVariable
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
from cdtime import comptime as CDTIMEcomptime
import cdutil
from genutil.statistics import correlation as GENUTILcorrelation
from genutil.statistics import linearregression as GENUTILlinearregression
from genutil.statistics import rms as GENUTILrms
from genutil.statistics import std as GENUTILstd
from MV2 import add as MV2add
from MV2 import arange as MV2arange
from MV2 import array as MV2array
from MV2 import average as MV2average
from MV2 import compress as MV2compress
from MV2 import concatenate as MV2concatenate
from MV2 import divide as MV2divide
from MV2 import masked_where as MV2masked_where
from MV2 import maximum as MV2maximum
from MV2 import minimum as MV2minimum
from MV2 import multiply as MV2multiply
from MV2 import ones as MV2ones
from MV2 import subtract as MV2subtract
from MV2 import sum as MV2sum
from MV2 import take as MV2take
from MV2 import where as MV2where
from MV2 import zeros as MV2zeros
from regrid2.horizontal import Horizontal as REGRID2horizontal__Horizontal


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of simple uvcdat functions used in EnsoMetricsLib.py
#
def ArrayOnes(tab, id='new_variable_ones'):
    """
    #################################################################################
    Description:
    Create a masked_array filled with ones with the same properties as tab (shape, axes, grid, mask)
    #################################################################################

    for more information:
    import MV2
    help(MV2.ones)
    """
    return CDMS2createVariable(MV2ones(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask, id=id)


def ArrayZeros(tab, id='new_variable_zeros'):
    """
    #################################################################################
    Description:
    Create a masked_array filled with zeros with the same properties as tab (shape, axes, grid, mask)
    #################################################################################

    for more information:
    import MV2
    help(MV2.zeros)
    """
    return CDMS2createVariable(MV2zeros(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask, id=id)


def AverageHorizontal(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 'xy' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    lat_num = get_num_axis(tab, 'latitude')
    lon_num = get_num_axis(tab, 'longitude')
    snum = str(lat_num) + str(lon_num)
    if areacell is None:
        try: averaged_tab = cdutil.averager(tab, axis='xy', weights='weighted', action='average')
        except:
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageHorizontal" + "\033[0m")
            print("\033[93m" + str().ljust(20) + "axes = " + str(snum) + "\033[0m")
            try: averaged_tab = cdutil.averager(tab, axis=snum, weights='weighted', action='average')
            except:
                if 'regridding' not in kwargs.keys() or isinstance(kwargs['regridding'], dict) is False:
                    kwargs2 = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                               'newgrid_name': 'generic_1x1deg'}
                else:
                    kwargs2 = kwargs['regridding']
                kwargs2["newgrid_name"] =\
                    closest_grid(region, len(tab.getAxis(lat_num)[:]), len(tab.getAxis(lon_num)[:]))
                print("\033[93m" + str().ljust(25) + "need to regrid to = " + str(kwargs2["newgrid_name"]) +
                      " to perform average \033[0m")
                tmp = Regrid(tab, None, region=region, **kwargs2)
                try: averaged_tab = cdutil.averager(tmp, axis=snum, weights='weighted', action='average')
                except:
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
                        str().ljust(5) + "cdutil.averager cannot perform horizontal average"
                    ]
                    EnsoErrorsWarnings.my_error(list_strings)
    else:
        if tab.getGrid().shape != areacell.getGrid().shape:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
                str().ljust(5) + "tab (" + str(tab.getGrid().shape) + ") and areacell ("
                + str(areacell.getGrid().shape) + ") are not on the same grid",
                str().ljust(5) + "cannot perform horizontal average"
            ]
            EnsoErrorsWarnings.my_error(list_strings)
        averaged_tab = MV2multiply(tab, areacell)
        for elt in snum[::-1]:
            averaged_tab = MV2sum(averaged_tab, axis=int(elt))
        averaged_tab = averaged_tab / float(MV2sum(areacell))
    return averaged_tab


def AverageMeridional(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 'y' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    lat_num = get_num_axis(tab, 'latitude')
    lon_num = get_num_axis(tab, 'longitude')
    snum = str(lat_num)
    if areacell is None:
        try: averaged_tab = cdutil.averager(tab, axis='y', weights='weighted', action='average')
        except:
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageMeridional" + "\033[0m")
            print("\033[93m" + str().ljust(20) + "axes = " + str(snum) + "\033[0m")
            try: averaged_tab = cdutil.averager(tab, axis=snum, weights='weighted', action='average')
            except:
                if 'regridding' not in kwargs.keys() or isinstance(kwargs['regridding'], dict) is False:
                    kwargs2 = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                               'newgrid_name': 'generic_1x1deg'}
                else:
                    kwargs2 = kwargs['regridding']
                kwargs2["newgrid_name"] = \
                    closest_grid(region, len(tab.getAxis(lat_num)[:]), len(tab.getAxis(lon_num)[:]))
                print("\033[93m" + str().ljust(25) + "need to regrid to = " + str(kwargs2["newgrid_name"]) +
                      " to perform average \033[0m")
                tmp = Regrid(tab, None, region=region, **kwargs2)
                try: averaged_tab = cdutil.averager(tmp, axis=snum, weights='weighted', action='average')
                except:
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional average",
                        str().ljust(5) + "cdutil.averager cannot perform meridional average"
                    ]
                    EnsoErrorsWarnings.my_error(list_strings)
    else:
        if tab.getGrid().shape != areacell.getGrid().shape:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
                str().ljust(5) + "tab (" + str(tab.getGrid().shape) + ") and areacell ("
                + str(areacell.getGrid().shape) + ") are not on the same grid",
                str().ljust(5) + "cannot perform meridional average"
            ]
            EnsoErrorsWarnings.my_error(list_strings)
        lat_num_area = get_num_axis(areacell, 'latitude')
        averaged_tab = MV2multiply(tab, areacell)
        averaged_tab = MV2sum(averaged_tab, axis=lat_num) / MV2sum(areacell, axis=lat_num_area)
    lon = tab.getLongitude()
    if len(lon.shape) > 1:
        lonn = CDMS2createAxis(MV2array(lon[0, :]), id='longitude')
        lonn.units = lon.units
        lon_num = get_num_axis(tab, 'longitude')
        try:
            averaged_tab.setAxis(lon_num, lonn)
        except:
            averaged_tab.setAxis(lon_num - 1, lonn)
    return averaged_tab


def AverageTemporal(tab, areacell=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 't' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    try: averaged_tab = cdutil.averager(tab, axis='t')
    except:
        time_num = get_num_axis(tab, 'time')
        try: averaged_tab = cdutil.averager(tab, axis=str(time_num))
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal average",
                str().ljust(5) + "cannot perform temporal average"
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    return averaged_tab


def AverageZonal(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 'x' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    lat_num = get_num_axis(tab, 'latitude')
    lon_num = get_num_axis(tab, 'longitude')
    snum = str(lon_num)
    if areacell is None:
        try: averaged_tab = cdutil.averager(tab, axis='x', weights='weighted', action='average')
        except:
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageZonal" + "\033[0m")
            print("\033[93m" + str().ljust(20) + "axes = " + str(snum) + "\033[0m")
            try: averaged_tab = cdutil.averager(tab, axis=snum, weights='weighted', action='average')
            except:
                if 'regridding' not in kwargs.keys() or isinstance(kwargs['regridding'], dict) is False:
                    kwargs2 = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                               'newgrid_name': 'generic_1x1deg'}
                else:
                    kwargs2 = kwargs['regridding']
                kwargs2["newgrid_name"] = \
                    closest_grid(region, len(tab.getAxis(lat_num)[:]), len(tab.getAxis(lon_num)[:]))
                print("\033[93m" + str().ljust(25) + "need to regrid to = " + str(kwargs2["newgrid_name"]) +
                      " to perform average \033[0m")
                tmp = Regrid(tab, None, region=region, **kwargs2)
                try: averaged_tab = cdutil.averager(tmp, axis=snum, weights='weighted', action='average')
                except:
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": zonal average",
                        str().ljust(5) + "cdutil.averager cannot perform zonal average"
                    ]
                    EnsoErrorsWarnings.my_error(list_strings)
    else:
        if tab.getGrid().shape != areacell.getGrid().shape:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
                str().ljust(5) + "tab (" + str(tab.getGrid().shape) + ") and areacell ("
                + str(areacell.getGrid().shape) + ") are not on the same grid",
                str().ljust(5) + "cannot perform zonal average"
            ]
            EnsoErrorsWarnings.my_error(list_strings)
        lon_num_area = get_num_axis(areacell, 'longitude')
        averaged_tab = MV2multiply(tab, areacell)
        averaged_tab = MV2sum(averaged_tab, axis=lon_num) / MV2sum(areacell, axis=lon_num_area)
    lat = tab.getLatitude()
    if len(lat.shape) > 1:
        latn = CDMS2createAxis(MV2array(lat[:, 0]), id='latitude')
        latn.units = lat.units
        lat_num = get_num_axis(tab, 'latitude')
        try:
            averaged_tab.setAxis(lat_num, latn)
        except:
            averaged_tab.setAxis(lat_num - 1, latn)
    return averaged_tab


# Dictionary of averaging methods
dict_average = {'horizontal': AverageHorizontal, 'meridional': AverageMeridional, 'time': AverageTemporal,
                'zonal': AverageZonal}


def Concatenate(tab1, tab2, events1=[], events2=[]):
    my_events = events1 + events2
    if len(my_events) > 0:
        my_events_sort = sorted(my_events)
        for yy in my_events_sort:
            try:
                tab_out
            except:
                if yy in events1:
                    tab_out = MV2array([tab1[events1.index(yy)]])
                else:
                    tab_out = MV2array([tab2[events2.index(yy)]])
            else:
                if yy in events1:
                    tab_out = MV2concatenate((tab_out, MV2array([tab1[events1.index(yy)]])))
                else:
                    tab_out = MV2concatenate((tab_out, MV2array([tab2[events2.index(yy)]])))
        axes = CDMS2createAxis(MV2array(my_events_sort, dtype='int32'), id='years')
        if len(events1):
            tmp = deepcopy(tab1)
        else:
            tmp = deepcopy(tab2)
        att = tmp.attributes
        if len(tmp.shape) > 1:
            mask = tmp[0].mask
            mask2 = MV2zeros(tab_out.shape)
            mask2[:] = mask
            dictvar = {"axes": [axes] + tab1[0].getAxisList(), "mask": mask2, "grid": tmp.getGrid(), "attributes": att}
        else:
            dictvar = {"axes": [axes], "attributes": att}
        tab_out = CDMS2createVariable(tab_out, **dictvar)
    else:
        tab_out = MyEmpty(tab1[:5, 0], time=True, time_id='years')
    return tab_out

def closest_grid(region, nlat, nlon):
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


def ComputeInterannualAnomalies(tab):
    """
    #################################################################################
    Description:
    Computes interannual anomalies
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.ANNUALCYCLE.departures)
    """
    return cdutil.ANNUALCYCLE.departures(tab)


def Correlation(tab, ref, weights=None, axis=0, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes correlation
    #################################################################################

    for more information:
    import genutil
    help(genutil.statistics.correlation)
    """
    return GENUTILcorrelation(tab, ref, weights=weights, axis=axis, centered=centered, biased=biased)


def OperationAdd(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Adds every elements of 'tab' by 'number_or_tab'
    If 'number_or_tab' is an array it must have the same shape as tab
    #################################################################################

    for more information:
    import MV2
    help(MV2.add)
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2add(tab, number_or_tab)


def OperationDivide(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Divides every elements of 'tab' by 'number_or_tab'
    #################################################################################

    for more information:
    import MV2
    help(MV2.divide)
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2divide(tab, number_or_tab)


def OperationMultiply(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Multiplies every elements of 'tab' by 'number_or_tab'
    #################################################################################

    for more information:
    import MV2
    help(MV2.multiply)
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    tab_out = MV2multiply(tab, number_or_tab)
    axes = tab.getAxisList()
    att = tab.attributes
    if len(tab.shape) > 1:
        mask = tab.mask
        dictvar = {"axes": axes, "mask": tab.mask, "grid": tab.getGrid(), "attributes": att}
    else:
        dictvar = {"axes": axes, "attributes": att}
    tab_out = CDMS2createVariable(tab_out, **dictvar)
    return tab_out


def OperationSubtract(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Subtracts every elements of 'tab' by 'number_or_tab'
    #################################################################################

    for more information:
    import MV2
    help(MV2.subtract)
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2subtract(tab, number_or_tab)


# Dictionary of operations
dict_operations = {'divide': OperationDivide, 'minus': OperationSubtract, 'multiply': OperationMultiply,
                   'plus': OperationAdd}


def RmsAxis(tab, ref, weights=None, axis=0, centered=0, biased=1):
    """
    #################################################################################
    Description:
    genutil.rms applied on two arrays

    Computes the root mean square difference between tab and ref on the given axis

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the modeled variable
    :param ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the observed variable
    :param weights: masked_array, optional
        weights applied to each grid point
        default value = None returns equally weighted statistic
        If you want to compute the weighted statistic, provide weights here
    :param axis: int or string, optional
        name ('t', 'x', 'y', 'xy',...) or number (0, 1, 2,...) of the axis over which you want to compute the rms
        default value = 0 returns statistic computed over the first axis
    :param centered: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :return rmse: float
        value of root mean square difference
    """
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, weights=weights, axis=axis, centered=centered, biased=biased)
    except:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": RMS over axis " + str(axis),
            str().ljust(5) + "cannot perform RMS along given axis",
            str().ljust(10) + "axes may not be in the same order in 'ref' and 'tab'",
            str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
            str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    return float(rmse)


def RmsHorizontal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    genutil.rms applied on two masked_arrays that are on the same grid

    Computes the root mean square difference between tab and ref on horizontal axes (lat and lon)

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the modeled variable
    :param ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the observed variable
    :param centered: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :return rmse: float
        value of root mean square difference
    """
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, weights='weighted', axis='xy', centered=centered, biased=biased)
    except:
        lat_num = get_num_axis(tab, 'latitude')
        lon_num = get_num_axis(tab, 'longitude')
        try:
            rmse = GENUTILrms(tab, ref, weights='weighted', axis=str(lat_num)+str(lon_num), centered=centered,
                              biased=biased)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal RMS",
                str().ljust(5) + "cannot perform horizontal RMS",
                str().ljust(10) + "either lat and lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat and lon are not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    return float(rmse)


def RmsMeridional(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    genutil.rms applied on two masked_arrays that are on the same grid

    Computes the root mean square difference between tab and ref on meridional axis (lat)

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the modeled variable
    :param ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the observed variable
    :param centered: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :return rmse: float
        value of root mean square difference
    """
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis='y', centered=centered, biased=biased)
    except:
        lat_num = get_num_axis(tab, 'latitude')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lat_num), centered=centered, biased=biased)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional RMS",
                str().ljust(5) + "cannot perform meridional RMS",
                str().ljust(10) + "lat cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    return float(rmse)


def RmsTemporal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    genutil.rms applied on two masked_arrays that are horizontally averaged

    Computes the root mean square difference between tab and ref along time axis

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the modeled variable
    :param ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the observed variable
    :param centered: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :return rmse: float
        value of root mean square difference
    """
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis='t', centered=centered, biased=biased)
    except:
        time_num = get_num_axis(tab, 'time')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(time_num), centered=centered, biased=biased)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal RMS",
                str().ljust(5) + "cannot perform temporal RMS",
                str().ljust(10) + "time cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or time is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    return float(rmse)


def RmsZonal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    genutil.rms applied on two masked_arrays that are on the same grid

    Computes the root mean square difference between tab and ref on zonal axis (lon)

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the modeled variable
    :param ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
        usually it is the observed variable
    :param centered: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :return rmse: float
        value of root mean square difference
    """
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis='x', centered=centered, biased=biased)
    except:
        lon_num = get_num_axis(tab, 'longitude')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lon_num), centered=centered, biased=biased)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": zonal RMS",
                str().ljust(5) + "cannot perform zonal RMS",
                str().ljust(10) + "lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lon is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList()),
            ]
            EnsoErrorsWarnings.my_error(list_strings)
        # Computes the root mean square difference
        try:
            rmse = GENUTILrms(tab, ref, axis='t', centered=centered, biased=1)
        except:
            time_num = get_num_axis(tab, 'time')
            try:
                rmse = GENUTILrms(tab, ref, axis=str(time_num), centered=centered, biased=1)
            except:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal RMS",
                    str().ljust(5) + "cannot perform temporal RMS",
                    str().ljust(10) + "time cannot be found in 'ref' / 'tab'",
                    str().ljust(10) + "or time is not in the same order in 'ref' and 'tab'",
                    str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                    str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
                ]
                EnsoErrorsWarnings.my_error(list_strings)
    return float(rmse)


# Dictionary of RMS methods
dict_rms = {'axis': RmsAxis, 'horizontal': RmsHorizontal, 'meridional': RmsMeridional, 'time': RmsTemporal,
            'zonal': RmsZonal}


def Std(tab, weights=None, axis=0, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes standard deviation
    #################################################################################

    for more information:
    import genutil
    help(genutil.statistics.std)
    """
    try:
        grid = tab.getGrid()
    except:
        pass
    tmp = GENUTILstd(tab, weights=weights, axis=axis, centered=centered, biased=biased)
    try:
        tmp.setGrid(grid)
    except:
        pass
    return tmp


def SumAxis(tab, axis=None, fill_value=0, dtype=None):
    """
    #################################################################################
    Description:
    MV2.sum applied on tab

    Sum of elements along a certain axis using fill_value for missing

    Uses CDAT
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param axis: int or string, optional
        name ('t', 'x', 'y', 'xy',...) or number (0, 1, 2,...) of the axis over which you want to compute the rms
        default value = 0 returns statistic computed over the first axis
    :param fill_value: int or float, optional
        fill_value used for missing (masked) data
    :param dtype : None or dtype, optional
        Determines the type of the returned array and of the accumulator where the elements are summed. If dtype has the
        value None and the type of tab is an integer type of precision less than the default platform integer, then the
        default platform integer precision is used. Otherwise, the dtype is the same as that of tab
    :return sum_along_axis: MaskedArray or scalar
        An array with the same shape as tab, with the specified axis removed. If tab is a 0-d array, or if axis is None,
        a scalar is returned.
    """
    try:
        sum_along_axis = MV2sum(tab, axis=axis, fill_value=fill_value, dtype=dtype)
    except:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": sum over axis " + str(axis),
            str().ljust(5) + "cannot perform sum along given axis",
            str().ljust(10) + "axes: " + str(tab.getAxisList()),
            str().ljust(15) + "axis = " + str(axis) + " ; fill_value = " + str(fill_value) + " ; dtype = " + str(dtype),
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    return sum_along_axis


def TimeBounds(tab):
    """
    #################################################################################
    Description:
    Finds first and last dates of tab's time axis, tab must be a uvcdat masked_array
    #################################################################################

    Returns a tuple of strings: e.g., ('1979-1-1 11:59:60.0', '2016-12-31 11:59:60.0')
    """
    time = tab.getTime().asComponentTime()
    return (str(time[0]), str(time[-1]))
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of more complex functions (based on uvcdat) used in EnsoMetricsLib.py
#
def annualcycle(tab):
    """
    #################################################################################
    Description:
    Computes the annual cycle (climatological value of each calendar month) of tab
    #################################################################################

    :param tab: masked_array
    :return: tab: array
        array of the monthly annual cycle
    """
    initorder = tab.getOrder()
    tab = tab.reorder('t...')
    axes = tab.getAxisList()
    time_ax = tab.getTime().asComponentTime()
    months = MV2array(list(tt.month for tt in time_ax))
    cyc = []
    for ii in range(12):
        ids = MV2compress(months == (ii + 1), range(len(tab)))
        tmp = MV2take(tab, ids, axis=0)
        # tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = MV2average(tmp, axis=0)
        cyc.append(tmp)
        del tmp
    time = CDMS2createAxis(range(12), id='time')
    moy = CDMS2createVariable(MV2array(cyc), axes=[time] + axes[1:], grid=tab.getGrid(), attributes=tab.attributes)
    moy = moy.reorder(initorder)
    time = CDMS2createAxis(range(12), id='months')
    moy.setAxis(get_num_axis(moy, 'time'), time)
    return moy


def ApplyLandmask(tab, landmask, maskland=True, maskocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if maskland is True, mask where landmask==100
        if maskocean is True, mask where landmask==0
    #################################################################################

    :param tab: masked_array
    :param landmask: masked_array
    :param maskland: boolean, optional
        masks land points
        default value is True
    :param maskocean: boolean, optional
        masks ocean points
        default value is False

    :return: tab: masked_array
        masked_array where land points and/or ocean points are masked
    """
    keyerror = None
    if maskland is True or maskocean is True:
        if tab.getGrid().shape != landmask.getGrid().shape:
            keyerror = "tab (" + str(tab.getGrid().shape) + ") and landmask (" + str(landmask.getGrid().shape) +\
                       ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": applying landmask",
                str().ljust(5) + keyerror,
                str().ljust(5) + "cannot apply landmask",
                str().ljust(5) + "this metric will be skipped"
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
        else:
            landmask_nd = MV2zeros(tab.shape)
            if landmask_nd.shape == landmask.shape:
                landmask_nd = deepcopy(landmask)
            else:
                try:
                    landmask_nd[:] = landmask
                except:
                    try:
                        landmask_nd[:, :] = landmask
                    except:
                        keyerror = "ApplyLandmask: tab must be more than 4D and this is not taken into account yet (" +\
                                   str(tab.shape) + ") and landmask (" + str(landmask.shape) + ")"
                        list_strings = [
                            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": landmask shape",
                            str().ljust(5) + keyerror,
                            str().ljust(5) + "cannot reshape landmask"
                        ]
                        EnsoErrorsWarnings.my_warning(list_strings)
            if keyerror is None:
                tab = MV2masked_where(landmask_nd.mask, tab)
                # if land = 100 instead of 1, divides landmask by 100
                if MV2minimum(landmask_nd) == 0 and MV2maximum(landmask_nd) == 100:
                    landmask_nd = landmask_nd / 100.
                if maskland is True:
                    tab = MV2masked_where(landmask_nd == 1, tab)
                if maskocean is True:
                    tab = MV2masked_where(landmask_nd == 0, tab)
    return tab, keyerror


def ApplyLandmaskToArea(area, landmask, maskland=True, maskocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if maskland is True, mask where landmask==1 and area=area*(1-landmask) (to weight island and coastal points)
        if maskocean is True, mask where landmask==0 and area=area*landmask (to weight island and coastal points)
    #################################################################################

    :param area: masked_array
        areacell
    :param landmask: masked_array
    :param maskland: boolean, optional
        masks land points and weights island and coastal points
        default value is True
    :param maskocean: boolean, optional
        masks ocean points and weights island and coastal points
        default value is False

    :return: tab: masked_array
        masked_array where land points and/or ocean points are masked
    """
    keyerror = None
    if maskland is True or maskocean is True:
        if area.getGrid().shape != landmask.getGrid().shape:
            keyerror = "ApplyLandmaskToArea: area (" + str(area.getGrid().shape) + ") and landmask (" +\
                       str(landmask.getGrid().shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": applying landmask to areacell",
                str().ljust(5) + keyerror,
                str().ljust(5) + "cannot apply landmask to areacell"
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
        if keyerror is None:
            # if land = 100 instead of 1, divides landmask by 100
            if MV2minimum(landmask) == 0 and MV2maximum(landmask) == 100:
                landmask = landmask / 100.
            area = MV2masked_where(landmask.mask, area)
            if maskland:
                area = MV2masked_where(landmask == 1, area)
                area = MV2multiply(area, 1-landmask)
            if maskocean:
                area = MV2masked_where(landmask == 0, area)
                area = MV2multiply(area, landmask)
    return area, keyerror


def ArrayListAx(tab, list1, ax_name_ax='', ax_long_name='', ax_ref=''):
    tab_out = MV2array(tab)
    ax = CDMS2createAxis(range(len(list1)), id=ax_name_ax)
    ax.regions = str(list1)
    if len(ax_long_name) > 0:
        ax.long_name = ax_long_name
    if len(ax_ref) > 0:
        ax.reference = ax_ref
    tab_out.setAxis(0, ax)
    return tab_out


def ArrayToList(tab):
    try:
        len(tab.mask)
    except:
        tmp_mask = [tab.mask]
    else:
        tmp_mask = tab.mask
    if all(ii is False for ii in tmp_mask) is True or all(ii == False for ii in tmp_mask) == True:
        tmp = NParray(tab)
    else:
        tmp = NParray(MV2where(tab.mask, 1e20, tab))
    if len(tab.shape) == 1:
        tab_out = list(tmp)
    elif len(tab.shape) == 2:
        tab_out = [list(tmp[ii]) for ii in range(len(tab))]
    else:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": bad shape",
            str().ljust(5) + "cannot transform this array to a list",
            str().ljust(10) + "the length (" + str(len(tab.shape)) + ") of the shape (" + str(tab.shape) +
            ") is too large",
            str().ljust(10) + "it is not programed yet",
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    return tab_out


def BasinMask(tab_in, region_mask, box=None, lat1=None, lat2=None, latkey='', lon1=None, lon2=None, lonkey='',
              debug=False):
    keyerror = None
    keys = ['between', 'outside']
    # temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # open file
    this_dir, this_filename = OSpath__split(__file__)
    # check basin file
    basin_generic_ncfile = OSpath__join(this_dir, '../share/EnsoMetrics/basin_generic_1x1deg.nc') 
    if not OSpath__isfile(basin_generic_ncfile):
        basin_generic_ncfile = OSpath__join(SYS_prefix, 'share', 'EnsoMetrics', 'basin_generic_1x1deg.nc')
    if debug is True:
        dict_debug = {'line1': '(path) ' + str(this_dir), 'line2': '(file) ' + str(this_filename),
                      'line3': '(basin) ' + str(basin_generic_ncfile)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'OSpath__split', 20, **dict_debug)
    ff = CDMS2open(basin_generic_ncfile)
    # read basins
    if box is not None:
        region_ref = ReferenceRegions(box)
        basin = ff('basin', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    else:
        basin = ff('basin')
    if debug is True:
        dict_debug = {'axes1': str([ax.id for ax in basin.getAxisList()]), 'shape1': str(basin.shape),
                      'line1': 'order = ' + str(basin.getOrder())}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'in BasinMask', 20, **dict_debug)
    # choose basin
    keybasin = {'atlantic': 1, 'pacific': 2, 'indian': 3, 'antarctic': 10, 'arctic': 11}
    mask = MV2zeros(basin.shape)
    if region_mask.lower() not in keybasin.keys():
        keyerror = "unknown region: " + region_mask + " (basin_generic_1x1deg.nc, regrided file from NOAA NODC" + \
                   "WOA09 Masks basin Data Files)"
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": region",
                        str().ljust(5) + keyerror,
                        str().ljust(5) + "https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA09/.Masks/.basin/"
                        + "datafiles.html"]
        EnsoErrorsWarnings.my_warning(list_strings)
    else:
        mask = MV2where(basin == keybasin[region_mask], 1, mask)
        mask = MV2where(basin.mask, 0, mask)
    # basin mask is selected only between or outside lat1 and lat2
    if latkey in keys and lat1 is not None and lat2 is not None:
        lat2d = MV2zeros(basin.shape)
        lat2d = lat2d.reorder('10')
        lat2d[:] = basin.getLatitude()
        lat2d = lat2d.reorder('10')
        tmp = MV2where(lat2d > lat1, 1, 0) + MV2where(lat2d < lat2, 1, 0)
        if latkey == 'between':
            mask = MV2where(tmp != 2, 0, mask)
        else:
            mask = MV2where(tmp == 2, 0, mask)
    # basin mask is selected only between or outside lon1 and lon2
    if lonkey in keys and lat1 is not None and lat2 is not None:
        lon2d = MV2zeros(basin.shape)
        lon2d[:] = basin.getLongitude()
        tmp = MV2where(lon2d > lon1, 1, 0) + MV2where(lon2d < lon2, 1, 0)
        if latkey == 'between':
            mask = MV2where(tmp != 2, 0, mask)
        else:
            mask = MV2where(tmp == 2, 0, mask)
    # apply mask
    tab_out = MV2masked_where(mask == 1, tab_in)
    tab_out = CDMS2createVariable(tab_out, axes=tab_in.getAxisList(), grid=tab_in.getGrid(), mask=tab_in.mask,
                                  attributes=tab_in.attributes, id=tab_in.id)
    return tab_out, keyerror


def CheckTime(tab1, tab2, frequency='monthly', min_time_steps=None, metric_name='', debug=False, **kwargs):
    """
    #################################################################################
    Description:
    Checks if tab1 and tab2 cover the same time period and adjust if not
    Checks if the minimum_length of the time period criterion if fulfilled
    #################################################################################

    :param tab1: masked_array
    :param tab2: masked_array
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
    :param metric_name: string, optional
        name of the metric calling the function
    :return:
    """
    if debug is True:
        # dict_debug = {'shape1': 'tab1.shape = ' + str(tab1.shape), 'shape2': 'tab2.shape = ' + str(tab2.shape),
        #               'time1': 'tab1.time = ' + str(TimeBounds(tab1)),
        #               'time2': 'tab2.time = ' + str(TimeBounds(tab2))}
        dict_debug = {'shape1': 'tab1.shape = ' + str(tab1.shape), 'shape2': 'tab2.shape = ' + str(tab2.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'in CheckTime (input)', 20, **dict_debug)
    # gets dates of the first and last the time steps of tab1
    stime1 = tab1.getTime().asComponentTime()[0]
    etime1 = tab1.getTime().asComponentTime()[-1]

    # gets dates of the first and last the time steps of tab2
    stime2 = tab2.getTime().asComponentTime()[0]
    etime2 = tab2.getTime().asComponentTime()[-1]

    # retains only the latest start date and the earliest end date
    if stime1.year > stime2.year:
        stime = stime1
    elif stime1.year < stime2.year:
        stime = stime2
    else:
        if stime1.month > stime2.month:
            stime = stime1
        elif stime1.month < stime2.month:
            stime = stime2
        else:
            if stime1.day > stime2.day:
                stime = stime1
            elif stime1.day < stime2.day:
                stime = stime2
            else:
                stime = max(stime1, stime2)
    if etime1.year < etime2.year:
        etime = etime1
    elif etime1.year > etime2.year:
        etime = etime2
    else:
        if etime1.month < etime2.month:
            etime = etime1
        elif etime1.month > etime2.month:
            etime = etime2
        else:
            if etime1.day < etime2.day:
                etime = etime1
            elif etime1.day > etime2.day:
                etime = etime2
            else:
                etime = min(etime1, etime2)
    # stime = max(stime1, stime2)
    # etime = min(etime1, etime2)

    # defines the period between the two dates
    if frequency == 'daily':
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime.day, 23, 59, 0)
    elif frequency == 'monthly':
        etime_day = monthrange(etime.year, etime.month)[-1]
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime_day, 23, 59, 0)
    elif frequency == 'yearly':
        stime_adjust = CDTIMEcomptime(stime.year, 1, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, 12, 31, 23, 59, 0)
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())

    # retains only the time-period common to both tab1 and tab2
    tab1_sliced = tab1(time=(stime_adjust, etime_adjust))
    tab2_sliced = tab2(time=(stime_adjust, etime_adjust))
    if debug is True:
        # dict_debug = {'shape1': 'tab1.shape = ' + str(tab1_sliced.shape),
        #               'shape2': 'tab2.shape = ' + str(tab2_sliced.shape),
        #               'time1': 'tab1.time = ' + str(TimeBounds(tab1_sliced)),
        #               'time2': 'tab1.time = ' + str(tab1_sliced.getTime().asComponentTime()[:]),
        #               'time3': 'tab2.time = ' + str(TimeBounds(tab2_sliced)),
        #               'time4': 'tab2.time = ' + str(tab2_sliced.getTime().asComponentTime()[:])}
        dict_debug = {'shape1': 'tab1.shape = ' + str(tab1_sliced.shape),
                      'shape2': 'tab2.shape = ' + str(tab2_sliced.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'in CheckTime (output)', 20, **dict_debug)
    if len(tab1_sliced.getTime()[:]) != len(tab2_sliced.getTime()[:]):
        keyerror1 = "missing time step within the given period"
    else:
        keyerror1 = None
        tab2_sliced.setAxis(0, tab1_sliced.getTime())

    # checks if the remaining time-period fulfills the minimum length criterion
    if min_time_steps is not None:
        if len(tab1_sliced) < min_time_steps or len(tab2_sliced) < min_time_steps:
            shortest = min(len(tab1_sliced), len(tab2_sliced))
            EnsoErrorsWarnings.too_short_time_period(metric_name, shortest, min_time_steps, INSPECTstack())
            keyerror2 = "too short time period (variable1:" + str(len(tab1_sliced)) + " ; variable2:" +\
                        str(len(tab2_sliced)) + ")"
        else:
            keyerror2 = None
    else:
        keyerror2 = None

    # errors
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
    return tab1_sliced, tab2_sliced, keyerror


def CheckUnits(tab, var_name, name_in_file, units, return_tab_only=True, **kwargs):
    """
    #################################################################################
    Description:
    Checks the units of the variable and changes it if necessary
    Works for current/wind velocities, heat flux, precipitation, pressure, temperature, wind stress

    Uses MV2 (uvcdat) to find the minimum value, to multiply and to subtract
    #################################################################################

    :param tab: array
        array containing 'var_name'
    :param var_name: string
        name of the variable included in 'tab'
    :param name_in_file: string
        name of the variable in the file (usually the short_name)
    :param units: string
        units of the variable included in 'tab'
    :param return_tab_only: boolean, optional
        default value = True, only the tab is returned
        True if you want only the tab, if you want the new units also pass anything but true
    :return tab: array
        array with new units (if applicable)
    """
    keyerror = None
    if var_name in ['temperature']:
        if units in ['K', 'Kelvin', 'Kelvins', 'degree K', 'degree Kelvin', 'degree Kelvins', 'degree_K',
                     'degree_Kelvin', 'degree_Kelvins', 'degreeK', 'degreeKelvin', 'degreeKelvins', 'degrees K',
                     'degrees Kelvin', 'degrees Kelvins', 'degrees_K', 'degrees_Kelvin', 'degrees_Kelvins', 'degreesK',
                     'degreesKelvin', 'degreesKelvins', 'deg K', 'deg Kelvin', 'deg Kelvins', 'deg_K', 'deg_Kelvin',
                     'deg_Kelvins', 'degK', 'degKelvin', 'degKelvins', 'deg. K', 'deg. Kelvin', 'deg. Kelvins']:
            # check if the temperature units is really K
            if float(MV2minimum(tab)) > 150:
                # unit change of the temperature: from K to degC
                tab = dict_operations['minus'](tab, 273.15)
            else:
                minmax = [MV2minimum(tab), MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, units, minmax, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        elif units in ['C', 'celsius', 'Celsius', 'degree C', 'degree celsius', 'degree Celsius', 'degree_C',
                       'degree_celsius', 'degree_Celsius', 'degreeC', 'degreecelsius', 'degreeCelsius', 'degrees C',
                       'degrees celsius', 'degrees Celsius', 'degrees_C', 'degrees_celsius', 'degrees_Celsius',
                       'degreesC', 'degreescelsius', 'degreesCelsius', 'deg C', 'deg celsius', 'deg Celsius', 'deg_C',
                       'deg_celsius', 'deg_Celsius', 'degC', 'degcelsius', 'degCelsius', 'deg. C', 'deg. celsius',
                       'deg. Celsius']:
            # check if the temperature units is really degC
            if float(MV2minimum(tab)) > 50:
                minmax = [MV2minimum(tab), MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, minmax, units, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "degC"
    elif var_name in ['precipitations']:
        if units in ["kg/m2/s", "kg m-2 s-1", "kg/m^2/s", "kg/m**2/s", "kg m**-2 s**-1", "Kg/m2/s", "Kg m-2 s-1",
                     "Kg/m^2/s", "Kg/m**2/s", "Kg m**-2 s**-1"]:
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = dict_operations['multiply'](tab, 86400)
        elif units in ["mm/day", "mm day-1", 'mm day**-1', "mm/d", "mm d-1", 'mm d**-1']:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "mm/day"
    elif var_name in ['wind stress']:
        if units not in ['N/m2', 'N m-2', 'N/m^2', 'N/m**2', 'N m**-2', 'Pa', 'pascal', 'pascals', 'Pascal', 'Pascals']:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "N/m2"
    elif var_name in ['velocity']:
        if units in ['cm/s', 'cm s-1', 'cm s**-1', 'cm/sec', 'cm sec-1', 'cm sec**-1']:
            # unit change of the velocity: from cm/s to m/s
            tab = dict_operations['multiply'](tab, 1e-2)
        elif units in ['m/s', 'm s-1', 'm s**-1', 'm/sec', 'm sec-1', 'm sec**-1']:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m/s"
    elif var_name in ['heat flux']:
        if units in ['W/m2', 'W m-2', 'W/m^2', 'W/m**2', 'W m**-2']:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "W/m2"
    elif var_name in ['pressure']:
        if units in ['N/m2', 'N m-2', 'N/m^2', 'N/m**2', 'N m**-2', 'Pa', 'pascal', 'pascals', 'Pascal', 'Pascals']:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "Pa"
    elif var_name in ['sea surface height']:
        if units in ['cm', 'centimeter']:
            # unit change of the sea surface height: from cm to m
            tab = dict_operations['multiply'](tab, 1e-2)
        elif units in ['m', 'meter']:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m"
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.my_warning(list_strings)
    if return_tab_only is True:
        return tab
    else:
        return tab, units, keyerror


def Event_selection(tab, frequency, nbr_years_window=None, list_event_years=[]):
    if frequency not in ['daily', 'monthly', 'yearly']:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    if len(list_event_years) == 0:
        tax = tab.getTime().asComponentTime()
        list_event_years = sorted(list(set([tax[ii].year for ii in range(len(tax))])))
    else:
        list_event_years = sorted(list_event_years)
    # function to fill array with masked value where the data is not available
    def fill_array(tab, units, freq):
        y1 = tab.getTime().asComponentTime()[0].year
        m1 = tab.getTime().asComponentTime()[0].month
        d1 = tab.getTime().asComponentTime()[0].day
        if len(tab.shape) == 1:
            tab_out = MV2zeros(nbr_years_window * 12)
        elif len(tab.shape) == 2:
            tab_out = MV2zeros((nbr_years_window * 12, tab.shape[1]))
        else:
            tab_out = MV2zeros((nbr_years_window * 12, tab.shape[1], tab.shape[2]))
        tab_out = MV2masked_where(tab_out == 0, tab_out)
        axis = CDMS2createAxis(range(len(tab_out)), id='time')
        axis.units = units
        tab_out.setAxis(0, axis)
        for ii in range(len(tab)):
            y2 = tab.getTime().asComponentTime()[ii].year
            m2 = tab.getTime().asComponentTime()[ii].month
            d2 = tab.getTime().asComponentTime()[ii].day
            if freq == 'yearly':
                if y2 == y1:
                    tab_out[ii:ii + len(tab)] = deepcopy(tab)
                    break
            elif freq == 'monthly':
                if y2 == y1 and m2 == m1:
                    tab_out[ii:ii + len(tab)] = deepcopy(tab)
                    break
            elif freq == 'daily':
                if y2 == y1 and m2 == m1 and d2 == d1:
                    tab_out[ii:ii + len(tab)] = deepcopy(tab)
                    break
        return tab_out
    # compute composite
    if nbr_years_window is not None:
        composite = list()
        for yy in list_event_years:
            # first and last years of the window
            yy1, yy2 = yy - ((nbr_years_window // 2) - 1), yy + (nbr_years_window // 2)
            # create time bounds from 'first and last years of the window'
            timebnds = (str(yy1) + '-01-01 00:00:00.0', str(yy2) + '-12-31  23:59:60.0')
            # select the right time period in the given tab
            tmp1 = tab(time=timebnds)
            # sometimes there is some errors with 'time=timebnds'
            # if the time slice selected has the right length: do nothing
            # else: fill the beginning / end of the time series by masked values (done by the function 'fill_array')
            if frequency == 'yearly':
                length = nbr_years_window
                units = 'years since ' + timebnds[0]
                units_out = 'years since 0001-07-02 12:00:00'
            elif frequency == 'monthly':
                length = nbr_years_window * 12
                units = 'months since ' + timebnds[0]
                units_out = 'months since 0001-01-15 12:00:00'
            elif frequency == 'daily':
                date1 = date(yy1, 01, 01)
                date2 = date(yy2, 12, 31)
                length = (date2 - date1).days
                units = 'days since ' + timebnds[0]
                units_out = 'days since 0001-01-01 12:00:00'
            if len(tmp1) == length:
                tmp2 = deepcopy(tmp1)
            else:
                tmp2 = fill_array(tmp1, units, frequency)
            # save the selected time slice
            composite.append(tmp2)
        composite = MV2array(composite)
        # axis list
        axis0 = CDMS2createAxis(MV2array(list_event_years, dtype='int32'), id='years')
        axis1 = CDMS2createAxis(range(len(composite[0])), id='months')
        axis1.units = units_out
        axes = [axis0, axis1]
        if tab.shape > 1:
            axes = axes + tab.getAxisList()[1:]
        composite.setAxisList(axes)
    else:
        time_ax = tab.getTime().asComponentTime()  # gets component time of tab
        list_years = [yy.year for yy in time_ax[:]]  # listing years in tab (from component time)
        indices = MV2arange(tab.size)
        # creates a tab of 'condition' where True is set when the event is found, False otherwise
        try:
            condition = [True if yy in list_event_years else False for yy in list_years]
        except:
            list_event_years = [str(yy) for yy in list_event_years]
            condition = [True if str(yy) in list_event_years else False for yy in list_years]
        ids = MV2compress(condition, indices)  # gets indices of events
        composite = MV2take(tab, ids, axis=0)  # gets events
        axis0 = CDMS2createAxis(MV2array(list_event_years, dtype='int32'), id='years')
        composite.setAxis(0, axis0)
    return composite


def Composite(tab, list_event_years, frequency, nbr_years_window=None):
    return MV2average(
        Event_selection(tab, frequency, nbr_years_window=nbr_years_window, list_event_years=list_event_years), axis=0)


def DetectEvents(tab, season, threshold, normalization=False, nino=True, compute_season=True):
    """
    #################################################################################
    Description:
    Detects Nina or Nino events
    These events are detected when 'tab' anomalies during 'season' are above (less) then 'threshold'
    The anomalies can be normalized

    Uses MV2 (uvcdat) to create an empty array, to create an array of indices, to define conditions, to select the
    indices depending on the conditions and to select the years of the events
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param season: string
        one month (e.g, 'DEC'), two months (e.g., 'DJ'), three months (e.g., 'NDJ'), four months (e.g., 'NDJF'), period
        when the events are detected
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param normalization: boolean, optional
        True if events are detected based on the standard deviation, if not pass anything but True
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :return list_of_years: list
        list of years including a detected event
    """
    # Seasonal mean and anomalies
    if compute_season is True:
        tab = SeasonalMean(tab, season, compute_anom=True)
    # Normalization ?
    if normalization is True:
        threshold = threshold * float(GENUTILstd(tab, weights=None, axis=0, centered=1, biased=1))
    # Initialization
    tab_threshold = MV2zeros(tab.shape)
    tab_threshold.fill(threshold)
    list_years = [tab.getTime().asComponentTime()[yy].year for yy in range(len(tab))]
    indices = MV2arange(len(tab))
    # Conditions
    if nino is True:
        condition = MV2where(tab > tab_threshold, True, False)
    else:
        condition = MV2where(tab < tab_threshold, True, False)
    # Indices of the events
    ids = MV2compress(condition, indices)
    # Events years
    return list(MV2take(list_years, ids, axis=0))


def Detrend(tab, info, axis=0, method='linear', bp=0):
    """
    #################################################################################
    Description:
    Removes trend along 'axis' from 'tab'
    #################################################################################

    :param tab: array
        tab of data to detrend
    :param info: string
        information about what is done to 'tab'
    :param axis: int, optional
        axis along which to detrend the data
        default value is the first axis (0)
    :param method: string, optional
        detrending method:
        'constant': only the mean of 'tab' is subtracted
        'linear':   the result of a linear least-squares fit to 'tab' is subtracted from 'tab'
    :param bp: array of integer, optional
        a sequence of break points. If given, an individual linear fit is performed for each part of 'tab' between two
        break points
        break points are specified as indices into 'tab'
    :return new_tab: array
        detrended data
    """
    if method not in ['linear', 'constant']:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": method",
            str().ljust(5) + "unknown method: " + str(method)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if method in ['linear', 'constant']:
        axes = tab.getAxisList()
        grid = tab.getGrid()
        mask = tab.mask
        mean = AverageTemporal(tab)
        new_tab = MV2array(SCIPYsignal_detrend(tab, axis=axis, type=method, bp=bp))
        new_tab = new_tab + mean
        new_tab = MV2masked_where(mask, new_tab)
        new_tab.setAxisList(axes)
        new_tab.setGrid(grid)
        if method == 'linear':
            info = info + ', time series are linearly detrended'
        else:
            info = info + ', the mean value of the time series is subtracted'
    return new_tab, info


def DurationAllEvent(tab, threshold, nino=True, debug=False):
    """
    #################################################################################
    Description:
    Duration of Nina or Nino events
    The duration is the number of consecutive timestep when tab < threshold for La Nina and tab > threshold for El Nino

    Uses CDMS2 (uvcdat) for axes
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return list_of_years: list
        list of years including a detected event
    """
    tmp = MV2array([DurationEvent(tab[tt], threshold, nino=nino, debug=debug) for tt in range(len(tab))])
    tmp.setAxis(0, tab.getAxis(0))
    return tmp


def DurationEvent(tab, threshold, nino=True, debug=False):
    """
    #################################################################################
    Description:
    Duration of Nina or Nino events
    The duration is the number of consecutive timestep when tab < threshold for La Nina and tab > threshold for El Nino

    Uses MV2 (uvcdat)
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return list_of_years: list
        list of years including a detected event
    """
    try:
        len(tab.mask)
    except:
        mask = [tab.mask]
    else:
        mask = tab.mask
    # if debug is True:
    #     dict_debug = {'line1': 'threshold = ' + str(threshold) + '  ;  nino = ' + str(nino)
    #                            + '  ;  len(tab) = ' + str(len(tab)),
    #                   'line2': 'tab = ' + str(tab) + '\nmask = ' + str(mask)
    #                   }
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    if all(ii is False for ii in mask) is True:
        pass
    else:
        if nino is True:
            tab = MV2where(tab.mask, -9999, tab)
        else:
            tab = MV2where(tab.mask, 9999, tab)
    # if debug is True:
    #     dict_debug = {'line1': 'after unmasking',
    #                   'line2': 'tab = ' + str(tab)}
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    tmp1 = list(reversed(tab[:len(tab)/2]))
    tmp2 = list(tab[len(tab)/2:])
    if nino is True:
        try:
            nbr_before = next(x[0] for x in enumerate(tmp1) if x[1] <= threshold)
        except:
            if all(ii == -9999 for ii in tmp1):
                nbr_before = 0
            elif all(ii > threshold for ii in tmp1):
                nbr_before = len(tmp1)
        try:
            nbr_after = next(x[0] for x in enumerate(tmp2) if x[1] <= threshold)
        except:
            if all(ii == -9999 for ii in tmp2):
                nbr_after = 0
            elif all(ii > threshold for ii in tmp2):
                nbr_after = len(tmp2)
    else:
        try:
            nbr_before = next(x[0] for x in enumerate(tmp1) if x[1] >= threshold)
        except:
            if all(ii == 9999 for ii in tmp1):
                nbr_before = 0
            elif all(ii < threshold for ii in tmp1):
                nbr_before = len(tmp1)
        try:
            nbr_after = next(x[0] for x in enumerate(tmp2) if x[1] >= threshold)
        except:
            if all(ii == 9999 for ii in tmp2):
                nbr_after = 0
            elif all(ii < threshold for ii in tmp2):
                nbr_after = len(tmp2)
    duration = nbr_before + nbr_after
    # if debug is True:
    #     dict_debug = {'line1': 'duration of the event = ' + str(duration)}
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    return duration


def get_num_axis(tab, name_axis):
    """
    #################################################################################
    Description:
    Finds the number of the axis named "name_axis"
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param name_axis: string
        name of an axis
        e.g., name_axis='latitude'
    :return number: int
        position of the axis named "name_axis"
    """
    num = None
    if name_axis == 'depth':
        axis_nick = 'lev'
        axis_nicks = ['z', 'Z', 'st_ocean', 'sw_ocean']
    if name_axis == 'latitude':
        axis_nick = 'lat'
        axis_nicks = ['j', 'y', 'Y', 'yt_ocean', 'yu_ocean']
    elif name_axis == 'longitude':
        axis_nick = 'lon'
        axis_nicks = ['i', 'x', 'X', 'xt_ocean', 'xu_ocean']
    elif name_axis == 'time':
        axis_nick = 'time'
        axis_nicks = ['t', 'T']
    for nn in range(len(tab.shape)):
        if axis_nick in tab.getAxisList()[nn].id:
            num = nn
            break
    if num is None:
        for nn in range(len(tab.shape)):
            for ax in axis_nicks:
                if ax == tab.getAxisList()[nn].id:
                    num = nn
                    break
    if num is None:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
            str().ljust(5) + "cannot find axis named: " + str(name_axis),
            str().ljust(5) + "axes: " + str(tab.getAxisList()),
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    return num


def get_year_by_year(tab, frequency='monthly'):
    """
    #################################################################################
    Description:
    Reshape array to a year by year array
    #################################################################################

    :param tab: masked_array
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :return: tab: array
        array of the year by year values
    """
    tab = tab.reorder('t...')
    time_ax = tab.getTime().asComponentTime()
    myshape = [1] + [ss for ss in tab.shape[1:]]
    zeros = MV2zeros(myshape)
    zeros = MV2masked_where(zeros == 0, zeros)
    if frequency == 'daily':
        days = MV2array(list(tt.day for tt in time_ax))
        months = MV2array(list(tt.month for tt in time_ax))
        months = MV2array([(mm*100)+dd for dd, mm in zip(days, months)])
        tmm = CDMS2createAxis(range(365), id='days')
        m1 = time_ax[0].day
        m2 = time_ax[-1].day
        t2 = 365
    elif frequency == 'monthly':
        months = MV2array(list(tt.month for tt in time_ax))
        tmm = CDMS2createAxis(range(12), id='months')
        m1 = time_ax[0].month
        m2 = time_ax[-1].month
        t2 = 12
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    years = sorted(set(MV2array(list(tt.year for tt in time_ax))))
    tyy = CDMS2createAxis(MV2array(years, dtype='int32'), id='years')
    axes = [tyy] + [tmm]
    val = sorted(set(months))
    tab_out = list()
    if frequency == 'daily':
        val.remove(129)
    for ii in val:
        tmp = tab.compress(months == ii, axis=0)
        if m1 != 1 and len(tmp) != len(years):
            tmp = MV2concatenate((zeros, tmp))
        if m2 != t2 and len(tmp) != len(years):
            tmp = MV2concatenate((tmp, zeros))
        tab_out.append(tmp)
    tab_out = MV2array(tab_out)
    tab_out = MV2masked_where(tab_out == 0, tab_out)
    tab_out = tab_out.reorder('10')
    if len(tab.shape) == 1:
        tab_out = CDMS2createVariable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
    else:
        axes = axes + tab.getAxisList()[1:]
        grid = tab[0].getGrid()
        mask = tab[0].mask
        mask_out = MV2zeros(tab_out.shape)
        mask_out[:, :] = mask
        tab_out = CDMS2createVariable(tab_out, axes=axes, grid=grid, mask=mask_out, attributes=tab.attributes,
                                      id=tab.id)
    return tab_out


def MinMax(tab):
    return [float(MV2minimum(tab)), float(MV2maximum(tab))]


def MyEmpty(tab, time=True, time_id=''):
    tab_out = ArrayZeros(tab)
    if time is True:
        axis = CDMS2createAxis(MV2array(len(tab_out), dtype='int32'), id=time_id)
        axes = [axis] + tab.getAxisList()[1:]
    else:
        axes = tab.getAxisList()
    tab_out.setAxisList(axes)
    return tab_out


def Normalize(tab, frequency):
    """
    #################################################################################
    Description:
    Removes trend along 'axis' from 'tab'
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :return tab: masked_array
        normalized data
    """
    axes = tab.getAxisList()
    if frequency == 'daily':
        time_steps_per_year = 365
    elif frequency == 'monthly':
        time_steps_per_year = 12
    elif frequency == 'yearly':
        time_steps_per_year = 1
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    if len(tab) % time_steps_per_year != 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": data length",
            str().ljust(5) + "the normalization function needs only full years: " +
            str(len(tab) // time_steps_per_year) +" years + " + str(len(tab) % time_steps_per_year),
            str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " + str(time_steps_per_year) +
            "), len(dataset) = " + str(len(tab)) + ", so " + str(len(tab) / float(time_steps_per_year)) + " years",
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    # reshape tab like [yy,nb]
    new_tab = list()
    for yy in range(len(tab)/time_steps_per_year):
        new_tab.append(tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year])
    new_tab = MV2array(new_tab)
    std = MV2zeros(new_tab[0].shape)
    for dd in range(time_steps_per_year):
        std[dd] = float(GENUTILstd(new_tab[:, dd], weights=None, axis=0, centered=1, biased=1))
    tab_out = deepcopy(tab)
    for yy in range(len(tab) // time_steps_per_year):
        tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
            tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] // std
    if len(tab.shape) == 1:
        tab_out = CDMS2createVariable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
    else:
        axes = axes + tab.getAxisList()[1:]
        grid = tab.getGrid()
        mask = tab.mask
        tab_out = CDMS2createVariable(tab_out, axes=axes, grid=grid, mask=mask, attributes=tab.attributes, id=tab.id)
    return tab_out


def ReadAndSelectRegion(filename, varname, box=None, time_bounds=None, frequency=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given 'varname' from the given 'filename' and selects the given 'box'

    Uses cdms2 (uvcdat) to read 'varname' from 'filename' and cdutil (uvcdat) to select the 'box'
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param varname: string
        name of the variable to read from 'filename'
    :param box: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None

    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if box is None:  # no box given
        if time_bounds is None:  # no time period given
            # read file
            tab = fi(varname)
        else:  # time period given by the user
            # read file
            tab = fi(varname, time=time_bounds)
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        if time_bounds is None:  # no time period given
            #  read file
            tab = fi(varname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        else:
            # read file
            tab = fi(varname, time=time_bounds, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    fi.close()
    # sign correction
    try:
        att1 = tab.attributes['standard_name'].lower().replace(" ", "_")
    except:
        att1 = ''
    try:
        att2 = tab.attributes['long_name'].lower().replace(" ", "_")
    except:
        att2 = ''
    if "latent_heat" in att1 or "latent_heat" in att2 or "sensible_heat" in att1 or "sensible_heat" in att2:
        if "upward" in att1 or "upward" in att2:
            # I need to be in the ocean point of view so the heat fluxes must be downwards
            tab = -1 * tab
    if time_bounds is not None:
        # sometimes the time boundaries are wrong, even with 'time=time_bounds'
        # this section checks if one time step has not been included by error at the beginning or the end of the time
        # series
        if isinstance(time_bounds[0], basestring):
            if str(tab.getTime().asComponentTime()[0]) < time_bounds[0]:
                tab = tab[1:]
            if str(tab.getTime().asComponentTime()[-1]) > time_bounds[1]:
                tab = tab[:-1]
    time_ax = tab.getTime()
    time_units = "days since " + str(time_ax.asComponentTime()[0].year) + "-01-01 12:00:00"
    time_ax.id = "time"
    time_ax.toRelativeTime(time_units)
    tab.setAxis(0, time_ax)
    if frequency is None:  # no frequency given
        pass
    elif frequency == 'daily':
        cdutil.setTimeBoundsDaily(tab)
    elif frequency == 'monthly':
        cdutil.setTimeBoundsMonthly(tab)
    elif frequency == 'yearly':
        cdutil.setTimeBoundsYearly(tab)
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    # remove axis 'level' if its length is 1
    if tab.getLevel():
        if len(tab.getLevel()) == 1:
            tab = tab(squeeze=1)
    # HadISST has -1000 values... mask them
    if 'HadISST' in filename or 'hadisst' in filename:
        tab = MV2masked_where(tab == -1000, tab)
    return tab


def ReadAreaSelectRegion(filename, areaname='', box=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given areacell from the given 'filename' and selects the given 'box'

    Uses cdms2 (uvcdat) to read areacell from 'filename' and cdutil (uvcdat) to select the 'box'
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param areaname: string, optional
        name of areacell (areacella, areacello,...) in 'filename'
    :param box: string, optional
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return area: masked_array
        masked_array containing areacell in 'box'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if box is None:  # no box given
        # read file
        try:
            areacell = fi(areaname)
        except:
            try:
                areacell = fi('areacell')
            except:
                try:
                    areacell = fi('areacella')
                except:
                    try:
                        areacell = fi('areacello')
                    except:
                        areacell = None
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        # read file
        try:
            areacell = fi(areaname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        except:
            try:
                areacell = fi('areacell', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
            except:
                try:
                    areacell = fi('areacella', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                except:
                    try:
                        areacell = fi('areacello', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                    except:
                        areacell = None
    fi.close()
    return areacell


def ReadLandmaskSelectRegion(tab, filename, landmaskname='', box=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given landmask from the given 'filename' and selects the given 'box'

    Uses cdms2 (uvcdat) to read areacell from 'filename' and cdutil (uvcdat) to select the 'box'
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param landmaskname: string, optional
        name of landmask (sftlf, lsmask, landmask,...) in 'filename'
    :param box: string, optional
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return area: masked_array
        masked_array containing landmask in 'box'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Get landmask
    if OSpath__isfile(filename):
        # Open file and get time dimension
        fi = CDMS2open(filename)
        if box is None:  # no box given
            # read file
            try:
                landmask = fi(landmaskname)
            except:
                try:
                    landmask = fi('landmask')
                except:
                    try:
                        landmask = fi('lsmask')
                    except:
                        try:
                            landmask = fi('sftlf')
                        except:
                            landmask = None
        else:  # box given by the user
            # define box
            region_ref = ReferenceRegions(box)
            # read file
            try:
                landmask = fi(landmaskname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
            except:
                try:
                    landmask = fi('landmask', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                except:
                    try:
                        landmask = fi('lsmask', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                    except:
                        try:
                            landmask = fi('sftlf', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                        except:
                            landmask = None
        fi.close()
    else:
        # Estimate landmask
        landmask = EstimateLandmask(tab)
        if box is not None:
            # define box
            region_ref = ReferenceRegions(box)
            # subset
            landmask = landmask(latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    # Return
    return landmask


def EstimateLandmask(d):
    """
    #################################################################################
    Description:
    Estimate landmask (when landmask was not given) 
    Uses cdutil (uvcdat) to create estimated landmask for model resolution
    #################################################################################

    :param d: array (CDMS)
        model variable 

    :return landmask: masked_array
        masked_array containing landmask
    """
    print('\033[93m' + str().ljust(25) + 'NOTE: Estimated landmask applied' + '\033[0m')
    n = 1
    sft = cdutil.generateLandSeaMask(d(*(slice(0, 1),) * n)) * 100.0
    sft[:] = sft.filled(100.0)
    lmsk = sft
    lmsk.setAxis(0, d.getAxis(1))
    lmsk.setAxis(1, d.getAxis(2))
    lmsk.id = 'sftlf'
    return lmsk


def Regrid(tab_to_regrid, newgrid, missing=None, order=None, mask=None, regridder='cdms', regridTool='esmf',
           regridMethod='linear', **kwargs):
    """
    #################################################################################
    Description:
    Regrids 'tab_to_regrid' to 'togrid'
    #################################################################################

    for more information:
    import cdms2
    help(cdms2.avariable)

    :param tab_to_regrid: masked_array
        masked_array to regrid (must include a CDMS grid!)
    :param newgrid: CDMS grid
        destination grid
    :param missing: float, optional
        missing values (missing data value, if any)
    :param order: string, optional
        axis order (form "tzyx", "tyx", etc.)
    :param mask: array of booleans, optional
        mask of the new grid (either 2-D or the same shape as togrid)
    :param regridder: string, optional
        regridders (either 'cdms', 'cdmsHorizontal')
        default value is 'cdms'
    :param regridTool: string, optional
        only if regrider is set to 'cdms'
        regridding tools (either 'regrid2', 'esmf', 'libcf')
        default value is 'esmf'
    :param regridMethod: string, optional
        regridding methods depend on regridder and regridTool
        'cdms'
            regridTool='regrid2' -> 'linear'
            regridTool='esmf'    -> 'conserve', 'linear', 'patch'
            regridTool='libcf'   -> 'linear'
        'cdmsHorizontal' -> None
        default value is 'linear'

    usual kwargs:
    :param newgrid_name: string, optional
        if 'newgrid' is not defined (is string) this will be used to create a grid:
            this string must contain two keywords: the grid type and the resolution (same resolution in lon and lat)
            grid type  -> 'equalarea', 'gaussian', 'generic', 'uniform'
            resolution -> '0.25x0.25deg', '0.5x0.5deg', '1x1deg', '2x2deg'
            e.g., newgrid_name='gaussian 1x1deg'
        default value is 'generic 1x1deg'
    :param region: string, optional
        if 'newgrid' is not defined (is string) this will be used to create a grid
        name of a region, domain where the grid will be defined, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return new_tab: masked_array
        tab_to_regrid regridded on newgrid
    """
    known_args = {"newgrid_name", "region"}
    extra_args = set(kwargs) - known_args
    if extra_args:
        EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
    # test given arguments
    known_regridder = ["cdms", "cdmsHorizontal"]
    if regridder not in known_regridder:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridder",
            str().ljust(5) + "unknown regridder: " + str(regridder),
            str().ljust(10) + "known regridder: " + str(known_regridder)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    elif regridder == "cdms":
        if regridTool in ["regrid2", "libcf"]:
            list_method = [None, "linear"]
        elif regridTool == "esmf":
            list_method = [None, "conserve", "linear", "patch"]
        if (regridTool is not None) and (regridTool not in ["regrid2", "esmf", "libcf"]):
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridTool",
                str().ljust(5) + "unknown regridTool: " + str(regridTool),
                str().ljust(10) + "known regridTool: " + str(["regrid2", "esmf", "libcf"])
            ]
            EnsoErrorsWarnings.my_error(list_strings)
        elif regridMethod not in list_method:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridMethod",
                str().ljust(5) + "unknown regridMethod (" + str(regridMethod) + ") for this regridTool ("
                + str(regridTool) + ")",
                str().ljust(10) + "known regridMethod: " + str(list_method)
            ]
            EnsoErrorsWarnings.my_error(list_strings)
    # test the given 'newgrid'
    if isinstance(newgrid, basestring) or newgrid is None:
        #
        # newgrid is not a grid, so a grid will be created
        # to do this, kwargs['newgrid_name'] and kwargs['region'] must be defined
        #
        # define the grid type
        for gtype in ["equalarea", "gaussian", "generic", "uniform"]:
            if gtype in kwargs['newgrid_name']:
                GridType = gtype
                break
        try:
            GridType
        except:
            GridType = "generic"
        # define resolution (same resolution in lon and lat)
        for res in ["0.25x0.25deg", "0.5x0.5deg", "0.75x0.75deg", "1x1deg", "1.25x1.25deg", "1.5x1.5deg",
                    "1.75x1.75deg", "2x2deg", "2.25x2.25deg", "2.5x2.5deg", "2.75x2.75deg"]:
            if res in kwargs['newgrid_name']:
                if res == "0.25x0.25deg":
                    GridRes = 0.25
                elif res == "0.5x0.5deg":
                    GridRes = 0.5
                elif res == "0.75x0.75deg":
                    GridRes = 0.75
                elif res == "1x1deg":
                    GridRes = 1.
                elif res == "1.25x1.25deg":
                    GridRes = 1.25
                elif res == "1.5x1.5deg":
                    GridRes = 1.5
                elif res == "1.75x1.75deg":
                    GridRes = 1.75
                elif res == "2x2deg":
                    GridRes = 2.
                elif res == "2.25x2.25deg":
                    GridRes = 2.25
                elif res == "2.5x2.5deg":
                    GridRes = 2.5
                else:
                    GridRes = 2.75
                break
        try:
            GridRes
        except:
            GridRes = 1.
        # define bounds of 'region'
        region_ref = ReferenceRegions(kwargs["region"])
        lat1, lat2 = region_ref["latitude"][0], region_ref["latitude"][1]
        lon1, lon2 = region_ref["longitude"][0], region_ref["longitude"][1]
        # create uniform axis
        nlat = lat2 - lat1
        lat = CDMS2createUniformLatitudeAxis(lat1 + (GridRes / 2.), nlat, GridRes)
        nlon = lon2 - lon1
        lon = CDMS2createUniformLongitudeAxis(lon1 + (GridRes / 2.), nlon, GridRes)
        # create grid
        newgrid = CDMS2createRectGrid(lat, lon, "yx", type=GridType, mask=None)
        newgrid.id = kwargs["newgrid_name"]
    #
    # regrid
    #
    if regridder == "cdms":
        axis = tab_to_regrid.getAxis(0)
        idname = deepcopy(axis.id)
        if len(tab_to_regrid.shape) == 3 and (axis.id == "months" or axis.id == "years"):
            axis.id = "time"
            tab_to_regrid.setAxis(0, axis)
        new_tab = tab_to_regrid.regrid(newgrid, missing=missing, order=order, mask=mask, regridTool=regridTool,
                                       regridMethod=regridMethod)
        axis = tab_to_regrid.getAxis(0)
        axis.id = idname
        tab_to_regrid.setAxis(0, axis)
        if tab_to_regrid.getGrid().shape == newgrid.shape:
            new_tab = MV2masked_where(tab_to_regrid.mask, new_tab)
    elif regridder == "cdmsHorizontal":
        regridFCT = REGRID2horizontal__Horizontal(tab_to_regrid.getGrid(), newgrid)
        new_tab = regridFCT(tab_to_regrid)
    return new_tab


def SaveNetcdf(netcdf_name, var1=None, var1_attributes={}, var1_name='', var1_time_name=None, var2=None,
               var2_attributes={}, var2_name='', var2_time_name=None, var3=None, var3_attributes={}, var3_name='',
               var3_time_name=None, var4=None, var4_attributes={}, var4_name='', var4_time_name=None, var5=None,
               var5_attributes={}, var5_name='', var5_time_name=None, var6=None, var6_attributes={}, var6_name='',
               var6_time_name=None, var7=None, var7_attributes={}, var7_name='', var7_time_name=None, var8=None,
               var8_attributes={}, var8_name='', var8_time_name=None, var9=None, var9_attributes={}, var9_name='',
               var9_time_name=None, var10=None, var10_attributes={}, var10_name='', var10_time_name=None, var11=None,
               var11_attributes={}, var11_name='', var11_time_name=None, var12=None, var12_attributes={}, var12_name='',
               var12_time_name=None, frequency="monthly", global_attributes={}, **kwargs):
    if OSpath_isdir(ntpath.dirname(netcdf_name)) is not True:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": given path does not exist",
            str().ljust(5) + "netcdf_name = " + str(netcdf_name)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if OSpath__isfile(netcdf_name) is True:
        o = CDMS2open(netcdf_name, "a")
    else:
        o = CDMS2open(netcdf_name, "w+")
    if var1 is not None:
        if var1_name == '':
            var1_name = var1.id
        if var1_time_name is not None:
            var1 = TimeButNotTime(var1, var1_time_name, frequency)
        o.write(var1, attributes=var1_attributes, dtype="float32", id=var1_name)
    if var2 is not None:
        if var2_name == '':
            var2_name = var2.id
        if var2_time_name is not None:
            var2 = TimeButNotTime(var2, var2_time_name, frequency)
        o.write(var2, attributes=var2_attributes, dtype="float32", id=var2_name)
    if var3 is not None:
        if var3_name == '':
            var3_name = var3.id
        if var3_time_name is not None:
            var3 = TimeButNotTime(var3, var3_time_name, frequency)
        o.write(var3, attributes=var3_attributes, dtype="float32", id=var3_name)
    if var4 is not None:
        if var4_name == '':
            var4_name = var4.id
        if var4_time_name is not None:
            var4 = TimeButNotTime(var4, var4_time_name, frequency)
        o.write(var4, attributes=var4_attributes, dtype="float32", id=var4_name)
    if var5 is not None:
        if var5_name == '':
            var5_name = var5.id
        if var5_time_name is not None:
            var5 = TimeButNotTime(var5, var5_time_name, frequency)
        o.write(var5, attributes=var5_attributes, dtype="float32", id=var5_name)
    if var6 is not None:
        if var6_name == '':
            var6_name = var6.id
        if var6_time_name is not None:
            var6 = TimeButNotTime(var6, var6_time_name, frequency)
        o.write(var6, attributes=var6_attributes, dtype="float32", id=var6_name)
    if var7 is not None:
        if var7_name == '':
            var7_name = var7.id
        if var7_time_name is not None:
            var7 = TimeButNotTime(var7, var7_time_name, frequency)
        o.write(var7, attributes=var7_attributes, dtype="float32", id=var7_name)
    if var8 is not None:
        if var8_name == '':
            var8_name = var8.id
        if var8_time_name is not None:
            var8 = TimeButNotTime(var8, var8_time_name, frequency)
        o.write(var8, attributes=var8_attributes, dtype="float32", id=var8_name)
    if var9 is not None:
        if var9_name == '':
            var9_name = var9.id
        if var9_time_name is not None:
            var9 = TimeButNotTime(var9, var9_time_name, frequency)
        o.write(var9, attributes=var9_attributes, dtype="float32", id=var9_name)
    if var10 is not None:
        if var10_name == '':
            var10_name = var10.id
        if var10_time_name is not None:
            var10 = TimeButNotTime(var10, var10_time_name, frequency)
        o.write(var10, attributes=var10_attributes, dtype="float32", id=var10_name)
    if var11 is not None:
        if var11_name == '':
            var11_name = var11.id
        if var11_time_name is not None:
            var11 = TimeButNotTime(var11, var11_time_name, frequency)
        o.write(var11, attributes=var11_attributes, dtype="float32", id=var11_name)
    if var12 is not None:
        if var12_name == '':
            var12_name = var12.id
        if var12_time_name is not None:
            var12 = TimeButNotTime(var12, var12_time_name, frequency)
        o.write(var12, attributes=var12_attributes, dtype="float32", id=var12_name)
    my_keys = sorted([key for key in kwargs.keys() if "var" in key and str(key.replace("var", "")).isdigit() is True],
                     key=lambda v: v.upper())
    for key in my_keys:
        if kwargs[key] is not None:
            if key + "_name" not in kwargs.keys() or (key + "_name" in kwargs.keys() and kwargs[key + "_name"] == ''):
                kwargs[key + "_name"] = kwargs[key].id
            if key + "_time_name" in kwargs.keys() and kwargs[key + "_time_name"] is not None:
                kwargs[key] = TimeButNotTime(kwargs[key], kwargs[key + "_time_name"], frequency)
            if key + "_attributes" not in kwargs.keys():
                kwargs[key + "_attributes"] = {}
            o.write(kwargs[key], attributes=kwargs[key + "_attributes"], dtype="float32", id=kwargs[key + "_name"])
    for att in global_attributes.keys():
        o.__setattr__(att, global_attributes[att])
    o.close()
    return


def SkewnessTemporal(tab):
    """
    #################################################################################
    Description:
    Computes the skewness along the time axis
    #################################################################################

    :param tab: masked_array
    :return: tab: array
        array of the temporal skewness
    """
    tab = tab.reorder('t...')
    if len(tab.shape) > 4:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many dimensions",
            str().ljust(5) + "tab.shape = " + str(tab.shape)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    elif len(tab.shape) == 1:
        skew = float(SCIPYstats__skew(tab))
    else:
        if len(tab.shape) == 2:
            skew = SCIPYstats__skew(tab, axis=0)
        else:
            # swith to numpy
            dataset = NPma__core__MaskedArray(tab)
            # masked values -> nan
            dataset = dataset.filled(fill_value=NPnan)
            # Store information about the shape/size of the input data
            time_ax = dataset.shape[0]
            spac_ax = dataset.shape[1:]
            channels = NPproduct(spac_ax)
            # Reshape to two dimensions (time, space) creating the design matrix
            dataset = dataset.reshape([time_ax, channels])
            # Find the indices of values that are not missing in one row. All the rows will have missing values in the same
            # places provided the array was centered. If it wasn't then it is possible that some missing values will be
            # missed and the singular value decomposition will produce not a number for everything.
            nonMissingIndex = NPwhere(NPisnan(dataset[0]) == False)[0]
            # Remove missing values from the design matrix.
            dataNoMissing = dataset[:, nonMissingIndex]
            new_dataset = SCIPYstats__skew(dataNoMissing, axis=0)
            flatE = NPones([channels], dtype=dataset.dtype) * NPnan
            flatE = flatE.astype(dataset.dtype)
            flatE[nonMissingIndex] = new_dataset
            skew = flatE.reshape(spac_ax)
            skew = MV2masked_where(NPisnan(skew), skew)
        skew = CDMS2createVariable(MV2array(skew), axes=tab.getAxisList()[1:], grid=tab.getGrid(), mask=tab[0].mask,
                                   attributes=tab.attributes, id='skewness')
    return skew


def SmoothGaussian(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using gaussian moving window average
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the triangle moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    # Reorder tab in order to put 'axis' in first position
    indices = range(len(tab.shape))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder + str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the gaussian weight array
    weightGauss = list()
    for ii in range(window):
        ii = ii - degree + 1
        frac = ii / float(window)
        gauss = float(1 / (NPexp((4 * frac) ** 2)))
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(gauss)
        weightGauss.append(ww)
    weightGauss = MV2array(weightGauss)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    sum_weight = MV2sum(weightGauss, axis=0)
    for ii in range(len(smoothed_tab)):
        smoothed_tab[ii] = MV2sum(MV2array(new_tab[ii:ii + window]) * weightGauss * weightGauss, axis=0) / sum_weight

    # Axes list
    axes0 = new_tab[(window // 2):len(new_tab) - (window // 2)].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


def SmoothSquare(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using square moving window average
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the square moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    # Reorder tab in order to put 'axis' in first position
    indices = range(len(tab.shape))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder+str(ii)
    new_tab = tab.reorder(newOrder)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        smoothed_tab[ii] = sum(new_tab[ii:ii + window]) / float(window)

    # Axes list
    axes0 = new_tab[(window/2):len(new_tab)-(window/2)].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


def SmoothTriangle(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using triangle moving window average
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the triangle moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    # Reorder tab in order to put 'axis' in first position
    indices = range(len(tab.shape))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder+str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the weight array (triangle)
    weight = list()
    for ii in range(0, (2 * degree)+1):
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(float(1 + degree - abs(degree - ii)))
        weight.append(ww)
    weight = MV2array(weight)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    sum_weight = MV2sum(weight, axis=0)
    for ii in range(len(smoothed_tab)):
        smoothed_tab[ii] = MV2sum(MV2array(new_tab[ii:ii + window]) * weight, axis=0) / sum_weight

    # Axes list
    axes0 = new_tab[(window // 2):len(new_tab) - (window // 2)].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)
    if tab.getGrid():
        try:
            smoothed_tab.setGrid(tab.getGrid())
        except:
            pass

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


# Dictionary of seasons
sea_dict = dict(JAN=cdutil.JAN, FEB=cdutil.FEB, MAR=cdutil.MAR, APR=cdutil.APR, MAY=cdutil.MAY, JUN=cdutil.JUN,
                JUL=cdutil.JUL, AUG=cdutil.AUG, SEP=cdutil.SEP, OCT=cdutil.OCT, NOV=cdutil.NOV, DEC=cdutil.DEC,
                JF=cdutil.times.Seasons('JF'), FM=cdutil.times.Seasons('FM'), MA=cdutil.times.Seasons('MA'),
                AM=cdutil.times.Seasons('AM'), MJ=cdutil.times.Seasons('MJ'), JJ=cdutil.times.Seasons('JJ'),
                JA=cdutil.times.Seasons('JA'), AS=cdutil.times.Seasons('AS'), SO=cdutil.times.Seasons('SO'),
                ON=cdutil.times.Seasons('ON'), ND=cdutil.times.Seasons('ND'), DJ=cdutil.times.Seasons('DJ'),
                JFM=cdutil.times.Seasons('JFM'), FMA=cdutil.times.Seasons('FMA'), MAM=cdutil.MAM,
                AMJ=cdutil.times.Seasons('AMJ'), MJJ=cdutil.times.Seasons('MJJ'), JJA=cdutil.JJA,
                JAS=cdutil.times.Seasons('JAS'), ASO=cdutil.times.Seasons('ASO'), SON=cdutil.SON,
                OND=cdutil.times.Seasons('OND'), NDJ=cdutil.times.Seasons('NDJ'), DJF=cdutil.DJF,
                JFMA=cdutil.times.Seasons('JFMA'),FMAM=cdutil.times.Seasons('FMAM'),MAMJ=cdutil.times.Seasons('MAMJ'),
                AMJJ=cdutil.times.Seasons('AMJJ'),MJJA=cdutil.times.Seasons('MJJA'),JJAS=cdutil.times.Seasons('JJAS'),
                JASO=cdutil.times.Seasons('JASO'),ASON=cdutil.times.Seasons('ASON'),SOND=cdutil.times.Seasons('SOND'),
                ONDJ=cdutil.times.Seasons('ONDJ'),NDJF=cdutil.times.Seasons('NDJF'),DJFM=cdutil.times.Seasons('DJFM'))


def SeasonalMean(tab, season, compute_anom=False):
    """
    #################################################################################
    Description:
    Creates a time series of the seasonal mean ('season') and computes the anomalies (difference from the mean value; if
    applicable)
    Improved cdutil seasonal mean (more seasons and incomplete seasons are removed)

    Uses cdutil (uvcdat) to select the 'season', to average it, and to compute the anomalies (if applicable)
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param season: string
        name of a season, must be defined in 'sea_dict'
    :param compute_anom: boolean, optional
        default value = True, computes anomalies (difference from the mean value)
        True if you want to compute anomalies, if you don't want to compute anomalies pass anything but true
    :return tab: masked_array
        time series of the seasonal mean ('season') anomalies (if applicable)
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Checks if the season has been defined
    try:
        sea_dict[season]
    except:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": season",
                        str().ljust(5) + "unknown season: " + str(season)]
        EnsoErrorsWarnings.my_error(list_strings)
    else:
        if season in ['DJ', 'NDJ', 'DJF', 'ONDJ', 'NDJF', 'NDJF']:
            # these 'seasons' are between two years
            # if I don't custom 'tab' cdutil will compute half season mean
            # (i.e., for NDJ the first element would be for J only and the last for ND only)
            time_ax_comp = tab.getTime().asComponentTime()
            ntime = len(time_ax_comp)
            ii, jj = 0, 0
            if season == 'DJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'NDJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 11: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'DJF':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 2: break
            elif season == 'ONDJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 10: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'NDJF':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 11: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 2: break
            elif season == 'DJFM':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 3: break
            tab = tab[ii:ntime - jj]
        if compute_anom:
            tab = sea_dict[season].departures(tab)  # extracts 'season' seasonal anomalies (from climatology)
        else:
            tab = sea_dict[season](tab)  # computes the 'season' climatology of a tab
    if season == 'DJF':
        time_ax = tab.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
        tab.setAxis(0, time_ax)
    return tab


# Dictionary of smoothing methods
dict_smooth = {'gaussian': SmoothGaussian, 'square': SmoothSquare, 'triangle': SmoothTriangle}


def Smoothing(tab, info, axis=0, window=5, method='triangle'):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using moving window average based on 'method'
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param info: string
        information about what was done on tab
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the moving window average
        default value is 5
    :param method: string, optional
        smoothing method:
            'gaussian': gaussian shaped window
            'square':   square shaped window
            'triangle': triangle shaped window
    :return: smoothed_tab: masked_array
        smoothed data
    """
    try: dict_smooth[method]
    except:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing method (running mean)",
            str().ljust(5) + "unkwown smoothing method: " + str(method),
            str().ljust(10) + "known smoothing method: " + str(sorted(dict_smooth.keys(), key=lambda v: v.upper())),
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    info = info + ', smoothing using a ' + str(method) + ' shaped window of ' + str(window) + ' points'
    return dict_smooth[method](tab, axis=axis, window=window), info


def SkewMonthly(tab):
    """
    #################################################################################
    Description:
    Computes the monthly standard deviation (value of each calendar month) of tab
    #################################################################################

    :param tab: masked_array
    :return: tab: array
        array of the monthly standard deviation
    """
    initorder = tab.getOrder()
    tab = tab.reorder('t...')
    axes = tab.getAxisList()
    time_ax = tab.getTime().asComponentTime()
    months = MV2array(list(tt.month for tt in time_ax))
    cyc = []
    for ii in range(12):
        tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = SCIPYstats__skew(tmp)
        cyc.append(tmp)
        del tmp
    time = CDMS2createAxis(range(12), id='time')
    skew = CDMS2createVariable(MV2array(cyc), axes=[time] + axes[1:], grid=tab.getGrid(), attributes=tab.attributes)
    skew = skew.reorder(initorder)
    time = CDMS2createAxis(range(12), id='months')
    skew.setAxis(get_num_axis(skew, 'time'), time)
    return skew


def StdMonthly(tab):
    """
    #################################################################################
    Description:
    Computes the monthly standard deviation (value of each calendar month) of tab
    #################################################################################

    :param tab: masked_array
    :return: tab: array
        array of the monthly standard deviation
    """
    initorder = tab.getOrder()
    tab = tab.reorder('t...')
    axes = tab.getAxisList()
    time_ax = tab.getTime().asComponentTime()
    months = MV2array(list(tt.month for tt in time_ax))
    cyc = []
    for ii in range(12):
        tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = Std(tmp, axis=0)
        cyc.append(tmp)
        del tmp
    time = CDMS2createAxis(range(12), id='time')
    std = CDMS2createVariable(MV2array(cyc), axes=[time] + axes[1:], grid=tab.getGrid(), attributes=tab.attributes)
    std = std.reorder(initorder)
    time = CDMS2createAxis(range(12), id='months')
    std.setAxis(get_num_axis(std, 'time'), time)
    return std


def TimeButNotTime(tab, new_time_name, frequency):
    tab_out = deepcopy(tab)
    time_num = get_num_axis(tab_out, 'time')
    timeax = tab_out.getAxis(time_num).asComponentTime()
    year1, month1, day1 = timeax[0].year, timeax[0].month, timeax[0].day
    if frequency == 'daily':
        freq = 'days'
    elif frequency == 'monthly':
        freq = 'months'
    elif frequency == 'yearly':
        freq = 'years'
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    axis = CDMS2createAxis(range(len(tab_out)), id=new_time_name)
    axis.units = freq + " since " + str(year1) + "-" + str(month1) + "-" + str(day1)
    axis.axis = freq
    tab_out.setAxis(time_num, axis)
    return tab_out
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of often used combinations of previous functions
#
def ComputePDF(tab, nbr_bins=10, interval=None, axis_name='axis'):
    """
    #################################################################################
    Description:
    Computes the PDF of tab based on numpy.histogram using nbr_bins in interval
    Returns the density (sum of pdf=1.) in each bin

    Uses uvcdat for array and axis
    #################################################################################

    :param tab: masked_array
        array for which you would like to compute the PDF
    :param nbr_bins: int, optional
        defines the number of equal-width bins in the given range
        default value = 10
    :param interval: [float,float], optional
        The lower and upper range of the bins. Values outside the range are ignored. The first element of the range must
        be less than or equal to the second.
        default value = [tab.min(), tab.max()]
    :param axis_name: string, optional
        name given to the axis of the pdf, e.g. 'longitude'
        default value = 'axis'

    :return: pdf: masked_array
        density (sum of pdf=1.) in each bin, with the center of each bin as axis
    """
    tmp = NPhistogram(tab, bins=nbr_bins, range=interval)
    axis = [(tmp[1][ii] + tmp[1][ii + 1]) / 2. for ii in range(len(tmp[1]) - 1)]
    pdf = MV2array(tmp[0]) / float(len(tab))
    axis = CDMS2createAxis(MV2array(axis, dtype='f'), id=axis_name)
    pdf.setAxis(0, axis)
    return pdf


def CustomLinearRegression(y, x, sign_x=0, return_stderr=True, return_intercept=True):
    """
    #################################################################################
    Description:
    Custom version of genutil.linearregression
    This function offers the possibility to compute the linear regression for all points, for values of x>=0, for values
    of x<=0

    Uses uvcdat
    #################################################################################

    :param y: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param x: masked_array
        masked_array (uvcdat cdms2) containing 'var_name', with many attributes attached (short_name, units,...)
    :param sign_x: int, optional
        default value = 0, computes the linear regression of y over x. You can pass -1 or 1 to compute the linear
        regression of y over x for x >=0 or x<=0 respectively
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :param return_intercept: boolean, optional
        default value = True, returns the the interception value of the linear regression
        True if you want the interception value, if you don't want it pass anything but true
    :return slope, stderr: floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    if sign_x != 0:
        try:
            len(y[0])
        except:
            slope, intercept, stderr = CustomLinearRegression1d(y, x, sign_x=sign_x)
        else:
            if x.shape != y.shape:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": array shape",
                    str().ljust(5) + "different array shape for x " + str(x.shape) + " and y " + str(y.shape),
                ]
                EnsoErrorsWarnings.my_error(list_strings)
            slope, intercept, stderr = MV2zeros(y[0].shape), MV2zeros(y[0].shape), MV2zeros(y[0].shape)
            for ii in range(len(y[0])):
                try:
                    len(y[0, ii])
                except:
                    slope[ii], intercept[ii], stderr[ii] = CustomLinearRegression1d(y[:, ii], x[:, ii], sign_x=sign_x)
                else:
                    for jj in range(len(y[0, ii])):
                        try:
                            len(y[0, ii, jj])
                        except:
                            slope[ii, jj], intercept[ii, jj], stderr[ii, jj] = \
                                CustomLinearRegression1d(y[:, ii, jj], x[:, ii, jj], sign_x=sign_x)
                        else:
                            for kk in range(len(y[0, ii, jj])):
                                try:
                                    len(y[0, ii, jj, kk])
                                except:
                                    slope[ii, jj, kk], intercept[ii, jj, kk], stderr[ii, jj, kk] = \
                                        CustomLinearRegression1d(y[:, ii, jj, kk], x[:, ii, jj, kk], sign_x=sign_x)
                                else:
                                    list_strings = [
                                        "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                                        ": array shape",
                                        str().ljust(5) + str(x.shape) + " too many dimensions (not programmed)",
                                        str().ljust(5) + "Please ckeck and modify the program if needed"
                                    ]
                                    EnsoErrorsWarnings.my_error(list_strings)

    else:
        results = GENUTILlinearregression(y, x=x, error=1, nointercept=None)
        slope, intercept, stderr = results[0][0], results[0][1], results[1][0]
        try:
            slope[0]
        except:
            slope, intercept, stderr = float(slope), float(intercept), float(stderr)
    try:
        len(slope)
    except:
        pass
    else:
        axes = y[0].getAxisList()
        grid = y[0].getGrid()
        mask = y[0].mask
        slope = CDMS2createVariable(MV2array(slope), mask=mask, grid=grid, axes=axes, id='slope')
        stderr = CDMS2createVariable(MV2array(stderr), mask=mask, grid=grid, axes=axes, id='standart_error')
        intercept = CDMS2createVariable(MV2array(intercept), mask=mask, grid=grid, axes=axes, id='intercept')
    if return_stderr is False and return_intercept is False:
        tab = deepcopy(slope)
    else:
        tab = [slope]
        if return_stderr is True:
            tab.append(stderr)
        if return_intercept is True:
            tab.append(intercept)
    return tab


def CustomLinearRegression1d(y, x, sign_x=1):
    x = NParray(x)
    y = NParray(y)
    if sign_x == 1:
        idx = NPnonzero(x > 0.)
    elif sign_x == -1:
        idx = NPnonzero(x < 0.)
    if len(idx[0]) == 0:
        slope, intercept, stderr = 0, 0, 0
    else:
        results = GENUTILlinearregression(y[idx], x=x[idx], error=1, nointercept=None)
        slope, intercept, stderr = float(results[0][0]), float(results[0][1]), float(results[1][0])
    return slope, intercept, stderr


def fill_dict_teleconnection(tab1, tab2, dataset1, dataset2, timebounds1, timebounds2, nyear1, nyear2, nbr, var_name,
                             add_name, units, centered_rmse=0, biased_rmse=1, dict_metric={}, dict_nc={}, ev_name=None,
                             events1=None, events2=None):
    # Metric 1
    rmse_dive = float(RmsAxis(tab1, tab2, axis="xy", centered=centered_rmse, biased=biased_rmse))
    rmse_error_dive = None
    # Metric 2
    corr_dive = float(Correlation(tab1, tab2, axis="xy", centered=1, biased=1))
    corr_error_dive = None
    # Metric 3
    std_mod_dive = Std(tab1, weights=None, axis="xy", centered=1, biased=1)
    std_obs_dive = Std(tab2, weights=None, axis="xy", centered=1, biased=1)
    std_dive = float(std_mod_dive) / float(std_obs_dive)
    std_error_dive = None
    list_met_name = ["RMSE_" + dataset2, "RMSE_error_" + dataset2, "CORR_" + dataset2, "CORR_error_" + dataset2,
                     "STD_" + dataset2, "STD_error_" + dataset2]
    list_metric_value = [rmse_dive, rmse_error_dive, corr_dive, corr_error_dive, std_dive, std_error_dive]
    for tmp1, tmp2 in zip(list_met_name, list_metric_value):
        dict_metric[tmp1 + "_" + add_name] = tmp2
    dict_nc["var" + str(nbr)] = tab1
    dict_dive = {"units": units, "number_of_years_used": nyear1, "time_period": str(timebounds1),
                 "spatialSTD_" + dataset1: std_mod_dive}
    if isinstance(events1, list) is True:
        dict_dive[ev_name + "_years"] = str(events1)
    dict_nc["var" + str(nbr) + "_attributes"] = dict_dive
    dict_nc["var" + str(nbr) + "_name"] = var_name + dataset1
    dict_dive = {"units": units, "number_of_years_used": nyear2, "time_period": str(timebounds2),
                 "spatialSTD_" + dataset2: std_obs_dive}
    if isinstance(events2, list) is True:
        dict_dive[ev_name + "_years"] = str(events2)
    dict_nc["var" + str(nbr + 1)] = tab2
    dict_nc["var" + str(nbr + 1) + "_attributes"] = dict_dive
    dict_nc["var" + str(nbr + 1) + "_name"] = var_name + dataset2
    return dict_metric, dict_nc


def FindXYMinMaxInTs(tab, return_val='both', smooth=False, axis=0, window=5, method='triangle'):
    """
    #################################################################################
    Description:
    Finds in tab in each time step the position (t,x,y,z) of the minimum (return_val='mini') or the maximum
    (return_val='maxi') or both values (if return_val is neither 'mini' nor 'maxi')
    Returned position(s) are not the position in tab but in the (t,x,y,z) space defined by tab axes

    Uses uvcdat for smoothing
    #################################################################################

    :param tab: masked_array
        array for which you would like to know the position (t,x,y,z) of the minimum and/or the maximum values
    :param return_val: string, optional
        'mini' to return the position of the minimum value
        'maxi' to return the position of the maximum value
        to return both minimum and maximum values, pass anything else
        default value = 'both', returns both minimum and maximum values
    :param smooth: boolean, optional
        True if you want to smooth tab, if you do not, pass anything but true
        default value = False, tab is not smoothed

    See function EnsoUvcdatToolsLib.Smoothing
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the moving window average
        default value is 5
    :param method: string, optional
        smoothing method:
            'gaussian': gaussian shaped window
            'square':   square shaped window
            'triangle': triangle shaped window

    :return: minimum/maximum position or both minimum and maximum positions, int, float or list
        position(s) in the (t,x,y,z) space defined by tab axes of the minimum and/or maximum values of tab
    """
    tab_ts = list()
    for tt in range(len(tab)):
        if smooth is True:
            tmp, unneeded = Smoothing(tab[tt], '', axis=axis, window=window, method=method)
        else:
            tmp = deepcopy(tab[tt])
        tab_ts.append(find_xy_min_max(tmp, return_val=return_val))
    tab_ts = MV2array(tab_ts)
    tab_ts.setAxis(0, tab.getAxis(0))
    return tab_ts


def MyDerive(project, internal_variable_name, dict_var):
    # test input parameters
    if not isinstance(project, basestring):
        EnsoErrorsWarnings.object_type_error('project', 'string', type(project), INSPECTstack())
    if not isinstance(internal_variable_name, basestring):
        EnsoErrorsWarnings.object_type_error('internal_variable_name', 'string', type(internal_variable_name),
                                             INSPECTstack())
    if not isinstance(dict_var, dict):
        EnsoErrorsWarnings.object_type_error('project', 'dictionary', type(dict_var), INSPECTstack())

    # get dictionary of observations
    dict_obs = ReferenceObservations()

    # wrong project?
    if 'CMIP' not in project and project not in dict_obs.keys():
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
            ": project",
            str().ljust(5) + "unknown 'project' (or observations dataset): " + str(project),
            str().ljust(10) + "it must be either a 'CMIP' project or an observations dataset defined in\
                    EnsoCollectionsLib.ReferenceObservations",
            str().ljust(10) + "known observations dataset: " + str(sorted(dict_obs.keys(), key=lambda v: v.upper()))
        ]
        EnsoErrorsWarnings.my_error(list_strings)

    # compute 'internal_variable_name' in 'CMIP' case
    if 'CMIP' in project:
        # get dictionary of CMIP
        dict_CMIP = CmipVariables()['variable_name_in_file']
        # test if 'internal_variable_name' is defined in EnsoCollectionsLib.CmipVariables
        if internal_variable_name in dict_CMIP.keys():
            list_var = dict_CMIP[internal_variable_name]['var_name']
            outvar = MyDeriveCompute(list_var, dict_var, dict_att=dict_CMIP, variable=internal_variable_name,
                                     project=project)
    # compute 'internal_variable_name' in 'obs' case
    else:
        # 'project' is defined in EnsoCollectionsLib.ReferenceObservations
        dict_obs_var = dict_obs[project]['variable_name_in_file']
        # test if 'internal_variable_name' is defined for this observations dataset
        if internal_variable_name in dict_obs_var.keys():
            list_var = dict_obs_var[internal_variable_name]['var_name']
            outvar = MyDeriveCompute(list_var, dict_var, dict_att=dict_obs_var, variable=internal_variable_name,
                                     isObs=True)
    return outvar


def MyDeriveCompute(list_var, dict_var, dict_att={}, variable='', isObs=False, project=''):
    # test if keys in list_var are in 'dict_var'
    string_in_dict(list_var, dict_var, INSPECTstack())
    if isinstance(list_var, basestring):
        # this 'internal_variable_name' is based on one variable
        outvar = dict_var[list_var]
    else:
        # this 'internal_variable_name' is based on several variables
        list_operator = dict_att[variable]['algebric_calculation']
        if len(list_operator) != len(list_var):
            if isObs is True:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    ": variable definition in EnsoCollectionsLib.ReferenceObservations(" + str(project) + ")",
                    str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                    str(variable) + " but " + str(len(list_operator)) + " operator(s) are given"
                ]
            else:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    ": variable definition in EnsoCollectionsLib.CmipVariables",
                    str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                    str(variable) + " but " + str(len(list_operator)) + " operator(s) are given"
                ]
            EnsoErrorsWarnings.my_error(list_strings)
        # compute the output variable
        if list_operator[0] == 'minus':
            outvar = -1 * dict_var[list_var[0]]
        else:
            outvar = dict_var[list_var[0]]
        for ii in range(1, len(list_var)):
            outvar = dict_operations[list_operator[ii]](outvar, dict_var[list_var[ii]])
        outvar.setAxisList(dict_var[list_var[0]].getAxisList())
        outvar = MV2masked_where(dict_var[list_var[0]].mask, outvar)
        outvar.setGrid(dict_var[list_var[0]].getGrid())
    return outvar


def LinearRegressionAndNonlinearity(y, x, return_stderr=True, return_intercept=True):
    """
    #################################################################################
    Description:
    CustomLinearRegression applied for all values of x, for values of x>=0, for values of x<=0

    Uses uvcdat
    #################################################################################

    :param y: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param x: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :param return_intercept: boolean, optional
        default value = True, returns the the interception value of the linear regression
        True if you want the interception value, if you don't want it pass anything but true
    :return [slope_all_values, stderr_all_values], [slope_positive_values, stderr_positive_values],
            [slope_negative_values, stderr_negative_values]: lists of floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    # all points
    all_values = CustomLinearRegression(y, x, 0, return_stderr=return_stderr, return_intercept=return_intercept)
    # positive SSTA = El Nino
    positive_values = CustomLinearRegression(y, x, 1, return_stderr=return_stderr, return_intercept=return_intercept)
    # negative SSTA = La Nina
    negative_values = CustomLinearRegression(y, x, -1, return_stderr=return_stderr, return_intercept=return_intercept)
    return all_values, positive_values, negative_values


def LinearRegressionTsAgainstMap(y, x, return_stderr=True):
    """
    #################################################################################
    Description:
    Custom version of genutil.linearregression
    This function offers the possibility to compute the linear regression of a time series against a map

    Uses uvcdat
    #################################################################################

    :param y: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param x: masked_array
        masked_array (uvcdat cdms2) containing 'var_name', with many attributes attached (short_name, units,...)
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :return slope, stderr: floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    tmp = MV2zeros(y.shape)
    for ii in range(len(tmp)):
        tmp[ii].fill(x[ii])
    tmp = CDMS2createVariable(tmp, mask=y.mask, grid=y.getGrid(), axes=y.getAxisList(), id=x.id)
    slope, stderr = GENUTILlinearregression(y, x=tmp, error=1, nointercept=1)
    if return_stderr:
        return slope, stderr
    else:
        return slope


def LinearRegressionTsAgainstTs(y, x, nbr_years_window, return_stderr=True, frequency=None, debug=False):
    """
    #################################################################################
    Description:
    Custom version of genutil.linearregression
    This function offers the possibility to compute the linear regression of a time series against a lead-lag time
    series

    Uses uvcdat
    #################################################################################

    :param y: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param x: masked_array
        masked_array (uvcdat cdms2) containing 'var_name', with many attributes attached (short_name, units,...)
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param debug: boolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return slope, stderr: floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    if frequency == 'daily':
        nbr_timestep = nbr_years_window * 365
    elif frequency == 'monthly':
        nbr_timestep = nbr_years_window * 12
    elif frequency == 'yearly':
        nbr_timestep = nbr_years_window
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    tab_yy_mm = Event_selection(y, frequency, nbr_years_window=nbr_years_window)
    myshape = [nbr_timestep] + [ss for ss in y.shape[1:]]
    tmp_ax = CDMS2createAxis(range(nbr_timestep), id='months')
    slope_out = MV2zeros(myshape)
    slope_out.setAxisList([tmp_ax] + y.getAxisList()[1:])
    stderr_out = MV2zeros(myshape)
    stderr_out.setAxisList([tmp_ax] + y.getAxisList()[1:])
    for ii in range(nbr_timestep):
        tmp1 = tab_yy_mm[:, ii]
        tmp2 = deepcopy(x)
        yy1 = tab_yy_mm.getAxis(0)[0]
        yy2 = tmp2.getTime().asComponentTime()[0].year
        if yy1 == yy2:
            tmp1 = tmp1[:len(tmp2)]
        elif yy1 < yy2:
            tmp1 = tmp1[yy2 - yy1:len(tmp2)]
        else:
            tmp2 = tmp2[yy2 - yy1:]
            tmp1 = tmp1[:len(x)]
        if len(tmp2) > len(tmp1):
            tmp2 = tmp2[:len(tmp1)]
        # if debug is True:
        #     yy1 = tmp1.getAxis(0)[0]
        #     yy2 = tmp2.getTime().asComponentTime()[0].year
        #     EnsoErrorsWarnings.DebugMode('\033[93m', "EnsoUvcdatToolsLib LinearRegressionTsAgainstTs", 20)
        #     dict_debug = {'axes1': str([ax.id for ax in tmp1.getAxisList()]), 'shape1': str(tmp1.shape),
        #                   'line1': "first year is " + str(yy1),
        #                   'axes2': str([ax.id for ax in tmp2.getAxisList()]), 'shape2': str(tmp2.shape),
        #                   'line2': "first year is " + str(yy2)}
        #     EnsoErrorsWarnings.DebugMode('\033[93m', str(x.id) + " regressed against " + str(y.id), 25, **dict_debug)
        if tmp2.shape == tmp1.shape:
            tmp3 = deepcopy(tmp2)
        else:
            tmp3 = MV2zeros(tmp1.shape)
            for jj in range(len(tmp3)):
                tmp3[jj].fill(tmp2[jj])
        tmp3 = CDMS2createVariable(tmp3, mask=tmp1.mask, grid=tmp1.getGrid(), axes=tmp1.getAxisList(), id=x.id)
        slope, stderr = GENUTILlinearregression(tmp1, x=tmp3, error=1, nointercept=1)
        slope_out[ii] = slope
        stderr_out[ii] = stderr
        del slope, stderr, tmp1, tmp2, tmp3, yy1, yy2
    if return_stderr:
        return slope_out, stderr_out
    else:
        return slope_out


def PreProcessTS(tab, info, areacell=None, average=False, compute_anom=False, compute_sea_cycle=False, debug=False,
                 region=None, **kwargs):
    # removes annual cycle (anomalies with respect to the annual cycle)
    if compute_anom is True:
        tab = ComputeInterannualAnomalies(tab)
    # Normalization of the anomalies
    if kwargs['normalization']:
        if kwargs['frequency'] is not None:
            tab = Normalize(tab, kwargs['frequency'])
            info = info + ', normalized'
    # Removing linear trend
    if isinstance(kwargs['detrending'], dict):
        known_args = {'axis', 'method', 'bp'}
        extra_args = set(kwargs['detrending']) - known_args
        if extra_args:
            EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
        tab, info = Detrend(tab, info, **kwargs['detrending'])
    # Smoothing time series
    if isinstance(kwargs['smoothing'], dict):
        known_args = {'axis', 'method', 'window'}
        extra_args = set(kwargs['smoothing']) - known_args
        if extra_args:
            EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
        tab, info = Smoothing(tab, info, **kwargs['smoothing'])
    # computes mean annual cycle
    if compute_sea_cycle is True:
        tab = annualcycle(tab)
    # average
    if average is not False:
        if debug is True:
            EnsoErrorsWarnings.debug_mode('\033[93m', "EnsoUvcdatToolsLib PreProcessTS", 20)
            dict_debug = {'axes1':  str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', "averaging to perform: " + str(average), 25, **dict_debug)
        if isinstance(average, basestring):
            try: dict_average[average]
            except:
                EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
            else:
                tab = dict_average[average](tab, areacell, region=region, **kwargs)
                if debug is True:
                    dict_debug = {'axes1': str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                    EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(average), 25, **dict_debug)
        elif isinstance(average, list):
            for av in average:
                try: dict_average[av]
                except:
                    EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
                else:
                    tab = dict_average[av](tab, areacell, region=region, **kwargs)
                    if debug is True:
                        dict_debug = {'axes1': str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                        EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(av), 25, **dict_debug)
        else:
            EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
    return tab, info


def ReadSelectRegionCheckUnits(filename, varname, varfamily, box=None, time_bounds=None, frequency=None, **keyarg):
    """
    #################################################################################
    Description:
    Combines ReadAndSelectRegion and CheckUnits
    Reads the given 'varname' from the given 'filename', selects the given 'box' and checks the 'varname''s units
    depending on 'vartype'

    Uses uvcdat
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param varname: string
        name of the variable to read from 'filename'
    :param varfamily: string
        family of variable encompassing 'varname' (temperature, velocity,...)
    :param box: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None

    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    tab = ReadAndSelectRegion(filename, varname, box=box, time_bounds=time_bounds, frequency=frequency)
    tab, units, keyerror = CheckUnits(tab, varfamily, varname, tab.units, return_tab_only=False)
    tab.name = varname
    tab.units = units
    return tab, keyerror


def Read_data_mask_area(file_data, name_data, type_data, metric, region, file_area='', name_area='', file_mask='',
                        name_mask='', maskland=False, maskocean=False, time_bounds=None, debug=False, **kwargs):
    keyerror1, keyerror2, keyerror3 = None, None, None
    # Read variable
    if debug is True:
        dict_debug = {'file1': '(' + type_data + ') ' + str(file_data), 'var1': '(' + type_data + ') ' + str(name_data)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'Files', 20, **dict_debug)
    variable, keyerror1 = ReadSelectRegionCheckUnits(file_data, name_data, type_data, box=region,
                                                     time_bounds=time_bounds, **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(' + type_data + ') ' + str([ax.id for ax in variable.getAxisList()]),
                      'shape1': '(' + type_data + ') ' + str(variable.shape),
                      'time1': '(' + type_data + ') ' + str(TimeBounds(variable))}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadSelectRegionCheckUnits', 20, **dict_debug)
    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        if len(variable) < kwargs['min_time_steps']:
            EnsoErrorsWarnings.too_short_time_period(metric, len(variable), kwargs['min_time_steps'], INSPECTstack())
            keyerror2 = "too short time period (" + str(len(variable)) + ")"
    # Read areacell & mask
    variable, areacell, keyerror3 = Read_mask_area(variable, file_data, type_data, region, file_area=file_area,
                                                   name_area=name_area, file_mask=file_mask, name_mask=name_mask,
                                                   maskland=maskland, maskocean=maskocean, debug=debug, **kwargs)
    if keyerror1 is not None or keyerror2 is not None or keyerror3 is not None:
        keyerror = ''
        if keyerror1 is not None:
            keyerror = keyerror1
        if len(keyerror) > 0 and keyerror2 is not None:
            keyerror += " ; "
        if keyerror2 is not None:
            keyerror += keyerror2
        if len(keyerror) > 0 and keyerror3 is not None:
            keyerror += " ; "
        if keyerror3 is not None:
            keyerror += keyerror3
    else:
        keyerror = None
    return variable, areacell, keyerror


def Read_data_mask_area_multifile(file_data, name_data, type_data, variable, metric, region, file_area='', name_area='',
                                  file_mask='', name_mask='', maskland=False, maskocean=False, debug=False,
                                  interpreter='', **kwargs):
    dict_area, dict_keye, dict_var = dict(), dict(), dict()
    if isinstance(file_data, basestring):
        tab, areacell, keyerror = \
            Read_data_mask_area(file_data, name_data, type_data, metric, region, file_area=file_area,
                                name_area=name_area, file_mask=file_mask, name_mask=name_mask, maskland=maskland,
                                maskocean=maskocean, debug=debug, **kwargs)
        dict_area[name_data], dict_keye[name_data], dict_var[name_data] = areacell, keyerror, tab
    else:
        for ii in range(len(file_data)):
            try:
                file_data[ii]
            except:
                ff1 = ''
            else:
                ff1 = file_data[ii]
            try:
                name_data[ii]
            except:
                nn1 = ''
            else:
                nn1 = name_data[ii]
            try:
                file_area[ii]
            except:
                fa1 = ''
            else:
                fa1 = file_area[ii]
            try:
                name_area[ii]
            except:
                an1 = ''
            else:
                an1 = name_area[ii]
            try:
                file_mask[ii]
            except:
                fl1 = ''
            else:
                fl1 = file_mask[ii]
            try:
                name_mask[ii]
            except:
                ln1 = ''
            else:
                ln1 = name_mask[ii]
            tab, areacell, keyerror = \
                Read_data_mask_area(ff1, nn1, type_data, metric, region, file_area=fa1, name_area=an1, file_mask=fl1,
                                    name_mask=ln1, maskland=maskland, maskocean=maskocean, debug=debug, **kwargs)
            dict_area[nn1], dict_keye[nn1], dict_var[nn1] = areacell, keyerror, tab
    list_var = sorted(dict_var.keys())
    if len(list_var) > 1:
        for ii in range(2):
            for var in list_var[1:]:
                dict_var[list_var[0]], dict_var[var], unneeded =\
                    CheckTime(dict_var[list_var[0]], dict_var[var], metric_name=metric, **kwargs)
    tab = MyDerive(kwargs[interpreter], variable, dict_var)
    areacell = dict_area[dict_area.keys()[0]]
    keyerror = ''
    for ii in dict_keye.keys():
        if len(keyerror) > 0 and dict_keye[ii] is not None:
            keyerror += " ; "
        if dict_keye[ii] is not None:
            keyerror += dict_keye[ii]
    if len(keyerror) == 0:
        keyerror = None
    return tab, areacell, keyerror


def Read_mask_area(tab, file_data, type_data, region, file_area='', name_area='', file_mask='', name_mask='',
                   maskland=False, maskocean=False, debug=False, **kwargs):
    tab_out = deepcopy(tab)
    keyerror1, keyerror2 = None, None
    # Read areacell
    if file_area:
        areacell = ReadAreaSelectRegion(file_area, areaname=name_area, box=region, **kwargs)
    else:
        areacell = ReadAreaSelectRegion(file_data, areaname=name_area, box=region, **kwargs)
    if debug is True:
        if areacell is not None:
            dict_debug = {'axes1': '(' + type_data + ') ' + str([ax.id for ax in areacell.getAxisList()]),
                          'shape1': '(' + type_data + ') ' + str(areacell.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'areacell is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
    # Read landmask
    if file_mask:
        landmask = ReadLandmaskSelectRegion(tab, file_mask, landmaskname=name_mask, box=region, **kwargs)
    else:
        landmask = ReadLandmaskSelectRegion(tab, file_data, landmaskname=name_mask, box=region, **kwargs)
    if debug is True:
        if landmask is not None:
            dict_debug = {'axes1': '(' + type_data + ') ' + str([ax.id for ax in landmask.getAxisList()]),
                          'shape1': '(' + type_data + ') ' + str(landmask.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'landmask is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
    # Apply landmask
    if landmask is not None:
        tab_out, keyerror1 = ApplyLandmask(tab_out, landmask, maskland=maskland, maskocean=maskocean)
        if keyerror1 is None:
            if areacell is None:
                areacell = ArrayOnes(landmask, id='areacell')
            areacell, keyerror2 = ApplyLandmaskToArea(areacell, landmask, maskland=maskland, maskocean=maskocean)
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = ''
        if keyerror1 is not None:
            keyerror = keyerror1
        if len(keyerror) > 0 and keyerror2 is not None:
            keyerror += " ; "
        if keyerror2 is not None:
            keyerror += keyerror2
    else:
        keyerror = None
    return tab_out, areacell, keyerror


def SlabOcean(tab1, tab2, month1, month2, events, frequency=None, debug=False):
    """
    #################################################################################
    Description:
    Compute a simple slab ocean by integrating the total heat fluxes over time

    Based on:
    Bayr, T., C. Wengel, M. Latif, D. Dommenget, J. Lübbecke, W. Park (2018) Error compensation of ENSO atmospheric
    feedbacks in climate models and its influence on simulated ENSO dynamics. Clim. Dyn., doi:10.1007/s00382-018-4575-7

    Uses CDAT
    #################################################################################

    :param tab1: masked_array
        masked_array (uvcdat cdms2) containing SSTA, with many attributes attached (short_name, units,...)
    :param tab2: masked_array
        masked_array (uvcdat cdms2) containing THFA, with many attributes attached (short_name, units,...)
    :param month1: string
        first month of integration (e.g., 'JUN')
    :param month2: string
        last month of integration (e.g., 'DEC')
    :param events: list of integer
        list of the years considered as ENSO events to be selected
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return dSST, dSSTthf, dSSToce: masked_array
        normalized cumulative SST change (from 0 to 1; in C/C)
        normalized cumulative heat flux-driven SST change (in C/C)
        normalized cumulative SST change by an anomalous ocean circulation (in C/C)
    """
    if debug is True:
        EnsoErrorsWarnings.debug_mode('\033[93m', "EnsoUvcdatToolsLib SlabOcean", 20)
    # months and associated position
    list_months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    if month1 in list_months and month2 in list_months:
        mm1 = list_months.index(month1)
        mm2 = list_months.index(month2)
    else:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": month"]
        if month1 not in list_months:
            list_strings.append(str().ljust(5) + "unknown month1 : " + str(month1))
        if month2 not in list_months:
            list_strings.append(str().ljust(5) + "unknown month2 : " + str(month2))
        EnsoErrorsWarnings.my_error(list_strings)
    # sea water constants
    cp = 4000   # J/(kg * K) (specific heat capacity at constant pressure of sea water)
    rho = 1024  # kg/m3      (average density of sea water)
    H = 50      # m          (depth of the slab ocean)
    fraction = 60 * 60 * 24 * 30.42 / (cp * rho * H)  # W/m2 to C
    # selecting events
    sstA = Event_selection(tab1, frequency, nbr_years_window=2, list_event_years=events)
    thfA = Event_selection(tab2, frequency, nbr_years_window=2, list_event_years=events)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sstA.getAxisList()]),
                      'axes2': '(thf) ' + str([ax.id for ax in thfA.getAxisList()]),
                      'shape1': '(sst) ' + str(sstA.shape), 'shape2': '(thf) ' + str(thfA.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after Event_selection', 25, **dict_debug)
    # cumulative anomalies
    myshape = [len(events), mm2-mm1+1] + [ss for ss in tab1.shape[1:]]
    dSST = MV2zeros(myshape)
    dSSTthf = MV2zeros(myshape)
    for ii in range(mm1, mm2):
        dSST[:, ii - mm1 + 1] = dSST[:, ii - mm1] + sstA[:, ii + 1] - sstA[:, ii]
        dSSTthf[:, ii - mm1 + 1] = dSSTthf[:, ii - mm1] + thfA[:, ii + 1]
    if debug is True:
        dict_debug = {'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after cumulative_anomalies', 25, **dict_debug)
    # normalized heat flux-driven SST change
    dt = MV2zeros(dSSTthf.shape)
    dt = dt.reorder('10')
    dt[:] = dSST[:, -1]
    dt = dt.reorder('10')
    dt = MV2masked_where(abs(dt) < 0.1, dt)
    dSSTthf[:] = fraction * dSSTthf[:] / dt
    # normalized SST change
    dSST[:] = dSST[:] / dt
    # normalized SST change by an anomalous ocean circulation
    dSSToce = dSST - dSSTthf
    # averaging across events
    dSST = MV2average(dSST, axis=0)
    dSSTthf = MV2average(dSSTthf, axis=0)
    dSSToce = MV2average(dSSToce, axis=0)
    # axes
    axes = [CDMS2createAxis(MV2array(range(12-len(dSST), 12)), id='months')]
    if debug is True:
        dict_debug = {'axes1': 'axes ' + str(axes[0]), 'axes2': 'axes[:] ' + str(axes[0][:]),
                      'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape),
                      'shape3': '(dSSToce) ' + str(dSSToce.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after mean dSST', 25, **dict_debug)
    if len(tab1.shape) > 1:
        axes = axes + tab1.getAxisList()[1:]
    dSST.setAxisList(axes)
    dSSTthf.setAxisList(axes)
    dSSToce.setAxisList(axes)
    if debug is True:
        dict_debug = {'axes1': '(dSST) ' + str([ax.id for ax in dSST.getAxisList()]),
                      'axes2': '(dSSTthf) ' + str([ax.id for ax in dSSTthf.getAxisList()]),
                      'axes3': '(dSSToce) ' + str([ax.id for ax in dSSToce.getAxisList()]),
                      'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape),
                      'shape3': '(dSSToce) ' + str(dSSToce.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'output', 25, **dict_debug)
    return dSST, dSSTthf, dSSToce


def TimeAnomaliesLinearRegressionAndNonlinearity(tab2, tab1, return_stderr=True):
    """
    #################################################################################
    Description:
    LinearRegressionAndNonlinearity applied on two 'raw' masked_arrays (i.e., the annual cycle is not removed and the
    spatial average is not computed)
    The linear regression of tab2 on tab1 is computed for all values of tab1, for values of tab1>=0, for values of
    tab1<=0

    Uses uvcdat
    #################################################################################

    :param tab2: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param tab1: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :return: [slope_all_values, stderr_all_values], [slope_positive_values, stderr_positive_values],
            [slope_negative_values, stderr_negative_values]: lists of floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    # horizontal average
    tab1 = dict_average['horizontal'](tab1)
    tab2 = dict_average['horizontal'](tab2)
    # removes annual cycle (anomalies with respect to the annual cycle)
    tab1 = cdutil.ANNUALCYCLE.departures(tab1)
    tab2 = cdutil.ANNUALCYCLE.departures(tab2)
    # computes linear regression of tab2 on tab1 for all values of tab1, for values of tab1>=0, for values of tab1<=0
    lr, lrpos, lrneg = LinearRegressionAndNonlinearity(tab2, tab1, return_stderr=return_stderr)
    return lr, lrpos, lrneg


def TimeAnomaliesStd(tab):
    """
    #################################################################################
    Description:
    Combines cdutil.averager and genutil.std
    Averages spatially and computes the standard deviation

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :return std: float
        standard deviation (one value) of the masked_array averaged spatially and with the annual cycle removed
    """
    # horizontal average
    tab = dict_average['horizontal'](tab)
    # computes standard deviation
    std = float(GENUTILstd(tab, weights=None, axis=0, centered=1, biased=1))
    return std


def TsToMap(tab, map_ref):
    """
    #################################################################################
    Description:
    Put a time series into a nD array according to the reference map

    Uses uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param map_ref: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :return map_out: masked_array
        tab (1D) values set into map_ref shape (nD)
    """
    if len(map_ref.shape) > 6:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many dimensions",
            str().ljust(5) + "map_ref.shape = " + str(map_ref.shape)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    map_out = MV2zeros(map_ref.shape)
    map_out = CDMS2createVariable(map_out, axes=map_ref.getAxisList(), grid=map_ref.getGrid(), mask=map_ref.mask,
                                  attributes=map_ref.attributes, id=tab.id)
    initorder = map_out.getOrder()
    map_out = map_out.reorder('...t')
    if len(map_ref.shape) == 2:
        map_out[:] = tab
    elif len(map_ref.shape) == 3:
        map_out[:, :] = tab
    elif len(map_ref.shape) == 4:
        map_out[:, :, :] = tab
    elif len(map_ref.shape) == 5:
        map_out[:, :, :, :] = tab
    else:
        map_out[:, :, :, :, :] = tab
    map_out = map_out.reorder(initorder)
    return map_out


def TwoVarRegrid(model, obs, info, region=None, model_orand_obs=0, newgrid=None, **keyarg):
    """
    #################################################################################
    Description:
    Regrids 'model', 'obs' or both

    Uses uvcdat
    #################################################################################

    :param model: masked_array
        model data
    :param obs: masked_array
        observations data
    :param info: string
        information about what was done to 'model' and 'obs'
    :param region: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param model_orand_obs: integer, optional
        0 if you want to regrid model data toward observations data
        1 if you want to regrid observations data toward model data
        2 if you want to regrid model AND observations data toward 'newgrid'
        default value = 0
    :param newgrid: CDMS grid
        grid toward which model data and observations data are regridded if model_orand_obs=2

    usual kwargs:
    :param newgrid_name: string, optional
        generates 'newgrid' depending on the given string using cdms2
        'newgrid_name' must contain a name of a grid type:
            'equalarea', 'gaussian', 'generic', 'uniform'
        and a name of a grid resolution:
            '0.25x0.25deg', '0.5x0.5deg', '1x1deg', '2x2deg'
        default value = 'generic 1x1deg'
        for more information:
        import cdms2
        help(cdms2.createUniformLatitudeAxis)
        help(cdms2.createUniformLongitudeAxis)
        help(cdms2.createRectGrid)
    see EnsoUvcdatToolsLib.Regrid for regridding options

    :return: model, obs, info
        model and obs on the same grid, and information about what has been done to 'model' and 'obs'
    """
    known_args = {'missing', 'order', 'mask', 'newgrid_name', 'regridder', 'regridTool', 'regridMethod'}
    extra_args = set(keyarg) - known_args
    if extra_args:
        EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
    grid_obs = obs.getGrid()
    grid_model = model.getGrid()
    # select case:
    if model_orand_obs == 0:
        model = Regrid(model, grid_obs, **keyarg)
        info = info + ', model regridded to observations'
    elif model_orand_obs == 1:
        obs = Regrid(obs, grid_model, **keyarg)
        info = info + ', observations regridded to model'
    elif model_orand_obs == 2:
        model = Regrid(model, newgrid, region=region, **keyarg)
        obs = Regrid(obs, newgrid, region=region, **keyarg)
        try: grid_name = newgrid.id
        except:
            try: grid_name = newgrid.name
            except:
                try: grid_name = keyarg['newgrid_name']
                except: grid_name = 'newgrid'
        info = info + ', observations and model regridded to ' + str(grid_name)
    else:
        info = info + ', observations and model NOT regridded'
    if model.shape == obs.shape:
        if model.mask.shape != ():
            mask = model.mask
            if obs.mask.shape != ():
                mask = MV2where(obs.mask, obs.mask, mask)
        else:
            if obs.mask.shape != ():
                mask = obs.mask
            else:
                mask = MV2where(MV2zeros(model.shape)==0, False, True)
        model = MV2masked_where(mask, model)
        obs = MV2masked_where(mask, obs)
    else:
        if obs[0].mask.shape != ():
            tab = MV2zeros(model.shape)
            for tt in range(len(tab)):
                tab[tt] = MV2masked_where(obs[0].mask, tab[tt])
            model = MV2masked_where(tab.mask, model)
        if model[0].mask != ():
            tab = MV2zeros(obs.shape)
            for tt in range(len(tab)):
                tab[tt] = MV2masked_where(model[0].mask, tab[tt])
            obs = MV2masked_where(tab.mask, obs)
    return model, obs, info
