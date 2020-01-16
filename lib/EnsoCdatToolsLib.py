# -*- coding:UTF-8 -*-
from __future__ import print_function
from calendar import monthrange
from copy import deepcopy
from inspect import stack as INSPECTstack
import ntpath
from numpy import exp as NPexp
from os.path import isdir as OSpath_isdir
from os.path import isfile as OSpath__isfile
from scipy.signal import detrend as SCIPYsignal_detrend

# ENSO_metrics package functions:
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings
from EnsoToolsLib import add_up_errors, closest_grid

# CDAT based functions:
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import createRectGrid as CDMS2createRectGrid
from cdms2 import createUniformLatitudeAxis as CDMS2createUniformLatitudeAxis
from cdms2 import createUniformLongitudeAxis as CDMS2createUniformLongitudeAxis
from cdms2 import createVariable as CDMS2createVariable
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
from cdtime import comptime as CDTIMEcomptime
import cdutil
from genutil.statistics import std as GENUTILstd
from MV2 import add as MV2add
from MV2 import array as MV2array
from MV2 import average as MV2average
from MV2 import compress as MV2compress
from MV2 import divide as MV2divide
from MV2 import masked_where as MV2masked_where
from MV2 import maximum as MV2maximum
from MV2 import minimum as MV2minimum
from MV2 import multiply as MV2multiply
from MV2 import ones as MV2ones
from MV2 import subtract as MV2subtract
from MV2 import sum as MV2sum
from MV2 import take as MV2take
from MV2 import zeros as MV2zeros
from regrid2.horizontal import Horizontal as REGRID2horizontal__Horizontal


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of simple CDAT functions used in EnsoMetricsLib.py
#
def array_ones(tab, var_id='new_variable_ones'):
    """
    #################################################################################
    Description:
    Create a masked_array filled with ones with the same properties as tab (shape, axes, grid, mask)
    #################################################################################
    for more information:
    import MV2
    help(MV2.ones)
    """
    return CDMS2createVariable(MV2ones(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask,
                               id=var_id)


def average_horizontal(tab, areacell=None, region=None, **kwargs):
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
            print("\033[93m" + str().ljust(15) + "EnsoCDATToolsLib AverageHorizontal" + "\033[0m")
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
                tmp = compute_regrid(tab, None, region=region, **kwargs2)
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


def average_meridional(tab, areacell=None, region=None, **kwargs):
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
            print("\033[93m" + str().ljust(15) + "EnsoCDATToolsLib AverageMeridional" + "\033[0m")
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
                tmp = compute_regrid(tab, None, region=region, **kwargs2)
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


def average_temporal(tab, areacell=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 't' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    averaged_tab = None
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


def average_zonal(tab, areacell=None, region=None, **kwargs):
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
            print("\033[93m" + str().ljust(15) + "EnsoCDATToolsLib AverageZonal" + "\033[0m")
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
                tmp = compute_regrid(tab, None, region=region, **kwargs2)
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
dict_average = {'horizontal': average_horizontal, 'meridional': average_meridional, 'time': average_temporal,
                'zonal': average_zonal}


def compute_interannual_anomalies(tab):
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


def operation_add(tab, number_or_tab):
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


def operation_divide(tab, number_or_tab):
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


def operation_multiply(tab, number_or_tab):
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


def operation_subtract(tab, number_or_tab):
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
dict_operations = {'divide': operation_divide, 'minus': operation_subtract, 'multiply': operation_multiply,
                   'plus': operation_add}


def compute_std(tab, weights=None, axis=0, centered=1, biased=1):
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


def get_time_bounds(tab):
    """
    #################################################################################
    Description:
    Finds first and last dates of tab's time axis, tab must be a CDAT masked_array
    #################################################################################

    Returns a tuple of strings: e.g., ('1979-1-1 11:59:60.0', '2016-12-31 11:59:60.0')
    """
    time = tab.getTime().asComponentTime()
    return str(time[0]), str(time[-1])
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of more complex functions (based on CDAT) used in EnsoMetricsLib.py
#
def apply_landmask(tab, landmask, mask_land=True, mask_ocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if mask_land is True, mask where landmask==100
        if mask_ocean is True, mask where landmask==0
    #################################################################################

    :param tab: masked_array
    :param landmask: masked_array
    :param mask_land: boolean, optional
        masks land points
        default value is True
    :param mask_ocean: boolean, optional
        masks ocean points
        default value is False

    :return tab: masked_array
        masked_array where land points and/or ocean points are masked
    :return keyerror: string or None
        name of the error encountered if applicable, else None
    """
    keyerror = None
    if mask_land is True or mask_ocean is True:
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
                if mask_land is True:
                    tab = MV2masked_where(landmask_nd == 1, tab)
                if mask_ocean is True:
                    tab = MV2masked_where(landmask_nd == 0, tab)
    return tab, keyerror


def apply_landmask_to_area(area, landmask, mask_land=True, mask_ocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if mask_land is True, mask where landmask==1 and area=area*(1-landmask) (to weight island and coastal points)
        if mask_ocean is True, mask where landmask==0 and area=area*landmask (to weight island and coastal points)
    #################################################################################

    :param area: masked_array
        areacell
    :param landmask: masked_array
    :param mask_land: boolean, optional
        masks land points and weights island and coastal points
        default value is True
    :param mask_ocean: boolean, optional
        masks ocean points and weights island and coastal points
        default value is False

    :return area: masked_array
        masked_array where land points and/or ocean points are masked
    :return keyerror: string or None
        name of the error encountered if applicable, else None
    """
    keyerror = None
    if mask_land is True or mask_ocean is True:
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
            if mask_land:
                area = MV2masked_where(landmask == 1, area)
                area = MV2multiply(area, 1-landmask)
            if mask_ocean:
                area = MV2masked_where(landmask == 0, area)
                area = MV2multiply(area, landmask)
    return area, keyerror


def check_time(tab1, tab2, frequency='monthly', min_time_steps=None, metric_name='', debug=False, **kwargs):
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
    :param debug: boolean, optional
        default value is False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)

    :return tab1_sliced: masked_array
        tab1 but only during the period common to tab1 and tab2
    :return tab2_sliced: masked_array
        tab2 but only during the period common to tab1 and tab2
    :return keyerror: string or None
        name of the error encountered if applicable, else None
    """
    if debug is True:
        dict_debug = {'shape1': 'tab1.shape = ' + str(tab1.shape),'shape2': 'tab2.shape = ' + str(tab2.shape),
                      'time1': 'tab1.time = ' + str(get_time_bounds(tab1)),
                      'time2': 'tab2.time = ' + str(get_time_bounds(tab2))}
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
        stime_adjust = 0
        etime_adjust = 1
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())

    # retains only the time-period common to both tab1 and tab2
    tab1_sliced = tab1(time=(stime_adjust, etime_adjust))
    tab2_sliced = tab2(time=(stime_adjust, etime_adjust))
    if debug is True:
        dict_debug = {'shape1': 'tab1.shape = ' + str(tab1_sliced.shape),
                      'shape2': 'tab2.shape = ' + str(tab2_sliced.shape),
                      'time1': 'tab1.time = ' + str(get_time_bounds(tab1_sliced)),
                      'time2': 'tab1.time = ' + str(tab1_sliced.getTime().asComponentTime()[:]),
                      'time3': 'tab2.time = ' + str(get_time_bounds(tab2_sliced)),
                      'time4': 'tab2.time = ' + str(tab2_sliced.getTime().asComponentTime()[:])}
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


def check_units(tab, var_name, name_in_file, units, return_tab_only=True, **kwargs):
    """
    #################################################################################
    Description:
    Checks the units of the variable and changes it if necessary
    Works for current/wind velocities, heat flux, precipitation, pressure, temperature, wind stress

    Uses MV2 (CDAT) to find the minimum value, to multiply and to subtract
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
        default value is True, only the tab is returned
        True if you want only the tab, if you want the new units also pass anything but true

    :return tab: array
        array with new units (if applicable)
    :return units: string
        name of the new units (if applicable)
    :return keyerror: string or None
        name of the error encountered if applicable, else None
    """
    keyerror = None
    if var_name in ['temperature']:
        if units == 'K':
            # check if the temperature units is really K
            if float(MV2minimum(tab)) > 150:
                # unit change of the temperature: from K to degC
                tab = dict_operations['minus'](tab, 273.15)
                units = "degC"
            else:
                minmax = [MV2minimum(tab),MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, units, minmax, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        elif units in ['C', 'degree_Celsius', 'deg_Celsius', 'deg. C', 'degCelsius', 'degree_C', 'deg_C', 'degC',
                       'degrees C']:
            # check if the temperature units is really degC
            if float(MV2minimum(tab)) < 50:
                units = "degC"
            else:
                minmax = [MV2minimum(tab), MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, minmax, units, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['precipitations']:
        if units == 'kg m-2 s-1':
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = dict_operations['multiply'](tab, 86400)
        elif units == 'mm/day':
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['wind stress']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "N/m2"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['velocity']:
        if units in ['cm s-1', 'cm/s', 'cm s**-1']:
            # unit change of the velocity: from cm/s to m/s
            tab = dict_operations['multiply'](tab, 1e-2)
            units = "m/s"
        elif units in ['m s-1', 'm/s', 'm s**-1', 'm/sec']:
            units = "m/s"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['heat flux']:
        if units in ['W/m2', 'W m-2', 'W/m^2']:
            units = "W/m2"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['pressure']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "Pa"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    elif var_name in ['sea surface height']:
        if units in ['cm', 'centimeter']:
            # unit change of the sea surface height: from cm to m
            tab = dict_operations['multiply'](tab, 1e-2)
            units = "m"
        elif units in ['m', 'meter']:
            units = "m"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.my_warning(list_strings)
    if return_tab_only is True:
       return tab
    else:
        return tab, units, keyerror


def compute_annual_cycle(tab):
    """
    #################################################################################
    Description:
    Computes the annual cycle (climatological value of each calendar month) of tab
    #################################################################################

    :param tab: masked_array
    :return moy: array
        array of the monthly annual cycle
    """
    initial_order = tab.getOrder()
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
    moy = moy.reorder(initial_order)
    time = CDMS2createAxis(range(12), id='months')
    moy.setAxis(get_num_axis(moy, 'time'), time)
    return moy


def compute_estimate_landmask(tab):
    """
    #################################################################################
    Description:
    Estimate landmask (when landmask was not given)
    Uses cdutil (CDAT) to create estimated landmask for model resolution
    #################################################################################

    :param tab: masked_array
        masked_array of data (sst, taux,...), used to estimate the landmask

    :return landmask: masked_array
        masked_array containing landmask
    """
    print('\033[93m' + str().ljust(25) + 'NOTE: Estimated landmask applied' + '\033[0m')
    n = 1
    sft = cdutil.generateLandSeaMask(tab(*(slice(0, 1),) * n)) * 100.0
    sft[:] = sft.filled(100.0)
    lmsk = sft
    lmsk.setAxis(0, tab.getAxis(1))
    lmsk.setAxis(1, tab.getAxis(2))
    lmsk.id = 'sftlf'
    return lmsk


def compute_regrid(tab_to_regrid, new_grid, missing=None, order=None, mask=None, regridder='cdms',
                   regridTool='esmf', regridMethod='linear', **kwargs):
    """
    #################################################################################
    Description:
    Regrids 'tab_to_regrid' to 'new_grid'
    #################################################################################

    for more information:
    import cdms2
    help(cdms2.avariable)

    :param tab_to_regrid: masked_array
        masked_array to regrid (must include a CDMS grid!)
    :param new_grid: CDMS grid
        destination grid
    :param missing: float, optional
        missing values (missing data value, if any)
    :param order: string, optional
        axis order (form "tzyx", "tyx", etc.)
    :param mask: array of booleans, optional
        mask of the new grid (either 2-D or the same shape as new_grid)
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
        if 'new_grid' is not defined (is string) this will be used to create a grid
        name of a region, domain where the grid will be defined, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return new_tab: masked_array
        tab_to_regrid regridded on new_grid
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
            str().ljust(10) + "known regridder: " + str(known_regridder)]
        EnsoErrorsWarnings.my_error(list_strings)
    elif regridder == "cdms":
        if regridTool in ["regrid2", "libcf"]:
            list_method = [None, "linear"]
        elif regridTool == "esmf":
            list_method = [None, "conserve", "linear", "patch"]
        if (regridTool is not None) and (regridTool not in ["regrid2", "esmf", "libcf"]):
            list_method = list()
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridTool",
                str().ljust(5) + "unknown regridTool: " + str(regridTool),
                str().ljust(10) + "known regridTool: " + str(["regrid2", "esmf", "libcf"])]
            EnsoErrorsWarnings.my_error(list_strings)
        elif regridMethod not in list_method:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridMethod",
                str().ljust(5) + "unknown regridMethod (" + str(regridMethod) + ") for this regridTool ("
                + str(regridTool) + ")",
                str().ljust(10) + "known regridMethod: " + str(list_method)]
            EnsoErrorsWarnings.my_error(list_strings)
    # test the given 'new_grid'
    if isinstance(new_grid, str) is True or isinstance(new_grid, basestring) is True or new_grid is None:
        #
        # new_grid is not a grid, so a grid will be created
        # to do this, kwargs['newgrid_name'] and kwargs['region'] must be defined
        #
        # define the grid type
        grid_type = "generic"
        for gtype in ["equalarea", "gaussian", "generic", "uniform"]:
            if gtype in kwargs['newgrid_name']:
                grid_type = gtype
                break
        try:
            grid_type
        except:
            grid_type = "generic"
        # define resolution (same resolution in lon and lat)
        grid_res = 1.
        for res in ["0.25x0.25deg", "0.5x0.5deg", "0.75x0.75deg", "1x1deg", "1.25x1.25deg", "1.5x1.5deg",
                    "1.75x1.75deg", "2x2deg", "2.25x2.25deg", "2.5x2.5deg", "2.75x2.75deg"]:
            if res in kwargs['newgrid_name']:
                if res == "0.25x0.25deg":
                    grid_res = 0.25
                elif res == "0.5x0.5deg":
                    grid_res = 0.5
                elif res == "0.75x0.75deg":
                    grid_res = 0.75
                elif res == "1x1deg":
                    grid_res = 1.
                elif res == "1.25x1.25deg":
                    grid_res = 1.25
                elif res == "1.5x1.5deg":
                    grid_res = 1.5
                elif res == "1.75x1.75deg":
                    grid_res = 1.75
                elif res == "2x2deg":
                    grid_res = 2.
                elif res == "2.25x2.25deg":
                    grid_res = 2.25
                elif res == "2.5x2.5deg":
                    grid_res = 2.5
                else:
                    grid_res = 2.75
                break
        # define bounds of 'region'
        region_ref = ReferenceRegions(kwargs["region"])
        lat1, lat2 = region_ref["latitude"][0], region_ref["latitude"][1]
        lon1, lon2 = region_ref["longitude"][0], region_ref["longitude"][1]
        # create uniform axis
        nlat = lat2 - lat1
        lat = CDMS2createUniformLatitudeAxis(lat1 + (grid_res / 2.), nlat, grid_res)
        nlon = lon2 - lon1
        lon = CDMS2createUniformLongitudeAxis(lon1 + (grid_res / 2.), nlon, grid_res)
        # create grid
        new_grid = CDMS2createRectGrid(lat, lon, "yx", type=grid_type, mask=None)
        new_grid.id = kwargs["newgrid_name"]
    #
    # regrid
    #
    if regridder == "cdms":
        axis = tab_to_regrid.getAxis(0)
        id_name = deepcopy(axis.id)
        if len(tab_to_regrid.shape) == 3 and (axis.id == "months" or axis.id == "years"):
            axis.id = "time"
            tab_to_regrid.setAxis(0, axis)
        tab_out = tab_to_regrid.regrid(new_grid, missing=missing, order=order, mask=mask, regridTool=regridTool,
                                       regridMethod=regridMethod)
        axis = tab_to_regrid.getAxis(0)
        axis.id = id_name
        tab_to_regrid.setAxis(0, axis)
        if tab_to_regrid.getGrid().shape == new_grid.shape:
            tab_out = MV2masked_where(tab_to_regrid.mask, tab_out)
    else:
        regrid_function = REGRID2horizontal__Horizontal(tab_to_regrid.getGrid(), new_grid)
        tab_out = regrid_function(tab_to_regrid)
    return tab_out


def create_time_axis_not_called_time(tab, new_time_name, frequency):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using moving window average based on 'method'
    #################################################################################

    :param tab: masked_array
        masked_array in which the time axis must be renamed
    :param new_time_name: string
        new name for the time axis
    :param frequency: string
        time frequency of the dataset
        e.g., frequency='monthly'

    :return tab_out: masked_array
        tab with the time axis renamed to new_time_name
    """
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


def detrend_time_series(tab, info, axis=0, method='linear', bp=0):
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

    :return tab_out: array
        detrended data
    :return info: string
        given info plus information about what is done to 'tab'
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
        mean = average_temporal(tab)
        tab_out = MV2array(SCIPYsignal_detrend(tab, axis=axis, type=method, bp=bp))
        tab_out = tab_out + mean
        tab_out = MV2masked_where(mask, tab_out)
        tab_out.setAxisList(axes)
        tab_out.setGrid(grid)
        if method == 'linear':
            info = info + ', time series are linearly detrended'
        else:
            info = info + ', the mean value of the time series is subtracted'
    return tab_out, info


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
    :return num: int
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
    else:
        axis_nick = "no_name"
        axis_nicks = "no_name"
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis name",
            str().ljust(5) + "unknown axis name: " + str(name_axis),
            str().ljust(10) + "known axis name(s): depth, latitude, longitude, time"]
        EnsoErrorsWarnings.my_error(list_strings)
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
            str().ljust(10) + "axes: " + str(tab.getAxisList())]
        EnsoErrorsWarnings.my_error(list_strings)
    return num


def min_and_max(tab):
    """
    #################################################################################
    Description:
    Finds the minimum and maximum values in the given tab
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :return min_max: list of floats
        minimum and maximum values in tab
    """
    return [float(MV2minimum(tab)), float(MV2maximum(tab))]


def normalize_time_series(tab, frequency):
    """
    #################################################################################
    Description:
    normalize time series by the temporal standard deviation (daily std, monthly std, yearly std)
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param frequency: string, optional
        time frequency of the dataset
        e.g., frequency='monthly'

    :return tab_out: masked_array
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
        time_steps_per_year = 1
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    if len(tab) % time_steps_per_year != 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": data length",
            str().ljust(5) + "the normalization function needs only full years: " +
            str(len(tab) // time_steps_per_year) + " years + " + str(len(tab) % time_steps_per_year),
            str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " + str(time_steps_per_year) +
            "), len(dataset) = " + str(len(tab)) + ", so " + str(len(tab) / float(time_steps_per_year)) + " years"]
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
    for yy in range(len(tab) / time_steps_per_year):
        tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
            tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] / std
    if len(tab.shape) == 1:
        tab_out = CDMS2createVariable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
    else:
        axes = axes + tab.getAxisList()[1:]
        grid = tab.getGrid()
        mask = tab.mask
        tab_out = CDMS2createVariable(tab_out, axes=axes, grid=grid, mask=mask, attributes=tab.attributes, id=tab.id)
    return tab_out


def read_and_select_region(filename, var_name, region=None, time_bounds=None, frequency=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given 'var_name' from the given 'filename' and selects the given 'region'

    Uses cdms2 (CDAT) to read 'var_name' from 'filename' and cdutil (CDAT) to select the 'region'
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param var_name: string
        name of the variable to read from 'filename'
    :param region: string
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
        masked_array containing 'var_name' in 'region'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if region is None:  # no region given
        if time_bounds is None: # no time period given
            # read file
            tab = fi(var_name)
        else:  # time period given by the user
            # read file
            tab = fi(var_name, time=time_bounds)
    else:  # region given by the user
        # define region
        region_ref = ReferenceRegions(region)
        if time_bounds is None:  # no time period given
            #  read file
            tab = fi(var_name, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        else:
            # read file
            tab = fi(var_name, time=time_bounds, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
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
        if isinstance(time_bounds[0], str) is True or isinstance(time_bounds[0], basestring) is True:
            if str(tab.getTime().asComponentTime()[0]) < time_bounds[0]:
                tab = tab[1:]
            if str(tab.getTime().asComponentTime()[-1]) > time_bounds[1]:
                tab = tab[:-1]
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


def read_area_and_select_region(filename, area_name='', region=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given areacell from the given 'filename' and selects the given 'region'

    Uses cdms2 (CDAT) to read areacell from 'filename' and cdutil (CDAT) to select the 'region'
    #################################################################################

    :param filename: string
        string of the path to the file and name of the file to read
    :param area_name: string, optional
        name of areacell (areacella, areacello,...) in 'filename'
    :param region: string, optional
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return area: masked_array
        masked_array containing areacell in 'region'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if region is None:  # no region given
        # read file
        try:
            areacell = fi(area_name)
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
    else:  # region given by the user
        # define region
        region_ref = ReferenceRegions(region)
        # read file
        try:
            areacell = fi(area_name, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
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


def read_landmask_and_select_region(tab, filename, landmask_name='', region=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given landmask from the given 'filename' and selects the given 'region'

    Uses cdms2 (CDAT) to read areacell from 'filename' and cdutil (CDAT) to select the 'region'
    #################################################################################

    :param tab: array
        tab of data (sst, taux,...), used to estimate the landmask if it cannot be found in the file
    :param filename: string
        string of the path to the file and name of the file to read
    :param landmask_name: string, optional
        name of landmask (sftlf, lsmask, landmask,...) in 'filename'
    :param region: string, optional
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return area: masked_array
        masked_array containing landmask in 'region'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Get landmask
    if OSpath__isfile(filename):
        # Open file and get time dimension
        fi = CDMS2open(filename)
        if region is None:  # no region given
            # read file
            try:
                landmask = fi(landmask_name)
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
        else:  # region given by the user
            # define region
            region_ref = ReferenceRegions(region)
            # read file
            try:
                landmask = fi(landmask_name, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
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
        landmask = compute_estimate_landmask(tab)
        if region is not None:
            # define region
            region_ref = ReferenceRegions(region)
            # subset
            landmask = landmask(latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    # Return
    return landmask


def save_netcdf(netcdf_name, var1=None, var1_attributes={}, var1_name='', var1_time_name=None, var2=None,
                var2_attributes={}, var2_name='', var2_time_name=None, var3=None, var3_attributes={}, var3_name='',
                var3_time_name=None, var4=None, var4_attributes={}, var4_name='', var4_time_name=None, var5=None,
                var5_attributes={}, var5_name='', var5_time_name=None, var6=None, var6_attributes={}, var6_name='',
                var6_time_name=None, var7=None, var7_attributes={}, var7_name='', var7_time_name=None, var8=None,
                var8_attributes={}, var8_name='', var8_time_name=None, var9=None, var9_attributes={}, var9_name='',
                var9_time_name=None, var10=None, var10_attributes={}, var10_name='', var10_time_name=None, var11=None,
                var11_attributes={}, var11_name='', var11_time_name=None, var12=None, var12_attributes={},
                var12_name='', var12_time_name=None, frequency='monthly', global_attributes={}):
    if OSpath_isdir(ntpath.dirname(netcdf_name)) is not True:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": given path does not exist",
            str().ljust(5) + "netcdf_name = " + str(netcdf_name)
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    if OSpath__isfile(netcdf_name) is True:
        o = CDMS2open(netcdf_name, 'a')
    else:
        o = CDMS2open(netcdf_name, 'w+')
    if var1 is not None:
        if var1_name == '':
            var1_name = var1.id
        if var1_time_name is not None:
            var1 = create_time_axis_not_called_time(var1, var1_time_name, frequency)
        o.write(var1, attributes=var1_attributes, dtype='float32', id=var1_name)
    if var2 is not None:
        if var2_name == '':
            var2_name = var2.id
        if var2_time_name is not None:
            var2 = create_time_axis_not_called_time(var2, var2_time_name, frequency)
        o.write(var2, attributes=var2_attributes, dtype='float32', id=var2_name)
    if var3 is not None:
        if var3_name == '':
            var3_name = var3.id
        if var3_time_name is not None:
            var3 = create_time_axis_not_called_time(var3, var3_time_name, frequency)
        o.write(var3, attributes=var3_attributes, dtype='float32', id=var3_name)
    if var4 is not None:
        if var4_name == '':
            var4_name = var4.id
        if var4_time_name is not None:
            var4 = create_time_axis_not_called_time(var4, var4_time_name, frequency)
        o.write(var4, attributes=var4_attributes, dtype='float32', id=var4_name)
    if var5 is not None:
        if var5_name == '':
            var5_name = var5.id
        if var5_time_name is not None:
            var5 = create_time_axis_not_called_time(var5, var5_time_name, frequency)
        o.write(var5, attributes=var5_attributes, dtype='float32', id=var5_name)
    if var6 is not None:
        if var6_name == '':
            var6_name = var6.id
        if var6_time_name is not None:
            var6 = create_time_axis_not_called_time(var6, var6_time_name, frequency)
        o.write(var6, attributes=var6_attributes, dtype='float32', id=var6_name)
    if var7 is not None:
        if var7_name == '':
            var7_name = var7.id
        if var7_time_name is not None:
            var7 = create_time_axis_not_called_time(var7, var7_time_name, frequency)
        o.write(var7, attributes=var7_attributes, dtype='float32', id=var7_name)
    if var8 is not None:
        if var8_name == '':
            var8_name = var8.id
        if var8_time_name is not None:
            var8 = create_time_axis_not_called_time(var8, var8_time_name, frequency)
        o.write(var8, attributes=var8_attributes, dtype='float32', id=var8_name)
    if var9 is not None:
        if var9_name == '':
            var9_name = var9.id
        if var9_time_name is not None:
            var9 = create_time_axis_not_called_time(var9, var9_time_name, frequency)
        o.write(var9, attributes=var9_attributes, dtype='float32', id=var9_name)
    if var10 is not None:
        if var10_name == '':
            var10_name = var10.id
        if var10_time_name is not None:
            var10 = create_time_axis_not_called_time(var10, var10_time_name, frequency)
        o.write(var10, attributes=var10_attributes, dtype='float32', id=var10_name)
    if var11 is not None:
        if var11_name == '':
            var11_name = var11.id
        if var11_time_name is not None:
            var11 = create_time_axis_not_called_time(var11, var11_time_name, frequency)
        o.write(var11, attributes=var11_attributes, dtype='float32', id=var11_name)
    if var12 is not None:
        if var12_name == '':
            var12_name = var12.id
        if var12_time_name is not None:
            var12 = create_time_axis_not_called_time(var12, var12_time_name, frequency)
        o.write(var12, attributes=var12_attributes, dtype='float32', id=var12_name)
    for att in global_attributes.keys():
        o.__setattr__(att, global_attributes[att])
    o.close()
    return


def smooth_gaussian(tab, axis=0, window=5):
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
    degree = window / 2

    # Create the gaussian weight array
    weightGauss = list()
    for ii in range(window):
        ii = ii - degree + 1
        frac = ii / float(window)
        gauss = float(1 / (NPexp((4 * (frac)) ** 2)))
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
    axes0 = new_tab[(window / 2):len(new_tab) - (window / 2)].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


def smooth_square(tab, axis=0, window=5):
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


def smooth_triangle(tab, axis=0, window=5):
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
    degree = window / 2

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
    axes0 = new_tab[(window / 2):len(new_tab) - (window / 2)].getAxisList()[0]
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


# Dictionary of smoothing methods
dict_smooth = {'gaussian': smooth_gaussian, 'square': smooth_square, 'triangle': smooth_triangle}


def smooth_data(tab, info, axis=0, window=5, method='triangle'):
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
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of often used combinations of previous functions
#
def pre_process_data(tab, info, areacell=None, average=False, interannual_anomalies=False, annual_cycle=False,
                     debug=False, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Combines compute_interannual_anomalies, normalize_time_series, detrend_time_series, smooth_data,
    compute_annual_cycle and dict_average (depending on the needed averaging method)
    Computes interannual anomalies if interannual_anomalies is True, annual cycle if annual_cycle is True, average if
    average is defined (string or list of strings), normalize time series (if applicable), smooth time series (if
    applicable)

    Uses CDAT
    #################################################################################

    :param tab: masked_array
        masked_array of data (sst, taux,...)
    :param info: string
        information about what is done to 'tab'
    :param areacell: masked_array, optional
        areacell defined where 'tab' is defined
        default value is None
    :param average: boolean or string or list, optional
        name or list of names of the average to perform: 'horizontal', 'meridional', 'time', 'zonal'
        default value is False, no average to perform
    :param interannual_anomalies: boolean, optional
        default value is False, interannual anomalies not computed
        If want to compute interannual anomalies set it to True
    :param annual_cycle: boolean, optional
        default value is False, annual cycle not computed
        If want to compute annual cycle set it to True
    :param debug: boolean, optional
        default value is False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param region: string, optional
        name of the region where 'tab' is defined, must be defined in EnsoCollectionsLib.ReferenceRegions

    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False

    :return tab: masked_array
        masked_array pre processed
    :return info: string
        given info plus information about what is done to 'tab'
    """
    # removes annual cycle (anomalies with respect to the annual cycle)
    if interannual_anomalies is True:
        tab = compute_interannual_anomalies(tab)
    # Normalization of the anomalies
    if kwargs['normalization'] is True:
        if kwargs['frequency'] is not None:
            tab = normalize_time_series(tab, kwargs['frequency'])
            info = info + ', normalized'
    # Removing linear trend
    if isinstance(kwargs['detrending'], dict) is True:
        known_args = {'axis', 'method', 'bp'}
        extra_args = set(kwargs['detrending']) - known_args
        if extra_args:
            EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
        tab, info = detrend_time_series(tab, info, **kwargs['detrending'])
    # Smoothing time series
    if isinstance(kwargs['smoothing'], dict) is True:
        known_args = {'axis', 'method', 'window'}
        extra_args = set(kwargs['smoothing']) - known_args
        if extra_args:
            EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
        tab, info = smooth_data(tab, info, **kwargs['smoothing'])
    # computes mean annual cycle
    if annual_cycle is True:
        tab = compute_annual_cycle(tab)
    # average
    if average is not False:
        if debug is True:
            EnsoErrorsWarnings.debug_mode('\033[93m', "EnsoCDATToolsLib PreProcessTS", 20)
            dict_debug = {'axes1':  str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', "averaging to perform: " + str(average), 25, **dict_debug)
        if isinstance(average, str) is True or isinstance(average, basestring) is True:  # only one average
            try: dict_average[average]
            except: EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
            else:
                tab = dict_average[average](tab, areacell, region=region, **kwargs)
                if debug is True:
                    dict_debug = {'axes1': str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                    EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(average), 25, **dict_debug)
        elif isinstance(average, list) is True:  # several averaging steps
            for av in average:
                try: dict_average[av]
                except: EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
                else:
                    tab = dict_average[av](tab, areacell, region=region, **kwargs)
                    if debug is True:
                        dict_debug = {'axes1': str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                        EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(av), 25, **dict_debug)
        else: EnsoErrorsWarnings.unknown_averaging(average, dict_average.keys(), INSPECTstack())
    return tab, info


def read_data(filename, var_name, var_family, metric, region, file_area='', name_area='', file_mask='', name_mask='',
              mask_land=False, mask_ocean=False, time_bounds=None, debug=False, **kwargs):
    """
    #################################################################################
    Description:
    Combines read_select_check and read_area_and_mask
    Reads the given 'var_name' from the given 'filename', selects the given 'region' and checks the 'var_name''s units
    depending on 'var_family', same for areacell, then masks the land if 'mask_land' is True and/or the ocean if
    'mask_ocean' is True on both tab_out (contaning 'var_name') and areacell

    Uses CDAT
    #################################################################################

    :param filename: string
        path_to/filename of the file (NetCDF) of 'var_name'
    :param var_name: string
        name of the variable to read from 'filename'
    :param var_family: string
        family of variable encompassing 'var_name' (temperature, velocity,...)
    :param metric: string
        name of the metric calling this function
    :param region: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param file_area: string
        path_to/filename of the file (NetCDF) of 'name_area'
    :param name_area: string
        name of areacell variable (e.g., areacell, areacella, areacello) in 'file_area'
    :param file_mask: string
        path_to/filename of the file (NetCDF) of 'name_mask'
    :param name_mask: string
        name of landmask variable (e.g., sftlf, lsmask, landmask) in 'file_mask'
    :param mask_land: boolean, optional
        masks land points
        default value is False
    :param mask_ocean: boolean, optional
        masks ocean points
        default value is False
    :param time_bounds: time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param debug: boolean, optional
        default value is False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    usual kwargs:
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None

    :return tab_out: masked_array
        masked_array of 'var_name' in 'region' with masked lands if mask_land is True and/or masked oceans if mask_ocean
        is True
    :return areacell: masked_array
        areacell in 'region' with masked lands if mask_land is True and/or masked oceans if mask_ocean is True
    :return keyerror: string or None
        name of the error(s) encountered if applicable, else None
    """
    keyerror1, keyerror2, keyerror3 = None, None, None
    # Read variable
    if debug is True:
        dict_debug = {'file1': '(' + var_family + ') ' + str(filename), 'var1': '(' + var_family + ') ' + str(var_name)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'Files', 20, **dict_debug)
    tab_out, keyerror1 = read_select_check(filename, var_name, var_family, region=region, time_bounds=time_bounds)
    if debug is True:
        dict_debug = {'axes1': '(' + var_family + ') ' + str([ax.id for ax in tab_out.getAxisList()]),
                      'shape1': '(' + var_family + ') ' + str(tab_out.shape),
                      'time1': '(' + var_family + ') ' + str(time_bounds(tab_out))}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadSelectRegionCheckUnits', 20, **dict_debug)
    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        if len(tab_out) < kwargs['min_time_steps']:
            EnsoErrorsWarnings.too_short_time_period(metric, len(tab_out), kwargs['min_time_steps'], INSPECTstack())
            keyerror2 = "too short time period (" + str(len(tab_out)) + ")"
    # Read areacell & mask
    tab_out, areacell, keyerror3 = read_area_and_mask(tab_out, filename, var_family, region, file_area=file_area,
                                                      name_area=name_area, file_mask=file_mask, name_mask=name_mask,
                                                      mask_land=mask_land, mask_ocean=mask_ocean, debug=debug)
    if keyerror1 is not None or keyerror2 is not None or keyerror3 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2, keyerror3])
    else:
        keyerror = None
    return tab_out, areacell, keyerror


def read_select_check(filename, var_name, var_family, region=None, time_bounds=None, frequency=None, **keyarg):
    """
    #################################################################################
    Description:
    Combines read_and_select_region and check_units
    Reads the given 'var_name' from the given 'filename', selects the given 'region' and checks the 'var_name''s units
    depending on 'var_family'

    Uses CDAT
    #################################################################################

    :param filename: string
        path_to/filename of the file (NetCDF) of 'var_name'
    :param var_name: string
        name of the variable to read from 'filename'
    :param var_family: string
        family of variable encompassing 'var_name' (temperature, velocity,...)
    :param region: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param frequency: string, optional
        time frequency of the dataset
        e.g., frequency='monthly'
        default value is None

    :return tab: masked_array
        masked_array containing 'var_name' in 'region'
    :return keyerror: string or None
        name of the error(s) encountered if applicable, else None
    """
    tab = read_and_select_region(filename, var_name, region=region, time_bounds=time_bounds, frequency=frequency)
    tab, units, keyerror = check_units(tab, var_family, var_name, tab.units, return_tab_only=False)
    tab.name = var_name
    tab.units = units
    return tab, keyerror


def read_area_and_mask(tab, file_data, var_family, region, file_area='', name_area='', file_mask='', name_mask='',
                       mask_land=False, mask_ocean=False, debug=False, **kwargs):
    """
    #################################################################################
    Description:
    Combines read_area_and_select_region, read_landmask_and_select_region, apply_landmask and apply_landmask_to_area
    Reads the given 'name_area' from the given 'file_area' and selects the given 'region', same for landmask, then masks
    the land if 'mask_land' is True and/or the ocean if 'mask_ocean' is True on both tab and areacell

    Uses CDAT
    #################################################################################

    :param tab: masked_array
        masked_array of data (sst, taux,...)
    :param file_data: string
        path_to/filename of the file (NetCDF) used to read 'tab'
    :param var_family: string
        family of variable encompassing 'var_name' (temperature, velocity,...)
    :param region: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param file_area: string
        path_to/filename of the file (NetCDF) of 'name_area'
    :param name_area: string
        name of areacell variable (e.g., areacell, areacella, areacello) in 'file_area'
    :param file_mask: string
        path_to/filename of the file (NetCDF) of 'name_mask'
    :param name_mask: string
        name of landmask variable (e.g., sftlf, lsmask, landmask) in 'file_mask'
    :param mask_land: boolean, optional
        masks land points
        default value is False
    :param mask_ocean: boolean, optional
        masks ocean points
        default value is False
    :param debug: boolean, optional
        default value is False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)

    :return tab_out: masked_array
        tab with masked lands if mask_land is True and/or masked oceans if mask_ocean is True
    :return areacell: masked_array
        areacell in 'region' with masked lands if mask_land is True and/or masked oceans if mask_ocean is True
    :return keyerror: string or None
        name of the error(s) encountered if applicable, else None
    """
    tab_out = deepcopy(tab)
    # Read areacell
    if file_area:
        areacell = read_area_and_select_region(file_area, area_name=name_area, region=region)
    else:
        areacell = read_area_and_select_region(file_data, area_name=name_area, region=region)
    if debug is True:
        if areacell is not None:
            dict_debug = {'axes1': '(' + var_family + ') ' + str([ax.id for ax in areacell.getAxisList()]),
                          'shape1': '(' + var_family + ') ' + str(areacell.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'areacell is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
    # Read landmask
    if file_mask:
        landmask = read_landmask_and_select_region(tab, file_mask, landmask_name=name_mask, region=region)
    else:
        landmask = read_landmask_and_select_region(tab, file_data, landmask_name=name_mask, region=region)
    if debug is True:
        if landmask is not None:
            dict_debug = {'axes1': '(' + var_family + ') ' + str([ax.id for ax in landmask.getAxisList()]),
                          'shape1': '(' + var_family + ') ' + str(landmask.shape)}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'landmask is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
    # Apply landmask
    keyerror1, keyerror2 = None, None
    if landmask is not None:
        tab_out, keyerror1 = apply_landmask(tab_out, landmask, mask_land=mask_land, mask_ocean=mask_ocean)
        if keyerror1 is None:
            if areacell is None:
                areacell = array_ones(landmask, var_id='areacell')
            areacell, keyerror2 = apply_landmask_to_area(areacell, landmask, mask_land=mask_land, mask_ocean=mask_ocean)
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
    return tab_out, areacell, keyerror
