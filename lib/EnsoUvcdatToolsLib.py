# -*- coding:UTF-8 -*-
from calendar import monthrange
from copy import deepcopy
from datetime import date
from inspect import stack as INSPECTstack
from numpy import array as NParray
from numpy import exp as NPexp
from numpy import nonzero as NPnonzero
from scipy.signal import detrend as SCIPYsignal_detrend

# ENSO_metrics package functions:
from EnsoCollectionsLib import CmipVariables
from EnsoCollectionsLib import ReferenceObservations
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings
from EnsoToolsLib import StringInDict

# uvcdat based functions:
import adamsregrid
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
from cdtime import comptime as CDTIMEcomptime
import cdutil
from genutil.statistics import linearregression as GENUTILlinearregression
from genutil.statistics import rms as GENUTILrms
from genutil.statistics import std as GENUTILstd
from MV2 import add as MV2add
from MV2 import arange as MV2arange
from MV2 import array as MV2array
from MV2 import average as MV2average
from MV2 import compress as MV2compress
from MV2 import divide as MV2divide
from MV2 import masked_where as MV2masked_where
from MV2 import minimum as MV2minimum
from MV2 import multiply as MV2multiply
from MV2 import subtract as MV2subtract
from MV2 import sum as MV2sum
from MV2 import take as MV2take
from MV2 import where as MV2where
from MV2 import zeros as MV2zeros


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of simple uvcdat functions used in EnsoMetricsLib.py
#
def AverageHorizontal(tab):
    """
    #################################################################################
    Description:
    Averages along 'xy' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    try: averaged_tab = cdutil.averager(tab, axis='xy')
    except:
        lat_num = get_num_axis(tab, 'latitude')
        lon_num = get_num_axis(tab, 'longitude')
        print " AverageHorizontal"
        print '\033[93m' + str().ljust(15) + "EnsoUvcdatToolsLib AverageHorizontal" + '\033[0m'
        print '\033[93m' + str().ljust(20) + "axes" + str(lat_num)+str(lon_num) + '\033[0m'
        try: averaged_tab = cdutil.averager(tab, axis=str(lat_num)+str(lon_num))
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": horizontal average",
                str().ljust(5) + "cannot perform horizontal average"
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return averaged_tab


def AverageMeridional(tab):
    """
    #################################################################################
    Description:
    Averages along 'y' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    try: averaged_tab = cdutil.averager(tab, axis='y')
    except:
        lat_num = get_num_axis(tab, 'latitude')
        try: averaged_tab = cdutil.averager(tab, axis=str(lat_num))
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": meridional average",
                str().ljust(5) + "cannot perform meridional average"
            ]
            EnsoErrorsWarnings.MyError(list_strings)
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


def AverageTemporal(tab):
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
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": temporal average",
                str().ljust(5) + "cannot perform temporal average"
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return averaged_tab


def AverageZonal(tab):
    """
    #################################################################################
    Description:
    Averages along 'x' axis
    #################################################################################

    for more information:
    import cdutil
    help(cdutil.averager)
    """
    try: averaged_tab = cdutil.averager(tab, axis='x')
    except:
        lon_num = get_num_axis(tab, 'longitude')
        try: averaged_tab = cdutil.averager(tab, axis=str(lon_num))
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": zonal average",
                str().ljust(5) + "cannot perform zonal average"
            ]
            EnsoErrorsWarnings.MyError(list_strings)
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
            EnsoErrorsWarnings.MismatchShapesError(tab, number_or_tab, INSPECTstack())
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
            EnsoErrorsWarnings.MismatchShapesError(tab, number_or_tab, INSPECTstack())
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
            EnsoErrorsWarnings.MismatchShapesError(tab, number_or_tab, INSPECTstack())
    return MV2multiply(tab, number_or_tab)


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
            EnsoErrorsWarnings.MismatchShapesError(tab, number_or_tab, INSPECTstack())
    return MV2subtract(tab, number_or_tab)


# Dictionary of operations
dict_operations = {'divide': OperationDivide, 'minus': OperationSubtract, 'multiply': OperationMultiply,
                   'plus': OperationAdd}


def RmsHorizontal(tab, ref, centered=0):
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
    try: rmse = GENUTILrms(tab, ref, weights='weighted', axis='xy', centered=centered, biased=1)
    except:
        lat_num = get_num_axis(tab, 'latitude')
        lon_num = get_num_axis(tab, 'longitude')
        try: rmse = GENUTILrms(tab, ref, weights='weighted', axis=str(lat_num)+str(lon_num), centered=centered, biased=1)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": horizontal RMS",
                str().ljust(5) + "cannot perform horizontal RMS",
                str().ljust(10) + "either lat and lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat and lon are not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return float(rmse)


def RmsMeridional(tab, ref, centered=0):
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
    try: rmse = GENUTILrms(tab, ref, axis='y', centered=centered, biased=1)
    except:
        lat_num = get_num_axis(tab, 'latitude')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lat_num), centered=centered, biased=1)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": meridional RMS",
                str().ljust(5) + "cannot perform meridional RMS",
                str().ljust(10) + "lat cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return float(rmse)


def RmsTemporal(tab, ref, centered=0):
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
        rmse = GENUTILrms(tab, ref, axis='t', centered=centered, biased=1)
    except:
        time_num = get_num_axis(tab, 'time')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(time_num), centered=centered, biased=1)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": temporal RMS",
                str().ljust(5) + "cannot perform temporal RMS",
                str().ljust(10) + "time cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or time is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return float(rmse)


def RmsZonal(tab, ref, centered=0):
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
        rmse = GENUTILrms(tab, ref, axis='x', centered=centered, biased=1)
    except:
        lon_num = get_num_axis(tab, 'longitude')
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lon_num), centered=centered, biased=1)
        except:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": zonal RMS",
                str().ljust(5) + "cannot perform zonal RMS",
                str().ljust(10) + "lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lon is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList()),
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    return float(rmse)


# Dictionary of RMS methods
dict_rms = {'horizontal': RmsHorizontal, 'meridional': RmsMeridional, 'time': RmsTemporal, 'zonal': RmsZonal}


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
    return float(GENUTILstd(tab, weights=weights, axis=axis, centered=centered, biased=biased))


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
def CheckTime(tab1, tab2, frequency='monthly', min_time_steps=None, metric_name='', **kwargs):
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
    # gets dates of the first and last the time steps of tab1
    stime1 = tab1.getTime().asComponentTime()[0]
    etime1 = tab1.getTime().asComponentTime()[-1]

    # gets dates of the first and last the time steps of tab2
    stime2 = tab2.getTime().asComponentTime()[0]
    etime2 = tab2.getTime().asComponentTime()[-1]

    # retains only the latest start date and the earliest end date
    tmp = stime1.cmp(stime2)
    # tmp is -1 when time1 < time2
    # tmp is  0 when time1 = time2
    # tmp is  1 when time1 > time2
    if tmp in [0, -1]:
        stime = stime2
    else:
        stime = stime1
    tmp = etime1.cmp(etime2)
    if tmp in [0, -1]:
        etime = etime1
    else:
        etime = etime2

    # defines the period between the two dates
    if frequency == 'daily':
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime.day, 23, 59, 60.0)
    elif frequency == 'monthly':
        etime_day = monthrange(etime.year, etime.month)[-1]
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime_day, 23, 59, 60.0)
    elif frequency == 'yearly':
        stime_adjust = CDTIMEcomptime(stime.year, 1, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, 12, 31, 23, 59, 60.0)
    else:
        EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())

    # retains only the time-period common to both tab1 and tab2
    tab1_sliced = tab1(time=(stime_adjust, etime_adjust))
    tab2_sliced = tab2(time=(stime_adjust, etime_adjust))

    # checks if the remaining time-period fulfills the minimum length criterion
    if min_time_steps is not None:
        if len(tab1_sliced)<min_time_steps or len(tab2_sliced)<min_time_steps:
            shortest = min(len(tab1_sliced), len(tab2_sliced))
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, shortest, min_time_steps, INSPECTstack())
    return tab1_sliced, tab2_sliced


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
    if var_name in ['temperature']:
        if units == 'K':
            # check if the temperature units is really K
            if float(MV2minimum(tab)) > 200:
                # unit change of the temperature: from K to degC
                tab = dict_operations['minus'](tab, 273.15)
                units = "degC"
            else:
                EnsoErrorsWarnings.UnlikelyUnits(var_name, name_in_file, units, INSPECTstack())
        elif units in ['C', 'degree_Celsius', 'deg_Celsius', 'deg. C', 'degCelsius', 'degree_C', 'deg_C', 'degC',
                       'degrees C']:
            # check if the temperature units is really degC
            if float(MV2minimum(tab)) < 50:
                units = "degC"
            else:
                EnsoErrorsWarnings.UnlikelyUnits(var_name, name_in_file, units, INSPECTstack())
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['precipitations']:
        if units == 'kg m-2 s-1':
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = dict_operations['multiply'](tab, 86400)
        elif units == 'mm/day':
            pass
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['wind stress']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "N/m2"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['velocity']:
        if units in ['cm s-1', 'cm/s', 'cm s**-1']:
            # unit change of the velocity: from cm/s to m/s
            tab = dict_operations['multiply'](tab, 1e-2)
            units = "m/s"
        elif units in ['m s-1', 'm/s', 'm s**-1', 'm/sec']:
            units = "m/s"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['heat flux']:
        if units in ['W/m2', 'W m-2', 'W/m^2']:
            units = "W/m2"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['pressure']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "Pa"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.MyWarning(list_strings)
    if return_tab_only:
       return tab
    else:
        return tab, units


def Composite(tab, list_event_years, frequency, nbr_years_window=None):
    # function to fill array with masked value where the data is not available
    def fill_array(tab, units, freq):
        y1 = tab.getTime().asComponentTime()[0].year
        m1 = tab.getTime().asComponentTime()[0].month
        d1 = tab.getTime().asComponentTime()[0].day
        tab_out = MV2zeros(nbr_years_window * 12)
        tab_out = MV2masked_where(tab_out == 0, tab_out)
        axis = CDMS2createAxis(MV2array(range(len(tab_out)), dtype='int32'), id='time')
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
            yy1, yy2 = yy - ((nbr_years_window / 2) - 1), yy + (nbr_years_window / 2)
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
        composite = MV2average(MV2array(composite), axis=0)
        # axis list
        axis0 = CDMS2createAxis(MV2array(range(len(composite)), dtype='int32'), id='time')
        axis0.units = units_out
        if tab.shape > 1:
            axes = [axis0] + tab.getAxisList()[1:]
        else:
            axes = [axis0]
        composite.setAxisList(axes)
    else:
        composite = MV2average(tab, axis=0)
    return composite


def DetectEvents(tab, season, threshold, normalization=False, nino=True):
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
    tab = SeasonalMean(tab, season, compute_anom=True)
    if season == 'DJF':
        time_ax = tab.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
    # Normalization ?
    if normalization:
        threshold = threshold * float(GENUTILstd(tab, weights=None, axis=0, centered=1, biased=1))
    # Initialization
    tab_threshold = MV2zeros(tab.shape)
    tab_threshold.fill(threshold)
    list_years = [tab.getTime().asComponentTime()[yy].year for yy in range(len(tab))]
    indices = MV2arange(len(tab))
    # Conditions
    if nino:
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
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": method",
            str().ljust(5) + "unknown method: " + str(method)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
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
        e.g., frequency='latitude'
    :return number: int
        position of the axis named "name_axis"
    """
    num = None
    if name_axis == 'depth':
        axis_nick = 'lev'
        axis_nicks = ['z','Z']
    if name_axis == 'latitude':
        axis_nick = 'lat'
        axis_nicks = ['j','y','Y']
    elif name_axis == 'longitude':
        axis_nick = 'lon'
        axis_nicks = ['i','x','X']
    elif name_axis == 'time':
        axis_nick = 'time'
        axis_nicks = ['t','T']
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
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "cannot find axis named: " + str(name_axis)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    return num


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
    if frequency == 'daily':
        time_steps_per_year = 365
    elif frequency == 'monthly':
        time_steps_per_year = 12
    elif frequency == 'yearly':
        time_steps_per_year = 1
    else:
        EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())
    if len(tab) % time_steps_per_year != 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": data length",
            str().ljust(5) + "the normalization function needs only full years: " +
            str(len(tab) // time_steps_per_year) +" years + " + str(len(tab) % time_steps_per_year),
            str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " + str(time_steps_per_year) +
            "), len(dataset) = " + str(len(tab)) + ", so " + str(len(tab) / float(time_steps_per_year)) + " years",
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    # reshape tab like [yy,nb]
    new_tab = list()
    for yy in range(len(tab)/time_steps_per_year):
        new_tab.append(tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year])
    new_tab = MV2array(new_tab)
    std = MV2zeros(new_tab[0].shape)
    for dd in range(time_steps_per_year):
        std[dd] = float(GENUTILstd(new_tab[:,dd], weights=None, axis=0, centered=1, biased=1))
    for yy in range(len(tab) / time_steps_per_year):
        tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
            tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] / std
    return tab


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
        if time_bounds is None: # no time period given
            # read file
            tab = fi(varname)
        else:  # time period given by the user
            # read file
            tab = fi(varname, time=time_bounds)
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        nbox = cdutil.region.domain(latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        if time_bounds is None:  # no time period given
            #  read file
#            tab = fi(varname, nbox)
            tab = fi(varname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        else:
            # read file
#            tab = fi(varname, nbox, time=time_bounds)
            tab = fi(varname, time=time_bounds, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    fi.close()
    if time_bounds is not None:
        # sometimes the time boundaries are wrong, even with 'time=time_bounds'
        # this section checks if one time step has not been included by error at the beginning or the end of the time
        # series
        if isinstance(time_bounds[0], basestring):
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
        EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())
    # remove axis 'level' if its length is 1
    if tab.getLevel():
        if len(tab.getLevel()) == 1:
            tab = tab(squeeze=1)
    return tab


def Regrid(tab_to_regrid, newgrid, missing=None, order=None, mask=None, regridTool='esmf', regridMethod='linear'):
    """
    #################################################################################
    Description:
    Regrids 'tab_to_regrid' to 'togrid'
    #################################################################################

    for more information:
    import cdms2
    help(cdms2.avariable)
    or (if the case of regridTool='adamsregrid'
    import adamsregrid
    To obtain a prescription for making an instance, type
        adamsregrid.help('Regrid')
    To acquire instructions on the use of the rgrd function, type
        adamsregrid.help('rgrd')
    To look at a general one dimensional example, type
        adamsregrid.help('OneDexample')
    To look at a general four dimensional example, type
        adamsregrid.help('FourDexample')

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
    :param regridTool: string, optional
        regridding tools (either 'regrid2', 'esmf', 'libcf', 'adamsregrid')
        default value is 'esmf'
    :param regridMethod: string, optional
        regridding methods:
            regridTool='regrid2'     -> 'linear'
            regridTool='esmf'        -> 'conserve', 'linear', 'patch'
            regridTool='libcf'       -> 'linear'
            regridTool='adamsregrid' -> 'linear', 'linearLog', 'cubic', 'cubicLog'
        default value is 'linear'
    :return new_tab: masked_array
        tab_to_regrid regridded on newgrid
    """
    if regridTool not in ['regrid2', 'esmf', 'libcf', 'adamsregrid']:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridTool",
            str().ljust(5) + "unknown regridTool: " + str(regridTool),
            str().ljust(10) + "known regridTool: " + str(['regrid2', 'esmf', 'libcf'])
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if (regridTool == 'esmf' and regridMethod not in ['conserve', 'patch', 'linear']) or \
            (regridTool in ['regrid2', 'libcf'] and regridMethod != 'linear'):
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridMethod",
            str().ljust(5) + "unknown regridMethod for this regridTool (" + str(regridMethod) + "): "
            + str(regridMethod),
            str().ljust(10) + "known regridTool: " + str(['conserve', 'patch', 'linear'])
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if (regridTool == 'adamsregrid' and regridMethod not in ['linear', 'linearLog', 'cubic', 'cubicLog']):
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridMethod",
            str().ljust(5) + "unknown regridMethod for this regridTool (" + str(regridMethod) + "): "
            + str(regridMethod),
            str().ljust(10) + "known regridTool: " + str(['conserve', 'patch', 'linear'])
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if regridTool in ['regrid2', 'esmf', 'libcf']:
        new_tab = tab_to_regrid.regrid(newgrid, missing=missing, order=order, mask=mask, regridTool=regridTool,
                                       regridMethod=regridMethod)
    else:
        axis = tab_to_regrid.getAxis(0)[:]
        if axis[0] > newgrid.getAxis(0)[0]:
            axis[0] = newgrid.getAxis(0)[0]
        if axis[-1] < newgrid.getAxis(0)[-1]:
            axis[-1] = newgrid.getAxis(0)[-1]
        if len(tab_to_regrid.shape) == 1:
            r = adamsregrid.Regrid(tab_to_regrid.getAxis(0)[:], newgrid.getAxis(0)[:], regridMethod, 0)
        elif len(tab_to_regrid.shape) == 2:
            r = adamsregrid.Regrid(tab_to_regrid.getAxis(0)[:], newgrid.getAxis(0)[:], regridMethod, 0,
                                   tab_to_regrid.getAxis(1)[:], newgrid.getAxis(1)[:], regridMethod, 1)
        elif len(tab_to_regrid.shape) == 3:
            r = adamsregrid.Regrid(tab_to_regrid.getAxis(0)[:], newgrid.getAxis(0)[:], regridMethod, 0,
                                   tab_to_regrid.getAxis(1)[:], newgrid.getAxis(1)[:], regridMethod, 1,
                                   tab_to_regrid.getAxis(2)[:], newgrid.getAxis(2)[:], regridMethod, 2)
        else:
            r = adamsregrid.Regrid(tab_to_regrid.getAxis(0)[:], newgrid.getAxis(0)[:], regridMethod, 0,
                                   tab_to_regrid.getAxis(1)[:], newgrid.getAxis(1)[:], regridMethod, 1,
                                   tab_to_regrid.getAxis(2)[:], newgrid.getAxis(2)[:], regridMethod, 2,
                                   tab_to_regrid.getAxis(3)[:], newgrid.getAxis(3)[:], regridMethod, 3)
        new_tab = MV2array(r.rgrd(tab_to_regrid))
        new_tab.setAxisList(newgrid.getAxisList())
    return new_tab


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
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
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
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
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
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        EnsoErrorsWarnings.MyError(list_strings)
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
        tab_sea = sea_dict[season]
    except:
        list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": season",
                        str().ljust(5) + "unknown season: " + str(season)]
        EnsoErrorsWarnings.MyError(list_strings)
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
    return tab


# Dictionary of smoothing methods
dict_smooth = {'gaussian': SmoothGaussian, 'square': SmoothSquare, 'triangle': SmoothTriangle}


def Smoothing(tab, info, axis=0, window=5, method='triangle'):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using moving window average based on 'method'
    #################################################################################

    :param tab:
    :param axis:
    :param window:
    :return smoothed_tab: masked_array
        smoothed data
    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the triangle moving window average
        default value is 5
    :param method: string, optional
        smoothing method:
            'gaussian': gaussian shaped window
            'square':   square shaped window
            'triangle': triangle shaped window
    :return:
    """
    try: dict_smooth[method]
    except:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing method (running mean)",
            str().ljust(5) + "unkwown smoothing method: " + str(method),
            str().ljust(10) + "known smoothing method: " + str(sorted(dict_smooth.keys())),
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    info = info + ', smoothing using a ' + str(method) + ' shaped window of ' + str(window) + ' points'
    return dict_smooth[method](tab, axis=axis, window=window), info
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Cet of often used combinations of previous functions
#
def CustomLinearRegression(y, x, sign_x=0, return_stderr=True):
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
    :return slope, stderr: floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    x = NParray(x)
    y = NParray(y)
    if sign_x == 1:
        idxpos = NPnonzero(x >= 0.)
        results = GENUTILlinearregression(y[idxpos], x=x[idxpos], error=1, nointercept=1)
    elif sign_x == -1:
        idxneg = NPnonzero(x <= 0.)
        results = GENUTILlinearregression(y[idxneg], x=x[idxneg], error=1, nointercept=1)
    else:
        results = GENUTILlinearregression(y, x=x, error=1, nointercept=1)
    slope, stderr = results
    if return_stderr:
        return float(slope), float(stderr)
    else:
        return float(slope)


def MyDerive(project, internal_variable_name, dict_var):
    # test input parameters
    if not isinstance(project, basestring):
        EnsoErrorsWarnings.ObjectTypeError('project', 'string', type(project), INSPECTstack())
    if not isinstance(internal_variable_name, basestring):
        EnsoErrorsWarnings.ObjectTypeError('internal_variable_name', 'string', type(internal_variable_name),
                                           INSPECTstack())
    if not isinstance(dict_var, dict):
        EnsoErrorsWarnings.ObjectTypeError('project', 'dictionary', type(dict_var), INSPECTstack())

    # get dictionary of observations
    dict_obs = ReferenceObservations(False)

    # compute 'internal_variable_name' in 'CMIP' case
    if 'CMIP' in project:
        # get dictionary of CMIP
        dict_CMIP = CmipVariables()['variable_name_in_file']
        # test if 'internal_variable_name' is defined in EnsoCollectionsLib.CmipVariables
        if internal_variable_name in dict_CMIP.keys():
            list_var = dict_CMIP[internal_variable_name]['var_name']
            # test if keys in list_var are in 'dict_var'
            StringInDict(list_var, dict_var, INSPECTstack())
            if isinstance(list_var, basestring):
                # this 'internal_variable_name' is based on one variable
                return dict_var[list_var]
            else:
                # this 'internal_variable_name' is based on several variables
                list_operator = dict_CMIP[internal_variable_name]['algebric_calculation']
                if len(list_operator) != len(list_var):
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                        ": variable definition in EnsoCollectionsLib.CmipVariables",
                        str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                        str(internal_variable_name) + " but " + str(len(list_operator)) + " operator(s) are given"
                    ]
                    EnsoErrorsWarnings.MyError(list_strings)
                # compute the output variable
                if list_operator[0] == 'minus':
                    outvar = -1 * dict_var[list_var[0]]
                else:
                    outvar = dict_var[list_var[0]]
                for ii in range(1,len(list_var)):
                    var, operator = list_var[ii], list_operator[ii]
                    outvar = dict_operations[operator](outvar, dict_var[var])
                outvar.setAxisList(dict_var[list_var[0]].getAxisList())
                outvar = MV2masked_where(dict_var[list_var[0]].mask, outvar)
                outvar.setGrid(dict_var[list_var[0]].getGrid())
                return outvar
    # compute 'internal_variable_name' in 'obs' case
    elif project in dict_obs.keys():
        # 'project' is defined in EnsoCollectionsLib.ReferenceObservations
        dict_obs_var = dict_obs[project]['variable_name_in_file']
        # test if 'internal_variable_name' is defined for this observations dataset
        if internal_variable_name in dict_obs_var.keys():
            list_var = dict_obs_var[internal_variable_name]['var_name']
            # test if keys in list_var are in 'dict_var'
            StringInDict(list_var, dict_var, INSPECTstack())
            if isinstance(list_var, basestring):
                # this 'internal_variable_name' is based on one variable
                return dict_var[list_var]
            else:
                # this 'internal_variable_name' is based on several variables
                list_operator = dict_obs_var[internal_variable_name]['algebric_calculation']
                if len(list_operator) != len(list_var):
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                        ": variable definition in EnsoCollectionsLib.ReferenceObservations(" + str(project) + ")",
                        str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                        str(internal_variable_name) + " but " + str(len(list_operator)) + " operator(s) are given"
                    ]
                    EnsoErrorsWarnings.MyError(list_strings)
                # compute the output variable
                if list_operator[0] == 'minus':
                    outvar = -1 * dict_var[list_var[0]]
                else:
                    outvar = dict_var[list_var[0]]
                for ii in range(1, len(list_var)):
                    var, operator = list_var[ii], list_operator[ii]
                    outvar = dict_operations[operator](outvar, dict_var[var])
                outvar.setAxisList(dict_var[list_var[0]].getAxisList())
                outvar = MV2masked_where(dict_var[list_var[0]].mask, outvar)
                outvar.setGrid(dict_var[list_var[0]].getGrid())
                return outvar
    else:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
            ": project",
            str().ljust(5) + "unknown 'project' (or observations dataset): " + str(project),
            str().ljust(10) + "it must be either a 'CMIP' project or an observations dataset defined in\
            EnsoCollectionsLib.ReferenceObservations",
            str().ljust(10) + "known observations dataset: " + str(sorted(dict_obs.keys()))
        ]
        EnsoErrorsWarnings.MyError(list_strings)


def LinearRegressionAndNonlinearity(y, x, return_stderr=True):
    """
    #################################################################################
    Description:
    UvcdatCustomLinearRegression applied for all values of x, for values of x>=0, for values of x<=0

    Uses uvcdat
    #################################################################################

    :param y: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param x: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :return [slope_all_values, stderr_all_values], [slope_positive_values, stderr_positive_values],
            [slope_negative_values, stderr_negative_values]: lists of floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    # all points
    all_values = CustomLinearRegression(y, x, 0, return_stderr=return_stderr)
    # positive SSTA = El Nino
    positive_values = CustomLinearRegression(y, x, 1, return_stderr=return_stderr)
    # negative SSTA = La Nina
    negative_values = CustomLinearRegression(y, x, -1, return_stderr=return_stderr)
    if return_stderr:
        return [float(all_values[0]), float(all_values[1])], [float(positive_values[0]), float(positive_values[1])],\
               [float(negative_values[0]), float(negative_values[1])]
    else:
        return [float(all_values)], [float(positive_values)], [float(negative_values)]


def PreProcessTS(tab, info, average=False, compute_anom=False, **kwargs):
    # removes annual cycle (anomalies with respect to the annual cycle)
    if compute_anom:
        tab = cdutil.ANNUALCYCLE.departures(tab)
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
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        tab, info = Detrend(tab, info, **kwargs['detrending'])
    # Smoothing time series
    if isinstance(kwargs['smoothing'], dict):
        known_args = {'axis', 'method', 'window'}
        extra_args = set(kwargs['smoothing']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        tab, info = Smoothing(tab, info, **kwargs['smoothing'])
    # average
    if average is not False:
        print '\033[93m' + str().ljust(15) + "EnsoUvcdatToolsLib PreProcessTS" + '\033[0m'
        print '\033[93m' + str().ljust(20) + "averaging to perform: " + str(average) + '\033[0m'
        print '\033[93m' + str().ljust(25) + "tab.shape = " + str(tab.shape) + '\033[0m'
        print '\033[93m' + str().ljust(25) + "tab.axes = " + str([ax.id for ax in tab.getAxisList()]) + '\033[0m'
        if isinstance(average, basestring):
            try: dict_average[average]
            except:
                EnsoErrorsWarnings.UnknownAveraging(average, dict_average.keys(), INSPECTstack())
            else:
                tab = dict_average[average](tab)
        elif isinstance(average, list):
            for av in average:
                print '\033[93m' + str().ljust(20) + "perform " + str(av) + '\033[0m'
                try: dict_average[av]
                except:
                    EnsoErrorsWarnings.UnknownAveraging(average, dict_average.keys(), INSPECTstack())
                else:
                    tab = dict_average[av](tab)
                    print '\033[93m' + str().ljust(25) + "tab.shape = " + str(tab.shape) + '\033[0m'
                    print '\033[93m' + str().ljust(25) + "tab.axes = "\
                          + str([ax.id for ax in tab.getAxisList()]) + '\033[0m'
        else:
            EnsoErrorsWarnings.UnknownAveraging(average, dict_average.keys(), INSPECTstack())
    return tab, info


def ReadSelectRegionCheckUnits(filename, varname, varfamily, box=None, time_bounds=None, frequency=None, **keyarg):
    """
    #################################################################################
    Description:
    Combines UvcdatReadAndSelectRegion and UvcdatCheckUnits
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
    tab, units = CheckUnits(tab, varfamily, varname, tab.units, return_tab_only=False)
    tab.name = varname
    tab.units = units
    return tab


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


def TwoVarRegrid(model, obs, info, model_to_obs=True, obs_to_model=False, model_and_obs_to_newgrid = False,
                 newgrid=None, **keyarg):
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
    :param model_to_obs: boolean, optional
        True if you want to regrid model data toward observations data, if you don't want it pass anything but True
    :param obs_to_model: boolean, optional
        True if you want to regrid observations data toward model data, if you don't want it pass anything but True
    :param model_and_obs_to_newgrid: boolean, optional
        True if you want to regrid model data and observations data toward 'newgrid', if you don't want it pass anything
        but True
    :param newgrid: CDMS grid
        grid toward which model data and observations data are regridded if model_and_obs_to_newgrid=True
    :param keyarg:
        see EnsoUvcdatToolsLib.Regrid for regridding options
    :return: model, obs, info
        model and obs on the same grid, and information about what has been done to 'model' and 'obs'
    """
    known_args = {'missing', 'order', 'mask', 'regridTool', 'regridMethod'}
    extra_args = set(keyarg) - known_args
    if extra_args:
        EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
    # if regridTool='adamsregrid', the grid is not sent but the masked_array
    if keyarg['regridTool'] == 'adamsregrid':
        grid_obs = deepcopy(obs)
        grid_model = deepcopy(model)
    else:
        grid_obs = obs.getGrid()
        grid_model = model.getGrid()
    # select case:
    if model_to_obs:
        model = Regrid(model, grid_obs, **keyarg)
        info = info + ', model regridded to observations'
    elif obs_to_model:
        obs = Regrid(obs, grid_model, **keyarg)
        info = info + ', observations regridded to model'
    elif model_and_obs_to_newgrid:
        model = Regrid(model, newgrid, **keyarg)
        obs = Regrid(obs, newgrid, **keyarg)
        try: grid_name = newgrid.id
        except:
            try: grid_name = newgrid.name
            except: grid_name = 'newgrid'
        info = info + ', observations and model regridded to ' + str(grid_name)
    if model.shape == obs.shape:
        mask = model.mask
        mask = MV2where(obs.mask, obs.mask, mask)
        model = MV2masked_where(mask, model)
        obs = MV2masked_where(mask, obs)
    else:
        tab = MV2zeros(model.shape)
        for tt in range(len(tab)):
            tab[tt] = MV2masked_where(obs[0].mask, tab[tt])
        model = MV2masked_where(tab.mask, model)
        tab = MV2zeros(obs.shape)
        for tt in range(len(tab)):
            tab[tt] = MV2masked_where(model[0].mask, tab[tt])
        obs = MV2masked_where(tab.mask, obs)
    return model, obs, info
