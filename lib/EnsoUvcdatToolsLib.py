# -*- coding:UTF-8 -*-
from calendar import monthrange
from inspect import stack as INSPECTstack
from numpy import array as NParray
from numpy import nonzero as NPnonzero

# ENSO_metrics package functions:
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings

# uvcdat based functions:
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
from cdtime import comptime as CDTIMEcomptime
import cdutil
from genutil.statistics import linearregression as GENUTILlinearregression
from genutil.statistics import rms as GENUTILrms
from genutil.statistics import std as GENUTILstd
from MV2 import minimum as MV2minimum
from MV2 import subtract as MV2subtract
from MV2 import multiply as MV2multiply


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of simple uvcdat functions used in EnsoMetricsLib.py
#
def Multiply(tab, number):
    """
    #################################################################################
    Description:
    Multiplies every elements of 'tab' by 'number'
    #################################################################################

    for more information:
    import MV2
    help(MV2.multiply)
    """
    return MV2multiply(tab, number)


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
# Set of more complex uvcdat functions used in EnsoMetricsLib.py
#
def CheckTime(tab1, tab2, frequency=None, minimum_length=None, metric_name=''):
    # gets dates of the first and last the time steps of tab1
    stime1 = tab1.getTime().asComponentTime()[0]
    etime1 = tab1.getTime().asComponentTime()[-1]

    # gets dates of the first and last the time steps of tab2
    stime2 = tab2.getTime().asComponentTime()[0]
    etime2 = tab2.getTime().asComponentTime()[-1]

    # retains only the latest start date and the earliest end date
    stime = max(stime1, stime2)
    etime = min(etime1, etime2)

    # defines the period between the two dates
    if frequency == 'daily':
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, 0, 0, 0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime.day, 23, 59, 59.99)
    elif frequency == 'monthly':
        etime_day = monthrange(etime.year, etime.month)[-1]
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, 0, 0, 0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime_day, 23, 59, 59.99)
    elif frequency == 'yearly':
        stime_adjust = CDTIMEcomptime(stime.year, 1, 1, 0, 0, 0)
        etime_adjust = CDTIMEcomptime(etime.year, 12, 31, 23, 59, 59.99)
    else:
        if frequency is None:
            pass
        else:
            EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, stime.hour, stime.minute, stime.second)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime.day, etime.hour, etime.minute, etime.second)

    # retains only the time-period common to both tab1 and tab2
    tab1_sliced = tab1(time=(stime_adjust, etime_adjust))
    tab2_sliced = tab2(time=(stime_adjust, etime_adjust))

    # checks if the remaining time-period fulfills the minimum length criterion
    if minimum_length is not None:
        if len(tab1_sliced)<minimum_length or len(tab2_sliced)<minimum_length:
            shortest = min(len(tab1_sliced), len(tab2_sliced))
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, shortest, minimum_length, INSPECTstack())
    return tab1_sliced, tab2_sliced


def UvcdatCheckUnits(tab, var_name, name_in_file, units, return_tab_only=True):
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
                tab = MV2subtract(tab, 273.15)
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
            tab = Multiply(tab, 86400)
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
            tab = Multiply(tab, 1e-2)
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


def UvcdatReadAndSelectRegion(filename, varname, box=None, time_bnds=None, frequency=None):
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
    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if box is None:  # no box given
        if time_bnds is None: # no time period given
            # read file
            tab = fi(varname)
        else:  # time period given by the user
            # read file
            tab = fi(varname, time=time_bnds)
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        nbox = cdutil.region.domain(latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        if time_bnds is None:  # no time period given
            #  read file
            tab = fi(varname, nbox)
        else:
            # read file
            tab = fi(varname, nbox, time=time_bnds)
    fi.close()
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
    return tab


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
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Cet of often used combinations of previous functions
#
def ReadSelectRegionCheckUnits(filename, varname, varfamily, box=None, time_bnds=None, frequency=None):
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
    :param box: string
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions
    :param varfamily: string
        family of variable encompassing 'varname' (temperature, velocity,...)
    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    tab = UvcdatReadAndSelectRegion(filename, varname, box=box, time_bnds=time_bnds, frequency=frequency)
    tab, units = UvcdatCheckUnits(tab, varfamily, varname, tab.units, return_tab_only=False)
    tab.name = varname
    tab.units = units
    return tab


def TimeAnomaliesStd(tab):
    """
    #################################################################################
    Description:
    Combines cdutil.averager, cdutil.ANNUALCYCLE.departures and genutil.std
    Averages spatially, removes the annual cycle (anomalies with respect to the annual cycle) and computes the standard
    deviation

    Uses  uvcdat
    #################################################################################

    :param tab: masked_array
        masked_array (uvcdat cdms2) containing a variable, with many attributes attached (short_name, units,...)
    :return std: float
        standard deviation (one value) of the masked_array averaged spatially and with the annual cycle removed
    """
    # horizontal average
    tab = cdutil.averager(tab, axis='xy')
    # removes annual cycle (anomalies with respect to the annual cycle)
    tab = cdutil.ANNUALCYCLE.departures(tab)
    # computes standard deviation
    std = float(GENUTILstd(tab, weights=None, axis=0, centered=1, biased=1))
    return std


def UvcdatCustomLinearRegression(y, x, sign_x=0, return_stderr=True):
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


def UvcdatLinearRegressionAndNonlinearity(y, x, return_stderr=True):
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
    all_values = UvcdatCustomLinearRegression(y, x, 0, return_stderr=return_stderr)
    # positive SSTA = El Nino
    positive_values = UvcdatCustomLinearRegression(y, x, 1, return_stderr=return_stderr)
    # negative SSTA = La Nina
    negative_values = UvcdatCustomLinearRegression(y, x, -1, return_stderr=return_stderr)
    if return_stderr:
        return [float(all_values[0]), float(all_values[1])], [float(positive_values[0]), float(positive_values[1])],\
               [float(negative_values[0]), float(negative_values[1])]
    else:
        return [float(all_values)], [float(positive_values)], [float(negative_values)]


def TimeAnomaliesLinearRegressionAndNonlinearity(tab2, tab1, return_stderr=True):
    """
    #################################################################################
    Description:
    UvcdatLinearRegressionAndNonlinearity applied on two 'raw' masked_arrays (i.e., the annual cycle is not removed and
    the spatial average is not computed)
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
    tab1 = cdutil.averager(tab1, axis='xy')
    tab2 = cdutil.averager(tab2, axis='xy')
    # removes annual cycle (anomalies with respect to the annual cycle)
    tab1 = cdutil.ANNUALCYCLE.departures(tab1)
    tab2 = cdutil.ANNUALCYCLE.departures(tab2)
    # computes linear regression of tab2 on tab1 for all values of tab1, for values of tab1>=0, for values of tab1<=0
    lr, lrpos, lrneg = UvcdatLinearRegressionAndNonlinearity(tab2, tab1, return_stderr=return_stderr)
    return lr, lrpos, lrneg


def SpatialRms(tab, ref, centered=0, regrid_tab_on_ref=True):
    """
    #################################################################################
    Description:
    UvcdatRms applied on two 'raw' masked_arrays (i.e., the annual cycle is not removed, the temporal average is not
    computed and the grids may be different)

    Averages temporally, regrids tab on ref grid (if applicable) and computes the root mean square difference between
    tab and ref

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
    :param regrid_tab_on_ref: boolean, optional
        default value = True, regrids tab on ref grid
        True if you want to regrid, if you don't want it pass anything but true
    :param return_stderr: boolean, optional
        default value = True, returns the the unadjusted standard error
        True if you want the unadjusted standard error, if you don't want it pass anything but true
    :return: [slope_all_values, stderr_all_values], [slope_positive_values, stderr_positive_values],
            [slope_negative_values, stderr_negative_values]: lists of floats
        slope of the linear regression of y over x
        unadjusted standard error of the linear regression of y over x (if return_stderr=True)
    """
    # Time average
    tab = cdutil.averager(tab, axis='t')
    ref = cdutil.averager(ref, axis='t')
    # Regrid model tab on ref grid
    if regrid_tab_on_ref:
        tab = tab.regrid(ref.getGrid(), regridTool='regrid2')
    # Computes the root mean square difference
    rmse = GENUTILrms(tab, ref, weights='weighted', axis='xy', centered=centered, biased=1)
    return float(rmse)
