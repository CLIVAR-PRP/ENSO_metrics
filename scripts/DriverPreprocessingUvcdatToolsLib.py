# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from scipy.signal import detrend as SCIPYsignal_detrend

# ENSO_metrics package functions:
from lib.EnsoCollectionsLib import ReferenceRegions
import lib.EnsoErrorsWarnings

# uvcdat based functions:
import adamsregrid
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import createUniformLatitudeAxis as CDMS2createUniformLatitudeAxis
from cdms2 import createUniformLongitudeAxis as CDMS2createUniformLongitudeAxis
from cdms2 import createRectGrid as CDMS2createRectGrid
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
import cdutil
from genutil.statistics import std as GENUTILstd
from MV2 import array as MV2array
from MV2 import masked_where as MV2masked_where
from MV2 import maximum as MV2maximum
from MV2 import minimum as MV2minimum
from MV2 import sum as MV2sum
from MV2 import zeros as MV2zeros


# ---------------------------------------------------------------------------------------------------------------------#
#
# uvcdat averaging functions
#
def AverageHorizontal(tab, info):
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
        lat_num = GetNumAxis(tab, 'latitude')
        lon_num = GetNumAxis(tab, 'longitude')
        try: averaged_tab = cdutil.averager(tab, axis=str(lat_num)+str(lon_num))
        except:
            list_strings = [
                "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": horizontal average",
                str().ljust(5) + "cannot perform horizontal average"
            ]
            lib.EnsoErrorsWarnings.MyError(list_strings)
    info = info + '\n' + str().ljust(5) + "horizontal average: cdutil.averager"
    return averaged_tab, info


def AverageMeridional(tab, info):
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
        lat_num = GetNumAxis(tab, 'latitude')
        try: averaged_tab = cdutil.averager(tab, axis=str(lat_num))
        except:
            list_strings = [
                "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": meridional average",
                str().ljust(5) + "cannot perform meridional average"
            ]
            lib.EnsoErrorsWarnings.MyError(list_strings)
    lon = tab.getLongitude()
    if len(lon.shape) > 1:
        lonn = CDMS2createAxis(MV2array(lon[0, :]), id='longitude')
        lonn.units = lon.units
        lon_num = GetNumAxis(tab, 'longitude')
        try:
            averaged_tab.setAxis(lon_num, lonn)
        except:
            averaged_tab.setAxis(lon_num - 1, lonn)
    info = info + '\n' + str().ljust(5) + "meridional average: cdutil.averager"
    return averaged_tab, info


def AverageTemporal(tab, info):
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
        time_num = GetNumAxis(tab, 'time')
        try: averaged_tab = cdutil.averager(tab, axis=str(time_num))
        except:
            list_strings = [
                "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": temporal average",
                str().ljust(5) + "cannot perform temporal average"
            ]
            lib.EnsoErrorsWarnings.MyError(list_strings)
    info = info + '\n' + str().ljust(5) + "temporal average: cdutil.averager"
    return averaged_tab, info


def AverageZonal(tab, info):
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
        lon_num = GetNumAxis(tab, 'longitude')
        try: averaged_tab = cdutil.averager(tab, axis=str(lon_num))
        except:
            list_strings = [
                "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": zonal average",
                str().ljust(5) + "cannot perform zonal average"
            ]
            lib.EnsoErrorsWarnings.MyError(list_strings)
    lat = tab.getLatitude()
    if len(lat.shape) > 1:
        latn = CDMS2createAxis(MV2array(lat[:, 0]), id='latitude')
        latn.units = lat.units
        lat_num = GetNumAxis(tab, 'latitude')
        try:
            averaged_tab.setAxis(lat_num, latn)
        except:
            averaged_tab.setAxis(lat_num - 1, latn)
    info = info + '\n' + str().ljust(5) + "zonal average: cdutil.averager"
    return averaged_tab, info


# Dictionary of averaging methods
dict_average = {'horizontal': AverageHorizontal, 'meridional': AverageMeridional, 'time': AverageTemporal,
                'zonal': AverageZonal}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# functions developed by Y.Y. Planton based on uvcdat
#
def ComputeAnomalies(tab, info, **kwargs):
    # removes annual cycle (anomalies with respect to the annual cycle)
    tab = cdutil.ANNUALCYCLE.departures(tab)
    info = info + '\n' + str().ljust(5)+ "anomalies with respect to the mean seasonal cycle: " + \
           "cdutil.ANNUALCYCLE.departures"
    return tab, info


def ComputeAverage(tab, info, average_dimension=[], **kwargs):
    if isinstance(average_dimension, basestring):
        try:
            dict_average[average_dimension]
        except:
            lib.EnsoErrorsWarnings.UnknownAveraging(average_dimension, dict_average.keys(), INSPECTstack())
        else:
            tab, info = dict_average[average_dimension](tab, info)
    elif isinstance(average_dimension, list):
        for av in average_dimension:
            try:
                dict_average[av]
            except:
                lib.EnsoErrorsWarnings.UnknownAveraging(average_dimension, dict_average.keys(), INSPECTstack())
            else:
                tab, info = dict_average[av](tab, info)
    else:
        lib.EnsoErrorsWarnings.UnknownAveraging(average_dimension, dict_average.keys(), INSPECTstack())
    return tab, info


def ComputeDetrend(tab, info, axis=0, method='linear', bp=0, **kwargs):
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
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": method",
            str().ljust(5) + "unknown method: " + str(method)
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    if method in ['linear', 'constant']:
        axes = tab.getAxisList()
        grid = tab.getGrid()
        mask = tab.mask
        mean = dict_average['time'](tab)
        num_axis = GetNumAxis(tab, axis)
        new_tab = MV2array(SCIPYsignal_detrend(tab, axis=num_axis, type=method, bp=bp))
        new_tab = new_tab + mean
        new_tab = MV2masked_where(mask, new_tab)
        new_tab.setAxisList(axes)
        new_tab.setGrid(grid)
        if method == 'linear':
            info = info + '\n' + str().ljust(5) + "time series linearly detrended: scipy.signal.detrend"
        else:
            info = info + '\n' + str().ljust(5) + "mean value of the time series subtracted: scipy.signal.detrend"
    return new_tab, info


def ComputeNormalize(tab, info, frequency='monthly', **kwargs):
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
        info = info + '\n' + str().ljust(5) + "time series normalized by the daily standard deviation: " + \
               "genutil.statistics.std"
    elif frequency == 'monthly':
        time_steps_per_year = 12
        info = info + '\n' + str().ljust(5) + "time series normalized by the monthly standard deviation: " + \
               "genutil.statistics.std"
    elif frequency == 'yearly':
        time_steps_per_year = 1
        info = info + '\n' + str().ljust(5) + "time series normalized by the yearly standard deviation: " + \
               "genutil.statistics.std"
    else:
        lib.EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())
    if len(tab) % time_steps_per_year != 0:
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": data length",
            str().ljust(5) + "the normalization function needs only full years: " +
            str(len(tab) // time_steps_per_year) + " years + " + str(len(tab) % time_steps_per_year),
            str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " + str(time_steps_per_year) +
            "), len(dataset) = " + str(len(tab)) + ", so " + str(len(tab) / float(time_steps_per_year)) + " years",
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    # reshape tab like [yy,nb]
    new_tab = list()
    for yy in range(len(tab) / time_steps_per_year):
        new_tab.append(tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year])
    new_tab = MV2array(new_tab)
    std = MV2zeros(new_tab[0].shape)
    for dd in range(time_steps_per_year):
        std[dd] = float(GENUTILstd(new_tab[:, dd], weights=None, axis=0, centered=1, biased=1))
    for yy in range(len(tab) / time_steps_per_year):
        tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
            tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] / std
    return tab, info


def ComputeRegrid(tab_to_regrid, newgrid, newgrid_name, info, missing=None, order=None, mask=None, regridTool='esmf',
                  regridMethod='linear', **kwargs):
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
        if newgrid type is string, a grid will be created using **kwargs
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
    # test if regridTool is a good keyword
    if regridTool not in ['regrid2', 'esmf', 'libcf', 'adamsregrid']:
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridTool",
            str().ljust(5) + "unknown regridTool: " + str(regridTool),
            str().ljust(10) + "known regridTool: " + str(['regrid2', 'esmf', 'libcf'])
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    # test if regridTool is a good keyword
    if (regridTool == 'esmf' and regridMethod not in ['conserve', 'patch', 'linear']) or \
            (regridTool in ['regrid2', 'libcf'] and regridMethod != 'linear'):
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridMethod",
            str().ljust(5) + "unknown regridMethod for this regridTool (" + str(regridMethod) + "): "
            + str(regridMethod),
            str().ljust(10) + "known regridTool: " + str(['conserve', 'patch', 'linear'])
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    if (regridTool == 'adamsregrid' and regridMethod not in ['linear', 'linearLog', 'cubic', 'cubicLog']):
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": regridMethod",
            str().ljust(5) + "unknown regridMethod for this regridTool (" + str(regridMethod) + "): "
            + str(regridMethod),
            str().ljust(10) + "known regridTool: " + str(['conserve', 'patch', 'linear'])
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    # test the given 'newgrid'
    if isinstance(newgrid, basestring):
        #
        # newgrid is not a grid, so a grid will be created
        # to do this, kwargs['newgrid_name'] and kwargs['region'] must be defined
        #
        # define the grid type
        for gtype in ['equalarea', 'gaussian', 'generic', 'uniform']:
            if gtype in kwargs['newgrid_name']:
                GridType = gtype
                break
        try: GridType
        except: GridType = 'generic'
        # define resolution (same resolution in lon and lat)
        for res in ['0.25x0.25deg', '0.5x0.5deg', '1x1deg', '2x2deg']:
            if res in kwargs['newgrid_name']:
                if res == '0.25x0.25deg':
                    GridRes = 0.25
                elif res == '0.5x0.5deg':
                    GridRes = 0.5
                elif res == '1x1deg':
                    GridRes = 1.
                else:
                    GridRes = 2.
                break
        try: GridRes
        except: GridRes = '1.'
        # define bounds of 'region'
        region_ref = ReferenceRegions(kwargs['region'])
        lat1, lat2 = region_ref['latitude'][0], region_ref['latitude'][1]
        lon1, lon2 = region_ref['longitude'][0], region_ref['longitude'][1]
        # create uniform axis
        nlat = lat2 - lat1
        lat = CDMS2createUniformLatitudeAxis(lat1+(GridRes/2.), nlat, GridRes)
        nlon = lon2 - lon1
        lon = CDMS2createUniformLongitudeAxis(lon1+(GridRes/2.), nlon, GridRes)
        # create grid
        newgrid = CDMS2createRectGrid(lat, lon, "yx", type=GridType, mask=None)
        newgrid.id = kwargs['newgrid_name']
    #
    # regrid
    #
    if regridTool in ['regrid2', 'esmf', 'libcf']:
        new_tab = tab_to_regrid.regrid(newgrid, missing=missing, order=order, mask=mask, regridTool=regridTool,
                                       regridMethod=regridMethod)
        info = info + '\n' + str().ljust(5) + "regridded toward " + str(newgrid_name) + " using " + str(regridTool) +\
               " " + str(regridMethod) + ": cdms2"
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
        info = info + '\n' + str().ljust(5) + "regridded toward " + str(newgrid_name) + " using " + str(regridTool) + \
               " " + str(regridMethod) + ": adamsregrid"
    return new_tab, info


def ComputeSmooth(tab, info, axis=0, window=5, method='triangle', **kwargs):
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
        smoothed tab
    """
    try: dict_smooth[method]
    except:
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing method (running mean)",
            str().ljust(5) + "unkwown smoothing method: " + str(method),
            str().ljust(10) + "known smoothing method: " + str(sorted(dict_smooth.keys())),
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    info = info + '\n' + str().ljust(5) + "time series smoothized using a " + str(method) + " shaped window of "\
           + str(window) + " points: MV2"
    axis_num = GetNumAxis(tab, axis)
    return dict_smooth[method](tab, axis=axis_num, window=window), info


def CheckUnits(tab, file_name, name_in_file, old_units, new_units):
    """
    #################################################################################
    Description:
    Checks the units of the variable and changes it if necessary
    Works for current/wind velocities, heat flux, precipitation, pressure, temperature, wind stress

    Uses MV2 (uvcdat) to find the minimum value, to multiply and to subtract
    #################################################################################

    :param tab: array
        array containing 'var_name'
    :param name_in_file: string
        name of the variable in the file (usually the short_name)
    :param units: string
        units of the variable included in 'tab'

    :return tab: array
        array with new units (if applicable)
    """
    heat_flux_units = ['W/m^2', 'W/m2', 'W m-2', 'W m**-2']
    precipitations_units = ['kg / (m^2 s^1)', 'kg / (m2 s)', 'kg m-2 s-1', 'kg m**-2 s**-1',
                            'mm/day', 'mm day-1', 'mm day**-1']
    pressure_units = ['N/m^2', 'N/m2', 'N m-2', 'N m**-2', 'Pa']
    temperature_units = ['C', 'degree_Celsius', 'deg_Celsius', 'deg. C', 'degCelsius', 'degree_C', 'deg_C', 'degC',
                         'degrees C', 'K']
    velocity_units = ['cm/s', 'cm s-1', 'cm s**-1', 'cm/sec', 'cm sec-1', 'cm sec**-1',
                      'm/s', 'm s-1', 'm s**-1', 'm/sec', 'm sec-1', 'm sec**-1']
    if new_units in temperature_units and old_units in temperature_units:
        if old_units == 'K':
            # check if the temperature units is really K
            if float(MV2minimum(tab)) > 200:
                # units seems to be K
                if new_units == 'K':
                    tab.units = "K"
                else:
                    # change from K to degC
                    tab = tab - 273.15
                    tab.units = "degC"
            else:
                lib.EnsoErrorsWarnings.UnlikelyUnits(file_name, name_in_file, old_units, INSPECTstack())
        else:
            # check if the temperature units is really degC
            if float(MV2minimum(tab)) < 50:
                # units seems to be degC
                if new_units in ['C', 'degree_Celsius', 'deg_Celsius', 'deg. C', 'degCelsius', 'degree_C', 'deg_C',
                                 'degC', 'degrees C']:
                    tab.units = "degC"
                else:
                    # change from degC to K
                    tab = tab + 273.15
                    tab.units = "K"
            else:
                lib.EnsoErrorsWarnings.UnlikelyUnits(file_name, name_in_file, old_units, INSPECTstack())
    elif new_units in precipitations_units and old_units in precipitations_units:
        if old_units in ['kg / (m^2 s^1)', 'kg / (m2 s)', 'kg m-2 s-1', 'kg m**-2 s**-1']:
            if new_units in ['kg / (m^2 s^1)', 'kg / (m2 s)', 'kg m-2 s-1', 'kg m**-2 s**-1']:
                tab.units = "kg m-2 s-1"
            else:
                # changes from kg/(m2.s) to mm/day
                # it must be divided by the density of water = 1000 kg/m3
                #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
                tab = tab * 86400.
                tab.units = "mm/day"
        else:
            if new_units in ['mm/day', 'mm day-1', 'mm day**-1']:
                tab.units = "mm/day"
            else:
                # changes from mm/day to kg/(m2.s)
                # it must be multiplied by the density of water = 1000 kg/m3
                #     and divided by 1000 (m to mm) and by 60*60*24 (s to day)
                tab = tab / 86400.
                tab.units = "kg m-2 s-1"
    elif new_units in velocity_units and old_units in velocity_units:
        if old_units in ['cm/s', 'cm s-1', 'cm s**-1', 'cm/sec', 'cm sec-1', 'cm sec**-1']:
            if new_units in ['cm/s', 'cm s-1', 'cm s**-1', 'cm/sec', 'cm sec-1', 'cm sec**-1']:
                tab.units = "cm/s"
            else:
                # change from cm/s to m/s
                tab = tab * 1e-2
                tab.units = "mm/s"
        else:
            if new_units in ['m/s', 'm s-1', 'm s**-1', 'm/sec', 'm sec-1', 'm sec**-1']:
                tab.units = "m/s"
            else:
                # change from m/s to cm/s
                tab = tab * 100.
                tab.units = "cm/s"
    elif new_units in heat_flux_units and old_units in heat_flux_units:
        tab.units = "W/m2"
    elif new_units in pressure_units and old_units in pressure_units:
        tab.units = "N/m2"
    else:
        # either input units ('old_units') or output units ('new_units') is unknown
        # or the change of units is impossible / weird / unlikely (like degC to mm/day)
        if old_units not in heat_flux_units and old_units not in precipitations_units \
            and old_units not in pressure_units and old_units not in temperature_units \
                and old_units not in velocity_units:
            list_strings = ["ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": units",
                            str().ljust(5) + "unknown input units: " + str(old_units)]
            lib.EnsoErrorsWarnings.MyError(list_strings)
        elif new_units not in heat_flux_units and new_units not in precipitations_units \
            and new_units not in pressure_units and new_units not in temperature_units \
                and new_units not in velocity_units:
            list_strings = ["ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": units",
                            str().ljust(5) + "unknown output units: " + str(new_units)]
            lib.EnsoErrorsWarnings.MyError(list_strings)
        else:
            list_strings = ["ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": units",
                            str().ljust(5) + "tries to change " + str(old_units) + " to " + str(new_units),
                            str().ljust(5) + "this seems impossible / weird / unlikely"]
            lib.EnsoErrorsWarnings.MyError(list_strings)
    return tab


def GetNumAxis(tab, name_axis, **kwargs):
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
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "cannot find axis named: " + str(name_axis)
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    return num


def ReadAndSelectRegion(filename, varname, info, box=None, time_bounds=None, frequency=None, units=None, **kwargs):
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
            tab = fi(varname, nbox)
#            tab = fi(varname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        else:
            # read file
            tab = fi(varname, nbox, time=time_bounds)
#            tab = fi(varname, time=time_bounds, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
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
        lib.EnsoErrorsWarnings.UnknownFrequency(frequency, INSPECTstack())
    # remove axis 'level' if its length is 1
    if tab.getLevel():
        if len(tab.getLevel()) == 1:
            tab = tab(squeeze=1)
    # set units to the given units
    if units is not None:
        tab = CheckUnits(tab, filename, varname, tab.units, units)
    # info
    info = info + '\n' + str().ljust(5) + "open file: cdms2.setAutoBounds('on') cdms2.open" + \
           '\n' + str().ljust(5) + "select region: cdutil.region.domain"
    return tab, info


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
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = [
            "ERROR" + lib.EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": axis",
            str().ljust(5) + "axis number too big: " + str(axis)
        ]
        lib.EnsoErrorsWarnings.MyError(list_strings)
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
dict_smooth = {'triangle': SmoothTriangle}
# ---------------------------------------------------------------------------------------------------------------------#
