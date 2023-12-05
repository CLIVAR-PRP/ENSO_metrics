# -*- coding:UTF-8 -*-

# basic python
import copy
from inspect import stack as inspect__stack
import ntpath
from os.path import isdir as os__path__isdir
from os.path import isfile as os__path__isfile


# matplotlib
from matplotlib.path import Path as matplotlib__path

# numpy
from numpy import arange as numpy__arange
from numpy import exp as numpy__exp
from numpy import meshgrid as numpy__meshgrid
from numpy import vstack as numpy__vstack
from numpy import where as numpy__where

# scipy
from scipy.signal import detrend as scipy__signal__detrend

# ENSO_metrics package functions:
from . import EnsoErrorsWarnings
from .EnsoCollectionsLib import ReferenceRegions
from .EnsoToolsLib import add_up_errors


# I need to convert
def CDMS2createAxis():
    return
def CDMS2createRectGrid():
    return
def CDMS2createUniformLatitudeAxis():
    return
def CDMS2createUniformLongitudeAxis():
    return
def CDMS2createVariable():
    return
def CDMS2open():
    return
def CDMS2setAutoBounds():
    return
def cdutil():
    #cdutil.ANNUALCYCLE.departures(tab)
    #cdutil.averager()
    #cdutil.generateLandSeaMask
    return
def GENUTILstd():
    return
def MV2add():
    return
def MV2array():
    return
def MV2average():
    return
def MV2compress():
    return
def MV2divide():
    return
def MV2masked_where():
    return
def MV2max():
    return
def MV2min():
    return
def MV2multiply():
    return
def MV2ones():
    return
def MV2subtract():
    return
def MV2sum():
    return
def MV2take():
    return
def MV2where():
    return
def MV2zeros():
    return
def REGRID2horizontal__Horizontal():
    return

# math
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
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, inspect__stack())
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
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, inspect__stack())
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
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, inspect__stack())
    tab_out = MV2multiply(tab, number_or_tab)
    axes = tab.getAxisList()
    att = tab.attributes
    if len(tab.shape) > 1:
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
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, inspect__stack())
    return MV2subtract(tab, number_or_tab)


dict_operations = {"divide": operation_divide, "minus": operation_subtract, "multiply": operation_multiply,
                   "plus": operation_add}


def compute_annual_cycle(tab):
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
    tab = tab.reorder("t...")
    axes = tab.getAxisList()
    time_ax = tab.getTime().asComponentTime()
    months = MV2array(list(tt.month for tt in time_ax))
    cyc = []
    for ii in range(12):
        ids = MV2compress(months == (ii + 1), list(range(len(tab))))
        tmp = MV2take(tab, ids, axis=0)
        # tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = MV2average(tmp, axis=0)
        cyc.append(tmp)
        del tmp
    time = CDMS2createAxis(list(range(12)), id="time")
    moy = CDMS2createVariable(MV2array(cyc), axes=[time] + axes[1:], grid=tab.getGrid(), attributes=tab.attributes)
    moy = moy.reorder(initorder)
    time = CDMS2createAxis(list(range(12)), id="months")
    moy.setAxis(get_num_axis(moy, "time"), time)
    return moy


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


def apply_landmask(tab, landmask, maskland=True, maskocean=False):
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
            keyerror = "tab (" + str(tab.getGrid().shape) + ") and landmask (" + str(landmask.getGrid().shape) + \
                       ") are not on the same grid"
            list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": applying landmask",
                            str().ljust(5) + keyerror, str().ljust(5) + "cannot apply landmask",
                            str().ljust(5) + "this metric will be skipped"]
            EnsoErrorsWarnings.my_warning(list_strings)
        else:
            landmask_nd = MV2zeros(tab.shape)
            if landmask_nd.shape == landmask.shape:
                landmask_nd = copy.copy(landmask)
            else:
                try:
                    landmask_nd[:] = landmask
                except:
                    try:
                        landmask_nd[:, :] = landmask
                    except:
                        keyerror = "ApplyLandmask: tab must be more than 4D and this is not taken into account yet (" + \
                                   str(tab.shape) + ") and landmask (" + str(landmask.shape) + ")"
                        list_strings = [
                            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": landmask shape",
                            str().ljust(5) + keyerror, str().ljust(5) + "cannot reshape landmask"]
                        EnsoErrorsWarnings.my_warning(list_strings)
            if keyerror is None:
                # if land = 100 instead of 1, divides landmask by 100
                if (MV2min(landmask) == 0 and MV2max(landmask) == 100) or \
                        ("units" in landmask.listattributes() and landmask.units == "%"):
                    landmask_nd = landmask_nd / 100.
                if maskland is True:
                    tab = MV2masked_where(landmask_nd >= 0.2, tab)
                if maskocean is True:
                    tab = MV2masked_where(landmask_nd <= 0.8, tab)
    return tab, keyerror


def apply_landmask_to_area(area, landmask, maskland=True, maskocean=False):
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
            keyerror = "ApplyLandmaskToArea: area (" + str(area.getGrid().shape) + ") and landmask (" + \
                       str(landmask.getGrid().shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": applying landmask to areacell",
                str().ljust(5) + keyerror, str().ljust(5) + "cannot apply landmask to areacell"]
            EnsoErrorsWarnings.my_warning(list_strings)
        if keyerror is None:
            # if land = 100 instead of 1, divides landmask by 100
            if (MV2min(landmask) == 0 and MV2max(landmask) == 100) or \
                    ("units" in landmask.listattributes() and landmask.units == "%"):
                landmask = landmask / 100.
            if maskland is True:
                area = MV2masked_where(landmask >= 0.2, area)
                area = MV2multiply(area, 1 - landmask)
            if maskocean is True:
                area = MV2masked_where(landmask <= 0.8, area)
                area = MV2multiply(area, landmask)
    return area, keyerror


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
    keyerror = None
    lat_num = get_num_axis(tab, "latitude")
    lon_num = get_num_axis(tab, "longitude")
    tmp = sorted([lat_num, lon_num])
    snum = str(tmp[0]) + str(tmp[1])
    if areacell is None or tab.getGrid().shape != areacell.getGrid().shape:
        print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib average_horizontal" + "\033[0m")
        if areacell is not None and tab.getGrid().shape != areacell.getGrid().shape:
            print("\033[93m" + str().ljust(25) + "tab.grid " + str(tab.getGrid().shape) +
                  " is not the same as areacell.grid " + str(areacell.getGrid().shape) + " \033[0m")
        try:
            averaged_tab = cdutil.averager(tab, axis="xy", weights="weighted", action="average")
        except:
            try:
                averaged_tab = cdutil.averager(tab, axis=snum, weights="weighted", action="average")
            except:
                if "regridding" not in list(kwargs.keys()) or isinstance(kwargs["regridding"], dict) is False:
                    kwargs2 = {"regridder": "cdms", "regridTool": "esmf", "regridMethod": "linear",
                               "newgrid_name": "generic_1x1deg"}
                else:
                    kwargs2 = kwargs["regridding"]
                kwargs2["newgrid_name"] = \
                    find_closest_grid(region, len(tab.getAxis(lat_num)[:]), len(tab.getAxis(lon_num)[:]))
                print("\033[93m" + str().ljust(25) + "need to regrid to = " + str(kwargs2["newgrid_name"]) +
                      " to perform average \033[0m")
                tmp = regridder(tab, None, region=region, **kwargs2)
                try:
                    averaged_tab = cdutil.averager(tmp, axis=snum, weights="weighted", action="average")
                except:
                    keyerror = "cannot perform horizontal average"
                    averaged_tab = None
                    list_strings = [
                        "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": horizontal average",
                        str().ljust(5) + "cdutil.averager cannot perform horizontal average"]
                    EnsoErrorsWarnings.my_warning(list_strings)
    else:
        averaged_tab = MV2multiply(tab, areacell)
        for elt in snum[::-1]:
            averaged_tab = MV2sum(averaged_tab, axis=int(elt))
        averaged_tab = averaged_tab / float(MV2sum(areacell))
    return averaged_tab, keyerror


def averager_polygon(larray, larea, lregion, input_region="global2", **kwargs):
    """
    #################################################################################
    Description:
    Adapt the given array to keep only the given region, even if the given region is a polygon
    #################################################################################

    :param larray: masked_array
        masked_array to average
    :param larea: masked_array
        masked_array with the area of each grid cell
    :param lregion: string
        name of a region (e.g., 'equatorial_pacific') to select and then average
        must be defined in .EnsoCollectionsLib.ReferenceRegions
    :param input_region: string, optional
        name of a region (e.g., 'global2') of the given larray (in case it needs to be regridded)
        must be defined in .EnsoCollectionsLib.ReferenceRegions
        default value is 'global2'
    :param kwargs: dict, optional
        dict with processing information, such as regridding if needed
        e.g.,
        kwargs = {
            "regridding": {"regridder": "cdms", "regridTool": "esmf", "regridMethod": "linear",
                           "newgrid_name": "generic_1x1deg"},
        }

    :return: larray_out: masked_array
        larray horizontally averaged in given lregion if the computation is successful, else, None
    :return: keyerror: string or None
        None if the computation is successful, else, a string defining the error
    """
    keyerror = None
    # get the definition of the region
    region_ref = ReferenceRegions(lregion)
    # if the region is a polygon, set unwanted area to 0 in unwanted region
    if "polygon" in list(region_ref.keys()) and region_ref["polygon"] is True:
        # if the grid is not regular, regrid
        if len(larray.getLatitude()[:].shape) > 1:
            # check if regridding is defined, if not, set a default
            if "regridding" not in list(kwargs.keys()) or isinstance(kwargs["regridding"], dict) is False:
                kwargs2 = {"regridder": "cdms", "regridTool": "esmf", "regridMethod": "linear",
                           "newgrid_name": "generic_1x1deg"}
            else:
                kwargs2 = kwargs["regridding"]
            # find generic horizontal resolution closest to the original grid
            lat_num = get_num_axis(larray, "latitude")
            lon_num = get_num_axis(larray, "longitude")
            kwargs2["newgrid_name"] = \
                find_closest_grid(input_region, len(larray.getAxis(lat_num)[:]), len(larray.getAxis(lon_num)[:]))
            print("\033[93m" + str().ljust(25) + "need to regrid to = " + str(kwargs2["newgrid_name"]) +
                  " to perform average \033[0m")
            # regrid given array
            larray_in = regridder(larray, None, region=input_region, **kwargs2)
            # if the given area does not correspond to the given array, create an area of ones, else, regrid given area
            if larea is None or larray.getGrid().shape != larea.getGrid().shape:
                larea_in = create_array_ones(larray_in[0], mid="areacell")
            else:
                larea_in = regridder(larea, None, region=lregion, **kwargs2)
            del kwargs2, lat_num, lon_num
        else:
            # no need to modify given array
            larray_in = copy.copy(larray)
            # if the given area does not correspond to the given array, create an area of ones, else, regrid given area
            if larea is None or larray.getGrid().shape != larea.getGrid().shape:
                larea_in = create_array_ones(larray_in[0], mid="areacell")
            else:
                larea_in = copy.copy(larea)
        # get lat and lon
        lat = copy.deepcopy(larea_in.getLatitude()[:])
        lon = copy.deepcopy(larea_in.getLongitude()[:])
        # map size
        nx, ny = len(lon), len(lat)
        # loop on the polygon's vertices
        poly_verts = list()
        for l1, l2 in zip(region_ref["latitude"], region_ref["longitude"]):
            # find the closest lat and lon grid point and get its coordinates (position in the grid)
            i1 = min(range(len(lat)), key=lambda ii: abs(lat[ii] - l1))
            i2 = min(range(len(lon)), key=lambda ii: abs(lon[ii] - l2))
            # define found coordinates as the vertex of the polygon
            poly_verts.append((i2, i1))
            del i1, i2
        # Create vertex coordinates for each grid cell...
        # (<0,0> is at the top left of the grid in this system)
        x, y = numpy__meshgrid(numpy__arange(nx), numpy__arange(ny))
        x, y = x.flatten(), y.flatten()
        points = numpy__vstack((x, y)).T
        path = matplotlib__path(poly_verts)
        mask1 = path.contains_points(points)
        mask1 = mask1.reshape((ny, nx))
        if min(region_ref["longitude"]) < 0 or max(region_ref["longitude"]) > 360:
            # loop on the polygon's vertices
            poly_verts = list()
            for l1, l2 in zip(region_ref["latitude"], region_ref["longitude"]):
                # change longitude
                if min(region_ref["longitude"]) < 0:
                    l3 = l2 + 360
                else:
                    l3 = l2 - 360
                # find the closest lat and lon grid point and get its coordinates (position in the grid)
                i1 = min(range(len(lat)), key=lambda ii: abs(lat[ii] - l1))
                i2 = min(range(len(lon)), key=lambda ii: abs(lon[ii] - l3))
                # define found coordinates as the vertex of the polygon
                poly_verts.append((i2, i1))
                del i1, i2, l3
            # Create vertex coordinates for each grid cell...
            x, y = numpy__meshgrid(numpy__arange(nx), numpy__arange(ny))
            x, y = x.flatten(), y.flatten()
            points = numpy__vstack((x, y)).T
            path = matplotlib__path(poly_verts)
            mask2 = path.contains_points(points)
            mask2 = mask2.reshape((ny, nx))
            # add mask2 to mask1
            mask1 = MV2where(mask2, mask2, mask1)
            del mask2
        # this new area is 0 outside the polygon, cell area inside the polygon
        larea_new = MV2where(mask1, larea_in, 0)
        larea_new = CDMS2createVariable(larea_new, axes=larea_in.getAxisList(), grid=larea_in.getGrid(),
                                        mask=larea_in.mask, id="areacell")
        # create a mask of the same dimension as larray_in
        mask_nd = MV2ones(larray_in.shape)
        if mask_nd.shape == larea_new.shape:
            mask_nd = MV2where(larea_new == 0, 0, mask_nd)
        elif mask_nd[0].shape == larea_new.shape:
            mask_nd[:] = MV2where(larea_new == 0, 0, mask_nd[0])
        elif mask_nd[0, 0].shape == larea_new.shape:
            mask_nd[:, :] = MV2where(larea_new == 0, 0, mask_nd[0, 0])
        elif mask_nd[0, 0, 0].shape == larea_new.shape:
            mask_nd[:, :, :] = MV2where(larea_new == 0, 0, mask_nd[0, 0, 0])
        elif mask_nd[0, 0, 0, 0].shape == larea_new.shape:
            mask_nd[:, :, :, :] = MV2where(larea_new == 0, 0, mask_nd[0, 0, 0, 0])
        else:
            keyerror = "averager_polygon, create mask: given array must be more than 4D and this is not taken into " + \
                       "account yet (" + str(larray_in.shape) + ") and area (" + str(larea_new.shape) + ")"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": region mask shape",
                str().ljust(5) + keyerror, str().ljust(5) + "cannot reshape region mask"]
            EnsoErrorsWarnings.my_warning(list_strings)
        # set larray to 0 outside the polygon
        if keyerror is None:
            laxes, lgrid, lmask = larray_in.getAxisList(), larray_in.getGrid(), larray_in.mask
            larray_in = MV2where(mask_nd == 0, 0, larray_in)
            larray_in = MV2masked_where(lmask, larray_in)
            larray_in.setAxisList(laxes)
            larray_in.setGrid(lgrid)
            del laxes, lgrid, lmask
        del larea_in, mask_nd, mask1, nx, ny, path, points, poly_verts, x, y
    else:
        # select region in given array
        larray_in = larray(latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        # select region in given areacell
        if larea is not None:
            larea_new = larea(latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        else:
            larea_new = copy.copy(larea)
    # use usual averager with modified array and/or area if applicable
    if keyerror is None:
        larray_out, keyerror = average_horizontal(larray_in, areacell=larea_new, region=lregion, **kwargs)
    else:
        larray_out = None
    return larray_out, keyerror


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
    keyerror = None
    try:
        averaged_tab = cdutil.averager(tab, axis="t")
    except:
        time_num = get_num_axis(tab, "time")
        try:
            averaged_tab = cdutil.averager(tab, axis=str(time_num))
        except:
            keyerror = "cannot perform temporal average"
            averaged_tab = None
            list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": temporal average",
                            str().ljust(5) + "cannot perform temporal average"]
            EnsoErrorsWarnings.my_warning(list_strings)
    return averaged_tab, keyerror


def check_units(tab, var_name, name_in_file, units, return_tab_only=True, **kwargs):
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
    if var_name in ["temperature"]:
        if units in ["K", "Kelvin", "Kelvins", "degree K", "degree Kelvin", "degree Kelvins", "degree_K",
                     "degree_Kelvin", "degree_Kelvins", "degreeK", "degreeKelvin", "degreeKelvins", "degrees K",
                     "degrees Kelvin", "degrees Kelvins", "degrees_K", "degrees_Kelvin", "degrees_Kelvins", "degreesK",
                     "degreesKelvin", "degreesKelvins", "deg K", "deg Kelvin", "deg Kelvins", "deg_K", "deg_Kelvin",
                     "deg_Kelvins", "degK", "degKelvin", "degKelvins", "deg. K", "deg. Kelvin", "deg. Kelvins"]:
            # check if the temperature units is really K
            if float(MV2min(tab)) > 150:
                # unit change of the temperature: from K to degC
                tab = dict_operations["minus"](tab, 273.15)
            else:
                minmax = get_min_and_max(tab)
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, units, minmax, inspect__stack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        elif units in ["C", "celsius", "Celsius", "degree C", "degree celsius", "degree Celsius", "degree_C",
                       "degree_celsius", "degree_Celsius", "degreeC", "degreecelsius", "degreeCelsius", "degrees C",
                       "degrees celsius", "degrees Celsius", "degrees_C", "degrees_celsius", "degrees_Celsius",
                       "degreesC", "degreescelsius", "degreesCelsius", "deg C", "deg celsius", "deg Celsius", "deg_C",
                       "deg_celsius", "deg_Celsius", "degC", "degcelsius", "degCelsius", "deg. C", "deg. celsius",
                       "deg. Celsius"]:
            # check if the temperature units is really degC
            if float(MV2min(tab)) > 50:
                minmax = get_min_and_max(tab)
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, minmax, units, inspect__stack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "degC"
    elif var_name in ["precipitation"]:
        if units in ["kg/m2/s", "kg/m^2/s", "kg/m**2/s", "kg m-2 s-1", "kg m^-2 s^-1", "kg m**-2 s**-1", "Kg/m2/s",
                     "Kg/m^2/s", "Kg/m**2/s", "Kg m-2 s-1", "Kg m^-2 s^-1", "Kg m**-2 s**-1"]:
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = dict_operations["multiply"](tab, 86400)
        elif units in ["mm/day", "mm day-1", "mm day^-1", "mm day**-1", "mm/d", "mm d-1", "mm d^-1", "mm d**-1"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "mm/day"
    elif var_name in ["wind stress"]:
        if units not in ["N/m2", "N/m^2", "N/m**2", "N m-2", "N m^-2", "N m**-2", "Pa", "pascal", "pascals", "Pascal",
                         "Pascals"]:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        else:
            tab = dict_operations["multiply"](tab, 1e3)
        units = "10-3 N/m2"
    elif var_name in ["velocity"]:
        if units in ["cm/s", "cm s-1", "cm s^-1", "cm s**-1", "cm/sec", "cm sec-1", "cm sec^-1", "cm sec**-1"]:
            # unit change of the velocity: from cm/s to m/s
            tab = dict_operations["multiply"](tab, 1e-2)
        elif units in ["m/s", "m s-1", "m s^-1", "m s**-1", "m/sec", "m sec-1", "m sec^-1", "m sec**-1"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m/s"
    elif var_name in ["heat flux"]:
        if units in ["W/m2", "W/m^2", "W/m**2", "W m-2", "W m^-2", "W m**-2", "Watt/m2", "Watt/m^2", "Watt/m**2",
                     "Watt m-2", "Watt m^-2", "Watt m**-2", "Watts/m2", "Watts/m^2", "Watts/m**2", "Watts m-2",
                     "Watts m^-2", "Watts m**-2"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "W/m2"
    elif var_name in ["pressure"]:
        if units in ["N/m2", "N/m^2", "N/m**2", "N m-2", "N m^-2", "N m**-2", "Pa", "pascal", "pascals", "Pascal",
                     "Pascals"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "Pa"
    elif var_name in ["depth"]:
        if units in ["mm", "millimetre", "millimetres"]:
            # unit change of the depth: from mm to m
            tab = dict_operations["multiply"](tab, 1e-3)
        elif units in ["cm", "centimeter", "centimeters"]:
            # unit change of the depth: from cm to m
            tab = dict_operations["multiply"](tab, 1e-2)
        elif units in ["m", "meter", "meters"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m"
    elif var_name in ["sea surface height"]:
        if units in ["mm", "millimetre", "millimetres"]:
            # unit change of the depth: from mm to cm
            tab = dict_operations["multiply"](tab, 1e-1)
        elif units in ["cm", "centimeter", "centimeters"]:
            pass
        elif units in ["m", "meter", "meters"]:
            # unit change of the sea surface height: from m to cm
            tab = dict_operations["multiply"](tab, 1e2)
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, inspect__stack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "cm"
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.my_warning(list_strings)
    if return_tab_only is True:
        return tab
    else:
        return tab, units, keyerror


def create_array_ones(tab, mid="new_variable_ones"):
    """
    #################################################################################
    Description:
    Create a masked_array filled with ones with the same properties as tab (shape, axes, grid, mask)
    #################################################################################

    for more information:
    import MV2
    help(MV2.ones)
    """
    return CDMS2createVariable(MV2ones(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask, id=mid)


def detrendder(tab, info, axis=0, method="linear", bp=0):
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
        "constant": only the mean of 'tab' is subtracted
        "linear":   the result of a linear least-squares fit to 'tab' is subtracted from 'tab'
    :param bp: array of integer, optional
        a sequence of break points. If given, an individual linear fit is performed for each part of 'tab' between two
        break points
        break points are specified as indices into 'tab'
    :return new_tab: array
        detrended data
    """
    if method not in ["linear", "constant"]:
        new_tab = None
        keyerror = "cannot detrend: unknown method"
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": method",
            str().ljust(5) + "unknown method: " + str(method)
        ]
        EnsoErrorsWarnings.my_warning(list_strings)
    else:
        axes = tab.getAxisList()
        grid = tab.getGrid()
        mask = tab.mask
        mean, keyerror = average_temporal(tab)
        new_tab = MV2array(scipy__signal__detrend(tab, axis=axis, type=method, bp=bp))
        new_tab = new_tab + mean
        new_tab = MV2masked_where(mask, new_tab)
        new_tab.setAxisList(axes)
        new_tab.setGrid(grid)
        if method == "linear":
            info = info + ", time series are linearly detrended"
        else:
            info = info + ", the mean value of the time series is subtracted"
    return new_tab, info, keyerror


def estimate_landmask(d):
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
    n = 1
    lmsk = cdutil.generateLandSeaMask(d(*(slice(0, 1),) * n)) * 100.0
    lmsk[:] = lmsk.filled(100.0)
    lmsk.setAxis(0, d.getAxis(1))
    lmsk.setAxis(1, d.getAxis(2))
    lmsk.id = "sftlf"
    return lmsk


def get_min_and_max(tab):
    return [float(MV2min(tab)), float(MV2max(tab))]


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
        e.g., name_axis="latitude"
    :return number: int
        position of the axis named "name_axis"
    """
    num = None
    if name_axis == "depth":
        axis_nick = "lev"
        axis_nicks = ["Level", "st_ocean", "sw_ocean", "z", "Z"]
    elif name_axis == "latitude":
        axis_nick = "lat"
        axis_nicks = ["j", "Latitude", "y", "Y", "yt_ocean", "yu_ocean"]
    elif name_axis == "longitude":
        axis_nick = "lon"
        axis_nicks = ["i", "Longitude", "x", "X", "xt_ocean", "xu_ocean"]
    elif name_axis == "time":
        axis_nick = "time"
        axis_nicks = ["t", "T", "Time"]
    else:
        axis_nick = None
        axis_nicks = [None]
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": axis",
                        str().ljust(5) + "unknown axis named: " + str(name_axis),
                        str().ljust(5) + "known names: depth, latitude, longitude, time"]
        EnsoErrorsWarnings.my_error(list_strings)
    for nn in range(len(tab.shape)):
        if axis_nick in str(tab.getAxisList()[nn].id).lower():
            num = nn
            break
    if num is None:
        for nn in range(len(tab.shape)):
            for ax in axis_nicks:
                if ax == str(tab.getAxisList()[nn].id).lower():
                    num = nn
                    break
    if num is None:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": axis",
                        str().ljust(5) + "cannot find axis named: " + str(name_axis),
                        str().ljust(5) + "axes: " + str(tab.getAxisList())]
        EnsoErrorsWarnings.my_error(list_strings)
    return num


def get_time_bounds(tab):
    """
    #################################################################################
    Description:
    Finds first and last dates of tab's time axis, tab must be a uvcdat masked_array
    #################################################################################

    Returns a tuple of strings: e.g., ('1979-1-1 11:59:60.0', '2016-12-31 11:59:60.0')
    """
    time = tab.getTime().asComponentTime()
    return str(time[0]), str(time[-1])


def find_closest_grid(region, nlat, nlon):
    res = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75]
    region_ref = ReferenceRegions(region)
    lats = region_ref["latitude"]
    dy = float(abs(max(lats) - min(lats))) / nlat
    lyy = [abs(dy - ii) for ii in res]
    lyy = res[lyy.index(min(lyy))]
    lons = region_ref["longitude"]
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


def normalizer(tab, frequency):
    """
    #################################################################################
    Description:
    Removes trend along 'axis' from 'tab'
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency="monthly"
    :return tab: masked_array
        normalized data
    """
    keyerror = None
    tab_out = None
    axes = tab.getAxisList()
    if frequency == "daily":
        time_steps_per_year = 365
    elif frequency == "monthly":
        time_steps_per_year = 12
    elif frequency == "yearly":
        time_steps_per_year = 1
    else:
        keyerror = "unknown frequency"
        time_steps_per_year = None
        EnsoErrorsWarnings.unknown_frequency(frequency, inspect__stack())
    if time_steps_per_year is not None:
        if len(tab) % time_steps_per_year != 0:
            tab_out = None
            keyerror = "cannot perform normalization: the function can only handle full years (len(tab) = " + \
                       str(len(tab)) + ")"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": data length",
                str().ljust(5) + "the normalization function can only handle full years: " +
                str(len(tab) // time_steps_per_year) + " years " + str(len(tab) % time_steps_per_year),
                str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " +
                str(time_steps_per_year) + "), len(dataset) = " + str(len(tab)) + ", so " +
                str(len(tab) / float(time_steps_per_year)) + " years",
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
    else:
        # reshape tab like [yy,nb]
        new_tab = list()
        for yy in range(len(tab) // time_steps_per_year):
            new_tab.append(tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year])
        new_tab = MV2array(new_tab)
        std = MV2zeros(new_tab[0].shape)
        for dd in range(time_steps_per_year):
            std[dd] = float(GENUTILstd(new_tab[:, dd], weights=None, axis=0, centered=1, biased=1))
        tab_out = copy.copy(tab)
        for yy in range(len(tab) // time_steps_per_year):
            tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
                tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] / std
        if len(tab.shape) == 1:
            tab_out = CDMS2createVariable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
        else:
            axes = axes + tab.getAxisList()[1:]
            grid = tab.getGrid()
            mask = tab.mask
            tab_out = CDMS2createVariable(tab_out, axes=axes, grid=grid, mask=mask, attributes=tab.attributes,
                                          id=tab.id)
    return tab_out, keyerror


def preprocess_ts_polygon(tab, info, areacell=None, average="horizontal", compute_anom=False, compute_sea_cycle=False,
                          region=None, **kwargs):
    keyerror = None
    # average
    if average == "horizontal" and isinstance(region, str) is True:
        tab, keyerror = averager_polygon(tab, areacell, region, **kwargs)
        if keyerror is None:
            region_ref = ReferenceRegions(region)
            llat = region_ref["latitude"]
            llat0 = int(round(llat[0], 0)) if int(llat[0]) == llat[0] else copy.deepcopy(llat[0])
            llat1 = int(round(llat[1], 0)) if int(llat[1]) == llat[1] else copy.deepcopy(llat[1])
            llat0 = str(abs(llat0)) + "S" if llat[0] < 0 else str(llat0) + "N"
            llat1 = str(abs(llat1)) + "S" if llat[1] < 0 else str(llat1) + "N"
            llon = region_ref["longitude"]
            llon0 = int(round(llon[0], 0)) if int(llon[0]) == llon[0] else copy.deepcopy(llon[0])
            llon1 = int(round(llon[1], 0)) if int(llon[1]) == llon[1] else copy.deepcopy(llon[1])
            llon0 = str(llon0) if llon[0] == 180 else (
                str(llon0) + "E" if llon[0] < 180 else str(abs(llon0 - 360)) + "W")
            if "E" in llon0 and int(llon0.replace("E", "")) < 0:
                llon0 = str(abs(int(llon0.replace("E", "")))) + "W"
            llon1 = str(llon1) if llon[1] == 180 else (
                str(llon1) + "E" if llon[1] < 180 else str(abs(llon1 - 360)) + "W")
            if region in ["global", "global_no_poles", "tropic"]:
                llon0, llon1 = "0E", "360E"
            info = str(info) + ", " + str(average) + " average in " + str(region) + " region [" + str(llat0) + " - "
            info += str(llat1) + "; " + str(llon0) + " - " + str(llon1) + "]"
    # continue preprocessing if no error happened yet
    if keyerror is None:
        # removing linear trend
        if isinstance(kwargs["detrending"], dict) is True:
            known_args = {"axis", "method", "bp"}
            extra_args = set(kwargs["detrending"]) - known_args
            if extra_args:
                EnsoErrorsWarnings.unknown_key_arg(extra_args, inspect__stack())
            tab, info, keyerror = detrendder(tab, info, **kwargs["detrending"])
        if keyerror is None:
            # computes mean annual cycle
            if compute_sea_cycle is True:
                tab = compute_annual_cycle(tab)
                info = info + ", seasonal cycle computed"
            # removes annual cycle (anomalies with respect to the annual cycle)
            if compute_anom is True:
                tab = compute_interannual_anomalies(tab)
                info = info + ", interannual anomalies computed"
            # normalization of the anomalies
            if kwargs["normalization"] is True:
                if kwargs["frequency"] is not None:
                    tab, keyerror = normalizer(tab, kwargs["frequency"])
                    info = info + ", normalized"
            # continue preprocessing if not error happened yet
            if keyerror is None:
                # smoothing time series
                if isinstance(kwargs["smoothing"], dict) is True:
                    known_args = {"axis", "method", "window"}
                    extra_args = set(kwargs["smoothing"]) - known_args
                    if extra_args:
                        EnsoErrorsWarnings.unknown_key_arg(extra_args, inspect__stack())
                    tab, info = smoother(tab, info, **kwargs["smoothing"])
            else:
                tab = None
        else:
            tab = None
    else:
        tab = None
    return tab, info, keyerror


def read_and_select_region(filename, varname, box=None, time_bounds=None, frequency=None, **kwargs):
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
        e.g., frequency="monthly"
        default value is None

    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds("on")
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
            tab = fi(varname, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        else:
            # read file
            tab = fi(varname, time=time_bounds, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
    # mask large values
    tab = MV2masked_where(tab >= 1e20, tab)
    tab = MV2masked_where(tab <= -1e20, tab)
    # sign correction
    try:
        att1 = tab.attributes["standard_name"].lower().replace(" ", "_")
    except:
        att1 = ""
    try:
        att2 = tab.attributes["long_name"].lower().replace(" ", "_")
    except:
        att2 = ""
    reversed_sign = False
    if "latent_heat" in att1 or "latent_heat" in att2 or "sensible_heat" in att1 or "sensible_heat" in att2 or \
            (varname in ["tauu", "tauuo", "tauv", "tauvo", "taux", "tauy", "uflx", "vflx"]):
        if "upward" in att1 or "upward" in att2 or \
                (varname in ["tauu", "tauuo", "tauv", "tauvo", "taux", "tauy", "uflx", "vflx"] and
                 ("in_air" in att1 or "in_air" in att2)):
            # I need to be in the ocean point of view so the heat fluxes must be downwards
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib ReadAndSelectRegion" + "\033[0m")
            print("\033[93m" + str().ljust(25) + varname + " sign reversed" + "\033[0m")
            print("\033[93m" + str().ljust(5) + "range old = " + "{0:+.2f}".format(round(float(MV2min(tab)), 2)) +
                  " to " + "{0:+.2f}".format(round(float(MV2max(tab)), 2)) + "\033[0m")
            tab = -1 * tab
            print("\033[93m" + str().ljust(5) + "range new = " + "{0:+.2f}".format(round(float(MV2min(tab)), 2)) +
                  " to " + "{0:+.2f}".format(round(float(MV2max(tab)), 2)) + "\033[0m")
            reversed_sign = True
    if time_bounds is not None:
        # sometimes the time boundaries are wrong, even with 'time=time_bounds'
        # this section checks if one time step has not been included by error at the beginning or the end of the time
        # series
        if isinstance(time_bounds[0], str):
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
    elif frequency == "daily":
        cdutil.setTimeBoundsDaily(tab)
    elif frequency == "monthly":
        cdutil.setTimeBoundsMonthly(tab)
    elif frequency == "yearly":
        cdutil.setTimeBoundsYearly(tab)
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, inspect__stack())
    # remove axis 'level' if its length is 1
    if tab.getLevel():
        if len(tab.getLevel()) == 1:
            tab = tab(squeeze=1)
    # HadISST has -1000 values... mask them
    if "HadISST" in filename or "hadisst" in filename:
        tab = MV2masked_where(tab == -1000, tab)
    # check if the mask is constant through time
    if len(numpy__where(tab[0].mask != tab[1:].mask)[0]) > 0:
        # the mask is not constant -> make it constant
        # sum mask through time
        mask = MV2sum(tab.mask.astype("f"), axis=0)
        # mask where at least one time step is masked
        mask = MV2where(mask > 0, True, False)
        # create a mask the same size as the original data
        mask_nd = MV2zeros(tab.shape)
        mask_nd[:] = mask
        # apply mask to original data
        tab = MV2masked_where(mask_nd, tab)
    # check wind stress sign
    if varname in ["taux", "tauu", "tauuo", "uflx", "tauy", "tauv", "tauvo", "vflx"] and reversed_sign is False:
        lreg = "nino4" if varname in ["taux", "tauu", "tauuo", "uflx"] else "itcz_se"
        # define box
        lregion_ref = ReferenceRegions(lreg)
        if time_bounds is None:  # no time period given
            #  read file
            ltau = fi(varname, latitude=lregion_ref["latitude"], longitude=lregion_ref["longitude"])
        else:
            # read file
            ltau = fi(varname, time=time_bounds, latitude=lregion_ref["latitude"], longitude=lregion_ref["longitude"])
        # horizontal average
        ltau, lkeyerror = average_horizontal(ltau, region=lreg)
        if lkeyerror is None:
            ltau, lkeyerror = average_temporal(ltau)
            if lkeyerror is None and ((varname in ["taux", "tauu", "tauuo", "uflx"] and float(ltau) > 0) or
                                      (varname in ["tauy", "tauv", "tauvo", "vflx"] and float(ltau) < 0)):
                print("\033[93m" + str().ljust(25) + "NOTE: wind stress sign reversed by the code (mean " + str(lreg) +
                      " = " + str(float(ltau)) + ")" + "\033[0m")
                tab = -1 * tab
        # delete
        del lkeyerror, lreg, lregion_ref, ltau
    fi.close()
    return tab


def read_area_select_region(filename, areaname="", box=None, **kwargs):
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
    CDMS2setAutoBounds("on")
    # Open file and get time dimension
    fi = CDMS2open(filename)
    if box is None:  # no box given
        # read file
        try:
            areacell = fi(areaname)
        except:
            try:
                areacell = fi("areacell")
            except:
                try:
                    areacell = fi("areacella")
                except:
                    try:
                        areacell = fi("areacello")
                    except:
                        areacell = None
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        # read file
        try:
            areacell = fi(areaname, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        except:
            try:
                areacell = fi('areacell', latitude=region_ref["latitude"], longitude=region_ref["longitude"])
            except:
                try:
                    areacell = fi('areacella', latitude=region_ref["latitude"], longitude=region_ref["longitude"])
                except:
                    try:
                        areacell = fi('areacello', latitude=region_ref["latitude"], longitude=region_ref["longitude"])
                    except:
                        areacell = None
    fi.close()
    return areacell


def read_select_region_check_units(filename, varname, varfamily, box=None, time_bounds=None, frequency=None, **keyarg):
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
        e.g., frequency="monthly"
        default value is None

    :return tab: masked_array
        masked_array containing 'varname' in 'box'
    """
    tab = read_and_select_region(filename, varname, box=box, time_bounds=time_bounds, frequency=frequency)
    tab, units, keyerror = check_units(tab, varfamily, varname, tab.units, return_tab_only=False)
    tab.name = varname
    tab.units = units
    return tab, keyerror


def read_data_mask_area(file_data, name_data, type_data, metric, region, file_area="", name_area="", file_mask="",
                        name_mask="", maskland=False, maskocean=False, time_bounds=None, debug=False, **kwargs):
    keyerror1, keyerror2, keyerror3 = None, None, None
    # Read variable
    if debug is True:
        dict_debug = {"file1": "(" + type_data + ") " + str(file_data), "var1": "(" + type_data + ") " + str(name_data)}
        EnsoErrorsWarnings.debug_mode("\033[93m", "Files", 20, **dict_debug)
    variable, keyerror1 = read_select_region_check_units(file_data, name_data, type_data, box=region,
                                                         time_bounds=time_bounds, **kwargs)
    if debug is True:
        dict_debug = {"axes1": "(" + type_data + ") " + str([ax.id for ax in variable.getAxisList()]),
                      "shape1": "(" + type_data + ") " + str(variable.shape),
                      "time1": "(" + type_data + ") " + str(get_time_bounds(variable))}
        EnsoErrorsWarnings.debug_mode("\033[93m", "after ReadSelectRegionCheckUnits", 20, **dict_debug)
    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs["min_time_steps"], int):
        if len(variable) < kwargs["min_time_steps"]:
            EnsoErrorsWarnings.too_short_time_period(metric, len(variable), kwargs["min_time_steps"], inspect__stack())
            keyerror2 = "too short time period (" + str(len(variable)) + ")"
    # Read areacell & mask
    variable, areacell, keyerror3 = read_mask_area(
        variable, name_data, file_data, type_data, region, file_area=file_area, name_area=name_area,
        file_mask=file_mask, name_mask=name_mask, maskland=maskland, maskocean=maskocean, debug=debug, **kwargs)
    if keyerror1 is not None or keyerror2 is not None or keyerror3 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2, keyerror3])
    else:
        keyerror = None
    return variable, areacell, keyerror


def read_landmask_select_region(tab, file_data, name_data, file_mask, landmaskname='', box=None, **kwargs):
    """
    #################################################################################
    Description:
    Reads the given landmask from the given 'filename' and selects the given 'box'

    Uses cdms2 (uvcdat) to read areacell from 'filename' and cdutil (uvcdat) to select the 'box'
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable
    :param file_data: string
        string of the path to the file and name of the file that was used to read 'tab'
    :param name_data: string
        name of the variable, read from 'file_data', stored in 'tab'
    :param file_mask: string
        string of the path to the file and name of the file containing the landmask
    :param landmaskname: string, optional
        name of landmask (sftlf, lsmask, landmask,...) in 'filename'
    :param box: string, optional
        name of a region to select, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return area: masked_array
        masked_array containing landmask in 'box'
    """
    landmask = None
    # corrections for cdms2 to find the right axis
    CDMS2setAutoBounds("on")
    # list files
    list_files = [ff for ff in [file_mask, file_data] if ff is not None and os__path__isfile(ff) is True]
    # list variables
    list_variables = [landmaskname, "landmask", "lsmask", "sftlf"]
    for filename in list_files:
        # Open file and get time dimension
        fi = CDMS2open(filename)
        # dict region
        dict_region = dict()
        if box is not None:
            # define box
            region_ref = ReferenceRegions(box)
            dict_region["latitude"] = region_ref["latitude"]
            dict_region["longitude"] = region_ref["longitude"]
        # read land mask
        for var in list_variables:
            try:
                landmask = fi(var, **dict_region)
            except:
                pass
            else:
                break
        fi.close()
        if landmask is not None:
            break
    if landmask is None or tab.getGrid().shape != landmask.getGrid().shape:
        # Estimate landmask
        try:
            landmask = estimate_landmask(tab)
        except:
            arr = read_and_select_region(file_data, name_data, box="global",  **kwargs)
            try:
                landmask = estimate_landmask(arr)
            except:
                pass
        if landmask is not None:
            print("\033[93m" + str().ljust(25) + "NOTE: Estimated landmask applied" + "\033[0m")
        if landmask is not None and box is not None:
            # define box
            region_ref = ReferenceRegions(box)
            # subset
            landmask = landmask(latitude=region_ref["latitude"], longitude=region_ref["longitude"])
    # Return
    return landmask


def read_mask_area(tab, name_data, file_data, type_data, region, file_area="", name_area="", file_mask="", name_mask="",
                   maskland=False, maskocean=False, debug=False, **kwargs):
    tab_out = copy.copy(tab)
    keyerror1, keyerror2 = None, None
    # Read areacell
    if file_area:
        areacell = read_area_select_region(file_area, areaname=name_area, box=region, **kwargs)
    else:
        areacell = read_area_select_region(file_data, areaname=name_area, box=region, **kwargs)
    if areacell is not None and tab.getGrid().shape != areacell.getGrid().shape:
        areacell = None
    if debug is True:
        if areacell is not None:
            dict_debug = {"axes1": "(" + type_data + ") " + str([ax.id for ax in areacell.getAxisList()]),
                          "shape1": "(" + type_data + ") " + str(areacell.shape)}
            EnsoErrorsWarnings.debug_mode("\033[93m", "after read_area_select_region", 20, **dict_debug)
        else:
            dict_debug = {"line1": "areacell is None"}
            EnsoErrorsWarnings.debug_mode("\033[93m", "after read_area_select_region", 20, **dict_debug)
    # Read landmask
    if name_data.lower() in [
        "adt", "eta_t", "evap_heat", "height", "latent_heatflux", "lhf", "lw_heat", "lwr", "meridional_wind_stress",
        "msla", "net_heating", "net_longwave_heatflux_downwards", "net_shortwave_heatflux_downwards",
        "net_surface_heatflux_downwards", "netflux", "sea_surface_height", "sea_surface_temperature", "sens_heat",
        "sensible_heatflux", "shf", "sla", "sohefldo", "solatent", "solongwa", "sometauy", "sosensib", "soshfldo",
        "sossheig", "sosstsst", "sozotaux", "ssh", "sshg", "sst", "swflx", "swr", "tau_x", "tau_y", "tauuo", "tauvo",
        "taux", "tauy", "thf", "thflx", "tmpsf", "tos", "zonal_wind_stress", "zos"] \
            and "_Amon_" not in file_data and "oisst" not in file_data.lower():
        landmask = None
    elif name_data.lower() in ["pr", "slp"] and "_Omon_" in file_data:
        landmask = None
    elif maskland is False and maskocean is False:
        landmask = None
    else:
        landmask = read_landmask_select_region(tab, file_data, name_data, file_mask, landmaskname=name_mask, box=region,
                                               **kwargs)
    if debug is True:
        if landmask is not None:
            dict_debug = {"axes1": "(" + type_data + ") " + str([ax.id for ax in landmask.getAxisList()]),
                          "shape1": "(" + type_data + ") " + str(landmask.shape)}
            EnsoErrorsWarnings.debug_mode("\033[93m", "after ReadLandmaskSelectRegion", 20, **dict_debug)
        else:
            dict_debug = {"line1": "landmask is None"}
            EnsoErrorsWarnings.debug_mode("\033[93m", "after ReadLandmaskSelectRegion", 20, **dict_debug)
    # Apply landmask
    if landmask is not None:
        tab_out, keyerror1 = apply_landmask(tab_out, landmask, maskland=maskland, maskocean=maskocean)
        if keyerror1 is None:
            if areacell is None:
                areacell = create_array_ones(landmask, mid="areacell")
            areacell, keyerror2 = apply_landmask_to_area(areacell, landmask, maskland=maskland, maskocean=maskocean)
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
    return tab_out, areacell, keyerror


def save_netcdf(netcdf_name, var1=None, var1_attributes={}, var1_name="", var1_time_name=None, var2=None,
               var2_attributes={}, var2_name="", var2_time_name=None, var3=None, var3_attributes={}, var3_name="",
               var3_time_name=None, var4=None, var4_attributes={}, var4_name="", var4_time_name=None, var5=None,
               var5_attributes={}, var5_name="", var5_time_name=None, var6=None, var6_attributes={}, var6_name="",
               var6_time_name=None, var7=None, var7_attributes={}, var7_name="", var7_time_name=None, var8=None,
               var8_attributes={}, var8_name="", var8_time_name=None, var9=None, var9_attributes={}, var9_name="",
               var9_time_name=None, var10=None, var10_attributes={}, var10_name="", var10_time_name=None, var11=None,
               var11_attributes={}, var11_name="", var11_time_name=None, var12=None, var12_attributes={}, var12_name="",
               var12_time_name=None, frequency="monthly", global_attributes={}, **kwargs):
    if os__path__isdir(ntpath.dirname(netcdf_name)) is not True:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": given path does not exist",
            str().ljust(5) + "netcdf_name = " + str(netcdf_name)]
        EnsoErrorsWarnings.my_error(list_strings)
    if os__path__isfile(netcdf_name) is True:
        o = CDMS2open(netcdf_name, "a")
    else:
        o = CDMS2open(netcdf_name, "w+")
    if var1 is not None:
        if var1_name == "":
            var1_name = var1.id
        if var1_time_name is not None:
            var1 = time_but_not_time(var1, var1_time_name, frequency)
        o.write(var1, attributes=var1_attributes, dtype="float32", id=var1_name)
    if var2 is not None:
        if var2_name == "":
            var2_name = var2.id
        if var2_time_name is not None:
            var2 = time_but_not_time(var2, var2_time_name, frequency)
        o.write(var2, attributes=var2_attributes, dtype="float32", id=var2_name)
    if var3 is not None:
        if var3_name == "":
            var3_name = var3.id
        if var3_time_name is not None:
            var3 = time_but_not_time(var3, var3_time_name, frequency)
        o.write(var3, attributes=var3_attributes, dtype="float32", id=var3_name)
    if var4 is not None:
        if var4_name == "":
            var4_name = var4.id
        if var4_time_name is not None:
            var4 = time_but_not_time(var4, var4_time_name, frequency)
        o.write(var4, attributes=var4_attributes, dtype="float32", id=var4_name)
    if var5 is not None:
        if var5_name == "":
            var5_name = var5.id
        if var5_time_name is not None:
            var5 = time_but_not_time(var5, var5_time_name, frequency)
        o.write(var5, attributes=var5_attributes, dtype="float32", id=var5_name)
    if var6 is not None:
        if var6_name == "":
            var6_name = var6.id
        if var6_time_name is not None:
            var6 = time_but_not_time(var6, var6_time_name, frequency)
        o.write(var6, attributes=var6_attributes, dtype="float32", id=var6_name)
    if var7 is not None:
        if var7_name == "":
            var7_name = var7.id
        if var7_time_name is not None:
            var7 = time_but_not_time(var7, var7_time_name, frequency)
        o.write(var7, attributes=var7_attributes, dtype="float32", id=var7_name)
    if var8 is not None:
        if var8_name == "":
            var8_name = var8.id
        if var8_time_name is not None:
            var8 = time_but_not_time(var8, var8_time_name, frequency)
        o.write(var8, attributes=var8_attributes, dtype="float32", id=var8_name)
    if var9 is not None:
        if var9_name == "":
            var9_name = var9.id
        if var9_time_name is not None:
            var9 = time_but_not_time(var9, var9_time_name, frequency)
        o.write(var9, attributes=var9_attributes, dtype="float32", id=var9_name)
    if var10 is not None:
        if var10_name == "":
            var10_name = var10.id
        if var10_time_name is not None:
            var10 = time_but_not_time(var10, var10_time_name, frequency)
        o.write(var10, attributes=var10_attributes, dtype="float32", id=var10_name)
    if var11 is not None:
        if var11_name == "":
            var11_name = var11.id
        if var11_time_name is not None:
            var11 = time_but_not_time(var11, var11_time_name, frequency)
        o.write(var11, attributes=var11_attributes, dtype="float32", id=var11_name)
    if var12 is not None:
        if var12_name == "":
            var12_name = var12.id
        if var12_time_name is not None:
            var12 = time_but_not_time(var12, var12_time_name, frequency)
        o.write(var12, attributes=var12_attributes, dtype="float32", id=var12_name)
    my_keys = sorted(
        [key for key in list(kwargs.keys()) if "var" in key and str(key.replace("var", "")).isdigit() is True],
        key=lambda v: v.upper())
    for key in my_keys:
        if kwargs[key] is not None:
            if key + "_name" not in list(kwargs.keys()) or (
                    key + "_name" in list(kwargs.keys()) and kwargs[key + "_name"] == ""):
                kwargs[key + "_name"] = kwargs[key].id
            if key + "_time_name" in list(kwargs.keys()) and kwargs[key + "_time_name"] is not None:
                kwargs[key] = time_but_not_time(kwargs[key], kwargs[key + "_time_name"], frequency)
            if key + "_attributes" not in list(kwargs.keys()):
                kwargs[key + "_attributes"] = {}
            o.write(kwargs[key], attributes=kwargs[key + "_attributes"], dtype="float32", id=kwargs[key + "_name"])
    for att in list(global_attributes.keys()):
        o.__setattr__(att, global_attributes[att])
    o.close()
    return


def regridder(tab_to_regrid, newgrid, missing=None, order=None, mask=None, regridder="cdms", regridTool="esmf",
           regridMethod="linear", **kwargs):
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
            regridTool='regrid2' -> "linear"
            regridTool='esmf'    -> 'conserve', "linear", 'patch'
            regridTool='libcf'   -> "linear"
        'cdmsHorizontal' -> None
        default value is "linear"

    usual kwargs:
    :param newgrid_name: string, optional
        if "newgrid" is not defined (is string) this will be used to create a grid:
            this string must contain two keywords: the grid type and the resolution (same resolution in lon and lat)
            grid type  -> 'equalarea', "gaussian", 'generic', 'uniform'
            resolution -> '0.25x0.25deg', '0.5x0.5deg', '1x1deg', '2x2deg'
            e.g., newgrid_name='gaussian 1x1deg'
        default value is 'generic 1x1deg'
    :param region: string, optional
        if "newgrid" is not defined (is string) this will be used to create a grid
        name of a region, domain where the grid will be defined, must be defined in EnsoCollectionsLib.ReferenceRegions

    :return new_tab: masked_array
        tab_to_regrid regridded on newgrid
    """
    known_args = {"newgrid_name", "region"}
    extra_args = set(kwargs) - known_args
    if extra_args:
        EnsoErrorsWarnings.unknown_key_arg(extra_args, inspect__stack())
    # test given arguments
    known_regridder = ["cdms", "cdmsHorizontal"]
    if regridder not in known_regridder:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": regridder",
            str().ljust(5) + "unknown regridder: " + str(regridder),
            str().ljust(10) + "known regridder: " + str(known_regridder)]
        EnsoErrorsWarnings.my_error(list_strings)
    elif regridder == "cdms":
        if regridTool in ["regrid2", "libcf"]:
            list_method = [None, "linear"]
        elif regridTool == "esmf":
            list_method = [None, "conserve", "linear", "patch"]
        if (regridTool is not None) and (regridTool not in ["regrid2", "esmf", "libcf"]):
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": regridTool",
                str().ljust(5) + "unknown regridTool: " + str(regridTool),
                str().ljust(10) + "known regridTool: " + str(["regrid2", "esmf", "libcf"])]
            EnsoErrorsWarnings.my_error(list_strings)
        elif regridMethod not in list_method:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": regridMethod",
                str().ljust(5) + "unknown regridMethod (" + str(regridMethod) + ") for this regridTool ("
                + str(regridTool) + ")",
                str().ljust(10) + "known regridMethod: " + str(list_method)]
            EnsoErrorsWarnings.my_error(list_strings)
    # test the given "newgrid"
    if isinstance(newgrid, str) or newgrid is None:
        #
        # newgrid is not a grid, so a grid will be created
        # to do this, kwargs["newgrid_name"] and kwargs['region'] must be defined
        #
        # define the grid type
        GridType = "generic"
        for gtype in ["equalarea", "gaussian", "generic", "uniform"]:
            if gtype in kwargs["newgrid_name"]:
                GridType = gtype
                break
        # define resolution (same resolution in lon and lat)
        GridRes = 1.
        for res in ["0.25x0.25deg", "0.5x0.5deg", "0.75x0.75deg", "1x1deg", "1.25x1.25deg", "1.5x1.5deg",
                    "1.75x1.75deg", "2x2deg", "2.25x2.25deg", "2.5x2.5deg", "2.75x2.75deg"]:
            if res in kwargs["newgrid_name"]:
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
        idname = copy.copy(axis.id)
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
    else:
        regridFCT = REGRID2horizontal__Horizontal(tab_to_regrid.getGrid(), newgrid)
        new_tab = regridFCT(tab_to_regrid)
    return new_tab


def smoother_gaussian(tab, axis=0, window=5):
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
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)
    # Reorder tab in order to put 'axis' in first position
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder + str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the gaussian weight array
    weight = list()
    for ii in range(window):
        ii = ii - degree + 1
        frac = ii / float(window)
        gauss = float(1. / numpy__exp((4 * frac) ** 2))
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(gauss)
        weight.append(ww)
        del frac, gauss, ww
    weight = MV2array(weight)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1 * tmp2, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0) / window > 0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
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


def smoother_square(tab, axis=0, window=5):
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
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)

    # Reorder tab in order to put 'axis' in first position
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder + str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the weight array (uniform)
    weight = MV2ones((window,) + new_tab.shape[1:])

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0) / window > 0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
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


def smoother_triangle(tab, axis=0, window=5):
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
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)

    # Reorder tab in order to put 'axis' in first position
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder + str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the weight array (triangle)
    weight = list()
    for ii in range(0, (2 * degree) + 1):
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(float(1 + degree - abs(degree - ii)))
        weight.append(ww)
        del ww
    weight = MV2array(weight)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    sum_weight = MV2sum(weight, axis=0)
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1 * tmp2, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0) / window > 0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
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


dict_smooth = {"gaussian": smoother_gaussian, "square": smoother_square, "triangle": smoother_triangle}


def smoother(tab, info, axis=0, window=5, method="triangle"):
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
            "gaussian": gaussian shaped window
            "square":   square shaped window
            "triangle": triangle shaped window
    :return: smoothed_tab: masked_array
        smoothed data
    """
    try:
        dict_smooth[method]
    except:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(inspect__stack()) + ": smoothing method (running mean)",
            str().ljust(5) + "unkwown smoothing method: " + str(method),
            str().ljust(10) + "known smoothing method: " + str(
                sorted(list(dict_smooth.keys()), key=lambda v: v.upper()))]
        EnsoErrorsWarnings.my_error(list_strings)
    info = info + ', smoothing using a ' + str(method) + ' shaped window of ' + str(window) + ' points'
    return dict_smooth[method](tab, axis=axis, window=window), info


def time_but_not_time(tab, new_time_name, frequency):
    tab_out = copy.copy(tab)
    time_num = get_num_axis(tab_out, "time")
    timeax = tab_out.getAxis(time_num).asComponentTime()
    year1, month1, day1 = timeax[0].year, timeax[0].month, timeax[0].day
    if frequency == "daily":
        freq = "days"
    elif frequency == "monthly":
        freq = "months"
    elif frequency == "yearly":
        freq = "years"
    else:
        freq = None
        EnsoErrorsWarnings.unknown_frequency(frequency, inspect__stack())
    axis = CDMS2createAxis(list(range(len(tab_out))), id=new_time_name)
    axis.units = freq + " since " + str(year1) + "-" + str(month1) + "-" + str(day1)
    axis.axis = freq
    tab_out.setAxis(time_num, axis)
    return tab_out
