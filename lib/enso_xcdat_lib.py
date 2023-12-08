# -*- coding:UTF-8 -*-

# basic python
from copy import deepcopy
from typing import Literal
# numpy
import numpy as np
# regionmask
import regionmask
# xarray
import xarray as xr
# xCDAT
import xcdat as xc
# ENSO_metrics
from . enso_xcdat_base import fct_array_ones, fct_averager, fct_get_axis_key, fct_get_latitude, \
    fct_get_longitude, fct_numpy_to_xarray, fct_sum, fct_standard_deviation, fct_uniform_grid, \
    fct_where_xarray, fct_where_dataarray, fct_xarray_to_numpy, \
    fct_horizontal_regrid, fct_attrs_update_global, fct_attrs_update_variable, fct_dataset_to_netcdf, \
    fct_max, fct_min, fct_get_attributes_list, fct_get_attribute


# ---------------------------------------------------------------------------------------------------------------------#
# xarray / xcdat data processing functions
# ---------------------------------------------------------------------------------------------------------------------#
def apply_landmask(ds: xr.Dataset, data_var: str, landmask: xr.DataArray, maskland : bool=True, maskocean: bool=False) -> xr.DataArray:
    da = ds[data_var].copy()
    # if land = 100 instead of 1, divides landmask by 100
    if (fct_min(landmask) == 0 and fct_max(landmask) == 100) or \
            ("units" in fct_get_attributes_list(landmask) and fct_get_attribute(landmask, "units") == "%"):
        landmask = landmask / 100.
    if maskland is True:
        da = fct_where_dataarray(da, landmask <= 0.8)
    if maskocean is True:
        da = fct_where_dataarray(da, landmask >= 0.2)       
    return da


def average_horizontal(ds: xr.Dataset, data_var: str, areacell: xr.DataArray = None) -> xr.Dataset:
    if areacell is None:
        weights = "generate"
    else:
        weights = generate_weights(ds, data_var, areacell, axis=["X", "Y"])
    return fct_averager(ds, data_var, ["X", "Y"], weights=weights)


def average_zonal(ds: xr.Dataset, data_var: str="ts", areacell: xr.DataArray=None) -> xr.Dataset:
    if areacell is None:
        weights = "generate"
    else:
        weights = generate_weights(ds, data_var, areacell, axis=["X"])
    return fct_averager(ds, data_var, ["X"], weights=weights)


def average_meridional(ds: xr.Dataset, data_var: str="ts", areacell: xr.DataArray=None) -> xr.Dataset:
    if areacell is None:
        weights = "generate"
    else:
        weights = generate_weights(ds, data_var, areacell, axis=["Y"])
    return fct_averager(ds, data_var, ["Y"], weights=weights)


def compute_standard_deviation(ds: xr.Dataset, data_var: str, areacell: xr.DataArray = None, ddof: int = 1,
                               axis: list[Literal["X", "Y", "T", "Z"]] = None):
    # get keys of each dimension
    if axis is None:
        # if axis is not given, compute the standard deviation on the time dimension
        axis_keys = [fct_get_axis_key(ds, "T")]
    else:
        axis_keys = [fct_get_axis_key(ds, k) for k in axis]
    # generate weights for given axes (dimensions)
    weights = generate_weights(ds, data_var, areacell, axis_keys)
    # multiply variable array with weights
    da = ds[data_var] * weights
    # compute the standard deviation along given axes
    return fct_standard_deviation(da, axis_keys, ddof)


def generate_land_sea_mask(ds: xr.Dataset, boolean: bool=False) -> xr.DataArray:
    """
    A function generates land sea mask (1: land, 0: sea) for given xarray Dataset,
    assuming the given xarray dataset and has latitude and longitude coordinates. 

    Parameters
    ----------
    ds : xr.Dataset
        A Dataset object.
    boolen : bool, optional
        Set mask value to True (land) or False (sea), by default False
        
    Returns
    -------
    xr.DataArray
        A DataArray of land sea mask (1: land, 0: sea)
    """
    # Create a land-sea mask using regionmask
    land_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110

    # Get the longitude and latitude from the xarray dataset  
    lon = fct_get_longitude(ds)
    lat = fct_get_latitude(ds)

    # Mask the land-sea mask to match the dataset's coordinates
    land_sea_mask = land_mask.mask(lon, lat)
    
    if not boolean:
        # Convert the land-sea mask to a boolean mask
        land_sea_mask = fct_where_xarray(land_sea_mask, 0, 1)  

    return land_sea_mask


def generate_uniform_grid(lat1: float, lat2: float, lon1: float, lon2: float,
                          target_grid: str = "1x1deg") -> xr.Dataset:
    # split resolution name
    resolution = target_grid.split("x")
    # remove 'deg'
    resolution = [k.replace("deg", "") for k in resolution]
    # str to float
    resolution = [float(k) for k in resolution]
    # get meridional and zonal resolutions
    lat_res = resolution[0]
    lon_res = resolution[1]
    # calculate the center of the first meridional and zonal grid points
    start_lat = lat1 + lat_res / 2.
    start_lon = lon1 + lon_res / 2.
    # generate uniform grid
    return fct_uniform_grid(start_lat, lat2, lat_res, start_lon, lon1, lon_res)


def generate_weights(ds: xr.Dataset, data_var: str, areacell: xr.DataArray = None,
                     axis: list[str] = None) -> xr.DataArray:
    if axis is None:
        axis = ["X", "Y"]
    if areacell is None or axis == ["T"]:
        # areacell is not provided so weights cannot be computed
        # or
        # time axis is used and so weights are not needed (each time step has the same weight)
        weights = fct_array_ones(ds[data_var])
    else:
        # get the name of the given axis (T, X, Y, Z)
        lon_key = fct_get_axis_key(areacell, "X")
        lat_key = fct_get_axis_key(areacell, "Y")
        if axis == ["X"] or axis == ["Y"]:
            # weights are generated for operations along a single axis
            if axis == ["X"]:
                # zonal operation: weights == 1 when summed along X axis
                total_area = fct_sum(areacell, lon_key)
                weights_numpy_2d = fct_xarray_to_numpy(areacell) / fct_xarray_to_numpy(total_area)[:, np.newaxis]
            else:
                # meridional operation: weights == 1 when summed along Y axis
                total_area = fct_sum(areacell, lat_key)
                weights_numpy_2d = areacell.to_numpy() / total_area.to_numpy()[np.newaxis, :]
            # create DataArray from numpy array
            weights = fct_numpy_to_xarray(weights_numpy_2d, (lat_key, lon_key),
                                            {lat_key: fct_get_latitude(ds), lon_key: fct_get_longitude(ds)})
        else:
            # horizontal operation: weights == 1 when summed
            total_area = fct_sum(areacell, (lat_key, lon_key))
            weights = areacell / total_area
    return weights


def save_netcdf(ds, output_file, list_of_variables: list[str]=None, list_of_variable_attributes: list[dict]=None, global_attributes: dict=None):
    if global_attributes is not None:
        fct_attrs_update_global(ds, global_attributes)
        
    if list_of_variables is not None:
        for (variable, variable_attributes) in zip(list_of_variables, list_of_variable_attributes):
            fct_attrs_update_variable(ds, variable, variable_attributes)
    fct_dataset_to_netcdf(ds, output_file, list_of_variables)



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
