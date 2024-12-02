# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools built over xarray and xcdat
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
from dataclasses import dataclass
from inspect import stack as inspect__stack
import logging
from re import split as re__split
from typing import Annotated, Any, Hashable, Literal, Union
# numpy
from numpy import array as numpy__array
from numpy import ndarray as numpy__ndarray
# regionmask
import regionmask

# local functions
from enso_metrics.tools.dictionary_tools import combine_dict_levels, put_in_dict, sort_dict
from enso_metrics.wrapper.xarray_base import array_wrapper, dataset_wrapper
from enso_metrics.wrapper import xarray_base
from enso_metrics.wrapper import xcdat_base
# ---------------------------------------------------#
log = logging.getLogger(__name__)
ds_error = " should be xarray.DataArray, xarray.Dataset"


# ---------------------------------------------------------------------------------------------------------------------#
# Classes
# ---------------------------------------------------------------------------------------------------------------------#
@dataclass
class ValueRange:
    min: float
    max: float
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def apply_mask(
        ds: Union[array_wrapper, dataset_wrapper],
        ds_mask: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        data_var_mask: str = None,
        mask_land: bool = True,
        mask_ocean: bool = False,
        mask_threshold: Annotated[float, ValueRange(0.0, 1.0)] = 0.8,
        **kwargs) -> array_wrapper:
    """
    Mask land and/or ocean cells of given xarray.DataArray's or xarray.Dataset[data_var]'s data.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset to mask
    :param ds_mask: xarray.DataArray or xarray.Dataset
        Land-sea mask
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    :param data_var_mask: str, optional
        Data variable in ‘ds_mask‘ if it is xarray.Dataset; e.g., data_var_mask = "ts".
        If ‘ds_mask‘ is xarray.Dataset, ‘data_var_mask‘ must be provided.
        Default is None
    :param mask_land: bool, optional
        If True, cells that are mostly land are masked; e.g., mask_land = True.
        Default is True
    :param mask_ocean: bool, optional
        If True, cells that are mostly ocean are masked; e.g., mask_ocean = False.
        Default is False
    :param mask_threshold: float, optional
        Threshold used to mask data, must be within [0, 1]; e.g., mask_threshold = 0.8.
        If ‘mask_land’ is True and ’mask_threshold’ is 0.8, ds is masked where ds_mask > 0.8.
        If ‘mask_ocean’ is True and ’mask_threshold’ is 0.8, ds is masked where ds_mask < 1 - 0.8.
        Default is 0.8
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray masked where desired.
    """
    # select data and copy it
    da_o = xarray_base.copy(ds, data_var=data_var)
    ls_o = xarray_base.copy(ds_mask, data_var=data_var_mask)
    # ds_mask maximum value
    ds_mask_max = max_global(ls_o)
    # land-sea mask units
    ds_mask_units = ""
    if "units" in xarray_base.get_attributes_keys(ls_o):
        ds_mask_units = xarray_base.get_attribute(ls_o, "units")
    # if land = 100 instead of 1, divides ds_mask by 100
    if ds_mask_max > 1 or ds_mask_units == "%":
        ls_o /= 100.
    # mask cells that are mostly land
    if mask_land is True:
        da_o = xarray_base.where(da_o, ds_mask <= mask_threshold)
    # msk cells that are mostly ocean
    if mask_ocean is True:
        da_o = xarray_base.where(da_o, ds_mask >= 1 - mask_threshold)
    return da_o


def average_spatiotemporal(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var_ds: str,
        cf_dim: list[Literal["T", "X", "Y"]],
        ds_area: Union[array_wrapper, dataset_wrapper, None] = None,
        data_var_area: str = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    ds_o = None
    # select case
    if isinstance(cf_dim, list) is True and len(cf_dim) == 1 and cf_dim[0] == "T":
        # temporal average
        ds_o = xcdat_base.average_temporal(ds, data_var_ds, weighted=False)
    elif isinstance(cf_dim, list) is True and (
            (len(cf_dim) == 1 and cf_dim[0] in ["X", "Y", "Z"]) or
            (len(cf_dim) == 2 and cf_dim == ["X", "Y"])):
        # generate areacell weights
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)) is True:
            weights = compute_weights(
                ds, cf_dim=cf_dim, data_var=data_var_ds, ds_area=ds_area, data_var_area=data_var_area)
        else:
            weights = "generate"
        # spatial average
        ds_o = xcdat_base.average_spatial(ds, data_var_ds, cf_dim=cf_dim, weights=weights)
    return ds_o


def check_time_bounds(
        ds: Union[array_wrapper, dataset_wrapper],
        dim: Union[Hashable, str],
        time_bounds: Union[tuple[str, str]],
        side: Literal["lower", "upper"],
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    fill_log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "dim": dim, "time_bounds": time_bounds, "side": side})
    # bound position (0 or 1)
    position = 0 if side == "lower" else 1
    # split time bound in ["year", "month", "day"]
    split_ta = re__split("-| |:", str(xarray_base.get_time_bounds(ds)[position]))
    split_td = re__split("-| |:", time_bounds[position])
    fill_log_debug(inspect__stack(), "available vs. desired time_bounds: " + str(side), details={
        "available time_bound": split_ta, "desired time_bound": split_td})
    while (side == "lower" and all([float(k1) >= float(k2) for k1, k2 in zip(split_ta, split_td)]) is False) or \
            (side == "upper" and all([float(k1) <= float(k2) for k1, k2 in zip(split_ta, split_td)]) is False):
        # if side == "lower" & one of available("year", "month", "day") is smaller than desired ("year", "month", "day")
        # or
        # if side == "upper" & one of available("year", "month", "day") is larger than desired ("year", "month", "day")
        # remove the first (last) time step
        time_slice = slice(1, int(1e20)) if side == "lower" else slice(0, -1)
        ds = xarray_base.select_index(ds, {dim: time_slice})
        split_ta = re__split("-| |:", str(xarray_base.get_time_bounds(ds)[position]))
        fill_log_debug(inspect__stack(), "xarray_base.select_index", details={
            "available time_bound": split_ta, "desired time_bound": split_td})
    fill_log_debug(inspect__stack(), "output", details={
        "available time_bound": split_ta, "desired time_bound": split_td})
    return ds


def compute_variable(
        dict_ds: dict[str, dataset_wrapper],
        dict_input: dict[str, Union[str, list[str], None, dict[str, Union[int, float]]]],
        list_variables: list[str],
        variable: str,
        variables_param: dict[str, str],
        dataset: str = "unknown",
        experiment: str = "unknown",
        member: str = "unknown",
        project: str = "unknown",
        **kwargs) -> (Union[dataset_wrapper, None], dict):
    attributes_global = ["grid_label", "source", "version"]
    attributes_variable = ["units"]
    arr, dict_a = None, {}
    # loop on variable names to combine
    for nn in list_variables:
        # get offset and scale_factor from input dictionary
        add_offset = dict_input["variable_offset"][nn]
        scale_factor = dict_input["variable_scaling"][nn]
        # get array
        if nn not in list(dict_ds.keys()):
            arr = None
            break
        arr_t = xarray_base.to_array(dict_ds[nn], nn)
        # adjust array
        arr_t = arr_t * scale_factor + add_offset
        # combine
        if arr is None:
            arr = arr_t
        else:
            arr += arr_t
        # get global attributes
        att = xarray_base.get_attributes(dict_ds[nn])
        for kk in attributes_global:
            val = att[kk] if kk in list(att.keys()) else (dataset if kk == "source" else "unknown")
            put_in_dict(dict_a, val, kk, nn)
        # get variable attributes
        att = xarray_base.get_attributes(dict_ds[nn], data_var=nn)
        for kk in attributes_variable:
            val = att[kk] if kk in list(att.keys()) else "unknown"
            put_in_dict(dict_a, val, kk, nn)
    metadata, ds = {}, None
    if arr is not None:
        # metadata
        dict_t = {"input": {}, "variable": {}}
        for k1 in list(dict_t.keys()):
            if k1 == "input":
                for k2 in list(dict_a.keys()):
                    # check if given metadata is the same for all variables
                    val = list(set(list(dict_a[k2].values())))
                    # keep only one value is possible
                    val = val[0] if len(val) == 1 else deepcopy(dict_a[k2])
                    put_in_dict(dict_t, val, k1, variable, k2)
                put_in_dict(dict_t, dict_input["variable_computation"], k1, variable, "computation")
                val = str(list_variables[0]) if len(list_variables) == 1 else deepcopy(list_variables)
                put_in_dict(dict_t, val, k1, variable, "out_name")
            else:
                # metadata 'variable'
                val = dict((k2, k3) for k2, k3 in variables_param.items() if k2 != "variable_type")
                put_in_dict(dict_t, val, k1, variable)
        # sort metadata dictionary
        metadata, _, _ = sort_dict(dict_t)
        # array to dataset
        ds = xarray_base.to_dataset(arr, variable)
        # set netCDF variable attributes
        dict_t, _, _ = combine_dict_levels(metadata)
        xarray_base.set_attributes_variable(ds, variable, **dict_t)
        # set netCDF global attributes
        dict_t = {
            "DISCLAIMER": "The results in this file were produced with the PMP 3.6.1 (https://github.com/PCMDI/" +
                          "pcmdi_metrics). They are for research purposes only. They are subject to ongoing quality " +
                          "control and change as the PMP software advances, interpolation methods are modified, " +
                          "observational data sets are updated, problems with model data are corrected, etc. Use of " +
                          "these results for research (presentation, publications, etc.) should reference: 1) Lee et " +
                          "al. 2024, https://doi.org/10.5194/gmd-17-3919-2024; 2) Planton et al. 2021, https://doi." +
                          "org/10.1175/BAMS-D-19-0337.1. If any problems are uncovered in using these results please " +
                          "contact the PMP development team at pcmdi-metrics@llnl.gov",
            "REFERENCE": "CLIVAR ENSO metrics package (https://github.com/CLIVAR-PRP/ENSO_metrics)",
        }
        xarray_base.set_attributes_global(ds, **dict_t)
    return ds, metadata


def compute_weights(
        ds: Union[array_wrapper, dataset_wrapper],
        cf_dim: list[Literal["T", "X", "Y"]] = None,
        data_var: str = None,
        ds_area: Union[array_wrapper, dataset_wrapper, None] = None,
        data_var_area: str = None,
        **kwargs) -> array_wrapper:
    """
    Compute weights (or try to) o average along given dimension(s).
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param cf_dim: list[{"T", "X", "Y"}], optional
        List of CF dimension(s) for which weights are needed; e.g., cf_dim = ["X"].
        Default is None (i.e., ["X", "Y"])
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    :param ds_area: xarray.DataArray or xarray.Dataset, optional
        DataArray or Dataset of areacell
    :param data_var_area: str, optional
        Data variable in ‘ds_area‘ if it is xarray.Dataset; e.g., data_var_area = "ts".
        If ‘ds_area‘ is xarray.Dataset, ‘data_var_area‘ must be provided.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray containing the weights to use during averaging.
    """
    if isinstance(cf_dim, list) is False:
        cf_dim = ["X", "Y"]
    # by defaults weights are equal to 1
    weights = create_array(ds, data_var=data_var, data_var_o="weights", value=1)
    # temporal or spatial weights?
    if isinstance(cf_dim, list) is True and len(cf_dim) == 1 and cf_dim[0] == "T":
        # try to generate time weights using xcdat
        try:
            weights = xcdat_base.weights_temporal(ds, data_var)
        except "weights to 1":
            pass
    else:
        # is area available?
        da_area = xarray_base.copy(ds_area, data_var=data_var_area)
        if isinstance(da_area, array_wrapper) is False:
            # try to generate spatial weights using xcdat
            try:
                weights = xcdat_base.weights_spatial(ds, data_var, cf_dim=cf_dim)
            except "force weights to 1":
                pass
        else:
            # get the name of the given axis (T, X, Y, Z)
            dim_lon = xcdat_base.get_axis_key(ds_area, "X")
            dim_lat = xcdat_base.get_axis_key(ds_area, "Y")
            if set(cf_dim) == {"X", "Y"}:
                # compute area summed along given dimensions
                total_area = xarray_base.sum_along_axis(da_area, dim=[dim_lat, dim_lon])
            else:
                # dimension to average
                dim_name = deepcopy(dim_lon if set(cf_dim) == {"X"} else dim_lat)
                # get dimension array and position in matrix (i.e., axis)
                dim_array = da_area[dim_name]
                dim_axis = xarray_base.get_dim_keys(da_area).index(dim_name)
                # compute area summed along given dimension
                total_area = xarray_base.sum_along_axis(da_area, dim=[dim_name])
                # expand array (i.e., recreate initial array shape)
                total_area = xarray_base.expand_dim(total_area, axis=dim_axis, dim={dim_name: dim_array})
            # operation: weights == 1 when summed along given dimension(s)
            weights = da_area / total_area
    return weights


def create_array(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        data_var_o: str = "",
        value: Union[float, int, None] = 0,
        **kwargs) -> array_wrapper:
    """
    Return a new DataArray of ’value’ with the same shape, axes, coordinates, attributes,... as input DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.zeros_like.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    :param data_var_o: str, optional
        Name of the output data variable.
        Default is ""
    :param value: float or int or None, optional
        Value to fill the array with; e.g., value = 1.
        If None, array filled with NaN.
        Default is 0
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray of ones with the same shape and type as ds.
    """
    # create an array of zeros with the same shape, axes, coordinates, attributes,... as input DataArray
    da_o = xarray_base.create_array_zero(ds, data_var)
    # name this new array
    da_o.name = data_var_o
    if value is None:
        da_o = xarray_base.where(da_o, da_o == 1)
    elif isinstance(value, (float, int)) is True and value != 0:
        da_o = xarray_base.where(da_o, da_o == value, other=value)
    return da_o


def create_land_sea_mask(
        ds: Union[array_wrapper, dataset_wrapper],
        mask_as_boolean: bool = False,
        **kwargs) -> array_wrapper:
    """
    Create a land-sea mask (1: land, 0: sea) for given dataset, using regionmask.
    https://regionmask.readthedocs.io/en/stable/defined_landmask.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset for which a land-sea mask must be created
    :param mask_as_boolean : bool, optional
        Define mask as (1: land, 0: sea) or (True: land, False: sea); e.g., mask_as_boolean = False.
        If True, define cells as boolean.
        If False, define cells as float.
        Default is False
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray of land-sea mask (1: land, 0: sea) or (True: land, False: sea).
    """
    # get latitude and longitude names from the input dataset or array
    dim_lon = xcdat_base.get_axis_key(ds, "X")
    dim_lat = xcdat_base.get_axis_key(ds, "Y")
    # get latitude and longitude arrays
    dim_lon_array = xarray_base.get_dim_array(ds, dim_lon)
    dim_lat_array = xarray_base.get_dim_array(ds, dim_lat)
    # create a land-sea mask using regionmask
    land_sea_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110
    # adapt land-sea mask to input dataset or array
    land_sea_mask = land_sea_mask.mask(dim_lon_array, dim_lat_array)
    if mask_as_boolean is True:
        # convert the land-sea mask to a boolean mask
        land_sea_mask = xarray_base.change_type(land_sea_mask, "bool")
    else:
        # convert masked cells (land) to 1
        land_sea_mask = xarray_base.where(land_sea_mask, land_sea_mask == 0, other=1)
    return land_sea_mask


def fill_log_debug(
        stack,
        message: str,
        adjust: int = None,
        data_var: str = None,
        details: dict[str, Any] = None,
        ds: Union[array_wrapper, dataset_wrapper, None] = None):
    if isinstance(adjust, int) is False:
        adjust = 20
    message = " function " + str(stack[0][3]) + " ; line " + str(stack[0][2]) + "\n" + str(message)
    tmp_dict = {}
    if isinstance(ds, dataset_wrapper) is True and isinstance(data_var, str) is True and \
            data_var in xarray_base.get_dataset_keys(ds):
        tmp_dict = {
            "variable": str(data_var), "dataset_keys": xarray_base.get_dataset_keys(ds),
            "dim_keys": xarray_base.get_dim_keys(ds), "shape": xarray_base.get_array_shape(ds, data_var),
            "min": min_global(ds, data_var), "max": max_global(ds, data_var)}
    if isinstance(ds, array_wrapper) is True:
        tmp_dict = {
            "array_name": xarray_base.get_array_name(ds), "dim_keys": xarray_base.get_dim_keys(ds),
            "shape": xarray_base.get_array_shape(ds), "min": min_global(ds), "max": max_global(ds)}
    for k1, k2 in tmp_dict.items():
        message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    if isinstance(details, dict) is True:
        for k1, k2 in details.items():
            message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    log.debug(message)


def fill_log_info(stack, message: str):
    log.info(" function " + str(stack[0][3]) + " ; line " + str(stack[0][2]) + "\n" + str(message))


def get_dim_array_latitude(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> array_wrapper:
    """
    Return latitude array of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Latitude DataArray.
    """
    # get latitude as named in xarray.DataArray or xarray.Dataset
    dim_name = xcdat_base.get_axis_key(ds, "Y")
    # get latitude array
    return xarray_base.get_dim_array(ds, dim_name)


def get_dim_array_longitude(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> array_wrapper:
    """
    Return longitude array of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Longitude DataArray.
    """
    # get longitude as named in xarray.DataArray or xarray.Dataset
    dim_name = xcdat_base.get_axis_key(ds, "X")
    # get longitude array
    return xarray_base.get_dim_array(ds, dim_name)


def get_dim_array_time(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> array_wrapper:
    """
    Return time array of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Time DataArray.
    """
    # get time as named in xarray.DataArray or xarray.Dataset
    dim_name = xcdat_base.get_axis_key(ds, "T")
    # get time array
    return xarray_base.get_dim_array(ds, dim_name)


def get_dim_array_vertical(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> array_wrapper:
    """
    Return vertical array of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Vertical DataArray.
    """
    # get vertical as named in xarray.DataArray or xarray.Dataset
    dim_name = xcdat_base.get_axis_key(ds, "Z")
    # get vertical array
    return xarray_base.get_dim_array(ds, dim_name)


def recreate_array(
        arr: numpy__ndarray,
        ds: Union[array_wrapper, dataset_wrapper],
        attrs_added: dict[str, str] = None,
        axis_added: Union[list[int], tuple[int], None] = None,
        coords_added: dict[str, Union[numpy__ndarray, array_wrapper]] = None,
        data_var: str = None,
        data_var_o: str = None,
        dim_added: Union[list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        dim_removed: Union[list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        **kwargs) -> array_wrapper:
    """
    Recreate output xarray.DataArray from input xarray.DataArray.
    
    Input:
    ------
    :param arr: numpy.ndarray
        Array derived from ‘ds‘ (e.g., a statistic was computed) that was transformed into a numpy.ndarray in the
        process
    :param ds: xarray.DataArray or xarray.Dataset
        Original xarray.DataArray or xarray.Dataset from which ‘arr‘ is derived
    :param attrs_added: dict[str, str] or None, optional
        Variable attributes to add to the DataArray; e.g, attrs_added = {"attr_name": "new attribute"}.
        If given, both ‘axis_added‘ and ‘dim_added‘ must be provided.
        Default is None (no dimension has been added)
    :param axis_added: list[int] or tuple[int] or None
        Position(s) of added dimension(s) (if any); e.g, axis_added = [0] or axis_added = [0, 1].
        If given, both ‘axis_added‘ and ‘dim_added‘ must be provided.
        Default is None (no dimension has been added)
    :param coords_added: dict[str, numpy.ndarray or xarray.DataArray] or None, optional
        Coordinates (tick labels) to use for indexing along each dimension.
        If ‘dim_added‘ is not in ‘coords_added‘, coordinates will be a sequence of numbers.
        Default is None
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    :param data_var_o: str, optional
        Name of the output data variable.
        Default is ""
    :param dim_added: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been added from ‘ds‘ to ‘arr‘ (if any);
        e.g, dim_added = ["x"] or dim_added = ["x", "y"].
        If given, both ‘axis_added‘ and ‘dim_added‘ must be provided.
        Default is None (no dimension has been added)
    :param dim_removed: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been removed from ‘ds‘ to ‘arr‘ (e.g., to compute a statistic);
        e.g., dim_removed = ["x"] or dim_removed = ["x", "y"].
        Default is None (no dimension was removed)
    
    Output:
    -------
    :return: xarray.DataArray
        Input ‘arr‘ wrapped in a xarray.DataArray.
    """
    # read variable from xarray.Dataset if needed
    ds = xarray_base.to_array(ds, data_var)
    # list dimensions
    dimensions = xarray_base.get_dim_keys(ds)
    # delete removed dimension(s)
    for k in dim_removed:
        if isinstance(k, str) is True and k in dimensions:
            dimensions.remove(dim_removed)
    # get coordinates corresponding to dimensions
    coordinates = dict((k, xarray_base.get_dim_array(ds, k)) for k in dimensions)
    # add given dimension(s)
    if isinstance(axis_added, (list, tuple)) is True and isinstance(dim_added, (list, tuple)) is True and \
            len(axis_added) == len(dim_added):
        for k1, k2 in zip(axis_added, dim_added):
            # add given dimension at given position
            dimensions.insert(k1, k2)
            # add coordinates in dictionary
            if isinstance(coords_added, dict) is True and k2 in coords_added.keys():
                coordinates[k2] = coords_added[k2]
            else:
                coordinates[k2] = numpy__array(list(range(arr.shape[k1])))
    # get input attributes
    attributes = xarray_base.get_attributes(ds)
    # add given attribute(s)
    attributes.update(attrs_added)
    attributes = dict((k, attributes[k]) for k in sorted(list(attributes.keys()), key=lambda v: v.lower()))
    # numpy.ndarray to xarray.DataArray
    return xarray_base.numpy_to_array(arr, coordinates, dimensions, attrs=attributes, name=data_var_o)


def max_global(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        **kwargs) -> float:
    """
    Get the maximum value of given object's array.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: float
        Maximum value.
    """
    # maximum value
    arr = xarray_base.maximum(ds, data_var=data_var)
    # to numpy
    arr = xarray_base.to_numpy(arr)
    # to float
    return float(arr)


def min_global(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        **kwargs) -> float:
    """
    Get the minimum value of given object's array.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: float
        Minimum value.
    """
    # minimum value
    arr = xarray_base.minimum(ds, data_var=data_var)
    # to numpy
    arr = xarray_base.to_numpy(arr)
    # to float
    return float(arr)


def min_max_global(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        **kwargs) -> tuple[float, float]:
    """
    Get the minimum and maximum values of given object's array.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: tuple[float, float]
        Minimum and maximum values.
    """
    # minimum and maximum values
    arr_min, arr_max = xarray_base.minimum(ds, data_var=data_var), xarray_base.maximum(ds, data_var=data_var)
    # to numpy
    arr_min, arr_max = xarray_base.to_numpy(arr_min), xarray_base.to_numpy(arr_max)
    # to float
    return float(arr_min), float(arr_max)


def processing_description(
        ds: dataset_wrapper,
        data_var: str,
        attribute_name: str,
        description: str,
        **kwargs) -> dataset_wrapper:
    # get all attributes
    attributes = xarray_base.get_attributes(ds, data_var=data_var)
    # select or create 'attribute_name'
    att_o = attributes[attribute_name] if attribute_name in list(attributes.keys()) else ""
    # processing number
    cc = 2
    if ";; 2) " in att_o:
        # split after ';; ' and then before ') ' to find the last processing number, then add 1
        cc = int(att_o.split(";; ")[-1].split(") ")[0]) + 1
    elif att_o == "":
        # att_o is an empty string so the processing number is 1
        cc = 1
    cc = str(cc) + ") " if cc == 1 else ";; " + str(cc) + ") "
    # update the processing description
    att_o += str(cc) + str(description)
    # update the processing description in the dataset
    xarray_base.set_attributes_variable(ds, attribute_name=att_o)
    return ds


def read_dataset(
        add_offset: Union[float, int] = 0,
        data_var: Union[str, None] = None,
        exponent: Union[float, int] = 1,
        filename: Union[str, list[str], None] = None,
        scale_factor: Union[float, int] = 1,
        time_bounds: Union[int, tuple[int], tuple[int, int], tuple[str, str], None] = None,
        kwargs_open_dataset: dict = None,
        kwargs_select_time: dict = None,
        **kwargs):
    if isinstance(kwargs_open_dataset, dict) is False:
        kwargs_open_dataset = {}
    fill_log_debug(inspect__stack(), "input", details={
        "add_offset": add_offset, "data_var": data_var, "exponent": exponent, "filename": filename,
        "scale_factor": scale_factor, **kwargs_open_dataset})
    ds = "error"
    # fake loop to be able to break out
    for _ in [0]:
        if isinstance(data_var, str) is False or isinstance(filename, (list, str)) is False:
            # file / variable to read not given
            details = {}
            if isinstance(data_var, str) is False:
                details["data_var"] = str(type(data_var)) + " should be str"
            if isinstance(filename, (list, str)) is False:
                details["filename"] = str(type(filename)) + " should be list[str] or str"
            fill_log_debug(inspect__stack(), "WARNING cannot read dataset", details=details)
            break
        # read dataset
        try:
            ds = xcdat_base.open_dataset(filename, data_var=data_var, **kwargs_open_dataset)
        except Exception as err:
            fill_log_debug(inspect__stack(), "ERROR cannot read dataset\n" + str(err))
            # no log: it is usual that the given dimension is not available and therefore squeeze is not applied
            break
        fill_log_debug(inspect__stack(), "xcdat_base.open_dataset", ds=ds)
        # select epoch
        if isinstance(time_bounds, (int, tuple)) is True:
            ds = select_time(ds, data_var=data_var, time_bounds=time_bounds, kwargs_select_time=kwargs_select_time)
            fill_log_debug(inspect__stack(), "select_time", ds=ds)
        # scale factor
        ds *= scale_factor
        fill_log_debug(inspect__stack(), "scale_factor", ds=ds)
        # add offset
        ds += add_offset
        fill_log_debug(inspect__stack(), "add_offset", ds=ds)
        # exponent
        ds **= exponent
        fill_log_debug(inspect__stack(), "exponent", ds=ds)
        # squeeze vertical dimension
        ds = squeeze_dimension(ds, "Z")
        fill_log_debug(inspect__stack(), "squeeze_dimension, Z", ds=ds)
        # mask values that are not within [-1e20, 1e20]
        ds = xarray_base.where(ds, ds > -1e20)
        fill_log_debug(inspect__stack(), "xarray_base.where, ds > -1e20", ds=ds)
        ds = xarray_base.where(ds, ds < 1e20)
        fill_log_debug(inspect__stack(), "xarray_base.where, ds < 1e20", ds=ds)
        # HadISST has -1000 values... mask values lower than 200
        if "hadisst" in filename.lower():
            ds = xarray_base.where(ds, ds > -200)
            fill_log_debug(inspect__stack(), "xarray_base.where, ds > -200", ds=ds)
    fill_log_debug(inspect__stack(), "output", data_var=data_var, ds=ds)
    return ds


def select_time(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: Union[str, None] = None,
        time_bounds: Union[int, tuple[int], tuple[int, int], tuple[str, str], None] = None,
        kwargs_select_time: dict = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Select epoch based on given time bounds.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str or None, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, log will have more details.
        Default is None
    :param time_bounds: int or tuple[int] or tuple[int, int] or tuple[str, str] or None, optional
        Time extent to select data; e.g., time_bounds = ("1980-01-01", "2014-12-31") or time_bounds = 120 or
        time_bounds = (0, 120).
        If ‘time_bounds‘ is int or tuple[int], the first ‘time_bounds‘ time steps will be selected.
        If ‘time_bounds‘ is tuple[int, int], time steps between the first int and second int will be selected.
        Default is None
    :param kwargs_select_time: dict or None, optional
        kwargs to be passed to xarray_base.select or xarray_base.select_index
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with data of each array selected within given ‘time_bounds‘.
    """
    if isinstance(kwargs_select_time, dict) is False:
        kwargs_select_time = {}
    fill_log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "data_var": data_var, "time_bounds": time_bounds, **kwargs_select_time})
    if time_bounds is None:
        # no time bounds to select data
        fill_log_debug(inspect__stack(), "time_bounds is None: NOT selected", details={"time_bounds": str(time_bounds)})
    elif isinstance(ds, (array_wrapper, dataset_wrapper)):
        # test if time dimension is available
        try:
            xcdat_base.get_axis_key(ds, "T")
        except (Exception,):
            fill_log_debug(inspect__stack(), "WARNING cannot select time_bounds", details={
                "dim_name": "time dimension not available", "dim_keys": xarray_base.get_dim_keys(ds)})
            pass
        else:
            # get dimension name
            dim_name = xcdat_base.get_axis_key(ds, "T")
            fill_log_debug(inspect__stack(), "xcdat_base.get_axis_key", data_var=data_var,
                           details={"dim_name": dim_name}, ds=ds)
            # check given time bounds type
            if isinstance(time_bounds, (list, tuple)) is True and len(time_bounds) == 2 and \
                    all([isinstance(k, str) for k in time_bounds]) is True:
                # select using time_bounds like ("1980-01-01", "2014-12-31")
                ds = xarray_base.select(ds, {dim_name: slice(*time_bounds)}, **kwargs_select_time)
                fill_log_debug(inspect__stack(), "xarray_base.select", data_var=data_var,
                               details={"available time_bounds": xarray_base.get_time_bounds(ds)}, ds=ds)
                # sometimes selecting time is slightly wrong
                # this section checks if one time step has not been included by error at the beginning or the end of the
                # time series
                # check lower time bound
                ds = check_time_bounds(ds, dim_name, time_bounds, "lower")
                # check upper time bound
                ds = check_time_bounds(ds, dim_name, time_bounds, "upper")
                fill_log_debug(inspect__stack(), "check_time_bounds", data_var=data_var,
                               details={"available time_bounds": xarray_base.get_time_bounds(ds)}, ds=ds)
            elif isinstance(time_bounds, int) is True or (
                    isinstance(time_bounds, (list, tuple)) is True and 0 < len(time_bounds) < 3 and
                    all([isinstance(k, int) for k in time_bounds]) is True):
                # select using time_bounds like (12, 24)
                tmp_time_bounds = deepcopy(time_bounds)
                if isinstance(time_bounds, int) is True:
                    tmp_time_bounds = (time_bounds,)
                ds = xarray_base.select_index(ds, {dim_name: slice(*tmp_time_bounds)}, **kwargs_select_time)
                fill_log_debug(inspect__stack(), "xarray_base.select_index", data_var=data_var, ds=ds)
            else:
                fill_log_debug(inspect__stack(), "WARNING cannot select time_bounds", details={
                    "time_bounds": str(time_bounds) + " should be int, tuple[str, str], tuple[int], tuple[int, int]"})
    else:
        # input is neither a dataset nor a dataarray
        fill_log_debug(inspect__stack(), "WARNING cannot select time_bounds", details={
                    "ds": str(type(ds)) + " should be xarray.DataArray, xarray.Dataset"})
    fill_log_debug(inspect__stack(), "output", data_var=data_var, ds=ds)
    return ds


def squeeze_dimension(
        ds: Union[array_wrapper, dataset_wrapper],
        cf_dim: Literal["T", "X", "Y", "Z"],
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    If given dimension is defined and is of length 1 in given object, squeeze it and remove associated bounds.
    Else, ignore.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param cf_dim: {"X", "Y", "T", "Z"}
        The CF axis (dimension) key
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        This object, but with given dimension, if defined and of length 1, removed.
    """
    fill_log_debug(inspect__stack(), "input", details={"ds.type": type(ds), "cf_dim": cf_dim})
    # fake loop to be able to break out
    for _ in [0]:
        if isinstance(ds, (array_wrapper, dataset_wrapper)) is False:
            # wrong input
            fill_log_debug(inspect__stack(), "WARNING cannot squeeze", details={"ds": str(type(ds)) + str(ds_error)})
            break
        # test is given dimension is available
        try:
            dim_name = xcdat_base.get_axis_key(ds, cf_dim)
        except (Exception,):
            # no log: it is usual that the given dimension is not available and therefore squeeze is not applied
            break
        else:
            fill_log_debug(inspect__stack(), "xcdat_base.get_axis_key", details={
                "dim_name": dim_name, "dim_keys": xarray_base.get_dim_keys(ds)})
            # test if given dimension can be squeezed
            if xarray_base.get_array_size(ds, dim_name) < 2:
                # squeeze dimension
                ds = xarray_base.squeeze(ds, dim=dim_name, drop=True)
                fill_log_debug(inspect__stack(), "xarray_base.squeeze", details={
                    "dim_keys": xarray_base.get_dim_keys(ds)})
                # remove bounds
                if isinstance(ds, dataset_wrapper) is True:
                    ds = xarray_base.drop_dataset_keys(
                        ds, [str(dim_name) + "_bounds", str(dim_name) + "_bnds"], errors="ignore")
                    fill_log_debug(inspect__stack(), "xarray_base.drop_dataset_keys", details={
                        "dataset_keys": xarray_base.get_dataset_keys(ds)})
    fill_log_debug(inspect__stack(), "output", details={"dim_keys": xarray_base.get_dim_keys(ds)})
    return ds
# ---------------------------------------------------------------------------------------------------------------------#
