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
from typing import Any, Hashable, Literal, Union
# numpy
from numpy import array as numpy__array
from numpy import ndarray as numpy__ndarray
# regionmask
import regionmask

# local functions
from enso_metrics.tools.default import set_instance
from enso_metrics.tools.dictionary_tools import combine_dict_levels, put_in_dict, sort_dict
from enso_metrics.wrapper import basics
from enso_metrics.wrapper import xarray_base
from enso_metrics.wrapper import xcdat_base
from enso_metrics.wrapper.xarray_base import array_wrapper, dataset_wrapper
# ---------------------------------------------------#
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
def average_spatiotemporal(
        ds: dataset_wrapper,
        data_var: str,
        cf_dim: list[Literal["T", "X", "Y"]],
        ds_area: Union[array_wrapper, dataset_wrapper, None] = None,
        data_var_area: str = None,
        kwargs_average_spatial: dict = None,
        kwargs_average_temporal: dict = None,
        **kwargs) -> dataset_wrapper:
    """
    Compute spatio-temporal average (as defined by ‘cf_dim’) using xcdat
    WARNING: xcdat only works with datasets

    :param ds: xarray.Dataset
    :param data_var: str
        Data variable in ‘ds’ on which to calculate averages; e.g., data_var = "ts"
    :param cf_dim: list[Literal["T", "X", "Y"]]
        List of cf_dim along which to calculate averages; e.g., cf_dim = ["X", "Y"]
    :param ds_area: xarray.Dataset, optional
        Area cell dataset
    :param data_var_area: str, optional
        Data variable in ‘ds_area’
    :param kwargs_average_spatial: dict, optional
        Key arguments to compute spatial average (see xcdat_base.average_spatial); e.g., kwargs_average_spatial = {}
        Default is None
    :param kwargs_average_temporal: dict, optional
        Key arguments to compute temporal average (see xcdat_base.average_temporal); e.g., kwargs_average_temporal = {}
        Default is None
    **kwargs - Discarded

    Output:
    -------
    :return: Dataset
        Input dataset with the average of the data variable.
    """
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds,
                   details={"ds.type": type(ds), "dim": cf_dim, "ds_area.type": type(ds_area)})
    ds_o = None
    # select case
    if isinstance(cf_dim, list) is True and "T" in cf_dim:
        # temporal average
        ds_o = xcdat_base.average_temporal(ds, data_var, **kwargs_average_temporal)
    if isinstance(cf_dim, list) is True and len(list(set(cf_dim) & {"X", "Y"})) > 0:
        cf_dim: list[Literal["X", "Y"]] = ["X"]
        if len(list({"X", "Y"} - set(cf_dim))) == 0:
            cf_dim = ["X", "Y"]
        elif len(list({"Y"} - set(cf_dim))) == 0:
            cf_dim = ["Y"]
        # generate areacell weights
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)) is True:
            weights = compute_weights(
                ds, cf_dim=cf_dim, data_var=data_var, ds_area=ds_area, data_var_area=data_var_area)
        else:
            weights = "generate"
        # spatial average
        ds_o = xcdat_base.average_spatial(ds, data_var, cf_dim=cf_dim, weights=weights, **kwargs_average_spatial)
    log_debug(inspect__stack(), "output", data_var=data_var, ds=ds_o)
    return ds_o


def check_multidimensional_coordinates(ds: Union[array_wrapper, dataset_wrapper], **kwargs) -> bool:
    """
    Check if input object uses multidimensional coordinates (e.g., longitude[y, x]).

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: bool
        True is input object uses multidimensional coordinates
    """
    dim_lat = xarray_base.convert_cf_dim_key(ds, "Y")
    dim_lon = xarray_base.convert_cf_dim_key(ds, "X")
    bool_o = False
    if (basics.is_dim(dim_lat) is True and len(ds[dim_lat].shape) > 1) or (
            basics.is_dim(dim_lon) is True and len(ds[dim_lon].shape) > 1):
        bool_o = True
    return bool_o


def check_time_bounds(
        ds: Union[array_wrapper, dataset_wrapper],
        dim: Union[Hashable, str],
        time_bounds: tuple[str, str],
        side: Literal["lower", "upper"],
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Check ‘ds’ time bounds and compare them to ‘time_bounds’
    Sometimes selecting time is slightly wrong so here it is checked if one time step has not been included by error at
    the beginning or the end of the time series

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    :param dim: Hashable, str
        Name of the time dimension; e.g., dim = "time"
    :param time_bounds: tuple[str, str]
        Desired time bounds; time_bounds = ("1980-01-01", "2014-12-31")
    :param side: Literal["lower", "upper"]
        Check "lower" or "upper" bound
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        Object (as input) with correct time bounds
    """
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "dim": dim, "time_bounds": time_bounds, "side": side})
    # bound position (0 or 1)
    position = 0 if side == "lower" else 1
    # split time bound in ["year", "month", "day"]
    split_ta = basics.split_time_bound(str(xarray_base.get_time_bounds(ds)[position]))
    split_td = basics.split_time_bound(time_bounds[position])
    log_debug(inspect__stack(), "available vs. desired time_bounds: " + str(side), details={
        "available time_bound": split_ta, "desired time_bound": split_td})
    while (side == "lower" and all([float(k1) >= float(k2) for k1, k2 in zip(split_ta, split_td)]) is False) or \
            (side == "upper" and all([float(k1) <= float(k2) for k1, k2 in zip(split_ta, split_td)]) is False):
        # if side == "lower" & one of available("year", "month", "day") is smaller than desired ("year", "month", "day")
        # or
        # if side == "upper" & one of available("year", "month", "day") is larger than desired ("year", "month", "day")
        # remove the first (last) time step
        time_slice = slice(1, int(1e20)) if side == "lower" else slice(0, -1)
        ds = xarray_base.select_index(ds, {dim: time_slice})
        split_ta = basics.split_time_bound(str(xarray_base.get_time_bounds(ds)[position]))
        log_debug(inspect__stack(), "xarray_base.select_index", details={
            "available time_bound": split_ta, "desired time_bound": split_td})
    log_debug(inspect__stack(), "output", details={
        "available time_bound": split_ta, "desired time_bound": split_td})
    return ds


def compute_mask(
        input_array: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_mask: Union[str, None] = None,
        **kwargs) -> (Union[dataset_wrapper, None], Union[dict, None]):
    """
    Wraps create_land_sea_mask to obtain a dataset

    Input:
    ------
    :param input_array: dict[str, Union[array_wrapper, dataset_wrapper, dict]]
        Dictionary containing the data of the main variable for which a mask is needed
    :param data_var: str, optional
        Data variable in ‘input_array['array']’ which will be used to create the mask; e.g., data_var = "ts"
    :param data_var_mask: str, optional
        Name of the data_var for the newly created mask
    **kwargs - Discarded

    Output:
    -------
    :return: (dataset_wrapper, dict)
        Mask dataset and associated metadata
    """
    basics.log_info(inspect__stack(), "")
    details = {}
    l1 = ["input_array", "data_var"]
    l2 = [input_array, data_var]
    for k1, k2 in zip(l1, l2):
        details[str(k1) + ".type"] = str(type(k2))
        if isinstance(k2, dict) is True:
            details[str(k1) + ".keys"] = ", ".join(sorted(list(k2.keys()), key=lambda v: v.lower()))
        elif isinstance(k2, (float, int, str)) is True or k2 is None:
            details[k1] = str(k2)
    log_debug(inspect__stack(), "input", adjust=5, details=details)
    # generate mask based on input array
    ds_mask, metadata = None, None
    if isinstance(input_array, dict) is True and "array" in list(input_array.keys()):
        # get input array
        ds_array = input_array["array"]
        # create DataArray mask
        da_mask = create_land_sea_mask(ds_array)
        # DataArray to dataset
        ds_mask = xarray_base.to_dataset(da_mask, data_var_mask)
        # set variable attributes
        metadata = {
            "comment": "created using regionmask (https://regionmask.readthedocs.io/en/stable/defined_landmask.html)",
            "description": ";; a) Grid-Cell Land Fraction for Atmospheric Variables estimated",
            "long_name": "Grid-Cell Land Fraction for Atmospheric Variables",
            "short_name": "lsmask",
            "units": "1",
        }
        xarray_base.set_attributes_variable(ds_mask, data_var=data_var_mask, **metadata)
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
        xarray_base.set_attributes_global(ds_mask, **dict_t)
    return ds_mask, metadata


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
        **kwargs) -> (Union[dataset_wrapper, Union[dict, None]], dict):
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", details={"variable": variable, "list_variables": list_variables})
    attributes_global = ["grid_label", "source", "version"]
    attributes_variable = ["units"]
    arr: Union[array_wrapper, numpy__ndarray, None] = None
    dict_a = {}
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
        dict_t = {}
        for k1 in ["input", "variable"]:
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
        if "variable" in list(dict_t.keys()) and variable in list(dict_t["variable"].keys()):
            for k1 in ["description", "long_name", "short_name", "units"]:
                if k1 in list(dict_t["variable"][variable].keys()):
                    dict_t[k1] = deepcopy(dict_t["variable"][variable][k1])
                    if k1 == "long_name":
                        txt = str(dict_t["variable"][variable][k1]) + " read"
                        dict_t["description"] = basics.description_writer("", txt)
        # sort metadata dictionary
        metadata, _, _ = sort_dict(dict_t)
        # array to dataset
        ds = xarray_base.to_dataset(arr, variable)
        # set variable attributes
        dict_t, _, _ = combine_dict_levels(metadata)
        xarray_base.set_attributes_variable(ds, data_var=variable, **dict_t)
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
        # add bounds from first input variable
        ds_i = dict_ds[list_variables[0]]
        for nn in list(ds_i.keys()):
            if nn == list_variables[0]:
                continue
            ds = xarray_base.assign_to_dataset(ds, variables_kwargs={nn: ds_i[nn]})
    log_debug(inspect__stack(), "output", data_var=variable, ds=ds)
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
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
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
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds, details={
        "ds.type": type(ds), "dim": cf_dim, "ds_area.type": type(ds_area)})
    if isinstance(cf_dim, list) is False:
        cf_dim = ["X", "Y"]
    # by defaults weights are equal to 1
    weights = create_array(ds, data_var=data_var, data_var_o="weights", value=1)
    # temporal or spatial weights?
    if isinstance(cf_dim, list) is True and len(cf_dim) == 1 and cf_dim[0] == "T":
        # try to generate time weights using xcdat
        try:
            weights = xcdat_base.weights_temporal(ds, data_var)
        except Exception as err:
            log_debug(inspect__stack(), "WARNING cannot generate time weights using xcdat\n" + str(err))
    else:
        # is area available?
        da_area = xarray_base.copy(ds_area, data_var=data_var_area)
        if isinstance(da_area, array_wrapper) is False:
            # try to generate spatial weights using xcdat
            try:
                weights = xcdat_base.weights_spatial(ds, data_var, cf_dim=cf_dim)
            except Exception as err:
                log_debug(inspect__stack(), "WARNING cannot generate spatial weights using xcdat\n" + str(err))
        else:
            # get the name of the given axis (X, Y, Z)
            dim_lon = xarray_base.convert_cf_dim_key(ds_area, "X")
            dim_lat = xarray_base.convert_cf_dim_key(ds_area, "Y")
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
    Return a new DataArray of ‘value’ with the same shape, axes, coordinates, attributes,... as input DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.zeros_like.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
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
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds, details={"ds.type": type(ds)})
    # create an array of zeros with the same shape, axes, coordinates, attributes,... as input DataArray
    da_o = xarray_base.create_array_zero(ds, data_var)
    # name this new array
    da_o.name = data_var_o
    if value is None:
        da_o = xarray_base.where(da_o, da_o == 1)
    elif isinstance(value, (float, int)) is True and value != 0:
        da_o = xarray_base.where(da_o, da_o == value, other=value)
    log_debug(inspect__stack(), "input", data_var=data_var_o, ds=da_o, details={"da_o.type": type(da_o)})
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
    dim_lon = xarray_base.convert_cf_dim_key(ds, "X")
    dim_lat = xarray_base.convert_cf_dim_key(ds, "Y")
    # get latitude and longitude arrays
    dim_lon_array = xarray_base.get_dim_array(ds, dim_lon)
    dim_lat_array = xarray_base.get_dim_array(ds, dim_lat)
    # create a land-sea mask using regionmask
    land_sea_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110
    # adapt land-sea mask to input dataset or array -> land-sea mask is 0 on land and NaN on ocean
    land_sea_mask = land_sea_mask.mask(dim_lon_array, dim_lat_array)
    # convert zeros to 1 (1: land) and NaN to 0 (0: sea)
    land_sea_mask = xarray_base.where(land_sea_mask, land_sea_mask.isnull(), other=1)
    land_sea_mask = xarray_base.where(land_sea_mask, land_sea_mask == 1, other=0)
    if mask_as_boolean is True:
        # convert the land-sea mask to a boolean mask
        land_sea_mask = xarray_base.change_type(land_sea_mask, "bool")
    # drop attributes
    xarray_base.drop_given_attributes(land_sea_mask, xarray_base.get_attributes_keys(land_sea_mask))
    return land_sea_mask


def log_debug(
        stack,
        message: str,
        adjust: int = None,
        data_var: str = None,
        details: dict[str, Any] = None,
        ds: Union[array_wrapper, dataset_wrapper, None] = None):
    if isinstance(adjust, int) is False:
        adjust = 0
    tmp_dict = {}
    if isinstance(ds, dataset_wrapper) is True and isinstance(data_var, str) is True and \
            data_var in xarray_base.get_dataset_keys(ds):
        tmp_dict = {
            "variable": str(data_var), "dataset_keys": xarray_base.get_dataset_keys(ds),
            "dim_keys": xarray_base.get_dim_keys(ds), "shape": xarray_base.get_array_shape(ds, data_var),
            "min": min_global(ds, data_var), "max": max_global(ds, data_var)}
    elif isinstance(ds, array_wrapper) is True:
        tmp_dict = {
            "array_name": xarray_base.get_array_name(ds), "dim_keys": xarray_base.get_dim_keys(ds),
            "shape": xarray_base.get_array_shape(ds), "min": min_global(ds), "max": max_global(ds)}
    elif ds is not None:
        tmp_dict["array.type"] = type(ds)
        if isinstance(ds, str) is True:
            tmp_dict["ds given"] = str(ds)
    for k1, k2 in tmp_dict.items():
        message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    if isinstance(details, dict) is True:
        for k1, k2 in details.items():
            message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    basics.log_debug(stack, message)


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
    dim_name = xarray_base.convert_cf_dim_key(ds, "Y")
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
    dim_name = xarray_base.convert_cf_dim_key(ds, "X")
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
    dim_name = xarray_base.convert_cf_dim_key(ds, "T")
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
    dim_name = xarray_base.convert_cf_dim_key(ds, "Z")
    # get vertical array
    return xarray_base.get_dim_array(ds, dim_name)


def recreate_array(
        arr: numpy__ndarray,
        ds: Union[array_wrapper, dataset_wrapper],
        attrs_added: dict[str, str] = None,
        axis_added: Union[list[int], tuple[int], None] = None,
        coords_added: dict[Union[Hashable, str], Union[numpy__ndarray, array_wrapper]] = None,
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
        Array derived from ‘ds’ (e.g., a statistic was computed) that was transformed into a numpy.ndarray in the
        process
    :param ds: xarray.DataArray or xarray.Dataset
        Original xarray.DataArray or xarray.Dataset from which ‘arr‘ is derived
    :param attrs_added: dict[str, str] or None, optional
        Variable attributes to add to the DataArray; e.g, attrs_added = {"attr_name": "new attribute"}.
        If given, both ‘axis_added’ and ‘dim_added’ must be provided.
        Default is None (no dimension has been added)
    :param axis_added: list[int] or tuple[int] or None
        Position(s) of added dimension(s) (if any); e.g, axis_added = [0] or axis_added = [0, 1].
        If given, both ‘axis_added’ and ‘dim_added’ must be provided.
        Default is None (no dimension has been added)
    :param coords_added: dict[str, numpy.ndarray or xarray.DataArray] or None, optional
        Coordinates (tick labels) to use for indexing along each dimension.
        If ‘dim_added’ is not in ‘coords_added’, coordinates will be a sequence of numbers.
        Default is None
    :param data_var: str, optional
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
        Default is None
    :param data_var_o: str, optional
        Name of the output data variable.
        Default is ""
    :param dim_added: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been added from ‘ds’ to ‘arr‘ (if any);
        e.g, dim_added = ["x"] or dim_added = ["x", "y"].
        If given, both ‘axis_added’ and ‘dim_added’ must be provided.
        Default is None (no dimension has been added)
    :param dim_removed: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been removed from ‘ds’ to ‘arr‘ (e.g., to compute a statistic);
        e.g., dim_removed = ["x"] or dim_removed = ["x", "y"].
        Default is None (no dimension was removed)
    
    Output:
    -------
    :return: xarray.DataArray
        Input ‘arr’ wrapped in a xarray.DataArray.
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
    coordinates: dict[Union[Hashable, str], Union[numpy__ndarray, array_wrapper]] = dict(
        (k, xarray_base.get_dim_array(ds, k)) for k in dimensions)
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
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
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
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
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
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset, ‘data_var’ must be provided.
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


def roll_longitude(
        ds: Union[array_wrapper, dataset_wrapper],
        new_lon_min: Union[float, int, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Roll longitude dimension, either to ensure that longitude ranges from 0 to 360E or to start from new_lon_min (e.g.,
    new_lon_min = -70, longitude = [-70; 290])

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    :param new_lon_min: float or int or None, optional
        New minimum longitude
        Default is None
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object (as input) with rolled longitude
    """
    dim_lon = xarray_base.convert_cf_dim_key(ds, "X")
    if basics.is_dim(dim_lon) is True:
        # update longitude
        if isinstance(new_lon_min, (float, int)) is True:
            # add minimum value to dataset's longitude to shift the dimension
            # e.g., initial longitude = [0; 360], new_lon_min = -70, new longitude = [-70; 290]
            coords_kwargs = {dim_lon: ds[dim_lon] + new_lon_min}
        else:
            # ensure that longitude ranges from 0 to 360E
            coords_kwargs = {dim_lon: (360 + (ds[dim_lon] % 360)) % 360}
        ds = xarray_base.assign_coords(ds, coords_kwargs=coords_kwargs)
        # roll so that the first longitude of the dimension is the minimum longitude
        if check_multidimensional_coordinates(ds) is False:
            # normal roll method
            shifts = {dim_lon: -ds[dim_lon].argmin().values}
        else:
            # for multidimensional coordinates (e.g., curvilinear grids)
            # average lon along Y
            arr_lon = ds[dim_lon]
            lon_x = xarray_base.to_numpy(arr_lon).mean(axis=0)
            # find minimum value
            min_x = lon_x.argmin()
            # shift the last dimension of longitude coordinate
            last_lon_dim = xarray_base.get_dim_keys(arr_lon)[-1]
            shifts = {last_lon_dim: -min_x}
        ds = xarray_base.roll(ds, roll_coords=True, shifts=shifts)
    return ds


def select_depth(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: Union[str, None] = None,
        depth_bounds: Union[tuple[Union[float, int], Union[float, int]], None] = None,
        kwargs_sel: dict = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper, None]:
    """
    Select epoch based on given depth bounds.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str or None, optional
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset and ‘data_var’ in ‘ds’, log will have more details.
        Default is None
    :param depth_bounds: tuple[float or int, float or int], optional
        Depth extent to select data; e.g., depth_bounds = (0, 300).
        Default is None
    :param kwargs_sel: dict or None, optional
        kwargs to be passed to xarray_base.select
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with data of each array selected within given ‘depth_bounds‘.
    """
    kwargs_sel = set_instance(kwargs_sel, dict, False, {})
    log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "data_var": data_var, "depth_bounds": depth_bounds, **kwargs_sel})
    # fake loop to be able to break out
    ds_o = None
    for _ in [0]:
        if depth_bounds is None:
            ds_o = ds
            # no depth bounds to select data
            log_debug(inspect__stack(), "depth_bounds is None: NOT selected",
                      details={"depth_bounds": str(depth_bounds)})
            break
        elif isinstance(ds, (array_wrapper, dataset_wrapper)) is False:
            # input is neither a dataset nor a dataarray
            log_debug(inspect__stack(), "WARNING cannot select horizontal_bounds", details={
                "ds": str(type(ds)) + " should be xarray.DataArray, xarray.Dataset"})
            break
        # check given depth_bounds
        elif isinstance(depth_bounds, (list, tuple)) is False or (
                isinstance(depth_bounds, (list, tuple)) is True and len(depth_bounds) != 2) or (
                isinstance(depth_bounds, (list, tuple)) is True and len(depth_bounds) == 2 and
                all([isinstance(k, (float, int)) for k in depth_bounds]) is False):
            # depth bounds can be:
            #    - tuple[float or int, float or int]: depth within ‘depth_bounds’ will be selected; e.g., (0, 300)
            # given ‘depth_bounds’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot select depth_bounds", details={
                "depth_bounds": str(depth_bounds) + " should be tuple[float or int, float or int]"})
            break
        # test if depth dimension is available
        try:
            xarray_base.convert_cf_dim_key(ds, "Z")
        except Exception as err:
            log_debug(inspect__stack(), "WARNING cannot select depth_bounds\n" + str(err), details={
                "dim_name": "depth dimension not available", "dim_keys": xarray_base.get_dim_keys(ds)})
            break
        # get dimension name
        dim_name = xarray_base.convert_cf_dim_key(ds, "Z")
        log_debug(inspect__stack(), "xarray_base.convert_cf_dim_key", data_var=data_var,
                  details={"dim_name": dim_name}, ds=ds)
        # select using depth_bounds like (0, 300)
        ds_o = xarray_base.select(ds, {dim_name: slice(*depth_bounds)}, **kwargs_sel)
        log_debug(inspect__stack(), "xarray_base.select", data_var=data_var, ds=ds)
    log_debug(inspect__stack(), "output", data_var=data_var, ds=ds_o)
    return ds_o


def select_horizontal(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: Union[str, None] = None,
        horizontal_bounds: Union[dict[Literal["X", "Y"], tuple[Union[float, int]]], None] = None,
        mask_only: bool = False,
        kwargs_sel: dict = None,
        kwargs_where: dict = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper, None]:
    kwargs_sel = set_instance(kwargs_sel, dict, False, {})
    kwargs_where = set_instance(kwargs_where, dict, False, {})
    log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "data_var": data_var, "horizontal_bounds": horizontal_bounds, "mask_only": mask_only,
        "kwargs_sel": kwargs_sel, "kwargs_where": kwargs_where})
    # fake loop to be able to break out
    ds_o = None
    for _ in [0]:
        if horizontal_bounds is None:
            ds_o = ds
            # no horizontal bounds to select data
            log_debug(inspect__stack(), "horizontal_bounds is None: NOT selected",
                      details={"horizontal_bounds": str(horizontal_bounds)})
            break
        elif isinstance(ds, (array_wrapper, dataset_wrapper)) is False:
            # input is neither a dataset nor a dataarray
            log_debug(inspect__stack(), "WARNING cannot select horizontal_bounds", details={
                "ds": str(type(ds)) + " should be xarray.DataArray, xarray.Dataset"})
            break
        # get region
        lats = horizontal_bounds["Y"] if "Y" in list(horizontal_bounds.keys()) else None
        lons = horizontal_bounds["X"] if "X" in list(horizontal_bounds.keys()) else None
        if (isinstance(lats, (list, tuple)) is True and len(lats) == 1) or (
                isinstance(lons, (list, tuple)) is True and len(lons) == 1) or (
                isinstance(lats, (list, tuple)) is True and len(lats) >= 2 and
                isinstance(lons, (list, tuple)) is True and len(lats) != len(lons)):
            # region definition is not good
            details = {}
            if isinstance(lats, (list, tuple)) is True and len(lats) == 1:
                details["lat_error"] = "region latitude should have at least 2 values (" + str(len(lats)) + " given)"
            if isinstance(lons, (list, tuple)) is True and len(lons) == 1:
                details["lon_error"] = "region longitude should have at least 2 values (" + str(len(lons)) + " given)"
            if isinstance(lats, (list, tuple)) is True and len(lats) >= 2 and \
                    isinstance(lons, (list, tuple)) is True and len(lats) != len(lons):
                details["vertices_error"] = "region should have the same number of values for latitude and " + \
                                            "longitude (lats = " + str(len(lats)) + " & lons = " + str(len(lons)) + ")"
            details["horizontal_bounds"] = horizontal_bounds
            log_debug(inspect__stack(), "WARNING cannot select horizontal_bounds", details=details)
            break
        # get latitude and longitude
        dim_lat, dim_lon = None, None
        try:
            dim_lat = xarray_base.convert_cf_dim_key(ds, "Y")
        except Exception as err:
            log_debug(inspect__stack(), "WARNING input data doesn't have latitude\n" + str(err), details={
                "dim_name": "latitude dimension not available", "dim_keys": xarray_base.get_dim_keys(ds)})
        try:
            dim_lon = xarray_base.convert_cf_dim_key(ds, "X")
        except Exception as err:
            log_debug(inspect__stack(), "WARNING input data doesn't have longitude\n" + str(err), details={
                "dim_name": "longitude dimension not available", "dim_keys": xarray_base.get_dim_keys(ds)})
        arr_lat = ds[dim_lat] if basics.is_dim(dim_lat) is True else None
        arr_lon = ds[dim_lon] if basics.is_dim(dim_lon) is True else None
        if (isinstance(lats, (list, tuple)) is True and arr_lat is None) or (
                isinstance(lons, (list, tuple)) is True and arr_lon is None):
            # lat and/or lon required but corresponding dimension(s) is not available
            details = {}
            if isinstance(lats, (list, tuple)) is True and arr_lat is None:
                details["lat_error"] = "latitude must be selected but latitude is not available"
            if isinstance(lons, (list, tuple)) is True and arr_lon is None:
                details["lon_error"] = "longitude must be selected but latitude is not longitude"
            details["dim_keys"] = xarray_base.get_dim_keys(ds)
            details["horizontal_bounds"] = horizontal_bounds
            log_debug(inspect__stack(), "WARNING cannot select horizontal_bounds", details=details)
            break
        # -- mask outside region
        if (isinstance(lats, (list, tuple)) is True and len(lats) == 2) or (
                isinstance(lons, (list, tuple)) is True and len(lons) == 2) or (
                isinstance(lats, (list, tuple)) is True and len(lats) == 2 and
                isinstance(lons, (list, tuple)) is True and len(lons) == 2):
            # -- rectangular region
            # create condition
            if isinstance(lats, (list, tuple)) is True and isinstance(lons, (list, tuple)) is False:
                cond = (min(lats) <= arr_lat) & (arr_lat <= max(lats))
            elif isinstance(lats, (list, tuple)) is False and isinstance(lons, (list, tuple)) is True:
                cond = (min(lons) <= arr_lon) & (arr_lon <= max(lons))
            else:
                cond = (min(lats) <= arr_lat) & (arr_lat <= max(lats)) & (min(lons) <= arr_lon) & (arr_lon <= max(lons))
            # mask data outside region
            ds_o = xarray_base.where(ds, cond, **kwargs_where)
        else:
            # -- polygonal region
            # create region using regionmask
            region = numpy__array([[lo, la] for la, lo in zip(lats, lons)])
            region = regionmask.Regions([region])
            mask = region.mask(arr_lon, arr_lat)
            # mask data outside region
            ds_o = xarray_base.where(ds, xarray_base.notnull(mask), **kwargs_where)
        # -- select region
        if basics.is_dim(dim_lon) is True and mask_only is False and isinstance(lons, (list, tuple)) is True:
            # -- roll longitude
            lon_min, lon_max = min(lons), max(lons)
            # desired longitudes are usually defined [0; 360], but the input may not be, roll longitude if needed
            # new minimum value for the longitude will be the minimum longitude of the given region
            new_lon_min = deepcopy(lon_min)
            if lon_max - new_lon_min < 360:
                # modify minimum to be slightly lower, a maximum of 10 degree lower
                new_lon_min -= min(10., 360 - (lon_max - lon_min) / 2)
            # update longitude and roll, i.e., add minimum value to dataset's longitude to shift the dimension
            # e.g., initial longitude = [0; 360], desired = [-60; 30], new_lon_min = -70, new longitude = [-70; 290]
            ds_o = roll_longitude(ds_o, new_lon_min=new_lon_min)
            # -- select region (i.e., reduce the shape of the input data)
            # create indexers
            indexers = {}
            if arr_lat.shape == 1 and arr_lon.shape == 1:
                # regular grid
                if isinstance(lats, (list, tuple)) is True:
                    indexers[dim_lat] = slice(*(min(lats), max(lats)))
                if isinstance(lons, (list, tuple)) is True:
                    indexers[dim_lon] = slice(*(min(lons), max(lons)))
            else:
                # non-regular grid -> dataarray's lat/lon must be j/i or y/x or something like that
                # create condition
                if isinstance(lats, (list, tuple)) is True and isinstance(lons, (list, tuple)) is False:
                    cond = (min(lats) <= arr_lat) & (arr_lat <= max(lats))
                elif isinstance(lats, (list, tuple)) is False and isinstance(lons, (list, tuple)) is True:
                    cond = (min(lons) <= arr_lon) & (arr_lon <= max(lons))
                else:
                    cond = (min(lats) <= arr_lat) & (arr_lat <= max(lats)) & (min(lons) <= arr_lon) & \
                           (arr_lon <= max(lons))
                # find lat/lon dimensions in dataarray's
                dim_y, dim_x = xarray_base.get_dim_keys(arr_lon)
                arr_y, arr_x = ds[dim_y], ds[dim_x]
                # change lat/lon condition to dataarray's dimensions condition
                if isinstance(lats, (list, tuple)) is True:
                    reg_y = xarray_base.where(arr_y, cond)
                    reg_y = (min_global(reg_y), max_global(reg_y))
                    indexers[dim_y] = slice(*(min(reg_y), max(reg_y)))
                if isinstance(lons, (list, tuple)) is True:
                    reg_x = xarray_base.where(arr_x, cond)
                    reg_x = (min_global(reg_x), max_global(reg_x))
                    indexers[dim_x] = slice(*(min(reg_x), max(reg_x)))
            # select region
            ds_o = xarray_base.select(ds_o, indexers, **kwargs_sel)
    log_debug(inspect__stack(), "output", data_var=data_var, ds=ds_o)
    return ds_o


def select_time(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: Union[str, None] = None,
        time_bounds: Union[int, tuple[int], tuple[int, int], tuple[str, str], None] = None,
        kwargs_sel: dict = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper, None]:
    """
    Select epoch based on given time bounds.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str or None, optional
        Data variable in ‘ds’ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds’ is xarray.Dataset and ‘data_var’ in ‘ds’, log will have more details.
        Default is None
    :param time_bounds: int or tuple[int] or tuple[int, int] or tuple[str, str] or None, optional
        Time extent to select data; e.g., time_bounds = ("1980-01-01", "2014-12-31") or time_bounds = 120 or
        time_bounds = (0, 120).
        If ‘time_bounds’ is int or tuple[int], the first ‘time_bounds’ time steps will be selected.
        If ‘time_bounds’ is tuple[int, int], time steps between the first int and second int will be selected.
        Default is None
    :param kwargs_sel: dict or None, optional
        kwargs to be passed to xarray_base.select or xarray_base.select_index
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with data of each array selected within given ‘time_bounds‘.
    """
    kwargs_sel = set_instance(kwargs_sel, dict, False, {})
    log_debug(inspect__stack(), "input", details={
        "ds.type": type(ds), "data_var": data_var, "time_bounds": time_bounds, **kwargs_sel})
    # fake loop to be able to break out
    ds_o = None
    for _ in [0]:
        if time_bounds is None:
            ds_o = ds
            # no time bounds to select data
            log_debug(inspect__stack(), "time_bounds is None: NOT selected",
                      details={"time_bounds": str(time_bounds)})
            break
        elif isinstance(ds, (array_wrapper, dataset_wrapper)) is False:
            # input is neither a dataset nor a dataarray
            log_debug(inspect__stack(), "WARNING cannot select horizontal_bounds", details={
                "ds": str(type(ds)) + " should be xarray.DataArray, xarray.Dataset"})
            break
        # check given time_bounds
        if isinstance(time_bounds, int) is True:
            time_bounds = (time_bounds,)
        elif any([isinstance(time_bounds, (list, tuple)) is True and 0 < len(time_bounds) < 3 and
                 all([isinstance(k, int) for k in time_bounds]) is True,
                  isinstance(time_bounds, (list, tuple)) is True and len(time_bounds) == 2 and \
                  all([isinstance(k, str) for k in time_bounds]) is True]) is False:
            # time bounds can be:
            #    - int: the first ‘time_bounds’ time steps will be selected; e.g., 120
            #    - tuple[int]: the first ‘time_bounds’ time steps will be selected; e.g., (120,)
            #    - tuple[int, int]: time steps within ‘time_bounds’ will be selected; e.g., (120, 240)
            #    - tuple[str, str]: time steps within ‘time_bounds’ will be selected; e.g., ("1980-01-01", "2014-12-31")
            # given ‘time_bounds’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot select time_bounds", details={
                "time_bounds": str(time_bounds) + " should be int, tuple[str, str], tuple[int], tuple[int, int]"})
            break
        # test if time dimension is available
        try:
            xarray_base.convert_cf_dim_key(ds, "T")
        except Exception as err:
            log_debug(inspect__stack(), "WARNING cannot select time_bounds\n" + str(err), details={
                "dim_name": "time dimension not available", "dim_keys": xarray_base.get_dim_keys(ds)})
            break
        # get dimension name
        dim_name = xarray_base.convert_cf_dim_key(ds, "T")
        log_debug(inspect__stack(), "xarray_base.convert_cf_dim_key", data_var=data_var,
                  details={"dim_name": dim_name}, ds=ds)
        # check given time bounds type
        if isinstance(time_bounds, (list, tuple)) is True and all([isinstance(k, str) for k in time_bounds]) is True:
            # select using time_bounds like ("1980-01-01", "2014-12-31")
            ds_o = xarray_base.select(ds, {dim_name: slice(*time_bounds)}, **kwargs_sel)
            log_debug(inspect__stack(), "xarray_base.select", data_var=data_var,
                      details={"available time_bounds": xarray_base.get_time_bounds(ds)}, ds=ds)
            # sometimes selecting time is slightly wrong
            # this section checks if one time step has not been included by error at the beginning or the end of the
            # time series
            # check lower time bound
            ds_o = check_time_bounds(ds_o, dim_name, time_bounds, "lower")
            # check upper time bound
            ds_o = check_time_bounds(ds_o, dim_name, time_bounds, "upper")
            log_debug(inspect__stack(), "check_time_bounds", data_var=data_var,
                      details={"available time_bounds": xarray_base.get_time_bounds(ds_o)}, ds=ds_o)
        else:
            # select using time_bounds like (12, 24)
            ds_o = xarray_base.select_index(ds, {dim_name: slice(*time_bounds)}, **kwargs_sel)
            log_debug(inspect__stack(), "xarray_base.select_index", data_var=data_var, ds=ds_o)
    log_debug(inspect__stack(), "output", data_var=data_var, ds=ds_o)
    return ds_o


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
    log_debug(inspect__stack(), "input", details={"ds.type": type(ds), "cf_dim": cf_dim})
    # fake loop to be able to break out
    for _ in [0]:
        if isinstance(ds, (array_wrapper, dataset_wrapper)) is False:
            # wrong input
            log_debug(inspect__stack(), "WARNING cannot squeeze", details={"ds": str(type(ds)) + str(ds_error)})
            break
        # test is given dimension is available
        try:
            dim_name = xarray_base.convert_cf_dim_key(ds, cf_dim)
        except (Exception,):
            # no log: it is usual that the given dimension is not available and therefore squeeze is not applied
            break
        else:
            log_debug(inspect__stack(), "xarray_base.convert_cf_dim_key", details={
                "dim_name": dim_name, "dim_keys": xarray_base.get_dim_keys(ds)})
            # test if given dimension can be squeezed
            if xarray_base.get_array_size(ds, dim_name) < 2:
                # squeeze dimension
                ds = xarray_base.squeeze(ds, dim=dim_name, drop=True)
                log_debug(inspect__stack(), "xarray_base.squeeze", details={
                    "dim_keys": xarray_base.get_dim_keys(ds)})
                # remove bounds
                if isinstance(ds, dataset_wrapper) is True:
                    ds = xarray_base.drop_dataset_keys(
                        ds, [str(dim_name) + "_bounds", str(dim_name) + "_bnds"], errors="ignore")
                    log_debug(inspect__stack(), "xarray_base.drop_dataset_keys", details={
                        "dataset_keys": xarray_base.get_dataset_keys(ds)})
    log_debug(inspect__stack(), "output", details={"dim_keys": xarray_base.get_dim_keys(ds)})
    return ds
# ---------------------------------------------------------------------------------------------------------------------#
