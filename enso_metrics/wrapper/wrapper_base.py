# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools built over xarray and xcdat
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from dataclasses import dataclass
from glob import iglob as glob__iglob
from inspect import stack as inspect__stack
from math import ceil as math__ceil
from math import floor as math__floor
from os.path import isdir as os__path__isdir
from os.path import isfile as os__path__isfile
from typing import Any, Hashable, Literal, Union
# numpy
from numpy import array as numpy__array
from numpy import cos as numpy__cos
from numpy import deg2rad as numpy__deg2rad
from numpy import float64 as numpy__float64
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
def average_spatial(
        ds: Union[array_wrapper, dataset_wrapper],
        cf_dim: Union[list[Literal["X", "Y"]], None] = None,
        data_var: Union[Hashable, str, None] = None,
        data_var_area: Union[Hashable, str, None] = None,
        ds_area: Union[array_wrapper, dataset_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Compute spatial average (as defined by ‘cf_dim’) using xarray (xcdat for weights if needed)

    :param ds: xarray.DataArray or xarray.Dataset
    :param cf_dim: list[Literal["X", "Y"]]
        List of cf_dim along which to calculate averages; e.g., cf_dim = ["X", "Y"]
    :param data_var: Hashable or str
        Data variable in ‘ds’ on which to calculate averages; e.g., data_var = "ts"
    :param data_var_area: str, optional
        Data variable in ‘ds_area’
    :param ds_area: xarray.DataArray or xarray.Dataset, optional
        Area cell dataset
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        Object (as input) with the average of the data variable.
    """
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds,
                   details={"ds.type": type(ds), "dim": cf_dim, "ds_area.type": type(ds_area)})
    # fake loop to be able to break out
    ds_o = None
    for _ in [0]:
        if cf_dim is None or (isinstance(cf_dim, (list, tuple)) is True and len(cf_dim) == 0):
            ds_o = ds
            # no depth bounds to select data
            log_debug(inspect__stack(), "cf_dim is None: NO average", details={"cf_dim": str(cf_dim)})
            break
        elif error_ds(ds, inspect__stack(), message="cannot average spatially"):
            break
        # check given depth_bounds
        elif not(isinstance(cf_dim, (list, tuple)) is True and
                 3 > len(cf_dim) > 0 == len(list(set(cf_dim) - {"X", "Y"}))):
            # cf_dim can be (list of tuple):
            #    - list["X"]: array will be averaged along the longitude
            #    - list["Y"]: array will be averaged along the latitude
            #    - list["X", "Y"]: array will be averaged along the latitude & longitude
            # given ‘depth_bounds’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot average spatially", details={
                "cf_dim": str(cf_dim) + " should be list[Literal['X', 'Y']] or tuple[Literal['X', 'Y']]"})
            break
        # generate weights
        weights = compute_weights(ds, cf_dim=cf_dim, data_var=data_var, ds_area=ds_area, data_var_area=data_var_area)
        if weights is None:
            log_debug(inspect__stack(), "WARNING cannot average spatially: cannot generate weights")
            break
        # from cf dims to array dims
        dims = []
        for k in cf_dim:
            # find dim in array
            da_dim = get_dim_latitude_array(ds) if k == "Y" else get_dim_longitude_array(ds)
            if not isinstance(da_dim, array_wrapper):
                log_debug(inspect__stack(), "WARNING cannot average spatially: cannot find cf dim " + str(k))
                break
            # are lat/lon multidimensional coordinates?
            if len(xarray_base.get_array_shape(da_dim)) > 1:
                dim_keys = xarray_base.get_dim_keys(da_dim)
                dim = dim_keys[0] if k == "Y" else dim_keys[1]
            else:
                dim = get_dim_latitude(ds) if k == "Y" else get_dim_longitude(ds)
            dims.append(dim)
        # spatial average
        ds_o = xarray_base.mean(ds, data_var=data_var, dim=dims, keep_attrs=True, skipna=True, weights=weights)
        # weighted mean returns a dataarray -> recreate dataset
        if isinstance(ds, dataset_wrapper) and isinstance(ds_o, array_wrapper):
            # bounds are a pain, in the input ds they are fct of all coordinates (e.g., time_bnds[time, bnds, lat, lon])
            # so they cannot be easily taken from input ds and must be averaged the same way as data_var
            dict_bnds = {}
            for k1 in xarray_base.get_dataset_keys(ds):
                if k1 in [data_var] or ("_bnd" in k1 and str(k1).split("_bnd")[0] in dims) or \
                        ("_bound" in k1 and str(k1).split("_bound")[0] in dims) or \
                        ("bnd_" in k1 and str(k1).split("bnd_")[-1] in dims) or \
                        ("bnds_" in k1 and str(k1).split("bnds_")[-1] in dims) or \
                        ("bound_" in k1 and str(k1).split("bound_")[-1] in dims) or \
                        ("bounds_" in k1 and str(k1).split("bounds_")[-1] in dims) or \
                        ("vertice" in k1 and "lat" in k1 and "Y" in cf_dim) or \
                        ("vertice" in k1 and "lon" in k1 and "X" in cf_dim):
                    continue
                # select averaged dimensions available in given bounds key
                da_dims = xarray_base.get_dim_keys(ds, data_var=k1)
                tmp_dims = [k2 for k2 in dims if k2 in da_dims]
                if len(tmp_dims) == 0:
                    dict_bnds[k1] = xarray_base.to_array(ds, data_var=k1)
                    continue
                dict_bnds[k1] = xarray_base.mean(
                    ds, data_var=k1, dim=dims, keep_attrs=True, skipna=True, weights=weights)
            # coordinates are a pain for curvilinear grids (e.g., variable = ssh[time, j, i], coordinates = j[j], i[i],
            # latitude[j, i], longitude[j, i]) so they must be averaged the same way as data_var
            dict_coords = {}
            if check_multidimensional_coordinates(ds) is True and len(list({"X", "Y"} - set(cf_dim))) == 1:
                for k in ["lat", "lon"]:
                    if "Y" not in cf_dim and k == "lat":
                        dim_name = get_dim_latitude(ds)
                    elif "X" not in cf_dim and k == "lon":
                        dim_name = get_dim_longitude(ds)
                    else:
                        continue
                    dict_coords[dim_name] = xarray_base.mean(xarray_base.get_dim_array(ds, dim_name), dim=dims,
                                                             keep_attrs=True, skipna=True, weights=weights)
            # recreate dataset with data_var and associated bounds
            ds_o = recreate_dataset(ds, ds_o, data_var, cf_dims_removed=cf_dim, dict_coords=dict_coords,
                                    dict_da=dict_bnds)
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
    # get latitude and longitude arrays
    da_lat, da_lon = get_dim_latitude_array(ds), get_dim_longitude_array(ds)
    # check if latitude and longitude arrays and 2D arrays
    bool_o = False
    if (isinstance(da_lat, (array_wrapper, numpy__ndarray)) and len(xarray_base.get_array_shape(da_lat)) > 1) or (
            isinstance(da_lon, (array_wrapper, numpy__ndarray)) and len(xarray_base.get_array_shape(da_lon)) > 1):
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
        if isinstance(k2, dict):
            details[str(k1) + ".keys"] = ", ".join(sorted(list(k2.keys()), key=lambda v: v.lower()))
        elif isinstance(k2, (float, int, str)) is True or k2 is None:
            details[k1] = str(k2)
    log_debug(inspect__stack(), "input", adjust=5, details=details)
    # generate mask based on input array
    ds_mask, metadata = None, None
    if isinstance(input_array, dict) and "array" in list(input_array.keys()):
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
                    val = val[0] if len(val) == 1 else copy__deepcopy(dict_a[k2])
                    put_in_dict(dict_t, val, k1, variable, k2)
                put_in_dict(dict_t, dict_input["variable_computation"], k1, variable, "computation")
                val = str(list_variables[0]) if len(list_variables) == 1 else copy__deepcopy(list_variables)
                put_in_dict(dict_t, val, k1, variable, "out_name")
            else:
                # metadata 'variable'
                val = dict((k2, k3) for k2, k3 in variables_param.items() if k2 != "variable_type")
                put_in_dict(dict_t, val, k1, variable)
        if "variable" in list(dict_t.keys()) and variable in list(dict_t["variable"].keys()):
            for k1 in ["description", "long_name", "short_name", "units"]:
                if k1 in list(dict_t["variable"][variable].keys()):
                    dict_t[k1] = copy__deepcopy(dict_t["variable"][variable][k1])
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
        for nn in xarray_base.get_dataset_keys(ds_i):
            if nn == list_variables[0]:
                continue
            ds = xarray_base.assign_to_dataset(ds, variables_kwargs={nn: xarray_base.to_array(ds_i, nn)})
        # reindex (reverse) latitude dimension if needed
        dim_y = xarray_base.convert_cf_dim_key(ds, "Y")
        if basics.is_dim(dim_y):
            da_y = ds[dim_y]
            if not check_multidimensional_coordinates(ds) and da_y.values[0] > da_y.values[-1]:
                ds = xarray_base.reindex(ds, {dim_y: da_y[::-1]})
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
    # fake loop to be able to break out
    weights = None
    for _ in [0]:
        if error_ds(ds, inspect__stack(), message="cannot compute weights"):
            break
        elif cf_dim is None or (isinstance(cf_dim, (list, tuple)) is True and len(cf_dim) == 0):
            # no cf_dim given
            log_debug(inspect__stack(), "cf_dim is None: NO weights", details={"cf_dim": str(cf_dim)})
            break
        elif not (isinstance(cf_dim, (list, tuple)) is True and 0 < len(cf_dim) < 3 and (
                (len(cf_dim) == 1 and len(list(set(cf_dim) - {"T", "X", "Y"})) == 0) or
                (len(cf_dim) == 2 and len(list(set(cf_dim) - {"X", "Y"})) == 0))):
            # cf_dim can be (list of tuple):
            #    - list["T"]: weights will be computed to average along the time
            #    - list["X"]: weights will be computed to average along the longitude
            #    - list["Y"]: weights will be computed to average along the latitude
            #    - list["X", "Y"]: weights will be computed to average along latitude & longitude
            # given ‘depth_bounds’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot compute weights", details={
                "cf_dim": str(cf_dim) + " should be list[Literal['T', 'X', 'Y']] or tuple[Literal['T', 'X', 'Y']]"})
            break
        # get array without time
        da = xarray_base.to_array(ds, data_var)
        dim_time = get_dim_time(ds)
        if basics.is_dim(dim_time):
            da = xarray_base.select_index(da, indexers={dim_time: 0}, drop=True)
        # temporal or spatial weights?
        if isinstance(cf_dim, (list, tuple)) is True and len(cf_dim) == 1 and cf_dim[0] == "T":
            # try to generate time weights using xcdat
            try:
                weights = xcdat_base.weights_temporal(ds, data_var)
            except Exception as err:
                log_debug(inspect__stack(), "WARNING cannot generate time weights using xcdat\n" + str(err))
                weights = create_array(ds, data_var=data_var, data_var_o="weights", value=1)
        else:
            # is area available?
            weights = xarray_base.to_array(ds_area, data_var=data_var_area)
            if not isinstance(weights, array_wrapper) and not check_multidimensional_coordinates(ds):
                # try to generate spatial weights using xcdat
                try:
                    weights = xcdat_base.weights_spatial(ds, data_var, cf_dim=cf_dim)
                except Exception as err:
                    log_debug(inspect__stack(), "WARNING cannot generate spatial weights using xcdat\n" + str(err))
            if not isinstance(weights, array_wrapper) and "Y" in cf_dim:
                # try to get latitude dimension
                dim_lat = get_dim_latitude(ds)
                if not basics.is_dim(dim_lat):
                    break
                else:
                    # get latitude array
                    da_lat = xarray_base.get_dim_array(ds, dim_lat)
                    # compute weights using cos(latitude) only if 'Y' in cf_dim
                    weights = numpy__cos(numpy__deg2rad(da_lat))
            elif not isinstance(weights, array_wrapper):
                weights = create_array(da, data_var_o="weights", value=1)
            if isinstance(weights, array_wrapper):
                # set to 0 (i.e., weight of NaN cells is 0)
                weights = xarray_base.fill_nan(weights, 0)
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
        **kwargs) -> Union[array_wrapper, None]:
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
    # fake loop to be able to break out
    land_sea_mask = None
    for _ in [0]:
        if error_ds(ds, inspect__stack(), message="cannot create land-sea mask"):
            break
        # get latitude and longitude arrays
        da_lat, da_lon = get_dim_latitude_array(ds), get_dim_longitude_array(ds)
        if isinstance(da_lat, array_wrapper) is False or isinstance(da_lon, array_wrapper) is False:
            break
        # create a land-sea mask using regionmask
        land_sea_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110
        # adapt land-sea mask to input dataset or array -> land-sea mask is 0 on land and NaN on ocean
        land_sea_mask = land_sea_mask.mask(da_lon, da_lat)
        # convert zeros to 1 (1: land) and NaN to 0 (0: sea)
        land_sea_mask = xarray_base.where(land_sea_mask, land_sea_mask.isnull(), other=1)
        land_sea_mask = xarray_base.where(land_sea_mask, land_sea_mask == 1, other=0)
        if mask_as_boolean:
            # convert the land-sea mask to a boolean mask
            land_sea_mask = xarray_base.change_type(land_sea_mask, "bool")
        # drop attributes
        xarray_base.drop_given_attributes(land_sea_mask, xarray_base.get_attributes_keys(land_sea_mask))
    return land_sea_mask


def error_ds(ds: Any, stack, is_error: bool = True, message: str = "", **kwargs) -> bool:
    """
    Test is given ‘ds’ is xarray.DataArray or xarray.Dataset

    Input:
    ------
    :param ds: Any
    :param stack:
        Must be the result of the inspect.stack() function to log properly
        https://docs.python.org/3/library/inspect.html#inspect.stack
    :param is_error: bool
        If ‘ds’ is neither xarray.DataArray nor xarray.Dataset, True to log as 'error' else log as warning
        Default is True
    :param message: str
        Test for the error log
    **kwargs - Discarded

    Output:
    -------
    :return: bool
        True if ‘ds’ is neither xarray.DataArray nor xarray.Dataset, else False
    """
    error = False
    if not isinstance(ds, (array_wrapper, dataset_wrapper)):
        error = True
        a1 = "ERROR " if isinstance(is_error, bool) and is_error else "WARNING "
        log_debug(stack, str(a1) + str(message),
                  details={"ds": str(type(ds)) + " should be xarray.DataArray, xarray.Dataset"})
    return error


def log_debug(
        stack,
        message: str,
        adjust: int = None,
        data_var: str = None,
        details: dict[str, Any] = None,
        ds: Union[array_wrapper, dataset_wrapper, None] = None):
    if not isinstance(adjust, int):
        adjust = 0
    tmp_dict = {}
    if isinstance(ds, dataset_wrapper) is True and isinstance(data_var, str) is True and \
            data_var in xarray_base.get_dataset_keys(ds):
        tmp_dict = {
            "variable": str(data_var), "dataset_keys": xarray_base.get_dataset_keys(ds),
            "dim_keys": xarray_base.get_dim_keys(ds), "shape": xarray_base.get_array_shape(ds, data_var),
            "min": min_global(ds, data_var), "max": max_global(ds, data_var)}
    elif isinstance(ds, array_wrapper):
        tmp_dict = {
            "array_name": xarray_base.get_array_name(ds), "dim_keys": xarray_base.get_dim_keys(ds),
            "shape": xarray_base.get_array_shape(ds), "min": min_global(ds), "max": max_global(ds)}
    elif ds is not None:
        tmp_dict["array.type"] = type(ds)
        if isinstance(ds, str):
            tmp_dict["ds given"] = str(ds)
    for k1, k2 in tmp_dict.items():
        message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    if isinstance(details, dict):
        for k1, k2 in details.items():
            message += "\n" + str(k1).rjust(adjust) + ": " + str(k2)
    basics.log_debug(stack, message)


def get_dim_latitude(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> Union[Hashable, str, None]:
    """
    Return latitude dimension name of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str
        Latitude dimension name.
    """
    # get latitude as named in xarray.DataArray or xarray.Dataset
    dim_name = xarray_base.convert_cf_dim_key(ds, "Y")
    if not basics.is_dim(dim_name):
        log_debug(inspect__stack(), "WARNING cannot find latitude dimension")
    return dim_name


def get_dim_latitude_array(
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
    dim = get_dim_latitude(ds)
    # get latitude array
    da = None
    if basics.is_dim(dim):
        da = xarray_base.get_dim_array(ds, dim)
    if not isinstance(ds, (array_wrapper, numpy__ndarray)):
        log_debug(inspect__stack(), "WARNING cannot find latitude coordinates array")
    return da


def get_dim_longitude(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> Union[Hashable, str, None]:
    """
    Return longitude dimension name of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str
        Longitude dimension name.
    """
    # get longitude as named in xarray.DataArray or xarray.Dataset
    dim_name = xarray_base.convert_cf_dim_key(ds, "X")
    if not basics.is_dim(dim_name):
        log_debug(inspect__stack(), "WARNING cannot find longitude dimension")
    return dim_name


def get_dim_longitude_array(
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
    dim = get_dim_longitude(ds)
    # get longitude array
    da = None
    if basics.is_dim(dim):
        da = xarray_base.get_dim_array(ds, dim)
    if not isinstance(ds, (array_wrapper, numpy__ndarray)):
        log_debug(inspect__stack(), "WARNING cannot find longitude coordinates array")
    return da


def get_dim_time(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> Union[Hashable, str, None]:
    """
    Return time dimension name of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str
        Time dimension name.
    """
    # get time as named in xarray.DataArray or xarray.Dataset
    dim_name = xarray_base.convert_cf_dim_key(ds, "T")
    if not basics.is_dim(dim_name):
        log_debug(inspect__stack(), "WARNING cannot find time dimension")
    return dim_name


def get_dim_time_array(
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
    dim = get_dim_time(ds)
    # get time array
    da = None
    if basics.is_dim(dim):
        da = xarray_base.get_dim_array(ds, dim)
    if not isinstance(ds, (array_wrapper, numpy__ndarray)):
        log_debug(inspect__stack(), "WARNING cannot find time coordinates array")
    return da


def get_dim_vertical(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> Union[Hashable, str, None]:
    """
    Return vertical (depth, height, pressure) dimension name of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str
        Vertical dimension name.
    """
    # get vertical as named in xarray.DataArray or xarray.Dataset
    dim_name = xarray_base.convert_cf_dim_key(ds, "X")
    if not basics.is_dim(dim_name):
        log_debug(inspect__stack(), "WARNING cannot find vertical dimension")
    return dim_name


def get_dim_vertical_array(
        ds: Union[array_wrapper, dataset_wrapper],
        **kwargs) -> array_wrapper:
    """
    Return vertical (depth, height, pressure) array of given xarray.DataArray or xarray.Dataset.

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
    dim = get_dim_vertical(ds)
    # get vertical array
    da = None
    if basics.is_dim(dim):
        da = xarray_base.get_dim_array(ds, dim)
    if not isinstance(ds, (array_wrapper, numpy__ndarray)):
        log_debug(inspect__stack(), "WARNING cannot find vertical coordinates array")
    return da


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
        Original xarray.DataArray or xarray.Dataset from which ‘arr’ is derived
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
        Dimension name(s) that has been added from ‘ds’ to ‘arr’ (if any);
        e.g, dim_added = ["x"] or dim_added = ["x", "y"].
        If given, both ‘axis_added’ and ‘dim_added’ must be provided.
        Default is None (no dimension has been added)
    :param dim_removed: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been removed from ‘ds’ to ‘arr’ (e.g., to compute a statistic);
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
            dimensions.remove(k)
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


def open_dataset(
        filename: Union[str, list[str]],
        data_var: str = None,
        package: Literal["xarray", "xcdat"] = "xcdat",
        kwargs_open_dataset: dict = None,
        **kwargs) -> Union[dataset_wrapper, None]:
    kwargs_open_dataset = set_instance(kwargs_open_dataset, dict, False, {})
    ds = None
    try:
        if package == "xcdat":
            ds = xcdat_base.open_dataset(filename, data_var=data_var, **kwargs_open_dataset)
        else:
            ds = xarray_base.open_dataset(filename, **kwargs_open_dataset)
        # ds = dict_packages[package](filename, data_var=data_var, **kwargs_open_dataset)
    except Exception as err:
        message = "can't read (" + str(data_var) + ")" + "\n" + str(err)
        basics.log_error(inspect__stack(), message)
        # WARNING: cannot read variable or file
        if isinstance(filename, str):
            filename = [filename]
        paths = ["/".join(k.split("/")[:-1]) for k in filename]
        paths_test = [os__path__isdir(k) for k in paths]
        files = []
        for k in filename:
            files += list(glob__iglob(k))
        files = sorted(list(set(files)), key=lambda s: s.lower())
        files_test = [os__path__isfile(k) for k in files]
        files_string = ""
        if len(files) > 0:
            for k in files:
                files_string += "\n" + str().ljust(5) + str(k)
        else:
            files_string = "no file matches this file pattern"
        details = {
            "directory": str(paths),
            "isdir": str(paths_test),
            "file": str(filename),
            "isfile": str(files_test),
            "list": str(files_string)}
        log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
    return ds


def recreate_dataset(
        ds: dataset_wrapper,
        da: array_wrapper,
        data_var: Union[Hashable, str],
        cf_dims_removed: list[Literal["T", "X", "Y", "Z"]] = None,
        dict_coords: dict[Union[Hashable, str], Union[array_wrapper, numpy__ndarray]] = None,
        dict_da: dict[Union[Hashable, str], array_wrapper] = None,
        **kwargs) -> dataset_wrapper:
    """
    Use input ‘ds’ (which is supposed to be the original dataset) to recreate a dataset for ‘da’ (which is supposed to
    be a dataarray from ‘ds’ on which an operation has been performed).

    Input:
    ------
    :param ds: xarray.Dataset
    :param da: xarray.DataArray
    :param data_var: Hashable or str
        Name of the variable in the dataset
    :param cf_dims_removed: list[Literal["T", "X", "Y", "Z"]], optional
        List of cf_dims removed
    :param dict_coords: dict[Union[Hashable, str], Union[array_wrapper, numpy__ndarray]], optional
        Array to add as coordinates
    :param dict_da: dict[Union[Hashable, str], array_wrapper], optional
        Other xarray.DataArray to add (usually bounds)
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.Dataset
    """
    cf_dims_removed = set_instance(cf_dims_removed, list, False, [])
    dict_coords = set_instance(dict_coords, dict, False, {})
    dict_da = set_instance(dict_da, dict, False, {})
    # input DataArray to Dataset
    ds_o = xarray_base.to_dataset(da, data_var)
    # add other DataArray to new Dataset
    if len(list(dict_da.keys())) > 0:
        ds_o = xarray_base.assign_to_dataset(ds_o, variables_kwargs=dict_da)
    # add coordinates
    if len(list(dict_coords.keys())) > 0:
        ds_o = xarray_base.assign_coords(ds_o, coords_kwargs=dict_coords)
    # get input dimensions
    dims = xarray_base.get_dim_keys(da)
    # get data_vars and bounds
    dict_var = {}
    for k in xarray_base.get_dataset_keys(ds):
        if k in [data_var] + list(dict_da.keys()) or ("_bnd" in k and str(k).split("_bnd")[0] not in dims) or \
                ("_bound" in k and str(k).split("_bound")[0] not in dims) or \
                    ("vertice" in k and "lat" in k and "Y" in cf_dims_removed) or \
                    ("vertice" in k and "lon" in k and "X" in cf_dims_removed):
            continue
        dict_var[k] = xarray_base.to_array(ds, k)
    # add them to output dataset
    if len(list(dict_var.keys())) > 0:
        ds_o = xarray_base.assign_to_dataset(ds_o, variables_kwargs=dict_var)
    # set global attributes from input dataset to output dataset
    xarray_base.set_attributes_global(ds_o, **xarray_base.get_attributes(ds))
    return ds_o


def redo_bounds(ds: dataset_wrapper, cf_dims: list[Literal["T", "X", "Y", "Z"]]) -> dataset_wrapper:
    # operations on ds changes bounds and xcdat doesn't like that: bounds must be deleted and recreated
    for cf_dim in cf_dims:
        if check_multidimensional_coordinates(ds) and cf_dim in ["X", "Y"]:
            continue
        # get dimension key
        dim = xarray_base.convert_cf_dim_key(ds, cf_dim)
        if not isinstance(dim, str):
            continue
        # get bounds key
        try:
            dim_bnds = xarray_base.get_attribute(ds, "bounds", data_var=dim)
        except (Exception,):
            # no log: sometimes bounds are not defined correctly
            continue
        if not isinstance(dim_bnds, str):
            continue
        # delete current time bounds
        if dim_bnds in list(ds.keys()):
            ds = xarray_base.drop_dataset_keys(ds, [dim_bnds])
        # use xcdat to set bounds
        ds = xcdat_base.set_auto_bounds(ds, [cf_dim])
        try:
            dim_bnds = xarray_base.get_attribute(ds, "bounds", data_var=dim)
        except (Exception,):
            # no log: sometimes bounds are not defined correctly
            continue
    return ds


def regrid_horizontal(
        ds: dataset_wrapper,
        data_var: Union[Hashable, str, None] = None,
        grid: Union[array_wrapper, dataset_wrapper, str, None] = None,
        method: Literal[
            "bilinear", "conservative", "conservative_normed", "patch", "nearest_s2d", "nearest_d2s"] = "conservative",
        tool: Literal["regrid2", "xesmf"] = "regrid2",
        kwargs_regridder_horizontal: dict = None,
        **kwargs) -> dataset_wrapper:
    kwargs_regridder_horizontal = set_instance(kwargs_regridder_horizontal, dict, False, {})
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds,
              details={"ds.type": type(ds), "grid.type": type(grid)})
    ds_o = None
    # fake loop to be able to break out when an error occurs
    for _ in [0]:
        if isinstance(grid, str) and "gaussian" in grid and "x" in grid:
            # grid should be like 'gaussian_1x1', the number of latitudes is computed using the first integer
            nlat = int(round(180. / int(grid.replace("gaussian_", "").split("x")[0]), 0))
            output_grid = xcdat_base.create_gaussian_grid(nlat)
        elif isinstance(grid, str) and "uniform" in grid and "x" in grid:
            # grid should be like 'gaussian_1x1', the number of latitudes is computed using the first integer
            lat_start, lat_stop, lon_start, lon_stop = -89.5, 89.5, 0.5, 359.5
            lat_delta = float(grid.replace("uniform_", "").split("x")[0])
            lon_delta = float(grid.replace("uniform_", "").split("x")[1])
            # # adapt lat_start, lat_stop based on input data
            # try:
            #     bounds_lat = xcdat_base.get_bounds(ds, "Y", data_var=data_var)
            # except (Exception,):
            #     # no log: sometimes bounds are not defined correctly
            #     # use latitudes to define lat_start and lat_stop
            #     da_lat = get_dim_latitude_array(ds)
            #     min_max = min_max_global(da_lat)
            #     arr_lat = xarray_base.to_numpy(da_lat)
            #     if check_multidimensional_coordinates(ds):
            #         # for multidimensional coordinates (e.g., curvilinear grids) average lat along X
            #         arr_lat = arr_lat.mean(axis=1)
            #     dy1, dy0 = (arr_lat[-1] - arr_lat[-2]) / 2, (arr_lat[1] - arr_lat[0]) / 2
            #     if min_max[1] + dy1 - min_max[0] - dy0 < 85:
            #         # not the entire globe is available: find new lat_start and lat_stop
            #         lat_start = max(-90, math__floor(min_max[0] - dy0))
            #         lat_stop = min(90, math__ceil(min_max[1] + dy1))
            # else:
            #     # use latitude bounds to define lat_start and lat_stop
            #     min_max = min_max_global(bounds_lat)
            #     if min_max[1] - min_max[0] < 85:
            #         # not the entire globe is available: find new lat_start and lat_stop
            #         lat_start = max(-90, math__floor(min_max[0]))
            #         lat_stop = min(90, math__ceil(min_max[1]))
            # # adapt lon_start, lon_stop based on input data
            # try:
            #     bounds_lon = xcdat_base.get_bounds(ds, "X", data_var=data_var)
            # except (Exception,):
            #     # no log: sometimes bounds are not defined correctly
            #     # use longitudes to define lon_start and lon_stop
            #     da_lon = get_dim_longitude_array(ds)
            #     min_max = min_max_global(da_lon)
            #     arr_lon = xarray_base.to_numpy(da_lon)
            #     if check_multidimensional_coordinates(ds):
            #         # for multidimensional coordinates (e.g., curvilinear grids) average lon along Y
            #         arr_lon = arr_lon.mean(axis=0)
            #     dx1, dx0 = (arr_lon[-1] - arr_lon[-2]) / 2, (arr_lon[1] - arr_lon[0]) / 2
            #     if min_max[1] + dx1 - min_max[0] - dx0 < 355:
            #         # not the entire globe is available: find new lon_start and lon_stop
            #         if (-360 <= min_max[0] <= 0 and -360 <= min_max[1] <= 0) or (
            #                 0 <= min_max[0] <= 360 and 0 <= min_max[1] <= 360):
            #             # lon_start and lon_stop can be defined between 0 and 360
            #             lon_start = max(0, math__floor(min_max[0] - dx0) % 360)
            #             lon_stop = min(360, math__ceil(min_max[1] + dx1) % 360)
            #         else:
            #             # lon_start and lon_stop cannot be defined between 0 and 360
            #             lon_start, lon_stop = math__floor(min_max[0] - dx0), math__ceil(min_max[1] + dx1)
            # else:
            #     # use longitude bounds to define lon_start and lon_stop
            #     min_max = min_max_global(bounds_lon)
            #     if min_max[1] - min_max[0] < 355:
            #         # not the entire globe is available: find new lon_start and lon_stop
            #         if (-360 <= min_max[0] <= 0 and -360 <= min_max[1] <= 0) or (
            #                 0 <= min_max[0] <= 360 and 0 <= min_max[1] <= 360):
            #             # lon_start and lon_stop can be defined between 0 and 360
            #             lon_start = max(0, math__floor(min_max[0]) % 360)
            #             lon_stop = min(360, math__ceil(min_max[1]) % 360)
            #         else:
            #             # lon_start and lon_stop cannot be defined between 0 and 360
            #             lon_start, lon_stop = math__floor(min_max[0]), math__ceil(min_max[1])
            # generate uniform grid
            output_grid = xcdat_base.create_uniform_grid(lat_start, lat_stop, lat_delta, lon_start, lon_stop, lon_delta)
        elif isinstance(grid, (array_wrapper, dataset_wrapper)):
            output_grid = grid
        else:
            # given ‘grid’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot regrid horizontally", details={
                "grid": str(grid) + " should be string like 'gaussian_latxlon' or 'uniform_latxlon' or array"})
            break
        # check keywords
        if check_multidimensional_coordinates(ds):
            tool: Literal["xesmf"] = "xesmf"
            method: Literal["bilinear"] = "bilinear"
        # regrid
        ds_o = xcdat_base.regridder_horizontal(ds, data_var, output_grid, tool=tool, method=method,
                                               **kwargs_regridder_horizontal)
    return ds_o


def regrid_vertical(
        ds: dataset_wrapper,
        data_var: Union[Hashable, str, None] = None,
        grid: Union[array_wrapper, dataset_wrapper, int, None] = None,
        tool: Literal["xgcm"] = "xgcm",
        kwargs_regridder_vertical: dict = None,
        **kwargs) -> dataset_wrapper:
    kwargs_regridder_vertical = set_instance(kwargs_regridder_vertical, dict, False, {})
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds,
              details={"ds.type": type(ds), "grid.type": type(grid)})
    ds_o = None
    # 50 and 60 levels grids are the same between 0 and 160m depth
    levels = {
        50: [5.00000000e+00, 1.50000000e+01, 2.50000000e+01, 3.50000000e+01, 4.50000000e+01, 5.50000000e+01,
             6.50000000e+01, 7.50000000e+01, 8.50000000e+01, 9.50000000e+01, 1.05000000e+02, 1.15000000e+02,
             1.25000000e+02, 1.35000000e+02, 1.45000000e+02, 1.55000000e+02, 1.65000000e+02, 1.75000000e+02,
             1.85000000e+02, 1.95000000e+02, 2.05000000e+02, 2.16846756e+02, 2.41349014e+02, 2.80780731e+02,
             3.43250458e+02, 4.27315552e+02, 5.36715637e+02, 6.65414124e+02, 8.12781616e+02, 9.69065125e+02,
             1.13093494e+03, 1.28960461e+03, 1.45577014e+03, 1.62292566e+03, 1.80155811e+03, 1.98485461e+03,
             2.18290479e+03, 2.38841748e+03, 2.61093506e+03, 2.84256445e+03, 3.09220483e+03, 3.35129468e+03,
             3.62805762e+03, 3.91326440e+03, 4.21449512e+03, 4.52191797e+03, 4.84256592e+03, 5.16612988e+03,
             5.49924512e+03, 5.83129443e+03],
        60: [5.00000000e+00, 1.50000000e+01, 2.50000000e+01, 3.50000000e+01, 4.50000000e+01, 5.50000000e+01,
             6.50000000e+01, 7.50000000e+01, 8.50000000e+01, 9.50000000e+01, 1.05000000e+02, 1.15000000e+02,
             1.25000000e+02, 1.35000000e+02, 1.45000000e+02, 1.55000000e+02, 1.65098398e+02, 1.75479043e+02,
             1.86291270e+02, 1.97660273e+02, 2.09711387e+02, 2.22578281e+02, 2.36408828e+02, 2.51370156e+02,
             2.67654199e+02, 2.85483652e+02, 3.05119219e+02, 3.26867988e+02, 3.51093477e+02, 3.78227617e+02,
             4.08784648e+02, 4.43377695e+02, 4.82736719e+02, 5.27728008e+02, 5.79372891e+02, 6.38862617e+02,
             7.07563281e+02, 7.87002500e+02, 8.78825234e+02, 9.84705859e+02, 1.10620422e+03, 1.24456688e+03,
             1.40049719e+03, 1.57394641e+03, 1.76400328e+03, 1.96894422e+03, 2.18645656e+03, 2.41397156e+03,
             2.64900125e+03, 2.88938469e+03, 3.13340469e+03, 3.37979344e+03, 3.62767031e+03, 3.87645188e+03,
             4.12576812e+03, 4.37539250e+03, 4.62519031e+03, 4.87508344e+03, 5.12502812e+03, 5.37500000e+03],
        75: [5.05760017e-01, 1.55585530e+00, 2.66768175e+00, 3.85627974e+00, 5.14036125e+00, 6.54303362e+00,
             8.09251839e+00, 9.82275043e+00, 1.17736795e+01, 1.39910380e+01, 1.65253215e+01, 1.94298028e+01,
             2.27576162e+01, 2.65583009e+01, 3.08745618e+01, 3.57402047e+01, 4.11800247e+01, 4.72118941e+01,
             5.38506372e+01, 6.11128402e+01, 6.90216839e+01, 7.76111618e+01, 8.69294254e+01, 9.70413126e+01,
             1.08030281e+02, 1.20000001e+02, 1.33075822e+02, 1.47406245e+02, 1.63164456e+02, 1.80549922e+02,
             1.99789960e+02, 2.21141180e+02, 2.44890622e+02, 2.71356387e+02, 3.00887515e+02, 3.33862834e+02,
             3.70688484e+02, 4.11793845e+02, 4.57625617e+02, 5.08639904e+02, 5.65292274e+02, 6.28025970e+02,
             6.97258648e+02, 7.73368259e+02, 8.56678942e+02, 9.47447897e+02, 1.04585430e+03, 1.15199125e+03,
             1.26586142e+03, 1.38737698e+03, 1.51636363e+03, 1.65256845e+03, 1.79567082e+03, 1.94529547e+03,
             2.10102652e+03, 2.26242161e+03, 2.42902521e+03, 2.60038049e+03, 2.77603935e+03, 2.95557038e+03,
             3.13856486e+03, 3.32464083e+03, 3.51344558e+03, 3.70465666e+03, 3.89798194e+03, 4.09315874e+03,
             4.28995243e+03, 4.48815461e+03, 4.68758110e+03, 4.88806979e+03, 5.08947856e+03, 5.29168316e+03,
             5.49457529e+03, 5.69806076e+03, 5.90205781e+03],
    }
    # fake loop to be able to break out when an error occurs
    for _ in [0]:
        if isinstance(grid, int) and grid in list(levels.keys()):
            output_grid = xcdat_base.create_grid(
                z=xcdat_base.create_axis("depth", numpy__array(levels[grid], dtype=numpy__float64)),
                attrs={"positive": "down", "units": "m"})
        elif isinstance(grid, (array_wrapper, dataset_wrapper)):
            output_grid = grid
        else:
            # given ‘grid’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot regrid vertically", details={
                "grid": str(grid) + " should be integer " + str(list(levels.keys())) + " or array"})
            break
        # regrid
        ds_o = xcdat_base.regridder_vertical(ds, data_var, output_grid, tool=tool, **kwargs_regridder_vertical)
    return ds_o


def remove_fit(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: Union[Hashable, str, None] = None,
        deg: int = 1,
        dim: Union[Hashable, str] = "T",
        kwargs_polyfit: dict = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    kwargs_polyfit = set_instance(kwargs_polyfit, dict, False, {})
    basics.log_info(inspect__stack(), "")
    log_debug(inspect__stack(), "input", data_var=data_var, ds=ds,
              details={"ds.type": type(ds), "dim": dim, "deg": deg})
    ds_o = None
    # fake loop to be able to break out when an error occurs
    for _ in [0]:
        # list dimensions
        list_dimensions = xarray_base.get_dim_keys(ds)
        # ensures that dim is a dimension in input ds
        if dim in ["T", "X", "Y", "Z"]:
            dim_name = xarray_base.convert_cf_dim_key(ds, dim)
        elif dim in list_dimensions:
            dim_name = copy__deepcopy(dim)
        else:
            # given ‘dim’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot remove fit", details={
                "dim": str(dim) + " should be cf_dim or dim in " + str(list_dimensions)})
            break
        # compute coefficient
        p = xarray_base.polyfit(ds, dim_name, deg, **kwargs_polyfit)
        # remove fit
        da_dim = xarray_base.get_dim_array(ds, dim_name)
        ds_o = xarray_base.copy(ds, deep=True) - xarray_base.polyval(da_dim, p[str(data_var) + "_polyfit_coefficients"])
    return ds_o


def rename_variable(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var_i: Union[Hashable, str],
        data_var_o: Union[Hashable, str],
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Rename data_var_i to data_var_o in input object.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    :param data_var_i: Hashable or str
        Data variable ro rename
    :param data_var_o: Hashable or str
        New name of data variable
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        Object (as input) with renamed variable
    """
    # get array
    da = xarray_base.to_array(ds, data_var_i)
    # rename variable
    if isinstance(da, array_wrapper):
        da = xarray_base.rename(da, data_var_o)
    if isinstance(da, dataset_wrapper):
        da = xarray_base.rename(da, {data_var_i: data_var_o})
    if isinstance(ds, dataset_wrapper):
        ds = xarray_base.set_array_in_place(ds, da, data_var_i)
        ds = xarray_base.rename_vars(ds, {data_var_i: data_var_o})
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
    dim_lon = get_dim_longitude(ds)
    if basics.is_dim(dim_lon):
        # --- Step 1: Update longitude
        # get longitude
        arr_lon = xarray_base.get_dim_array(ds, dim_lon)
        if isinstance(new_lon_min, (float, int)):
            # add minimum value to dataset's longitude to shift the dimension
            # e.g., initial longitude = [0; 359], new_lon_min = -70, new longitude = [-70; 289]
            if new_lon_min >= 0:
                arr_lon = xarray_base.where(arr_lon, arr_lon >= new_lon_min, other=arr_lon + 360)
            else:
                arr_lon = xarray_base.where(arr_lon, arr_lon < 360 - new_lon_min, other=arr_lon - 360)
            coords_kwargs = {dim_lon: arr_lon}
        else:
            # ensure that longitude ranges from 0 to 360E
            coords_kwargs = {dim_lon: (360 + (arr_lon % 360)) % 360}
        # update longitude
        ds = xarray_base.assign_coords(ds, coords_kwargs=coords_kwargs)
        # --- Step 2: Roll so that the first longitude of the dimension is the minimum longitude
        if not check_multidimensional_coordinates(ds):
            # normal roll method
            shifts = {dim_lon: -int(xarray_base.to_numpy(arr_lon).argmin())}
        else:
            # for multidimensional coordinates (e.g., curvilinear grids)
            # average lon along Y
            lon_x = xarray_base.to_numpy(arr_lon).mean(axis=0)
            # find minimum value
            min_x = lon_x.argmin()
            # shift the last dimension of longitude coordinate
            last_lon_dim = xarray_base.get_dim_keys(arr_lon)[-1]
            shifts = {last_lon_dim: -min_x}
        # roll
        ds = xarray_base.roll(ds, roll_coords=True, shifts=shifts)
        # --- Step 3: Update longitude bounds (if applicable)
        arr_bnds, dim_bnds = None, None
        # find bounds
        attrs = xarray_base.get_attributes(arr_lon)
        if "bounds" in list(attrs.keys()):
            # bounds named in attributes
            dim_bnds = attrs["bounds"]
            # check if it exists
            if dim_bnds in list(ds.keys()):
                arr_bnds = ds[dim_bnds]
        if isinstance(arr_bnds, array_wrapper):
            # check if any variable looks like relevant bounds
            for k in list(ds.keys()):
                if str(dim_lon) + "_bnd" in k or str(dim_lon) + "_bound" in k:
                    arr_bnds = ds[k]
                    dim_bnds = copy__deepcopy(k)
                    break
        # update longitude bounds
        if isinstance(arr_bnds, array_wrapper):
            # remove mean bounds and add longitude
            # e.g., initial longitude = [0; 359], new_lon_min = -70, new longitude = [-70; 289]
            # rolled bounds are still bnds = [[289.5, 290.5] ... [288.5, 289.5]]
            # remove mean: bnds = [[-0.5, 0.5] ... [-0.5, 0.5]]
            # add longitude: bnds = [[-70.5, -69.5] ... [288.5, 289.5]]
            arr_bnds = arr_bnds - xarray_base.mean(arr_bnds, dim="B") + ds[dim_lon]
            ds = xarray_base.set_array_in_place(ds, arr_bnds, data_var=dim_bnds)
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
        New object with data of each array selected within given ‘depth_bounds’.
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
        elif error_ds(ds, inspect__stack(), message="cannot select depth_bounds"):
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
        # get dimension name
        dim_vertical = get_dim_vertical(ds)
        log_debug(inspect__stack(), "get_dim_vertical", data_var=data_var,
                  details={"dim_vertical": dim_vertical}, ds=ds)
        if not basics.is_dim(dim_vertical):
            break
        # select using depth_bounds like (0, 300)
        ds_o = xarray_base.select(ds, {dim_vertical: slice(*depth_bounds)}, **kwargs_sel)
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
        elif error_ds(ds, inspect__stack(), message="cannot select horizontal_bounds"):
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
        dim_lat, dim_lon = get_dim_latitude(ds), get_dim_longitude(ds)
        da_lat, da_lon = get_dim_latitude_array(ds), get_dim_longitude_array(ds)
        if (isinstance(lats, (list, tuple)) is True and isinstance(da_lat, array_wrapper) is False) or (
                isinstance(lons, (list, tuple)) is True and isinstance(da_lon, array_wrapper) is False):
            # lat and/or lon required but corresponding dimension(s) is not available
            details = {}
            if isinstance(lats, (list, tuple)) is True and isinstance(da_lat, array_wrapper) is False:
                details["lat_error"] = "latitude must be selected but latitude is not available"
            if isinstance(lons, (list, tuple)) is True and isinstance(da_lon, array_wrapper) is False:
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
                cond = (min(lats) <= da_lat) & (da_lat <= max(lats))
            elif isinstance(lats, (list, tuple)) is False and isinstance(lons, (list, tuple)) is True:
                cond = (min(lons) <= da_lon) & (da_lon <= max(lons))
            else:
                cond = (min(lats) <= da_lat) & (da_lat <= max(lats)) & (min(lons) <= da_lon) & (da_lon <= max(lons))
            # mask data outside region
            da = xarray_base.where(ds, cond, data_var=data_var, **kwargs_where)
            ds_o = xarray_base.set_array_in_place(ds, da, data_var=data_var)
        else:
            # -- polygonal region
            # create region using regionmask
            region = numpy__array([[lo, la] for la, lo in zip(lats, lons)])
            region = regionmask.Regions([region])
            mask = region.mask(da_lon, da_lat)
            # mask data outside region
            da = xarray_base.where(ds, xarray_base.notnull(mask), data_var=data_var, **kwargs_where)
            ds_o = xarray_base.set_array_in_place(ds, da, data_var=data_var)
        # -- select region
        if basics.is_dim(dim_lon) and isinstance(mask_only, bool) and not mask_only and isinstance(lons, (list, tuple)):
            # -- roll longitude
            lon_min, lon_max = min(lons), max(lons)
            # desired longitudes are usually defined [0; 360], but the input may not be, roll longitude if needed
            # new minimum value for the longitude will be the minimum longitude of the given region
            new_lon_min = copy__deepcopy(lon_min)
            if lon_max - new_lon_min < 360:
                # modify minimum to be slightly lower, a maximum of 10 degree lower
                new_lon_min -= min(10., 360 - (lon_max - lon_min) / 2)
            # update longitude and roll, i.e., add minimum value to dataset's longitude to shift the dimension
            # e.g., initial longitude = [0; 360], desired = [-60; 30], new_lon_min = -70, new longitude = [-70; 290]
            ds_o = roll_longitude(ds_o, new_lon_min=new_lon_min)
            # -- select region (i.e., reduce the shape of the input data)
            # create indexers
            indexers = {}
            if len(xarray_base.get_array_shape(da_lat)) == 1 and len(xarray_base.get_array_shape(da_lon)) == 1:
                # regular grid
                if isinstance(lats, (list, tuple)):
                    indexers[dim_lat] = slice(*(min(lats), max(lats)))
                if isinstance(lons, (list, tuple)):
                    indexers[dim_lon] = slice(*(min(lons), max(lons)))
            else:
                # non-regular grid -> dataarray's lat/lon must be j/i or y/x or something like that
                # create condition
                if isinstance(lats, (list, tuple)) is True and isinstance(lons, (list, tuple)) is False:
                    cond = (min(lats) <= da_lat) & (da_lat <= max(lats))
                elif isinstance(lats, (list, tuple)) is False and isinstance(lons, (list, tuple)) is True:
                    cond = (min(lons) <= da_lon) & (da_lon <= max(lons))
                else:
                    cond = (min(lats) <= da_lat) & (da_lat <= max(lats)) & (min(lons) <= da_lon) & \
                           (da_lon <= max(lons))
                # find lat/lon dimensions in dataarray's
                dim_y, dim_x = xarray_base.get_dim_keys(da_lon)
                arr_y, arr_x = xarray_base.get_dim_array(ds, dim_y), xarray_base.get_dim_array(ds, dim_x)
                # change lat/lon condition to dataarray's dimensions condition
                if isinstance(lats, (list, tuple)):
                    reg_y = xarray_base.where(arr_y, cond)
                    reg_y = (min_global(reg_y), max_global(reg_y))
                    indexers[dim_y] = slice(*(min(reg_y), max(reg_y)))
                if isinstance(lons, (list, tuple)):
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
        elif error_ds(ds, inspect__stack(), message="cannot select time_bounds"):
            break
        # check given time_bounds
        if isinstance(time_bounds, int):
            time_bounds = (time_bounds,)
        elif not any([isinstance(time_bounds, (list, tuple)) is True and 0 < len(time_bounds) < 3 and
                      all([isinstance(k, int) for k in time_bounds]) is True,
                      isinstance(time_bounds, (list, tuple)) is True and len(time_bounds) == 2 and \
                      all([isinstance(k, str) for k in time_bounds]) is True]):
            # time bounds can be:
            #    - int: the first ‘time_bounds’ time steps will be selected; e.g., 120
            #    - tuple[int]: the first ‘time_bounds’ time steps will be selected; e.g., (120,)
            #    - tuple[int, int]: time steps within ‘time_bounds’ will be selected; e.g., (120, 240)
            #    - tuple[str, str]: time steps within ‘time_bounds’ will be selected; e.g., ("1980-01-01", "2014-12-31")
            # given ‘time_bounds’ format is wrong
            log_debug(inspect__stack(), "WARNING cannot select time_bounds", details={
                "time_bounds": str(time_bounds) + " should be int, tuple[str, str], tuple[int], tuple[int, int]"})
            break
        # get time dimension name
        dim_time = get_dim_time(ds)
        log_debug(inspect__stack(), "get_dim_time", data_var=data_var, details={"dim_time": dim_time}, ds=ds)
        if not basics.is_dim(dim_time):
            break
        # check given time bounds type
        if isinstance(time_bounds, (list, tuple)) is True and len(time_bounds) == 2 and \
                all([isinstance(k, str) for k in time_bounds]) is True:
            # select using time_bounds like ("1980-01-01", "2014-12-31")
            ds_o = xarray_base.select(ds, {dim_time: slice(*time_bounds)}, **kwargs_sel)
            log_debug(inspect__stack(), "xarray_base.select", data_var=data_var,
                      details={"available time_bounds": xarray_base.get_time_bounds(ds)}, ds=ds)
            # sometimes selecting time is slightly wrong
            # this section checks if one time step has not been included by error at the beginning or the end of the
            # time series
            # check lower time bound
            ds_o = check_time_bounds(ds_o, dim_time, time_bounds, "lower")
            # check upper time bound
            ds_o = check_time_bounds(ds_o, dim_time, time_bounds, "upper")
            log_debug(inspect__stack(), "check_time_bounds", data_var=data_var,
                      details={"available time_bounds": xarray_base.get_time_bounds(ds_o)}, ds=ds_o)
        else:
            # select using time_bounds like (12, 24)
            ds_o = xarray_base.select_index(ds, {dim_time: slice(*time_bounds)}, **kwargs_sel)
            log_debug(inspect__stack(), "xarray_base.select_index", data_var=data_var, ds=ds_o)
    log_debug(inspect__stack(), "output", data_var=data_var, ds=ds_o)
    return ds_o


def smooth_along_dimension(
        ds: Union[array_wrapper, dataset_wrapper],
        cf_dim: Literal["T", "X", "Y", "Z"],
        method: Literal["triangular", "uniform"],
        window: int,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    # here is an example of solution
    # https://stackoverflow.com/questions/48510784/xarray-rolling-mean-with-weights
    # https://docs.xarray.dev/en/stable/generated/xarray.computation.rolling.DataArrayRolling.construct.html
    # step 1: compute weights
    # triangle: 121
    # uniform: 111
    # get time / lat / lon weights
    # combine dimension weights and window weights
    # step 2: average
    # ds2 = ds * weights
    # https://docs.xarray.dev/en/stable/generated/xarray.DataArray.rolling.html
    # ds2.rolling(dim={dim: window}, min_periods=min_periods, center=True).mean()
    # weights.rolling(dim={dim: window}, min_periods=min_periods, center=True).mean()
    # ds2 / weights
    # da.rolling(dim_0=3, center=True).construct('window').dot(weight)
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
    log_debug(inspect__stack(), "input", details={"ds.type": type(ds), "cf_dim": cf_dim})
    # fake loop to be able to break out
    for _ in [0]:
        if error_ds(ds, inspect__stack(), message="cannot squeeze"):
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
                if isinstance(ds, dataset_wrapper):
                    ds = xarray_base.drop_dataset_keys(
                        ds, [str(dim_name) + "_bounds", str(dim_name) + "_bnds"], errors="ignore")
                    log_debug(inspect__stack(), "xarray_base.drop_dataset_keys", details={
                        "dataset_keys": xarray_base.get_dataset_keys(ds)})
    log_debug(inspect__stack(), "output", details={"dim_keys": xarray_base.get_dim_keys(ds)})
    return ds
# ---------------------------------------------------------------------------------------------------------------------#

