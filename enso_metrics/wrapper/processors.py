# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Processors built over xarray and xcdat
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from glob import iglob as glob__iglob
from inspect import stack as inspect__stack
from json import dumps as json__dumps
from os.path import isdir as os__path__isdir
from os.path import isfile as os__path__isfile
from typing import Literal, Union, Hashable

# local functions
from enso_metrics.tools.default import set_instance
from enso_metrics.wrapper import basics
from enso_metrics.wrapper import wrapper_base as wb
from enso_metrics.wrapper import xarray_base as xab
from enso_metrics.wrapper import xcdat_base as xcb
from enso_metrics.wrapper.xarray_base import array_wrapper, dataset_wrapper
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: processors
# ---------------------------------------------------------------------------------------------------------------------#
def averager(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        cf_dims: list[Literal["T", "X", "Y"]] = None,
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_average_spatial: Union[dict, None] = None,
        kwargs_average_temporal: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "cf_dims", "kwargs_average_spatial",
          "kwargs_average_temporal"]
    l2 = [input_array, input_area, data_var, data_var_area, cf_dims, kwargs_average_spatial, kwargs_average_temporal]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # is cf_dims defined
        if cf_dims is None or (isinstance(cf_dims, (list, tuple)) is True and len(cf_dims) == 0):
            output_array, output_area = input_array, input_area
            # WARNING: user calls averager without providing a dimension to average
            message = "WARNING: user calls averager without providing a dimension to average"
            details = basics.log_details(["cf_dims"], [cf_dims])
            wb.log_debug(inspect__stack(), message, adjust=5, details=details)
            break
        # check desired average
        if isinstance(cf_dims, (list, tuple)) is False or (
                isinstance(cf_dims, (list, tuple)) is True and all([isinstance(k, str) for k in cf_dims]) is False) or (
                isinstance(cf_dims, (list, tuple)) is True and all([isinstance(k, str) for k in cf_dims]) is True and
                all([k in ["T", "X", "Y"] for k in cf_dims]) is False):
            # given ‘cf_dims’ format is wrong
            wb.log_debug(inspect__stack(), "WARNING cannot perform averager", details={
                "cf_dims": str(cf_dims) + " should be list[Literal['T', 'X', 'Y']]"})
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform averager") is True:
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) is True and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # perform average
        if "T" in cf_dims:
            try:
                ds_array = xcb.average_temporal(ds_array, data_var, **kwargs_average_temporal)
            except Exception as err:
                wb.log_debug(inspect__stack(), "WARNING cannot perform temporal average using xcdat\n" + str(err))
                break
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            text = "temporal average computed"
            metadata["description"] = basics.description_writer(description, text)
        if "X" in cf_dims or "Y" in cf_dims:
            # cf dims
            cf_dims_tmp = [k for k in cf_dims if k in ["X", "Y"]]
            # spatial average
            ds_array = wb.average_spatial(ds_array, cf_dim=cf_dims_tmp, data_var=data_var, data_var_area=data_var_area,
                                          ds_area=ds_area)
            if ds_array is None:
                wb.log_debug(inspect__stack(), "WARNING cannot perform spatial average")
                break
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            text = "zonal" if len(list(set(cf_dims) - {"X"})) == 0 else (
                "meridional" if len(list(set(cf_dims) - {"Y"})) == 0 else "horizontal")
            text += " average computed"
            metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)) is True:
            output_array = {"array": ds_array, "metadata": metadata}
        # TO DO: average area?
    return output_array, output_area


def detrender(
        ds: dataset_wrapper,
        degree: int,
        **kwargs) -> dataset_wrapper:
    basics.log_info(inspect__stack(), "")
    return ds


def masker(
        input_array: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        data_var_mask: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        input_mask: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        maskland: bool = False,
        maskocean: bool = False,
        region: Union[str, None] = None,
        regions_param: dict[str, dict[str, Union[bool, str, tuple]]] = None,
        tolerance: Union[float, int] = 0.,
        kwargs_where: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], Union[dict, None]):
    basics.log_info(inspect__stack(), "")
    # update input data using region
    if isinstance(region, str) is True and isinstance(regions_param, dict) is True and \
            region in list(regions_param.keys()):
        maskland, maskocean = regions_param[region]["maskland"], regions_param[region]["maskocean"]
    l1 = ["input_array", "input_area", "input_mask", "data_var", "data_var_area", "data_var_mask", "maskland",
          "maskocean", "region", "regions_param", "tolerance", "kwargs_where"]
    l2 = [input_array, input_area, input_mask, data_var, data_var_area, data_var_mask, maskland, maskocean, region,
          regions_param, tolerance, kwargs_where]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to xarray.Dataset.where (used to mask data): drop, other.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_where'.
    kwargs_where = set_instance(kwargs_where, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # mask land or ocean
        if (maskland is True and maskocean is True) or (maskland is False and maskocean is False):
            output_array, output_area = input_array, input_area
            # WARNING: user asked to mask everything or nothing
            message = "WARNING: user asked to mask "
            message += "everything" if maskland is True and maskocean is True else "nothing"
            wb.log_debug(inspect__stack(), message, adjust=5)
            break
        # -- get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform masker") is True:
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) is True and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        ds_mask, metadata_mask = None, None
        if isinstance(input_mask, dict) is True and "array" in list(input_mask.keys()) and \
                "metadata" in list(input_mask.keys()):
            ds_mask, metadata_mask = input_mask["array"], input_mask["metadata"]
        # -- get mask DataArray
        if isinstance(ds_mask, (array_wrapper, dataset_wrapper)) is False or (
                isinstance(ds_mask, dataset_wrapper) is True and isinstance(data_var_mask, (Hashable, str)) is False):
            # mask DataArray must be available to mask data
            wb.log_debug(inspect__stack(), "WARNING cannot perform masker", details={
                "ds_mask": str(type(ds_mask)) + " should be xarray.DataArray, or xarray.Dataset and data_var_mask " +
                           "given",
                "data_var_mask": str(data_var_mask)})
            break
        # get array
        da_mask = xab.to_array(ds_mask, data_var=data_var_mask)
        if isinstance(da_mask, array_wrapper) is False:
            # mask DataArray must be available to mask data
            details = {
                "da_mask": str(type(da_mask)) + " should be xarray.DataArray",
                "ds_mask": str(type(ds_mask)),
                "data_var_mask": str(data_var_mask)
            }
            if isinstance(ds_mask, dataset_wrapper) is True:
                details["ds_mask.keys"] = str(xab.get_dataset_keys(ds_mask))
            wb.log_debug(inspect__stack(), "WARNING cannot perform masker", details=details)
        # maximum value
        da_mask_max = wb.max_global(da_mask)
        # units
        da_mask_units = ""
        if "units" in list(metadata_mask.keys()):
            da_mask_units = metadata_mask["units"]
        # if land = 100 instead of 1, divide ds_mask by 100
        if da_mask_max > 1 or da_mask_units == "%":
            da_mask /= 100.
        # apply mask
        wb.log_debug(inspect__stack(), "input array", adjust=5, ds=ds_array, data_var=data_var)
        wb.log_debug(inspect__stack(), "landmask", adjust=5, ds=da_mask)
        # by definition: mask is 1 on land & 0 on ocean
        # reverse it if land must be masked (as ocean is kept, weights must be 1 on ocean and 0 on land)
        if maskland is True:
            da_mask = 1 - da_mask
        # now mask is 1 where we want to keep data
        # where function keeps data where the given condition is True
        # if tolerance is 0, only cells that are 100% land or ocean are kept
        # else some mix cells will we kept and if tolerance is 1 only opposite cells are masked
        # e.g., if maskland = True and tolerance = 1, only cells that are 100% land are masked
        condition = da_mask != 1 if tolerance == 0 else da_mask > 1 - tolerance
        ds_array = xab.where(ds_array, condition, **kwargs_where)
        if isinstance(ds_area, array_wrapper) is True or (
                isinstance(ds_area, (array_wrapper, dataset_wrapper)) and
                isinstance(data_var_area, (Hashable, str)) is True):
            wb.log_debug(inspect__stack(), "areacell", adjust=5, ds=ds_area, data_var=data_var_area)
            # for the area dataset, values are not masked but set to 0 (xarray.DataArray.weighted does not accept
            # missing values)
            tmp_kwargs_where = copy__deepcopy(kwargs_where)
            if "other" in list(tmp_kwargs_where.keys()):
                del tmp_kwargs_where["other"]
            ds_area = xab.where(ds_area, condition, other=0, **kwargs_where)
            # multiply area by mask so that area take into account the fraction of land or ocean
            da_area = xab.to_array(ds_area, data_var=data_var_area)
            da_area = da_area * da_mask
            ds_area = xab.assign_to_dataset(ds_area, variables_kwargs={data_var_area: da_area})
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        text = "land masked (mask " if maskland is True else "ocean masked (mask "
        if maskland is True:
            text += "!= 0)" if tolerance == 0 else ">= " + str(tolerance) + ")"
        else:
            text += "!= 1)" if tolerance == 0 else "<= " + str(1 - tolerance) + ")"
        if data_var_mask == "estimate":
            text = text.replace(" masked (mask ", " masked (estimated mask ")
        metadata["description"] = basics.description_writer(description, text)
        wb.log_debug(inspect__stack(), "output array (xarray_base.where)", adjust=5, ds=ds_array,
                     data_var=data_var, details={"description": metadata["description"]})
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)) is True:
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)) is True:
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


def selector(
        input_array: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        depth_bounds: Union[tuple[Union[float, int], Union[float, int]], None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        region: Union[str, None] = None,
        regions_param: dict[str, dict[str, Union[bool, str, tuple[Union[float, int]]]]] = None,
        time_bounds: Union[tuple[str, str], None] = None,
        kwargs_select_depth: Union[dict, None] = None,
        kwargs_select_horizontal: Union[dict, None] = None,
        kwargs_select_time: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], Union[dict, None]):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "depth_bounds", "region", "regions_param",
          "time_bounds", "kwargs_select_depth", "kwargs_select_horizontal", "kwargs_select_time"]
    l2 = [input_array, input_area, data_var, data_var_area, depth_bounds, region, regions_param,
          time_bounds, kwargs_select_depth, kwargs_select_horizontal, kwargs_select_time]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to wrapper_base.select_horizontal (used to select region):
    #      - mask_only: True to only mask data (for spatial average), False to select region
    #      - kwargs_sel (to select on regular grid): drop, method, tolerance (see xarray_base.select)
    #      - kwargs_where (to mask non-regular grid): drop (see xarray_base.where)
    # If desired they must be defined in a dictionary under the keyword 'kwargs_select_horizontal'.
    kwargs_select_horizontal = set_instance(kwargs_select_horizontal, dict, False, {})
    # One argument (a dictionary containing more arguments) can be specified to wrapper_base.select_depth (used to
    # select depth) and to wrapper_base.select_time (used to select time): kwargs_sel.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_select_depth' or 'kwargs_select_time'.
    kwargs_select_depth = set_instance(kwargs_select_depth, dict, False, {})
    kwargs_select_time = set_instance(kwargs_select_time, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # check if depth_bounds is defined
        do_z = True
        if isinstance(depth_bounds, (list, tuple)) is False or (
                isinstance(depth_bounds, (list, tuple)) is True and (
                len(depth_bounds) != 2 or all([isinstance(k, (float, int)) for k in depth_bounds]) is False)):
            do_z = False
        # check if time_bounds is defined
        do_t = True
        if isinstance(time_bounds, (list, tuple)) is False or (
                isinstance(time_bounds, (list, tuple)) is True and (
                len(time_bounds) != 2 or all([isinstance(k, str) for k in time_bounds]) is False)):
            do_t = False
        # check if region is defined
        do_xy = True
        if isinstance(region, str) is False or isinstance(regions_param, dict) is False or (
                isinstance(regions_param, dict) is True and region not in list(regions_param.keys())):
            do_xy = False
        if do_t is False and do_xy is False and do_z is False:
            # WARNING: user calls selector without providing something to select
            message = "WARNING: user calls selector without providing something to select"
            l1 = ["depth_bounds", "region", "regions_param", "time_bounds"]
            l2 = [depth_bounds, region, regions_param, time_bounds]
            details = basics.log_details(l1, l2)
            wb.log_debug(inspect__stack(), message, adjust=5, details=details)
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) is True and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # select time
        if do_t is True:
            ds_array = wb.select_time(ds_array, data_var=data_var, time_bounds=time_bounds, **kwargs_select_time)
            if ds_array is None:
                break
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            text = "time selected (" + str(time_bounds) + ")"
            metadata["description"] = basics.description_writer(description, text)
        # select depth
        if do_z is True:
            ds_array = wb.select_depth(ds_array, data_var=data_var, depth_bounds=time_bounds, **kwargs_select_depth)
            if ds_array is None:
                break
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            text = "time selected (" + str(kwargs_select_time) + ")"
            metadata["description"] = basics.description_writer(description, text)
        # select region
        if do_xy is True:
            # region
            horizontal_bounds = {}
            if "latitude" in list(regions_param[region].keys()) and \
                    isinstance(regions_param[region]["latitude"], (list, tuple)) is True and \
                    len(regions_param[region]["latitude"]) >= 2:
                horizontal_bounds["Y"] = regions_param[region]["latitude"]
            if "longitude" in list(regions_param[region].keys()) and \
                    isinstance(regions_param[region]["longitude"], (list, tuple)) is True and \
                    len(regions_param[region]["longitude"]) >= 2:
                horizontal_bounds["X"] = regions_param[region]["longitude"]
            # mask or select region
            ds_array = wb.select_horizontal(ds_array, data_var=data_var, horizontal_bounds=horizontal_bounds,
                                            **kwargs_select_horizontal)
            if ds_array is None:
                break
            if isinstance(ds_area, (array_wrapper, dataset_wrapper)) is True:
                ds_area = wb.select_horizontal(ds_area, data_var=data_var_area, horizontal_bounds=horizontal_bounds,
                                               **kwargs_select_horizontal)
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            n_sho, n_lon = copy__deepcopy(region), ""
            if "short_name" in list(regions_param[region].keys()):
                n_sho = regions_param[region]["short_name"]
            if "long_name" in list(regions_param[region].keys()):
                n_lon = " (" + str(regions_param[region]["long_name"]) + ")"
            text = "selected in " + str(n_sho) + str(n_lon)
            if "Y" in list(horizontal_bounds.keys()) and "X" in list(horizontal_bounds.keys()) and \
                    len(horizontal_bounds["Y"]) == 2 and len(horizontal_bounds["X"]) == 2:
                text += " " + str(basics.write_coordinates(*horizontal_bounds["Y"], *horizontal_bounds["X"]))
            metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)) is True:
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)) is True:
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: loop on processors
# ---------------------------------------------------------------------------------------------------------------------#
dict_processors = {}
for key in list(locals().keys()):
    if callable(locals()[key]) and locals()[key].__module__ == __name__:
        dict_processors[key] = locals()[key]
list_processors = sorted(list(dict_processors.keys()), key=lambda s: s.lower())


def loop(processors, input_dataset, input_param, **kwargs) -> Union[dict, None]:
    basics.log_info(inspect__stack(), "")
    print("processors", list(processors.keys()))
    print("list_processors", list_processors)
    # fx variables
    fx = dict((k1, d1) for k1, d1 in input_dataset.items() if k1 in ["areacell", "areacella", "areacello", "landmask"])
    # loop on variables to process
    dict_o = {}
    for k1 in list(processors.keys()):
        print(k1)
        # input and output variable names
        variable_i = processors[k1]["variable"]
        region_i = processors[k1]["region"]
        variable_o = copy__deepcopy(k1)
        # check if given variable is available
        if variable_i not in list(input_dataset.keys()) or isinstance(input_param, dict) is False or \
                variable_i not in list(input_param.keys()):
            # WARNING: variable must be defined
            details = {"variable": str(variable_i),
                       "in input_dataset": str(variable_i in list(input_dataset.keys())),
                       "input_param.type": str(type(input_param))}
            if isinstance(input_param, dict) is True:
                details["input_param.keys"] = ", ".join(list(input_param.keys()))
            wb.log_debug(inspect__stack(), "WARNING: variable must be defined", adjust=5, details=details)
            break
        print(variable_i, variable_o)
        # get param given variable as well as area and mask names related to given variable
        param = input_param[variable_i]
        data_var_area = param["area"] if "area" in list(param.keys()) else None
        data_var_mask = param["mask"] if "mask" in list(param.keys()) else None
        print("area", data_var_area, "mask", data_var_mask)
        # get variable, area and mask dictionaries (i.e., {"array": xarray.Dataset, "metadata": {}})
        dict_array = input_dataset[variable_i]
        dict_area, dict_mask = None, None
        for n1, n2 in zip(["area", "mask"], [data_var_area, data_var_mask]):
            if isinstance(n2, str) is True and (n2 not in list(input_dataset.keys()) or
                                                isinstance(input_param, dict) is False or
                                                n2 not in list(input_param.keys())):
                # WARNING: given variable must be defined
                details = {"variable": str(n2),
                           "in input_dataset": str(n2 in list(input_dataset.keys())),
                           "input_param.type": str(type(input_param))}
                if isinstance(input_param, dict) is True:
                    details["input_param.keys"] = ", ".join(list(input_param.keys()))
                message = "WARNING: variable " + str(n1) + " must be defined"
                print(message)
                wb.log_debug(inspect__stack(), message, adjust=5, details=details)
                break
            else:
                if n1 == "area" and n2 is not None:
                    dict_area = input_dataset[n2]
                elif n1 == "mask" and n2 is not None:
                    dict_mask = input_dataset[n2]
        if isinstance(dict_array, dict) and "array" in list(dict_array.keys()):
            print("array", type(dict_array["array"]))
        if isinstance(dict_area, dict) and "array" in list(dict_area.keys()):
            print("area", type(dict_area["array"]))
        # loop on processors to apply to given variable
        print("processors", list(processors[k1]["to_do"].keys()))
        for k2 in list(processors[k1]["to_do"].keys()):
            print(k2)
            process = k2.split("__")[-1]
            if process in list_processors:
                # call processor
                print(k2, "call processor")
                if process == "masker" and dict_mask is None:
                    wb.log_debug(inspect__stack(), "WARNING " + str(variable_i) + " not masked as mask not provided")
                    continue
                local_kwargs = processors[k1]["to_do"][k2]
                if isinstance(dict_array, dict) and "array" in list(dict_array.keys()):
                    print("array", type(dict_array["array"]))
                if isinstance(dict_area, dict) and "array" in list(dict_area.keys()):
                    print("area", type(dict_area["array"]))
                print("local_kwargs", list(local_kwargs.keys()))
                dict_array, dict_area = dict_processors[process](
                    dict_array, data_var=variable_i, data_var_area=data_var_area, data_var_mask=data_var_mask,
                    input_area=dict_area, input_mask=dict_mask, region=region_i, **local_kwargs, **kwargs)
                if dict_array is None:
                    print("dict_array is None -> must break")
                    break
        if dict_array is None:
            break
        dict_o[variable_o] = dict_array
    if len(dict_o.keys()) != len(processors.keys()):
        dict_o = None
    return dict_o


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: reader (to read netCDF files)
# ---------------------------------------------------------------------------------------------------------------------#
def reader(
        input_param: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, None]]]]],
        variables: Union[str, list[str]],
        kwargs_reader: dict = None,
        variables_param: dict[str, dict[str, str]] = None,
        **kwargs) -> dict:
    basics.log_info(inspect__stack(), "")
    # Several arguments can be specified to open_dataset.
    # Some are explicitly named (see enso_metrics.wrapper.xcdat_base):
    # add_bounds, center_times, data_vars, decode_times, lon_orient, preprocess
    # but more can be given.
    # All these arguments are optionals.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_reader'
    # Here add_bounds and decode_times are taken out as we impose a default value, and these two arguments may have to
    # be altered if time is not defined in a netCDF (e.g., for cmip variables like areacella, areacello, landmask)
    kwargs_reader = set_instance(kwargs_reader, dict, False, {})
    add_bounds: list[Literal["T", "X", "Y", "Z"]] = ["T", "X", "Y"]
    if "add_bounds" in list(kwargs_reader.keys()):
        add_bounds = copy__deepcopy(kwargs_reader["add_bounds"])
        del kwargs_reader["add_bounds"]
    decode_times: bool = True
    if "decode_times" in list(kwargs_reader.keys()):
        add_bounds = copy__deepcopy(kwargs_reader["decode_times"])
        del kwargs_reader["decode_times"]
    # variables to list
    if isinstance(variables, str) is True:
        variables = [variables]
    # read variable (as in file)
    dict_t = {}
    for kk in variables + ["areacella", "areacello", "landmask"]:
        if kk not in list(input_param.keys()) or kk not in list(variables_param.keys()):
            continue
        # loop on netCDF files-variables to read all required data
        for ff, nn in zip(input_param[kk]["file_name"], input_param[kk]["variable"]):
            # remove time dimension from add_bounds and decode_times if not defined in netCDF ('fx' frequency)
            ab, dt = copy__deepcopy(add_bounds), copy__deepcopy(decode_times)
            if kk in ["areacella", "areacello", "landmask"]:
                if "T" in ab:
                    while "T" in ab:
                        ab.remove("T")
                if dt is True:
                    dt = False
            # try to open_dataset and save Dataset in a dictionary using netCDF variables (names) as keys
            try:
                dict_t[nn] = xcb.open_dataset(
                    ff, add_bounds=ab, data_var=nn, decode_times=dt, **kwargs_reader)
            except Exception as err:
                message = "can't read (" + str(nn) + ") " + str(ff) + "\n" + str(err)
                basics.log_error(inspect__stack(), message)
                # WARNING: cannot read variable or file
                path = "/".join(ff.split("/")[:-1])
                files = sorted(list(glob__iglob(ff)), key=lambda s: s.lower())
                files_string = ""
                if len(files) > 0:
                    for k in files:
                        files_string += "\n" + str().ljust(5) + str(k)
                else:
                    files_string = "no file matches this file pattern"
                details = {
                    "directory": str(path),
                    "isdir": str(os__path__isdir(path)),
                    "file": str(ff),
                    "isfile": str(os__path__isfile(ff)),
                    "list": str(files_string)}
                wb.log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
    # compute variables
    dict_output = {}
    for kk in variables + ["areacella", "areacello", "landmask"]:
        if kk not in list(input_param.keys()) or kk not in list(variables_param.keys()):
            continue
        arr, metadata = wb.compute_variable(
            dict_t, input_param[kk], input_param[kk]["variable"], kk, variables_param[kk], **kwargs)
        dict_output[kk] = {"array": arr, "metadata": metadata}
    # check if landmask must be estimated
    for k1 in variables:
        if k1 not in list(dict_output.keys()) or k1 not in list(input_param.keys()) or \
                isinstance(input_param[k1], dict) is False or "mask" not in list(input_param[k1].keys()):
            continue
        data_var_mask = input_param[k1]["mask"]
        if data_var_mask is not None and data_var_mask not in list(dict_output.keys()):
            data_var_mask = "estimate"
            input_param[k1]["mask"] = copy__deepcopy(data_var_mask)
            arr, metadata = wb.compute_mask(dict_output[k1], data_var=k1, data_var_mask=data_var_mask)
            dict_output[data_var_mask] = {"array": arr, "metadata": metadata}
    if len([k for k in list(dict_output.keys()) if k in variables]) != len(variables):
        dict_output = None
    return dict_output
# ---------------------------------------------------------------------------------------------------------------------#
