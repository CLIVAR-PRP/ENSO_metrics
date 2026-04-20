# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Processors built over xarray and xcdat
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from dataclasses import dataclass as dataclasses__dataclass
from inspect import stack as inspect__stack
from os import remove as os__remove
from typing import Annotated, Literal, Union, Hashable

# local functions
from enso_metrics.tools.default import set_instance
from enso_metrics.wrapper import basics
from enso_metrics.wrapper import wrapper_base as wb
from enso_metrics.wrapper import xarray_base as xab
from enso_metrics.wrapper import xcdat_base as xcb
from enso_metrics.wrapper.xarray_base import array_wrapper, dataset_wrapper
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Classes: range
# ---------------------------------------------------------------------------------------------------------------------#
@dataclasses__dataclass
class IntRange:
    min: int
    max: int
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: processors
# ---------------------------------------------------------------------------------------------------------------------#
def anomaler(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_anomalies: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "kwargs_anomalies"]
    l2 = [input_array, input_area, data_var, data_var_area, kwargs_anomalies]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to xcdat_base.interannual_anomalies: frequency, keep_weights, reference_period,
    # season_config, skipna, weighted.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_anomalies'.
    kwargs_anomalies = set_instance(kwargs_anomalies, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform anomaler"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # perform interannual anomalies
        try:
            ds_array = xcb.interannual_anomalies(ds_array, data_var, **kwargs_anomalies)
        except Exception as err:
            print("WARNING cannot perform interannual_anomalies using xcdat")
            print(err)
            wb.log_debug(inspect__stack(), "WARNING cannot perform interannual_anomalies using xcdat\n" + str(err))
            break
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        metadata["description"] = basics.description_writer(description, "seasonal cycle removed")
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


def averager(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        cf_dims: list[Literal["T", "X", "XY", "Y"]] = None,
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_average_temporal: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "cf_dims", "kwargs_average_temporal"]
    l2 = [input_array, input_area, data_var, data_var_area, cf_dims, kwargs_average_temporal]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to xcdat_base.average_temporal: keep_weights, skipna, weighted.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_average_temporal'.
    kwargs_average_temporal = set_instance(kwargs_average_temporal, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # is cf_dims defined
        if cf_dims is None or (isinstance(cf_dims, (list, tuple)) and len(cf_dims) == 0):
            output_array, output_area = input_array, input_area
            # WARNING: user calls averager without providing a dimension to average
            message = "WARNING: user calls averager without providing a dimension to average"
            details = basics.log_details(["cf_dims"], [cf_dims])
            wb.log_debug(inspect__stack(), message, adjust=5, details=details)
            break
        # check desired average
        if not isinstance(cf_dims, (list, tuple)) or (
                isinstance(cf_dims, (list, tuple)) and not all([isinstance(k, str) for k in cf_dims])) or (
                isinstance(cf_dims, (list, tuple)) and all([isinstance(k, str) for k in cf_dims]) and
                not all([k in ["T", "X", "XY", "Y"] for k in cf_dims])):
            # given ‘cf_dims’ format is wrong
            wb.log_debug(inspect__stack(), "WARNING cannot perform averager", details={
                "cf_dims": str(cf_dims) + " should be list[Literal['T', 'X', 'XY', 'Y']]"})
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform averager"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # perform average
        for cf_dim in cf_dims:
            if cf_dim == "T":
                try:
                    ds_array = xcb.average_temporal(ds_array, data_var, **kwargs_average_temporal)
                except Exception as err:
                    wb.log_debug(inspect__stack(), "WARNING cannot perform temporal average using xcdat\n" + str(err))
                    break
            else:
                # cf dims
                cf_dims_tmp: list[Literal["X", "Y"]] = ["X", "Y"] if cf_dim == "XY" else [cf_dim]
                # spatial average
                ds_array = wb.average_spatial(ds_array, cf_dim=cf_dims_tmp, data_var=data_var,
                                              data_var_area=data_var_area, ds_area=ds_area)
                if ds_array is None:
                    wb.log_debug(inspect__stack(), "WARNING cannot perform spatial average")
                    break
            # redo bounds (operations on ds changes bounds and xcdat doesn't like that)
            # cf_dim_to_redo = [k for k in ["T", "X", "Y", "Z"] if k != cf_dim or (k in ["X", "Y"] and cf_dim != "XY")]
            # ds_array = wb.redo_bounds(ds_array, cf_dim_to_redo)
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            if cf_dim == "T":
                text = "temporal average computed"
            else:
                text = "zonal" if cf_dim == "X" else ("meridional" if cf_dim == "Y" else "horizontal")
                text += " average computed"
            metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        # TO DO: average area?
    return output_array, output_area


def detrender(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        degree: Annotated[int, IntRange(0, 3)] = 1,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_detrend: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "degree", "kwargs_detrend"]
    l2 = [input_array, input_area, data_var, data_var_area, degree, kwargs_detrend]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to wrapper_base.remove_fit: deg, dim, kwargs_polyfit.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_detrend'.
    kwargs_detrend = set_instance(kwargs_detrend, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # check polynomial degree
        if not isinstance(degree, int) or not (0 <= degree <= 3):
            # given ‘degree’ format is wrong
            wb.log_debug(inspect__stack(), "WARNING cannot perform detrender", details={
                "degree": str(degree) + " should be 0 <= degree <= 3"})
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform detrender"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # perform detrend
        ds_array = wb.remove_fit(ds_array, data_var=data_var, deg=degree, dim="T", kwargs_polyfit=kwargs_detrend)
        if ds_array is None:
            wb.log_debug(inspect__stack(), "WARNING cannot perform detrender")
            break
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        text = "time mean value removed"
        if degree > 0:
            text = "linearly" if degree == 1 else ("quadratically" if degree == 2 else "cubically")
            text = "time series " + str(text) + " detrended"
        metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


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
    if isinstance(region, str) and isinstance(regions_param, dict) and region in list(regions_param.keys()):
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
        if isinstance(maskland, bool) and isinstance(maskocean, bool) and maskland == maskocean:
            output_array, output_area = input_array, input_area
            # WARNING: user asked to mask everything or nothing
            message = "WARNING: user asked to mask "
            message += "everything" if maskland and maskocean else "nothing"
            wb.log_debug(inspect__stack(), message, adjust=5)
            break
        # -- get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform masker"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        ds_mask, metadata_mask = None, None
        if isinstance(input_mask, dict) and "array" in list(input_mask.keys()) and \
                "metadata" in list(input_mask.keys()):
            ds_mask, metadata_mask = input_mask["array"], input_mask["metadata"]
        # -- get mask DataArray
        if not isinstance(ds_mask, (array_wrapper, dataset_wrapper)) or (
                isinstance(ds_mask, dataset_wrapper) and not isinstance(data_var_mask, (Hashable, str))):
            if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
                output_array = {"array": ds_array, "metadata": metadata}
            if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
                output_area = {"array": ds_area, "metadata": metadata_area}
            # mask DataArray must be available to mask data
            wb.log_debug(inspect__stack(), "WARNING cannot perform masker", details={
                "ds_mask": str(type(ds_mask)) + " should be xarray.DataArray, or xarray.Dataset and data_var_mask " +
                           "given",
                "data_var_mask": str(data_var_mask)})
            break
        # get array
        da_mask = xab.to_array(ds_mask, data_var=data_var_mask)
        if not isinstance(da_mask, array_wrapper):
            # mask DataArray must be available to mask data
            details = {
                "da_mask": str(type(da_mask)) + " should be xarray.DataArray",
                "ds_mask": str(type(ds_mask)),
                "data_var_mask": str(data_var_mask)
            }
            if isinstance(ds_mask, dataset_wrapper):
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
        if isinstance(maskland, bool) and maskland:
            da_mask = 1 - da_mask
        # now mask is 1 where we want to keep data
        # where function keeps data where the given condition is True
        # if tolerance is 0, only cells that are 100% land or ocean are kept
        # else some mix cells will we kept and if tolerance is 1 only opposite cells are masked
        # e.g., if maskland = True and tolerance = 1, only cells that are 100% land are masked
        condition = da_mask == 1 if tolerance == 0 else da_mask > 1 - tolerance
        da_array = xab.where(ds_array, condition, data_var=data_var, **kwargs_where)
        ds_array = xab.set_array_in_place(ds_array, da_array, data_var=data_var)
        if isinstance(ds_area, array_wrapper) or (
                isinstance(ds_area, dataset_wrapper) and isinstance(data_var_area, (Hashable, str))):
            wb.log_debug(inspect__stack(), "areacell", adjust=5, ds=ds_area, data_var=data_var_area)
            # for the area dataset, values are not masked but set to 0 (xarray.DataArray.weighted does not accept
            # missing values)
            tmp_kwargs_where = copy__deepcopy(kwargs_where)
            if "other" in list(tmp_kwargs_where.keys()):
                del tmp_kwargs_where["other"]
            da_area = xab.where(ds_area, condition, data_var=data_var_area, **kwargs_where)
            # multiply area by mask so that area take into account the fraction of land or ocean
            da_area = da_area * da_mask
            ds_area = xab.set_array_in_place(ds_area, da_area, data_var=data_var_area)
            # ds_area = xab.assign_to_dataset(ds_area, variables_kwargs={data_var_area: da_area})
            # ds_area = wb.redo_bounds(ds_area, ["X", "Y"])
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        text = "land masked (mask " if isinstance(maskland, bool) and maskland else "ocean masked (mask "
        if isinstance(maskland, bool) and maskland:
            text += "!= 0)" if tolerance == 0 else ">= " + str(tolerance) + ")"
        else:
            text += "!= 1)" if tolerance == 0 else "<= " + str(1 - tolerance) + ")"
        if data_var_mask == "estimate":
            text = text.replace(" masked (mask ", " masked (estimated mask ")
        metadata["description"] = basics.description_writer(description, text)
        wb.log_debug(inspect__stack(), "output array (xarray_base.where)", adjust=5, ds=ds_array,
                     data_var=data_var, details={"description": metadata["description"]})
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


def regridder(
        input_array: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None],
        cf_dims: Union[list[Literal["X", "XY", "Y", "Z"]], Literal["X", "XY", "Y", "Z"], None] = None,
        data_var: Union[str, None] = None,
        grid_xy: Union[array_wrapper, dataset_wrapper, str, None] = None,
        grid_z: Union[array_wrapper, dataset_wrapper, str, None] = None,
        kwargs_regrid_xy: Union[dict, None] = None,
        kwargs_regrid_z: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], Union[dict, None]):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "data_var", "grid_xy", "grid_z", "kwargs_regrid_xy", "kwargs_regrid_z"]
    l2 = [input_array, data_var, grid_xy, grid_z, kwargs_regrid_xy, kwargs_regrid_z]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to wrapper_base.regrid (used to regrid horizontally):
    #      - output_grid: xr.Dataset grid or standard grid name to transform inputs to
    #      - method: regridding method (see xcdat_base.regrid_horizontal)
    #      - tool: name of the tool to use.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_regrid_xy'.
    kwargs_regrid_xy = set_instance(kwargs_regrid_xy, dict, False, {})
    # Several arguments can be specified to wrapper_base.regrid (used to regrid vertically):
    #      - output_grid: xr.Dataset grid or standard number of depth levels to transform inputs to
    #      - tool: name of the tool to use (only 'xgcm' at the moment)
    #      - **options: options passed directly to the tool (see specific regridder for available options)
    kwargs_regrid_z = set_instance(kwargs_regrid_z, dict, False, {})
    # Dimension(s) along which to apply
    if isinstance(cf_dims, str):
        cf_dims = [cf_dims]
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # check if grid_h is defined
        do_xy = False
        if isinstance(cf_dims, list) and "XY" in cf_dims and (
                isinstance(grid_xy, (array_wrapper, dataset_wrapper)) or
                (isinstance(grid_xy, str) and ("gaussian" in grid_xy or "uniform" in grid_xy))):
            do_xy = True
        # check if grid_z is defined
        do_z = False
        if isinstance(cf_dims, list) and "Z" in cf_dims and (
                isinstance(grid_z, (array_wrapper, dataset_wrapper)) or
                (isinstance(grid_z, int) and grid_z in [50, 60, 75])):
            grid_z: Union[array_wrapper, dataset_wrapper, str, None]
            do_z = True
        if not do_xy and not do_z:
            # WARNING: user calls selector without providing something to select
            message = "WARNING: user calls regridder without providing something to regrid"
            l1 = ["cf_dims", "grid_xy", "grid_z"]
            l2 = [cf_dims, grid_xy, grid_z]
            details = basics.log_details(l1, l2)
            wb.log_debug(inspect__stack(), message, adjust=5, details=details)
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        # regrid
        for cf_dim in cf_dims:
            if (cf_dim == "XY" and not do_xy) or (cf_dim == "Z" and not do_z):
                continue
            # regrid
            if cf_dim == "XY":
                ds_array = wb.regrid_horizontal(ds_array, data_var=data_var, grid=grid_xy, **kwargs_regrid_xy)
            else:
                ds_array = wb.regrid_vertical(ds_array, data_var=data_var, grid=grid_z, **kwargs_regrid_z)
            if ds_array is None:
                break
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            if cf_dim == "XY":
                text = "horizontally reggrided ("
                if isinstance(grid_xy, str) and ("gaussian" in grid_xy or "uniform" in grid_xy):
                    text += str(grid_xy) + ")"
                else:
                    text += "provided grid)"
            else:
                text = "vertically reggrided ("
                if isinstance(grid_z, int) and grid_z in [50, 60, 75]:
                    text += "standard_" + str(grid_z) + "_depth_level)"
                else:
                    text += "provided grid)"
            metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
    return output_array, output_area


def remover(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_averager: Union[dict, None] = None,
        kwargs_selector: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "kwargs_averager", "kwargs_selector"]
    l2 = [input_array, input_area, data_var, data_var_area, kwargs_averager, kwargs_selector]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to processors.averager: cf_dims, kwargs_average_temporal.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_averager'.
    kwargs_averager = set_instance(kwargs_averager, dict, False, {})
    # Several arguments can be specified to processors.selector: bounds_t, bounds_z, region, regions_param,
    # kwargs_select_t, kwargs_select_xy, kwargs_select_z.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_selector'.
    kwargs_selector = set_instance(kwargs_selector, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform remover"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # -- remove regional mean
        # select region
        tmp_kwargs = {k: v for k, v in kwargs.items() if k not in list(kwargs_selector.keys())}
        reg_array, reg_area = selector(input_array, data_var=data_var, data_var_area=data_var_area,
                                       input_area=input_area, **kwargs_selector, **tmp_kwargs)
        if reg_array is None:
            print("reg_array is None after remover.selector -> must break")
            break
        # spatial average
        tmp_kwargs = {k: v for k, v in kwargs.items() if k not in list(kwargs_averager.keys())}
        reg_array, reg_area = averager(reg_array, data_var=data_var, data_var_area=data_var_area,
                                       input_area=reg_area, **kwargs_averager, **tmp_kwargs)
        if reg_array is None:
            print("reg_array is None after remover.averager -> must break")
            break
        # remove region averaged array from original array (e.g., compute relative SST)
        da_array = xab.to_array(ds_array, data_var) - xab.to_array(reg_array["array"], data_var)
        ds_array = xab.set_array_in_place(ds_array, da_array, data_var=data_var)
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        text = ""
        for k in reg_array["metadata"]["description"].split(";; ")[-2:]:
            if isinstance(text, str) and len(text) > 0:
                text += ", "
            text += k.replace(str(k.split(") ")[0] + ") "), "")
        text += ", removed from array"
        metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


def seasonal_cycler(
        input_array: Union[dict[str, Union[dataset_wrapper, dict]], None],
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        kwargs_annual_cycle: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], None):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "kwargs_annual_cycle"]
    l2 = [input_array, input_area, data_var, data_var_area, kwargs_annual_cycle]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to xcdat_base.annual_cycle: frequency, keep_weights, reference_period,
    # season_config, skipna, weighted.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_annual_cycle'.
    kwargs_annual_cycle = set_instance(kwargs_annual_cycle, dict, False, {})
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        if wb.error_ds(ds_array, inspect__stack(), message="cannot perform seasonal_cycler"):
            break
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # perform interannual anomalies
        try:
            ds_array = xcb.annual_cycle(ds_array, data_var, **kwargs_annual_cycle)
        except Exception as err:
            print("WARNING cannot perform annual_cycle using xcdat")
            print(err)
            wb.log_debug(inspect__stack(), "WARNING cannot perform annual_cycle using xcdat\n" + str(err))
            break
        # adapt metadata
        description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
        metadata["description"] = basics.description_writer(description, "seasonal cycle computed")
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
            output_area = {"array": ds_area, "metadata": metadata_area}
    return output_array, output_area


def selector(
        input_array: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None],
        cf_dims: Union[list[Literal["T", "X", "XY", "Y", "Z"]], Literal["T", "X", "XY", "Y", "Z"], None] = None,
        bounds_t: Union[tuple[str, str], None] = None,
        bounds_z: Union[tuple[Union[float, int], Union[float, int]], None] = None,
        data_var: Union[str, None] = None,
        data_var_area: Union[str, None] = None,
        input_area: Union[dict[str, Union[array_wrapper, dataset_wrapper, dict]], None] = None,
        region: Union[str, None] = None,
        regions_param: dict[str, dict[str, Union[bool, str, tuple[Union[float, int]]]]] = None,
        kwargs_select_t: Union[dict, None] = None,
        kwargs_select_xy: Union[dict, None] = None,
        kwargs_select_z: Union[dict, None] = None,
        **kwargs) -> (Union[dict, None], Union[dict, None]):
    basics.log_info(inspect__stack(), "")
    l1 = ["input_array", "input_area", "data_var", "data_var_area", "bounds_t", "region", "regions_param",
          "bounds_z", "kwargs_select_t", "kwargs_select_xy", "kwargs_select_z"]
    l2 = [input_array, input_area, data_var, data_var_area, bounds_t, region, regions_param, bounds_z, kwargs_select_t,
          kwargs_select_xy, kwargs_select_z]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to wrapper_base.select_horizontal (used to select region):
    #      - mask_only: True to only mask data (for spatial average), False to select region
    #      - kwargs_sel (to select on regular grid): drop, method, tolerance (see xarray_base.select)
    #      - kwargs_where (to mask non-regular grid): drop (see xarray_base.where)
    # If desired they must be defined in a dictionary under the keyword 'kwargs_select_xy'.
    kwargs_select_xy = set_instance(kwargs_select_xy, dict, False, {})
    # One argument (a dictionary containing more arguments) can be specified to wrapper_base.select_depth (used to
    # select depth) and to wrapper_base.select_time (used to select time): kwargs_sel.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_select_t' or 'kwargs_select_z'.
    kwargs_select_t = set_instance(kwargs_select_t, dict, False, {})
    kwargs_select_z = set_instance(kwargs_select_z, dict, False, {})
    # Select depth, time and/or space
    if isinstance(cf_dims, str):
        cf_dims = [cf_dims]
    # fake loop to be able to break out if an error occurs
    output_array, output_area = None, None
    for _ in range(1):
        # check if time_bounds is defined
        do_t = False
        if isinstance(cf_dims, list) and "T" in cf_dims and isinstance(bounds_t, tuple) and len(bounds_t) == 2 and \
                all(isinstance(k, str) for k in bounds_t):
            do_t = True
        # check if region is defined
        do_xy = False
        if isinstance(cf_dims, list) and "XY" in cf_dims and isinstance(region, str) and \
                isinstance(regions_param, dict) and region in list(regions_param.keys()):
            do_xy = True
        # check if depth_bounds is defined
        do_z = False
        if isinstance(cf_dims, list) and "Z" in cf_dims and isinstance(bounds_z, tuple) and len(bounds_z) == 2 and \
                all(isinstance(k, (float, int)) for k in bounds_z):
            do_z = True
        if not do_t and not do_xy and not do_z:
            # WARNING: user calls selector without providing something to select
            message = "WARNING: user calls selector without providing something to select"
            l1 = ["bounds_t", "region", "regions_param", "bounds_z"]
            l2 = [bounds_t, region, regions_param, bounds_z]
            details = basics.log_details(l1, l2)
            wb.log_debug(inspect__stack(), message, adjust=5, details=details)
            break
        # get array and related input param
        ds_array = input_array["array"]
        metadata = copy__deepcopy(input_array["metadata"])
        ds_area, metadata_area = None, None
        if isinstance(input_area, dict) and "array" in list(input_area.keys()) and \
                "metadata" in list(input_area.keys()):
            ds_area, metadata_area = input_area["array"], input_area["metadata"]
        # select
        for cf_dim in cf_dims:
            if (cf_dim == "T" and not do_t) or (cf_dim == "XY" and not do_xy) or (cf_dim == "Z" and not do_z):
                continue
            bounds_xy = {}
            if cf_dim == "T":
                ds_array = wb.select_time(ds_array, data_var=data_var, time_bounds=bounds_t, **kwargs_select_t)
            elif cf_dim == "XY":
                if "latitude" in list(regions_param[region].keys()) and \
                        isinstance(regions_param[region]["latitude"], (list, tuple)) and \
                        len(regions_param[region]["latitude"]) >= 2:
                    bounds_xy["Y"] = regions_param[region]["latitude"]
                if "longitude" in list(regions_param[region].keys()) and \
                        isinstance(regions_param[region]["longitude"], (list, tuple)) and \
                        len(regions_param[region]["longitude"]) >= 2:
                    bounds_xy["X"] = regions_param[region]["longitude"]
                ds_array = wb.select_horizontal(ds_array, data_var=data_var, horizontal_bounds=bounds_xy,
                                                **kwargs_select_xy)
                if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
                    ds_area = wb.select_horizontal(ds_area, data_var=data_var_area, horizontal_bounds=bounds_xy,
                                                   **kwargs_select_xy)
            else:
                ds_array = wb.select_depth(ds_array, data_var=data_var, depth_bounds=bounds_z, **kwargs_select_z)
            if ds_array is None:
                break
            # redo bounds (operations on ds changes bounds and xcdat doesn't like that)
            # ds_array = wb.redo_bounds(ds_array, ["T", "X", "Y", "Z"])
            # adapt metadata
            description = copy__deepcopy(metadata["description"]) if "description" in list(metadata.keys()) else ""
            if cf_dim == "T":
                text = "time selected " + str(bounds_t)
            elif cf_dim == "XY":
                n_sho, n_lon = copy__deepcopy(region), ""
                if "short_name" in list(regions_param[region].keys()):
                    n_sho = regions_param[region]["short_name"]
                if "long_name" in list(regions_param[region].keys()):
                    n_lon = " (" + str(regions_param[region]["long_name"]) + ")"
                text = "selected in " + str(n_sho) + str(n_lon)
                if "Y" in list(bounds_xy.keys()) and "X" in list(bounds_xy.keys()) and len(bounds_xy["Y"]) == 2 and \
                        len(bounds_xy["X"]) == 2:
                    text += " " + str(basics.write_coordinates(*bounds_xy["Y"], *bounds_xy["X"]))
            else:
                text = "depth selected (" + str(bounds_z) + ")"
            metadata["description"] = basics.description_writer(description, text)
        # prepare output
        if isinstance(ds_array, (array_wrapper, dataset_wrapper)):
            output_array = {"array": ds_array, "metadata": metadata}
        if isinstance(ds_area, (array_wrapper, dataset_wrapper)):
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
    # print("processors", list(processors.keys()))
    # print("list_processors", list_processors)
    # fx variables
    fx = dict((k1, d1) for k1, d1 in input_dataset.items() if k1 in ["areacell", "areacella", "areacello", "landmask"])
    # loop on variables to process
    dict_o = {}
    for k1 in list(processors.keys()):
        # print(k1)
        # k1 should be the name of the new output variable and d1 is a dictionary.
        # E.g., d1 = {
        #     "region": "region_name",
        #     "variable": "variable_name",
        #     "to_do": {
        #         "1__masker": {arguments},
        #         "2__selector": {arguments},
        #         "3__???": {arguments},
        #         ...}}
        # input and output variable names
        variable_i = processors[k1]["variable"]
        region_i = processors[k1]["region"]
        variable_o = copy__deepcopy(k1)
        # check if given variable is available
        if variable_i not in list(input_dataset.keys()) or not isinstance(input_param, dict) or \
                variable_i not in list(input_param.keys()):
            # WARNING: variable must be defined
            details = {"variable": str(variable_i),
                       "in input_dataset": str(variable_i in list(input_dataset.keys())),
                       "input_param.type": str(type(input_param))}
            if isinstance(input_param, dict):
                details["input_param.keys"] = ", ".join(list(input_param.keys()))
            wb.log_debug(inspect__stack(), "WARNING: variable must be defined", adjust=5, details=details)
            break
        # get param given variable as well as area and mask names related to given variable
        param = input_param[variable_i]
        data_var_area = param["area"] if "area" in list(param.keys()) else None
        data_var_mask = param["mask"] if "mask" in list(param.keys()) else None
        # get variable, area and mask dictionaries (i.e., {"array": xarray.Dataset, "metadata": {}})
        dict_array = input_dataset[variable_i]
        dict_area, dict_mask = None, None
        for n1, n2 in zip(["area", "mask"], [data_var_area, data_var_mask]):
            if isinstance(n2, str) and (n2 not in list(input_dataset.keys()) or not isinstance(input_param, dict) or
                                        n2 not in list(input_param.keys())):
                # WARNING: given variable must be defined
                details = {"variable": str(n2),
                           "in input_dataset": str(n2 in list(input_dataset.keys())),
                           "input_param.type": str(type(input_param))}
                if isinstance(input_param, dict):
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
        # loop on processors to apply to given variable
        # print("processors", list(processors[k1]["to_do"].keys()))
        for k2 in list(processors[k1]["to_do"].keys()):
            # print(k2)
            process = k2.split("__")[-1]
            if process in list_processors:
                # call processor
                # print(k2, "call processor")
                if process == "masker" and dict_mask is None:
                    wb.log_debug(inspect__stack(), "WARNING " + str(variable_i) + " not masked as mask not provided")
                    continue
                local_kwargs = processors[k1]["to_do"][k2]
                # print("local_kwargs", list(local_kwargs.keys()))
                dict_array, dict_area = dict_processors[process](
                    dict_array, data_var=variable_i, data_var_area=data_var_area, data_var_mask=data_var_mask,
                    input_area=dict_area, input_mask=dict_mask, region=region_i, **local_kwargs, **kwargs)
                if dict_array is None:
                    print("dict_array is None -> must break")
                    break
        if dict_array is None:
            break
        # change data_var name in dict_array["array"]
        # print("loop", type(dict_array["array"]))
        # for k in list(dict_array["array"].keys()):
        #     print(str(k).rjust(20), dict_array["array"][k].shape, dict_array["array"][k].dims)
        if isinstance(variable_o, str) and len(variable_o) > 0 and variable_i != variable_o:
            dict_array["array"] = wb.rename_variable(dict_array["array"], variable_i, variable_o)
        # store dict_array["array"] in output dictionary
        dict_o[variable_o] = dict_array
    if len(dict_o.keys()) != len(processors.keys()):
        dict_o = None
    return dict_o


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: reader (to read netCDF files) & saver (to save in a netCDF file)
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
    if isinstance(variables, str):
        variables = [variables]
    # read variable (as in file)
    dict_t = {}
    for kk in variables + ["areacella", "areacello", "landmask"]:
        if kk not in list(input_param.keys()) or kk not in list(variables_param.keys()):
            continue
        # loop on netCDF files-variables to read all required dataf
        for ff, nn in zip(input_param[kk]["file_name"], input_param[kk]["variable"]):
            # remove time dimension from add_bounds and decode_times if not defined in netCDF ('fx' frequency)
            ab, dt = copy__deepcopy(add_bounds), copy__deepcopy(decode_times)
            if kk in ["areacella", "areacello", "landmask"]:
                if "T" in ab:
                    while "T" in ab:
                        ab.remove("T")
                if dt:
                    dt = False
            # try to open_dataset and save Dataset in a dictionary using netCDF variables (names) as keys
            kwargs_open_dataset = {"add_bounds": ab, "decode_times": dt, **kwargs_reader}
            ds = wb.open_dataset(ff, data_var=nn, package="xcdat", kwargs_open_dataset=kwargs_open_dataset)
            if isinstance(ds, dataset_wrapper):
                dict_t[nn] = ds
            # try:
            #     dict_t[nn] = xcb.open_dataset(
            #         ff, add_bounds=ab, data_var=nn, decode_times=dt, **kwargs_reader)
            # except Exception as err:
            #     message = "can't read (" + str(nn) + ") " + str(ff) + "\n" + str(err)
            #     basics.log_error(inspect__stack(), message)
            #     # WARNING: cannot read variable or file
            #     path = "/".join(ff.split("/")[:-1])
            #     files = sorted(list(glob__iglob(ff)), key=lambda s: s.lower())
            #     files_string = ""
            #     if len(files) > 0:
            #         for k in files:
            #             files_string += "\n" + str().ljust(5) + str(k)
            #     else:
            #         files_string = "no file matches this file pattern"
            #     details = {
            #         "directory": str(path),
            #         "isdir": str(os__path__isdir(path)),
            #         "file": str(ff),
            #         "isfile": str(os__path__isfile(ff)),
            #         "list": str(files_string)}
            #     wb.log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
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
                not isinstance(input_param[k1], dict) or "mask" not in list(input_param[k1].keys()):
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


def saver(
        processed_data,
        dataset: str = None,
        experiment: str = None,
        extension: str = None,
        filename_netcdf: str = None,
        member: str = None,
        path: str = None,
        project: str = None,
        recipe: str = None,
        recipe_add: str = None,
        kwargs_generate_filename: dict = None,
        kwargs_merge: dict = None,
        kwargs_to_netcdf: dict = None,
        **kwargs):
    basics.log_info(inspect__stack(), "")
    l1 = ["processed_data", "filename_netcdf", "path", "project", "dataset", "experiment", "member", "recipe",
          "kwargs_merge", "kwargs_to_netcdf"]
    l2 = [processed_data, filename_netcdf, path, project, dataset, experiment, member, recipe, kwargs_merge,
          kwargs_to_netcdf]
    details = basics.log_details(l1, l2)
    wb.log_debug(inspect__stack(), "input", adjust=5, details=details)
    # Several arguments can be specified to xcdat_base.merge: combine_attrs, compat, fill_value, join.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_merge'.
    kwargs_merge = set_instance(kwargs_merge, dict, False, {})
    # Several arguments can be specified to xcdat_base.to_netcdf: file_format, list_of_variables, mode.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_to_netcdf'.
    kwargs_to_netcdf = set_instance(kwargs_merge, dict, False, {})
    # Extra keywords can be provided to basics.generate_filename to generate the output filename.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_generate_filename'.
    kwargs_generate_filename = set_instance(kwargs_generate_filename, dict, False, {})
    # --- Step 1: check dataset and assign variable attributes
    dict_dataset, dict_dim = {}, {}
    for k1, d1 in processed_data.items():
        # k1 should be the name of the new output variable and d1 is a dictionary
        ds_array, metadata = None, None
        # get dataset
        if isinstance(d1, dict) and "array" in list(d1.keys()):
            ds_array = d1["array"]
        if not isinstance(ds_array, dataset_wrapper):
            continue
        # get variable metadata and set to variable
        if isinstance(d1, dict) and "metadata" in list(d1.keys()):
            metadata = basics.flatten_dict(d1["metadata"])
            xab.drop_given_attributes(ds_array, xab.get_attributes_keys(ds_array, data_var=k1), data_var=k1)
            xab.set_attributes_variable(ds_array, k1, **metadata)
        # put ds_array in temporary dictionary
        dict_dataset[k1] = ds_array
    # --- Step 2: merge
    # remove bounds (not easy to merge them)
    list_dat = [xab.drop_dataset_keys(d1, [k2 for k2 in list(d1.keys()) if "bound_" in k2 or "bounds_" in k2 or
                                           "bnd_" in k2 or "bnds_" in k2 or "_bound" in k2 or "_bnd" in k2])
                for k1, d1 in dict_dataset.items()]
    # merge datasets
    ds_o = xab.merge(list_dat, **kwargs_merge)
    # ds_o = xab.merge([k for k in dict_dataset.values() if not("_bnd" in k or "_bound" in k)], **kwargs_merge)
    # merge global attributes
    metadata = basics.merge_metadata({k: xab.get_attributes(v) for k, v in dict_dataset.items()})
    xab.drop_given_attributes(ds_o, xab.get_attributes_keys(ds_o))
    xab.set_attributes_global(ds_o, **metadata)
    # --- Step 3: save
    # output file name
    fo = basics.generate_filename(dataset=dataset, experiment=experiment, extension=extension,
                                  filename_netcdf=filename_netcdf, member=member, path=path, project=project,
                                  recipe=recipe, recipe_add=recipe_add, **kwargs_generate_filename)
    # save netCDF
    xab.to_netcdf(ds_o, fo, **kwargs_to_netcdf)
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions: comparator (to compute metrics)
# ---------------------------------------------------------------------------------------------------------------------#
def comparator(
        input_dataset: str,
        input_reference: dict[str, str],
        metric_method: str,
        data_var: str = None,
        metric_dictionary: dict = None,
        metric_dictionary_keys: tuple[str] = None,
        **kwargs) -> dict:
    metric_dictionary = set_instance(metric_dictionary, dict, False, {})
    # fake loop to be able to break out
    known_methods = {
        "difference": "difference",
        "difference_relative": "relative difference",
        "relative_difference_relative_absolute": "absolute relative difference",
        "correlation": "correlation",
        "rmse": "rmse"}
    ds = None
    dt = {}
    for _ in [0]:
        if metric_method not in list(known_methods.keys()):
            message = "can't compare (" + str(metric_method) + "): unknow metric method"
            basics.log_error(inspect__stack(), message)
            # WARNING: cannot compute metric
            details = {
                "file": str(input_dataset),
                "data_var": str(data_var),
                "metric_method": str(metric_method),
                "known method(s)": ", ".join(list(known_methods.keys()))}
            wb.log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
        # --- Step 1: open input dataset (to compare to references)
        # open
        ds = wb.open_dataset(input_dataset, package="xarray")
        if not isinstance(ds, dataset_wrapper):
            break
        # read data_var
        if not (isinstance(data_var, str) and data_var in list(ds.keys())):
            data_var = [k for k in list(ds.keys()) if str(k).split("_")[0] == "d1"]
            if len(data_var) != 1:
                message = "can't compare (" + str(data_var) + "): cannot find the diagnostic"
                basics.log_error(inspect__stack(), message)
                # WARNING: cannot compute metric
                details = {
                    "file": str(input_dataset),
                    "data_var": str(data_var),
                    "ds.keys": str(list(ds.keys()))}
                wb.log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
                break
            data_var = data_var[0]
        da = xab.to_array(ds, data_var)
        # --- Step 2: open input dataset (to compare to references)
        for n, f in input_reference.items():
            # open
            ds_r = wb.open_dataset(f, package="xarray")
            if not isinstance(ds, dataset_wrapper):
                continue
            # read data_var
            da_r = xab.to_array(ds_r, data_var)
            # -- Step 3: compute the metric
            if isinstance(da, array_wrapper) and len(xab.get_array_shape(da)) == 0 and xab.get_array_size(da) == 1 and \
                    isinstance(da_r, array_wrapper) and len(xab.get_array_shape(da_r)) == 0 and \
                    xab.get_array_size(da_r) == 1 and \
                    metric_method in ["difference", "difference_relative", "relative_difference_relative_absolute"]:
                mv = da - da_r
                if metric_method == "difference_relative":
                    mv = mv * 100 / da_r
                if metric_method == "relative_difference_relative_absolute":
                    mv = abs(mv)
                mv = float(mv.values)
            elif isinstance(da, array_wrapper) and len(xab.get_array_shape(da)) > 0 and \
                    isinstance(da_r, array_wrapper) and len(xab.get_array_shape(da_r)) > 0 and \
                    metric_method in ["correlation", "rmse"]:
                dim_t, dim_x, dim_y = wb.get_dim_time(da), wb.get_dim_longitude(da), wb.get_dim_latitude(da)
                dims: list[Literal["T", "X", "Y"]] = ["T", "X", "Y"]
                dims = [k1 for k1, k2 in zip(dims, [dim_t, dim_x, dim_y]) if basics.is_dim(k2)]
                weights = wb.compute_weights(da, cf_dim=dims)
                if metric_method == "correlation":
                    mv = xab.correlation(da, da_r, dim=dims, weights=weights)
                else:
                    mv = xab.mean((da - da_r) ** 2, dim=dims, skipna=True, weights=weights)
                mv = float(mv.values)
            else:
                message = "can't compare: type error"
                basics.log_error(inspect__stack(), message)
                # WARNING: cannot compute metric
                details = {
                    "file": str(input_dataset),
                    "data_var": str(data_var),
                    "metric_method": str(metric_method)}
                for k1, k2 in zip(["da_mod", "da_" + str(n)], [da, da_r]):
                    details[str(k1) + ".type"] = str(type(k2))
                    if isinstance(da, array_wrapper):
                        details[str(k1) + ".shape"] = str(xab.get_array_shape(da))
                        details[str(k1) + ".dims"] = str(xab.get_dim_keys(da))
                        details[str(k1) + ".size"] = str(xab.get_array_size(da))
                wb.log_debug(inspect__stack(), "WARNING: " + str(message), adjust=5, details=details)
                continue
            # save metric in a dictionary
            dt[n] = mv
        if len(list(dt.keys())) == 0:
            continue
        # --- Step 4: add to metadata and resave netCDF
        # get attributes
        attrs = xab.get_attributes(da)
        # adapt metadata
        description = copy__deepcopy(attrs["description"]) if "description" in list(attrs.keys()) else ""
        attrs["description"] = basics.description_writer(description, str(known_methods[metric_method]) + " computed")
        attrs["metric_method"] = copy__deepcopy(metric_method)
        for k, v in dt.items():
            attrs["metric_value__" + str(k)] = copy__deepcopy(v)
        attrs = dict(sorted(attrs.items()))
        # drop attrs and set attrs
        da = xab.drop_attrs(da)
        xab.set_attributes_variable(da, **attrs)
        # place in input dataset
        ds = xab.set_array_in_place(ds, da, data_var=data_var)
        # save netCDF
        os__remove(input_dataset)
        xab.to_netcdf(ds, input_dataset, mode="a")
        # fill output dictionary
        if isinstance(metric_dictionary_keys, tuple) and len(metric_dictionary_keys) > 0:
            metric_dictionary = basics.set_nested_value(metric_dictionary, dt, *metric_dictionary_keys)
        else:
            metric_dictionary = copy__deepcopy(dt)
    return metric_dictionary
# ---------------------------------------------------------------------------------------------------------------------#
