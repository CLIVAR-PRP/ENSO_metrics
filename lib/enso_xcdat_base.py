# -*- coding:UTF-8 -*-

# basic python package
from copy import deepcopy
from typing import Literal, Union, Tuple
# numpy
import numpy as np
# xarray
import xarray as xr
# xCDAT
import xcdat as xc


# ---------------------------------------------------------------------------------------------------------------------#
# xarray / xcdat base functions
# ---------------------------------------------------------------------------------------------------------------------#
def fct_array_ones(da: xr.DataArray, array_name: str = "my_name") -> xr.DataArray:
    da_ones = xr.ones_like(da)
    da_ones.name = array_name
    return da_ones


def fct_array_zeros(da: xr.DataArray, array_name: str = "my_name") -> xr.DataArray:
    da_zeros = xr.zeros_like(da)
    da_zeros.name = array_name
    return da_zeros


def fct_averager(ds: xr.Dataset, data_var: str, axis: list[str], weights: Union[str, xr.DataArray]) -> xr.Dataset:
    return ds.spatial.average(data_var, axis=axis, weights=weights)


def fct_get_axis_key(da: Union[xr.DataArray, xr.Dataset], axis: Literal["X", "Y", "T", "Z"]) -> str:
    return xc.axis.get_dim_keys(da, axis=axis)


def fct_get_axis_list(da: Union[xr.DataArray, xr.Dataset]) -> list:
    return list(da.coords.keys())


def fct_get_dimensions(da: xr.DataArray) -> list:
    return list(da.dims)


def fct_get_latitude(da: Union[xr.DataArray, xr.Dataset]) -> xr.DataArray:
    key_lat = fct_get_axis_key(da, "Y")
    return da[key_lat]


def fct_get_longitude(da: Union[xr.DataArray, xr.Dataset]) -> xr.DataArray:
    key_lon = fct_get_axis_key(da, "X")
    return da[key_lon]


def fct_get_time_bounds(da: Union[xr.DataArray, xr.Dataset]) -> list[str]:
    key_time = fct_get_axis_key(da, "T")
    time1 = da[key_time][0].values.item().isoformat().split("T")[0]
    time2 = da[key_time][-1].values.item().isoformat().split("T")[0]
    return [time1, time2]


def fct_interannual_anomalies(ds: xr.Dataset, data_var: str,
                              freq: Literal["year", "season", "month", "day", "hour"] = "month") -> xr.Dataset:
    return ds.temporal.departures(data_var, freq=freq, weighted=True)


def fct_annual_cycle(ds: xr.Dataset, data_var: str, 
                     freq: Literal["season", "month", "day"] = "month") -> xr.Dataset:
    return ds.temporal.climatology(data_var, freq=freq, weighted=True)


def fct_max(da: xr.DataArray) -> float:
    return float(da.max().values)


def fct_min(da: xr.DataArray) -> float:
    return float(da.min().values)


def fct_numpy_to_xarray(arr: np.ndarray, axes: tuple, axes_coordinates: dict) -> xr.DataArray:
    return xr.DataArray(arr, dims=axes, coords=axes_coordinates)


def fct_open_dataset(data_path: str) -> xr.Dataset:
    return xc.open_mfdataset(data_path)


def fct_set_auto_bounds(ds: xr.Dataset) -> xr.Dataset:
    return ds.bounds.add_missing_bounds()


def fct_sum(da: Union[xr.DataArray, xr.Dataset], axis: Union[str, tuple]) -> Union[xr.DataArray, xr.Dataset]:
    return da.sum(dim=axis)


def fct_standard_deviation(da: Union[xr.DataArray, xr.Dataset], axis: Union[list, str, tuple],
                             ddof: int) -> Union[xr.DataArray, xr.Dataset]:
    return da.std(dim=axis, ddof=ddof)


def fct_transpose(da: xr.DataArray, axis_list: list[str]):
    return da.transpose(*axis_list)


def fct_uniform_grid(lat_start: float, lat_end: float, lat_res: float, lon_start: float, lon_end: float,
                       lon_res: float) -> xr.Dataset:
    return xc.create_uniform_grid(lat_start, lat_end, lat_res, lon_start, lon_end, lon_res)


def fct_xarray_to_numpy(arr: xr.DataArray) -> np.ndarray:
    return arr.to_numpy()


def fct_where_xarray(cond: xr.DataArray, x, y) -> xr.DataArray:
    """_summary_
    See https://docs.xarray.dev/en/stable/generated/xarray.where.html for details

    Parameters
    ----------
    cond : xr.DataArray
        When True, return values from x, otherwise returns values from y.
    x : _type_
        values to choose from where cond is True
    y : _type_
        values to choose from where cond is False

    Returns
    -------
    xr.DataArray
    """
    return xr.where(cond, x, y)


def fct_where_dataarray(arr: xr.DataArray, cond: xr.DataArray) -> xr.DataArray:
    """_summary_
    See https://docs.xarray.dev/en/stable/generated/xarray.DataArray.where.html for details

    Parameters
    ----------
    cond : xr.DataArray
        Locations at which to preserve this object’s values. dtype must be bool. If a callable, the callable is passed this object, and the result is used as the value for cond.

    Returns
    -------
    xr.DataArray
    """
    return arr.where(cond)


def fct_horizontal_regrid(ds: xr.Dataset, data_var: str, target_grid: xr.DataArray, regrid_tool: str="xesmf", regrid_method :str="bilinear") -> xr.Dataset:
    return ds.regridder.horizontal(data_var, target_grid, tool=regrid_tool, method=regrid_method, unmapped_to_nan=True)


def fct_attrs_update_global(ds: xr.Dataset, global_attributes: dict):
    ds.attrs.update(global_attributes)

    
def fct_attrs_update_variable(ds: xr.Dataset, data_var: str, variable_attributes: dict):
    ds[data_var].attrs.update(variable_attributes)
    

def fct_dataset_to_netcdf(ds: xr.Dataset, filename: str, list_of_variables=list[str]):
    if list_of_variables is None:
        ds.to_netcdf(filename)
    else:
        ds[list_of_variables].to_netcdf(filename)
        

def fct_get_attributes_list(arr: xr.DataArray) -> list[str]:
    return list(arr.attrs.keys())


def fct_get_attribute(arr: xr.DataArray, attribute_name: str) -> str:
    return arr.attrs[attribute_name]

