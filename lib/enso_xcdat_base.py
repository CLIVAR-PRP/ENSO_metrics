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


def fct_averager(ds: xr.Dataset, data_var: str, weights: Union[str, xr.DataArray]) -> xr.Dataset:
    return ds.spatial.average(data_var, weights=weights)


def fct_get_axis_key(da: Union[xr.DataArray, xr.Dataset], axis: Literal["X", "Y", "T", "Z"]) -> str:
    return xc.axis.get_dim_keys(da, axis=axis)


def fct_get_axis_list(ds: xr.Dataset) -> list:
    return list(ds.coords.keys())


def fct_get_longitude(ds: xr.Dataset) -> xr.DataArray:
    key_lon = fct_get_axis_key(ds, "X")
    return ds[key_lon]


def fct_get_latitude(ds: xr.Dataset) -> xr.DataArray:
    key_lat = fct_get_axis_key(ds, "Y")
    return ds[key_lat]


def fct_interannual_anomalies(ds: xr.Dataset, data_var: str,
                                freq: Literal["year", "season", "month", "day", "hour"] = "month") -> xr.Dataset:
    return ds.temporal.departures(data_var, freq=freq, weighted=True)


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


def fct_uniform_grid(lat_start: float, lat_end: float, lat_res: float, lon_start: float, lon_end: float,
                       lon_res: float) -> xr.Dataset:
    return xc.create_uniform_grid(lat_start, lat_end, lat_res, lon_start, lon_end, lon_res)


def fct_xarray_to_numpy(arr: xr.DataArray) -> np.ndarray:
    return arr.to_numpy()
