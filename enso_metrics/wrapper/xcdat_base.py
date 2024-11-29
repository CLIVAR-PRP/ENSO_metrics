# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# N-dimensional array with labeled coordinates and dimensions are needed for the package.
# xCDAT is an extension of xarray for climate data analysis on structured grids. It serves as a modern successor to the
# Community Data Analysis Tools (CDAT) library.
# This file regroups used xCDAT functions
# https://xcdat.readthedocs.io/en/latest/
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from typing import Callable, Literal, Union
# xarray
from xarray import DataArray as xarray__DataArray
from xarray import Dataset as xarray__Dataset
# xCDAT
import xcdat
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def annual_cycle(
        ds: xarray__Dataset,
        data_var: str,
        frequency: Literal["day", "month", "season"] = "month",
        keep_weights: bool = False,
        reference_period: Union[tuple[str, str], None] = None,
        weighted: bool = True) -> xarray__Dataset:
    """
    Returns a Dataset with the climatology of a data variable.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.temporal.climatology.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param frequency: {"day", "month", "season"}, optional
        The time frequency to group by; e.g., frequency = "month".
            - “day”: groups by (month, day) for the daily cycle climatology. If the CF calendar type is "gregorian",
                     "proleptic_gregorian", or "standard", leap days (if present) are dropped to avoid inconsistencies
                     when calculating climatologies.
            - “month”: groups by month for the annual cycle climatology.
            - “season”: groups by season for the seasonal cycle climatology.
        Default is "month"
    :param keep_weights: bool, optional
        If calculating averages using weights, keep the weights in the final dataset output; e.g., keep_weights = False.
        Default is False
    :param reference_period: Tuple[str, str], None, optional
        The climatological reference period, which is a subset of the entire time series. This parameter must be
        tuple of strings in the format ‘yyyy-mm-dd’; e.g., reference_period = ('1850-01-01', '1899-12-31').
        If no value is provided, the climatological reference period will be the full period covered by the dataset.
        Default is None
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Default is True
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the climatology of given data variable.
    """
    tmp_kwargs = {"keep_weights": keep_weights, "weighted": weighted}
    if reference_period is not None:
        tmp_kwargs["reference_period"] = reference_period
    return ds.temporal.climatology(data_var, frequency, **tmp_kwargs)


def average_spatial(
        ds: xarray__Dataset,
        data_var: str,
        cf_dim: list[Literal["X", "Y"]] = None,
        keep_weights: bool = False,
        weights: Union[str, xarray__DataArray] = "generate") -> xarray__Dataset:
    """
    Return a Dataset with the average of a data variable and the given spatial dimension(s) removed.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.spatial.average.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param cf_dim: list[{"X", "Y"}]
        List of axis dimensions to average over, valid axis keys include 'X' and 'Y'; e.g., cf_axis = ["X", "Y"].
        Default is None (i.e., ['X', 'Y'])
    :param keep_weights: bool, optional
        If calculating averages using weights, keep the weights in the final dataset output; e.g., keep_weights = False.
        Default is False
    :param weights: Union["generate", xr.DataArray]
        If "generate", then weights are generated, otherwise, DataArray must contain the regional weights used for
        weighted averaging.
        Default is "generate"
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the average of a data variable and the given spatial dimension(s) removed.
    """
    if isinstance(cf_dim, list) is False:
        cf_dim = ["X", "Y"]
    # spatial average
    return ds.spatial.average(data_var, axis=cf_dim, keep_weights=keep_weights, weights=weights)


def average_temporal(
        ds: xarray__Dataset,
        data_var: str,
        keep_weights: bool = False,
        weighted: bool = True) -> xarray__Dataset:
    """
    Return a Dataset with the average of a data variable and the time dimension removed.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.temporal.average.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param keep_weights: bool, optional
        If calculating averages using weights, keep the weights in the final dataset output; e.g., keep_weights = False.
        Default is False
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Default is True
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the average of a data variable and the time dimension removed.
    """
    # temporal average
    return ds.temporal.average(data_var, keep_weights=keep_weights, weighted=weighted)


def create_uniform_grid(lat_start: float, lat_stop: float, lat_delta: float, lon_start: float, lon_stop: float,
                        lon_delta: float) -> xarray__Dataset:
    """
    Create a uniform rectilinear grid and sets appropriate the attributes for the lat/lon axis.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.create_uniform_grid.html
    
    Inputs:
    -------
    :param lat_start: float
        First latitude
    :param lat_stop: float
        Last latitude
    :param lat_delta: float
        Difference between two points of axis
    :param lon_start: float
        First longitude
    :param lon_stop: float
        Last longitude
    :param lon_delta: float
        Difference between two points of axis
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with uniform lat/lon grid.
    """
    return xcdat.create_uniform_grid(lat_start, lat_stop, lat_delta, lon_start, lon_stop, lon_delta)


def get_axis_key(
        ds: Union[xarray__DataArray, xarray__Dataset],
        cf_dim: Literal["X", "Y", "T", "Z"]) -> Union[str, list[str]]:
    """
    Gets the dimension key(s) for an axis.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.get_dim_keys.html
    
    Inputs:
    -------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param cf_dim: {"X", "Y", "T", "Z"}
        The CF axis (dimension) key
    
    Output:
    -------
    :return: str or list[str]
        The dimension string or a list of dimensions strings for an axis.
    """
    return xcdat.axis.get_dim_keys(ds, axis=cf_dim)


def interannual_anomalies(
        ds: xarray__Dataset,
        data_var: str,
        frequency: Literal["day", "month", "season"] = "month",
        keep_weights: bool = False,
        reference_period: Union[tuple[str, str], None] = None,
        weighted: bool = True) -> xarray__Dataset:
    """
    Returns a Dataset with the climatological departures (anomalies) for a data variable.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.temporal.departures.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param frequency: {"day", "month", "season"}, optional
        The frequency of time to group by; e.g., frequency = "month".
            - “day”: groups by (month, day) for the daily cycle departures. If the CF calendar type is "gregorian",
                     "proleptic_gregorian", or "standard", leap days (if present) are dropped to avoid inconsistencies
                     when calculating climatologies.
            - “month”: groups by month for the annual cycle departures.
            - “season”: groups by season for the seasonal cycle departures.
        Default is "month"
    :param keep_weights: bool, optional
        If calculating averages using weights, keep the weights in the final dataset output; e.g., keep_weights = False.
        Default is False
    :param reference_period: Tuple[str, str], None, optional
        The climatological reference period, which is a subset of the entire time series. This parameter must be
        tuple of strings in the format ‘yyyy-mm-dd’; e.g., reference_period = ('1850-01-01', '1899-12-31').
        If no value is provided, the climatological reference period will be the full period covered by the dataset.
        Default is None
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Default is True
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the climatological departures (anomalies) for a data variable.
    """
    tmp_kwargs = {"keep_weights": keep_weights, "weighted": weighted}
    if reference_period is not None:
        tmp_kwargs["reference_period"] = reference_period
    return ds.temporal.departures(data_var, frequency, **tmp_kwargs)


def open_dataset(
        paths: Union[str, list[str]],
        add_bounds: Union[list[Literal["T", "X", "Y", "Z"]], bool, None] = None,
        center_times: bool = False,
        data_var: str = None,
        data_vars: Union[Literal["minimal", "different", "all"], list[str]] = "minimal",
        decode_times: bool = True,
        lon_orient: Union[tuple[float, float], tuple[int, int], None] = None,
        preprocess: Union[Callable, None] = None,
        **kwargs) -> xarray__Dataset:
    """
    Open multiple files as a single dataset (wraps xarray.open_mfdataset() with post-processing options).
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.open_mfdataset.html
    
    Inputs:
    -------
    :param paths: str or list[str]
        Paths to dataset files:
            - Directory path (e.g., "path/to/files"), which is converted to a string glob of *.nc files
            - String glob (e.g., "path/to/files/*.nc"), which is expanded to a 1-dimensional list of file paths
            - File path to dataset (e.g., "path/to/files/file1.nc")
            - List of file paths (e.g., ["path/to/files/file1.nc", ...])
        If multiple files, concatenation along the time dimension recommended
    :param add_bounds: list[{"T", "X", "Y", "Z"}] or None or bool
        List of CF axes to try to add bounds for (if missing); e.g., add_bounds = ["X", "Y"].
        This parameter calls xarray.Dataset.bounds.add_missing_bounds().
        Default is None (i.e., ["T", "X", "Y"])
    :param center_times: bool, optional
        If True, attempt to center time coordinates using the midpoint between its upper and lower bounds. Otherwise,
        use the provided time coordinates; e.g., center_times = False.
        Default is False
    :param data_var: str, optional
        The key of the data variable to keep in the Dataset; e.g., data_var = "ts".
        Default is None
    :param data_vars: {"minimal", "different", "all" or list of str}, optional
        These data variables will be concatenated together; e.g., data_vars = "minimal".
            - “minimal”: Only data variables in which the dimension already appears are included.
            - “different”: Data variables which are not equal (ignoring attributes) across all datasets are also
                           concatenated (as well as all for which dimension already appears). Beware: this option may
                           load the data payload of data variables into memory if they are not already loaded.
            - “all”: All data variables will be concatenated.
            - list of str: The listed data variables will be concatenated, in addition to the “minimal” data variables.
        Default is "minimal"
    :param decode_times: bool, optional
        If True, attempt to decode times encoded in the standard NetCDF datetime format into cftime.datetime objects.
        Otherwise, leave them encoded as numbers. This keyword may not be supported by all the backends;
        e.g., decode_times = True.
        Default is True
    :param lon_orient: Tuple[float, float] or None, optional
        Orientation to use for the Dataset’s longitude axis (if it exists); e.g., lon_orient = (-180, 180).
            - None: use the current orientation (if the longitude axis exists).
            - (-180, 180): represents [-180, 180] in math notation.
            - (0, 360): represents [0, 360] in math notation.
        Default is None
    :param preprocess: Callable, optional
        If provided, call this function on each dataset prior to concatenation. You can find the file-name from which
        each dataset was loaded in ds.encoding["source"].
        Default is None
    **kwargs – Additional keyword arguments passed on to xarray.open_mfdataset.

    Output:
    -------
    :return: xarray.Dataset
        Newly created dataset.
    """
    if add_bounds is None:
        add_bounds = ["T", "X", "Y"]
    tmp_kwargs = {"add_bounds": add_bounds, "center_times": center_times, "data_var": data_var, "data_vars": data_vars,
                  "decode_times": decode_times, "lon_orient": lon_orient, "preprocess": preprocess, **kwargs}
    return xcdat.open_mfdataset(paths, **tmp_kwargs)


def regrid_horizontal(
        ds: xarray__Dataset,
        data_var: str,
        output_grid: Union[xarray__DataArray, xarray__Dataset],
        method: Literal[
            "bilinear", "conservative", "conservative_normed", "patch", "nearest_s2d", "nearest_d2s"] = "conservative",
        tool: Literal["regrid2", "xesmf"] = "xesmf",
        unmapped_to_nan: bool = True,
        **kwargs) -> xarray__Dataset:
    """
    Regrid data_var to output_grid.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.regridder.horizontal.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param output_grid: xarray.DataArray or xarray.Dataset
        Grid to transform data_var to
    :param tool: {"regrid2", "xesmf"}, optional
        Name of the tool to use; e.g., tool = "regrid2"
    :param method: {"bilinear", "conservative", "conservative_normed", "patch", "nearest_s2d", "nearest_d2s"}, optional
        Regridding method to apply; e.g., method = "conservative".
        If tool is "regrid2": "conservative".
        If tool is "xesmf": "bilinear", "conservative", "conservative_normed", "patch", "nearest_s2d", "nearest_d2s".
        Default is "conservative"
    :param unmapped_to_nan: bool, optional
        Sets values of unmapped points to numpy.nan instead of 0; e.g., unmapped_to_nan = True.
        Default is True
    **kwargs – Additional keyword arguments passed on to the regridder.
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the data_var transformed to the output_grid.
    """
    tmp_kwargs = {"method": method, "tool": tool, "unmapped_to_nan": unmapped_to_nan, **kwargs}
    return ds.regridder.horizontal(data_var, output_grid, **tmp_kwargs)


def regrid_vertical(
        ds: xarray__Dataset,
        data_var: str,
        output_grid: Union[xarray__DataArray, xarray__Dataset],
        method: Literal["linear", "conservative", "log"] = "linear",
        tool: Literal["xgcm"] = "xgcm",
        **kwargs) -> xarray__Dataset:
    """
    Regrid data_var to output_grid.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.regridder.horizontal.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param output_grid: xarray.DataArray or xarray.Dataset
        Grid to transform data_var to
    :param method: {"linear", "conservative", "log"}, optional
        Regridding method to apply; e.g., method = "conservative".
        If tool is "xgcm": "linear", "conservative", "log".
        Default is "linear"
    :param tool: {"xgcm"}, optional
        Name of the tool to use; e.g., tool = "xgcm"
    **kwargs – Additional keyword arguments passed on to the regridder.
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with the data_var transformed to the output_grid.
    """
    return ds.regridder.horizontal(data_var, output_grid, method=method, tool=tool, **kwargs)


def set_auto_bounds(ds: xarray__Dataset, cf_dim: list[Literal["T", "X", "Y", "Z"]] = None) -> xarray__Dataset:
    """
    Adds missing coordinate bounds for supported axes in the Dataset.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.bounds.add_missing_bounds.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param cf_dim: list[{"T", "X", "Y", "Z"}], optional
        List of CF axes that function should operate on. Supported CF axes include “X”, “Y”, “Z”, and “T”.
        Default is None (i.e., ["X", "Y", "Z"])
    
    Output:
    -------
    :return: xarray.Dataset
        Dataset with new bounds where missing.
    """
    if isinstance(cf_dim, list) is False:
        cf_dim = ["X", "Y", "Z"]
    return ds.bounds.add_missing_bounds(axes=cf_dim)


def weights_spatial(ds: xarray__Dataset, data_var: str, cf_dim: list[Literal["X", "Y"]] = None) -> xarray__DataArray:
    """
    Return a DataArray with area weights for specified ‘cf_dim’.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.spatial.SpatialAccessor.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param cf_dim: list[{"X", "Y"}]
        List of axis dimensions to average over, valid axis keys include 'X' and 'Y'; e.g., cf_axis = ["X", "Y"].
        Default is None (i.e., ['X', 'Y'])
    
    Output:
    -------
    :return: xarray.DataArray
        DataArray containing the area weights to use during averaging. Weights are per given dimension in ‘cf_dim’.
    """
    if isinstance(cf_dim, list) is False:
        cf_dim = ["X", "Y"]
    # area weights
    return ds.spatial.get_weights(axis=cf_dim, data_var=data_var)


def weights_temporal(ds: xarray__Dataset, data_var: str) -> xarray__DataArray:
    """
    Return a DataArray with time weights based on a specified frequency.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.temporal.TemporalAccessor.html
    
    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    
    Output:
    -------
    :return: xarray.DataArray
        DataArray containing the time weights to use during averaging.
    """
    # get time axis name
    cf_time = xcdat.axis.get_dim_keys(ds, axis="T")
    # get frequency
    frequency = xcdat.temporal._infer_freq(ds[cf_time])
    # set attributes (necessary to compute weights using temporal._get_weights)
    ds.temporal._set_data_var_attrs(data_var)
    ds.temporal._set_arg_attrs("average", frequency, True)
    # add missing bounds if necessary (necessary to compute weights using temporal._get_weights)
    ds.bounds.add_missing_bounds(axes=["T"])
    # get time bounds
    time_bounds = ds.bounds.get_bounds("T", var_key=data_var)
    # time weights
    return ds.temporal._get_weights(time_bounds)
# ---------------------------------------------------------------------------------------------------------------------#
