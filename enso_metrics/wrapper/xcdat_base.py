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
# numpy
from numpy import ndarray as numpy__ndarray
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
        season_config: Union[dict[list, Union[bool, str]], None] = None,
        skipna: Union[bool, None] = None,
        weighted: bool = True,
        **kwargs) -> xarray__Dataset:
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
    :param reference_period: tuple[str, str], None, optional
        The climatological reference period, which is a subset of the entire time series. This parameter must be
        tuple of strings in the format ‘yyyy-mm-dd’; e.g., reference_period = ('1850-01-01', '1899-12-31').
        If no value is provided, the climatological reference period will be the full period covered by the dataset.
        Default is None
    :param season_config: dict[list, Union[bool, str]], None, optional
        A dictionary for “season” frequency configurations. If configs for predefined seasons are passed, configs for
        custom seasons are ignored and vice versa.
            - “drop_incomplete_seasons” (bool, by default False)
              Seasons are considered incomplete if they do not have all the required months to form the season.
            - “dec_mode” (Literal[“DJF”, “JFD”], by default “DJF”)
              The mode for the season that includes December in the list of list of pre-defined seasons (“DJF”/“JFD”,
              “MAM”, “JJA”, “SON”). This config is ignored if the custom_seasons config is set.
            - “custom_seasons” ([list[list[str]]], by default None)
              List of sublists containing month strings, with each sublist representing a custom season. Month strings
              must be in the three letter format (e.g., ‘Jan’). Order of the months in each custom season does not
              matter. Custom seasons can vary in length.
        Default is None
    :param skipna: bool, None, optional
        If True, skip missing values (as marked by NaN). By default, only skips missing values for float dtypes; other
        dtypes either do not have a sentinel missing value (int) or skipna=True has not been implemented
        (object, datetime64 or timedelta64).
        Default is None
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Default is True
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with the climatology of given data variable.
    """
    tmp_kwargs: dict[str, Union[bool, tuple[str, str]]] = {"keep_weights": keep_weights, "weighted": weighted}
    for k1, k2 in zip(["reference_period", "season_config", "skipna"], [reference_period, season_config, skipna]):
        if k2 is not None:
            tmp_kwargs[k1] = k2
    # operations on ds changes bounds and xcdat doesn't like that: time bounds must be deleted and recreated
    # get time dimension key
    dim_time = xcdat.get_dim_keys(ds, "T")
    # get time bounds key
    dim_time_bnds = ds[dim_time].attrs["bounds"]
    # delete current time bounds
    ds = ds.drop_vars([dim_time_bnds])
    # use xcdat to set bounds
    ds = ds.bounds.add_missing_bounds(axes=("T",))
    return ds.temporal.climatology(data_var, frequency, **tmp_kwargs)


def average_temporal(
        ds: xarray__Dataset,
        data_var: str,
        keep_weights: bool = False,
        skipna: Union[bool, None] = None,
        weighted: bool = True,
        **kwargs) -> xarray__Dataset:
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
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN); e.g., skipna = None.
        Only skips missing values for float dtypes; other dtypes either do not have a sentinel missing value (int) or
        skipna=True has not been implemented (object, datetime64 or timedelta64).
        Default is None
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Weights are calculated by first determining the length of time for each coordinate point using the difference of
        its upper and lower bounds. The time lengths are grouped, then each time length is divided by the total sum of
        the time lengths to get the weight of each coordinate point. The weight of masked (missing) data is excluded
        when averages are taken. This is the same as giving them a weight of 0.
        Default is True
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with the average of a data variable and the time dimension removed.
    """
    # temporal average
    return ds.temporal.average(data_var, keep_weights=keep_weights, skipna=skipna, weighted=weighted)


def create_axis(
        name: str,
        data: Union[list[Union[int, float]], numpy__ndarray],
        bounds: Union[list[Union[int, float]], numpy__ndarray, None] = None,
        generate_bounds: bool = True,
        attrs: Union[dict[str, str], None] = None,
        **kwargs) -> xarray__DataArray:
    """
    Creates an axis and optional bounds.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.create_axis.html

    Inputs:
    -------
    :param name: str
        The CF standard name for the axis (e.g., “longitude”, “latitude”, “height”). xCDAT also accepts additional names
        such as “lon”, “lat”, and “lev”. Refer to xcdat.axis.VAR_NAME_MAP for accepted names.
    :param data:
        1-D axis data consisting of integers or floats.
    :param bounds: list[int | float] or numpy__ndarray
        2-D axis bounds data consisting of integers or floats, defaults to None. Must have a shape of n x 2, where n is the length of data.
    :param generate_bounds: list[int | float] or numpy__ndarray or None, optional
        Generate bounds for the axis if bounds is None, by default True.
    :param attrs: dict[str, str] or None, optional
        Custom attributes to be added to the generated xr.DataArray axis, by default None.
        User provided attrs will be merged with a set of default attributes.
        Default attributes (“axis”, “coordinate”, “bnds”) cannot be overwritten. The default “units” attribute is the
        only default that can be overwritten.
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.Dataset
        New DataArray containing the axis data and optional bounds.
    """
    return xcdat.create_axis(name, data, bounds=bounds, generate_bounds=generate_bounds, attrs=attrs)


def create_gaussian_grid(nlats: int, **kwargs) -> xarray__Dataset:
    """
    Create a grid with Gaussian latitudes and uniform longitudes.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.create_gaussian_grid.html

    Input:
    ------
    :param nlats: int
        Number of latitudes.
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.Dataset
        New Dataset with new grid, containing Gaussian latitudes.
    """
    return xcdat.create_gaussian_grid(nlats)


def create_grid(
        x: Union[xarray__DataArray, tuple[xarray__DataArray, xarray__DataArray, None], None] = None,
        y: Union[xarray__DataArray, tuple[xarray__DataArray, xarray__DataArray, None], None] = None,
        z: Union[xarray__DataArray, tuple[xarray__DataArray, xarray__DataArray, None], None] = None,
        attrs: Union[dict[str, str], None] = None,
        **kwargs) -> xarray__Dataset:
    """
    Creates a grid dataset using the specified axes.
    https://xcdat.readthedocs.io/en/latest/generated/xcdat.create_grid.html

    Inputs:
    -------
    :param x: xarray.DataArray or tuple[xarray.DataArray, xarray.DataArray, None] or None, optional
        An optional dataarray or tuple of a datarray with optional bounds to use for the “X” axis, by default None.
    :param y: xarray.DataArray or tuple[xarray.DataArray, xarray.DataArray, None] or None, optional
        An optional dataarray or tuple of a datarray with optional bounds to use for the “Y” axis, by default None.
    :param z: xarray.DataArray or tuple[xarray.DataArray, xarray.DataArray, None] or None, optional
        An optional dataarray or tuple of a datarray with optional bounds to use for the “Z” axis, by default None.
    :param attrs: dict[str, str] or None, optional
        Custom attributes to be added to the generated xarray.Dataset.
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.Dataset
        New Dataset with grid axes
    """
    return xcdat.create_grid(x=x, y=y, z=z, attrs=attrs)


def create_uniform_grid(
        lat_start: float,
        lat_stop: float,
        lat_delta: float,
        lon_start: float,
        lon_stop: float,
        lon_delta: float,
        **kwargs) -> xarray__Dataset:
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
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.Dataset
        New Dataset with uniform lat/lon grid.
    """
    return xcdat.create_uniform_grid(lat_start, lat_stop, lat_delta, lon_start, lon_stop, lon_delta)


def get_bounds(
        ds: xarray__Dataset,
        cf_dim: Literal["T", "X", "Y", "Z"],
        data_var: str = None,
        **kwargs) -> Union[xarray__DataArray, xarray__Dataset]:
    """
    Gets coordinate bounds.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.bounds.get_bounds.html

    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param cf_dim: {"T", "X", "Y", "Z"}
        CF axis that function should operate on. Supported CF axes include “X”, “Y”, “Z”, and “T”.
    :param data_var: str, optional
        Data variable in ds; e.g., data_var = "ts"
        The data variable to get axis bounds for. This parameter is useful if you only want the single bounds DataArray
        related to the axis on the variable (e.g., “ts” has a “lat” dimension, and you want “lat_bnds”).
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        A Dataset of N bounds variables, or a single bounds variable DataArray.
    """
    return ds.bounds.get_bounds(cf_dim, var_key=data_var)


def interannual_anomalies(
        ds: xarray__Dataset,
        data_var: str,
        frequency: Literal["day", "month", "season"] = "month",
        keep_weights: bool = False,
        reference_period: Union[tuple[str, str], None] = None,
        season_config: Union[dict[list, Union[bool, str]], None] = None,
        skipna: Union[bool, None] = None,
        weighted: bool = True,
        **kwargs) -> xarray__Dataset:
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
    :param season_config: dict[list, Union[bool, str]], None, optional
        A dictionary for “season” frequency configurations. If configs for predefined seasons are passed, configs for
        custom seasons are ignored and vice versa.
            - “drop_incomplete_seasons” (bool, by default False)
              Seasons are considered incomplete if they do not have all the required months to form the season.
            - “dec_mode” (Literal[“DJF”, “JFD”], by default “DJF”)
              The mode for the season that includes December in the list of list of pre-defined seasons (“DJF”/“JFD”,
              “MAM”, “JJA”, “SON”). This config is ignored if the custom_seasons config is set.
            - “custom_seasons” ([list[list[str]]], by default None)
              List of sublists containing month strings, with each sublist representing a custom season. Month strings
              must be in the three letter format (e.g., ‘Jan’). Order of the months in each custom season does not
              matter. Custom seasons can vary in length.
        Default is None
    :param skipna: bool, None, optional
        If True, skip missing values (as marked by NaN). By default, only skips missing values for float dtypes; other
        dtypes either do not have a sentinel missing value (int) or skipna=True has not been implemented
        (object, datetime64 or timedelta64).
        Default is None
    :param weighted: bool, optional
        Calculate averages using weights; e.g., weighted = True.
        Default is True
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with the climatological departures (anomalies) for a data variable.
    """
    tmp_kwargs: dict[str, Union[bool, tuple[str, str]]] = {"keep_weights": keep_weights, "weighted": weighted}
    for k1, k2 in zip(["reference_period", "season_config", "skipna"], [reference_period, season_config, skipna]):
        if k2 is not None:
            tmp_kwargs[k1] = k2
    # operations on ds changes bounds and xcdat doesn't like that: time bounds must be deleted and recreated
    # get time dimension key
    dim_time = xcdat.get_dim_keys(ds, "T")
    # get time bounds key
    dim_time_bnds = ds[dim_time].attrs["bounds"]
    # delete current time bounds
    ds = ds.drop_vars([dim_time_bnds])
    # use xcdat to set bounds
    ds = ds.bounds.add_missing_bounds(axes=("T",))
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
        Default is (0, 360)
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


def regridder_horizontal(
        ds: xarray__Dataset,
        data_var: str,
        output_grid: Union[xarray__DataArray, xarray__Dataset],
        method: Literal[
            "bilinear", "conservative", "conservative_normed", "patch", "nearest_s2d", "nearest_d2s"] = "conservative",
        tool: Literal["regrid2", "xesmf"] = "regrid2",
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
    **kwargs – Additional keyword arguments passed on to the regridder.
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with the data_var transformed to the output_grid.
    """
    tmp_kwargs = {"method": method, "tool": tool, **kwargs}
    return ds.regridder.horizontal(data_var, output_grid, **tmp_kwargs)


def regridder_vertical(
        ds: xarray__Dataset,
        data_var: str,
        output_grid: Union[xarray__DataArray, xarray__Dataset],
        tool: Literal["xgcm"] = "xgcm",
        **kwargs) -> xarray__Dataset:
    """
    Regrid data_var to output_grid.
    https://xcdat.readthedocs.io/en/latest/generated/xarray.Dataset.regridder.vertical.html

    Inputs:
    -------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param data_var: str
        Data variable in ds; e.g., data_var = "ts"
    :param output_grid: xarray.DataArray or xarray.Dataset
        Grid to transform data_var to
    :param tool: {"xgcm"}, optional
        Name of the tool to use; e.g., tool = "xgcm"
    **kwargs – Additional keyword arguments passed on to the regridder.
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with the data_var transformed to the output_grid.
    """
    return ds.regridder.vertical(data_var, output_grid, tool=tool, **kwargs)


def set_auto_bounds(
        ds: xarray__Dataset,
        cf_dim: list[Literal["T", "X", "Y", "Z"]] = None,
        **kwargs) -> xarray__Dataset:
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
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        Input object with new bounds where missing.
    """
    if not isinstance(cf_dim, list):
        cf_dim = ["X", "Y", "Z"]
    return ds.bounds.add_missing_bounds(axes=cf_dim)


def weights_spatial(
        ds: xarray__Dataset,
        data_var: str,
        cf_dim: list[Literal["X", "Y"]] = None,
        **kwargs) -> xarray__DataArray:
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
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray containing the area weights to use during averaging. Weights are per given dimension in ‘cf_dim’.
    """
    if not isinstance(cf_dim, list):
        cf_dim = ["X", "Y"]
    # area weights
    return ds.spatial.get_weights(axis=cf_dim, data_var=data_var)


def weights_temporal(
        ds: xarray__Dataset,
        data_var: str,
        **kwargs) -> xarray__DataArray:
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
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        New DataArray containing the time weights to use during averaging.
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
