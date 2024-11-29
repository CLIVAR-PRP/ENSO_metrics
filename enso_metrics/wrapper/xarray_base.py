# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# N-dimensional array with labeled coordinates and dimensions are needed for the package.
# In xarray:
#     - A dataset resembles an in-memory representation of a NetCDF file, and consists of variables, coordinates and
#       attributes which together form a self describing dataset.
#     - A DataArray provides a wrapper around numpy ndarrays that uses labeled dimensions and coordinates to support
#       metadata aware operations.
# This file regroups xarray functions.
# https://docs.xarray.dev/en/latest/
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from inspect import stack as inspect__stack
from typing import Any, Callable, Hashable, Literal, Union

# numpy
from numpy import dtype as numpy__dtype
from numpy import ndarray as numpy__ndarray
# xarray
import xarray
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
array_wrapper = xarray.DataArray
dataset_wrapper = xarray.Dataset


def change_type(
        ds: Union[array_wrapper, dataset_wrapper],
        dtype: Union[str, numpy__dtype],
        casting: Union[Literal["no", "equiv", "safe", "same_kind", "unsafe"], None] = None,
        copy: Union[bool, None] = None,
        data_var: str = None,
        keep_attrs: bool = True,
        order: Union[Literal["C", "F", "A", "K"], None] = None,
        subok: bool = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Copy of the xarray object, with data cast to a specified type. Leaves coordinate dtype unchanged.
    Return a copy of this xarray.DataArray's or xarray.Dataset[data_var]'s data.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.astype.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.astype.html
    https://numpy.org/doc/stable/reference/generated/numpy.ndarray.astype.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dtype: str or dtype
        Typecode or data-type to which the array is cast
    :param casting: {"no", "equiv", "safe", "same_kind", "unsafe"}, optional
        Controls what kind of data casting may occur.
            - ‘no’ means the data types should not be cast at all.
            - ‘equiv’ means only byte-order changes are allowed.
            - ‘safe’ means only casts which can preserve values are allowed.
            - ‘same_kind’ means only safe casts or casts within a kind, like float64 to float32, are allowed.
            - ‘unsafe’ means any data conversions may be done.
        Default is None ("unsafe" according to ndarray.astype)
    :param copy: bool, optional
        Copy or modify input array; e.g., copy = True.
        If True, returns a newly allocated array.
        If False and the dtype requirement is satisfied, the input array is returned.
        Default is None (True according to ndarray.astype)
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param keep_attrs: bool or None, optional
        Keep or not input attributes; e.g., keep_attrs = True.
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is True
    :param order: {"C", "F", "A", "K"}, optional
        Controls the memory layout order of the result.
            - ‘C’ means C order.
            - ‘F’ means Fortran order.
            - ‘A’ means ‘F’ order if all the arrays are Fortran contiguous, ‘C’ order otherwise.
            - ‘K’ means as close to the order the array elements appear in memory as possible.
        Default is None ("K" according to ndarray.astype)
    :param subok: bool, optional
        Use sub-classes or base-class array; e.g., subok = True.
        If True, then sub-classes will be passed-through, otherwise the returned array will be forced to be a
        base-class array.
        Default is None (True according to ndarray.astype)
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray__array or xarray.Dataset
        New object with data cast to the specified type.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # cast object to new type
    return ds.astype(dtype, casting=casting, copy=copy, keep_attrs=keep_attrs, order=order, subok=subok)


def convert_cf_dim_key(ds: Union[array_wrapper, dataset_wrapper], cf_dim: Literal["T", "X", "Y", "Z"], **kwargs) -> str:
    """
    Return dimension name corresponding to CF dimension name (T is time, X is longitude, Y is latitude, Z is depth or
    level).
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param cf_dim: {"X", "Y", "T", "Z"}
        Name of a CF dimension
    **kwargs - Discarded
    
    Output:
    -------
    :return: str
        Name of given CF dimension in input xarray.DataArray or xarray.Dataset.
    """
    # list dimension
    list_dim = list(ds.dims)
    # dimension to find
    dim_to_find = {
        "T": ["time"],
        "X": ["lon"],
        "Y": ["lat"],
        "Z": ["depth", "height", "level", "pressure", "vertical", "lev"]}
    # find the name in the dataset
    dim_o = ""
    for k1 in dim_to_find[cf_dim]:
        for k2 in list_dim:
            if k1 in k2:
                dim_o = copy__deepcopy(k2)
                break
        if dim_o != "":
            break
    if dim_o == "":
        stack = inspect__stack()
        error = "ERROR: file " + str(stack[0][1]) + " ; fct " + str(stack[0][3]) + " ; line " + str(stack[0][2])
        error += "\n" + str().ljust(5) + "cannot find " + str(cf_dim)
        error += "\n" + str().ljust(5) + "dimension(s): " + ", ".join([repr(k) for k in list_dim])
        raise ValueError(error)
    return dim_o


def convert_dim_keys(
        ds: Union[array_wrapper, dataset_wrapper],
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None],
        **kwargs) -> Union[str, list[str]]:
    """
    Return dimension name(s) from input xarray.DataArray or xarray.Dataset.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None
        Name(s) of dimension or CF dimension
    **kwargs - Discarded
    
    Output:
    -------
    :return: str
        Name(s) of dimension in input xarray.DataArray or xarray.Dataset.
    """
    dim_o = None
    if dim is not None:
        # input dimension to list
        dimensions_asked = copy__deepcopy(dim)
        if isinstance(dim, (Hashable, str)) is True:
            dimensions_asked = [dim]
        # dimensions in input xarray.DataArray or xarray.Dataset
        dimensions_available = list(ds.dims)
        # match asked and available dimensions
        dim_o = []
        for k in dimensions_asked:
            if k in dimensions_available:
                dim_o.append(k)
            elif isinstance(k, (Hashable, str)) is True and k in ["T", "X", "Y", "Z"]:
                dim_o.append(convert_cf_dim_key(ds, k))
            else:
                stack = inspect__stack()
                error = "ERROR: file " + str(stack[0][1]) + " ; fct " + str(stack[0][3]) + " ; line " + str(stack[0][2])
                error += "\n" + str().ljust(5) + "cannot find " + str(k)
                error += "\n" + str().ljust(5) + "dimension(s): " + ", ".join([repr(k) for k in dimensions_available])
                raise ValueError(error)
        # list to str if needed
        if isinstance(dim, (Hashable, str)) is True and len(dim_o) == 1:
            dim_o = dim_o[0]
    return dim_o


def copy(
        ds: Union[array_wrapper, dataset_wrapper],
        data: Union[numpy__ndarray, array_wrapper, None] = None,
        data_var: str = None,
        deep: bool = True,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Return a copy of this xarray.DataArray's or xarray.Dataset[data_var]'s data.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.copy.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.copy.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data: numpy.ndarray or xarray.DataArray or None, optional
        Data to use in the new object. Must have the same shape as original. When data is used, deep is ignored for all
        data variables, and only used for coords.
        Default is None
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param deep: bool, optional
        Whether the data array and its coordinates are loaded into memory and copied onto the new object.
        Default is True
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, and optionally data copied from original.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # copy object
    return ds.copy(data=data, deep=deep)


def correlation(
        ds_a: Union[array_wrapper, dataset_wrapper],
        ds_b: Union[array_wrapper, dataset_wrapper],
        data_var_a: str = None,
        data_var_b: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> array_wrapper:
    """
    Compute the Pearson correlation coefficient between two DataArray objects along shared dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.corr.html

    Input:
    ------
    :param ds_a: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param ds_b: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var_a: str, optional
        Data variable in ‘ds_a‘ if it is xarray.Dataset; e.g., data_var_a = "ts".
        If ‘ds_a‘ is xarray.Dataset, ‘data_var_a‘ must be provided.
        Default is None
    :param data_var_b: str, optional
        Data variable in ‘ds_b‘ if it is xarray.Dataset; e.g., data_var_b = "ts".
        If ‘ds_b‘ is xarray.Dataset, ‘data_var_b‘ must be provided.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension[s] along which to apply var; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param weights: xarray.DataArray, optional
        Array of weights.
        Default is None
     **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        DataArray with the correlation and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds_a = to_array(ds_a, data_var_a)
    ds_b = to_array(ds_b, data_var_b)
    # get dimension(s) as named in xarray.DataArray
    dim_name = convert_dim_keys(ds_a, dim)
    # correlation value
    return xarray.corr(ds_a, ds_b, dim=dim_name, weights=weights)


def create_array_zero(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> array_wrapper:
    """
    Return a new DataArray of zero with the same shape, axes, coordinates, attributes,... as input DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.zeros_like.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object of zeros with the same shape and type as ‘ds‘.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # create an array of zeros with the same shape, axes, coordinates, attributes,... as input DataArray
    return xarray.zeros_like(ds)


def drop_dataset_keys(
        ds: dataset_wrapper,
        names: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str]],
        errors: Literal["raise", "ignore"] = "raise",
        **kwargs) -> dataset_wrapper:
    """
    Drop variables from this dataset.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.drop_vars.html#xarray.Dataset.drop_vars
    
    Input:
    ------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    :param names: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str]
        Name(s) of variables to drop; e.g., "ts" or ["ts", "zos"]
    :param errors: {"raise", "ignore"}, optional
        Controls how to handle errors; e.g., errors = "raise".
            - raise’, raises a ValueError error if any of the variable passed are not in the dataset.
            - ‘ignore’, any given names that are in the dataset are dropped and no error is raised.
        Default is "raise"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object of zeros with the same shape and type as ‘ds‘.
    """
    return ds.drop_vars(names, errors=errors)


def expand_dim(
        ds: Union[array_wrapper, dataset_wrapper],
        axis: Union[int, list[int], tuple[int], None] = None,
        data_var: str = None,
        dim: Union[dict[str, Union[int, numpy__ndarray, array_wrapper]], str, list[str], tuple[str] or None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Return a new object with an additional axis (or axes) inserted at the corresponding position in the array shape.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.expand_dims.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.expand_dims.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param axis: int or list[int] or tuple[int], optional
        Axis position(s) where new axis is to be inserted (position(s) on the result array); e.g., axis=1 or axis=[1,2].
        If a sequence of integers is passed, multiple axes are inserted. In this case, dim arguments should be same
        length list.
        If axis=None is passed, all the axes will be inserted to the start of the result array.
        Default is None
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: dict[str, int or numpy.ndarray or xarray.DataArray] or str or list[str] or tuple[str] or None, optional
        Dimensions to include on the new variable; e.g., dim = {"latitude": array_latitude}.
        If provided as a dict, then the keys are the new dimensions and the values are either integers (giving the
        length of the new dimensions) or ndarray/DataArray (giving the coordinates of the new dimensions).
        If provided as str or sequence of str, then dimensions are inserted with length 1.
        If None, array filled with NaN.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return da_ones: xarray.DataArray or xarray.Dataset
        Input object, but with additional dimension(s).
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # expand array
    return ds.expand_dims(axis=axis, dim=dim)


def get_array_name(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> str:
    """
    Return name of given xarray.DataArray or xarray.Dataset[data_var].
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.name.html
    
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
    :return: str
        The name of this array.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get array name
    return str(ds.name)


def get_array_shape(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> tuple[int, ...]:
    """
    Return shape of given xarray.DataArray or xarray.Dataset[data_var].
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.shape.html
    
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
    :return: str
        Tuple of array dimensions.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get array shape
    return ds.shape


def get_array_size(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> int:
    """
    Return shape of given xarray.DataArray or xarray.Dataset[data_var].
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.size.html
    
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
    :return: str
        Number of elements in the array. Equal to numpy.prod(a.shape), i.e., the product of the array's dimensions.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # read given attribute
    return ds.size


def get_attribute(
        ds: Union[array_wrapper, dataset_wrapper],
        attribute_name: str,
        data_var: str = None,
        **kwargs) -> str:
    """
    Return the given global attribute of given xarray.Dataset or the given variable attribute of given xarray.DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.attrs.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.attrs.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param attribute_name: str
        Name of desired attribute; e.g., attribute_name = "units"
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, ‘attribute_name‘ must be a variable attribute.
        Else, ‘attribute_name‘ must be a global attribute.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: str
        Global attribute or variable attribute.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get given attribute
    return ds.attrs[attribute_name]


def get_attributes(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> dict[str, str]:
    """
    Return a dictionary of global attributes of given xarray.Dataset or variable attributes of given xarray.DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.attrs.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.attrs.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, a dictionary of variable attributes is returned.
        Else, a dictionary of global attributes is returned.
        Default is None
    **kwargs - Discarded

    Output:
    -------
    :return: dict[str, str]
        Global attributes or variable attributes.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get attributes dictionary
    return ds.attrs


def get_attributes_keys(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> list[str]:
    """
    Return the list of global attribute names of given xarray.Dataset or the list of variable attribute names of given
    xarray.DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.attrs.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.attrs.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, a list of variable attribute names is returned.
        Else, a list of global attribute names is returned.
        Default is None
    **kwargs - Discarded

    Output:
    -------
    :return: list[str]
        List of global attribute names or list of variable attribute names.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # list attribute names
    return sorted(list(ds.attrs.keys()), key=lambda v: v.lower())


def get_dataset_keys(ds: dataset_wrapper, **kwargs) -> list[Hashable]:
    """
    List variables in Dataset.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.keys.html
    
    Input:
    ------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    **kwargs - Discarded
    
    Output:
    -------
    :return: list[Hashable]
        List of variables names in the dataset.
    """
    return sorted(list(ds.keys()), key=lambda v: v.lower())


def get_dim_array(
        ds: Union[array_wrapper, dataset_wrapper],
        dim: Union[Hashable, str], **kwargs) -> array_wrapper:
    """
    Return dimension array of given xarray.DataArray or xarray.Dataset corresponding to given dimension name.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dim: Hashable or str
        Name of dimension; e.g., dim="x"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        Dimension DataArray.
    """
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # get DataArray
    return ds[dim_name]


def get_dim_keys(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> list[Hashable]:
    """
    Return dimension name(s) from input xarray.DataArray xarray.Dataset.
    If ‘ds‘ is xarray.DataArray (or xarray.Dataset and ‘data_var‘ in ‘ds‘), dimension names are ordered as in DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.coords.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.dims.html
    
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
    :return: list[Hashable]
        Dimension name(s) ordered as in input xarray.DataArray or xarray.Dataset.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension keys
    if isinstance(ds, dataset_wrapper) is True:
        dim_o = list(ds.coords.keys())
    else:
        dim_o = list(ds.dims)
    return dim_o


def get_time_bounds(ds: Union[array_wrapper, dataset_wrapper], **kwargs) -> list[str]:
    """
    Return first and last time values of given xarray.DataArray or xarray.Dataset.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    **kwargs - Discarded

    Output:
    -------
    :return: list[str]
        First and last time values.
    """
    dim_name = convert_cf_dim_key(ds, "T")
    time_initial = ds[dim_name][0].values.item().isoformat().split("T")[0]
    time_final = ds[dim_name][-1].values.item().isoformat().split("T")[0]
    return [time_initial, time_final]


def maximum(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying max along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.max.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.max.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply max; e.g., dim = "X".
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating max on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with max applied to its data and the
        indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray
    dim_name = convert_dim_keys(ds, dim)
    # maximum value
    return ds.max(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def mean(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying mean or weighted mean along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.mean.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.mean.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.mean.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply mean; e.g., dim = "X".
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).mean() is computed.
        Else, ds.mean() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating mean on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with mean or weighted mean applied to its
        data and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # mean value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).mean(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.mean(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def median(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying median or weighted quantile along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.median.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.median.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.quantile.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply median; e.g., dim = "X".
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).quantile(50) is computed.
        Else, ds.median() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating median on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with median or weighted quantile applied to
        its data and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # median value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).quantile(50, dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.median(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def minimum(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying min along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.min.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.min.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply min; e.g., dim = "X".
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating min on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with min applied to its data and the
        indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # minimum value
    return ds.min(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def numpy_to_array(
        arr: numpy__ndarray,
        coordinates: dict,
        dimensions: Union[list[Hashable], list[str], tuple[str], tuple[Hashable]],
        attrs: Union[dict[str, Union[float, int, str]], None] = None,
        name: Union[str, None] = None,
        **kwargs) -> array_wrapper:
    """
    Return an N-dimensional array with labeled coordinates and dimensions.
    https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html
    
    Input:
    ------
    :param arr: array_like
        Values for this array. Must be a numpy.ndarray, ndarray like, or castable to a ndarray
    :param coordinates: dict[str, array_like]
        Coordinates to use for indexing along each dimension
    :param dimensions: list[Hashable] or list[str] or tuple[Hashable] or tuple[str]
        Name(s) of the data dimension(s)
    :param attrs: dict[str, float or int or str] or None, optional
        Attributes to assign to the new instance.
        Default is None
    :param name: str or None, optional
        Name of this array.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        N-dimensional array with labeled coordinates and dimensions.
    """
    return array_wrapper(attrs=attrs, coords=coordinates, data=arr, dims=dimensions, name=name)


def open_dataset(
        paths: Union[str, list[str]],
        attrs_file: Union[str, list[str], None] = None,
        chunks: Union[Literal["auto"], int, dict, None] = None,
        combine: Literal["by_coords", "nested"] = "by_coords",
        combine_attrs: Union[
            Literal["drop", "identical", "no_conflicts", "drop_conflicts", "override"], Callable] = "override",
        compat: Literal["identical", "equals", "broadcast_equals", "no_conflicts", "override"] = "no_conflicts",
        concat_dim: Union[str, array_wrapper, None] = None,
        coords: Union[Literal["minimal", "different", "all"], list[str]] = "different",
        data_vars: Union[Literal["minimal", "different", "all"], list[str]] = "all",
        engine: Literal["netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None] = None,
        join: Literal["outer", "inner", "left", "right", "exact", "override"] = "outer",
        parallel: bool = False,
        preprocess: Union[Callable, None] = None,
        **kwargs) -> dataset_wrapper:
    """
    Open multiple files as a single dataset.
    https://docs.xarray.dev/en/latest/generated/xarray.open_mfdataset.html
    
    Input:
    ------
    :param paths: str or list[str]
        Paths to dataset files:
            - Directory path (e.g., "path/to/files"), which is converted to a string glob of *.nc files
            - String glob (e.g., "path/to/files/*.nc"), which is expanded to a 1-dimensional list of file paths
            - File path to dataset (e.g., "path/to/files/file1.nc")
            - List of file paths (e.g., ["path/to/files/file1.nc", ...])
        If multiple files, concatenation along the time dimension recommended
    :param attrs_file: str or list[str], optional
        Path of the file used to read global attributes from.
        If None, global attributes are read from the first file provided, with wildcard matches sorted by filename.
        Default is None
    :param chunks: int or dict or "auto" or None, optional
        Dictionary with keys given by dimension names and values given by chunk sizes. In general, these should divide
        the dimensions of each dataset.
        If int, chunk each dimension by chunks.
        If None, chunks will be chosen to load entire input files into memory at once.
        Default is None
    :param combine: {"by_coords", "nested"}, optional
        Whether xarray.combine_by_coords or xarray.combine_nested is used to combine all the data.
        Default is "by_coords"
    :param combine_attrs: {"drop", "identical", "no_conflicts", "drop_conflicts", "override"} or Callable, optional
        A callable or a string indicating how to combine attrs of the objects being merged:
            - “drop”: empty attrs on returned Dataset.
            - “identical”: all attrs must be the same on every object.
            - “no_conflicts”: attrs from all objects are combined, any that have the same name must also have the same
                              value.
            - “drop_conflicts”: attrs from all objects are combined, any that have the same name but different values
                                are dropped.
            - “override”: skip comparing and copy attrs from the first dataset to the result.
        If a callable, it must expect a sequence of attrs dicts and a context object as its only parameters.
        Default is "override"
    :param compat: {"identical", "equals", "broadcast_equals", "no_conflicts", "override"}, optional
        String indicating how to compare variables of the same name for potential conflicts when merging:
            - “broadcast_equals”: all values must be equal when variables are broadcast against each other to ensure
                                  common dimensions.
            - “equals”: all values and dimensions must be the same.
            - “identical”: all values, dimensions and attributes must be the same.
            - “no_conflicts”: only values which are not null in both datasets must be equal. The returned dataset then
                              contains the combination of all non-null values.
            - “override”: skip comparing and pick variable from first dataset
        Default is "no_conflicts"
    :param concat_dim: str or DataArray or None, optional
        Dimensions to concatenate files along. You only need to provide this argument if combine="nested", and if any of
        the dimensions along which you want to concatenate is not a dimension in the original datasets, e.g., if you
        want to stack a collection of 2D arrays along a third dimension. Set concat_dim=[..., None, ...] explicitly to
        disable concatenation along a particular dimension.
        If None, which for a 1D list of filepaths is equivalent to opening the files separately and then merging them
        with xarray.merge.
        Default is None
    :param coords: {"minimal", "different", "all"} or list of str, optional
        These coordinate variables will be concatenated together:
            - “minimal”: Only coordinates in which the dimension already appears are included.
            - “different”: Coordinates which are not equal (ignoring attributes) across all datasets are also
                           concatenated (as well as all for which dimension already appears). Beware: this option may
                           load the data payload of coordinate variables into memory if they are not already loaded.
            - “all”: All coordinate variables will be concatenated, except those corresponding to other dimensions.
            - list of str: The listed coordinate variables will be concatenated, in addition the “minimal” coordinates.
        Default is "different"
    :param data_vars: {"minimal", "different", "all"} or list of str
        These data variables will be concatenated together:
            - “minimal”: Only data variables in which the dimension already appears are included.
            - “different”: Data variables which are not equal (ignoring attributes) across all datasets are also
                           concatenated (as well as all for which dimension already appears). Beware: this option may
                           load the data payload of data variables into memory if they are not already loaded.
            - “all”: All data variables will be concatenated.
            - list of str: The listed data variables will be concatenated, in addition to the “minimal” data variables.
        Default is "all"
    :param engine: {"netcdf4", "scipy", "pydap", "h5netcdf", "zarr", None}, optional
        Engine to use when reading files. If not provided, the default engine is chosen based on available dependencies,
        with a preference for “netcdf4”.
        Default is None
    :param join: {"outer", "inner", "left", "right", "exact", "override"}, optional
        String indicating how to combine differing indexes (excluding concat_dim) in objects:
            - “outer”: use the union of object indexes.
            - “inner”: use the intersection of object indexes.
            - “left”: use indexes from the first object with each dimension.
            - “right”: use indexes from the last object with each dimension.
            - “exact”: instead of aligning, raise ValueError when indexes to be aligned are not equal.
            - “override”: if indexes are of same size, rewrite indexes to be those of the first object with that
                          dimension. Indexes for the same dimension must have the same size in all objects.
        Default is "outer"
    :param parallel: bool, optional
        If True, the open and preprocess steps of this function will be performed in parallel using dask.delayed.
        Default is False
    :param preprocess: Callable, optional
        If provided, call this function on each dataset prior to concatenation. You can find the file-name from which
        each dataset was loaded in ds.encoding["source"].
        Default is None
    **kwargs – Additional keyword arguments passed on to xarray.open_dataset(). For an overview of some possible
               options, see the documentation of xarray.open_dataset()
               https://docs.xarray.dev/en/latest/generated/xarray.open_dataset.html

    Output:
    -------
    :return: xarray.Dataset
        Newly created dataset.
    """
    tmp_kwargs = {"attrs_file": attrs_file, "chunks": chunks, "combine": combine, "combine_attrs": combine_attrs,
                  "compat": compat, "concat_dim": concat_dim, "coords": coords, "data_vars": data_vars,
                  "engine": engine, "join": join, "parallel": parallel, "preprocess": preprocess, **kwargs}
    return xarray.open_mfdataset(paths, **tmp_kwargs)


def polyfit(
        ds: Union[array_wrapper, dataset_wrapper],
        dim: Union[Hashable, str],
        deg: int,
        cov: Union[bool, Literal["unscaled"]] = False,
        full: bool = False,
        rcond: Union[float, None] = None,
        skipna: Union[bool, None] = None,
        w: Union[Hashable, numpy__ndarray, Any] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Least squares polynomial fit.
    This replicates the behaviour of numpy.polyfit but differs by skipping invalid values when skipna = True.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.polyfit.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.polyfit.html
    
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dim: Hashable or str
        Name of dimension along which to apply polyfit; e.g., dim = "x"
    :param deg: int
        Degree of the fitting polynomial; e.g., deg = 1
    :param cov: bool or {"unscaled"}, optional
        Whether to return to the covariance matrix in addition to the coefficients; e.g., cov = False.
        The matrix is not scaled if cov = "unscaled".
        Default is False
    :param full: bool, optional
        Whether to return the residuals, matrix rank and singular values in addition to the coefficients;
        e.g., full = False.
        Default is False
    :param rcond: float or None, optional
        Relative condition number to the fit; e.g., rcond = None.
        Default is None
    :param skipna: bool or None, optional
        If True, removes all invalid values (as marked by NaN) before fitting each 1D slices of the array.
        Default is None (i.e., True if data is stored in a dask.array or if there is any invalid values, False
        otherwise)
    :param w: Hashable or numpy.ndarray or Any, optional
        Weights to apply to the y-coordinate of the sample points. Can be an array-like object or the name of a
        coordinate in the dataset.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        A single dataset which contains:
            - “polyfit_coefficients”: The coefficients of the best fit.
            - “polyfit_residuals”:    The residuals of the least-square computation (only included if full=True). When
                                      the matrix rank is deficient, numpy.nan is returned.
            - “[dim]_matrix_rank”:    The effective rank of the scaled Vandermonde coefficient matrix (only included if
                                      full=True).
            - “[dim]_singular_value”: The singular values of the scaled Vandermonde coefficient matrix (only included if
                                      full=True).
            - “polyfit_covariance”:   The covariance matrix of the polynomial coefficient estimates (only included if
                                      full=False and cov=True).
        If ‘ds‘ is xarray.Dataset, “[var]_” is added at the beginning of “polyfit_coefficients”, “polyfit_residuals” and
        “polyfit_covariance” (for each “var” in the input dataset).
    """
    return ds.polyfit(dim, deg, cov=cov, full=full, rcond=rcond, skipna=skipna, w=w)


def polyval(
        coords: Union[array_wrapper, dataset_wrapper],
        coeffs: Union[array_wrapper, dataset_wrapper],
        dim_degree: str = "degree",
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Evaluate a polynomial at specific values.
    https://docs.xarray.dev/en/stable/generated/xarray.polyval.html
    
    Input:
    ------
    :param coords: xarray.DataArray or xarray.Dataset
        DataArray or Dataset of values at which to evaluate the polynomial
    :param coeffs: xarray.DataArray or xarray.Dataset
        DataArray or Dataset of coefficients of the polynomial
    :param dim_degree: str, optional
        Name of the polynomial degree dimension in ‘coeffs‘; e.g. dim_degree = "degree"
        Default is "degree"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        This object with evaluated polynomial.
    """
    return xarray.polyval(coords, coeffs, degree_dim=dim_degree)


def quantile(
        ds: Union[array_wrapper, dataset_wrapper],
        q: Union[float, int, list[float], list[int], tuple[float], tuple[int]],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        method: Literal["inverted_cdf", "averaged_inverted_cdf", "closest_observation", "interpolated_inverted_cdf",
                        "hazen", "weibull", "linear", "median_unbiased", "normal_unbiased"] = "linear",
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying quantile or weighted quantile along given dimension(s).
    If weighted quantile, the only interpolation method supported corresponds to method = "linear".
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.quantile.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.quantile.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.quantile.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param q: float or int or list[float] or list[int] or tuple[float] or tuple[int]
        Quantile(s) to compute, which must be between 0 and 1 inclusive; e.g., q = 0.5 or q = [0.25, 0.75]
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply quantile; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param method: {"inverted_cdf", "averaged_inverted_cdf", "closest_observation", "interpolated_inverted_cdf",
                    "hazen", "weibull", "linear", "median_unbiased", "normal_unbiased"}, optional
        Interpolation method to use when the desired quantile lies between two data points; e.g., method = "linear":
            - “inverted_cdf”, “averaged_inverted_cdf”, “closest_observation”, “interpolated_inverted_cdf”, “hazen”,
              “weibull”, “linear”, “median_unbiased”, “normal_unbiased”
        Default is "linear"
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).quantile() is computed.
        Else, ds.quantile() is computed.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        If q is a single quantile, then a new object with quantile applied to its data and the indicated dimension(s)
        removed.
        If multiple percentiles are given, then first axis of the new object corresponds to the quantile, a quantile
        dimension is added and the indicated dimension(s).
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # quantile value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).quantile(q, dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.quantile(q, dim=dim_name, keep_attrs=keep_attrs, method=method, skipna=skipna)


def rename(
        ds: Union[array_wrapper, dataset_wrapper],
        name_dict: Union[dict[str, str], str],
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Returns a new object with renamed variables, coordinates and dimensions.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.rename.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.rename.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param name_dict: str or dict[str, str]
        Dictionary whose keys are current variable, coordinate or dimension names and whose values are the desired
        names.
        If ‘ds‘ is xarray.DataArray and ‘name_dict‘ is str, it as the new name for this array.
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        This object with renamed variables, coordinates and dimensions.
    """
    return ds.rename(name_dict)


def select(
        ds: Union[array_wrapper, dataset_wrapper],
        indexers: dict,
        drop: bool = False,
        method: Union[Literal["nearest", "pad", "ffill", "backfill", "bfill"], None] = None,
        tolerance=None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Returns a new object with each array indexed by tick labels along the specified dimension(s).
    In contrast to isel, indexers for this method should use labels instead of integers.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.sel.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.sel.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param indexers: dict
        A dict with keys matching dimensions and values given by scalars, slices or arrays of tick labels.
        For dimensions with multi-index, the indexer may also be a dict-like object with keys matching index level
        names. If DataArrays are passed as indexers, xarray-style indexing will be carried out.
        See Indexing and selecting data for the details
    :param drop: bool, optional
        If drop=True, drop coordinates variables in indexers instead of making them scalar; e.g., drop = False.
        Default is False
    :param method: {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
        Method to use for inexact matches; method = None.
            - None: only exact matches
            - “pad” / “ffill”: propagate last valid index value forward
            - “backfill” / “bfill”: propagate next valid index value backward
            - “nearest”: use nearest valid index value
        Default is None
    :param tolerance:
        Maximum distance between original and new labels for inexact matches. The values of the index at the matching
        locations must satisfy the equation abs(index[indexer] - target) <= tolerance
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with the same contents as ‘ds‘, except each variable and dimension is indexed by the appropriate
        indexers. If indexer DataArrays have coordinates that do not conflict with this object, then these coordinates
        will be attached. In general, each array’s data will be a view of the array’s data in ‘ds‘, unless vectorized
        indexing was triggered by using an array indexer, in which case the data will be a copy.
    """
    return ds.sel(drop=drop, indexers=indexers, method=method, tolerance=tolerance)


def select_index(
        ds: Union[array_wrapper, dataset_wrapper],
        indexers: dict,
        drop: bool = False,
        missing_dims: Literal["raise", "warn", "ignore"] = "raise",
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Returns a new object with each array indexed along the specified dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.isel.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.isel.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param indexers: dict
        A dict with keys matching dimensions and values given by integers, slice objects or arrays. indexer can be an
        integer, slice, array-like or DataArray. If DataArrays are passed as indexers, xarray-style indexing will be
        carried out. See Indexing and selecting data for the details
    :param drop: bool, optional
        If drop=True, drop coordinates variables indexed by integers instead of making them scalar; e.g., drop = False.
        Default is False
    :param missing_dims: {"raise", "warn", "ignore"}, optional
        What to do if dimensions that should be selected from are not present in the object:
            - “raise”: raise an exception
            - “warn”: raise a warning, and ignore the missing dimensions
            - “ignore”: ignore the missing dimensions
        Default is "raise"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with the same contents as ‘ds‘, except each array and dimension is indexed by the appropriate
        indexers. If indexer DataArrays have coordinates that do not conflict with this object, then these coordinates
        will be attached. In general, each array’s data will be a view of the array’s data in ‘ds‘, unless vectorized
        indexing was triggered by using an array indexer, in which case the data will be a copy.
    """
    return ds.isel(drop=drop, indexers=indexers, missing_dims=missing_dims)


def set_attributes(
        ds: Union[array_wrapper, dataset_wrapper],
        *args, **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Assign new attrs to this object.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.assign_attrs.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.assign_attrs.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    *args – positional arguments passed into attrs.update.
    **kwargs – keyword arguments passed into attrs.update.

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object equivalent to self.attrs.update(*args, **kwargs).
    """
    return ds.assign_attrs(*args, **kwargs)


def set_attributes_global(ds: dataset_wrapper, *args, **kwargs):
    """
    Update global attributes.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.attrs.html

    Input:
    ------
    :param ds: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    *args – positional arguments passed into attrs.update.
    **kwargs – keyword arguments passed into attrs.update.
    """
    ds.attrs.update(*args, **kwargs)

    
def set_attributes_variable(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        *args, **kwargs):
    """
    Update variable attributes.
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.attrs.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    *args – positional arguments passed into attrs.update.
    **kwargs – keyword arguments passed into attrs.update.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # update attributes
    ds.attrs.update(*args, **kwargs)


def squeeze(
        ds: Union[array_wrapper, dataset_wrapper],
        axis: Union[int, list[int], tuple[int], None] = None,
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        drop: bool = False,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Return a new object with squeezed data.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.squeeze.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.squeeze.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param axis: int or list[int] or tuple[int] or None, optional
        Like dim, but positional; e.g., axis = 0 or axis = [0, 1]
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Selects a subset of the length one dimensions; e.g., dim="x" or dim=["x", "y"].
        If a dimension is selected with length greater than one, an error is raised.
        If None, all length one dimensions are squeezed.
        Default is None
    :param drop: bool, optional
        If drop=True, drop squeezed coordinates instead of making them scalar; e.g., drop = False.
        Default is False
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        This object, but with all or a subset of the dimensions of length 1 removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # squeeze
    return ds.squeeze(axis=axis, dim=dim_name, drop=drop)


def sum_along_axis(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        min_count: Union[int, None] = None,
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying sum or weighted sum along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.sum.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.sum.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.sum.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply sum; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param min_count: int or None, optional
        The required number of valid values to perform the operation. If fewer than min_count non-NA values are present
        the result will be NA. Only used if skipna is set to True or defaults to True for the array's dtype.
        Default is None
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).sum() is computed.
        Else, ds.sum() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating sum on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with sum or weighted sum applied to its
        data and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # sum value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).sum(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.sum(dim=dim_name, keep_attrs=keep_attrs, min_count=min_count, skipna=skipna, **kwargs)


def standard_deviation(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        ddof: int = 0,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying std or weighted std along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.std.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.std.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.std.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param ddof: int
        Delta Degrees of Freedom: the divisor used in the calculation is N - ddof, where N represents the number of
        elements.
        Default is 0
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply std; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).std() is computed.
        Else, ds.std() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating std on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with std or weighted std applied to its
        data and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # standard deviation value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).std(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.std(dim=dim_name, ddof=ddof, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def to_array(ds: Union[array_wrapper, dataset_wrapper], data_var: str, **kwargs) -> array_wrapper:
    """
    Return xarray.DataArray from input ‘ds’, i.e., ‘ds’ if it is xarray.DataArray or ds[data_var] if it is
    xarray.Dataset.
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str
        Data variable in ds if it is xarray.Dataset; e.g., data_var = "ts"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray
        N-dimensional array with labeled coordinates and dimensions.
    """
    if isinstance(ds, dataset_wrapper) is True and data_var in ds.keys():
        return ds[data_var]
    else:
        return ds


def to_dataset(ds: Union[array_wrapper, dataset_wrapper], data_var: str, **kwargs) -> dataset_wrapper:
    """
    Return xarray.Dataset from input ‘ds’, i.e., ‘ds’ if it is xarray.Dataset or ds.to_dataset() if it is
    xarray.DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.to_dataset.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str
        Data variable to name ‘ds’ in a new xarray.Dataset; e.g., data_var = "ts"
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.Dataset
        An in-memory representation of a NetCDF file, and consists of variables, coordinates and attributes which
        together form a self describing dataset
    """
    if isinstance(ds, dataset_wrapper) is True:
        return ds
    else:
        return ds.to_dataset(name=data_var)


def to_netcdf(
        ds: Union[array_wrapper, dataset_wrapper],
        filename: str,
        file_format: Literal["NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", "NETCDF3_CLASSIC"] = "NETCDF4",
        list_of_variables: list[str] = None,
        mode: Literal["a", "w"] = "a",
        **kwargs):
    """
    Write dataset contents to a netCDF file.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.to_netcdf.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param filename: str
        Path and filename to save this dataset
    :param file_format: {"NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", "NETCDF3_CLASSIC"}, optional
        File format for the resulting netCDF file:
            - NETCDF4: Data is stored in an HDF5 file, using netCDF4 API features.
            - NETCDF4_CLASSIC: Data is stored in an HDF5 file, using only netCDF 3 compatible API features.
            - NETCDF3_64BIT: 64-bit offset version of the netCDF 3 file format, which fully supports 2+ GB files, but is
                             only compatible with clients linked against netCDF version 3.6.0 or later.
            - NETCDF3_CLASSIC: The classic netCDF 3 file format. It does not handle 2+ GB files very well.
        All formats are supported by the netCDF4-python library. scipy.io.netcdf only supports the last two formats.
        Default is "NETCDF4"
    :param list_of_variables: list[str]
        List of data variable(s) in ds to save
    :param mode: {"w", "a"}, optional
        Write (‘w’) or append (‘a’) mode.
        If mode=’w’, any existing file at this location will be overwritten.
        If mode=’a’, existing variables will be overwritten.
        Default is "a"
    **kwargs - Discarded
    """
    # Write dataset contents to a netCDF file
    if isinstance(ds, array_wrapper) is True and isinstance(list_of_variables, list) is True and \
            len(list_of_variables) > 0:
        ds[list_of_variables].to_netcdf(filename, format=file_format, mode=mode)
    else:
        ds.to_netcdf(filename, format=file_format, mode=mode)


def to_numpy(ds: Union[array_wrapper, dataset_wrapper], data_var: str = None, **kwargs) -> numpy__ndarray:
    """
    Coerces wrapped data to numpy and returns a numpy.ndarray.
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.to_numpy.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        xarray.DataArray or xarray.Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset, ‘data_var‘ must be provided.
        Default is None
    **kwargs - Discarded
    
    Output:
    -------
    :return: numpy.ndarray
        Input DataArray as ndarray.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # DataArray to numpy
    return ds.to_numpy()


def transpose(
        ds: Union[array_wrapper, dataset_wrapper],
        dimensions: Union[list[Hashable], list[str], tuple[Hashable], tuple[str]],
        data_var: str = None,
        missing_dims: Literal["raise", "warn", "ignore"] = "raise",
        transpose_coords: bool = True,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Transpose dimensions of this object's data.
    https://docs.xarray.dev/en/stable/generated/xarray.Dataset.transpose.html
    https://docs.xarray.dev/en/stable/generated/xarray.DataArray.transpose.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dimensions: list[Hashable] or list[str] or tuple[Hashable] or tuple[str]
        List of dimension names with the desired order
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param missing_dims: {"raise", "warn", "ignore"}, optional
        What to do if dimensions that should be selected from are not present in the DataArray;
        e.g., missing_dims = "raise".
            - “raise”: raise an exception.
            - “warn”: raise a warning, and ignore the missing dimensions.
            - “ignore”: ignore the missing dimensions.
        Default is "raise"
    :param transpose_coords: bool, optional
        If True, also transpose the coordinates of this DataArray; transpose_coords = True.
        Used only if ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘).
        Default is True
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        Input dataset with each array (including) coordinates transposed to the given order.
        Or Input DataArray with array transposed to the given order.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dimensions_name = convert_dim_keys(ds, dimensions)
    # transpose
    tmp_kwargs = {"missing_dims": missing_dims}
    if isinstance(ds, array_wrapper) is True:
        tmp_kwargs["transpose_coords"] = transpose_coords
    return ds.transpose(*dimensions_name, missing_dims=missing_dims, transpose_coords=transpose_coords)


def variance(
        ds: Union[array_wrapper, dataset_wrapper],
        data_var: str = None,
        ddof: int = 0,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        keep_attrs: Union[bool, None] = False,
        skipna: Union[bool, None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Reduce this object's data by applying var or weighted var along given dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.var.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.var.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/latest/generated/xarray.core.weighted.DataArrayWeighted.var.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param ddof: int
        Delta Degrees of Freedom: the divisor used in the calculation is N - ddof, where N represents the number of
        elements.
        Default is 0
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) along which to apply var; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is False
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If given and ‘ds‘ is xarray.DataArray (or ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘),
        ds.weighted(weights).var() is computed.
        Else, ds.var() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating var on this
               object's data. These could include dask-specific kwargs like split_every.
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with var or weighted var applied to its
        data and the indicated dimension(s) removed.
    """
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # get dimension(s) as named in xarray.DataArray or xarray.Dataset
    dim_name = convert_dim_keys(ds, dim)
    # variance value
    if isinstance(ds, array_wrapper) is True and isinstance(weights, array_wrapper) is True:
        return ds.weighted(weights).var(dim=dim_name, keep_attrs=keep_attrs, skipna=skipna)
    else:
        return ds.var(dim=dim_name, ddof=ddof, keep_attrs=keep_attrs, skipna=skipna, **kwargs)


def where(
        ds: Union[array_wrapper, dataset_wrapper],
        cond: Union[numpy__ndarray, array_wrapper, dataset_wrapper],
        data_var: str = None,
        drop: bool = False,
        other: Union[bool, float, int, numpy__ndarray, array_wrapper, dataset_wrapper, None] = None,
        **kwargs) -> Union[array_wrapper, dataset_wrapper]:
    """
    Filter elements from this object according to a condition.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.where.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.where.html
    
    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param cond: numpy.ndarray or xarray.DataArray or xarray.Dataset
        Locations at which to preserve this object's values. dtype must be bool
    :param data_var: str, optional
        Data variable in ‘ds‘ if it is xarray.Dataset; e.g., data_var = "ts".
        If ‘ds‘ is xarray.Dataset and ‘data_var‘ in ‘ds‘, output will be xarray.DataArray.
        Default is None
    :param drop: bool, optional
        If True, coordinate labels that only correspond to False values of the condition are dropped from the result.
        Default is False
    :param other: bool or float or int or numpy.ndarray or xarray.DataArray or xarray.Dataset or None, optional
        Value to use for locations in this object where ‘cond’ is False.
        Default is None (NA where ‘cond’ is False)
    **kwargs - Discarded
    
    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset
        New object with dimensions, attributes, coordinates, name, encoding, with elements from ‘ds’ where ‘cond’ is
        True, otherwise fill in ‘other’.
    """
    tmp_kwargs = {"drop": drop}
    if other is not None:
        tmp_kwargs["other"] = other
    # read variable from xarray.Dataset if needed
    ds = to_array(ds, data_var)
    # where
    return ds.where(cond, **tmp_kwargs)
# ---------------------------------------------------------------------------------------------------------------------#
