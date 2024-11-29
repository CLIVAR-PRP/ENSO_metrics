# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Processors built over xarray and xcdat
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# local functions
from enso_metrics.tools.default import set_instance
from . tool_base import *
from enso_metrics.wrapper.xarray_base import array_wrapper, dataset_wrapper
from enso_metrics.wrapper import xarray_base
from enso_metrics.wrapper import xcdat_base
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def averager(
        ds: dataset_wrapper,
        data_var_ds: str,
        ds_area: Union[array_wrapper, dataset_wrapper, None] = None,
        data_var_area: str = None,
        method_average: Literal["horizontal", "meridional", "temporal", "zonal"] = True,
        **kwargs) -> dataset_wrapper:
    # define CF dimension
    cf_dim: list[Literal["T", "X", "Y"]] = []
    if method_average in ["horizontal"]:
        cf_dim += ["X", "Y"]
    elif method_average in ["meridional", "zonal"]:
        cf_dim += ["Y"] if method_average == "meridional" else ["X"]
    else:
        cf_dim += ["T"]
    # perform average
    ds_o = average_spatiotemporal(ds, data_var_ds, cf_dim, ds_area=ds_area, data_var_area=data_var_area)
    # metadata
    description = str(method_average) + " average computed"
    ds_o = processing_description(ds_o, data_var_ds, "method_processing", description)
    return ds_o


def detrender(
        ds: dataset_wrapper,
        degree: int,
        **kwargs) -> dataset_wrapper:
    return ds


def masker(
        input_dataset: dict,
        input_param: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, None]]]]] = None,
        kwargs_masker: dict = None,
        variable: str = None,
        variable_new: str = None,
        **kwargs) -> Union[dict, None]:
    # Several arguments can be specified to where (used to mask data): drop, other.
    # If desired they must be defined in a dictionary under the keyword 'kwargs_masker'.
    # Other arguments are used in this function to mask data:
    #     - maskland: bool, True to mask over land
    #     - maskocean: bool, True to mask over ocean
    #     - tolerance: float, > 0 to give a tolerance and land / ocean values (e.g., land is where mask = 1, but users
    #                  may want to mask as land where mask >= 1 - tolerance)
    kwargs_masker = set_instance(kwargs_masker, dict, False, {})
    mask_land, mask_ocean, tolerance = False, False, 0
    if "maskland" in list(kwargs_masker.keys()) and isinstance(kwargs_masker["maskland"], bool) is True:
        mask_land = deepcopy(kwargs_masker["maskland"])
    if "maskocean" in list(kwargs_masker.keys()) and isinstance(kwargs_masker["maskocean"], bool) is True:
        mask_ocean = deepcopy(kwargs_masker["maskocean"])
    if "mask_tolerance" in list(kwargs_masker.keys()) and isinstance(kwargs_masker["mask_tolerance"], float) is True:
        tolerance = deepcopy(kwargs_masker["mask_tolerance"])
    # fake loop to be able to break out if an error occurs
    dict_o = None
    for _ in range(1):
        # check if given variables is available
        if variable is None or variable not in list(input_dataset.keys()) or isinstance(input_param, dict) is False or \
                variable not in list(input_param.keys()):
            # error: variable must be defined
            break
        # get array and related input param
        array = input_dataset[variable]["array"]
        metadata = deepcopy(input_dataset[variable]["metadata"])
        # get mask array
        mask = input_param[variable]["mask"]
        if mask in list(input_dataset.keys()):
            mask = input_dataset[mask]["array"]
        elif mask == "estimate":
            mask = create_land_sea_mask(array)
        # mask land or ocean
        if mask_land is True and mask_ocean is True:
            # error: user asked to mask everything
            break
        elif (mask_land is True or mask_ocean is True) and isinstance(array, (array_wrapper, dataset_wrapper)) is True \
                and isinstance(mask, (array_wrapper, dataset_wrapper)) is True:
            if mask_land is True:
                array = xarray_base.where(array, mask < 1 - tolerance, **kwargs_masker)
            else:
                array = xarray_base.where(array, mask > 0 + tolerance, **kwargs_masker)
            # adapt metadata
            # !!!! TO DO !!!!
        # prepare output
        dict_o = {"array": array, "metadata": metadata}
    return dict_o


def reader(
        input_param: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, None]]]]],
        variables: Union[str, list[str]],
        kwargs_reader: dict = None,
        variables_param: dict[str, dict[str, str]] = None,
        **kwargs) -> dict:
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
        add_bounds = deepcopy(kwargs_reader["add_bounds"])
        del kwargs_reader["add_bounds"]
    decode_times: bool = True
    if "decode_times" in list(kwargs_reader.keys()):
        add_bounds = deepcopy(kwargs_reader["decode_times"])
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
            ab, dt = deepcopy(add_bounds), deepcopy(decode_times)
            if kk in ["areacella", "areacello", "landmask"]:
                if "T" in ab:
                    while "T" in ab:
                        ab.remove("T")
                if dt is True:
                    dt = False
            # try to open_dataset and save Dataset in a dictionary using netCDF variables (names) as keys
            try:
                dict_t[nn] = xcdat_base.open_dataset(
                    ff, add_bounds=ab, data_var=nn, decode_times=dt, **kwargs_reader)
            except "can't read (" + str(nn) + ") " + str(ff):
                pass
    # compute variables
    dict_output = {}
    for kk in variables + ["areacella", "areacello", "landmask"]:
        if kk not in list(input_param.keys()) or kk not in list(variables_param.keys()):
            continue
        arr, metadata = compute_variable(
            dict_t, input_param[kk], input_param[kk]["variable"], kk, variables_param[kk], **kwargs)
        import json
        print(json.dumps(metadata, indent=4))
        stop
        dict_output[kk] = {"array": arr, "metadata": metadata}
    if len([k for k in list(dict_output.keys()) if k in variables]) != len(variables):
        dict_output = None
    return dict_output
# ---------------------------------------------------------------------------------------------------------------------#
