# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack

# ENSO_metrics package functions:
from .EnsoCollectionsLib import ReferenceRegions, reference_variables
from . import EnsoErrorsWarnings
from .EnsoPlotLib import metric_variable_names
from .EnsoToolsLib import add_up_errors, simple_stats
from .EnsoUvcdatToolsLib import AverageMeridional, fill_dict_axis, preprocess_ts_polygon, PreProcessTS, \
    Read_data_mask_area, RmsZonal, SaveNetcdf, TimeBounds, TwoVarRegrid
from .KeyArgLib import default_arg_values


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO recipes
# These functions have file names and variable names as inputs and metric as output
#
def stat_box(var_file, var_name, var_areafile, var_areaname, var_landmaskfile, var_landmaskname, var_box, dataset="",
             debug=False, internal_variable_name="ts", metname="", netcdf=False, netcdf_name="", statistic="average",
             **kwargs):
    """
    The stat_box() function computes the 'statistic' of 'var_box' averaged 'internal_variable_name' (lhf, lwr, pr, shf,
    slp, ssh, sst, swr, taux, tauy, thf, ts)

    Author: Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Sep 08 2022

    Inputs:
    ------
    :param var_file: string
        path_to/filename of the file (NetCDF)
    :param var_name: string
        name of a variable in 'var_file'
            lhf:  evap_heat, hfls, latent_heatflux, lhf, lhtfl, slhf, solatent
            lwr:  lw_heat, lwr (dlwrf - ulwrf or rlds - rlus), net_longwave_heatflux_downwards, rls, solongwa, str
            pr:   mtpr, pr, prate, precip
            shf:  hfss, shf, shtfl, sens_heat, sensible_heatflux, sosensib, sshf
            slp:  msl, prmsl, psl, slp
            ssh:  adt, eta_t, height, sea_surface_height, SLA, sossheig, ssh, sshg, zos
            sst:  sea_surface_temperature, skt, sosstsst, sst, tmpsf, tos, ts
            swr:  net_shortwave_heatflux_downwards, soshfldo, rss, ssr, swflx, swr (dswrf - uswrf or rsds - rsus)
            taux: ewss, sozotaux, tau_x, tauu, tauuo, taux, uflx, zonal_wind_stress
            tauy: nsss, meridional_wind_stress, sometauy, tau_y, tauv, tauvo, tauy, vflx
            thf:  net_heating, net_surface_heatflux_downwards, netflux, sohefldo, thf (lhf + lwr + shf + swr), thflx
            ts:   skt, ts
    :param var_areafile: string
        path_to/filename of the file (NetCDF) of the areacell for the given variable
    :param var_areaname: string
        name of areacell variable (areacella, areacello) in 'var_areafile'
    :param var_landmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for the given variable
    :param var_landmaskname: string
        name of landmask variable (land, lsm, lsmask, landmask, mask, sftlf) in 'var_landmaskfile'
    :param var_box: string
        name of box (e.g. 'global') for the given variable
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param internal_variable_name: string, optional
        name of an internal variable name (lhf, lwr, pr, shf, slp, ssh, sst, swr, taux, tauy, thf, ts)
        default value = 'ts' (surface temperature)
    :param metname: string, optional
        default value = '' metric name is not changed
        e.g., metname = 'stat_box_2'
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
    :param statistic: string, optional
        name of statistic to compute (average, skewness, standard deviation, variance)
        default value = 'average' (time averaged value)
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param regridding: dict, optional
        see EnsoUvcdatToolsLib.TwoVarRegrid and EnsoUvcdatToolsLib.Regrid for options
        the aim if to specify if the model is regridded toward the observations or vice versa, of if both model and
        observations are regridded toward another grid
        interpolation tool and method can be specified
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return metric_output: dict
        dive_down_diag, keyerror, method, method_detail, name, nyears, ref, time_frequency, time_period, units, value,
        value_error

    Method:
    -------
        uses tools from cdat library

    """
    # Test given kwargs
    needed_kwarg = ["detrending", "frequency", "min_time_steps", "normalization", "smoothing", "time_bounds"]
    for arg in needed_kwarg:
        if arg not in list(kwargs.keys()):
            kwargs[arg] = default_arg_values(arg)
    # reference variable
    ref_var = deepcopy(reference_variables(internal_variable_name))
    var_ln = ref_var["long_name"]
    var_sn = ref_var["short_name"]
    # reference region
    reg_di = {var_box: deepcopy(ReferenceRegions(var_box))}
    reg_ln = reg_di[var_box]["long_name"]
    reg_sn = reg_di[var_box]["short_name"]
    reg_s2 = ""
    if "regions" in list(reg_di[var_box].keys()):
        reg_list = reg_di[var_box]["regions"]
        for reg in reg_list:
            reg_di[reg] = deepcopy(ReferenceRegions(reg))
            reg_s2 += reg_di[reg]["short_name"]
            if reg != reg_list[-1]:
                reg_s2 += " minus "
    else:
        reg_list = [var_box]
        reg_s2 = deepcopy(reg_sn)
    if internal_variable_name == "ssh" and var_box != "global":
        reg_di["global"] = deepcopy(ReferenceRegions("global"))
        reg_list += ["global"]
        reg_s2 += " GMSL removed"
    # reference statistic
    if statistic in ["average", "skewness", "standard deviation", "variance"]:
        sta_ano = "" if statistic == "average" else "A"
        sta_sh2 = "mean" if statistic == "average" else ("SKE" if statistic == "skewness" else (
                  "STD" if statistic == "standard deviation" else "VAR"))
        sta_sh1 = "ave" if statistic == "average" else str(sta_sh2).lower()
    else:
        sta_ano, sta_sh1, sta_sh2 = None, None, None

    # Define metric attributes
    name = str(sta_sh2) + " of " + str(reg_sn) + " " + str(var_sn) + str(sta_ano)
    units1 = deepcopy(ref_var["units"])
    units2 = deepcopy(units1)
    if statistic == "skewness":
        units2 = ""
    elif statistic == "variance":
        units2 = "degC2" if units2 == "degC" else ("cm2" if units2 == "cm" else ("m4" if units2 == "m2" else (
                 "mm2/day2" if units2 == "mm/day" else ("1e-6 N2/m4" if units2 == "1e-3 N/m2" else (
                 "Pa2" if units2 == "Pa" else ("W2/m4" if units2 == "W/m2" else deepcopy(units2)))))))
    method = str(statistic) + " of " + str(reg_ln) + " (" + str(reg_s2) + ") averaged " + str(var_ln) + " (" + \
             str(var_sn) + ")"
    if sta_ano == "A":
        method += " anomalies"
    method_var = str(var_sn) + ": "
    ref = "Using CDAT averaging"
    if sta_ano == "A":
        ref += " and computing interannual anomalies"
    metric = "stat_box_" + str(sta_sh1) + "_" + str(var_sn).lower() + "_" + str(var_box)
    if metname == "":
        metname = deepcopy(metric)

    # Values in case of error (failure)
    mv, mv_error, dive_down_diag = None, None, {"value": None, "axis": None}
    nbr_year, timebounds, keyerror = None, None, None

    # Create a fake loop to be able to break out if an error occur
    for break_loop in range(1):
        counter = 1
        # ------------------------------------------------
        # 1. Read files
        # ------------------------------------------------
        r1 = reg_list[0]
        if debug is True:
            EnsoErrorsWarnings.debug_mode("\033[92m", metric, 10)
        # 1.1 read file and select the 'r1' region
        mask_l = reg_di[r1]["maskland"] if "maskland" in list(reg_di[r1].keys()) else False
        mask_o = reg_di[r1]["maskocean"] if "maskocean" in list(reg_di[r1].keys()) else False
        if mask_l is True and mask_o is True:
            keyerror = "both land and ocean are masked (" + str(r1) + ": land = " + str(mask_l) + ", ocean = " + \
                       str(mask_o) + "), it should not be the case"
            break
        arr, areacell, keyerror = Read_data_mask_area(
            var_file, var_name, ref_var["variable_type"], metric, r1, file_area=var_areafile, name_area=var_areaname,
            file_mask=var_landmaskfile, name_mask=var_landmaskname, maskland=mask_l, maskocean=mask_o, debug=debug,
            **kwargs)
        if keyerror is not None:
            break
        method_var += str(counter) + ") read in " + str(r1)
        if mask_l is True:
            method_var += " & land masked"
        if mask_o is True:
            method_var += " & ocean masked"
        # 1.2 compute number of years
        nbr_year = arr.shape[0] // 12
        # 1.3 read time period used
        timebounds = TimeBounds(arr)

        # ------------------------------------------------
        # 2. Preprocess
        # ------------------------------------------------
        my_det = deepcopy(kwargs["detrending"]) if "detrending" in list(kwargs.keys()) else False
        kwargs["detrending"] = False
        my_nor = deepcopy(kwargs["normalization"]) if "normalization" in list(kwargs.keys()) else False
        kwargs["normalization"] = False
        my_smo = deepcopy(kwargs["smoothing"]) if "smoothing" in list(kwargs.keys()) else False
        kwargs["smoothing"] = False
        # 2.1 horizontal average of 'internal_variable_name' in 'r1'
        arr, tmp, keyerror = preprocess_ts_polygon(
            arr, "", areacell=areacell, average="horizontal", region=r1, **kwargs)
        if keyerror is not None:
            break
        for k1 in tmp.split(", "):
            method_var += ";; " + str(counter) + ") " + str(k1)
            counter += 1
        if debug is True:
            dict_debug = {"axes1": str([ax.id for ax in arr.getAxisList()]), "shape1": str(arr.shape),
                          "time1": str(TimeBounds(arr))}
            EnsoErrorsWarnings.debug_mode("\033[92m", "after preprocess_ts_polygon", 15, **dict_debug)
        # 2.2 read, average and remove from 'arr'
        for r2 in reg_list[1:]:
            # 2.2.1 read file and select the 'r2' region
            mask_l = reg_di[r2]["maskland"] if "maskland" in list(reg_di[r2].keys()) else False
            mask_o = reg_di[r2]["maskocean"] if "maskocean" in list(reg_di[r2].keys()) else False
            if mask_l is True and mask_o is True:
                keyerror = "both land and ocean are masked (" + str(r2) + ": land = " + str(mask_l) + ", ocean = " + \
                           str(mask_o) + "), it should not be the case"
                break
            ar2, areacell, keyerror = Read_data_mask_area(
                var_file, var_name, ref_var["variable_type"], metric, r2, file_area=var_areafile,
                name_area=var_areaname, file_mask=var_landmaskfile, name_mask=var_landmaskname, maskland=mask_l,
                maskocean=mask_o, debug=debug, **kwargs)
            if keyerror is not None:
                break
            method_var += ";; " + str(counter) + ") read in " + str(r2)
            if mask_l is True:
                method_var += " & land masked"
            if mask_o is True:
                method_var += " & ocean masked"
            # 2.2.2 horizontal average of 'internal_variable_name' in 'r2'
            ar2, tmp, keyerror = preprocess_ts_polygon(
                ar2, "", areacell=areacell, average="horizontal", region=r2, **kwargs)
            if keyerror is not None:
                break
            if debug is True:
                dict_debug = {"axes1": str([ax.id for ax in ar2.getAxisList()]), "shape1": str(ar2.shape),
                              "time1": str(TimeBounds(ar2))}
                EnsoErrorsWarnings.debug_mode("\033[92m", "after preprocess_ts_polygon " + str(r2), 15, **dict_debug)
            for k1 in tmp.split(", "):
                method_var += " & " + str(k1)
            # 2.2.3 remove ar2 from arr
            arr -= ar2
            method_var += " & removed at each time step"
            counter += 1
            if debug is True:
                dict_debug = {"axes1": str([ax.id for ax in ar2.getAxisList()]), "shape1": str(arr.shape),
                              "time1": str(TimeBounds(arr))}
                EnsoErrorsWarnings.debug_mode("\033[92m", "after arr -= ar2", 15, **dict_debug)
        if keyerror is not None:
            break
        # 2.3 compute anomalies / normalization / detrend / smooth if applicable
        kwargs["detrending"] = deepcopy(my_det)
        kwargs["normalization"] = deepcopy(my_nor)
        kwargs["smoothing"] = deepcopy(my_smo)
        compute_anom = True if sta_ano == "A" else False
        arr, tmp, keyerror = preprocess_ts_polygon(arr, "", compute_anom=compute_anom, **kwargs)
        if keyerror is not None:
            break
        for k1 in tmp.split(", "):
            method_var += ";; " + str(counter) + ") " + str(k1)
            counter += 1
        if debug is True:
            dict_debug = {"axes1": str([ax.id for ax in arr.getAxisList()]), "shape1": str(arr.shape),
                          "time1": str(TimeBounds(arr))}
            EnsoErrorsWarnings.debug_mode("\033[92m", "after preprocess_ts_polygon (ano, det, ...)", 15, **dict_debug)

        # ------------------------------------------------
        # 3. Metric value
        # ------------------------------------------------
        # 3.1 compute metric
        ave, ske, std, var = simple_stats(arr, axis=0)
        # 3.2 compute metric error
        nn = len(arr)
        if statistic == "average":
            mv, mv_error = deepcopy(ave), std / nn ** (1 / 2)
        elif statistic == "skewness":
            mv, mv_error = deepcopy(ske), (6 * nn * (nn - 1) / ((nn - 2) * (nn + 1) * (nn + 3))) ** (1 / 2)
        elif statistic == "standard deviation":
            mv, mv_error = deepcopy(std), std / (2 * (nn - 1)) ** (1 / 2)
        else:
            mv, mv_error = deepcopy(var), var * (2 / (nn - 1)) ** (1 / 2)
        method_var += ";; " + str(counter) + ") compute " + str(statistic)

        # ------------------------------------------------
        # 4. Supplementary dive down diagnostics
        # ------------------------------------------------
        if netcdf is True:
            # supplementary recipes
            dict_nc = {"var1": arr}
            dict_nc["var1_attributes"] = {
                "arrayAVE": ave, "arraySKE": ske, "arraySTD": std, "arrayVAR": var,
                "description": "time series of " + str(var_box) + " averaged " + str(var_sn) + str(sta_ano),
                "number_of_years_used": nbr_year, "time_period": str(timebounds), "units": units1}
            dict_nc["var1_name"] = str(var_sn).lower() + str(sta_ano) + "_" + str(var_box) + "__" + str(dataset)
            dict_nc["var1_time_name"] = "months_" + str(dataset)
            # save netCDF
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + str(metname) + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + str(metname) + ".nc"
            dict1 = {"computation_steps": method_var, "frequency": kwargs["frequency"], "metric_method": method,
                     "metric_name": name, "metric_reference": ref}
            SaveNetcdf(file_name, global_attributes=dict1, **dict_nc)
    # Metric value
    if debug is True:
        dict_debug = {"line1": "metric value: " + str(mv), "line2": "metric value_error: " + str(mv_error)}
        EnsoErrorsWarnings.debug_mode("\033[92m", "end of " + metric, 10, **dict_debug)
    # Create output
    metric_output = {
        "dive_down_diag": dive_down_diag, "keyerror": keyerror, "method": method, "method_detail": method_var,
        "name": name, "nyears": nbr_year, "ref": ref, "time_frequency": kwargs["frequency"], "time_period": timebounds,
        "units": units2, "value": mv, "value_error": mv_error}
    return metric_output
# ---------------------------------------------------------------------------------------------------------------------#
