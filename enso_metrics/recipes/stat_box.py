# -*- coding:UTF-8 -*-


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
from typing import Union

# local functions
from enso_metrics.tools.default import default_arg_values, input_dictionary_formater, set_default_str
from enso_metrics.tools import prints
from enso_metrics.wrapper import processors
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
def stat_box(
        dict_input: dict[
            str, dict[
                str, Union[int, float, str, list[str], None, dict[
                    str, Union[int, float, str, list[str], None]]]]],
        dataset: str = "unknown",
        experiment: str = "unknown",
        member: str = "unknown",
        netcdf: bool = False,
        netcdf_filename: str = None,
        project: str = "unknown",
        region: str = None,
        statistic: str = None,
        variable: str = None,
        **kwargs):
    # get default recipe, regions and variables parameters if they are not given
    # needed_kwarg = ["recipe_param", "regions_param", "variables_param"]
    needed_kwarg = ["regions_param", "variables_param"]
    for k in needed_kwarg:
        if k not in list(kwargs.keys()):
            kwargs[k] = default_arg_values(k)
    # check region, statistic and variable values as they are required
    region = set_default_str(region, list(kwargs["regions_param"].keys()), "nino3.4")
    list_statistics = ["average", "skewness", "standard_deviation", "variance"]
    statistic = set_default_str(statistic, list_statistics, "average")
    variable = set_default_str(variable, list(kwargs["variables_param"].keys()), "ts")
    # check / format dict_input
    dict_input = input_dictionary_formater(dict_input)
    # fake loop to be able to break out if an error occurs
    for _ in range(1):
        # ------------------------------------------------
        # 1. Read files
        # ------------------------------------------------
        # 1.1 Create a dictionary with variables and files
        dict_data = processors.reader(dict_input, variable, **kwargs)
        print(list(dict_data.keys()))
        # processors
        dict_processors = {
            "1__variable_name": {
                "name_out": "new_name",
                "to_do": {
                    "1__masker": "do",
                    "2__selecter": "do",
                    "3__averager": "do",
                    "4__detrender": "do",
                    "5__anomaler": "do",
                    "6__seasonal_cycler": "do",
                    "7__normalizer": "do",
                    "8__smoother": "do",
                },
            },
        }
        stop
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
        r1 = reg_list[0]
        
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
        lwr_mod, areacell_mod, keyerror_mod = Read_data_mask_area_multifile(
            lwrfilemod, lwrnamemod, "heat flux", "lwr", metric, box, file_area=lwrareafilemod, name_area=lwrareanamemod,
            file_mask=lwrlandmaskfilemod, name_mask=lwrlandmasknamemod, maskland=True, maskocean=False,
            time_bounds=kwargs["time_bounds_mod"], debug=debug, interpreter="project_interpreter_mod_var1", **kwargs)
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
            prints.debug_mode("\033[92m", "after preprocess_ts_polygon", 15, **dict_debug)
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
                prints.debug_mode("\033[92m", "after preprocess_ts_polygon " + str(r2), 15, **dict_debug)
            for k1 in tmp.split(", "):
                method_var += " & " + str(k1)
            # 2.2.3 remove ar2 from arr
            arr -= ar2
            method_var += " & removed at each time step"
            counter += 1
            if debug is True:
                dict_debug = {"axes1": str([ax.id for ax in ar2.getAxisList()]), "shape1": str(arr.shape),
                              "time1": str(TimeBounds(arr))}
                prints.debug_mode("\033[92m", "after arr -= ar2", 15, **dict_debug)
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
            prints.debug_mode("\033[92m", "after preprocess_ts_polygon (ano, det, ...)", 15, **dict_debug)

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
        prints.debug_mode("\033[92m", "end of " + metric, 10, **dict_debug)
    # Create output
    metric_output = {
        "dive_down_diag": dive_down_diag, "keyerror": keyerror, "method": method, "method_detail": method_var,
        "name": name, "nyears": nbr_year, "ref": ref, "time_frequency": kwargs["frequency"], "time_period": timebounds,
        "units": units2, "value": mv, "value_error": mv_error}
    return metric_output
# ---------------------------------------------------------------------------------------------------------------------#
