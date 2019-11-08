# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from math import ceil as MATHceil
from math import floor as MATHfloor
from numpy import arange as NUMPYarange
from numpy import around as NUMPYaround
from numpy import array as NUMPYarray
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
# xarray based functions
from xarray import open_dataset
# ENSO_metrics functions
from EnsoCollectionsLib import ReferenceObservations
import EnsoErrorsWarnings


calendar_months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
observations = sorted(ReferenceObservations().keys(), key=lambda v: v.upper())


def create_labels(label_name, label_ticks):
    if label_name == "months":
        if len(label_ticks) > 40:
            mult = 6
        elif len(label_ticks) > 10:
            mult = 4
        else:
            mult = 3
        label_ticks = [ii for ii in label_ticks if ii % mult == 0]
        label = [calendar_months[ii % 12] for ii in label_ticks]
    elif label_name == "latitude":
        if len(label_ticks) < 40:
            mult = 10
        else:
            mult = 20
        label_ticks = [ii for ii in label_ticks if ii % mult == 0]
        if min(label_ticks) < 0 and max(label_ticks) > 0 and 0 not in label_ticks:
            label_ticks = NUMPYarray(label_ticks)
            while 0 not in label_ticks:
                label_ticks = label_ticks + 1
        label = [str(abs(int(ii))) + '$^\circ$S' if ii < 0 else (str(abs(int(ii))) + '$^\circ$N' if ii > 0 else 'eq')
                 for ii in label_ticks]
    elif label_name == "longitude":
        if len(label_ticks) < 200:
            mult = 40
        else:
            mult = 90
        label_ticks = [ii for ii in label_ticks if ii % mult == 0]
        if min(label_ticks) < 180 and max(label_ticks) > 180 and 180 not in label_ticks:
            label_ticks = NUMPYarray(label_ticks)
            while 180 not in label_ticks:
                label_ticks = label_ticks + 10
        label = [str(int(ii)) + "$^\circ$E" if ii < 180 else (
            str(abs(int(ii) - 360)) + "$^\circ$W" if ii > 180 else "180$^\circ$") for ii in label_ticks]
    return label_ticks, label


def create_levels(labelbar):
    diff = round(float(labelbar[1] - labelbar[0]), 2)
    if diff in [0.3, 0.6, 0.9] or diff % 3 == 0:
        mult = 3
    elif diff in [0.1, 0.2, 0.4, 0.8, 1.0, 2, 4, 8, 10, 20, 40, 60, 80, 100]:
        mult = 4
    elif diff in [0.5, 5, 25]:
        mult = 5
    else:
        mult = 6
    delta = float(diff) / mult
    return [round(kk + jj * delta, 2) for kk in labelbar[:-1] for jj in range(mult)] + [labelbar[-1]]


def format_metric(metric_type, metric_value, metric_units):
    if metric_type in ["CORR", "RMSE"]:
        mytext = deepcopy(metric_type)
    else:
        if metric_type == "difference":
            mytext = "model-ref"
        elif metric_type == "ratio":
            mytext = r"$\frac{model}{ref}$"
        else:
            mytext = r"$\frac{model-ref}{ref}$"
    return mytext + ": " + "{0:.2f}".format(metric_value) + " " + metric_units


def minmax_plot(tab, metric=False):
    # define minimum and maximum
    tmp = [NUMPYma__masked_invalid(NUMPYarray(tmp)) for tmp in tab]
    tmp = [tt.min() for tt in tmp] + [tt.max() for tt in tmp]
    mini, maxi = min(tmp), max(tmp)
    if mini < 0 and maxi > 0:
        locmaxi = max([abs(mini), abs(maxi)])
        locmini = -deepcopy(locmaxi)
    else:
        locmini, locmaxi = deepcopy(mini), deepcopy(maxi)
    # find the power of ten to get an interval between 1 and 10
    mult = pow(10, int(str("%e" % abs(locmaxi - locmini)).split('e')[1]))
    locmini, locmaxi = int(MATHfloor(float(locmini) / mult)), int(MATHceil(float(locmaxi) / mult))
    if locmaxi == 2 and maxi < 15 and mult == 10 and abs(locmini) != locmaxi:
        locmini, locmaxi = 0, 15
        mult = 1.
    scalmini, scalemaxi = mini / mult, maxi / mult
    interval = locmaxi - locmini
    listbase = list(NUMPYaround([ii*10**exp for exp in range(-1, 1) for ii in range(1, 6)], decimals=1))
    listbase = listbase + listbase
    listmult = [3] * int(len(listbase)/2) + [4] * int(len(listbase)/2)
    list1 = list(NUMPYaround([listbase[ii] * listmult[ii] for ii in range(len(listbase))], decimals=1))
    list2 = list(NUMPYaround([abs(ii - interval) for ii in list1], decimals=1))
    interval = list1[list2.index(min(list2))]
    base = listbase[list1.index(interval)]
    if base * 4.5 < interval:
        ii = 1
        tmp = sorted(list2)
        while base * 4.5 < interval:
            interval = list1[list2.index(tmp[ii])]
            base = listbase[list1.index(interval)]
            ii += 1
    if abs(locmini) == locmaxi:
        maxi_out = 2 * base
        while maxi_out - base > locmaxi:
            maxi_out -= base
        if metric is True and maxi_out < scalemaxi + base * 0.4:
            maxi_out += base
        mini_out = -maxi_out
    else:
        if locmini < 0 and locmaxi <= 0:
            locmini, locmaxi = abs(locmaxi), abs(locmini)
            sign = -1
        else:
            sign = 1
        half_int = int(round(interval / 2.))
        tmp_middle = locmini + half_int
        mini_out = max([0, tmp_middle - half_int])
        while mini_out > locmini:
            mini_out -= base
        while mini_out + base < locmini:
            mini_out += base
        maxi_out = mini_out + 2 * base
        while maxi_out < locmaxi:
            maxi_out += base
        while maxi_out - base > locmaxi:
            maxi_out -= base
        minmax = list(NUMPYaround(NUMPYarray([mini_out, maxi_out]) * sign, decimals=0).astype(int))
        mini_out, maxi_out = min(minmax), max(minmax)
        if metric is True:
            if maxi_out < scalemaxi + base * 0.4:
                maxi_out += base
    tick_labels = NUMPYarange(mini_out, maxi_out + base / 2., base)
    tick_labels = list(NUMPYaround(tick_labels * mult, decimals=4))
    if all([True if int(ii) == float(ii) else False for ii in tick_labels]):
        tick_labels = list(NUMPYaround(tick_labels, decimals=0).astype(int))
    if len(tick_labels) > 6:
        list_strings = [
            "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too many ticks for axis",
            str().ljust(5) + str(len(tick_labels)) + " ticks: " + str(tick_labels),
            str().ljust(5) + "there should not be more than 6"
        ]
        EnsoErrorsWarnings.MyWarning(list_strings)
    if min(tick_labels) > mini or max(tick_labels) < maxi:
        list_strings = ["WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) +
                        ": wrong bounds in ticks for axis"]
        if min(tick_labels) > mini:
            list_strings += [str().ljust(5) + "ticks minimum (" + str(min(tick_labels)) + ") > tab minimum (" +
                             str(mini) + ")"]
        if max(tick_labels) < maxi:
            list_strings += [str().ljust(5) + "ticks maximum (" + str(max(tick_labels)) + ") > tab maximum (" +
                             str(maxi) + ")"]
        EnsoErrorsWarnings.MyWarning(list_strings)
    return tick_labels


def read_diag(dict_diag, dict_metric, model, reference, metric_variables):
    if isinstance(model, str):
        diag_mod = dict_diag[model]
    else:
        diag_mod = [dict_diag[mod] for mod in model]
    if reference in dict_diag.keys():
        obs = deepcopy(reference)
    else:
        if len(metric_variables) == 1:
            for obs in observations:
                if obs in dict_diag.keys():
                    break
        else:
            for obs1 in observations:
                for obs2 in observations:
                    obs = obs1 + "_" + obs2
                    if obs in dict_diag.keys():
                        break
    diag_obs = dict_diag[obs]
    if isinstance(model, str):
        metric_value = dict_metric[obs][model]
    else:
        metric_value = [dict_metric[obs][mod] for mod in model]
    # metric_value = dict_metric[obs]  # ["ref_" + obs]
    return diag_mod, diag_obs, metric_value, obs


def read_obs(xml, variables_in_xml, metric_variables, varname, dict_metric, model):
    if len(metric_variables) == 1:
        for obs in observations:
            newvar = varname + obs
            if newvar in variables_in_xml:
                break
    else:
        for obs1 in observations:
            for obs2 in observations:
                obs = obs1 + "_" + obs2
                newvar = varname + obs
                if newvar in variables_in_xml:
                    break
    # if "_lon__" in newvar or "_hov__" in newvar:
    #     list_strings = [
    #         "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": reader",
    #         str().ljust(5) + str(newvar) + " trick: read only between 140E and 96W",
    #         str().ljust(5) + "this should not stay like that!!!"
    #     ]
    #     EnsoErrorsWarnings.MyWarning(list_strings)
    #     tab_out = xml[newvar].sel(longitude=slice(140, 264))
    # else:
    #     tab_out = xml[newvar]
    tab_out = xml[newvar]
    metric_value = dict_metric[obs][model]#["ref_" + obs]
    return tab_out, metric_value, obs


def reader(filename_nc, model, reference, var_to_read, metric_variables, dict_metric):
    ff = open_dataset(filename_nc, decode_times=False)
    variables_in_file = sorted([var for var in ff.keys()], key=lambda v: v.upper())
    # read model
    tab_mod = list()
    for var in var_to_read:
        # if "_lon__" in var or "_hov__" in var:
        #     list_strings = [
        #         "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": reader",
        #         str().ljust(5) + str(var) + " trick: read only between 140E and 96W",
        #         str().ljust(5) + "this should not stay like that!!!"
        #     ]
        #     EnsoErrorsWarnings.MyWarning(list_strings)
        #     tab_mod.append(ff[var + model].sel(longitude=slice(140, 264)))
        # else:
        #     tab_mod.append(ff[var + model])
        tab_mod.append(ff[var + model])
    # reab obs
    tab_obs = list()
    for var in var_to_read:
        varobs = var + reference
        if varobs in variables_in_file:
            # if "_lon__" in varobs or "_hov__" in varobs:
            #     list_strings = [
            #         "WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": reader",
            #         str().ljust(5) + str(varobs) + " trick: read only between 140E and 96W",
            #         str().ljust(5) + "this should not stay like that!!!"
            #     ]
            #     EnsoErrorsWarnings.MyWarning(list_strings)
            #     tab = ff[varobs].sel(longitude=slice(140, 264))
            # else:
            #     tab = ff[varobs]
            tab = ff[varobs]
            metval = dict_metric[reference][model]  # ["ref_" + reference]
            obs = deepcopy(reference)
        else:
            tab, metval, obs = read_obs(ff, variables_in_file, metric_variables, var, dict_metric, model)
        tab_obs.append(tab)
    ff.close()
    return tab_mod, tab_obs, metval, obs


def read_var(var_to_read, filename_nc, model, reference, metric_variables, dict_metric):
    if isinstance(var_to_read, str):
        var_to_read = [var_to_read]
    if isinstance(filename_nc, str):
        tab_mod, tab_obs, metval, obs = reader(filename_nc, model, reference, var_to_read, metric_variables,
                                               dict_metric)
    else:
        tab_mod, metval, obs = list(), list(), list()
        for ii in range(len(filename_nc)):
            tmp1, tab_obs, tmp2, tmp3 = reader(filename_nc[ii], model[ii], reference, var_to_read, metric_variables,
                                               dict_metric)
            tab_mod.append(tmp1)
            metval.append(tmp2)
            obs.append(tmp3)
        obs = list(set(obs))
        if len(obs) > 1:
            list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too many obs",
                            str().ljust(5) + "var_to_read = "+str(var_to_read),
                            str().ljust(5) + "filename_nc = " + str(filename_nc),
                            str().ljust(5) + "model = " + str(model),
                            str().ljust(5) + "reference = " + str(reference)]
            EnsoErrorsWarnings.MyError(list_strings)
        else:
            obs = obs[0]
    return tab_mod, tab_obs, metval, obs

