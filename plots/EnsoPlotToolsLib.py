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
                             str(min(mini)) + ")"]
        if max(tick_labels) < maxi:
            list_strings += [str().ljust(5) + "ticks maximum (" + str(max(tick_labels)) + ") > tab maximum (" +
                             str(min(maxi)) + ")"]
        EnsoErrorsWarnings.MyWarning(list_strings)
    return tick_labels


def read_diag(dict_diag, dict_metric, model, reference, metric_variables):
    diag_mod = dict_diag[model]
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
    metric_value = dict_metric[obs]  # ["ref_" + obs]
    return diag_mod, diag_obs, metric_value, obs


def read_obs(xml, variables_in_xml, metric_variables, varname, dict_metric):
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
    tab_out = xml[newvar]
    metric_value = dict_metric[obs]#["ref_" + obs]
    return tab_out, metric_value, obs


def read_var(var_to_read, filename_nc, model, reference, metric_variables, dict_metric):
    if isinstance(var_to_read, str):
        var_to_read = [var_to_read]
    ff = open_dataset(filename_nc, decode_times=False)
    variables_in_file = sorted([var for var in ff.keys()], key=lambda v: v.upper())
    # read model
    tab_mod = list()
    for var in var_to_read:
        tab_mod.append(ff[var + model])
    # reab obs
    tab_obs = list()
    for var in var_to_read:
        varobs = var + reference
        if varobs in variables_in_file:
            tab = ff[varobs]
            metval = dict_metric[reference]#["ref_" + reference]
            obs = deepcopy(reference)
        else:
            tab, metval, obs = read_obs(ff, variables_in_file, metric_variables, var, dict_metric)
        tab_obs.append(tab)
    ff.close()
    return tab_mod, tab_obs, metval, obs

