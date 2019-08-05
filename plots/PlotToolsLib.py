# -*- coding:UTF-8 -*-

# ---------------------------------------------------#
# Import the right package
# ---------------------------------------------------#
from copy import deepcopy
from inspect import stack as INSPECTstack
from math import ceil as MATHceil
from math import floor as MATHfloor
from numpy import arange as NUMPYarange
from numpy import around as NUMPYaround
from numpy import array as NUMPYarray
# ENSO_metrics functions
import EnsoErrorsWarnings
from EnsoCollectionsLib import ReferenceObservations

observations = sorted(ReferenceObservations().keys(), key=lambda v: v.upper())


def format_scale(scale):
    if scale < 0.1 or scale > 10:
        scale = "{0:.0e}".format(scale)
    elif scale == 0.1:
        scale = "{0:.1f}".format(scale)
    else:
        scale = str(int(scale))
    return scale


def minmax_plot(tab, metric=False):
    # define minimum and maximum
    mini, maxi = min(tab), max(tab)
    if realmini < 0 and maxi > 0:
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


def read_obs(xml, variables_in_xml, metrics_variables, varname, model, dict_metric):
    if len(metrics_variables) == 1:
        for obs in observations:
            newvar = varname.replace(model, obs)
            if newvar in variables_in_xml:
                break
    else:
        for obs1 in observations:
            for obs2 in observations:
                obs = obs1 + "_" + obs2
                newvar = varname.replace(model, obs)
                if newvar in variables_in_xml:
                    break
    tab_out = xml[newvar]
    metric_value = dict_metric["ref_" + obs]["value"]
    return tab_out, metric_value