# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from math import ceil as MATHceil
from math import floor as MATHfloor
from numpy import arange as NUMPYarange
from numpy import around as NUMPYaround
from numpy import array as NUMPYarray
from numpy import isnan as NUMPYisnan
from numpy import mean as NUMPYmean
from numpy import nan as NUMPYnan
from numpy import sort as NUMPYsort
from numpy import where as NUMPYwhere
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.random import randint as NUMPYrandom__randint
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile

# xarray based functions
from xarray import open_dataset
from xarray import where as XARRAYwhere

# ENSO_metrics functions
from EnsoMetrics.EnsoCollectionsLib import ReferenceObservations
from EnsoMetrics.EnsoPlotLib import plot_param
from EnsoMetrics import EnsoErrorsWarnings


calendar_months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
observations = sorted(list(ReferenceObservations().keys()), key=lambda v: v.upper())

# metrics order
metrics_background = [
    "BiasPrLatRmse", "BiasPrLonRmse", "BiasSshLatRmse", "BiasSshLonRmse", "BiasSstLatRmse", "BiasSstLonRmse",
    "BiasTauxLatRmse", "BiasTauxLonRmse", "SeasonalPrLatRmse", "SeasonalPrLonRmse", "SeasonalSshLatRmse",
    "SeasonalSshLonRmse", "SeasonalSstLatRmse", "SeasonalSstLonRmse", "SeasonalTauxLatRmse", "SeasonalTauxLonRmse"]
metrics_basic = [
    "EnsoSstLonRmse", "EnsoPrTsRmse", "EnsoSstTsRmse", "EnsoTauxTsRmse", "EnsoAmpl", "EnsoSeasonality", "EnsoSstSkew",
    "EnsoDuration", "EnsoSstDiversity", "EnsoSstDiversity_1", "EnsoSstDiversity_2", "NinoSstDiversity",
    "NinoSstDiversity_1", "NinoSstDiversity_2"]
metrics_teleconnection = [
    "EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoPrMapStd", "EnsoPrMapDjfCorr", "EnsoPrMapDjfRmse", "EnsoPrMapDjfStd",
    "EnsoPrMapJjaCorr", "EnsoPrMapJjaRmse", "EnsoPrMapJjaStd", "EnsoSlpMapCorr", "EnsoSlpMapRmse", "EnsoSlpMapStd",
    "EnsoSlpMapDjfCorr", "EnsoSlpMapDjfRmse", "EnsoSlpMapDjfStd", "EnsoSlpMapJjaCorr", "EnsoSlpMapJjaRmse",
    "EnsoSlpMapJjaStd", "EnsoSstMapCorr", "EnsoSstMapRmse", "EnsoSstMapStd", "EnsoSstMapDjfCorr", "EnsoSstMapDjfRmse",
    "EnsoSstMapDjfStd", "EnsoSstMapJjaCorr", "EnsoSstMapJjaRmse", "EnsoSstMapJjaStd"]
metrics_process = [
    "EnsoFbSstTaux", "EnsoFbTauxSsh", "EnsoFbSshSst", "EnsoFbSstThf", "EnsoFbSstSwr", "EnsoFbSstLhf", "EnsoFbSstLwr",
    "EnsoFbSstShf", "EnsodSstOce", "EnsodSstOce_1", "EnsodSstOce_2"]
# models order
models_order = [
    "ACCESS1-0", "ACCESS1-3", "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM1-1", "BCC-CSM1-1-M", "BCC-CSM2-MR", "BCC-ESM1",
    "BNU-ESM", "CAMS-CSM1-0", "CanCM4", "CanESM2", "CanESM5", "CanESM5-CanOE", "CCSM4", "CESM1-BGC", "CESM1-CAM5",
    "CESM2", "CESM2-FV2", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2-WACCM", "CESM2-WACCM-FV2", "CMCC-CESM", "CMCC-CM",
    "CMCC-CMS", "CNRM-CM5", "CNRM-CM5-2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "CSIRO-Mk3-6-0",
    "CSIRO-Mk3L-1-2", "E3SM-1-0", "E3SM-1-1", "EC-EARTH", "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", "FGOALS-g2",
    "FGOALS-s2", "FIO-ESM", "GFDL-CM2p1", "GFDL-CM3", "GFDL-CM4", "GFDL-ESM2G", "GFDL-ESM2M", "GFDL-ESM4",
    "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-1-H", "GISS-E2-R", "GISS-E2-R-CC", "HadCM3",
    "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "HadGEM3-GC31-LL", "INMCM4", "INM-CM4-8", "INM-CM5-0", "IPSL-CM5A-LR",
    "IPSL-CM5A-MR", "IPSL-CM5B-LR", "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC4h", "MIROC5", "MIROC6", "MIROC-ESM",
    "MIROC-ESM-CHEM", "MIROC-ES2L", "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P", "MPI-ESM-1-2-HAM", "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR", "MRI-CGCM3", "MRI-ESM1", "MRI-ESM2-0", "NESM3", "NorESM1-M", "NorESM1-ME", "NorCPM1", "NorESM2-LM",
    "NorESM2-MM", "SAM0-UNICON", "TaiESM1", "UKESM1-0-LL"]


def bootstrap(tab, num_samples=1000000, alpha=0.05, nech=None, statistic=NUMPYmean):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    n = len(tab)
    if nech is None:
        nech = deepcopy(n)
    idx = NUMPYrandom__randint(0, n, (num_samples, nech))
    samples = tab[idx]
    stat = NUMPYsort(statistic(samples, 1))
    return [stat[int(round((alpha / 2.)*num_samples))], stat[int(round((1 - alpha / 2.)*num_samples))]]


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
        if len(label_ticks) > 80:
            label_ticks = [-60, 0, 60]
        else:
            label_ticks = [ii for ii in label_ticks if ii % mult == 0]
        if min(label_ticks) < 0 and max(label_ticks) > 0 and 0 not in label_ticks:
            label_ticks = NUMPYarray(label_ticks)
            while 0 not in label_ticks:
                label_ticks = label_ticks + 1
        label = [str(abs(int(ii))) + "S" if ii < 0 else (str(abs(int(ii))) + "N" if ii > 0 else "eq")
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
        label = [str(int(ii)) + "E" if ii < 180 else (
            str(abs(int(ii) - 360)) + "W" if ii > 180 else "180") for ii in label_ticks]
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


def create_lines(line_i, line_o=[], threshold=50):
    if len(line_i) > 1:
        words = line_i.split(" ")
        tmp = ""
        for ii, elt in enumerate(words):
            tmp += elt
            if ii != len(words) - 1:
                if len(tmp + " " + str(words[ii + 1])) > threshold:
                    line_o.append(tmp)
                    tmp = ""
                else:
                    tmp += " "
            else:
                line_o.append(tmp)
        del tmp, words
    return line_o


def create_round_string(number):
    if number < 1.:
        nbr_out = "{0:.2f}".format(round(number, 2))
    elif number < 10.:
        nbr_out = "{0:.1f}".format(round(number, 1))
    else:
        nbr_out = str(round(number, 0))
    return nbr_out


def format_metric(metric_type, metric_value, metric_units):
    if metric_type in ["CORR", "RMSE"]:
        mytext = deepcopy(metric_type)
    else:
        if metric_type == "difference":
            mytext = "model-ref"
        elif metric_type == "ratio":
            mytext = r"$\frac{model}{ref}$"
        elif metric_type == "relative_difference":
            mytext = r"$\frac{model-ref}{ref}$"
        else:
            mytext = r"$abs\left(\frac{model-ref}{ref}\right)$"
    if metric_value is not None:
        return mytext + ": " + "{0:.2f}".format(metric_value) + " " + metric_units
    else:
        return None


def get_reference(metric_collection, metric):
    if metric_collection in ["ENSO_tel"] and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    return plot_param(metric_collection, my_met)['metric_reference']


def minimaxi(tab):
    tmp = [my_mask(tmp, remove_masked=True) for tmp in tab]
    tmp = [tt.min() for tt in tmp] + [tt.max() for tt in tmp]
    return min(tmp), max(tmp)


def minmax_plot(tab, metric=False):
    # define minimum and maximum
    mini, maxi = minimaxi(tab)
    if mini == maxi or abs(maxi-mini)/float(abs(maxi+mini))<1e-2:
        tt = max(abs(mini), abs(maxi)) / 10.
        tmp = int(str("%.e" % tt)[3:])
        if mini == maxi or (mini < 0 and maxi > 0):
            tmp = 10**-tmp if tt < 1 else 10**tmp
        elif mini >= 0:
            tmp = 10**-tmp if tt < 1 else 10**tmp
        else:
            tmp = 10**-tmp if tt < 1 else 10**tmp
        mini = 0 if mini > 0 and (mini-tmp)<0 else mini - tmp
        maxi = 0 if maxi < 0 and (maxi+tmp)>0 else maxi + tmp
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
    listmult = [3] * int(round(len(listbase) / 2.)) + [4] * int(round(len(listbase) / 2.))
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
    if len(tick_labels) > 7:
        list_strings = [
            "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many ticks for axis",
            str().ljust(5) + str(len(tick_labels)) + " ticks: " + str(tick_labels),
            str().ljust(5) + "there should not be more than 7"
        ]
        EnsoErrorsWarnings.my_warning(list_strings)
    if min(tick_labels) > mini or max(tick_labels) < maxi:
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                        ": wrong bounds in ticks for axis"]
        if min(tick_labels) > mini:
            list_strings += [str().ljust(5) + "ticks minimum (" + str(min(tick_labels)) + ") > tab minimum (" +
                             str(mini) + ")"]
        if max(tick_labels) < maxi:
            list_strings += [str().ljust(5) + "ticks maximum (" + str(max(tick_labels)) + ") > tab maximum (" +
                             str(maxi) + ")"]
        EnsoErrorsWarnings.my_warning(list_strings)
    return tick_labels


def my_average(tab, axis=None, remove_masked=False):
    tmp = my_mask(tab, remove_masked=remove_masked)
    return NUMPYmean(tmp, axis=axis)


def my_bootstrap(tab1, tab2):
    mea1 = float(NUMPYarray(tab1).mean())
    mea2 = float(NUMPYarray(tab2).mean())
    bst1 = bootstrap(NUMPYarray(tab1), nech=len(tab2))
    bst2 = bootstrap(NUMPYarray(tab2), nech=len(tab1))
    return bst1, bst2, mea1, mea2


def my_legend(modname, obsname, filename_nc, models2=None, member=None, plot_metric=True, shading=False):
    legend = ["ref: " + obsname]
    if isinstance(filename_nc, str) is True or isinstance(filename_nc, str) is True:
        if isinstance(modname, str) is True or isinstance(modname, str) is True:
            if isinstance(member, str) is True or isinstance(member, str) is True:
                legend += [modname + " " + member]
            else:
                legend += [modname]
        else:
            legend += ["model"]
    else:
        if shading is True and (isinstance(filename_nc, dict) is True or plot_metric is False):
            if isinstance(member, list) is True and len(member) == len(modname):
                legend += [mod.upper() + mem + " (" + str(len(models2[mod])) + ")"
                           for mod, mem in zip(modname, member)]
            else:
                legend += [mod.upper() + "(" + str(len(models2[mod])) + ")" for mod in modname]
        elif shading is True and plot_metric is True:
            if isinstance(member, list) is True and len(member) == len(modname):
                legend += [mod.upper() + mem for mod, mem in zip(modname, member)]
            else:
                legend += [mod.upper() for mod in modname]
        else:
            if isinstance(member, list) is True and len(member) == len(modname):
                legend += [mod + mem for mod, mem in zip(modname, member)]
            else:
                legend += modname
    return legend


def my_mask(tab, remove_masked=False):
    tmp = NUMPYarray(tab, dtype=float)
    tmp = NUMPYwhere(tmp == None, NUMPYnan, tmp)
    tmp = NUMPYma__masked_invalid(tmp)
    if remove_masked is True:
        # tmp = tmp[~tmp.mask]
        tmp = tmp.compressed()
    return tmp


def my_mask_map(tab_ref, tab_mod):
    return XARRAYwhere(NUMPYisnan(tab_ref), NUMPYnan, tab_mod)


def read_diag(dict_diag, dict_metric, model, reference, metric_variables, shading=False, member=None):
    if member is not None:
        modelKeyName = model + "_" + member
    else:
        modelKeyName = model
    if isinstance(model, str):
        diag_mod = dict_diag[modelKeyName]
    else:
        if shading is True:
            diag_mod =\
                [[dict_diag[mod][mm] for mm in sorted(list(dict_diag[mod].keys()), key=lambda v: v.upper())] for mod in modelKeyName]
        else:
            diag_mod = [dict_diag[mod] for mod in modelKeyName]
    if shading is True:
        my_ref = list(dict_diag["obs"].keys())
    else:
        my_ref = list(dict_diag.keys())
    if reference in my_ref:
        obs = deepcopy(reference)
    else:
        if len(metric_variables) == 1:
            for obs1 in observations:
                if obs1 in my_ref:
                    obs = deepcopy(obs1)
                    break
        else:
            for obs1 in observations:
                for obs2 in observations:
                    obs3 = obs1 + "_" + obs2
                    if obs3 in my_ref:
                        obs = deepcopy(obs3)
                        break
        try: obs
        except: obs = sorted(my_ref)[0]
    if shading is True:
        diag_obs = dict_diag["obs"][obs]
    else:
        diag_obs = dict_diag[obs]
    if isinstance(model, str):
        metric_value = dict_metric[obs][model]
    else:
        if shading is True:
            metric_value = [
                my_average([dict_metric[mod][obs][mm]
                            for mm in sorted(list(dict_metric[mod][obs].keys()), key=lambda v: v.upper())],
                           remove_masked=True)
                for mod in model]
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
        my_break = False
        for obs1 in observations:
            for obs2 in observations:
                obs = obs1 + "_" + obs2
                newvar = varname + obs
                if newvar in variables_in_xml:
                    my_break = True
                    break
            if my_break is True:
                break
    # if "_lon__" in newvar or "_hov__" in newvar:
    #     list_strings = [
    #         "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": reader",
    #         str().ljust(5) + str(newvar) + " trick: read only between 140E and 96W",
    #         str().ljust(5) + "this should not stay like that!!!"
    #     ]
    #     EnsoErrorsWarnings.my_warning(list_strings)
    #     tab_out = xml[newvar].sel(longitude=slice(140, 264))
    # else:
    #     tab_out = xml[newvar]
    tab_out = xml[newvar]
    metric_value = dict_metric[obs][model]#["ref_" + obs]
    return tab_out, metric_value, obs


def reader(filename_nc, model, reference, var_to_read, metric_variables, dict_metric, member=None, met_in_file=False,
           met_type=None, met_pattern=""):
    ff = open_dataset(filename_nc, decode_times=False)
    variables_in_file = sorted([var for var in list(ff.keys())], key=lambda v: v.upper())
    # read model
    tab_mod = list()
    for var in var_to_read:
        # if "_lon__" in var or "_hov__" in var:
        #     list_strings = [
        #         "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": reader",
        #         str().ljust(5) + str(var) + " trick: read only between 140E and 96W",
        #         str().ljust(5) + "this should not stay like that!!!"
        #     ]
        #     EnsoErrorsWarnings.my_warning(list_strings)
        #     tab_mod.append(ff[var + model].sel(longitude=slice(140, 264)))
        # else:
        #     tab_mod.append(ff[var + model])
        varName_in_nc = var + model
        if member is not None:
            varName_in_nc += "_" + member 
        #tab_mod.append(ff[var + model])
        tab_mod.append(ff[varName_in_nc])
    # reab obs
    tab_obs = list()
    for var in var_to_read:
        varobs = var + reference
        if varobs in variables_in_file:
            # if "_lon__" in varobs or "_hov__" in varobs:
            #     list_strings = [
            #         "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": reader",
            #         str().ljust(5) + str(varobs) + " trick: read only between 140E and 96W",
            #         str().ljust(5) + "this should not stay like that!!!"
            #     ]
            #     EnsoErrorsWarnings.my_warning(list_strings)
            #     tab = ff[varobs].sel(longitude=slice(140, 264))
            # else:
            #     tab = ff[varobs]
            tab = ff[varobs]
            metval = dict_metric[reference][model]  # ["ref_" + reference]
            obs = deepcopy(reference)
        else:
            tab, metval, obs = read_obs(ff, variables_in_file, metric_variables, var, dict_metric, model)
        tab_obs.append(tab)
    if isinstance(var_to_read, list) is True and len(var_to_read) == 1:
        if met_in_file is True:
            if isinstance(met_type, str):
                for key in list(ff.attrs.keys()):
                    if met_type + "_" + obs + "_" + met_pattern == key:
                        val = ff.attrs[key]
                try: val
                except: val = None
                metval = deepcopy(val)
                del val
            elif isinstance(met_type, list):
                metval = list()
                for mety in met_type:
                    for key in list(ff.attrs.keys()):
                        if mety + "_" + obs + "_" + met_pattern == key or (met_pattern == "" and mety + "_" + obs == key):
                            val = ff.attrs[key]
                    try: val
                    except: val = None
                    metval.append(val)
                    del val
    elif isinstance(var_to_read, list) is True and len(var_to_read) == 2 and\
            ("nina" in var_to_read[0] or "nino" in var_to_read[0]):
        metval = list()
        for var in var_to_read:
            add = "nina" if "nina" in var else "nino"
            if met_in_file is True:
                if isinstance(met_type, str):
                    for key in list(ff.attrs.keys()):
                        if met_type + "_" + obs + "_" + add + "_" + met_pattern == key:
                            val = ff.attrs[key]
                    try:    val
                    except: val = None
                    metval.append(val)
                    del val
                elif isinstance(met_type, list):
                    tmpval = list()
                    for mety in met_type:
                        for key in list(ff.attrs.keys()):
                            if mety + "_" + obs + "_" + add + "_" + met_pattern == key or (
                                    met_pattern == "" and mety + "_" + obs + "_" + add == key):
                                val = ff.attrs[key]
                        try: val
                        except: val = None
                        tmpval.append(val)
                        del val
                    metval.append(tmpval)
                    del tmpval
            del add
    ff.close()
    return tab_mod, tab_obs, metval, obs


def read_var(var_to_read, filename_nc, model, reference, metric_variables, dict_metric, models2=None, member=None,
             shading=False, met_in_file=False, met_type=None, met_pattern=""):
    if isinstance(var_to_read, str):
        var_to_read = [var_to_read]
    if isinstance(filename_nc, str):
        tab_mod, tab_obs, metval, obs = reader(filename_nc, model, reference, var_to_read, metric_variables,
                                               dict_metric, member=member, met_in_file=met_in_file, met_type=met_type,
                                               met_pattern=met_pattern)
    else:
        tab_mod, metval, obs = list(), list(), list()
        if shading is True:
            for jj in range(len(model)):
                tmp1, tmp2 = list(), list()
                for ii in range(len(filename_nc[model[jj]])):
                    ttt1, tab_obs, ttt2, tmp3 = reader(
                        filename_nc[model[jj]][ii], models2[model[jj]][ii], reference, var_to_read, metric_variables,
                        dict_metric[model[jj]], member=member, met_in_file=met_in_file, met_type=met_type, met_pattern=met_pattern)
                    tmp1.append(ttt1)
                    tmp2.append(ttt2)
                    obs.append(tmp3)
                tab_mod.append(tmp1)
                metval.append(tmp2)
        else:
            for ii in range(len(filename_nc)):
                tmp1, tab_obs, tmp2, tmp3 = reader(filename_nc[ii], model[ii], reference, var_to_read, metric_variables,
                                                   dict_metric, member=member)
                tab_mod.append(tmp1)
                metval.append(tmp2)
                obs.append(tmp3)
        obs = list(set(obs))
        if len(obs) > 1:
            list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many obs",
                            str().ljust(5) + "var_to_read = "+str(var_to_read),
                            str().ljust(5) + "filename_nc = " + str(filename_nc),
                            str().ljust(5) + "model = " + str(model),
                            str().ljust(5) + "reference = " + str(reference)]
            EnsoErrorsWarnings.my_error(list_strings)
        else:
            obs = obs[0]
    return tab_mod, tab_obs, metval, obs


def shading_levels(tab, lev=[5, 25, 75, 95], axis=None):
    return [SCIPYstats__scoreatpercentile(tab, ll, axis=axis) for ll in lev] + [my_average(tab, axis=axis)]


def remove_metrics(metrics_in, metric_collection, reduced_set=False, portraitplot=False):
    """
    Removes some metrics from given list

    Inputs:
    ------
    :param metrics_in: list of string
        List of metrics.
    :param metric_collection: string
        Name of a metric collection.
    **Optional arguments:**
    :param reduced_set: boolean, optional
        True to remove extra metrics that are not in the final set chosen by CLIVAR PRP.
        If set to False it removes metrics that are in more than one metric collection.
        Default value is False
    :param portraitplot: boolean, optional
        True to remove extra metrics that are not in the final set chosen by CLIVAR PRP but keep metrics that are in
        more than one metric collection.
        If set to False it removes metrics that are in more than one metric collection.
        Default value is False

    Output:
    ------
    :return metrics_out: list of string
        Input list of metrics minus some metrics depending on given metric collection.
    """
    metrics_out = deepcopy(metrics_in)
    if reduced_set is True:
        if portraitplot is True:
            if metric_collection == "ENSO_perf":
                to_remove = [
                    'BiasSshLatRmse', 'BiasSshLonRmse', 'BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse',
                    'EnsoSstDiversity_1', 'EnsoTauxTsRmse', 'NinaSstDur', 'NinaSstDur_1',
                    'NinaSstDur_2', 'NinaSstLonRmse', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse',
                    'NinaSstTsRmse_1', 'NinaSstTsRmse_2', 'NinoSstDiversity', 'NinoSstDiversity_1',
                    'NinoSstDiversity_2', 'NinoSstDur', 'NinoSstDur_1', 'NinoSstDur_2', 'NinoSstLonRmse',
                    'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse', 'NinoSstTsRmse_1', 'NinoSstTsRmse_2',
                    "SeasonalSshLatRmse", "SeasonalSshLonRmse", "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
            elif metric_collection == "ENSO_proc":
                to_remove = [
                    'BiasSshLonRmse', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf', 'EnsoFbSstSwr']
            else:
                to_remove = [
                    'EnsoPrMapCorr', 'EnsoPrMapRmse', 'EnsoPrMapStd',
                    'EnsoPrMapDjfCorr', 'EnsoPrMapDjfStd', 'EnsoPrMapJjaCorr', 'EnsoPrMapJjaStd', 'EnsoSlpMapCorr',
                    'EnsoSlpMapRmse', 'EnsoSlpMapStd', 'EnsoSlpMapDjfCorr', 'EnsoSlpMapDjfRmse', 'EnsoSlpMapDjfStd',
                    'EnsoSlpMapJjaCorr', 'EnsoSlpMapJjaRmse', 'EnsoSlpMapJjaStd', 'EnsoSstMapCorr', 'EnsoSstMapRmse',
                    'EnsoSstMapStd', 'EnsoSstMapDjfCorr', 'EnsoSstMapDjfStd', 'EnsoSstMapJjaCorr', 'EnsoSstMapJjaStd',
                    'NinaPrMapCorr', 'NinaPrMap_1Corr', 'NinaPrMap_2Corr', 'NinaPrMapRmse', 'NinaPrMap_1Rmse',
                    'NinaPrMap_2Rmse', 'NinaPrMapStd', 'NinaPrMap_1Std', 'NinaPrMap_2Std', 'NinaSlpMapCorr',
                    'NinaSlpMap_1Corr', 'NinaSlpMap_2Corr', 'NinaSlpMapRmse', 'NinaSlpMap_1Rmse', 'NinaSlpMap_2Rmse',
                    'NinaSlpMapStd', 'NinaSlpMap_1Std', 'NinaSlpMap_2Std', 'NinaSstLonRmse', 'NinaSstLonRmse_1',
                    'NinaSstLonRmse_2', 'NinaSstMapCorr', 'NinaSstMap_1Corr', 'NinaSstMap_2Corr', 'NinaSstMapRmse',
                    'NinaSstMap_1Rmse', 'NinaSstMap_2Rmse', 'NinaSstMapStd', 'NinaSstMap_1Std', 'NinaSstMap_2Std',
                    'NinoPrMapCorr', 'NinoPrMap_1Corr', 'NinoPrMap_2Corr', 'NinoPrMapRmse', 'NinoPrMap_1Rmse',
                    'NinoPrMap_2Rmse', 'NinoPrMapStd', 'NinoPrMap_1Std', 'NinoPrMap_2Std', 'NinoSlpMapCorr',
                    'NinoSlpMap_1Corr', 'NinoSlpMap_2Corr', 'NinoSlpMap_1Rmse', 'NinoSlpMap_2Rmse', 'NinoSlpMapStd',
                    'NinoSlpMap_1Std', 'NinoSlpMap_2Std', 'NinoSstLonRmse', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                    'NinoSstMapCorr', 'NinoSstMap_1Corr', 'NinoSstMap_2Corr', 'NinoSstMapRmse', 'NinoSstMap_1Rmse',
                    'NinoSstMap_2Rmse', 'NinoSstMapStd', 'NinoSstMap_1Std', 'NinoSstMap_2Std']
        else:
            if metric_collection == "ENSO_perf":
                to_remove = [
                    'BiasSshLatRmse', 'BiasSshLonRmse', 'BiasSstLatRmse', 'BiasTauxLatRmse', 'EnsoPrTsRmse',
                    'EnsoSstDiversity_1', 'EnsoTauxTsRmse', 'NinaSstDur', 'NinaSstDur_1',
                    'NinaSstDur_2', 'NinaSstLonRmse', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse',
                    'NinaSstTsRmse_1', 'NinaSstTsRmse_2', 'NinoSstDiversity', 'NinoSstDiversity_1',
                    'NinoSstDiversity_2', 'NinoSstDur', 'NinoSstDur_1', 'NinoSstDur_2', 'NinoSstLonRmse',
                    'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse', 'NinoSstTsRmse_1', 'NinoSstTsRmse_2',
                    "SeasonalSshLatRmse", "SeasonalSshLonRmse", "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
            elif metric_collection == "ENSO_proc":
                to_remove = [
                    'BiasSshLonRmse', 'BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality',
                    'EnsoSstLonRmse', 'EnsoSstSkew', 'EnsodSstOce_1', 'EnsoFbSstLhf', 'EnsoFbSstLwr', 'EnsoFbSstShf',
                    'EnsoFbSstSwr']
            else:
                to_remove = [
                    'EnsoAmpl', 'EnsoSeasonality', 'EnsoSstLonRmse', 'EnsoPrMapCorr', 'EnsoPrMapRmse', 'EnsoPrMapStd',
                    'EnsoPrMapDjfCorr', 'EnsoPrMapDjfStd', 'EnsoPrMapJjaCorr', 'EnsoPrMapJjaStd', 'EnsoSlpMapCorr',
                    'EnsoSlpMapRmse', 'EnsoSlpMapStd', 'EnsoSlpMapDjfCorr', 'EnsoSlpMapDjfRmse', 'EnsoSlpMapDjfStd',
                    'EnsoSlpMapJjaCorr', 'EnsoSlpMapJjaRmse', 'EnsoSlpMapJjaStd', 'EnsoSstMapCorr', 'EnsoSstMapRmse',
                    'EnsoSstMapStd', 'EnsoSstMapDjfCorr', 'EnsoSstMapDjfStd', 'EnsoSstMapJjaCorr', 'EnsoSstMapJjaStd',
                    'NinaPrMapCorr', 'NinaPrMap_1Corr', 'NinaPrMap_2Corr', 'NinaPrMapRmse', 'NinaPrMap_1Rmse',
                    'NinaPrMap_2Rmse', 'NinaPrMapStd', 'NinaPrMap_1Std', 'NinaPrMap_2Std', 'NinaSlpMapCorr',
                    'NinaSlpMap_1Corr', 'NinaSlpMap_2Corr', 'NinaSlpMapRmse', 'NinaSlpMap_1Rmse', 'NinaSlpMap_2Rmse',
                    'NinaSlpMapStd', 'NinaSlpMap_1Std', 'NinaSlpMap_2Std', 'NinaSstLonRmse', 'NinaSstLonRmse_1',
                    'NinaSstLonRmse_2', 'NinaSstMapCorr', 'NinaSstMap_1Corr', 'NinaSstMap_2Corr', 'NinaSstMapRmse',
                    'NinaSstMap_1Rmse', 'NinaSstMap_2Rmse', 'NinaSstMapStd', 'NinaSstMap_1Std', 'NinaSstMap_2Std',
                    'NinoPrMapCorr', 'NinoPrMap_1Corr', 'NinoPrMap_2Corr', 'NinoPrMapRmse', 'NinoPrMap_1Rmse',
                    'NinoPrMap_2Rmse', 'NinoPrMapStd', 'NinoPrMap_1Std', 'NinoPrMap_2Std', 'NinoSlpMapCorr',
                    'NinoSlpMap_1Corr', 'NinoSlpMap_2Corr', 'NinoSlpMap_1Rmse', 'NinoSlpMap_2Rmse', 'NinoSlpMapStd',
                    'NinoSlpMap_1Std', 'NinoSlpMap_2Std', 'NinoSstLonRmse', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2',
                    'NinoSstMapCorr', 'NinoSstMap_1Corr', 'NinoSstMap_2Corr', 'NinoSstMapRmse', 'NinoSstMap_1Rmse',
                    'NinoSstMap_2Rmse', 'NinoSstMapStd', 'NinoSstMap_1Std', 'NinoSstMap_2Std']
    else:
        if portraitplot is True:
            to_remove = []
        else:
            if metric_collection == "ENSO_perf":
                to_remove = ['BiasSshLatRmse', 'BiasSshLonRmse', "SeasonalSshLatRmse", "SeasonalSshLonRmse"]
            elif metric_collection == "ENSO_proc":
                to_remove = ['BiasSshLonRmse', 'BiasSstLonRmse', 'BiasTauxLonRmse', 'EnsoAmpl', 'EnsoSeasonality',
                             'EnsoSstLonRmse', 'EnsoSstSkew']
            else:
                to_remove = ['EnsoAmpl', 'EnsoSeasonality', 'EnsoSstLonRmse', 'NinaSstLonRmse', 'NinaSstLonRmse_1',
                             'NinaSstLonRmse_2', 'NinoSstLonRmse', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2']
    for met in to_remove:
        while met in metrics_out:
            metrics_out.remove(met)
    return metrics_out


def find_first_member(members_in):
    """
    Finds first member name

    Inputs:
    ------
    :param members_in: list of string
        list of member names (e.g., "r1i1p1", "r1i1p2")

    Output:
    ------
    :return member_out: string
        first member of the given list
    """
    if "r1i1p1" in members_in:
        member_out = "r1i1p1"
    elif "r1i1p1f1" in members_in:
        member_out = "r1i1p1f1"
    elif "r1i1p1f2" in members_in:
        member_out = "r1i1p1f2"
    else:
        member_out = sort_members(members_in)[0]
    return member_out


def return_metrics_type():
    return metrics_background, metrics_basic, metrics_teleconnection, metrics_process


def sort_members(members_in):
    """
    Finds first member name

    Inputs:
    ------
    :param members_in: list of string
        list of member names (e.g., "r1i1p1", "r1i1p2")

    Output:
    ------
    :return members_out: list of string
        given list of member names sorted
    """
    members_tmp = list()
    for mem in members_in:
        mem2 = mem.replace("r1i", "r01i").replace("r2i", "r02i").replace("r3i", "r03i").replace("r4i", "r04i")
        mem2 = mem2.replace("r5i", "r05i").replace("r6i", "r06i").replace("r7i", "r07i").replace("r8i", "r08i")
        mem2 = mem2.replace("r9i", "r09i")
        mem2 = mem2.replace("i1p", "i01p").replace("i2p", "i02p").replace("i3p", "i03p").replace("i4p", "i04p")
        mem2 = mem2.replace("i5p", "i05p").replace("i6p", "i06p").replace("i7p", "i07p").replace("i8p", "i08p")
        mem2 = mem2.replace("i9p", "i09p")
        if "f" in mem2:
            mem2 = mem2.replace("p1f", "p01f").replace("p2f", "p02f").replace("p3f", "p03f").replace("p4f", "p04f")
            mem2 = mem2.replace("p5f", "p05f").replace("p6f", "p06f").replace("p7f", "p07f").replace("p8f", "p08f")
            mem2 = mem2.replace("p9f", "p09f")
            mem2 = mem2.replace("f1", "f01").replace("f2", "f02").replace("f3", "f03").replace("f4", "f04")
            mem2 = mem2.replace("f5", "f05").replace("f6", "f06").replace("f7", "f07").replace("f8", "f08")
            mem2 = mem2.replace("f9", "f09")
        else:
            mem2 = mem2.replace("p1", "p01").replace("p2", "p02").replace("p3", "p03").replace("p4", "p04")
            mem2 = mem2.replace("p5", "p05").replace("p6", "p06").replace("p7", "p07").replace("p8", "p08")
            mem2 = mem2.replace("p9", "p09")
        members_tmp.append(mem2)
    members_tmp = sorted(list(set(members_tmp)), key=lambda v: v.upper())
    members_out = list()
    for mem in members_tmp:
        mem2 = mem.replace("r0", "r").replace("i0", "i").replace("p0", "p").replace("f0", "f")
        members_out.append(mem2)
    return members_out


def sort_metrics(metrics_in):
    """
    Puts given list of models in a certain order (putting together model generations)

    Input:
    -----
    :param metrics_in: list of string
        List of metrics.

    Output:
    ------
    :return metrics_out: list of string
        Input list of metrics reordered
    """
    # metric type and order
    met_order = metrics_background + metrics_basic + metrics_teleconnection + metrics_process
    metrics_out = sorted(list(set(metrics_in)), key=lambda v: v.upper())
    metrics_out = [met for met in met_order if met in metrics_out]
    metrics_out += sorted(list(set(metrics_in) - set(metrics_out)), key=lambda v: v.upper())
    return metrics_out


def sort_models(models_in):
    """
    Puts given list of models in a certain order (putting together model generations)

    Input:
    -----
    :param models_in: list of string
        List of models.

    Output:
    ------
    :return models_out: list of string
        Input list of models reordered
    """
    # model order
    models_out = sorted(list(set(models_in)), key=lambda v: v.upper())
    models_out = [mod for mod in models_order if mod in models_out]
    models_out += sorted(list(set(models_in) - set(models_out)), key=lambda v: v.upper())
    return models_out
