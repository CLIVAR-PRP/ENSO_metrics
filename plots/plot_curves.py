# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
# plot curves for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right package and initialization of CVS
# ---------------------------------------------------#
from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from math import ceil as MATHceil
from math import floor as MATHfloor
import matplotlib.pyplot as plt
from numpy import array as NUMPYarray
from os.path import join as OSpath__join
from re import search as REsearch
from sys import argv as SYSargv
from xarray import open_dataset
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection, ReferenceObservations


# ---------------------------------------------------#
# colors for printing
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ---------------------------------------------------#


# ---------------------------------------------------#
# Variables initialization
# ---------------------------------------------------#
print(bcolors.OKGREEN + """ 
Needed argument(s):
    1) metric collection: ENSO_perf, ENSO_tel, MC1,...
    2) experiment:        historical, piControl,...
    3) metric:            one of the metric of the given metric collection
Optional argument(s):
By default:
python -i zzz_resave_netcdf.py ENSO_perf historical EnsoAmpl
""" + bcolors.ENDC)
print('')
# ---------------------------------------------------#
# axes label
seasons_1m = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
# references
ref_colors = {'model': 'r', 'obs': 'k'}
ref_obs = {'ssh': 'AVISO', 'pr': 'GPCPv2.3', 'sst': 'HadISST', 'lhf': 'Tropflux', 'lwr': 'Tropflux', 'shf': 'Tropflux',
           'swr': 'Tropflux', 'taux': 'Tropflux', 'thf': 'Tropflux'}
# metric collections
list_MC = sorted(list(defCollection().keys()), key=lambda v: v.upper())
list_obs = sorted(list(ReferenceObservations().keys()), key=lambda v: v.upper())
# model
model = "CNRM-CM5"
# path
path_main = '/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_05_report/Data/v20190515'
path_plot = '/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_05_report/Wiki'
# ---------------------------------------------------#
# get arguments
try:
    SYSargv[1]
except:
    print('     no argument given')
    print('')
    metric_col = 'ENSO_perf'  # 'CM6_perf'
    experiment = 'historical'
else:
    metric_col = str(SYSargv[1])
    experiment = str(SYSargv[2])
# ---------------------------------------------------#
# print work to do
print(bcolors.OKGREEN + '%%%%%     -----     %%%%%')
print(str().ljust(5) + "create curve plots for:")
print(str().ljust(10) + 'metric collection: ' + str(metric_col))
print(str().ljust(10) + 'experiment: ' + str(experiment))
print('%%%%%     -----     %%%%%' + bcolors.ENDC)
for ii in range(3): print('')


# ---------------------------------------------------#
# functions
def create_minmax_plot(tab):
    mini, maxi = min(tab), max(tab)
    print(mini, maxi)
    if mini < 0 and maxi > 0:
        locmaxi = max([abs(mini), abs(maxi)])
        locmini = -deepcopy(locmaxi)
    else:
        locmaxi = deepcopy(maxi)
        locmini = deepcopy(mini)
    interval = locmaxi - locmini
    tmp = str("%e" % abs(interval))
    exp = int(tmp.split('e')[1])
    mult = pow(10, exp)
    locmini, locmaxi = float(locmini) / mult, float(locmaxi) / mult
    interval = float(interval) / mult
    listbase = [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 4, 3, 5]
    list1 = [round(base * 4, 1) if base < 1 else int(round(base * 4, 0)) for base in listbase]
    print(list1)
    list2 = [abs(ii - interval) for ii in list1]
    tmp = sorted(list2)
    interval = list1[list2.index(tmp[0])]
    base = listbase[list1.index(interval)]
    ii = 0
    while base * 5 < locmaxi - locmini:
        interval = list1[list2.index(tmp[ii])]
        base = listbase[list1.index(interval)]
        ii += 1
    print(interval, base)
    if abs(locmini) == locmaxi:
        mini_out = -((interval / 2.) + base)
        maxi_out = (interval / 2.) + base
        for ii in range(3):
            if mini_out < locmini and mini_out + base < locmini:
                mini_out += base
                maxi_out -= base
    else:
        tmp_middle = locmini + (locmaxi - locmini) / 2.
        print("middle = " + str(tmp_middle))
        tmp2 = str("%e" % abs(tmp_middle))
        exp2 = int(tmp2.split('e')[1])
        mul2 = pow(10, exp2)
        print(tmp2, exp2, mul2)
        if mul2 >= 10:
            mul2 = mul2 / 10.
        tmp_middle = round(tmp_middle / mul2) * mul2
        print("middle = " + str(tmp_middle))
        if mini >= 0:
            mini_out = max([0, tmp_middle - interval / 2. - base])
        else:
            mini_out = tmp_middle - interval / 2. - base
        for ii in range(4):
            if mini_out > locmini:
                mini_out -= base
        for ii in range(4):
            if mini_out < locmini and mini_out + base < locmini:
                mini_out += base
        if mini >= 0:
            if mini_out == 0:
                maxi_out = mini_out + 2 * base
            else:
                maxi_out = tmp_middle + interval / 2. + base
        else:
            maxi_out = min([0, tmp_middle + interval / 2. + base])
        print(maxi_out)
        for ii in range(4):
            if maxi_out < locmaxi:
                maxi_out += base
                print(maxi_out)
        for ii in range(4):
            if maxi_out > locmaxi and maxi_out - base > locmaxi:
                maxi_out -= base
                print(maxi_out)
        if maxi_out == 0:
            mini_out = maxi_out - 2 * base
            for ii in range(4):
                if mini_out > locmini:
                    mini_out -= base
            for ii in range(4):
                if mini_out < locmini and mini_out + base < locmini:
                    mini_out += base
        if maxi_out < locmaxi + base * 0.75:
            maxi_out += base
    mini_out = round(mini_out, 5)
    maxi_out = round(maxi_out, 5)
    if abs(mini_out) < 1 and abs(maxi_out) < 1:
        mini_out = mini_out * 10
        maxi_out = maxi_out * 10
        mult = mult / 10
    if abs(mini_out) > 10 and abs(maxi_out) > 10:
        mini_out = mini_out / 10
        maxi_out = maxi_out / 10
        mult = mult * 10
    if mini >= 0 and maxi < 40 and mult > 1:
        mini_out = mini_out * 10
        maxi_out = maxi_out * 10
        mult = mult / 10
    if ((mini >= 0 and maxi > 1 and maxi < 10) or mini > -5 and maxi < 5) and mult < 1:
        mini_out = mini_out / 10
        maxi_out = maxi_out / 10
        mult = mult * 10
    print([mini_out, maxi_out], mult)
    return [mini_out, maxi_out], mult


def format_scale(scale):
    if scale < 0.1 or scale > 10:
        scale = "{0:.0e}".format(scale)
    elif scale == 0.1:
        scale = "{0:.1f}".format(scale)
    else:
        scale = str(int(scale))
    return scale


def plot_curves(tab_mod, tab_obs, name_png, title='', xname='', yname='', units='', mytext=[]):
    axisname = list(tab_mod.coords.keys())[0]
    axis = list(NUMPYarray(tab_mod.coords[axisname]))
    tab_mod = NUMPYarray(tab_mod)
    tab_obs = NUMPYarray(tab_obs)
    # figure initialization
    if (axisname == "month" or axisname == "months" or axisname == "time") and len(tab_mod) == 72:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
    # title
    ax.set_title(title, fontsize=20, y=1.01, loc='left')
    # x axis
    tmp = list(range(int(MATHfloor(min(axis))), int(MATHceil(max(axis))) + 1))
    plt.xlim(min(axis), max(axis))
    if axisname == "month" or axisname == "months" or axisname == "time":
        if len(tmp) > 40:
            mult = 6
        else:
            mult = 4
        tmp = [ii for ii in tmp if ii % mult == 0]
        label = [seasons_1m[ii % 12] for ii in tmp]
        if xname == '':
            xname = 'months'
    elif "lat" in axisname:
        if len(tmp) < 40:
            mult = 10
        else:
            mult = 20
        tmp = [ii for ii in tmp if ii % mult == 0]
        if min(tmp) < 0 and max(tmp) > 0 and 0 not in tmp:
            tmp = NUMPYarray(tmp)
            while 0 not in tmp:
                tmp = tmp + 1
        label = [str(abs(int(ii))) + '$^\circ$S' if ii < 0 else (str(abs(int(ii))) + '$^\circ$N' if ii > 0 else 'eq')
                 for ii in tmp]
        if xname == '':
            xname = 'latitude'
    elif "lon" in axisname:
        if len(tmp) < 200:
            mult = 40
        else:
            mult = 90
        tmp = [ii for ii in tmp if ii % mult == 0]
        if min(tmp) < 180 and max(tmp) > 180 and 180 not in tmp:
            tmp = NUMPYarray(tmp)
            while 180 not in tmp:
                tmp = tmp + 10
        label = [str(int(ii)) + '$^\circ$E' if ii < 180 else (
            str(abs(int(ii) - 360)) + '$^\circ$W' if ii > 180 else '180$^\circ$') for ii in tmp]
        if xname == '':
            xname = 'longitude'
    plt.xticks(tmp, label)
    ax.set_xlabel(xname, fontsize=20)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    # y axis
    minmax, mult = create_minmax_plot([tab_mod.min(), tab_mod.max(), tab_obs.min(), tab_obs.max()])
    if mult != 1.:
        tab_mod = tab_mod / mult
        tab_obs = tab_obs / mult
        tmp = units.split(" ")
        if len(tmp) == 1:
            scale = format_scale(mult)
            name_units = deepcopy(units)
        else:
            scale = format_scale(float(tmp[0]) * mult)
            name_units = deepcopy(tmp[1])
        ylabel = yname.replace(units, scale + " " + name_units)
    else:
        ylabel = deepcopy(yname)
    plt.ylim(min(minmax), max(minmax))
    plt.locator_params(axis='y', nbins=4)
    ax.set_ylabel(ylabel, fontsize=20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    if min(minmax) < 0 and max(minmax) > 0:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    # plot curves
    ax.plot(axis, list(tab_mod), lw=4, label='model', color=ref_colors['model'])
    ax.plot(axis, list(tab_obs), lw=4, label='obs', color=ref_colors['obs'])
    # relative space
    x1 = min(axis);
    x2 = max(axis)
    dx = (x2 - x1) / 100.
    y1 = min(minmax);
    y2 = max(minmax)
    dy = (y2 - y1) / 100.
    # legend
    leg_text = ['obs', 'model']
    for ii, txt in enumerate(leg_text):
        plt.text((2 * dx) + x1, y2 - (((ii + 1) * 6) * dy), txt, fontsize=15, color=ref_colors[txt],
                 horizontalalignment='left', verticalalignment='center')
    # my text
    for ii, txt in enumerate(mytext):
        plt.text(x2 - (2 * dx), y2 - (((ii + 1) * 6) * dy), txt, fontsize=15, color='k', horizontalalignment='right',
                 verticalalignment='center')
    plt.grid(linestyle='--', linewidth=1, which='major')
    plt.savefig(name_png, bbox_inches='tight')
    plt.close()
    return


def read_obs(xml, variables_in_xml, metrics_variables, varname, dict_metric):
    if len(metrics_variables) == 1:
        for obs in list_obs:
            newvar = varname.replace(model, obs)
            if newvar in variables_in_xml:
                break
    else:
        for obs1 in list_obs:
            for obs2 in list_obs:
                obs = obs1 + "_" + obs2
                newvar = varname.replace(model, obs)
                if newvar in variables_in_xml:
                    break
    tab_out = xml[newvar]
    metric_value = dict_metric["ref_" + obs]["value"]
    return tab_out, metric_value


# ---------------------------------------------------#


# ---------------------------------------------------#
# main
# ---------------------------------------------------#
no_var = ['months', 'bounds_months', 'bounds_time', 'latitude', 'bounds_latitude', 'bounds_latitude_a',
          'bounds_latitude_b', 'longitude', 'bounds_longitude', 'bounds_longitude_a', 'bounds_longitude_b',
          'bounds_years', 'bounds_years_a', 'bounds_years_b', 'bounds_years_c', 'bounds_years_d', 'bounds_years_e']
path_js = OSpath__join(path_main, "JSONs")
path_nc = OSpath__join(path_main, metric_col + "_cmip5")
dict_MC = defCollection(metric_col)['metrics_list']
# ---------------------------------------------------#
print('curves')
# read data
dict1 = dict()
if ' ':
    # json file
    file_js = sorted(
        list(GLOBiglob(OSpath__join(path_js, "cmip?_" + experiment + "_" + metric_col + "_v2019????.json"))),
        key=lambda v: v.upper())
    if "cmip5" in path_nc:
        file_js = file_js[0]
    else:
        file_js = file_js[1]
    with open(file_js) as ff:
        data_json = json.load(ff)
    data_json = data_json["RESULTS"]["model"][model]
    ff.close()
    del ff, file_js
    files_nc = sorted(list(
        GLOBiglob(OSpath__join(path_nc, "cmip?_" + experiment + "_" + metric_col + "_v2019????_" + model + "_*.nc"))),
                      key=lambda v: v.upper())
    for file1 in files_nc:
        # find metric name
        metric = REsearch("_v2019(.+?)_" + model + "_(.+?).nc", file1)
        metric = metric.group(2)
        if 'SstDiv' in metric or 'SstDur' in metric:
            pass
        else:
            # find metric observations and variable(s)
            met_obs = dict_MC[metric]["obs_name"]
            met_var = dict_MC[metric]["variables"]
            print(str().ljust(10) + metric)
            # open file
            ff = open_dataset(file1, decode_times=False)
            # variables to read
            var_in = sorted([var for var in list(ff.keys()) if var not in no_var and ff[var].ndim == 1],
                            key=lambda v: v.upper())
            var_to_read = [var for var in var_in if model in var]
            for nbr, var in enumerate(var_to_read):
                # read model
                tab1 = ff[var]
                # reab obs
                if len(met_var) == 1:
                    obsname = ref_obs[met_var[0]]
                    yname = met_var[0].upper()
                else:
                    obsname = ref_obs[met_var[0]] + "_" + ref_obs[met_var[1]]
                    yname = ''
                varobs = var.replace(model, obsname)
                if varobs in var_in:
                    tab2 = ff[varobs]
                    metval = data_json["value"][metric]["metric"]["ref_" + obsname]["value"]
                else:
                    tab2, metval = read_obs(ff, var_in, met_var, var, data_json["value"][metric]["metric"])
                # plot
                metuni = data_json["metadata"]["metrics"][metric]["metric"]["units"].replace("C", "$^\circ$C")
                if "Corr" in metric:
                    mytext = "CORR"
                elif "Rmse" in metric:
                    mytext = "RMSE"
                elif "STD" in metric:
                    mytext = r"$\frac{model}{obs}$"
                else:
                    mettyp = dict_MC[metric]['metric_computation']
                    if mettyp == "difference":
                        mytext = "model-obs"
                    elif mettyp == "ratio":
                        mytext = r"$\frac{model}{obs}$"
                    else:
                        mytext = r"$\frac{model-obs}{obs}$"
                mytext = [mytext + ": " + "{0:.2f}".format(metval)]
                if "Corr" not in metric:
                    mytext = [mytext[0] + " " + metuni]
                units = tab1.units.replace("C", "$^\circ$C")
                yname = yname + " (" + units + ")"
                name_png = OSpath__join(path_plot, metric + "_" + model + "_" + str(nbr + 1).zfill(2))
                plot_curves(tab1, tab2, name_png, title=metric, xname='', yname=yname, units=units, mytext=mytext)
                del metuni, metval, mytext, name_png, obsname, tab1, tab2, units, varobs, yname
            ff.close()
            del ff, met_obs, met_var, var_in, var_to_read
        del metric
    del data_json, files_nc
