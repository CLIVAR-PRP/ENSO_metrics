# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Library for the main driver to plot outputs of the ENSO_metrics package
#---------------------------------------------------#


#---------------------------------------------------#
# import python packages
# usual python package
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
import math
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from numpy import array as NUMPYarray
from numpy import arange as NUMPYarange
# CDAT package
from MV2 import array as MV2array
# ENSO_metrics package
from EnsoCollectionsLib import ReferenceObservations, ReferenceRegions
from EnsoErrorsWarnings import MyError
#---------------------------------------------------#


dict_colors = {'ERA-Interim': 'green', 'ERSSTv5': 'skyblue', 'GPCPv2.3': 'orange', 'HadISST': 'orange',
               'Tropflux': 'purple', 'model': 'black'}


def isobs(name):
    # list observations
    list_observations = sorted(ReferenceObservations().keys())
    # test if the given name contains the name of an observational dataset
    isobs = False
    for obs in list_observations:
        if obs in name:
            isobs = True
            break
    return isobs


def create_box(region):
    dict_reg = ReferenceRegions(region)
    lat, lon = dict_reg['latitude'], dict_reg['longitude']
    maxlat, minlat, maxlon, minlon = max(lat), min(lat), max(lon), min(lon)
    li_type, li_width, li_x1x2, li_y1y2 = ['solid' for ii in range(4)], [4 for ii in range(4)], [[minlon,maxlon],[minlon,maxlon],[minlon,minlon],[maxlon,maxlon]], [[maxlat,maxlat],[minlat,minlat],[minlat,maxlat],[minlat,maxlat]]
    return li_type, li_width, li_x1x2, li_y1y2


def create_dom(tab):
    min1 = min(tab.getAxis(0)[:])
    max1 = max(tab.getAxis(0)[:])
    min2 = min(tab.getAxis(1)[:])
    max2 = max(tab.getAxis(1)[:])
    return [min1,max1,min2,max2]


def create_minmax_plot(tab):
    locmini, locmaxi = deepcopy(min(tab)), deepcopy(max(tab))
    print "in = " + str([locmini, locmaxi])
    if locmini < 0 and locmaxi > 0:
        locmaxi = max([abs(locmini), abs(locmaxi)])
        locmini = -deepcopy(locmaxi)
    interval = locmaxi - locmini
    tmp = str("%e" % abs(interval))
    exp = int(tmp.split('e')[1])
    mult = pow(10, exp)
    locmini, locmaxi = float(locmini) / mult, float(locmaxi) / mult
    interval = float(interval) / mult
    listbase = [0.1, 0.2, 0.4, 0.5, 1, 2, 4, 5]
    list1 = [round(base * 6, 1) if base < 1 else int(round(base * 6, 0)) for base in listbase]
    list2 = [abs(ii - interval) for ii in list1]
    interval = list1[list2.index(min(list2))]
    base = listbase[list1.index(interval)]
    print "base = " +str(base)
    if abs(locmini) == locmaxi:
        mini_out = (interval / 2.) - base
        maxi_out = (interval / 2.) + base
    else:
        mini_out = max([0, (round(locmini + (locmaxi - locmini) / 2.) - interval / 2.) - base])
        maxi_out = (round(locmini + (locmaxi - locmini) / 2.) + interval / 2.) + base
    if exp in [-1,0,1]:
        mini_out = mini_out * mult
        maxi_out = maxi_out * mult
        mult = 1.
    print "out = " + str([mini_out, maxi_out])
    return mini_out, maxi_out, mult

def create_label(mini, maxi, ratio=1.):
    """
    Computes range and create label for plots

    :param minmax: list of float
        minimum and maximum values of the array

    return: label: list of float
    """
    if mini<0 and maxi>0:
        maxi = max([abs(mini),abs(maxi)])
        mini = -deepcopy(maxi)
    interval = (maxi-mini)*ratio
    tmp = str("%e" % abs(interval))
    exp = int(tmp.split('e')[1])
    mult = pow(10,exp)
    locmini, locmaxi = float(mini)/mult, float(maxi)/mult
    interval = float(interval)/mult
    listbase = [0.1,0.2,0.4,0.5,1,2,4,5]
    list1 = [round(base*6,1) if base<1 else int(round(base*6,0)) for base in listbase]
    list2 = [abs(ii-interval) for ii in list1]
    interval = list1[list2.index(min(list2))]
    base = listbase[list1.index(interval)]
    if abs(locmini)==locmaxi:
        label = [ii*base-(interval/2.) for ii in range(7)]
    else:
        label = [ii*base+(round(locmini+(locmaxi-locmini)/2.)-interval/2.) for ii in range(6)]
    if exp in [-1,0,1]:
        base = float(base)*mult
        label = NUMPYarray(label)*mult
        mult = 1.
    tmp2 = str("%e" % base)
    exp2 = int(tmp2.split('e')[1])
    label = [round(ii,2) if exp2==-2 else (round(ii,1) if exp2==-1 else int(ii)) for ii in label]
    return label, mult


def create_label_and_color(name, minmax, ratio=1.):
    if 'corr' in name.lower():
        colormap = 'bl_to_darkred'
        label, mult = MV2array([round(ii,1) for ii in NUMPYarange(-1.0,1.0+0.4,0.4)]), 1.
    elif 'std' in name.lower():
        colormap = 'ltbl_to_drkbl'
        label, mult = create_label(0, max(minmax), ratio=ratio)
    else:
        colormap = 'bl_to_darkred'
        label, mult = create_label(min(minmax), max(minmax), ratio=ratio)
    return label, colormap, mult


def myround(x, base=5, side='up', delta=None):
    """ Rounds x using the given base and side (up or down)
        e.g., 13 rounded up   on base 5 is 15
              13 rounded down on base 5 is 10
              -0.11 rounded up   on base 0.05 is -0.1
              -0.11 rounded down on base 0.05 is -0.15
        Returns: x rounded

        x     : float, number to round
        base  : (optional) float, base on which x is rounded
        side  : (optional) string, side on which x is rounded (up or down)
        delta : (optional) float, a better round is adapted on delta (for intervals)
    """
    if delta is None: delta = x
    jj = 10 ** (math.floor(math.log10(abs(delta))))
    x = x / jj
    if x < 0:
        sign, x2 = -1, -x
    else:
        sign, x2 = 1, x
    if (side == 'up' and sign > 0) or (side == 'down' and sign < 0):
        return sign * jj * round(((x2 // (base / jj)) + 1) * (base / jj), 1)
    elif (side == 'up' and sign < 0) or (side == 'down' and sign > 0):
        return sign * jj * round(((x2 // (base / jj))) * (base / jj), 1)


def obs_or_mod(name, list_observations, list_models):
    for obs in list_observations:
        if obs in name:
            name_out = 'obs_' + str(obs)
            break
    for mod in list_models + ['model']:
        if mod in name:
            name_out = 'mod_' + str(mod)
            break
    try:
        name_out
    except:
        list_strings = ['ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
                        'unknown name (obs or model): ' + str(name), 'known observations: ' + str(list_observations),
                        'known models: ' + str(list_models)]
        MyError(list_strings)
    return name_out


def replace_in_obsname(name):
    return deepcopy(name).split('__')[1]


def is_nina_and_nino(list_var):
    tmp = ''
    for var in list_var:
        if 'nina' in var.lower():
            tmp += 'nina'
        elif 'nino' in var.lower():
            tmp += 'nino'
    if 'nina' in tmp and 'nino' in tmp:
        value = True
    else:
        value = False
    return value


# ---------------------------------------------------------------------------------------------------------------------#
# boxplot
def get_nbr_years(name):
    return int(name.split('years_')[0].split('_')[-1])


def read_json(json_file, key_val):
    with open(json_file) as ff:
        data = json.load(ff)
    list_ref = sorted(data[key_val].keys())
    metric = dict()
    for ref in list_ref:
        val = data[key_val][ref]
        if len(val) == 1:
            metric[ref] = val[0]
        else:
            metric[ref] = val
        del val
    return metric


def metric_boxplot(json_pattern, key_val, output_name, title, yname):
    # list all files
    list_files = sorted(list(GLOBiglob(json_pattern + ".json")))
    # nbr years used
    list_nbr = [get_nbr_years(file1) for file1 in list_files]
    # read json
    files_val = dict()
    for nbr, file1 in zip(list_nbr, list_files):
        files_val[nbr] = read_json(file1, key_val)
    # reorganize
    list_ref = sorted(files_val[list_nbr[0]].keys())
    if 'metric' in key_val:
        dict_out = dict((ref, [files_val[nbr][ref] for nbr in list_nbr]) for ref in list_ref)
    else:
        dict_out = dict((ref, [files_val[nbr][ref] for nbr in list_nbr]) if isobs(ref) is False else
                        (ref, files_val[list_nbr[0]][ref]) for ref in list_ref)
    # plot
    alignment = {'horizontalalignment': 'right', 'verticalalignment': 'baseline'}
    font0 = FontProperties()
    font1 = font0.copy()
    font1.set_size('large')
    boxprops = dict(linestyle='-', linewidth=2, color='k')
    capprops = dict(linestyle='-', linewidth=2, color='k')
    marprops = dict(markerfacecolor='k', marker='o', markersize=2.0)
    meaprops = dict(markerfacecolor='r', marker='D', markersize=8.0)
    medprops = dict(linestyle='-', linewidth=2, color='k')
    wisprops = dict(linestyle='-', linewidth=2, color='k')
    tab = list()
    for ref in list_ref:
        if isinstance(dict_out[ref], float):
            tab.append(dict_out[ref])
        else:
            for elt in dict_out[ref]:
                if isinstance(elt, float):
                    tab.append(elt)
                else:
                    tab += deepcopy(elt)
    mini, maxi, mult = create_minmax_plot(tab)
    if mult != 1.:
        title = title + "(*" + str("%e" % mult) + ")"
        for ref in list_ref:
            if isinstance(dict_out[ref], float):
                dict_out[ref] = dict_out[ref] / mult
            else:
                dict_out[ref] = NUMPYarray(dict_out[ref]) / mult
    if 'metric' in key_val:
        for ref in list_ref:
            fig1, ax1 = plt.subplots()
            plt.ylim(mini, maxi)
            plt.title(title, fontsize=20)
            plt.xlabel('number of simulated years', fontsize=15)
            plt.ylabel(yname, fontsize=15)
            for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontsize(15)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(15)
            plt.text(0.98, 0.95, ref, fontproperties=font1, transform=ax1.transAxes, **alignment)
            ax1.boxplot(dict_out[ref], whis=[5, 95], labels=list_nbr, showmeans=True, boxprops=boxprops,
                        capprops=capprops, flierprops=marprops, meanprops=meaprops, medianprops=medprops,
                        whiskerprops=wisprops)
            plt.savefig(output_name + '_' + ref)
            plt.close()
    else:
        dict_ref = dict()
        for ref in list_ref:
            if isinstance(dict_out[ref], float):
                dict_ref[ref] = [dict_out[ref], dict_out[ref]]
            else:
                tab_to_plot = dict_out[ref]
        fig1, ax1 = plt.subplots()
        plt.ylim(mini, maxi)
        plt.title(title, fontsize=20)
        plt.xlabel('number of simulated years', fontsize=15)
        plt.ylabel(yname, fontsize=15)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label.set_fontsize(15)
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(15)
        ax1.boxplot(tab_to_plot, whis=[5, 95], labels=list_nbr, showmeans=True, boxprops=boxprops, capprops=capprops,
                    flierprops=marprops, meanprops=meaprops, medianprops=medprops, whiskerprops=wisprops)
        nref = sorted(dict_ref.keys())
        if len(nref) > 0:
            lines, names = list(), list()
            for ref in nref:
                plt.plot([-100, 100], dict_ref[ref], dict_colors[ref], lw=2)
                lines.append(Line2D([0], [0], color=dict_colors[ref], lw=2))
                names.append(ref)
            plt.legend(lines, names)
        plt.savefig(output_name)
        plt.close()
    return
# ---------------------------------------------------------------------------------------------------------------------#
