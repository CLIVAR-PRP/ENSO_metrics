# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Library for the main driver to plot outputs of the ENSO_metrics package
#---------------------------------------------------#


#---------------------------------------------------#
# import python packages
# usual python package
import calendar
import collections
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
import locale
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from numpy import array as NUMPYarray
from numpy import arange as NUMPYarange
from numpy import isnan as NUMPYisnan
from numpy import meshgrid as NUMPYmeshgrid
from numpy import nan as NUMPYnan
from numpy import ones as NUMPYones
from numpy import product as NUMPYproduct
from numpy import where as NUMPYwhere
from numpy.ma.core import MaskedArray as NUMPYma__core__MaskedArray
from scipy.stats import scoreatpercentile as SCIPYstats__scoreatpercentile
# CDAT package
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import createVariable as CDMS2createVariable
from cdms2 import open as CDMS2open
from MV2 import array as MV2array
from MV2 import average as MV2average
from MV2 import masked_where as MV2masked_where
from MV2 import maximum as MV2maximum
from MV2 import minimum as MV2minimum
from MV2 import zeros as MV2zeros
# ENSO_metrics package
from EnsoCollectionsLib import ReferenceObservations, ReferenceRegions
from EnsoErrorsWarnings import MyError
#---------------------------------------------------#


dict_colors = {'ERA-Interim': 'green', 'ERSSTv5': 'skyblue', 'GPCPv2.3': 'orange', 'HadISST': 'orange',
               'Tropflux': 'purple', 'model': 'black'}

dimensions = ['bounds_latitude', 'bounds_latitude_a', 'bounds_latitude_b', 'bounds_longitude', 'bounds_longitude_a',
              'bounds_longitude_b', 'bounds_months', 'bounds_time', 'bounds_years', 'bounds_years_a', 'bounds_years_b',
              'bounds_years_c', 'bounds_years_d',' bounds_years_e', 'latitude', 'latitude_a', 'latitude_b', 'longitude',
              'longitude_a', 'longitude_b', 'months', 'years', 'years_a', 'years_b', 'years_c', 'years_d', 'years_e']


def return_dim_names():
    return dimensions


def add_units_and_scale(units, scale):
    uni = str(units)
    name_out = ""
    if scale != 1.:
        tmp = " *" + str("%e" % scale)
        if units != '':
            name_out += " (" + uni + tmp + ")"
        else:
            name_out += " (" + tmp + ")"
    else:
        if units != '':
            name_out += " (" + uni + ")"
    return name_out

def my_format(tab):
    interval = tab[1] - tab[0]
    tmp = str("%e" % abs(interval))
    exp = int(tmp.split('e')[1])
    if exp == -1:
        f_out = "{0:.1f}"
    elif exp == -2:
        f_out = "{0:.2f}"
    elif exp == -3:
        f_out = "{0:.3f}"
    elif exp < -3 or exp > 3:
        f_out = "{0:.1e}"
    else:
        f_out = "{0:d}"
    return f_out


def axes_to_time(axis, axis_name):
    axis_out = list(axis)
    labe_out = ''
    name_out = deepcopy(axis_name)
    if len(axis) == 12:
        axis_out = list(range(12))
        if 'month' in axis_name or 'time' in axis_name:
            locale.setlocale(category=locale.LC_ALL, locale="en_US")
            labe_out = calendar.month_abbr[1:13]
            name_out = 'months'
    return axis_out, labe_out, name_out


def compute_1d_function_on_nd(tab, myfct, myfctval='', myname=''):
    # swith to numpy
    dataset = NUMPYma__core__MaskedArray(tab)
    # masked values -> nan
    dataset = dataset.filled(fill_value=NUMPYnan)
    # Store information about the shape/size of the input data
    time_ax = dataset.shape[0]
    spac_ax = dataset.shape[1:]
    channels = NUMPYproduct(spac_ax)
    # Reshape to two dimensions (time, space) creating the design matrix
    dataset = dataset.reshape([time_ax, channels])
    # Find the indices of values that are not missing in one row. All the rows will have missing values in the same
    # places provided the array was centered. If it wasn't then it is possible that some missing values will be
    # missed and the singular value decomposition will produce not a number for everything.
    nonMissingIndex = NUMPYwhere(NUMPYisnan(dataset[0]) == False)[0]
    # Remove missing values from the design matrix.
    dataNoMissing = dataset[:, nonMissingIndex]
    if myfctval != '':
        new_dataset = myfct(dataNoMissing, myfctval, axis=0)
    else:
        new_dataset = myfct(dataNoMissing, axis=0)
    flatE = NUMPYones([channels], dtype=dataset.dtype) * NUMPYnan
    flatE = flatE.astype(dataset.dtype)
    flatE[nonMissingIndex] = new_dataset
    tab_out = flatE.reshape(spac_ax)
    tab_out = MV2array(MV2masked_where(NUMPYisnan(tab_out), tab_out))
    tab_out = CDMS2createVariable(tab_out, axes=tab.getAxisList()[1:], grid=tab.getGrid(), mask=tab[0].mask,
                                  attributes=tab.attributes, id=myname)
    return tab_out


def get_nbr_years(name):
    return int(name.split('years_')[0].split('_')[-1])


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


def find_var_category(name_cat, name_var, tab, tab_units):
    compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
    dict_type = dict()
    for name in name_cat:
        tmp1, tmp2, units = list(), list(), list()
        axes, axna, nbax = list(), list(), list()
        for var, val, uni in zip(name_var, tab, tab_units):
            if name in var:
                tmpa, tmpn = list(), list()
                for ii in range(len(val.shape)):
                    tmpa.append(list(val.getAxis(ii)[:]))
                    tmpn.append(val.getAxis(ii).id)
                axes.append(tmpa)
                axna.append(tmpn)
                nbax.append(len(val.shape))
                tmp1.append(val)
                tmp2.append(var)
                units.append(uni)
                del tmpa, tmpn
        # number of axes
        nbax = list(set(nbax))
        if len(nbax) == 1:
            nbax = nbax[0]
        else:
            list_strings = [
                'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
                'the number of axes is not the same for all variables: ' + str(axes)]
            MyError(list_strings)
        # axis values uniqueness
        axes_out = list()
        for ax in range(nbax):
            tmp = [compare(axes[0][ax], axes[ii+1][ax]) for ii in range(len(axes[1:]))]
            if all(ii is True for ii in tmp) is False:
                list_strings = [
                    'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
                    'multiple axe' + str(ax) + ' values for variable ' + str(name) + ': ' + str(axes[:][0][ax])]
                MyError(list_strings)
            else:
                axes_out.append(axes[0][ax])
        # axis names uniqueness
        axna_out = list()
        for ax in range(nbax):
            tmp = list(set([axna[ii][ax] for ii in range(len(axna))]))
            if len(tmp) != 1:
                list_strings = [
                    'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
                    'multiple axe' + str(ax) + ' names for variable ' + str(name) + ': ' + str(axna[:][ax])]
                MyError(list_strings)
            else:
                axna_out.append(axna[0][ax])
        # units uniqueness
        units = list(set(units))
        if len(units) > 1:
            list_strings = [
                'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
                'multiple units for variable ' + str(name) + ': ' + str(units)]
            MyError(list_strings)
        elif len(units) == 0:
            units = ''
        else:
            units = units[0]
        dict_type[name] = {'axes': axes_out, 'axes_name': axna_out, 'units': units, 'value': tmp1, 'variable': tmp2}
        del axes, axes_out, axna, axna_out, tmp1, tmp2, units
    return dict_type


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
    if abs(locmini) == locmaxi:
        mini_out = (interval / 2.) - base
        maxi_out = (interval / 2.) + base
    else:
        tmp_middle = locmini + (locmaxi - locmini) / 2.
        tmp2 = str("%e" % abs(tmp_middle))
        exp2 = int(tmp2.split('e')[1])
        mul2 = pow(10, exp2)
        tmp_middle = round(tmp_middle / mul2) * mul2
        mini_out = max([0, tmp_middle - interval / 2. - base])
        maxi_out = tmp_middle + interval / 2. + base
    if exp in [-1, 0, 1]:
        mini_out = mini_out * mult
        maxi_out = maxi_out * mult
        mult = 1.
    print "out = " + str([mini_out, maxi_out])
    return [mini_out, maxi_out], mult


def create_label(tab, nbr_sca=6, ratio=1.):
    """
    Computes range and create label for plots

    :param minmax: list of float
        minimum and maximum values of the array

    return: label: list of float
    """
    mini, maxi = deepcopy(min(tab)), deepcopy(max(tab))
    if mini<0 and maxi>0:
        maxi = max([abs(mini), abs(maxi)])
        mini = -deepcopy(maxi)
    interval = (maxi - mini) * ratio
    tmp = str("%e" % abs(interval))
    exp = int(tmp.split('e')[1])
    mult = pow(10,exp)
    locmini, locmaxi = float(mini) / mult, float(maxi) / mult
    interval = float(interval) / mult
    listbase = [0.1, 0.2, 0.4, 0.5, 1, 2, 4, 5]
    list1 = [round(base * nbr_sca, 1) if base < 1 else int(round(base * nbr_sca, 0)) for base in listbase]
    list2 = [abs(ii-interval) for ii in list1]
    interval = list1[list2.index(min(list2))]
    base = listbase[list1.index(interval)]
    if abs(locmini) == locmaxi:
        label = [ii * base - (interval/2.) for ii in range(nbr_sca + 1)]
    else:
        tmp_middle = locmini + (locmaxi - locmini) / 2.
        tmp2 = str("%e" % abs(tmp_middle))
        exp2 = int(tmp2.split('e')[1])
        mul2 = pow(10, exp2)
        tmp_middle = round(tmp_middle / mul2) * mul2
        label = [ii * base + (tmp_middle - interval / 2.) for ii in range(nbr_sca)]
        del exp2, mul2, tmp_middle, tmp2
    if exp in [-1,0,1]:
        base = float(base) * mult
        label = NUMPYarray(label) * mult
        mult = 1.
    tmp2 = str("%e" % base)
    exp2 = int(tmp2.split('e')[1])
    label = [round(ii, 2) if exp2 == -2 else (round(ii, 1) if exp2 == -1 else int(ii)) for ii in label]
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


# ---------------------------------------------------------------------------------------------------------------------#
# boxplot
def read_json(json_file, key_val):
    with open(json_file) as ff:
        data = json.load(ff)
    list_ref = sorted(data[key_val].keys())
    dict_out = dict()
    for ref in list_ref:
        val = data[key_val][ref]
        if len(val) == 1:
            dict_out[ref] = val[0]
        else:
            dict_out[ref] = val
        del val
    metric = data['metadata']['metrics'].keys()[0]
    if 'metric' in key_val:
        units = data['metadata']['metrics'][metric]['metric']['units']
    else:
        units = data['metadata']['metrics'][metric]['diagnostic']['units']
    return dict_out, units


def metric_boxplot(json_pattern, key_val, output_name, title, yname):
    # list all files
    list_files = sorted(list(GLOBiglob(json_pattern + ".json")))
    # nbr years used
    list_nbr = [get_nbr_years(file1) for file1 in list_files]
    # read json
    units = list()
    files_val = dict()
    for nbr, file1 in zip(list_nbr, list_files):
        met, uni = read_json(file1, key_val)
        files_val[nbr] = met
        units.append(uni)
        del met, uni
    units = list(set(units))
    if len(units) == 0:
        units = ''
    elif len(units) == 1:
        units = units[0]
    else:
        list_strings = [
            'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2]),
            'multiple units for variable ' + str(yname) + ': ' + str(units)]
        MyError(list_strings)
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
    marprops = dict(marker='o', markersize=2.0, markeredgecolor='k', markerfacecolor='k', markeredgewidth=0)
    meaprops = dict(marker='D', markersize=8.0, markeredgecolor='r', markerfacecolor='r', markeredgewidth=0)
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
    minmax, mult = create_minmax_plot(tab)
    add_to = add_units_and_scale(units, mult)
    yname_plot = yname + add_to
    if mult != 1.:
        for ref in list_ref:
            if isinstance(dict_out[ref], float):
                dict_out[ref] = dict_out[ref] / mult
            else:
                dict_out[ref] = NUMPYarray(dict_out[ref]) / mult
    if 'metric' in key_val:
        for ref in list_ref:
            fig1, ax1 = plt.subplots()
            plt.ylim(min(minmax), max(minmax))
            plt.title(title, fontsize=20)
            plt.xlabel('number of simulated years', fontsize=15)
            plt.ylabel(yname_plot, fontsize=15)
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
        plt.ylim(min(minmax), max(minmax))
        plt.title(title, fontsize=20, y=1.02)
        plt.xlabel('number of simulated years', fontsize=15)
        plt.ylabel(yname_plot, fontsize=15)
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


# ---------------------------------------------------------------------------------------------------------------------#
# curve or map plot
def read_nc_bydim(nc_file, metric, nbr_dim):
    # open file
    ff = CDMS2open(nc_file)
    # read 1d variable names
    var_in = [var for var in ff.listvariables() if var not in dimensions]
    var_in.sort(key=lambda v: v.lower())
    tab, tabu = list(), list()
    for var in var_in:
        # read variable in the given time, latitude, longitude window
        if nbr_dim == "1d" and\
                (metric == 'EnsoDiversity' or metric == 'NinaSstDiv' or metric == 'NinoSstDiv' or 'SstDur' in metric):
            tmp = ff(var)
            axe = CDMS2createAxis(MV2array(range(len(tmp)), dtype='f'), id='years')
            tmp.setAxis(0, axe)
            tab.append(tmp)
        else:
            tab.append(ff(var))
        list_att = ff.listattribute(var)
        for att in ['_FillValue', 'coordinates', 'history', 'missing_value', 'var_desc']:
            while att in list_att:
                list_att.remove(att)
        if 'units' in ff.listattribute(var):
            tabu.append(ff.getattribute(var, 'units'))
    return var_in, tab, tabu


def metric_curveplot(name_plot, title_plot, x_axis, x_label, x_name, y_name, y_range, dict_to_plot):
    model = dict_to_plot['model']
    mean = MV2average(model, axis=0)
    c05 = MV2array(SCIPYstats__scoreatpercentile(model, 5, axis=0))
    c25 = MV2array(SCIPYstats__scoreatpercentile(model, 25, axis=0))
    c75 = MV2array(SCIPYstats__scoreatpercentile(model, 75, axis=0))
    c95 = MV2array(SCIPYstats__scoreatpercentile(model, 95, axis=0))
    fig1, ax1 = plt.subplots()
    plt.ylim(min(y_range), max(y_range))
    plt.title(title_plot, fontsize=20, y=1.02)
    plt.xlabel(x_name, fontsize=15)
    plt.ylabel(y_name, fontsize=15)
    if x_label != '':
        plt.xticks(x_axis, x_label)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    ax1.plot(x_axis, list(mean), lw=4, label='model', color=dict_colors['model'])
    ax1.fill_between(x_axis, list(c05), list(c95), facecolor=dict_colors['model'], alpha=0.3)
    ax1.fill_between(x_axis, list(c25), list(c75), facecolor=dict_colors['model'], alpha=0.2)
    lines, names = [Line2D([0], [0], lw=2, color=dict_colors['model'])], ['model']
    datkeys = sorted(dict_to_plot.keys())
    for dat in datkeys:
        if dat != 'model':
            ax1.plot(x_axis, list(dict_to_plot[dat]), lw=4, label=dat, color=dict_colors[dat])
            lines.append(Line2D([0], [0], color=dict_colors[dat], lw=2))
            names.append(dat)
    plt.legend(lines, names)
    plt.savefig(name_plot)
    plt.close()
    return


def metric_mapplot(name_plot, title_plot, lats, lons, z1_range, z2_range, units, dict_to_plot, color_bar='Oranges',
                   z2_name=''):
    for dat in dict_to_plot.keys():
        if dat != 'model':
            mean = deepcopy(dict_to_plot[dat])
        else:
            tab = deepcopy(dict_to_plot[dat])
    for dat in dict_to_plot.keys():
        if dat == 'model':
            tab = deepcopy(dict_to_plot[dat])
            mean = MV2average(tab, axis=0)
            mean = CDMS2createVariable(mean, axes=tab.getAxisList()[1:], grid=tab.getGrid(), mask=tab[0].mask,
                                       id='mean')
            mean.setAxisList(tab.getAxisList()[1:])
            mean.setGrid(tab.getGrid())
            c05 = compute_1d_function_on_nd(tab, SCIPYstats__scoreatpercentile, myfctval=5, myname='c05')
            c95 = compute_1d_function_on_nd(tab, SCIPYstats__scoreatpercentile, myfctval= 95, myname='c95')
            spread = c95 - c05
            del c05, c95
        else:
            mean = deepcopy(dict_to_plot[dat])
            spread = MV2zeros(mean.shape)
            if color_bar in ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdBu', 'bwr', 'seismic']:
                spread.fill((z2_range[0]+z2_range[-1])/2.)
            elif min(z2_range) == 0:
                spread.fill(-1)
        spread = CDMS2createVariable(MV2array(spread), axes=mean.getAxisList(), grid=mean.getGrid(), mask=mean.mask,
                                     id=z2_name)
        xx, yy = NUMPYmeshgrid(lons, lats)
        map = Basemap(projection='cyl', llcrnrlat=lats[0]-5, urcrnrlat=lats[-1]+5, llcrnrlon=lons[0],
                      urcrnrlon=lons[-1])
        # draw coastlines
        map.drawcoastlines()
        # draw countries
        # map.drawcountries()
        # fill continents
        map.fillcontinents(color='gainsboro')
        # draw parallels
        parallels = NUMPYarange(-80., 80, 10.)
        map.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=15, dashes=[6,900], linewidth=1)
        # draw meridians
        meridians = NUMPYarange(0., 360., 30.)
        map.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=15, dashes=[10,900], linewidth=1)
        # draw mesh
        cs = map.pcolormesh(xx, yy, spread, vmin=z2_range[0], vmax=z2_range[-1], cmap=color_bar)#, cmap='bwr')
        if color_bar not in ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdBu', 'bwr', 'seismic']:
            cs.cmap.set_under('white')
        # add colorbar
        cbar = map.colorbar(cs, location='bottom', ticks=z2_range[::2], pad=0.35, size=0.15, extend='both')
        cbar.set_label(z2_name + " (" + units + ")", fontsize=15)
        cbar.ax.tick_params(labelsize=15)
        # contour data over the map
        cs = map.contour(xx, yy, mean, z1_range[::2], linewidths=1.5, colors=['black'])
        my_f = my_format(z1_range)
        # Recast levels to new class
        cs.levels = [my_f.format(ii) for ii in cs.levels]
        plt.clabel(cs, cs.levels, inline=True, fontsize=10, use_clabeltext=True)
        plt.title(title_plot + " " + dat, fontsize=20, y=1.02)
        plt.savefig(name_plot + "_" + dat)
        plt.close()
    return


def myplot(nc_pattern, output_name, metric, title, plot_type, z_name=''):
    # list all files
    list_files = sorted(list(GLOBiglob(nc_pattern + "_" + str(plot_type) + ".nc")))
    # nbr years used
    list_nbr = [get_nbr_years(file1) for file1 in list_files]
    # loop on files
    dict_range = dict()
    for nbr, file1 in zip(list_nbr, list_files):
        print str().ljust(5) + str(nbr).zfill(4) + 'years file: ' + file1
        var_names, tab, tabunits = read_nc_bydim(file1, metric, plot_type)
        name_type = list(set([var.split('__')[0] for var in var_names]))
        dict_type = find_var_category(name_type, var_names, tab, tabunits)
        del tab, tabunits, var_names
        for var_typ in name_type:
            axis, axis_label, axis_name = list(), list(), list()
            for ii in range(len(dict_type[var_typ]['axes'])):
                tmp1, tmp2, tmp3 = axes_to_time(dict_type[var_typ]['axes'][ii], dict_type[var_typ]['axes_name'][ii])
                axis.append(tmp1)
                axis_label.append(tmp2)
                axis_name.append(tmp3)
                del tmp1, tmp2, tmp3
            units = dict_type[var_typ]['units']
            dict_data = dict()
            for ii, var in enumerate(dict_type[var_typ]['variable']):
                if isobs(var) is True:
                    dict_data[var.split('__')[1]] = dict_type[var_typ]['value'][ii]
                else:
                    try:
                        dict_data['model']
                    except:
                        dict_data['model'] = [dict_type[var_typ]['value'][ii]]
                    else:
                        dict_data['model'] += [dict_type[var_typ]['value'][ii]]
            dict_data['model'] = MV2array(dict_data['model'])
            try:
                dict_range[var_typ]
            except:
                tab = list()
                for dat in dict_data.keys():
                    tab += [MV2minimum(dict_data[dat]), MV2maximum(dict_data[dat])]
                if plot_type == '1d':
                    minmax, mult = create_minmax_plot(tab)
                else:
                    minmax, mult = create_label(tab, nbr_sca=12, ratio=0.75)
                dict_range[var_typ] = {'minmax': minmax, 'mul': mult}
                del minmax, mult, tab
            add_to = add_units_and_scale(units, dict_range[var_typ]['mul'])
            if dict_range[var_typ]['mul'] != 1.:
                for dat in dict_data.keys():
                    dict_data[dat] = dict_data[dat] / dict_range[var_typ]['mul']
            if plot_type == "1d":
                name_png = output_name + '_' + var_typ + '_' + str(nbr).zfill(4)
                title_plot = title + " (" + str(nbr).zfill(4) + " years)"
                metric_curveplot(name_png, title_plot, axis[0], axis_label[0], axis_name[0], var_typ + add_to,
                                 dict_range[var_typ]['minmax'], dict_data)
            else:
                name_png = output_name + '_' + var_typ + '_' + str(nbr).zfill(4)
                title_plot = title + " (" + str(nbr).zfill(4) + " years)" + add_to
                minmax2 = NUMPYarray(dict_range[var_typ]['minmax']) - dict_range[var_typ]['minmax'][0]
                minmax2 = minmax2 * 0.75
                metric_mapplot(name_png, title_plot, axis[0], axis[1], dict_range[var_typ]['minmax'],
                               minmax2, units, dict_data, z2_name=z_name)
            del name_png, title_plot
            del add_to, axis, axis_label, axis_name, dict_data, units
        del dict_type, name_type
    return
# ---------------------------------------------------------------------------------------------------------------------#

