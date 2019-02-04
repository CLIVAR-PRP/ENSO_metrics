# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Library for the main driver to plot outputs of the ENSO_metrics package
#---------------------------------------------------#


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
        print bcolors.FAIL + '%%%%%     -----     %%%%%'
        print 'unknown name (obs or model): ' + str(name)
        print 'known observations: ' + str(list_observations)
        print 'known models: ' + str(list_models)
        print '%%%%%     -----     %%%%%' + bcolors.ENDC
        SYSexit('')
    return name_out


def replace_in_obsname(name):
    new_name = deepcopy(name).replace('sst_map__', '').replace('pr__', '').replace('pr_lat__', '').replace('pr_lon__',
                                                                                                           '').replace(
        'pr_map__', '').replace('sst__', '').replace('sst_lat__', '').replace('sst_lon__', '').replace('sst_map__',
                                                                                                       '').replace(
        'sst_ts__', '').replace('sstComp_lon__', '').replace('sstComp_map__', '').replace('sstStd_MAM_map__',
                                                                                          '').replace(
        'sstStd_NDJ_map__', '').replace('sstSke_lon__', '').replace('sstSke_map__', '').replace('sstSke_monthly__',
                                                                                                '').replace(
        'sstStd_map__', '').replace('sstStd_monthly__', '').replace('taux__', '').replace('taux_lat__', '').replace(
        'taux_lon__', '').replace('taux_map__', '').replace('Nina_duration__', '').replace('Nina_lon_pos_minSSTA__',
                                                                                           '').replace(
        'Nino_duration__', '').replace('Nino_lon_pos_maxSSTA__', '').replace('pdf__', '')
    return new_name


def names_to_colors(name, dict_colors):
    if 'GPCPv2.3' in name:
        col = dict_colors['GPCPv2.3']
    elif 'HadISST' in name:
        col = dict_colors['HadISST']
    elif 'Tropflux' in name:
        col = dict_colors['Tropflux']
    elif 'IPSL-CM5A-LR' in name:
        col = dict_colors['IPSL-CM5A-LR']
    elif 'IPSL-CM5A-MR' in name:
        col = dict_colors['IPSL-CM5A-MR']
    elif 'IPSL-CM5B-LR' in name:
        col = dict_colors['IPSL-CM5B-LR']
    else:
        col = dict_colors['IPSL-CM6A-LR']
    return col


def both_nina_and_nino(list_var):
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
