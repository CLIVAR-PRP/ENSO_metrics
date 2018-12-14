from cdms2 import open as CDMS2open
from os.path import join as join_path
from os import environ
from sys import exit

#from EnsoMetrics.EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
#from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import ComputeCollection

xmldir = environ['XMLDIR']

def find_xml(name, frequency, variable, project='', experiment='', ensemble='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_name, file_area, file_land = find_xml_obs(name, frequency, variable)
    else:
        file_name, file_area, file_land = find_xml_cmip(name, project, experiment, ensemble, frequency, realm, variable)
    return file_name, file_area, file_land

def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                          '_glob_' + str(frequency) + '_' + str(realm) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        if realm == 'O':
            new_realm = 'A'
        elif realm == 'A':
            new_realm = 'O'
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere)
        file_name = join_path(xmldir, str(model) + '_' + str(project) + '_' + str(experiment) + '_' + str(ensemble) +
                              '_glob_' + str(frequency) + '_' + str(new_realm) + '.xml')
        xml = CDMS2open(file_name)
        listvar2 = sorted(xml.listvariables())
        if variable not in listvar2:
            print '\033[95m' + str().ljust(5) + "CMIP var " + str(variable) + " cannot be found (realm A and O)"\
                  + '\033[0m'
            print '\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m'
            print '\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m'
            print '\033[95m' + str().ljust(10) + "AND" + '\033[0m'
            print '\033[95m' + str().ljust(10) + "variables = " + str(listvar2) + '\033[0m'
            exit(1)
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=new_realm)
    else:
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=realm)
    return file_name, file_area, file_land

def find_xml_fx(name, project='', experiment='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_area = join_path(xmldir, 'obs_' + str(name) + '_glob_fx_O_areacell.xml')
        file_land = join_path(xmldir, 'obs_' + str(name) + '_glob_fx_O_landmask.xml')
    else:
        file_area = join_path(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_'
                                   + str(realm) + '_areacell.xml')
        file_land = join_path(xmldir, str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_'
                                   + str(realm) + '_landmask.xml')
    try: xml = CDMS2open(file_area)
    except: file_area = None
    try: xml = CDMS2open(file_land)
    except: file_land = None
    return file_area, file_land

def find_xml_obs(obs, frequency, variable):
    file_name = join_path(xmldir, 'obs_' + str(obs) + '_glob_' + str(frequency) + '_O.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print '\033[95m' + str().ljust(5) + "obs var " + str(variable) + " cannot be found" + '\033[0m'
        print '\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m'
        print '\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m'
        exit(1)
    file_area, file_land = find_xml_fx(obs)
    return file_name, file_area, file_land


# metric collection
mc_name = 'ENSO_perf'#'ENSO_tel'#'MC1'#
dict_mc = defCollection(mc_name)
list_metric = sorted(dict_mc['metrics_list'].keys())

# parameters
project = 'CMIP5'
experiment = 'historical'
ensemble = 'r1i1p1'
frequency = 'mon'
realm = 'A'

# list of variables
list_variables = list()
for metric in list_metric:
    listvar = dict_mc['metrics_list'][metric]['variables']
    for var in listvar:
        if var not in list_variables:
            list_variables.append(var)
list_variables = sorted(list_variables)
print '\033[95m' + str(list_variables) + '\033[0m'

# list of observations
list_obs = list()
for metric in list_metric:
    dict_var_obs = dict_mc['metrics_list'][metric]['obs_name']
    for var in dict_var_obs.keys():
        for obs in dict_var_obs[var]:
            if obs not in list_obs:
                list_obs.append(obs)
list_obs = sorted(list_obs)
if mc_name == 'MC1':
    list_obs = ['Tropflux']
elif mc_name == 'ENSO_perf':
    list_obs = ['Tropflux']#['HadISST']#['HadISST','Tropflux','GPCPv2.3']#['Tropflux','GPCPv2.3']
elif mc_name == 'ENSO_tel':
    list_obs = ['HadISST','GPCPv2.3']
print '\033[95m' + str(list_obs) + '\033[0m'


#
# finding file and variable name in file for each observations dataset
#
dict_obs = dict()
for obs in list_obs:
    # @jiwoo: be sure to add your datasets to EnsoCollectionsLib.ReferenceObservations if needed
    dict_var = ReferenceObservations(obs)['variable_name_in_file']
    dict_obs[obs] = dict()
    for var in list_variables:
        #
        # finding variable name in file
        #
        # @jiwoo: correct / adapt the 'varname' in
        # EnsoCollectionsLib.ReferenceObservations(obs)['variable_name_in_file'][var] if it is not correct or if you
        # changed a name in the xml
        # I usually alias the variable names from observations and models in the xml in order to have the same name
        # for sst (or any other variable) in every xml. This way I don not need to go through this function to know the
        # variable name in file
        try: var_in_file = dict_var[var]['var_name']
        except:
            print '\033[95m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m'
        else:
            if isinstance(var_in_file, list):
                var0 = var_in_file[0]
            else:
                var0 = var_in_file
            #
            # finding file for 'obs', 'var'
            #
            # @jiwoo: pretty easy as I have all variables in one file
            file_name, file_areacell, file_landmask = find_xml(obs, frequency, var0)
            try:
                areacell_in_file = dict_var['areacell']['var_name']
            except:
                areacell_in_file = None
            try:
                landmask_in_file = dict_var['landmask']['var_name']
            except:
                landmask_in_file = None
            # if var_in_file is a list (like for thf) all variables should be read from the same realm
            if isinstance(var_in_file, list):
                list_files = [file_name for var1 in var_in_file]
                list_areacell = [file_areacell for var1 in var_in_file]
                list_name_area = [areacell_in_file for var1 in var_in_file]
                list_landmask = [file_landmask for var1 in var_in_file]
                list_name_land = [landmask_in_file for var1 in var_in_file]
            else:
                list_files = file_name
                list_areacell = file_areacell
                list_name_area = areacell_in_file
                list_landmask = file_landmask
                list_name_land = landmask_in_file
            dict_obs[obs][var] = {'path + filename': list_files, 'varname': var_in_file,
                                  'path + filename_area': list_areacell, 'areaname': list_name_area,
                                  'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}

# models
list_models = ['IPSL-CM5B-LR']#['CNRM-CM5']#['IPSL-CM5B-LR']#['CNRM-CM5','IPSL-CM5B-LR']#
#
# finding file and variable name in file for each observations dataset
#
dict_metric, dict_dive = dict(), dict()
dict_var = CmipVariables()['variable_name_in_file']
for mod in list_models:
    dict_mod = {mod: {}}
    # ------------------------------------------------
    # @jiwoo: between these dash the program is a bit ad hoc...
    # it works well for me because I am looking for sst and taux on the ocean grid, and fluxes [lhf, lwr, swr, shf, thf]
    # on the atmosphere grid
    # if you want to use atmosphere only, do not use this or create your own way to find the equivalent between the
    # variable name in the program and the variable name in the file
    for var in list_variables:
        #
        # finding variable name in file
        #
        var_in_file = dict_var[var]['var_name']
        if isinstance(var_in_file, list):
            var0 = var_in_file[0]
        else:
            var0 = var_in_file
        #
        # finding file for 'mod', 'var'
        #
        # @jiwoo: first try in the realm 'O' (for ocean)
        file_name, file_areacell, file_landmask = find_xml(mod, frequency, var0, project=project, experiment=experiment,
                                                           ensemble=ensemble, realm=realm)
        try:
            areacell_in_file = dict_var['areacell']['var_name']
        except:
            areacell_in_file = None
        try:
            landmask_in_file = dict_var['landmask']['var_name']
        except:
            landmask_in_file = None
        if isinstance(var_in_file, list):
            list_files = [file_name for var1 in var_in_file]
            list_areacell = [file_areacell for var1 in var_in_file]
            list_name_area = [areacell_in_file for var1 in var_in_file]
            list_landmask = [file_landmask for var1 in var_in_file]
            list_name_land = [landmask_in_file for var1 in var_in_file]
        else:
            list_files = file_name
            list_areacell = file_areacell
            list_name_area = areacell_in_file
            list_landmask = file_landmask
            list_name_land = landmask_in_file
        dict_mod[mod][var] = {'path + filename': list_files, 'varname': var_in_file,
                              'path + filename_area': list_areacell, 'areaname': list_name_area,
                              'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}
    # dictionary needed by nsoMetrics.ComputeMetricsLib.ComputeCollection
    # @jiwoo the ComputeCollection function it still on development and it does not read the observations requirement
    # defined in the metric collection, i.e., defCollection(mc_name)['metrics_list']['<metric name>']['obs_name']
    # so the function does not take a specific obs to compute the metric so for every obs in 'dict_obs' we must include
    # every variables needed by the metric collection [lhf, lwr, swr, shf, sst, taux, thf] even if its coming from
    # another dataset
    dictDatasets = {'model': dict_mod, 'observations': dict_obs}
    # regridding dictionary (only if you want to specify the regridding)
    dict_regrid = {}
    # dict_regrid = {
    #     'regridding': {
    #         'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
    #         'newgrid_name': 'generic 1x1deg'},
    # }
    # Computes the metric collection
    #dict_metric[mod] = ComputeCollection(mc_name, dictDatasets, user_regridding=dict_regrid, debug=False)
    #dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, user_regridding=dict_regrid,
    #                                                     debug=False, dive_down=True)
    netcdf_path = '/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Test'
    netcdf_name = '20181210_YANN_PLANTON_' + mc_name + '_' + mod
    dict_metric[mod], dict_dive[mod] = ComputeCollection(mc_name, dictDatasets, netcdf=True, netcdf_path=netcdf_path,
                                                         netcdf_name=netcdf_name, debug=True)
    stop
    # Prints the metrics values
    for ii in range(3): print ''
    print '\033[95m' + str().ljust(5) + str(mod) + '\033[0m'
    list_metric = dict_metric[mod]['value'].keys()
    for metric in list_metric:
        print '\033[95m' + str().ljust(10) + str(metric) + '\033[0m'
        metric_dict = dict_metric[mod]['value'][metric]['metric']
        for ref in metric_dict.keys():
            print '\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + str(metric_dict[ref]['value'])\
                  + ', error = ' + str(metric_dict[ref]['value_error']) + '\033[0m'
            if 'value2' in metric_dict[ref].keys():
                print '\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' +\
                      str(metric_dict[ref]['value2']) + ', error = ' + str(metric_dict[ref]['value_error2']) + '\033[0m'
            if 'value3' in metric_dict[ref].keys():
                print '\033[95m' + str().ljust(15) + 'metric: ' + str(ref) + ' value = ' + \
                      str(metric_dict[ref]['value3']) + ', error = ' + str(metric_dict[ref]['value_error3']) + '\033[0m'
# Plot
#stop
#if ' ':
for mod in list_models:
    from numpy import arange as NUMPYarange
    from cdms2 import createAxis as CDMS2createAxis
    from MV2 import array as MV2array
    from MV2 import masked_where as MV2masked_where
    from MV2 import maximum as MV2maximum
    from MV2 import minimum as MV2minimum
    import plot_frame as PFRAME
    import plot_functions as PF
    path_plot = '/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Plots'
#if ' ':
    if ' ':
        for metric in list_metric:
            print '\033[95m' + str().ljust(10) + str(metric) + '\033[0m'
            metric_dict = dict_metric[mod]['value'][metric]['metric']
            # metric
            dict_m1, dict_m2, dict_m3 = dict(), dict(), dict()
            for ref in metric_dict.keys():
                dict_m1[ref] = metric_dict[ref]['value']
                if 'value2' in metric_dict[ref].keys():
                    dict_m2[ref] = metric_dict[ref]['value2']
                if 'value3' in metric_dict[ref].keys():
                    dict_m3[ref] = metric_dict[ref]['value3']
            # dive down
            dive_model = dict_dive[mod]['value'][metric]['model']
            if metric in ['EnsoPrMap', 'EnsoSstMap']:
                tmp_dive, tmp_axis = dict(), dict()
                for ref in dict_dive[mod]['value'][metric].keys():
                    if ref != 'model':
                        tmp_dive['ref_' + ref] = dict_dive[mod]['value'][metric][ref]
                        tmp1 = dict_dive[mod]['metadata']['metrics'][metric][ref]['axisLat']
                        axis1 = CDMS2createAxis(MV2array(tmp1), id='latitude')
                        tmp2 = dict_dive[mod]['metadata']['metrics'][metric][ref]['axisLon']
                        axis2 = CDMS2createAxis(MV2array(tmp2), id='longitude')
                        tmp_axis['ref_' + ref] = [axis1, axis2]
                        del axis1, axis2, tmp1, tmp2
                # plot
                x_axis, inc = [0, 360], 60
                x_dict = dict((ii, str(ii) + 'E') if ii <= 180 else (ii, str(abs(ii-360)) + 'W') for ii in
                              range(x_axis[0], x_axis[1] + inc, inc))
                y_axis, inc = [-60, 60], 20
                y_dict = dict((ii, str(abs(ii)) + 'S') if ii < 0 else ((ii, str(ii) + 'N') if ii>0 else (ii, 'Eq')) for
                              ii in range(y_axis[0], y_axis[1] + inc, inc))
                dom = (y_axis[0], y_axis[1], x_axis[0], x_axis[1])
                if metric in ['EnsoPrMap']:
                    label_col = MV2array(range(-3, 3 + 1, 1))
                elif metric in ['EnsoSstMap']:
                    label_col = MV2array([round(ii, 1) for ii in NUMPYarange(-1.2, 1.2 + 0.4, 0.4)])
                for ref in dict_m1.keys():
                    tab1 = MV2array(dive_model)
                    tab1.setAxisList(tmp_axis[ref])
                    m1 = 'Metric 1: ' + str("%.2f" % round(dict_m1[ref], 2))
                    m2 = 'Metric 2: ' + str("%.2f" % round(dict_m2[ref], 2))
                    m3 = 'Metric 3: ' + str("%.2f" % round(dict_m3[ref], 2))
                    tab2 = MV2array(tmp_dive[ref])
                    tab2.setAxisList(tmp_axis[ref])
                    print str().ljust(10) + 'range = ' + str("%.2f" % round(min(MV2minimum(tab1),MV2minimum(tab2)), 2))\
                          + ' ' + str("%.2f" % round(max(MV2maximum(tab1),MV2maximum(tab2)), 2))
                    name = metric + ' in Historical (' + mod + ')'
                    name_png = path_plot + '/' + metric + '_' + mod
                    PFRAME.plot_my_map(tab1, label_col, dom, white_zero=0, x_dico=x_dict, y_dico=y_dict, name=name,
                                       path_plus_name_png=name_png, bg=1)
                    name = metric + ' in Historical (' + ref + ')'
                    name_png = path_plot + '/' + metric + '_' + ref
                    PFRAME.plot_my_map(tab2, label_col, dom, white_zero=0, x_dico=x_dict, y_dico=y_dict, name=name,
                                       path_plus_name_png=name_png, bg=1)
                    del m1, m2, m3, name, name_png, tab1, tab2
            elif metric in ['EnsoPrJjaTel', 'EnsoPrNdjTel']:
                tmp_dive, tmp_axis = dict(), dict()
                for ref in dict_dive[mod]['value'][metric].keys():
                    if ref != 'model':
                        tmp_dive['ref_' + ref] = dict_dive[mod]['value'][metric][ref]
                        tmp_axis['ref_' + ref] = dict_dive[mod]['metadata']['metrics'][metric][ref]['axis']
                # plot
                y_axis, inc = [-2.0, 6.2], 0.2
                y_dict = dict((round(elt, 1), "{0:.1f}".format(round(elt, 1)))
                              if (round(elt, 1) * 10) % round(5 * round(inc, 1) * 10, 1) == 0
                              else (round(elt, 1), '') for elt in NUMPYarange(y_axis[0], y_axis[1] + inc, inc))
                for ref in dict_m1.keys():
                    axis = CDMS2createAxis(MV2array(range(len(tmp_axis[ref])), dtype='int32'), id='regions')
                    x_dict = dict((elt, tmp_axis[ref][elt]) for elt in range(len(tmp_axis[ref])))
                    x_axis = [-1.0, len(tmp_axis[ref])]
                    tab1 = MV2array(dive_model)
                    tab1.setAxisList([axis])
                    m1 = 'Metric 1: ' + str("%.2f" % round(dict_m1[ref], 2))
                    m2 = 'Metric 2: ' + str("%.1f" % round(dict_m2[ref]*100, 1))
                    tab2 = MV2array(tmp_dive[ref])
                    tab2.setAxisList([axis])
                    print str().ljust(10) + 'range = ' + str("%.2f" % round(min(MV2minimum(tab1),MV2minimum(tab2)), 2))\
                          + ' ' + str("%.2f" % round(max(MV2maximum(tab1),MV2maximum(tab2)), 2))
                    list_curve = [tab1, tab2]
                    list_col = ['black', 'red']
                    # strings to write
                    l_w = [m1, m2]
                    l_w_xy = [[97, 100 - (ii + 1) * 6] for ii in range(len(l_w))]
                    l_w_si = [30 for ii in range(len(l_w))]
                    l_w_ha = ['right' for ii in range(len(l_w))]
                    # lines to plot
                    lines_y1y2 = [[round(ii, 1), round(ii, 1)] for ii in y_dict.keys() if y_dict[ii] != '' and
                                  round(ii, 1) != 0 and round(ii, 1) not in y_axis]
                    lines_x1x2 = [x_axis for ii in range(len(lines_y1y2))]
                    lines_colo = ['grey' for ii in range(len(lines_y1y2))]
                    name = metric + ' metric in Historical (' + mod + ')'
                    yname = 'El Nino (PR) minus La Nina (PR)'
                    name_png = path_plot + '/' + metric + '_' + mod + '_ ' + ref
                    PFRAME.curves_plot(list_curve, list_col=list_col, x_axis=x_axis, x_dico=x_dict, y_axis=y_axis,
                                       y_dico=y_dict, name_in_xlabel=True, name=name, xname='', yname=yname,
                                       list_writings=l_w, list_writings_pos_xy=l_w_xy, list_writings_size=l_w_si,
                                       list_writings_halign=l_w_ha, plot_lines=True, lines_x1x2=lines_x1x2,
                                       lines_y1y2=lines_y1y2, lines_color=lines_colo, path_plus_name_png=name_png,
                                       draw_white_background=True, save_ps=False, bg=1)
                    del l_w, l_w_ha, l_w_si, l_w_xy, lines_colo, lines_x1x2, lines_y1y2, list_curve, m1,\
                        m2, name, name_png, yname
            elif metric in ['BiasSstLonRmse', 'BiasSstSkLonRmse', 'SeasonalSstLonRmse', 'NinaSstTsRmse',
                            'NinoSstTsRmse']:
                tmp_dive, tmp_axis = dict(), dict()
                for ref in dict_dive[mod]['value'][metric].keys():
                    if ref != 'model':
                        tmp_dive['ref_' + ref] = dict_dive[mod]['value'][metric][ref]
                        tmp_axis['ref_' + ref] = dict_dive[mod]['metadata']['metrics'][metric][ref]['axis']
                # plot
                for ref in dict_m1.keys():
                    axis = CDMS2createAxis(MV2array(tmp_axis[ref], dtype='float32'), id='axis')
                    tab1 = MV2array(dive_model)
                    tab1.setAxisList([axis])
                    tab1 = MV2masked_where(tab1>=1e20, tab1)
                    if metric in ['BiasSstLonRmse', 'BiasSstSkLonRmse', 'NinoSstLonRmse', 'SeasonalSstLonRmse']:
                        inc = 30
                        if min(axis[:])<0:
                            x_axis = [-250, -70]
                            tmp = [x_axis[0]+10, x_axis[1]-10]
                            x_dict = dict((ii, str(ii + 360) + 'E') if ii < -180 else (ii, str(abs(ii)) + 'W') for ii in
                                    range(tmp[0], tmp[1] + inc, inc))
                        else:
                            x_axis = [110, 290]
                            tmp = [x_axis[0] + 10, x_axis[1] - 10]
                            x_dict = dict((ii, str(ii) + 'E') if ii < 180 else (ii, str(abs(ii - 360)) + 'W') for ii in
                                          range(tmp[0], tmp[1] + inc, inc))
                    elif metric in ['NinoSstTsRmse']:
                        x_axis, inc = [-1, len(axis)], 1
                        tmp = ['M', 'J', 'S', 'D']
                        x_dict = dict((ii, tmp[(((ii + 1) / 3) % 4) - 1]) if (ii + 1) % 3 == 0 else (ii, '') for ii in
                                      range(x_axis[0], x_axis[1] + inc, inc))
                    m1 = 'Metric: ' + str("%.2f" % round(dict_m1[ref], 2))
                    tab2 = MV2array(tmp_dive[ref])
                    tab2.setAxisList([axis])
                    tab2 = MV2masked_where(tab2 >= 1e20, tab2)
                    print str().ljust(10) + 'range = ' + str("%.2f" % round(min(MV2minimum(tab1),MV2minimum(tab2)), 2))\
                          + ' ' + str("%.2f" % round(max(MV2maximum(tab1),MV2maximum(tab2)), 2))
                    y_axis, y_dict = PF.create_dico([min(MV2minimum(tab1),MV2minimum(tab2)),
                                                    max(MV2maximum(tab1),MV2maximum(tab2))])
                    list_curve = [tab1, tab2]
                    list_col = ['black', 'red']
                    # strings to write
                    l_w = [m1]
                    l_w_xy = [[97, 100 - (ii + 1) * 6] for ii in range(len(l_w))]
                    l_w_si = [30 for ii in range(len(l_w))]
                    l_w_ha = ['right' for ii in range(len(l_w))]
                    # lines to plot
                    lines_y1y2 = [[round(ii, 1), round(ii, 1)] for ii in y_dict.keys() if y_dict[ii] != '' and
                                  round(ii, 1) != 0 and round(ii, 1) not in y_axis]
                    lines_x1x2 = [x_axis for ii in range(len(lines_y1y2))]
                    if metric in ['BiasSstLonRmse', 'BiasSstSkLonRmse', 'NinoSstLonRmse', 'SeasonalSstLonRmse']:
                        xname = 'longitude'
                        lines_x1x2 = lines_x1x2 + [[ii, ii] for ii in x_dict.keys() if x_dict[ii] != '' and ii != 0
                                                   and ii not in x_axis]
                        lines_y1y2 = lines_y1y2 + [y_axis for ii in x_dict.keys() if x_dict[ii] != '' and ii != 0
                                                   and ii not in x_axis]
                    elif metric in ['NinoSstTsRmse']:
                        xname = 'time'
                        lines_x1x2 = lines_x1x2 + [[ii, ii] for ii in x_dict.keys() if (ii + 1) % 12 == 0 and ii != 0
                                                   and ii not in x_axis]
                        lines_y1y2 = lines_y1y2 + [y_axis for ii in x_dict.keys() if (ii + 1) % 12 and ii != 0
                                                   and ii not in x_axis]
                    lines_colo = ['grey' for ii in range(len(lines_y1y2))]
                    name = metric + ' metric (' + mod + ')'
                    print metric, mod, ref
                    name_png = path_plot + '/' + metric + '_' + mod + '_ ' + ref
                    if metric in ['NinoSstLonRmse', 'NinoSstTsRmse', 'SeasonalSstLonRmse']:
                        yname = 'SSTA (degC)'
                    elif metric in ['BiasSstLonRmse', 'BiasSstSkLonRmse']:
                        yname = 'SST (degC)'
                    PFRAME.curves_plot(list_curve, list_col=list_col, x_axis=x_axis, x_dico=x_dict, y_axis=y_axis,
                                       y_dico=y_dict, name_in_xlabel=False, name=name, xname=xname, yname=yname,
                                       list_writings=l_w, list_writings_pos_xy=l_w_xy, list_writings_size=l_w_si,
                                       list_writings_halign=l_w_ha, plot_lines=True, lines_x1x2=lines_x1x2,
                                       lines_y1y2=lines_y1y2, lines_color=lines_colo, path_plus_name_png=name_png,
                                       draw_white_background=True, save_ps=False, bg=1)
                    del l_w, l_w_ha, l_w_si, l_w_xy, lines_colo, lines_x1x2, lines_y1y2, list_curve, m1,\
                        name, name_png, xname, yname


