# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Library for the main driver to compute the ENSO_metrics package
#---------------------------------------------------#


#---------------------------------------------------#
# import python packages
# usual python package
from copy import deepcopy
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
from math import floor as MATHfloor
from os import remove as OSremove
from os.path import join as OSpath__join
import sys
# CDAT package
from cdms2 import open as CDMS2open
# ENSO_metrics package
from EnsoCollectionsLib import CmipVariables, defCollection, ReferenceObservations
from EnsoComputeMetricsLib import InternCompute
# My (YYP) package
# set new path where to find programs
sys.path.insert(0, "/home/yplanton/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import find_path_and_files
from getfiles_sh_to_py import get_ensembles
from getfiles_sh_to_py import get_time_size
#---------------------------------------------------#


#---------------------------------------------------#
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
#---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# define path and dimensions
xmldir = '/home/yplanton/New_XMLDIR'
# CMIP variable names
dict_CMIPvar = CmipVariables()['variable_name_in_file']
dimensions = ['bounds_latitude', 'bounds_latitude_a', 'bounds_latitude_b', 'bounds_longitude', 'bounds_longitude_a',
              'bounds_longitude_b', 'bounds_months', 'bounds_time', 'bounds_years', 'bounds_years_a', 'bounds_years_b',
              'bounds_years_c', 'bounds_years_d',' bounds_years_e', 'latitude', 'latitude_a', 'latitude_b', 'longitude',
              'longitude_a', 'longitude_b', 'months', 'years', 'years_a', 'years_b', 'years_c', 'years_d', 'years_e']


def return_dim_names():
    return dimensions
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# find files
def file_model(experiment, ensemble, frequency, model, project, realm, list_var):
    dict_mod = {model: {}}
    for var in list_var:
        # find variable name in file
        var_in_file = dict_CMIPvar[var]['var_name']
        if isinstance(var_in_file, list):
            var0 = var_in_file[0]
        else:
            var0 = var_in_file
        # find file for 'mod', 'var'
        file_name, file_areacell, file_landmask =\
            find_xml(model, frequency, var0, project=project, experiment=experiment, ensemble=ensemble, realm=realm)
        try:
            areacell_in_file = dict_CMIPvar['areacell']['var_name']
        except:
            areacell_in_file = None
        try:
            landmask_in_file = dict_CMIPvar['landmask']['var_name']
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
        dict_mod[model][var] = {'path + filename': list_files, 'varname': var_in_file,
                                'path + filename_area': list_areacell, 'areaname': list_name_area,
                                'path + filename_landmask': list_landmask, 'landmaskname': list_name_land}
    return dict_mod


def file_obs(list_obs, list_var, frequency):
    # find file and variable name in file for each observations dataset
    dict_obs = dict()
    for obs in list_obs:
        # get variable name in file, defined by the ENSO_metric package
        dict_var = ReferenceObservations(obs)['variable_name_in_file']
        dict_obs[obs] = dict()
        for var in list_var:
            # find variable name in file
            try:
                var_in_file = dict_var[var]['var_name']
            except:
                print '\033[93m' + str(var) + " is not available for " + str(obs) + " or unscripted" + '\033[0m'
            else:
                if isinstance(var_in_file, list):
                    var0 = var_in_file[0]
                else:
                    var0 = var_in_file
                # find file for 'obs', 'var'
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
    return dict_obs


def find_fx(model, project='', experiment='', ensemble='', realm=''):
    if project in ['CMIP5', 'CMIP6']:
        if project in ['CMIP5']:
            my_ens = 'r0i0p0'
        else:
            my_ens = deepcopy(ensemble)
        if realm == 'A':
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='areacella')
            fland1, fland2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='sftlf')
            file_land = OSpath__join(fland1, fland2[0])
        else:
            farea1, farea2 = find_path_and_files(ens=my_ens, exp=experiment, fre='fx', mod=model, pro=project,
                                                 rea=realm, var='areacello')
            file_land = None
        file_area = OSpath__join(farea1, farea2[0])
    else:
        file_area, file_land = find_xml_fx(model, project=project, experiment=experiment, realm=realm)
    try:
        CDMS2open(file_area)
    except:
        file_area = None
    try:
        CDMS2open(file_land)
    except:
        file_land = None
    return file_area, file_land


def find_xml(name, frequency, variable, project='', experiment='', ensemble='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_name, file_area, file_land = find_xml_obs(name, variable)
    else:
        file_name, file_area, file_land = find_xml_cmip(name, project, experiment, ensemble, frequency, realm, variable)
    return file_name, file_area, file_land


def find_xml_cmip(model, project, experiment, ensemble, frequency, realm, variable):
    try:
        pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                             rea=realm, var=variable)
    except:
        if realm == 'O':
            new_realm = 'A'
        elif realm == 'A':
            new_realm = 'O'
        # if var is not in realm 'O' (for ocean), look for it in realm 'A' (for atmosphere), and conversely
        try:
            pathnc, filenc = find_path_and_files(ens=ensemble, exp=experiment, fre=frequency, mod=model, pro=project,
                                                 rea=new_realm, var=variable)
        except:
            # given variable is neither in realm 'A' nor 'O'
            print bcolors.FAIL + '%%%%%     -----     %%%%%'
            print 'ERROR: function: ' + str(INSPECTstack()[0][3]) + ', line: ' + str(INSPECTstack()[0][2])
            print 'given variable cannot be found in either realm A or O: ' + str(variable)
            print 'param: ' + str(model) + ', ' + str(project) + ', ' + str(experiment) + ', ' + str(ensemble) + \
                  ', ' + str(frequency) + ', ' + str(realm)
            print '%%%%%     -----     %%%%%' + bcolors.ENDC
            sys.exit('')
        file_area, file_land = find_fx(model, project=project, experiment=experiment, ensemble=ensemble,
                                       realm=new_realm)
    else:
        file_area, file_land = find_fx(model, project=project, experiment=experiment, ensemble=ensemble, realm=realm)
    file_name = OSpath__join(pathnc, filenc[0])
    return file_name, file_area, file_land


def find_xml_fx(name, project='', experiment='', realm=''):
    list_obs = ReferenceObservations().keys()
    if name in list_obs:
        file_area = OSpath__join(xmldir, 'obs_' + str(name) + '_areacell.xml')
        file_land = OSpath__join(xmldir, 'obs_' + str(name) + '_landmask.xml')
    else:
        filename = str(name) + '_' + str(project) + '_' + str(experiment) + '_r0i0p0_glob_fx_' + str(realm)
        file_area = OSpath__join(xmldir, filename + '_areacell.xml')
        file_land = OSpath__join(xmldir, filename + '_landmask.xml')
    return file_area, file_land


def find_xml_obs(obs, variable):
    if obs == 'HadISST':
        file_name = OSpath__join(xmldir, 'obs_' + str(obs) + 'v1.1.xml')
    else:
        file_name = OSpath__join(xmldir, 'obs_' + str(obs) + '.xml')
    xml = CDMS2open(file_name)
    listvar1 = sorted(xml.listvariables())
    if variable not in listvar1:
        print '\033[95m' + str().ljust(5) + "obs var " + str(variable) + " cannot be found" + '\033[0m'
        print '\033[95m' + str().ljust(10) + "file_name = " + str(file_name) + '\033[0m'
        print '\033[95m' + str().ljust(10) + "variables = " + str(listvar1) + '\033[0m'
        sys.exit('')
    file_area, file_land = find_fx(obs)
    return file_name, file_area, file_land
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# functions to find the number of simulated years available for a model
def nbryear_from_filename(filename):
    length = get_time_size(filename)
    nbryear = length / 12.
    if nbryear == int(nbryear):
        nbryear = int(nbryear)
    else:
        nbryear = MATHfloor(nbryear)
    return nbryear


def nbryear_from_model(experiment, ensemble, frequency, model, project, realm, variable):
    dict_mod = file_model(experiment, ensemble, frequency, model, project, realm, [variable])
    filename = dict_mod[dict_mod.keys()[0]][dict_mod[dict_mod.keys()[0]].keys()[0]]['path + filename']
    nbryear = nbryear_from_filename(filename)
    return nbryear
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# additional functions to (re)save outputs of the ENSO_metrics package
def attributes_global(xml, name, att_dict={}):
    list_glob = xml.listglobal()
    for att in ['cdms_filemap', 'history']:
        while att in list_glob:
            list_glob.remove(att)
    for att in list_glob:
        if 'metric_value' in att:
            att2 = att + name
        else:
            att2 = deepcopy(att)
        att_dict[att2] = xml.attributes[att]
    return


def attributes_variable(xml, variable_name):
    list_att = xml.listattribute(variable_name)
    for att in ['_FillValue','coordinates','history','missing_value','var_desc']:
        while att in list_att:
            list_att.remove(att)
    attributes = dict((att,xml.getattribute(variable_name, att)) for att in list_att)
    return attributes


def save_json(dict_in, json_name):
    list_sli = sorted(dict_in.keys())
    list_ref = sorted(dict_in[list_sli[0]]['value']['metric'].keys())
    list_dat = sorted(dict_in[list_sli[0]]['value']['diagnostic'].keys())
    dict_metric, dict_metric_err = dict(), dict()
    for ref in list_ref:
        tab1, tab2 = list(), list()
        for sli in list_sli:
            tab1.append(dict_in[sli]['value']['metric'][ref]['value'])
            tab2.append(dict_in[sli]['value']['metric'][ref]['value_error'])
        dict_metric[ref], dict_metric_err[ref] = tab1, tab2
    dict_diag, dict_diag_err = dict(), dict()
    for dat in list_dat:
        tab1, tab2 = list(), list()
        for sli in list_sli:
            tab1.append(dict_in[sli]['value']['diagnostic'][dat]['value'])
            tab2.append(dict_in[sli]['value']['diagnostic'][dat]['value_error'])
        dict_diag[dat], dict_diag_err[dat] = tab1, tab2
    dict_metadata = dict_in[list_sli[0]]['metadata']
    dict_out = {'metric': dict_metric, 'metric_err': dict_metric_err, 'diagnostic': dict_diag,
                'diagnostic_err': dict_diag_err, 'metadata': dict_metadata}
    # save as json file
    with open(json_name + '.json', 'w') as outfile:
        json.dump(dict_out, outfile)
    return


def save_netcdf(netcdf_name, metric):
    # list all files
    list_files = sorted(list(GLOBiglob(netcdf_name + '_*' + metric + '.nc')))
    # loop on files
    dict_att_glob, dict_att_var, dict_var = dict(), dict(), dict()
    for file1 in list_files:
        # find slice number
        snbr = file1.split('_slice_')[1].slpit('_to_')[0]
        # open current file
        ff = CDMS2open(file1)
        # read global attributes
        attributes_global(ff, snbr, att_dict=dict_att_glob)
        # list all variables in current file that are not dimensions
        list_variables = [var for var in ff.listvariables() if var not in dimensions]
        list_variables.sort(key=lambda v: v.lower())
        # loop on variables
        for var in list_variables:
            # read variable in the current file
            dict_var[var + '__' + str(snbr)] = ff(var)
            # read variable attributes
            dict_att_var[var + '__' + str(snbr)] = attributes_variable(ff, var)
        ff.close()
        del ff, list_variables, snbr
    # save files
    o = CDMS2open(netcdf_name + '_' + metric + '.nc', 'w+')
    for var in dict_var.keys():
        o.write(dict_var[var], attributes=dict_att_var[var], dtype='float32', id=var)
    for att in dict_att_glob.keys():
        o.__setattr__(att, dict_att_glob[att])
    o.close()
    # delete selected files
    for file1 in list_files:
        OSremove(file1)
    del dict_att_glob, dict_att_var, dict_var, list_files, o
    return
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# main function to compute ENSO_metrics
def main_compute(metricCollection, metric, nbr_years, path, file_name, experiment, frequency, model, project, realm):
    # metric collection dictionary for the given metric collection (mc_name)
    dict_mc = defCollection(metricCollection)
    # list of variables for the given metric
    list_variables = sorted(list(set(dict_mc['metrics_list'][metric]['variables'])))
    print '\033[95m' + 'variables: ' + str(list_variables) + '\033[0m'
    # list of observations
    list_observations = ['GPCPv2.3', 'HadISST', 'Tropflux']
    print '\033[95m' + 'observations: ' + str(list_observations) + '\033[0m'
    # list of ensembles for the given experiment/frequency/model/project/realm
    list_ensembles = sorted(get_ensembles(exp=experiment, fre=frequency, mod=model, pro=project, rea=realm))
    print '\033[95m' + 'ensembles: ' + str(list_ensembles) + '\033[0m'
    # dictionary for observations
    dict_obs = file_obs(list_observations, list_variables, frequency)
    for ens in list_ensembles:
        dict_mod = file_model(experiment, ens, frequency, model, project, realm, list_variables)
        model_file_name = dict_mod[dict_mod.keys()[0]][dict_mod[dict_mod.keys()[0]].keys()[0]]['path + filename']
        model_nbr_years = nbryear_from_filename(model_file_name)
        final_name_out = file_name + '_' + ens + '_' + '_' + str(nbr_years.zfill(3)) + 'years'
        dictDatasets = {'model': dict_mod, 'observations': dict_obs}
        dict1 = dict()
        for ystart in range(1, model_nbr_years-nbr_years):
            yend = ystart+nbr_years
            period = slice(ystart*12, (yend)*12)
            # file names
            file_out = final_name_out + '_slice_' + str(str(ystart).zfill(3)) + '_to_' + str(str(yend).zfill(3))
            file_out = OSpath__join(path, file_out)
            dict1[str(ystart).zfill(3)] = InternCompute(metricCollection, metric, dictDatasets, debug=False,
                                                        netcdf=True, netcdf_name=file_out, period=period)
            del file_out, period, yend
        # save json
        save_json(dict1, final_name_out)
        # resave NetCDF
        save_netcdf(final_name_out, metric)
        del dict_mod, dict1, dictDatasets, model_file_name, model_nbr_years
    del dict_obs
    return
# ---------------------------------------------------------------------------------------------------------------------#
