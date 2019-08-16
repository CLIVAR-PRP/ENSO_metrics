# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoMetricsLib import BiasPrLatRmse, BiasPrLonRmse, BiasPrRmse, BiasSstLonRmse, BiasSstLatRmse, BiasSstSkLonRmse,\
    BiasSstRmse, BiasTauxLatRmse, BiasTauxLonRmse, BiasTauxRmse, EnsoAmpl, EnsoDiversity, EnsodSstOce, EnsoDuration,\
    EnsoFbSshSst, EnsoFbSstLhf, EnsoFbSstLwr, EnsoFbSstShf, EnsoFbSstSwr, EnsoFbSstTaux, EnsoFbSstThf, EnsoFbTauxSsh,\
    EnsoPrMap, EnsoPrJjaTel, EnsoPrNdjTel, EnsoPrTsRmse, EnsoSeasonality, EnsoSlpMap, EnsoSstLonRmse, EnsoSstMap,\
    EnsoSstSkew, EnsoSstTsRmse, EnsoTauxTsRmse, NinaPrJjaTel, NinaPrNdjTel, NinaPrMap, NinaSlpMap, NinaSstDiv,\
    NinaSstDivRmse, NinaSstDur, NinaSstLonRmse, NinaSstMap, NinaSstTsRmse, NinoPrJjaTel, NinoPrNdjTel, NinoPrMap,\
    NinoSlpMap, NinoSstDiv, NinoSstDiversity, NinoSstDivRmse, NinoSstDur, NinoSstLonRmse, NinoSstMap, NinoSstTsRmse,\
    SeasonalPrLatRmse, SeasonalPrLonRmse, SeasonalSstLatRmse, SeasonalSstLonRmse, SeasonalTauxLatRmse,\
    SeasonalTauxLonRmse
from KeyArgLib import DefaultArgValues


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric collection
#
def ComputeCollection(metricCollection, dictDatasets, modelName, user_regridding={}, debug=False, dive_down=False,
                      netcdf=False, netcdf_name=''):
    """
    The ComputeCollection() function computes all the diagnostics / metrics associated with the given Metric Collection

    Inputs:
    ------
    :param MetricCollection: string
        name of a Metric Collection, must be defined in EnsoCollectionsLib.defCollection()
    :param dictDatasets: dict
        dictionary containing all information needed to compute the Metric Collection for one model and observations
        it must be like:
        dictDatasets = {
            'model': {
                'modelName': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                },
            },
            'observations': {
                # obsName1 can be 'HadISST' if all variables come from HadISST, or it could be 'ERA-Interim + HadISST'
                # if the variables come from ERA-Interim AND HadISST,...
                'obsName1': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                    'variable n': {
                        'path + filename': ['path_to_file1/filename1', 'path_to_file2/filename2', ...],
                        'varname': ['variable_name_in_file1', 'variable_name_in_file2', ...],
                        # this last one is not compulsory if 'obsName1' is defined in
                        # EnsoCollectionsLib.ReferenceObservations()
                        'algebric_calculation': [+1,-1,..., '/'],

                    },
                },
                'obsName2': {
                    'variable1': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    'variable2': {'path + filename': 'path_to_file/filename', 'varname': 'variable_name_in_file'},
                    ...
                },
            },
        }
    :param user_regridding: dict, optional
        regridding parameters selected by the user and not the "climate experts" (i.e. the program)
        to see parameters:
        help(EnsoUvcdatToolsLib.TwoVarRegrid)
        help(EnsoUvcdatToolsLib.Regrid)
        e.g.:
        user_regridding = {
            'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                'newgrid_name': 'generic 1x1deg'},
        }
    :param debug: boolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param dive_down: boolean, optional
        default value = False dive_down are not saved in a dictionary
        If you want to save the dive down diagnostics set it to True
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' root name of the saved NetCDFs
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='USER_DATE_METRICCOLLECTION_MODEL'

    :return: MCvalues: dict
        name of the Metric Collection, Metrics, value, value_error, units, ...
        MCvalues = {
            'MetricCollection': {
                'information about the MetricCollection': 'descriptions',
                'metrics': {
                    'metric1': {
                        'metric_values': {
                            'obsName1': {
                                'value': 'value of the metric',
                                'value_error': 'estimation of the error on the metric',
                            },
                            'obsName2': {'value': val, 'value_error': err},
                            ...
                        },
                        'information about the metric': 'description of how this metric if computed from model and
                                                        observations values',
                        'raw_values': {
                            'model': {
                                'value': 'model value of the diagnostic',
                                'value_error': 'estimation of the error on the diagnostic',
                                'nyears': 'number of years used for the computation of the diagnostic',
                                'time_period': 'period used for the computation of the diagnostic',
                            },
                            'observations': {
                                'obsName1': {'value': val, 'value_error': err, 'nyears': ny, 'time_period': bnds},
                                'obsName2': {'value': val, 'value_error': err, 'nyears': ny, 'time_period': bnds},
                                ...
                            },
                            'units': 'units of the diagnostic', 'method': 'method used to compute the diagnostic',
                            'time_frequency': 'data frequency used to compute the diagnostic',
                            'ref': 'reference paper', ...
                        },
                    },
                    'metric2': {
                        ...
                    },
                    ...
                },
            },
        }
    """
    dict_mc = defCollection(metricCollection)
    dict_col_meta = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_col_dd_meta = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_col_valu = dict()
    dict_col_dd_valu = dict()
    dict_m = dict_mc['metrics_list']
    list_metrics = sorted(dict_m.keys(), key=lambda v: v.upper())
    for metric in list_metrics:
        print '\033[94m' + str().ljust(5) + "ComputeCollection: metric = " + str(metric) + '\033[0m'
        # sets arguments for this metric
        list_variables = dict_m[metric]['variables']
        dict_regions = dict_m[metric]['regions']
        # model name, file, variable name in file
        try:
            modelFile1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename']
        except:
            modelFile1 = ''
        try:
            modelVarName1 = dictDatasets['model'][modelName][list_variables[0]]['varname']
        except:
            modelVarName1 = ''
        try:
            modelFileArea1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename_area']
        except:
            modelFileArea1, modelAreaName1 = None, None
        else:
            modelAreaName1 = dictDatasets['model'][modelName][list_variables[0]]['areaname']
        try:
            modelFileLandmask1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename_landmask']
        except:
            modelFileLandmask1, modelLandmaskName1 = None, None
        else:
            modelLandmaskName1 = dictDatasets['model'][modelName][list_variables[0]]['landmaskname']
        # observations name(s), file(s), variable(s) name in file(s)
        obsNameVar1, obsFile1, obsVarName1, obsFileArea1, obsAreaName1 = list(), list(), list(), list(), list()
        obsFileLandmask1, obsLandmaskName1 = list(), list()
        for obs in dictDatasets['observations'].keys():
            try: dictDatasets['observations'][obs][list_variables[0]]
            except: pass
            else:
                obsNameVar1.append(obs)
                obsFile1.append(dictDatasets['observations'][obs][list_variables[0]]['path + filename'])
                obsVarName1.append(dictDatasets['observations'][obs][list_variables[0]]['varname'])
                try:
                    obsFileArea1.append(dictDatasets['observations'][obs][list_variables[0]]['path + filename_area'])
                except:
                    obsFileArea1.append(None)
                    obsAreaName1.append(None)
                else:
                    obsAreaName1.append(dictDatasets['observations'][obs][list_variables[0]]['areaname'])
                try:
                    obsFileLandmask1.append(
                        dictDatasets['observations'][obs][list_variables[0]]['path + filename_landmask'])
                except:
                    obsFileLandmask1.append(None)
                    obsLandmaskName1.append(None)
                else:
                    obsLandmaskName1.append(dictDatasets['observations'][obs][list_variables[0]]['landmaskname'])
        # same if a second variable is needed
        # this time in the form of a keyarg dictionary
        arg_var2 = {'modelFileArea1': modelFileArea1, 'modelAreaName1': modelAreaName1,
                    'modelFileLandmask1': modelFileLandmask1, 'modelLandmaskName1': modelLandmaskName1,
                    'obsFileArea1': obsFileArea1, 'obsAreaName1': obsAreaName1, 'obsFileLandmask1': obsFileLandmask1,
                    'obsLandmaskName1': obsLandmaskName1}
        if len(list_variables) > 1:
            try:
                arg_var2['modelFile2'] = dictDatasets['model'][modelName][list_variables[1]]['path + filename']
            except:
                arg_var2['modelFile2'] = ''
            try:
                arg_var2['modelVarName2'] = dictDatasets['model'][modelName][list_variables[1]]['varname']
            except:
                arg_var2['modelVarName2'] = ''
            arg_var2['regionVar2'] = dict_regions[list_variables[1]]
            try:
                arg_var2['modelFileArea2'] = dictDatasets['model'][modelName][list_variables[1]]['path + filename_area']
            except:
                arg_var2['modelFileArea2'], arg_var2['modelAreaName2'] = None, None
            else:
                arg_var2['modelAreaName2'] = dictDatasets['model'][modelName][list_variables[1]]['areaname']
            try:
                arg_var2['modelFileLandmask2'] = \
                    dictDatasets['model'][modelName][list_variables[1]]['path + filename_landmask']
            except:
                arg_var2['modelFileLandmask2'], arg_var2['modelLandmaskName2'] = None, None
            else:
                arg_var2['modelLandmaskName2'] = dictDatasets['model'][modelName][list_variables[1]]['landmaskname']
            obsNameVar2, obsFile2, obsVarName2, obsFileArea2, obsAreaName2 = list(), list(), list(), list(), list()
            obsFileLandmask2, obsLandmaskName2 = list(), list()
            for obs in dictDatasets['observations'].keys():
                try: dictDatasets['observations'][obs][list_variables[1]]
                except: pass
                else:
                    obsNameVar2.append(obs)
                    obsFile2.append(dictDatasets['observations'][obs][list_variables[1]]['path + filename'])
                    obsVarName2.append(dictDatasets['observations'][obs][list_variables[1]]['varname'])
                    try:
                        obsFileArea2.append(
                            dictDatasets['observations'][obs][list_variables[1]]['path + filename_area'])
                    except:
                        obsFileArea2.append(None)
                        obsAreaName2.append(None)
                    else:
                        obsAreaName2.append(dictDatasets['observations'][obs][list_variables[1]]['areaname'])
                    try:
                        obsFileLandmask2.append(
                            dictDatasets['observations'][obs][list_variables[1]]['path + filename_landmask'])
                    except:
                        obsFileLandmask2.append(None)
                        obsLandmaskName2.append(None)
                    else:
                        obsLandmaskName2.append(dictDatasets['observations'][obs][list_variables[1]]['landmaskname'])
            arg_var2['obsNameVar2'] = obsNameVar2
            arg_var2['obsFile2'] = obsFile2
            arg_var2['obsVarName2'] = obsVarName2
            arg_var2['obsFileArea2'] = obsFileArea2
            arg_var2['obsAreaName2'] = obsAreaName2
            arg_var2['obsFileLandmask2'] = obsFileLandmask2
            arg_var2['obsLandmaskName2'] = obsLandmaskName2
        # computes the metric
        if len(modelFile1) == 0 or (len(list_variables) > 1 and len(arg_var2['modelFile2']) == 0):
            print '\033[94m' + str().ljust(5) + "ComputeCollection: " + str(metricCollection) + ", metric " \
                  + str(metric) + " not computed" + '\033[0m'
            print '\033[94m' + str().ljust(10) + "reason(s):" + '\033[0m'
            if len(modelFile1) == 0:
                print '\033[94m' + str().ljust(11) + "no modeled " + list_variables[0] + " given" + '\033[0m'
            if len(list_variables) > 1 and len(arg_var2['modelFile2']) == 0:
                print '\033[94m' + str().ljust(11) + "no modeled " + list_variables[1] + " given" + '\033[0m'
        elif len(obsFile1) == 0 or (len(list_variables) > 1 and len(arg_var2['obsFile2']) == 0):
            print '\033[94m' + str().ljust(5) + "ComputeCollection: " + str(metricCollection) + ", metric "\
                  + str(metric) + " not computed" + '\033[0m'
            print '\033[94m' + str().ljust(10) + "reason(s):" + '\033[0m'
            if len(obsFile1) == 0:
                print '\033[94m' + str().ljust(11) + "no observed " + list_variables[0] + " given" + '\033[0m'
            if len(list_variables) > 1 and len(arg_var2['obsFile2']) == 0:
                print '\033[94m' + str().ljust(11) + "no observed " + list_variables[1] + " given" + '\033[0m'
        else:
            valu, vame, dive, dime = ComputeMetric(
                metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                dict_regions[list_variables[0]], user_regridding=user_regridding, debug=debug, netcdf=netcdf,
                netcdf_name=netcdf_name, **arg_var2)
            keys1 = valu.keys()
            keys2 = list(set([kk.replace('value', '').replace('__', '').replace('_error', '')
                              for ll in valu[keys1[0]].keys() for kk in valu[keys1[0]][ll].keys()]))
            if len(keys2) > 1:
                for kk in keys2:
                    mm, dd = dict(), dict()
                    keys3 = valu['metric'].keys()
                    for ll in keys3:
                        mm[ll] = {'value': valu['metric'][ll][kk + '__value'],
                                  'value_error': valu['metric'][ll][kk + '__value_error']}
                    keys3 = valu['diagnostic'].keys()
                    for ll in keys3:
                        dd[ll] = {'value': valu['diagnostic'][ll][kk + '__value'],
                                  'value_error': valu['diagnostic'][ll][kk + '__value_error']}
                    dict_col_valu[metric + kk] = {'metric': mm, 'diagnostic': dd}
                    mm = dict((ii, vame['metric'][ii]) for ii in vame['metric'].keys() if 'units' not in ii)
                    mm['units'] = vame['metric'][kk + '__units']
                    dict_col_meta['metrics'][metric + kk] = {'metric': mm, 'diagnostic': vame['diagnostic']}
                    dict_col_dd_valu[metric + kk], dict_col_dd_meta['metrics'][metric + kk] = dive, dime
                    del mm, dd
            else:
                dict_col_valu[metric], dict_col_meta['metrics'][metric] = valu, vame
                dict_col_dd_valu[metric], dict_col_dd_meta['metrics'][metric] = dive, dime
    if dive_down is True:
        return {'value': dict_col_valu, 'metadata': dict_col_meta},\
               {'value': dict_col_dd_valu, 'metadata': dict_col_dd_meta}
    else:
        return {'value': dict_col_valu, 'metadata': dict_col_meta}, {}
# ---------------------------------------------------------------------------------------------------------------------#



# ---------------------------------------------------------------------------------------------------------------------#
#
# Mathematical definition of ENSO metrics
#
def MathMetriComputation(model, model_err, obs=None, obs_err=None, keyword='difference'):
    """

    :param model:
    :param obs:
    :param keyword_metric_computation:
    :return:
    """
    if keyword not in ['difference', 'ratio', 'relative_difference']:
        metric, metric_err, description_metric =\
            None, None, "unknown keyword for the mathematical computation of the metric: " + str(keyword)
        list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": keyword",
                        str().ljust(5) + description_metric]
        EnsoErrorsWarnings.MyWarning(list_strings)
    else:
        if model is not None and obs is not None:
            if keyword == 'difference':
                description_metric =\
                    "The metric is the difference between model and observations values (M = model - obs)"
                metric = model - obs
            elif keyword == 'ratio':
                description_metric = "The metric is the ratio between model and observations values (M = model / obs)"
                metric = model / obs
            else:
                description_metric = \
                    "The metric is the relative difference between model and observations values (M = [model-obs] / obs)"
                metric = (model - obs) / obs
        else:
            metric, description_metric = None, ''
        if model_err is not None or obs_err is not None:
            if keyword == 'difference':
                # mathematical definition of the error on addition / subtraction
                metric_err = model_err + obs_err
            elif keyword == 'ratio':
                # mathematical definition of the error on division
                if model is not None and obs is not None:
                    metric_err = float((obs * model_err + model * obs_err) / NUMPYsquare(obs))
                else:
                    metric_err = None
            else:
                # mathematical definition of the error on division
                if model is not None and obs is not None:
                    metric_err = float((obs * (model_err + obs_err) + (model-obs) * obs_err) / NUMPYsquare(obs))
                else:
                    metric_err = None
        else:
            metric_err = None
    return metric, metric_err, description_metric
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric
#
dict_oneVar_modelAndObs = {
    'BiasPrLatRmse': BiasPrLatRmse, 'BiasPrLonRmse': BiasPrLonRmse, 'BiasPrRmse': BiasPrRmse,
    'BiasSstLatRmse': BiasSstLatRmse, 'BiasSstLonRmse': BiasSstLonRmse, 'BiasSstRmse': BiasSstRmse,
    'BiasTauxLatRmse': BiasTauxLatRmse, 'BiasTauxLonRmse': BiasTauxLonRmse, 'BiasTauxRmse': BiasTauxRmse,
    'BiasSstSkLonRmse': BiasSstSkLonRmse,
    'EnsoSstLonRmse': EnsoSstLonRmse, 'NinaSstLonRmse': NinaSstLonRmse, 'NinoSstTsRmse': NinoSstTsRmse,
    'EnsoSstMap': EnsoSstMap, 'NinaSstMap': NinaSstMap, 'NinoSstMap': NinoSstMap,
    'EnsoSstTsRmse': EnsoSstTsRmse, 'NinaSstTsRmse': NinaSstTsRmse, 'NinoSstLonRmse': NinoSstLonRmse,
    'NinaSstDivRmse': NinaSstDivRmse, 'NinoSstDivRmse': NinoSstDivRmse,
    'SeasonalPrLatRmse': SeasonalPrLatRmse, 'SeasonalPrLonRmse': SeasonalPrLonRmse,
    'SeasonalSstLatRmse': SeasonalSstLatRmse, 'SeasonalSstLonRmse': SeasonalSstLonRmse,
    'SeasonalTauxLatRmse': SeasonalTauxLatRmse, 'SeasonalTauxLonRmse': SeasonalTauxLonRmse,
}

dict_twoVar_modelAndObs = {
    'EnsoPrMap': EnsoPrMap, 'EnsoPrJjaTel': EnsoPrJjaTel, 'EnsoPrNdjTel': EnsoPrNdjTel, 'EnsoSlpMap': EnsoSlpMap,
    'NinaPrMap': NinaPrMap, 'NinaPrJjaTel': NinaPrJjaTel, 'NinaPrNdjTel': NinaPrNdjTel, 'NinaSlpMap': NinaSlpMap,
    'NinoPrMap': NinoPrMap, 'NinoPrJjaTel': NinoPrJjaTel, 'NinoPrNdjTel': NinoPrNdjTel, 'NinoSlpMap': NinoSlpMap,
    'EnsoPrTsRmse': EnsoPrTsRmse, 'EnsoTauxTsRmse': EnsoTauxTsRmse,
}

dict_oneVar = {
    'EnsoAmpl': EnsoAmpl, 'EnsoDuration': EnsoDuration, 'EnsoDiversity': EnsoDiversity,
    'EnsoSeasonality': EnsoSeasonality, 'EnsoSstSkew': EnsoSstSkew, 'NinaSstDiv': NinaSstDiv, 'NinaSstDur': NinaSstDur,
    'NinoSstDiv': NinoSstDiv, 'NinoSstDiversity': NinoSstDiversity, 'NinoSstDur': NinoSstDur,
}

dict_twoVar = {
    'EnsoFbSshSst': EnsoFbSshSst, 'EnsoFbSstLhf': EnsoFbSstLhf, 'EnsoFbSstLwr': EnsoFbSstLwr,
    'EnsoFbSstShf': EnsoFbSstShf, 'EnsoFbSstSwr': EnsoFbSstSwr, 'EnsoFbSstTaux': EnsoFbSstTaux,
    'EnsoFbSstThf': EnsoFbSstThf, 'EnsoFbTauxSsh': EnsoFbTauxSsh, 'EnsodSstOce': EnsodSstOce,
}


def ComputeMetric(metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                  regionVar1, modelFileArea1='', modelAreaName1='', modelFileLandmask1='', modelLandmaskName1='',
                  obsFileArea1='', obsAreaName1='', obsFileLandmask1='', obsLandmaskName1='',
                  modelFile2='', modelVarName2='', modelFileArea2='', modelAreaName2='', modelFileLandmask2='',
                  modelLandmaskName2='', obsNameVar2='', obsFile2='', obsVarName2='', obsFileArea2='', obsAreaName2='',
                  obsFileLandmask2='', obsLandmaskName2='', regionVar2='', user_regridding={}, debug=False,
                  netcdf=False, netcdf_name=''):
    """
    :param metricCollection: string
        name of a Metric Collection, must be defined in EnsoCollectionsLib.defCollection()
    :param metric: string
        name of a Metric, must be defined in EnsoCollectionsLib.defCollection()
    :param modelName: string
        path_to/filename of a model
    :param modelFile1: string or list of strings
        model file number 1
        e.g.: 'model_file.nc' or ['model_file1.nc','model_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param modelVarName1: string or list of strings
        model variable number 1, must be defined in EnsoCollectionsLib.CmipVariables()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param obsNameVar1: string
        path_to/filename of the observations for variable number 1, must be defined in
        EnsoCollectionsLib.ReferenceObservations()
    :param obsFile1: string or list of strings
        observations file number 1
        e.g.: 'observations_file.nc' or ['observations_file1.nc','observations_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param obsVarName1: string or list of strings
        observations variable number 1, must be defined in EnsoCollectionsLib.ReferenceObservations()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param regionVar1: string
        name of box ('tropical_pacific') for variable number 1, must be defined in EnsoCollectionsLib.ReferenceRegions()
    :param modelFileArea1: string, optional
        path_to/filename of the model areacell for variable number 1
    :param modelAreaName1: string, optional
        name the model areacell (e.g. 'areacella' or 'areacello') for variable number 1
    :param modelFileLandmask1: string, optional
        path_to/filename of the model landmask for variable number 1
    :param modelLandmaskName1: string, optional
        name the model landmask (e.g. 'landmask' or 'sftlf') for variable number 1
    :param obsFileArea1: string, optional
        path_to/filename of the observations areacell for variable number 1
    :param obsAreaName1: string, optional
        name the observations areacell (e.g. 'areacella' or 'areacello') for variable number 1
    :param obsFileLandmask1: string, optional
        path_to/filename of the observations landmask for variable number 1
    :param obsLandmaskName1: string, optional
        name the observations landmask (e.g. 'landmask' or 'sftlf') for variable number 1
    :param modelFile2: string or list of strings, optional
        model file number 2
        e.g.: 'model_file.nc' or ['model_file1.nc','model_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param modelVarName2: string or list of strings, optional
        model variable number 2, must be defined in EnsoCollectionsLib.CmipVariables()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param modelFileArea2: string, optional
        path_to/filename of the model areacell for variable number 2
    :param modelAreaName2: string, optional
        name the model areacell (e.g. 'areacella' or 'areacello') for variable number 2
    :param modelFileLandmask2: string, optional
        path_to/filename of the model landmask for variable number 2
    :param modelLandmaskName2: string, optional
        name the model landmask (e.g. 'landmask' or 'sftlf') for variable number 2
    :param obsNameVar2: string, optional
        path_to/filename of the observations for variable number 2, must be defined in
        EnsoCollectionsLib.ReferenceObservations()
    :param obsFile2: string or list of strings, optional
        observations file number 2
        e.g.: 'observations_file.nc' or ['observations_file1.nc','observations_file2.nc']
        a list of files is given if a variable is the result of more than one CMIP variable (lwr = rlds - rlus)
    :param obsVarName2: string or list of strings, optional
        observations variable number 2, must be defined in EnsoCollectionsLib.ReferenceObservations()
        e.g.: 'lwr' or ['rlds', 'rlus']
    :param obsFileArea2: string, optional
        path_to/filename of the observations areacell for variable number 2
    :param obsAreaName2: string, optional
        name the observations areacell (e.g. 'areacella' or 'areacello') for variable number 2
    :param obsFileLandmask2: string, optional
        path_to/filename of the observations landmask for variable number 2
    :param obsLandmaskName2: string, optional
        name the observations landmask (e.g. 'landmask' or 'sftlf') for variable number 2
    :param regionVar2: string
        name of box ('tropical_pacific') for variable number 2, must be defined in EnsoCollectionsLib.ReferenceRegions()
    :param user_regridding: dict, optional
        regridding parameters selected by the user and not the "climate experts" (i.e. the program)
        to see parameters:
        help(EnsoUvcdatToolsLib.TwoVarRegrid)
        help(EnsoUvcdatToolsLib.Regrid)
        e.g.:
        user_regridding = {
            'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                'newgrid_name': 'generic 1x1deg'},
        }
    :param debug: boolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param dive_down: boolean, optional
        default value = False dive_down are not saved in a dictionary
        If you want to save the dive down diagnostics set it to True
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' root name of the saved NetCDFs
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='USER_DATE_METRICCOLLECTION_MODEL'

    :return:
    """
    tmp_metric = deepcopy(metric)
    metric = metric.replace('_1', '').replace('_2', '').replace('_3', '').replace('_4', '').replace('_5', '')
    # retrieving keyargs from EnsoCollectionsLib.defCollection
    dict_mc = defCollection(metricCollection)

    # common_collection_parameters
    keyarg = dict()
    for arg in dict_mc['common_collection_parameters'].keys():
        keyarg[arg] = dict_mc['common_collection_parameters'][arg]
    for arg in dict_mc['metrics_list'][tmp_metric].keys():
        keyarg[arg] = dict_mc['metrics_list'][tmp_metric][arg]
    # if 'metric_computation' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its
    # default value
    try:
        keyarg['metric_computation']
    except:
        keyarg['metric_computation'] = DefaultArgValues('metric_computation')
    # if 'modeled_period' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg['time_bounds_mod'] = keyarg['modeled_period']
    except:
        keyarg['time_bounds_mod'] = DefaultArgValues('time_bounds_mod')
    # if 'modeled_period' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg['time_bounds_obs'] = keyarg['observed_period']
    except:
        keyarg['time_bounds_obs'] = DefaultArgValues('time_bounds_obs')
    # if the user gave a specific regridding Tool / method, use it
    if tmp_metric in user_regridding.keys():
        keyarg['regridding'] = user_regridding[tmp_metric]
    elif 'regridding' in user_regridding.keys():
        keyarg['regridding'] = user_regridding['regridding']

    # obsName could be a list if the user wants to compare the model with a set of observations
    # if obsName is just a name (string) it is put in a list
    if isinstance(obsNameVar1, basestring):
        obsNameVar1 = [obsNameVar1]
    if isinstance(obsFile1, basestring):
        obsFile1 = [obsFile1]
    if isinstance(obsVarName1, basestring):
        obsVarName1 = [obsVarName1]
    # Var2, do the same as for Var1
    if isinstance(obsNameVar2, basestring):
        obsNameVar2 = [obsNameVar2]
    if isinstance(obsFile2, basestring):
        obsFile2 = [obsFile2]
    if isinstance(obsVarName2, basestring):
        obsVarName2 = [obsVarName2]

    dict_metric_val = dict()
    dict_diagnostic = dict()
    dict_diagnostic_metadata = dict()
    dict_dive_down = dict()
    dict_dive_down_metadata = dict()

    multimetric = False

    if metric in dict_oneVar_modelAndObs.keys() or metric in dict_twoVar_modelAndObs.keys():
        #
        # this part regroups all diagnostics comparing model and obs (rmse)
        # so the diagnostic is the metric
        #
        description_metric = "The metric is the statistical value between the model and the observations"
        diagnostic1 = dict()
        for ii in range(len(obsNameVar1)):
            # computes the diagnostic/metric
            if metric in dict_oneVar_modelAndObs.keys():
                output_name = obsNameVar1[ii]
                print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarRMSmetric, " + metric + " = " + modelName\
                      + " and " + output_name + '\033[0m'
                diagnostic1[output_name] = dict_oneVar_modelAndObs[metric](
                    modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                    obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                    obsLandmaskName1[ii], regionVar1, dataset1=modelName, dataset2=output_name, debug=debug,
                    netcdf=netcdf, netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
            elif metric in dict_twoVar_modelAndObs.keys():
                for jj in range(len(obsNameVar2)):
                    output_name = obsNameVar1[ii] + '_' + obsNameVar2[jj]
                    print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarRMSmetric, " + metric + " = " + modelName\
                          + " and " + output_name + '\033[0m'
                    diagnostic1[output_name] = dict_twoVar_modelAndObs[metric](
                        modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1,
                        modelLandmaskName1, modelFile2, modelVarName2, modelFileArea2, modelAreaName2,
                        modelFileLandmask2, modelLandmaskName2, obsFile1[ii], obsVarName1[ii], obsFileArea1[ii],
                        obsAreaName1[ii], obsFileLandmask1[ii], obsLandmaskName1[ii], obsFile2[jj], obsVarName2[jj],
                        obsFileArea2[jj], obsAreaName2[jj], obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar1,
                        regionVar2, dataset1=modelName, dataset2=output_name, debug=debug, netcdf=netcdf,
                        netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
        for obs in diagnostic1.keys():
            # puts metric values in its proper dictionary
            if 'value' in diagnostic1[obs].keys():
                dict_metric_val[obs] = {
                    'value': diagnostic1[obs]['value'], 'value_error': diagnostic1[obs]['value_error']}
                dict_diagnostic[modelName] = {'value': None, 'value_error': None}
                dict_diagnostic[obs] = {'value': None, 'value_error': None}
            else:
                multimetric = True
                lkeys = list(set([key.split('__')[0] for key in diagnostic1[obs].keys() if 'value' in key]))
                for key in lkeys:
                    try:
                        dict_metric_val[obs]
                    except:
                        dict_metric_val[obs] = {key + '__value': diagnostic1[obs][key + '__value'],
                                                         key + '__value_error': diagnostic1[obs][key + '__value_error']}
                        dict_diagnostic[modelName] = {key + '__value': None, key + '__value_error': None}
                        dict_diagnostic[obs] = {key + '__value': None, key + '__value_error': None}
                    else:
                        dict_metric_val[obs][key + '__value'] = diagnostic1[obs][key + '__value']
                        dict_metric_val[obs][key + '__value_error'] = diagnostic1[obs][key + '__value_error']
                        dict_diagnostic[modelName][key + '__value'] = None
                        dict_diagnostic[modelName][key + '__value_error'] = None
                        dict_diagnostic[obs][key + '__value'] = None
                        dict_diagnostic[obs][key + '__value_error'] = None
            if 'dive_down_diag' in diagnostic1[obs].keys():
                dict_dive_down[modelName] = diagnostic1[obs]['dive_down_diag']['model']
                dict_dive_down[obs] = diagnostic1[obs]['dive_down_diag']['observations']
                dict1 = {}
                for elt in diagnostic1[obs]['dive_down_diag'].keys():
                    if elt not in ['model', 'observations']:
                        dict1[elt] = diagnostic1[obs]['dive_down_diag'][elt]
                dict_dive_down_metadata[obs] = dict1
                del dict1
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata[modelName] = {'name': modelName, 'nyears': diagnostic1[obs]['nyears_model'],
                                                   'time_period': diagnostic1[obs]['time_period_model']}
            dict_diagnostic_metadata[obs] = {'name': obs, 'nyears': diagnostic1[obs]['nyears_observations'],
                                             'time_period': diagnostic1[obs]['time_period_observations']}
            if 'events_model' in diagnostic1[obs].keys():
                dict_diagnostic_metadata[modelName]['events'] = diagnostic1[obs]['events_model']
                dict_diagnostic_metadata[obs]['events'] = diagnostic1[obs]['events_observations']
        units = diagnostic1[obs]['units']
        diagnostic1 = diagnostic1[obs]
    else:
        #
        # model diagnostic
        #
        keyarg['time_bounds'] = keyarg['time_bounds_mod']
        keyarg['project_interpreter_var1'] = keyarg['project_interpreter']
        if metric in dict_oneVar.keys():
            # computes diagnostic that needs only one variable
            print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(modelName) + '\033[0m'
            diagnostic1 = dict_oneVar[metric](
                modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                regionVar1, dataset=modelName, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name, metname=tmp_metric,
                **keyarg)
        elif metric in dict_twoVar.keys():
            # computes diagnostic that needs two variables
            print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(modelName) + '\033[0m'
            keyarg['project_interpreter_var2'] = keyarg['project_interpreter']
            diagnostic1 = dict_twoVar[metric](
                modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                regionVar1, modelFile2, modelVarName2, modelFileArea2, modelAreaName2, modelFileLandmask2,
                modelLandmaskName2, regionVar2, dataset=modelName, debug=debug, netcdf=netcdf, netcdf_name=netcdf_name,
                metname=tmp_metric, **keyarg)
        else:
            diagnostic1 = None
            list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": metric", str().ljust(5) +
                            "unknown metric name: " + str(metric)]
            EnsoErrorsWarnings.MyError(list_strings)
        # puts metric / diagnostic values in its proper dictionary
        dict_diagnostic[modelName] = {'value': diagnostic1['value'], 'value_error': diagnostic1['value_error']}
        if 'nonlinearity' in diagnostic1.keys():
            dict_diagnostic[modelName]['nonlinearity'] = diagnostic1['nonlinearity']
            dict_diagnostic[modelName]['nonlinearity_error'] = diagnostic1['nonlinearity_error']
        if 'dive_down_diag' in diagnostic1.keys():
            dict_dive_down[modelName] = diagnostic1['dive_down_diag']['value']
            for elt in diagnostic1['dive_down_diag'].keys():
                if elt not in ['value']:
                    try:
                        dict_dive_down_metadata[modelName]
                    except:
                        dict_dive_down_metadata[modelName] = {elt: diagnostic1['dive_down_diag'][elt]}
                    else:
                        dict_dive_down_metadata[modelName][elt] = diagnostic1['dive_down_diag'][elt]
        # puts diagnostic metadata in its proper dictionary
        dict_diagnostic_metadata[modelName] = {
            'name': modelName, 'nyears': diagnostic1['nyears'], 'time_period': diagnostic1['time_period'],
        }
        if 'events_model' in diagnostic1.keys():
            dict_diagnostic_metadata[modelName]['events'] = diagnostic1['events']
        #
        # observations diag
        #
        diag_obs = dict()
        keyarg['time_bounds'] = keyarg['time_bounds_obs']
        for ii in range(len(obsNameVar1)):
            # sets observations
            obs1 = obsNameVar1[ii]
            keyarg['project_interpreter_var1'] = obs1
            #            keyarg['project_interpreter'] = 'CMIP'
            if metric in dict_oneVar.keys():
                print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(obs1) + '\033[0m'
                output_name = obs1
                diag_obs[output_name] = dict_oneVar[metric](
                    obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                    obsLandmaskName1[ii], regionVar1, dataset=output_name, debug=debug, netcdf=netcdf,
                    netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
            elif metric in dict_twoVar.keys():
                for jj in range(len(obsNameVar2)):
                    obs2 = obsNameVar2[jj]
                    output_name = obs1 + '_' + obs2
                    keyarg['project_interpreter_var2'] = obs2
                    print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(output_name) + '\033[0m'
                    diag_obs[output_name] = dict_twoVar[metric](
                        obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                        obsLandmaskName1[ii], regionVar1, obsFile2[jj], obsVarName2[jj], obsFileArea2[jj],
                        obsAreaName2[jj], obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar2, dataset=output_name,
                        debug=debug, netcdf=netcdf, netcdf_name=netcdf_name, metname=tmp_metric, **keyarg)
        for obs in diag_obs.keys():
            # computes the metric
            metric_val, metric_err, description_metric = MathMetriComputation(
                diagnostic1['value'], diagnostic1['value_error'], obs=diag_obs[obs]['value'],
                obs_err=diag_obs[obs]['value_error'], keyword=keyarg['metric_computation'])
            dict_metric_val[obs] = {'value': metric_val, 'value_error': metric_err}
            # puts metric / diagnostic values in its proper dictionary
            dict_diagnostic[obs] = {'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error']}
            if 'nonlinearity' in diag_obs[obs].keys():
                dict_diagnostic[obs]['nonlinearity'] = diag_obs[obs]['nonlinearity']
                dict_diagnostic[obs]['nonlinearity_error'] = diag_obs[obs]['nonlinearity_error']
            if 'dive_down_diag' in diag_obs[obs].keys():
                dict_dive_down[obs] = diag_obs[obs]['dive_down_diag']['value']
                for elt in diag_obs[obs]['dive_down_diag'].keys():
                    if elt not in ['value']:
                        try:
                            dict_dive_down_metadata[obs]
                        except:
                            dict_dive_down_metadata[obs] = {elt: diag_obs[obs]['dive_down_diag'][elt]}
                        else:
                            dict_dive_down_metadata[obs][elt] = diag_obs[obs]['dive_down_diag'][elt]
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata[obs] = {
                'name': obs, 'nyears': diag_obs[obs]['nyears'], 'time_period': diag_obs[obs]['time_period'],
            }
            if 'events_model' in diag_obs[obs].keys():
                dict_diagnostic_metadata[obs]['events'] = diag_obs[obs]['events']
        if keyarg['metric_computation'] in ['ratio', 'relative_difference']:
            units = diagnostic1['units'] + ' / ' + diagnostic1['units']
        else:
            units = diagnostic1['units']
    # finishes to fill the diagnostic dictionary
    list_keys = ['method', 'name', 'ref', 'time_frequency', 'units']
    for key in list_keys:
        if key == 'method':
            dict_diagnostic_metadata[key] = diagnostic1['method']
            dict_dive_down_metadata[key] = diagnostic1['method']
            try:
                diagnostic1['nonlinearity']
            except:
                pass
            else:
                dict_diagnostic_metadata['method_nonlinearity'] = diagnostic1['method_nonlinearity']
        else:
            dict_diagnostic_metadata[key] = diagnostic1[key]
            dict_dive_down_metadata[key] = diagnostic1[key]
    # creates the output dictionaries
    dict_metrics = {
        'metric': dict_metric_val, 'diagnostic': dict_diagnostic,
    }
    datasets = modelName + "; "
    for ii in range(len(obsNameVar1)):
        datasets = datasets + obsNameVar1[ii] + "'s " + obsVarName1[ii]
        if ii != len(obsNameVar1) - 1:
            datasets = datasets + " & "
        else:
            datasets = datasets + "; "
    if len(obsNameVar2) > 0:
        for ii in range(len(obsNameVar2)):
            if isinstance(obsVarName2[ii], list):
                datasets = datasets + obsNameVar2[ii] + "'s "
                for jj in range(len(obsVarName2[ii])):
                    datasets += obsVarName2[ii][jj]
                    if jj != len(obsVarName2[ii])-1:
                        datasets += " & "
            else:
                datasets = datasets + obsNameVar2[ii] + "'s " + obsVarName2[ii]
            if ii != len(obsNameVar2) - 1:
                datasets = datasets + " & "
    if multimetric is True:
        tmp = {'name': metric, 'method': description_metric, 'datasets': datasets}
        for key in lkeys:
            tmp[key + '__units'] = diagnostic1[key + '__units']
        dict_metadata = {'metric': tmp, 'diagnostic': dict_diagnostic_metadata}
    else:
        dict_metadata = {
            'metric': {'name': metric, 'method': description_metric, 'datasets': datasets, 'units': units},
            'diagnostic': dict_diagnostic_metadata,
        }
    return dict_metrics, dict_metadata, dict_dive_down, dict_dive_down_metadata
# ---------------------------------------------------------------------------------------------------------------------#
