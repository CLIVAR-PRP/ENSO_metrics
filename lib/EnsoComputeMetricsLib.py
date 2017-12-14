# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import sqrt as NUMPYsqrt
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoMetricsLib import EnsoAlphaLhf, EnsoAlphaLwr, EnsoAlphaSwr, EnsoAlphaThf, EnsoAmpl, EnsoMu, EnsoRMSE,\
    EnsoSeasonality
from EnsoUvcdatToolsLib import CheckTime, Multiply, SeasonalMean, TimeBounds,\
    TimeAnomaliesLinearRegressionAndNonlinearity, TimeAnomaliesStd, ReadSelectRegionCheckUnits,\
    SpatialRms


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
    if keyword == 'difference':
        description_metric = "The metric is the difference between model and observations values (M = model - obs)"
        metric = model-obs
        # mathematical definition of the error on addition / subtraction
        metric_err = model_err + obs_err
    elif keyword == 'ratio':
        description_metric = "The metric is the ratio between model and observations values (M = model / obs)"
        metric = model/obs
        # mathematical definition of the error on division
        metric_err = float((obs * model_err + model * obs_err) / NUMPYsquare(obs))
    elif keyword == 'relative_difference':
        description_metric = "The metric is the relative difference between model and observations values (M = [model-obs] / obs)"
        metric = (model-obs)/obs
        # mathematical definition of the error on division
        metric_err = float((obs * (model_err + obs_err) + (model-obs) * obs_err) / NUMPYsquare(obs))
    else:
        metric, metric_err, description_metric = None, None, ''
        list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": keyword", str().ljust(5) +
                        "unknown keyword for the mathematical computation of the metric: " + str(keyword)]
        EnsoErrorsWarnings.MyError(list_strings)
    return metric, metric_err, description_metric
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric
#
dict_oneVar_modelAndObs = {'EnsoRMSE': EnsoRMSE,}

dict_oneVar = {'EnsoAmpl': EnsoAmpl, 'EnsoSeasonality': EnsoSeasonality,}

dict_twoVar = {
    'EnsoAlphaThf': EnsoAlphaThf, 'EnsoAlphaLhf': EnsoAlphaLhf, 'EnsoAlphaLwr': EnsoAlphaLwr,
    'EnsoAlphaSwr': EnsoAlphaSwr, 'EnsoMu': EnsoMu,
}


def ComputeMetric(metricCollection, metric, modelName, modelFile1, modelVarName1, obsName, obsFile1, obsVarName1,
                  regionVar1='', modelFile2='', modelVarName2='', obsFile2='', obsVarName2='', regionVar2=''):
    """
    The ComputeMetric() function computes the given metric for the given model and observations

    :param metricCollection:
    :param metric:
    :param modelName:
    :param modelFile1:
    :param modelVarName1:
    :param obsName:
    :param obsFile1:
    :param obsVarName1:
    :param regionVar1:
    :param regionVar2:
    :param modelFile2:
    :param modelVarName2:
    :param obsFile2:
    :param obsVarName2:
    :return:
    """
    frequency = defCollection(metricCollection)['common_collection_parameters']['frequency']
    mintimesteps = defCollection(metricCollection)['common_collection_parameters']['minimum_number_of_time_steps']
    try:
        observed_period = defCollection(metricCollection)['common_collection_parameters']['observed_period']
    except:
        observed_period = ''
    try:
        modeled_period = defCollection(metricCollection)['common_collection_parameters']['modeled_period']
    except:
        modeled_period = ''

    list_variables = defCollection(metricCollection)['metrics_list'][metric]['variables']
    dict_regions = defCollection(metricCollection)['metrics_list'][metric]['regions']
    try: keyword_metric_computation = defCollection(metricCollection)['metrics_list'][metric]['metric_computation']
    except: keyword_metric_computation = 'difference'

    if not regionVar1:
        regionVar1 = dict_regions[list_variables[0]]
    if len(list_variables) > 1:
        if not regionVar2:
            regionVar2 = dict_regions[list_variables[1]]

    if not isinstance(obsName, list):
        obsName = [obsName]

    dict_metric_val = dict()
    if metric in dict_oneVar_modelAndObs.keys():
        description_metric = "The metric is the statistical value between the model and the observations"
        for ii in range(len(obsName)):
            obs, tmp_obsFile1, tmp_obsVarName1 = obsName[ii], obsFile1[ii], obsVarName1[ii]
            diagnostic1 = dict_oneVar_modelAndObs[metric](modelFile1, modelVarName1, tmp_obsFile1, tmp_obsVarName1,
                                                          regionVar1, timeboundsmodel=modeled_period,
                                                          timeboundsobs=observed_period, centered_rmse=0,
                                                          frequency=frequency, mintimesteps=mintimesteps,
                                                          regrid_model_on_obs=True)
            dict_metric_val['ref_' + obs] = {'value': diagnostic1['value'], 'value_error': diagnostic1['value_error']}
            try:
                dict_diagnostic
            except:
                dict_diagnostic = {
                    'model': {
                        'name': modelName, 'value': None, 'value_error': None, 'nyears': diagnostic1['nyears_model'],
                        'time_period': diagnostic1['time_period_model'],
                    },
                    'observations': {
                        obs: {
                            'name': obs, 'value': None, 'value_error': None,
                            'nyears': diagnostic1['nyears_observations'],
                            'time_period': diagnostic1['time_period_observations'],
                        },
                    },
                }
            else:
                dict_diagnostic['observations'][obs] = {
                    'name': obs, 'value': None, 'value_error': None, 'nyears': diagnostic1['nyears_observations'],
                    'time_period': diagnostic1['time_period_observations'],
                }
    else:
        diag_obs = dict()
        if metric in dict_oneVar.keys():
            diagnostic1 = dict_oneVar[metric](modelFile1, modelVarName1, regionVar1)
        elif metric in dict_twoVar.keys():
            diagnostic1 = dict_twoVar[metric](modelFile1, modelFile2, modelVarName1, modelVarName2, regionVar1,
                                              regionVar2)
        else:
            diagnostic1 = None
            list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": metric", str().ljust(5) +
                            "unknown metric name: " + str(metric)]
            EnsoErrorsWarnings.MyError(list_strings)
        dict_diagnostic = {
            'model': {
                'name': modelName, 'value': diagnostic1['value'], 'value_error': diagnostic1['value_error'],
                'nyears': diagnostic1['nyears'], 'time_period': diagnostic1['time_period'],
            },
        }
        try:
            diagnostic1['nonlinearity']
        except:
            pass
        else:
            dict_diagnostic['model']['nonlinearity'] = diagnostic1['nonlinearity']
            dict_diagnostic['model']['nonlinearity_error'] = diagnostic1['nonlinearity_error']
        for ii in range(len(obsName)):
            obs = obsName[ii]
            if metric in dict_oneVar.keys():
                tmp_obsFile1, tmp_obsVarName1 = obsFile1[ii], obsVarName1[ii]
                diag_obs[obs] = dict_oneVar[metric](tmp_obsFile1, tmp_obsVarName1, regionVar1)
            elif metric in dict_twoVar.keys():
                tmp_obsFile1, tmp_obsVarName1 = obsFile1[ii], obsVarName1[ii]
                tmp_obsFile2, tmp_obsVarName2 = obsFile2[ii], obsVarName2[ii]
                diag_obs[obs] = dict_twoVar[metric](tmp_obsFile1, tmp_obsFile2, tmp_obsVarName1, tmp_obsVarName2,
                                                    regionVar1, regionVar2)
            try:
                dict_diagnostic['observations']
            except:
                dict_diagnostic['observations'] = {
                    obs: {
                        'name': obs, 'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error'],
                        'nyears': diagnostic1['nyears'], 'time_period': diagnostic1['time_period'],
                    },
                }
            else:
                dict_diagnostic['observations'][obs] = {
                    'name': obs, 'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error'],
                    'nyears': diagnostic1['nyears'], 'time_period': diagnostic1['time_period'],
                }
            metric_val, metric_err, descript = MathMetriComputation(diagnostic1['value'],
                                                                              diagnostic1['value_error'],
                                                                              obs=diag_obs[obs]['value'],
                                                                              obs_err=diag_obs[obs]['value_error'],
                                                                              keyword=keyword_metric_computation)
            description_metric = descript
            dict_metric_val['ref_' + obs] = {'value': metric_val, 'value_error': metric_err}
            try:
                diagnostic1['nonlinearity']
            except:
                pass
            else:
                dict_diagnostic['observations'][obs]['nonlinearity'] = diag_obs[obs]['nonlinearity']
                dict_diagnostic['observations'][obs]['nonlinearity_error'] = diag_obs[obs]['nonlinearity_error']
                metric_val, metric_err, string = MathMetriComputation(diagnostic1['nonlinearity'],
                                                                                  diagnostic1['nonlinearity_error'],
                                                                                  obs=diag_obs[obs]['nonlinearity'],
                                                                                  obs_err=diag_obs[obs]['nonlinearity_error'],
                                                                                  keyword=keyword_metric_computation)
                dict_metric_val['ref_' + obs]['nonlinearity'] = metric_val
                dict_metric_val['ref_' + obs]['nonlinearity_error'] = metric_err

    list_keys = ['name','units', 'time_frequency', 'method_to_compute_diagnostic', 'ref']
    for key in list_keys:
        if key == 'method_to_compute_diagnostic':
            dict_diagnostic[key] = diagnostic1['method']
            try:
                diagnostic1['nonlinearity']
            except:
                pass
            else:
                dict_diagnostic['method_to_compute_nonlinearity'] = diagnostic1['method_nonlinearity']
        else:
            dict_diagnostic[key] = diagnostic1[key]

    dict_metrics = {
        'name': metric, 'datasets':[modelName]+obsName, 'metric_values': dict_metric_val,
        'method_to_compute_metric': description_metric, 'diagnostic': dict_diagnostic,
    }
    return dict_metrics
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric collection
#
def ComputeCollection(metricCollection, dictDatasets):
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
    :return: MCvalues: dict
        name of the Metric Collection, Metrics, value, value_error, units, ...
        MCvalues = {
            'MetricCollection': {
                'information about the MetricCollection': 'descriptions',
                'metrics': {
                    'metric1': {
                        'metric_values': {
                            'ref_obsName1': {
                                'value': 'value of the metric',
                                'value_error': 'estimation of the error on the metric',
                            },
                            'ref_obsName2': {'value': val, 'value_error': err},
                            ...
                        },
                        'information about the metric': 'description of how this metric if computed from model and
                                                        observations values',
                        'diagnostic': {
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
    dict_collection = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_m = dict_mc['metrics_list']
    list_metrics = dict_m.keys()
    for metric in list_metrics:
        list_variables = dict_m[metric]['variables']
        dict_regions = dict_m[metric]['regions']
        modelName = dictDatasets['model'].keys()[0]
        modelFile1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename']
        modelVarName1 = dictDatasets['model'][modelName][list_variables[0]]['varname']
        if len(list_variables) > 1:
            modelFile2 = dictDatasets['model'][modelName][list_variables[1]]['path + filename']
            modelVarName2 = dictDatasets['model'][modelName][list_variables[1]]['varname']
            obsFile2, obsVarName2 = list(), list()
        obsName, obsFile1, obsVarName1 = list(), list(), list()
        for obs in dictDatasets['observations'].keys():
            obsName.append(obs)
            obsFile1.append(dictDatasets['observations'][obs][list_variables[0]]['path + filename'])
            obsVarName1.append(dictDatasets['observations'][obs][list_variables[0]]['varname'])
            if len(list_variables) > 1:
                obsFile2.append(dictDatasets['observations'][obs][list_variables[1]]['path + filename'])
                obsVarName2.append(dictDatasets['observations'][obs][list_variables[1]]['varname'])
        if len(list_variables) == 1:
            dict_collection['metrics'][metric] = ComputeMetric(metricCollection, metric, modelName, modelFile1,
                                                               modelVarName1, obsName, obsFile1, obsVarName1,
                                                               regionVar1=dict_regions[list_variables[0]])
        else:
            dict_collection['metrics'][metric] = ComputeMetric(metricCollection, metric, modelName, modelFile1,
                                                               modelVarName1, obsName, obsFile1, obsVarName1,
                                                               regionVar1=dict_regions[list_variables[0]],
                                                               modelFile2=modelFile2, modelVarName2=modelVarName2,
                                                               obsFile2=obsFile2, obsVarName2=obsVarName2,
                                                               regionVar2=dict_regions[list_variables[1]])
    return dict_collection
# ---------------------------------------------------------------------------------------------------------------------#
