# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoMetricsLib import EnsoRmse
from KeyArgLib import DefaultArgValues


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
    dict_collection = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_m = dict_mc['metrics_list']
    list_metrics = sorted(dict_m.keys())
    for metric in list_metrics:
        # sets arguments for this metric
        list_variables = dict_m[metric]['variables']
        dict_regions = dict_m[metric]['regions']
        # model name, file, variable name in file
        modelName = dictDatasets['model'].keys()[0]
        modelFile1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename']
        modelVarName1 = dictDatasets['model'][modelName][list_variables[0]]['varname']
        # observations name(s), file(s), variable(s) name in file(s)
        obsName, obsFile1, obsVarName1 = list(), list(), list()
        for obs in dictDatasets['observations'].keys():
            obsName.append(obs)
            obsFile1.append(dictDatasets['observations'][obs][list_variables[0]]['path + filename'])
            obsVarName1.append(dictDatasets['observations'][obs][list_variables[0]]['varname'])
        # same if a second variable is needed
        # this time in the form of a keyarg dictionary
        arg_var2 = dict()
        if len(list_variables) > 1:
            arg_var2['modelFile2'] = dictDatasets['model'][modelName][list_variables[1]]['path + filename']
            arg_var2['modelVarName2'] = dictDatasets['model'][modelName][list_variables[1]]['varname']
            arg_var2['regionVar2'] = dict_regions[list_variables[1]]
            obsFile2, obsVarName2 = list(), list()
            for obs in dictDatasets['observations'].keys():
                obsFile2.append(dictDatasets['observations'][obs][list_variables[1]]['path + filename'])
                obsVarName2.append(dictDatasets['observations'][obs][list_variables[1]]['varname'])
            arg_var2['obsFile2'] = obsFile2
            arg_var2['obsVarName2'] = obsVarName2
        # computes the metric
        dict_collection['metrics'][metric] = ComputeMetric(
            metricCollection, metric, modelName, modelFile1, modelVarName1, obsName, obsFile1, obsVarName1,
            dict_regions[list_variables[0]], **arg_var2)
    return dict_collection
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
        description_metric = \
            "The metric is the relative difference between model and observations values (M = [model-obs] / obs)"
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
dict_oneVar_modelAndObs = {'EnsoRmse': EnsoRmse}

dict_oneVar = {}

dict_twoVar = {}


def ComputeMetric(metricCollection, metric, modelName, modelVar1, obsNameVar1, obsVar1, **kwargs):
    """
    The ComputeMetric() function computes the given metric for the given model and observations

    :param metricCollection: string
        name of the metric collection (from which the given 'metric' is from)
    :param metric: string
        name of the metric
    :param modelName: string
        name of the model (for which the metric is computed)
    :param modelVar1: masked_array
        masked_array of the first modeled variable of the given metric
    :param obsNameVar1: string or list
        name of the observations used for this first variable
        if several observations are used for this variable, 'obsNameVar1' is a list of observations' name
    :param obsVar1: masked_array or list of masked_array
        masked_array of the first observed variable of the given metric
        if several observations are used for this variable, 'obsNameVar1' is a list of observations' masked_array
    usual kwargs:
    :param modelVar2: masked_array
        masked_array of the second modeled variable of the given metric
    :param obsNameVar2: string or list
        name of the observations used for this second variable
        if several observations are used for this variable, 'obsNameVar1' is a list of observations' name
    :param obsVar2: masked_array or list of masked_array
        masked_array of the second observed variable of the given metric
        if several observations are used for this variable, 'obsNameVar1' is a list of observations' masked_array

    :return: dict()
        dict_metrics = {'name', 'datasets', 'metric_values', 'method_to_compute_metric', 'raw_values', dict_diagnostic}
    """
    # retrieving metric collection dictionary: EnsoCollectionsLib.defCollection
    dict_mc = defCollection(metricCollection)
    # read the list of variables for the given metric
    list_variables = dict_mc['metrics_list'][metric]['variables']
    # read the dictionary of parameters
    dict_param = defCollection(metricCollection)['common_collection_parameters']

    # common_collection_parameters
    kwargs['metric_collection'] = metricCollection
    for arg in dict_param.keys():
        kwargs[arg] = dict_param[arg]
    # if 'metric_computation' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its
    # default value
    try:
        kwargs['metric_computation']
    except:
        kwargs['metric_computation'] = DefaultArgValues('metric_computation')

    # obsName could be a list if the user wants to compare the model with a set of observations
    # if obsName is just a name (string) it is put in a list
    if isinstance(obsNameVar1, basestring):
        obsNameVar1 = [obsNameVar1]
        obsVar1 = [obsVar1]
    try: kwargs['obsNameVar2']
    except:
        if isinstance(kwargs['obsNameVar2'], basestring):
            kwargs['obsNameVar2'] = [kwargs['obsNameVar2']]
            kwargs['obsVar2'] = [kwargs['obsVar2']]

    dict_metric_val = dict()
    if metric in dict_oneVar_modelAndObs.keys():
        #
        # this part regroups all diagnostics comparing model and obs (rmse)
        # so the diagnostic is the metric
        #
        description_metric = "The metric is the statistical value between the model and the observations"
        for ii in range(len(obsNameVar1)):
            # sets observations
            obs, tab_obs = obsNameVar1[ii], obsVar1[ii]
            # computes the diagnostic/metric
            diagnostic1 = dict_oneVar_modelAndObs[metric](modelVar1, modelName, obsVar1, obs, **kwargs)
            # puts metric values in its proper dictionary
            dict_metric_val['ref_' + obs] = {'value': diagnostic1['value'], 'value_error': diagnostic1['value_error']}
            # puts diagnostic values in its proper dictionary
            try:
                dict_diagnostic
            except:
                # first observations (first round in the obs loop)
                # as the diagnostic is the metric, only the information about the model and the observations are saved
                dict_diagnostic = {
                    'model': {
                        'name': modelName, 'value': None, 'value_error': None, 'nyears': diagnostic1['nyears_model'],
                    },
                    'observations': {
                        obs: {
                            'name': obs, 'value': None, 'value_error': None,
                            'nyears': diagnostic1['nyears_observations'],
                        },
                    },
                }
                if 'events_model' in diagnostic1.keys():
                    dict_diagnostic['model']['events'] = diagnostic1['events_model']
                if 'events_observations' in diagnostic1.keys():
                    dict_diagnostic['observations'][obs]['events'] = diagnostic1['events_observations']
            else:
                # next rounds in the obs loop
                # information about the model have already been saved: information about the observations are saved
                dict_diagnostic['observations'][obs] = {
                    'name': obs, 'value': None, 'value_error': None, 'nyears': diagnostic1['nyears_observations'],
                }
                if 'events_observations' in diagnostic1.keys():
                    dict_diagnostic['observations'][obs]['events'] = diagnostic1['events_observations']
    else:
        #
        # this part regroups all diagnostics that are computed separately for model and obs
        #
        # model diagnostic
        #
        if metric in dict_oneVar.keys():
            # computes diagnostic that needs only one variable
            diagnostic1 = dict_oneVar[metric](modelVar1, modelName, **kwargs)
        elif metric in dict_twoVar.keys():
            # computes diagnostic that needs two variables
            diagnostic1 = dict_twoVar[metric](modelVar1, modelName, kwargs['modelVar2'], **kwargs)
        else:
            diagnostic1 = None
            list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": metric", str().ljust(5) +
                            "unknown metric name: " + str(metric)]
            EnsoErrorsWarnings.MyError(list_strings)
        # saves the model diagnostic and the information about the model
        dict_diagnostic = {
            'model': {
                'name': modelName, 'value': diagnostic1['value'], 'value_error': diagnostic1['value_error'],
                'nyears': diagnostic1['nyears'],
            },
        }
        try:
            diagnostic1['nonlinearity']
        except:
            pass
        else:
            dict_diagnostic['model']['nonlinearity'] = diagnostic1['nonlinearity']
            dict_diagnostic['model']['nonlinearity_error'] = diagnostic1['nonlinearity_error']
        #
        # observations diagnostic
        #
        diag_obs = dict()
        for ii in range(len(obsNameVar1)):
            # sets observations
            obs1, tab_obs1 = obsNameVar1[ii], obsVar1[ii]
            if metric in dict_oneVar.keys():
                diag_obs[obs1] = dict_oneVar[metric](tab_obs1, obs1, **kwargs)
            elif metric in dict_twoVar.keys():
                for jj in range(len(kwargs['obsNameVar2'])):
                    # sets observations
                    obs2, tab_obs2 = kwargs['obsNameVar2'][ii], kwargs['obsVar2'][ii]
                    output_name = list_variables[0] + '_' + str(obs1) + '__' + list_variables[1] + '_' + str(obs2)
                    diag_obs[output_name] = dict_twoVar[metric](tab_obs1, obs1, tab_obs2, obs2, **kwargs)
        # saves the model diagnostic and the information about the model
        for obs in diag_obs.keys():
            try:
                dict_diagnostic['observations']
            except:
                # first observations (first round in the obs loop
                dict_diagnostic['observations'] = {
                    obs: {
                        'name': obs, 'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error'],
                        'nyears': diag_obs[obs]['nyears'],
                    },
                }
            else:
                # next rounds in the obs loop
                dict_diagnostic['observations'][obs] = {
                    'name': obs, 'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error'],
                    'nyears': diagnostic1['nyears'],
                }
            # computes the metric
            metric_val, metric_err, descript = MathMetriComputation(
                diagnostic1['value'], diagnostic1['value_error'], obs=diag_obs[obs]['value'],
                obs_err=diag_obs[obs]['value_error'], keyword=kwargs['metric_computation'])
            description_metric = descript
            dict_metric_val['ref_' + str(obs)] = {'value': metric_val, 'value_error': metric_err}
            try:
                diagnostic1['nonlinearity']
            except:
                pass
            else:
                # more computation if the diagnostic includes a nonlinearity measurement
                dict_diagnostic['observations'][obs]['nonlinearity'] = diag_obs[obs]['nonlinearity']
                dict_diagnostic['observations'][obs]['nonlinearity_error'] = diag_obs[obs]['nonlinearity_error']
                metric_val, metric_err, string = MathMetriComputation(
                    diagnostic1['nonlinearity'], diagnostic1['nonlinearity_error'], obs=diag_obs[obs]['nonlinearity'],
                    obs_err=diag_obs[obs]['nonlinearity_error'], keyword=kwargs['metric_computation'])
                dict_metric_val['ref_' + str(obs)]['nonlinearity'] = metric_val
                dict_metric_val['ref_' + str(obs)]['nonlinearity_error'] = metric_err
    # finishes to fill the diagnostic dictionary
    list_keys = ['name', 'units', 'time_frequency', 'method_to_compute_diagnostic', 'ref']
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
    # creates the output dictionary
    datasets = [modelName] + obsNameVar1
    if len(list_variables) > 1:
        datasets = datasets + kwargs['obsNameVar2']
    dict_metrics = {
        'name': metric, 'datasets': datasets, 'metric_values': dict_metric_val,
        'method_to_compute_metric': description_metric, 'raw_values': dict_diagnostic,
    }
    return dict_metrics
# ---------------------------------------------------------------------------------------------------------------------#
