# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoMetricsLib import BiasPrLatRmse, BiasPrLonRmse, BiasPrRmse, BiasSstLonRmse, BiasSstLatRmse, BiasSstRmse, \
    BiasTauxLatRmse, BiasTauxLonRmse, BiasTauxRmse, EnsoAlphaLhf, EnsoAlphaLwr, EnsoAlphaShf, EnsoAlphaSwr, \
    EnsoAlphaThf, EnsoAmpl, EnsoMu, EnsoSeasonality, NinaSstLonRmse, NinaSstTsRmse, NinoSstLonRmse, NinoSstTsRmse, \
    SeasonalPrLatRmse, SeasonalPrLonRmse, SeasonalSstLatRmse, SeasonalSstLonRmse
from KeyArgLib import DefaultArgValues



# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric collection
#
def ComputeCollection(metricCollection, dictDatasets, user_regridding={}):
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
    dict_collection_metadata = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_collection_value = dict()
    dict_m = dict_mc['metrics_list']
    list_metrics = sorted(dict_m.keys())
    for metric in list_metrics:
        print '\033[94m' + str().ljust(5) + "ComputeCollection: metric = " + str(metric) + '\033[0m'
        # sets arguments for this metric
        list_variables = dict_m[metric]['variables']
        dict_regions = dict_m[metric]['regions']
        # model name, file, variable name in file
        modelName = dictDatasets['model'].keys()[0]
        modelFile1 = dictDatasets['model'][modelName][list_variables[0]]['path + filename']
        modelVarName1 = dictDatasets['model'][modelName][list_variables[0]]['varname']
        # observations name(s), file(s), variable(s) name in file(s)
        obsNameVar1, obsFile1, obsVarName1 = list(), list(), list()
        for obs in dictDatasets['observations'].keys():
            try: dictDatasets['observations'][obs][list_variables[0]]
            except: pass
            else:
                obsNameVar1.append(obs)
                obsFile1.append(dictDatasets['observations'][obs][list_variables[0]]['path + filename'])
                obsVarName1.append(dictDatasets['observations'][obs][list_variables[0]]['varname'])
        # same if a second variable is needed
        # this time in the form of a keyarg dictionary
        arg_var2 = {'obsFile2': []}
        if len(list_variables) > 1:
            arg_var2['modelFile2'] = dictDatasets['model'][modelName][list_variables[1]]['path + filename']
            arg_var2['modelVarName2'] = dictDatasets['model'][modelName][list_variables[1]]['varname']
            arg_var2['regionVar2'] = dict_regions[list_variables[1]]
            obsNameVar2, obsFile2, obsVarName2 = list(), list(), list()
            for obs in dictDatasets['observations'].keys():
                try: dictDatasets['observations'][obs][list_variables[1]]
                except: pass
                else:
                    obsNameVar2.append(obs)
                    obsFile2.append(dictDatasets['observations'][obs][list_variables[1]]['path + filename'])
                    obsVarName2.append(dictDatasets['observations'][obs][list_variables[1]]['varname'])
            arg_var2['obsNameVar2'] = obsNameVar2
            arg_var2['obsFile2'] = obsFile2
            arg_var2['obsVarName2'] = obsVarName2
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
            dict_collection_value[metric], dict_collection_metadata['metrics'][metric] = ComputeMetric(
                metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                dict_regions[list_variables[0]], user_regridding=user_regridding, **arg_var2)

    return {'value': dict_collection_value, 'metadata': dict_collection_metadata}
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
dict_oneVar_modelAndObs = {
    'BiasPrLatRmse': BiasPrLatRmse, 'BiasPrLonRmse': BiasPrLonRmse, 'BiasPrRmse': BiasPrRmse,
    'BiasSstLatRmse': BiasSstLatRmse, 'BiasSstLonRmse': BiasSstLonRmse, 'BiasSstRmse': BiasSstRmse,
    'BiasTauxLatRmse': BiasTauxLatRmse, 'BiasTauxLonRmse': BiasTauxLonRmse, 'BiasTauxRmse': BiasTauxRmse,
    'NinaSstTsRmse': NinaSstTsRmse, 'NinaSstLonRmse': NinaSstLonRmse,
    'NinoSstTsRmse': NinoSstTsRmse, 'NinoSstLonRmse': NinoSstLonRmse,
    'SeasonalPrLatRmse': SeasonalPrLatRmse, 'SeasonalPrLonRmse': SeasonalPrLonRmse,
    'SeasonalSstLatRmse': SeasonalSstLatRmse, 'SeasonalSstLonRmse': SeasonalSstLonRmse,
}

dict_oneVar = {'EnsoAmpl': EnsoAmpl, 'EnsoSeasonality': EnsoSeasonality,}

dict_twoVar = {
    'EnsoAlphaLhf': EnsoAlphaLhf, 'EnsoAlphaLwr': EnsoAlphaLwr, 'EnsoAlphaShf': EnsoAlphaShf,
    'EnsoAlphaSwr': EnsoAlphaSwr, 'EnsoAlphaThf': EnsoAlphaThf, 'EnsoMu': EnsoMu,
}


def ComputeMetric(metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                  regionVar1, modelFile2='', modelVarName2='', obsNameVar2='', obsFile2='', obsVarName2='',
                  regionVar2='', user_regridding={}):
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
    :param user_regridding:
    :return:
    """
    # retrieving keyargs from EnsoCollectionsLib.defCollection
    dict_mc = defCollection(metricCollection)
    # read the list of variables for the given metric
    list_variables = dict_mc['metrics_list'][metric]['variables']

    # common_collection_parameters
    keyarg = dict()
    for arg in dict_mc['common_collection_parameters'].keys():
        keyarg[arg] = dict_mc['common_collection_parameters'][arg]
    for arg in dict_mc['metrics_list'][metric].keys():
        keyarg[arg] = dict_mc['metrics_list'][metric][arg]
    # if 'metric_computation' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its
    # default value
    try:
        keyarg['metric_computation']
    except:
        keyarg['metric_computation'] = DefaultArgValues('metric_computation')
    # if 'modeled_period' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg['time_bounds_model'] = keyarg['modeled_period']
    except:
        keyarg['time_bounds_model'] = DefaultArgValues('time_bounds_model')
    # if 'modeled_period' is not defined for this metric (in EnsoCollectionsLib.defCollection), sets it to its default
    # value
    try:
        keyarg['time_bounds_obs'] = keyarg['observed_period']
    except:
        keyarg['time_bounds_obs'] = DefaultArgValues('time_bounds_obs')
    # if the user gave a specific regridding Tool / method, use it
    if metric in user_regridding.keys():
        keyarg['regridding'] = user_regridding[metric]
    elif 'regridding'in user_regridding.keys():
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
    dict_diagnostic = dict()#{'model': {}}
    dict_diagnostic_metadata = dict()#{'model': {}}
    dict_dive_down = dict()#{'model': {}}
    dict_dive_down_metadata = dict()
    # for obs1 in obsNameVar1:
    #     if len(obsNameVar2) > 0:
    #         for obs2 in obsNameVar2:
    #             dict_diagnostic[obs1 + '_' + obs2] = {}
    #             dict_diagnostic_metadata[obs1 + '_' + obs2] = {}
    #             dict_dive_down[obs1 + '_' + obs2] = {}
    #     else:
    #         dict_diagnostic[obs1] = {}
    #         dict_diagnostic_metadata[obs1 + '_' + obs2] = {}
    #         dict_dive_down[obs1 + '_' + obs2] = {}

    if metric in dict_oneVar_modelAndObs.keys():
        #
        # this part regroups all diagnostics comparing model and obs (rmse)
        # so the diagnostic is the metric
        #
        description_metric = "The metric is the statistical value between the model and the observations"
        for ii in range(len(obsNameVar1)):
            # sets observations
            obs, tmp_obsFile1, tmp_obsVarName1 = obsNameVar1[ii], obsFile1[ii], obsVarName1[ii]
            print '\033[94m' + str().ljust(5) + "ComputeMetric: RMSmetric = " + str(modelName) + " and " + str(obs)\
                  + '\033[0m'
            # computes the diagnostic/metric
            diagnostic1 = dict_oneVar_modelAndObs[metric](
                modelFile1, modelVarName1, tmp_obsFile1, tmp_obsVarName1, box=regionVar1, **keyarg)
            # puts metric / diagnostic values in its proper dictionary
            dict_metric_val[obs] = {'value': diagnostic1['value'], 'value_error': diagnostic1['value_error']}
            dict_diagnostic['model'] = {'value': None, 'value_error': None}
            dict_diagnostic[obs] = {'value': None, 'value_error': None}
            if 'dive_down_diag' in diagnostic1.keys():
                dict_dive_down['model'] = diagnostic1['dive_down_diag']['model']
                dict_dive_down[obs] = diagnostic1['dive_down_diag']['observations']
                for elt in diagnostic1['dive_down_diag'].keys():
                    if elt not in ['model', 'observations']:
                        dict_dive_down_metadata[elt] = diagnostic1['dive_down_diag'][elt]
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata['model'] = {
                'name': modelName, 'nyears': diagnostic1['nyears_observations'],
                'time_period': diagnostic1['time_period_model'],
            }
            dict_diagnostic_metadata[obs] = {
                'name': obs, 'nyears': diagnostic1['nyears_model'],
                'time_period': diagnostic1['time_period_observations'],
            }
            if 'events_model' in diagnostic1.keys():
                dict_diagnostic_metadata['model']['events'] = diagnostic1['events_model']
                dict_diagnostic_metadata[obs]['events'] = diagnostic1['events_observations']
        units = diagnostic1['units']
    else:
        #
        # model diagnostic
        #
        keyarg['time_bounds'] = keyarg['time_bounds_model']
        keyarg['project_interpreter_var1'] = keyarg['project_interpreter']
        if metric in dict_oneVar.keys():
            # computes diagnostic that needs only one variable
            print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(modelName) + '\033[0m'
            diagnostic1 = dict_oneVar[metric](modelFile1, modelVarName1, regionVar1, **keyarg)
        elif metric in dict_twoVar.keys():
            # computes diagnostic that needs two variables
            print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(modelName) + '\033[0m'
            keyarg['project_interpreter_var2'] = keyarg['project_interpreter']
            diagnostic1 = dict_twoVar[metric](modelFile1, modelFile2, modelVarName1, modelVarName2, regionVar1,
                                              regionVar2, **keyarg)
        else:
            diagnostic1 = None
            list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": metric", str().ljust(5) +
                            "unknown metric name: " + str(metric)]
            EnsoErrorsWarnings.MyError(list_strings)
        # puts metric / diagnostic values in its proper dictionary
        dict_diagnostic['model'] = {'value': diagnostic1['value'], 'value_error': diagnostic1['value_error']}
        if 'nonlinearity' in diagnostic1.keys():
            dict_diagnostic['model']['nonlinearity'] = diagnostic1['nonlinearity']
            dict_diagnostic['model']['nonlinearity_error'] = diagnostic1['nonlinearity_error']
        if 'dive_down_diag' in diagnostic1.keys():
            dict_dive_down['model'] = diagnostic1['dive_down_diag']['value']
            for elt in diagnostic1['dive_down_diag'].keys():
                if elt not in ['value']:
                    dict_dive_down_metadata[elt] = diagnostic1['dive_down_diag'][elt]
        # puts diagnostic metadata in its proper dictionary
        dict_diagnostic_metadata['model'] = {
            'name': modelName, 'nyears': diagnostic1['nyears'], 'time_period': diagnostic1['time_period'],
        }
        if 'events_model' in diagnostic1.keys():
            dict_diagnostic_metadata['model']['events'] = diagnostic1['events']
        #
        # observations diag
        #
        diag_obs = dict()
        keyarg['time_bounds'] = keyarg['time_bounds_obs']
        for ii in range(len(obsNameVar1)):
            # sets observations
            obs1, tmp_obsFile1, tmp_obsVarName1 = obsNameVar1[ii], obsFile1[ii], obsVarName1[ii]
            keyarg['project_interpreter_var1'] = obs1
#            keyarg['project_interpreter'] = 'CMIP'
            if metric in dict_oneVar.keys():
                print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(obs1) + '\033[0m'
                output_name = obs1
                diag_obs[output_name] = dict_oneVar[metric](tmp_obsFile1, tmp_obsVarName1, regionVar1, **keyarg)
            elif metric in dict_twoVar.keys():
                for jj in range(len(obsNameVar2)):
                    obs2, tmp_obsFile2, tmp_obsVarName2 = obsNameVar2[jj], obsFile2[jj], obsVarName2[jj]
                    output_name = obs1 + '_' + obs2
                    keyarg['project_interpreter_var2'] = obs2
                    print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(output_name) + '\033[0m'
                    diag_obs[output_name] = dict_twoVar[metric](tmp_obsFile1, tmp_obsFile2, tmp_obsVarName1,
                                                                tmp_obsVarName2, regionVar1, regionVar2, **keyarg)
        for obs in diag_obs.keys():
            # puts metric / diagnostic values in its proper dictionary
            dict_diagnostic[obs] = {'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error']}
            if 'nonlinearity' in diag_obs[obs].keys():
                dict_diagnostic[obs]['nonlinearity'] = diag_obs[obs]['nonlinearity']
                dict_diagnostic[obs]['nonlinearity_error'] = diag_obs[obs]['nonlinearity_error']
            if 'dive_down_diag' in diag_obs[obs].keys():
                dict_dive_down[obs] = diag_obs[obs]['dive_down_diag']['value']
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata[obs] = {
                'name': modelName, 'nyears': diag_obs[obs]['nyears'], 'time_period': diag_obs[obs]['time_period'],
            }
            if 'events_model' in diag_obs[obs].keys():
                dict_diagnostic_metadata['model']['events'] = diag_obs[obs]['events']
        if keyarg['metric_computation'] in ['ratio', 'relative_difference']:
            units = ''
        else:
            units = diagnostic1['units']
    # finishes to fill the diagnostic dictionary
    list_keys = ['method', 'name', 'ref', 'time_frequency', 'units']
    for key in list_keys:
        if key == 'method':
            dict_diagnostic_metadata[key] = diagnostic1['method']
            try:
                diagnostic1['nonlinearity']
            except:
                pass
            else:
                dict_diagnostic_metadata['method_nonlinearity'] = diagnostic1['method_nonlinearity']
        else:
            dict_diagnostic_metadata[key] = diagnostic1[key]
    # creates the output dictionaries
    dict_metrics = {
        'metric': dict_metric_val, 'diagnostic': dict_diagnostic,
    }
    datasets = modelName+'; '
    for ii in range(len(obsNameVar1)):
        datasets = datasets + obsNameVar1[ii] + obsVarName1[ii] + "'s"
        if ii != len(obsNameVar1) - 1:
            datasets = datasets + ' & '
        else:
            datasets = datasets + '; '
    if len(obsNameVar2) > 0:
        for ii in range(len(obsNameVar2)):
            datasets = datasets + obsNameVar2[ii] + obsVarName2[ii] + "'s"
            if ii != len(obsNameVar2) - 1:
                datasets = datasets + ' & '
    dict_metadata = {
        'metric': {
            'name': metric, 'method': description_metric, 'datasets': datasets, 'units': units,
        },
        'diagnostic': dict_diagnostic_metadata,
    }
    if 'dive_down_diag' in diagnostic1.keys():
        dict_metrics['dive_down_diagnostic'] = dict_dive_down
        dict_metadata['dive_down_diagnostic'] = dict_dive_down_metadata
    return dict_metrics, dict_metadata
# ---------------------------------------------------------------------------------------------------------------------#
