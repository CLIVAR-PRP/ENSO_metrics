# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoMetricsLib import BiasPrLatRmse, BiasPrLonRmse, BiasPrRmse, BiasSstLonRmse, BiasSstLatRmse, BiasSstRmse, \
    BiasTauxLatRmse, BiasTauxLonRmse, BiasTauxRmse, EnsoAlphaLhf, EnsoAlphaLwr, EnsoAlphaShf, EnsoAlphaSwr, \
    EnsoAlphaThf, EnsoAmpl, EnsoMu, EnsoPrJjaTel, EnsoSeasonality, NinaSstLonRmse, NinaSstTsRmse, NinoSstLonRmse, \
    NinoSstTsRmse, SeasonalPrLatRmse, SeasonalPrLonRmse, SeasonalSstLatRmse, SeasonalSstLonRmse
from KeyArgLib import DefaultArgValues



# ---------------------------------------------------------------------------------------------------------------------#
#
# Computation of the metric collection
#
def ComputeCollection(metricCollection, dictDatasets, user_regridding={}, degug=False, dive_down=False):
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
    dict_col_meta = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_col_dd_meta = {
        'name': dict_mc['long_name'], 'description_of_the_collection': dict_mc['description'], 'metrics': {},
    }
    dict_col_valu = dict()
    dict_col_dd_valu = dict()
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
            arg_var2['modelFile2'] = dictDatasets['model'][modelName][list_variables[1]]['path + filename']
            arg_var2['modelVarName2'] = dictDatasets['model'][modelName][list_variables[1]]['varname']
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
        else:dict_col_valu[metric], dict_col_meta['metrics'][metric], dict_col_dd_valu[metric], \
            dict_col_dd_meta['metrics'][metric] = ComputeMetric(
                metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                dict_regions[list_variables[0]], user_regridding=user_regridding, degug=degug, **arg_var2)
    if dive_down is True:
        return {'value': dict_col_valu, 'metadata': dict_col_meta},\
               {'value': dict_col_dd_valu, 'metadata': dict_col_dd_meta}
    else:
        return {'value': dict_col_valu, 'metadata': dict_col_meta}
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

dict_twoVar_modelAndObs = {'EnsoPrJjaTel': EnsoPrJjaTel}

dict_oneVar = {'EnsoAmpl': EnsoAmpl, 'EnsoSeasonality': EnsoSeasonality}

dict_twoVar = {
    'EnsoAlphaLhf': EnsoAlphaLhf, 'EnsoAlphaLwr': EnsoAlphaLwr, 'EnsoAlphaShf': EnsoAlphaShf,
    'EnsoAlphaSwr': EnsoAlphaSwr, 'EnsoAlphaThf': EnsoAlphaThf, 'EnsoMu': EnsoMu,
}


def ComputeMetric(metricCollection, metric, modelName, modelFile1, modelVarName1, obsNameVar1, obsFile1, obsVarName1,
                  regionVar1, modelFileArea1='', modelAreaName1='', modelFileLandmask1='', modelLandmaskName1='',
                  obsFileArea1='', obsAreaName1='', obsFileLandmask1='', obsLandmaskName1='',
                  modelFile2='', modelVarName2='', modelFileArea2='', modelAreaName2='', modelFileLandmask2='',
                  modelLandmaskName2='', obsNameVar2='', obsFile2='', obsVarName2='', obsFileArea2='', obsAreaName2='',
                  obsFileLandmask2='', obsLandmaskName2='', regionVar2='', degug=False, user_regridding={}):
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

    if metric in dict_oneVar_modelAndObs.keys():
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
                print '\033[94m' + str().ljust(5) + "ComputeMetric: RMSmetric, " + metric + " = " + modelName + " and "\
                      + output_name + '\033[0m'
                diagnostic1[output_name] = dict_oneVar_modelAndObs[metric](
                    modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                    obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                    obsLandmaskName1[ii], box=regionVar1, degug=degug, **keyarg)
            elif metric in dict_twoVar_modelAndObs.keys():
                for jj in range(len(obsNameVar2)):
                    output_name = obsNameVar1[ii] + '_' + obsNameVar2[jj]
                    print '\033[94m' + str().ljust(5) + "ComputeMetric: RMSmetric, " + metric + " = " + modelName\
                          + " and " + output_name + '\033[0m'
                    diagnostic1[output_name] = dict_twoVar_modelAndObs[metric](
                        modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1,
                        modelLandmaskName1, modelFile2, modelVarName2, modelFileArea2, modelAreaName2,
                        modelFileLandmask2, modelLandmaskName2, obsFile1[ii], obsVarName1[ii], obsFileArea1[ii],
                        obsAreaName1[ii], obsFileLandmask1[ii], obsLandmaskName1[ii], obsFile2[jj], obsVarName2[jj],
                        obsFileArea2[jj], obsAreaName2[jj], obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar1,
                        regionVar2, degug=degug, **keyarg)
        for obs in diagnostic1.keys():
            # puts metric values in its proper dictionary
            dict_metric_val['ref_' + obs] = {
                'value': diagnostic1[obs]['value'], 'value_error': diagnostic1[obs]['value_error']}
            if 'value2' in diagnostic1.keys():
                dict_metric_val['ref_' + obs]['value2'] = diagnostic1[obs]['value2']
                dict_metric_val['ref_' + obs]['value_error2'] = diagnostic1[obs]['value_error2']
            dict_diagnostic['model'] = {'value': None, 'value_error': None}
            dict_diagnostic[obs] = {'value': None, 'value_error': None}
            if 'dive_down_diag' in diagnostic1[obs].keys():
                dict_dive_down['model'] = diagnostic1['dive_down_diag']['model']
                dict_dive_down[obs] = diagnostic1['dive_down_diag']['observations']
                dict1 = {}
                for elt in diagnostic1[obs]['dive_down_diag'].keys():
                    if elt not in ['model', 'observations']:
                        dict1[elt] = diagnostic1[obs]['dive_down_diag'][elt]
                dict_dive_down_metadata[obs] = dict1
                del dict1
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata['model'] = {'name': modelName, 'nyears': diagnostic1[obs]['nyears_model'],
                                                 'time_period': diagnostic1[obs]['time_period_model']},
            dict_diagnostic_metadata[obs] = {'name': obs, 'nyears': diagnostic1[obs]['nyears_model'],
                                             'time_period': diagnostic1[obs]['time_period_observations']},
            if 'events_model' in diagnostic1[obs].keys():
                dict_diagnostic_metadata['model']['events'] = diagnostic1[obs]['events_model']
                dict_diagnostic_metadata[obs]['events'] = diagnostic1[obs]['events_observations']
        del diagnostic1
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
            diagnostic1 = dict_oneVar[metric](modelFile1, modelVarName1, modelFileArea1, modelAreaName1,
                                              modelFileLandmask1, modelLandmaskName1, regionVar1, degug=degug, **keyarg)
        elif metric in dict_twoVar.keys():
            # computes diagnostic that needs two variables
            print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(modelName) + '\033[0m'
            keyarg['project_interpreter_var2'] = keyarg['project_interpreter']
            diagnostic1 = dict_twoVar[metric](
                modelFile1, modelVarName1, modelFileArea1, modelAreaName1, modelFileLandmask1, modelLandmaskName1,
                modelFile2, modelVarName2, modelFileArea2, modelAreaName2, modelFileLandmask2, modelLandmaskName2,
                regionVar1, regionVar2, degug=degug, **keyarg)
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
            obs1 = obsNameVar1[ii]
            keyarg['project_interpreter_var1'] = obs1
            #            keyarg['project_interpreter'] = 'CMIP'
            if metric in dict_oneVar.keys():
                print '\033[94m' + str().ljust(5) + "ComputeMetric: oneVarmetric = " + str(obs1) + '\033[0m'
                output_name = obs1
                diag_obs[output_name] = dict_oneVar[metric](
                    obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                    obsLandmaskName1[ii], regionVar1, degug=degug, **keyarg)
            elif metric in dict_twoVar.keys():
                for jj in range(len(obsNameVar2)):
                    obs2 = obsNameVar2[jj]
                    output_name = obs1 + '_' + obs2
                    keyarg['project_interpreter_var2'] = obs2
                    print '\033[94m' + str().ljust(5) + "ComputeMetric: twoVarmetric = " + str(output_name) + '\033[0m'
                    diag_obs[output_name] = dict_twoVar[metric](
                        obsFile1[ii], obsVarName1[ii], obsFileArea1[ii], obsAreaName1[ii], obsFileLandmask1[ii],
                        obsLandmaskName1[ii], obsFile2[jj], obsVarName2[jj], obsFileArea2[jj], obsAreaName2[jj],
                        obsFileLandmask2[jj], obsLandmaskName2[jj], regionVar1, regionVar2, **keyarg)
        for obs in diag_obs.keys():
            # computes the metric
            metric_val, metric_err, description_metric = MathMetriComputation(
                diagnostic1['value'], diagnostic1['value_error'], obs=diag_obs[obs]['value'],
                obs_err=diag_obs[obs]['value_error'], keyword=keyarg['metric_computation'])
            dict_metric_val['ref_' + obs] = {'value': metric_val, 'value_error': metric_err}
            # puts metric / diagnostic values in its proper dictionary
            dict_diagnostic[obs] = {'value': diag_obs[obs]['value'], 'value_error': diag_obs[obs]['value_error']}
            if 'nonlinearity' in diag_obs[obs].keys():
                dict_diagnostic[obs]['nonlinearity'] = diag_obs[obs]['nonlinearity']
                dict_diagnostic[obs]['nonlinearity_error'] = diag_obs[obs]['nonlinearity_error']
            if 'dive_down_diag' in diag_obs[obs].keys():
                dict_dive_down[obs] = diag_obs[obs]['dive_down_diag']['value']
            # puts diagnostic metadata in its proper dictionary
            dict_diagnostic_metadata[obs] = {
                'name': obs, 'nyears': diag_obs[obs]['nyears'], 'time_period': diag_obs[obs]['time_period'],
            }
            if 'events_model' in diag_obs[obs].keys():
                dict_diagnostic_metadata[obs]['events'] = diag_obs[obs]['events']
        if keyarg['metric_computation'] in ['ratio', 'relative_difference']:
            units = ''
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
    datasets = modelName + '; '
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
        return dict_metrics, dict_metadata, dict_dive_down, dict_dive_down_metadata
    else:
        return dict_metrics, dict_metadata
# ---------------------------------------------------------------------------------------------------------------------#
