# -*- coding:UTF-8 -*-
# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
from EnsoToolsLib import ComputeMetric

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
