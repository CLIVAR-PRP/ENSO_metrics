# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
import numpy

# ENSO_metrics package functions:
import EnsoErrorsWarnings
from EnsoToolsLib import __rms


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO metrics
# These functions have file names and variable names as inputs and metric as output
#
def EnsoRmse(modelVar1, modelNameVar1, obsVar1, obsNameVar1, centered_rmse=0, **kwargs):
    """
    The EnsoRmse() function computes the SST spatial root mean square error (RMSE) between model and observations
        Modeled and observed data must be on the same grid

    Inputs:
    ------
    :param modelVar1: masked_array
        masked_array of the modeled SST (already preprocessed)
    :param modelNameVar1: string
        name of the model for the first variable
    :param obsVar1: masked_array
        masked_array of the observed SST (already preprocessed)
    :param obsNameVar1: string
        name of the observations for the first variable
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    usual kwargs:
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_number_of_years: int, optional
        minimum number of years for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None

    Output:
    ------
    :return rmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, ref

    Method:
    -------
        uses tools from numpy library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    metric_name = "EnsoRmse"
    # get
    mod_tab, mod_axes, mod_weights, mod_Nyears = modelVar1['array'], modelVar1['axes'], modelVar1['weights'], \
                                                 modelVar1['Nyears']
    obs_tab, obs_axes, obs_weights, obs_Nyears = obsVar1['array'], obsVar1['axes'], obsVar1['weights'], \
                                                 obsVar1['Nyears']
    region = kwargs['regionVar1']
    # Define metric attributes
    Name = 'ENSO RMSE'
    Units = 'C'
    Method = 'Spatial root mean square error of ' + region + ' sst'
    if centered_rmse == 1:
        rms_centered = 'centered'
    else:
        rms_centered = 'uncentered'
    Ref = 'Using numpy (v' + str(numpy.__version__) + ') rms (' + str(rms_centered) + ' and biased) calculation'

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_number_of_years'], int):
        mini = kwargs['min_number_of_years']
        if mod_Nyears < mini:
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, mod_Nyears, mini, INSPECTstack(), info=modelNameVar1)
        if obs_Nyears < mini:
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, obs_Nyears, mini, INSPECTstack(), info=obsNameVar1)

    # Computes the root mean square difference
    sstRmse = __rms(mod_tab.flatten(), obs_tab.flatten(), weights=mod_weights.flatten(),
                    centered=centered_rmse, biased=1)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': mod_Nyears, 'nyears_observations': obs_Nyears, 'time_frequency': kwargs['frequency'],
        'ref': Ref,
    }
    return rmseMetric


def EnsoLatRmse(modelVar1, modelNameVar1, obsVar1, obsNameVar1, centered_rmse=0, **kwargs):
    """
    The EnsoLatRmse() function computes the SST meridional (latitude) root mean square error (RMSE) between model and
        observations
        Modeled and observed data must be on the same grid

    Inputs:
    ------
    :param modelVar1: masked_array
        masked_array of the modeled SST (already preprocessed)
    :param modelNameVar1: string
        name of the model for the first variable
    :param obsVar1: masked_array
        masked_array of the observed SST (already preprocessed)
    :param obsNameVar1: string
        name of the observations for the first variable
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    usual kwargs:
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_number_of_years: int, optional
        minimum number of years for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None

    Output:
    ------
    :return LatRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, ref

    Method:
    -------
        uses tools from numpy library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    metric_name = "EnsoLatRmse"
    # get
    mod_tab, mod_axes, mod_weights, mod_Nyears = modelVar1['array'], modelVar1['axes'], modelVar1['weights'], \
                                                 modelVar1['Nyears']
    obs_tab, obs_axes, obs_weights, obs_Nyears = obsVar1['array'], obsVar1['axes'], obsVar1['weights'], \
                                                 obsVar1['Nyears']
    region = kwargs['regionVar1']
    # Define metric attributes
    Name = 'ENSO Meridional RMSE'
    Units = 'C'
    Method = 'Meridional root mean square error of ' + region + ' sst'
    if centered_rmse == 1:
        rms_centered = 'centered'
    else:
        rms_centered = 'uncentered'
    Ref = 'Using numpy (v' + str(numpy.__version__) + ') rms (' + str(rms_centered) + ' and biased) calculation'

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_number_of_years'], int):
        mini = kwargs['min_number_of_years']
        if mod_Nyears < mini:
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, mod_Nyears, mini, INSPECTstack(), info=modelNameVar1)
        if obs_Nyears < mini:
            EnsoErrorsWarnings.TooShortTimePeriod(metric_name, obs_Nyears, mini, INSPECTstack(), info=obsNameVar1)

    # Computes the root mean square difference
    sstRmse = __rms(mod_tab.flatten(), obs_tab.flatten(), weights=mod_weights.flatten(),
                    centered=centered_rmse, biased=1)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': mod_Nyears, 'nyears_observations': obs_Nyears, 'time_frequency': kwargs['frequency'],
        'ref': Ref,
    }
    return rmseMetric
# ---------------------------------------------------------------------------------------------------------------------#

