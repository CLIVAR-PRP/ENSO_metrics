# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import sqrt as NUMPYsqrt
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
import EnsoErrorsWarnings
from EnsoUvcdatToolsLib import CheckTime, DetectEvents, Detrend, LinearRegressionAndNonlinearity, MyDerive,\
    OperationMultiply, PrePocessTS, ReadSelectRegionCheckUnits, SeasonalMean, SpatialRms, Std, TimeBounds, TwoVarRegrid
from KeyArgLib import DefaultArgValues


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO metrics
# These functions have file names and variable names as inputs and metric as output
#
def EnsoAlphaLhf(sstfile, lhffile, sstname, lhfname, sstbox, lhfbox, **kwargs):
    """
    The EnsoAlphaLhf() function computes the regression of 'lhfbox' lhfA (latent heat flux anomalies) over 'sstbox' sstA
    (usually the regression of nino3 lhfA over nino3 sstA)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param lhffile: string
        path_to/filename of the file (NetCDF) of the LHF dataset
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param lhfname: string
        name of LHF variable (lhf, hfls) in 'lhffile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param lhfbox: string
        name of box ('nino3') for LHF
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return alphaLhfMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Latent feedback (alpha_lh)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + lhfbox + ' lhfA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    lhf = ReadSelectRegionCheckUnits(lhffile, lhfname, 'heat flux', box=lhfbox, **kwargs)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lhf = CheckTime(sst, lhf, metric_name='EnsoAlphaLhf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    lhf, unneeded = PrePocessTS(lhf, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    alphaLhf, alphaLhfPos, alphaLhfNeg = LinearRegressionAndNonlinearity(lhf, sst, return_stderr=True)

    # Create output
    alphaLhfMetric = {
        'name': Name, 'value': alphaLhf[0], 'value_error': alphaLhf[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaLhfNeg[0] - alphaLhfPos[0],
        'nonlinearity_error': alphaLhfNeg[1] + alphaLhfPos[1],
    }
    return alphaLhfMetric


def EnsoAlphaLwr(sstfile, lwrfile, sstname, lwrname, sstbox, lwrbox, **kwargs):
    """
    The EnsoAlphaLwr() function computes the regression of 'lwrbox' lwrA (net surface longwave radiation anomalies) over
    'sstbox' sstA (usually the regression of nino3 lwrA over nino3 sstA)

    The net surface longwave radiation is not a CMIP variable.
    Either the user computes it and sends the filename and the varname or he feeds into lwrfile and lwrname of this
    function a list() of downward and upward radiations files and variable names (CMIP: rlds and rlus)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param lwrfile: string
        path_to/filename of the file (NetCDF) of the LWR dataset (may be a list of files)
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param lwrname: string
        name of LWR variable (lwr, rlds - rlus) (may be a list of variables) in 'lwrfile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param lwrbox: string
        name of box ('nino3') for LWR
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return alphaLwrMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Longwave feedback (alpha_lwr)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + lwrbox + ' lwrA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(lwrfile, basestring):
        dict_var[lwrname] = ReadSelectRegionCheckUnits(lwrfile, lwrname, 'heat flux', box=lwrbox, **kwargs)
    else:
        for ii in range(len(lwrfile)):
            filename, varname = lwrfile[ii], lwrname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=lwrbox, **kwargs)
    lwr = MyDerive(kwargs['project_interpreter'], 'lwr', dict_var)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lwr = CheckTime(sst, lwr, metric_name='EnsoAlphaLwr', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    lwr, unneeded = PrePocessTS(lwr, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    alphaLwr, alphaLwrPos, alphaLwrNeg = LinearRegressionAndNonlinearity(lwr, sst, return_stderr=True)

    # Create output
    alphaLwrMetric = {
        'name': Name, 'value': alphaLwr[0], 'value_error': alphaLwr[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaLwrNeg[0] - alphaLwrPos[0],
        'nonlinearity_error': alphaLwrNeg[1] + alphaLwrPos[1],
    }
    return alphaLwrMetric


def EnsoAlphaShf(sstfile, shffile, sstname, shfname, sstbox, shfbox, **kwargs):
    """
    The EnsoAlphaShf() function computes the regression of 'shfbox' shfA (sensible heat flux anomalies) over 'sstbox'
    sstA (usually the regression of nino3 shfA over nino3 sstA)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param shffile: string
        path_to/filename of the file (NetCDF) of the SHF dataset
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param shfname: string
        name of SHF variable (shf, hfss) in 'shffile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param shfbox: string
        name of box ('nino3') for SHF
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return alphaShfMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Latent feedback (alpha_lh)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + shfbox + ' shfA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    shf = ReadSelectRegionCheckUnits(shffile, shfname, 'heat flux', box=shfbox, **kwargs)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, shf = CheckTime(sst, shf, metric_name='EnsoAlphaShf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    shf, unneeded = PrePocessTS(shf, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    alphaShf, alphaShfPos, alphaShfNeg = LinearRegressionAndNonlinearity(shf, sst, return_stderr=True)

    # Create output
    alphaShfMetric = {
        'name': Name, 'value': alphaShf[0], 'value_error': alphaShf[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaShfNeg[0] - alphaShfPos[0],
        'nonlinearity_error': alphaShfNeg[1] + alphaShfPos[1],
    }
    return alphaShfMetric


def EnsoAlphaSwr(sstfile, swrfile, sstname, swrname, sstbox, swrbox, **kwargs):
    """
    The EnsoAlphaSwr() function computes the regression of 'swrbox' swrA (net surface shortwave radiation anomalies)
    over 'sstbox' sstA (usually the regression of nino3 swrA over nino3 sstA)

    The net surface shortwave radiation is not a CMIP variable.
    Either the user computes it and sends the filename and the varname or he feeds into swrfile and swrname of this
    function a list() of downward and upward radiations files and variable names (CMIP: rsds and rsus)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param swrfile: string
        path_to/filename of the file (NetCDF) of the SWR dataset (may be a list of files)
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param swrname: string
        name of SWR variable (swr, rsds - rsus) (may be a list of variables) in 'swrfile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param swrbox: string
        name of box ('nino3') for SWR
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return alphaSwrMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Longwave feedback (alpha_swr)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + swrbox + ' swrA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(swrfile, basestring):
        dict_var[swrname] = ReadSelectRegionCheckUnits(swrfile, swrname, 'heat flux', box=swrbox, **kwargs)
    else:
        for ii in range(len(swrfile)):
            filename, varname = swrfile[ii], swrname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=swrbox, **kwargs)
    swr = MyDerive(kwargs['project_interpreter'], 'swr', dict_var)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, swr = CheckTime(sst, swr, metric_name='EnsoAlphaSwr', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    swr, unneeded = PrePocessTS(swr, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    alphaSwr, alphaSwrPos, alphaSwrNeg = LinearRegressionAndNonlinearity(swr, sst, return_stderr=True)

    # Create output
    alphaSwrMetric = {
        'name': Name, 'value': alphaSwr[0], 'value_error': alphaSwr[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaSwrNeg[0] - alphaSwrPos[0],
        'nonlinearity_error': alphaSwrNeg[1] + alphaSwrPos[1],
    }
    return alphaSwrMetric


def EnsoAlphaThf(sstfile, thffile, sstname, thfname, sstbox, thfbox, **kwargs):
    """
    The EnsoAlphaThf() function computes the regression of 'thfbox' thfA (total heat flux anomalies) over 'sstbox' sstA
    (usually the regression of nino3 thfA over nino3 sstA)
    The total heat flux is the sum of four term:
         - net surface shortwave radiation,
         - net surface longwave radiation,
         - latent heat flux,
         - sensible heat flux

    The total heat flux is not always available is models or observations.
    Either the user computes it and sends the filename and the varname or he feeds into thffile and thfname of this
    function a list() of the four needed files and variable names (CMIP: rsds-rsus, rlds-rlus, hfls, hfss)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param thffile: string
        path_to/filename of the file (NetCDF) of the THF dataset (may be a list of files)
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param thfname: string
        name of THF variable (thf, netflux, thflx, swr + lwr + lhf + shf) (may be a list of variables) in 'thffile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param thfbox: string
        name of box ('nino3') for THF
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return alphaMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Heat flux feedback (alpha)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + thfbox + ' thfA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(thffile, basestring):
        dict_var[thfname] = ReadSelectRegionCheckUnits(thffile, thfname, 'heat flux', box=thfbox, **kwargs)
    else:
        for ii in range(len(thffile)):
            filename, varname = thffile[ii], thfname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=thfbox, **kwargs)
    thf = MyDerive(kwargs['project_interpreter'], 'thf', dict_var)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, thf = CheckTime(sst, thf, metric_name='EnsoAlphaThf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    thf, unneeded = PrePocessTS(thf, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    alpha, alphaPos, alphaNeg = LinearRegressionAndNonlinearity(thf, sst, return_stderr=True)

    # Create output
    alphaMetric = {
        'name': Name, 'value': alpha[0], 'value_error': alpha[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaNeg[0] - alphaPos[0],
        'nonlinearity_error': alphaNeg[1] + alphaPos[1],
    }
    return alphaMetric


def EnsoAmpl(sstfile, sstname, sstbox, **kwargs):
    """
    The EnsoAmpl() function computes the standard deviation of 'sstbox' sstA (usually the standard deviation of nino3
    sstA)

    Author:	Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Yann Planton : yann.planton@locean-ipsl.upmc.fr

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstbox: string
        name of box ('nino3', 'nino3.4', 'nino4') for SST
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return amplMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO amplitude'
    Units = 'C'
    Method = 'Standard deviation of ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        if len(sst) < kwargs['min_time_steps']:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoAmpl', len(sst), kwargs['min_time_steps'], INSPECTstack())

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the standard deviation
    sstStd = Std(sst)

    # Standard Error of the Standard Deviation (function of nyears)
    sstStdErr = sstStd / NUMPYsqrt(yearN)

    # Create output
    amplMetric = {
        'name': Name, 'value': sstStd, 'value_error': sstStdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref,
    }
    return amplMetric


def EnsoMu(sstfile, tauxfile, sstname, tauxname, sstbox, tauxbox, **kwargs):
    """
    The EnsoMu() function computes the regression of 'tauxbox' tauxA (surface downward zonal stress anomalies) over
    'sstbox' sstA (usually the regression of nino4 tauxA over nino3 sstA)

    Author:	Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Yann Planton : yann.planton@locean-ipsl.upmc.fr

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param tauxfile: string
        path_to/filename of the file (NetCDF) of the TAUX dataset
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param tauxname: string
        name of TAUX variable (taux, tauu) in 'tauxfile'
    :param sstbox: string
        name of box ('nino3') for SST
    :param tauxbox: string
        name of box ('nino4') for TAUX
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return muMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Bjerknes feedback (mu)'
    Units = '10-3 N/m2/C'
    Method = 'Regression of ' + tauxbox + ' tauxA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    taux = ReadSelectRegionCheckUnits(tauxfile, tauxname, 'wind stress', box=tauxbox, **kwargs)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, taux = CheckTime(sst, taux, metric_name='EnsoMu', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)
    taux, unneeded = PrePocessTS(taux, Method, average='horizontal', compute_anom=True, **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    mu, muPos, muNeg = LinearRegressionAndNonlinearity(taux, sst, return_stderr=True)

    # Change units
    mu = [OperationMultiply(mu[0], 1000.), OperationMultiply(mu[1], 1000.)]
    muPos = [OperationMultiply(muPos[0], 1000.), OperationMultiply(muPos[1], 1000.)]
    muNeg = [OperationMultiply(muNeg[0], 1000.), OperationMultiply(muNeg[1], 1000.)]

    # Create output
    muMetric = {
        'name': Name, 'value': mu[0], 'value_error': mu[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': muNeg[0] - muPos[0],
        'nonlinearity_error': muNeg[1] + muPos[1],
    }
    return muMetric


def EnsoRMSE(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoRMSE() function computes the SST spatial root mean square error (RMSE) in a 'box' (usually the tropical
    Pacific)

    Inputs:
    ------
    :param sstfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemodel: string
        name of SST variable (tos, ts) in 'sstfilemodel'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param box: string
        name of box ('tropical_pacific') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param regridding: dict, optional
        see EnsoUvcdatToolsLib.TwoVarRegrid and EnsoUvcdatToolsLib.Regrid for options
        the aim if to specify if the model is regridded toward the observations or vice versa, of if both model and
        observations are regridded toward another grid
        interpolation tool and method can be specified
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds_model: tuple, optional
        tuple of the first and last dates to extract from the modeled SST file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed SST file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return rmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_model', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO RMSE'
    Units = 'C'
    Method = 'Spatial root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bnds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bnds=kwargs['time_bounds_obs'], **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "the modeled time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "the observed time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    sst_model, Method = PrePocessTS(sst_model, Method, average='time', compute_anom=False, **kwargs)
    sst_obs, unneeded = PrePocessTS(sst_obs, Method, average='time', compute_anom=False, **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_to_obs', 'obs_to_model', 'model_and_obs_to_newgrid', 'newgrid', 'missing', 'order',
                      'mask', 'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, **kwargs['regridding'])

    # Averages temporally (if not performed earlier) and computes the root mean square difference
    sstRmse = SpatialRms(sst_model, sst_obs, centered=centered_rmse)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return rmseMetric


def EnsoSeasonality(sstfile, sstname, box, **kwargs):
    """
    The EnsoSeasonality() function computes ratio between the November-December-January (NDJ) and March-April-May (MAM)
    average standard deviation of 'box' sstA (usually nino3 sstA)

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST dataset
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param box: string
        name of box ('nino3') for SST
    usual kwargs:
    :param detrending: dict, optional
        see EnsoUvcdatToolsLib.Detrend for options
        the aim if to specify if the trend must be removed
        detrending method can be specified
        default value is False
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
        default value is None
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
        default value is None
    :param normalization: boolean, optional
        True to normalize by the standard deviation (needs the frequency to be defined), if you don't want it pass
        anything but true
        default value is False
    :param smoothing: dict, optional
        see EnsoUvcdatToolsLib.Smoothing for options
        the aim if to specify if variables are smoothed (running mean)
        smoothing axis, window and method can be specified
        default value is False
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return SeaMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO Seasonality'
    Units = ''
    Method = 'Ratio between NDJ and MAM standard deviation ' + box + ' sstA'
    Ref = 'Using CDAT std dev calculation'
    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst) < mini:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoSeasonality', len(sst), mini, INSPECTstack())

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=False, **kwargs)

    # Seasonal mean
    sst_NDJ = SeasonalMean(sst, 'NDJ', compute_anom=True)
    sst_MAM = SeasonalMean(sst, 'MAM', compute_anom=True)

    # Compute std dev and ratio
    sst_NDJ_std = Std(sst_NDJ)
    sst_MAM_std = Std(sst_MAM)
    ratioStd = float(sst_NDJ_std / sst_MAM_std)

    # Standard Error of the Standard Deviation (function of nyears)
    sst_NDJ_std_err = sst_NDJ_std / NUMPYsqrt(yearN - 1)
    sst_MAM_std_err = sst_MAM_std / NUMPYsqrt(yearN)

    # The error (dy) on ratio ('y = x/z'): dy = (z*dx + x*dz) / z2
    ratio_std_err = float((sst_MAM_std * sst_NDJ_std_err + sst_NDJ_std * sst_MAM_std_err) / NUMPYsquare(sst_MAM_std))

    # Create output
    seaMetric = {
        'name': Name, 'value': ratioStd, 'value_error': ratio_std_err, 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref,
    }
    return seaMetric


def EnsoCompositeTS(sstfile, sstname, box, season, threshold, nbr_years_window, nino=True, **kwargs):
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    if nino:
        event = 'Nino'
    else:
        event = 'Nina'
    Name = event + ' Composite Time Series'
    Method = event + ' events detected using ' + box + ' sstA during ' + season + ', time series of ' + \
             nbr_years_window + ' years (centered on events)'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT std dev calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=box, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst) < mini:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoCompositeTS', len(sst), mini, INSPECTstack())

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst, Method = PrePocessTS(sst, Method, average='horizontal', compute_anom=True, **kwargs)

    # Lists event years
    list_event_years = DetectEvents(sst, season, threshold, normalization=kwargs['normalization'], nino=nino)
    # Removes years that are too close to the time bounds to select the proper 'nbr_years_window'
    tmp_years = deepcopy(list_event_years)
    for yy in tmp_years:
        if yy < int(actualtimebounds[0].split('-')[0]) + (nbr_years_window / 2) -1:
            list_event_years.remove(yy)
        if yy > int(actualtimebounds[-1].split('-')[0]) - (nbr_years_window / 2):
            list_event_years.remove(yy)

    # composites
    composite = list()
    for yy in list_event_years:
        timebnds = (str(yy-2) + '-01-01 00:00:00.0', str(yy+3) + '-12-31  23:59:60.0')
        composite.append(sst(time=timebnds))
    return
# ---------------------------------------------------------------------------------------------------------------------#




