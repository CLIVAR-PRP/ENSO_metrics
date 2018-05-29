# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import sqrt as NUMPYsqrt
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings
from EnsoUvcdatToolsLib import AverageMeridional, AverageZonal, CheckTime, Composite, DetectEvents,\
    LinearRegressionAndNonlinearity, MyDerive, PreProcessTS, ReadAreaSelectRegion, ReadSelectRegionCheckUnits,\
    RmsHorizontal, RmsMeridional, RmsTemporal, RmsZonal, SeasonalMean, Std, TimeBounds, TwoVarRegrid
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
    print '\033[92m' + str().ljust(10) + "EnsoAlphaLhf" + '\033[0m'
    print '\033[92m' + str().ljust(15) + "model var is " + str(sstname) + ", file is " + str(sstfile) + '\033[0m'
    print '\033[92m' + str().ljust(15) + "obs var is " + str(lhfname) + ", file is " + str(lhffile) + '\033[0m'
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    lhf = ReadSelectRegionCheckUnits(lhffile, lhfname, 'heat flux', box=lhfbox, **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadSelectRegionCheckUnits" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.shape = " + str(sst.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.timebounds = " + str(TimeBounds(sst)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.axes = " + str([ax.id for ax in sst.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.grid = " + str(sst.getGrid()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.shape = " + str(lhf.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.timebounds = " + str(TimeBounds(lhf)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.axes = " + str([ax.id for ax in lhf.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.grid = " + str(lhf.getGrid()) + '\033[0m'
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    lhf_areacell = ReadAreaSelectRegion(lhffile, box=lhfbox, **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadAreaSelectRegion" + '\033[0m'
    if sst_areacell is None:
        print '\033[92m' + str().ljust(20) + "No sst_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "sst_areacell.shape = " + str(sst_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "sst_areacell.axes = " + str([ax.id for ax in sst_areacell.getAxisList()])\
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "sst_areacell.grid = " + str(sst_areacell.getGrid()) + '\033[0m'
    if lhf_areacell is None:
        print '\033[92m' + str().ljust(20) + "No lhf_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "lhf_areacell.shape = " + str(lhf_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "lhf_areacell.axes = " + str([ax.id for ax in lhf_areacell.getAxisList()])\
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "lhf_areacell.grid = " + str(lhf_areacell.getGrid()) + '\033[0m'

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lhf = CheckTime(sst, lhf, metric_name='EnsoAlphaLhf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    print '\033[92m' + str().ljust(15) + "after CheckTime" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.shape = " + str(sst.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.timebounds = " + str(TimeBounds(sst)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.axes = " + str([ax.id for ax in sst.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.grid = " + str(sst.getGrid()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.shape = " + str(lhf.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.timebounds = " + str(TimeBounds(lhf)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.axes = " + str([ax.id for ax in lhf.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lhf.grid = " + str(lhf.getGrid()) + '\033[0m'
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    lhf, unneeded = PreProcessTS(lhf, Method, areacell=lhf_areacell, average='horizontal', compute_anom=True, **kwargs)

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
    print '\033[92m' + str().ljust(10) + "EnsoAlphaLwr" + '\033[0m'
    print '\033[92m' + str().ljust(15) + "sst var is " + str(sstname) + ", file is " + str(sstfile) + '\033[0m'
    print '\033[92m' + str().ljust(15) + "lwr var is " + str(lwrname) + ", file is " + str(lwrfile) + '\033[0m'
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(lwrfile, basestring):
        dict_var[lwrname] = ReadSelectRegionCheckUnits(lwrfile, lwrname, 'heat flux', box=lwrbox, **kwargs)
    else:
        for ii in range(len(lwrfile)):
            filename, varname = lwrfile[ii], lwrname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=lwrbox, **kwargs)
    lwr = MyDerive(kwargs['project_interpreter_var2'], 'lwr', dict_var)
    print '\033[92m' + str().ljust(15) + "after ReadSelectRegionCheckUnits" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.shape = " + str(sst.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.timebounds = " + str(TimeBounds(sst)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.axes = " + str([ax.id for ax in sst.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.shape = " + str(lwr.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.timebounds = " + str(TimeBounds(lwr)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.axes = " + str([ax.id for ax in lwr.getAxisList()]) + '\033[0m'
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    if isinstance(lwrfile, basestring):
        lwr_areacell = ReadAreaSelectRegion(lwrfile, box=lwrbox, **kwargs)
    else:
        for ii in range(len(lwrfile)):
            lwr_areacell = ReadAreaSelectRegion(lwrfile[ii], box=lwrbox, **kwargs)
            if lwr_areacell is not None:
               break
    print '\033[92m' + str().ljust(15) + "after ReadAreaSelectRegion" + '\033[0m'
    if sst_areacell is None:
        print '\033[92m' + str().ljust(20) + "No sst_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "sst_areacell.shape = " + str(sst_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "sst_areacell.axes = " + str([ax.id for ax in sst_areacell.getAxisList()]) \
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "sst_areacell.grid = " + str(sst_areacell.getGrid()) + '\033[0m'
    if lwr_areacell is None:
        print '\033[92m' + str().ljust(20) + "No lwr_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "lwr_areacell.shape = " + str(lwr_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "lwr_areacell.axes = " + str([ax.id for ax in lwr_areacell.getAxisList()]) \
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "lwr_areacell.grid = " + str(lwr_areacell.getGrid()) + '\033[0m'

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lwr = CheckTime(sst, lwr, metric_name='EnsoAlphaLwr', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)
    print '\033[92m' + str().ljust(15) + "after CheckTime" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.shape = " + str(sst.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.timebounds = " + str(actualtimebounds) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.axes = " + str([ax.id for ax in sst.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.shape = " + str(lwr.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.timebounds = " + str(TimeBounds(lwr)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.axes = " + str([ax.id for ax in lwr.getAxisList()]) + '\033[0m'

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    lwr, unneeded = PreProcessTS(lwr, Method, areacell=lwr_areacell, average='horizontal', compute_anom=True, **kwargs)
    print '\033[92m' + str().ljust(15) + "after PreProcessTS" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.shape = " + str(sst.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.timebounds = " + str(TimeBounds(sst)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "sst.axes = " + str([ax.id for ax in sst.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.shape = " + str(lwr.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.timebounds = " + str(TimeBounds(lwr)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "lwr.axes = " + str([ax.id for ax in lwr.getAxisList()]) + '\033[0m'

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
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    shf_areacell = ReadAreaSelectRegion(shffile, box=shfbox, **kwargs)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, shf = CheckTime(sst, shf, metric_name='EnsoAlphaShf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    shf, unneeded = PreProcessTS(shf, Method, areacell=shf_areacell, average='horizontal', compute_anom=True, **kwargs)

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
    swr = MyDerive(kwargs['project_interpreter_var2'], 'swr', dict_var)
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    if isinstance(swrfile, basestring):
        swr_areacell = ReadAreaSelectRegion(swrfile, box=swrbox, **kwargs)
    else:
        for ii in range(len(swrfile)):
            swr_areacell = ReadAreaSelectRegion(swrfile[ii], box=swrbox, **kwargs)
            if swr_areacell is not None:
                break

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, swr = CheckTime(sst, swr, metric_name='EnsoAlphaSwr', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    swr, unneeded = PreProcessTS(swr, Method, areacell=swr_areacell, average='horizontal', compute_anom=True, **kwargs)

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
    thf = MyDerive(kwargs['project_interpreter_var2'], 'thf', dict_var)
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    if isinstance(thffile, basestring):
        thf_areacell = ReadAreaSelectRegion(thffile, box=thfbox, **kwargs)
    else:
        for ii in range(len(thffile)):
            thf_areacell = ReadAreaSelectRegion(thffile[ii], box=thfbox, **kwargs)
            if thf_areacell is not None:
               break

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, thf = CheckTime(sst, thf, metric_name='EnsoAlphaThf', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    thf, unneeded = PreProcessTS(thf, Method, areacell=thf_areacell, average='horizontal', compute_anom=True, **kwargs)

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
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        if len(sst) < kwargs['min_time_steps']:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoAmpl', len(sst), kwargs['min_time_steps'], INSPECTstack())

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)

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
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=sstbox, **kwargs)
    taux_areacell = ReadAreaSelectRegion(tauxfile, box=tauxbox, **kwargs)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, taux = CheckTime(sst, taux, metric_name='EnsoMu', **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
    taux, unneeded = PreProcessTS(taux, Method, areacell=taux_areacell, average='horizontal', compute_anom=True,
                                  **kwargs)

    # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
    mu, muPos, muNeg = LinearRegressionAndNonlinearity(taux, sst, return_stderr=True)

    # Change units
    mu = [mu[0] * 1000., mu[1] * 1000.]
    muPos = [muPos[0] * 1000., muPos[1] * 1000.]
    muNeg = [muNeg[0] * 1000., muNeg[1] * 1000.]

    # Create output
    muMetric = {
        'name': Name, 'value': mu[0], 'value_error': mu[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': muNeg[0] - muPos[0],
        'nonlinearity_error': muNeg[1] + muPos[1],
    }
    return muMetric


def EnsoRmse(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoRmse() function computes the SST spatial root mean square error (RMSE) in a 'box' (usually the tropical
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
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoRmse: the modeled time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoRmse: the observed time-period is too short: " + str(len(sst_obs))
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
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average='time', compute_anom=False,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                     **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=box, **kwargs['regridding'])

    # Computes the root mean square difference
    sstRmse = RmsHorizontal(sst_model, sst_obs, centered=centered_rmse)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return rmseMetric


def EnsoLatRmse(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoLatRmse() function computes the SST meridional (latitude) root mean square error (RMSE) in a 'box'
    (usually 'nino3_reduced')

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
        name of box ('nino3_reduced') for SST
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
    :return LatRmseMetric: dict
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
    Name = 'ENSO Meridional RMSE'
    Units = 'C'
    Method = 'Meridional root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    print '\033[92m' + str().ljust(10) + "EnsoLatRmse" + '\033[0m'
    print '\033[92m' + str().ljust(15) + "model var is " + str(sstnamemodel) + ", file is "\
          + str(sstfilemodel) + '\033[0m'
    print '\033[92m' + str().ljust(15) + "obs var is " + str(sstnameobs) + ", file is " + str(sstfileobs) + '\033[0m'
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadSelectRegionCheckUnits" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.timebounds = " + str(TimeBounds(sst_model)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.latitude = " + str(sst_model.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.longitude = " + str(sst_model.getLongitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.timebounds = " + str(TimeBounds(sst_obs)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.latitude = " + str(sst_obs.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.longitude = " + str(sst_obs.getLongitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=box, **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadAreaSelectRegion" + '\033[0m'
    if model_areacell is None:
        print '\033[92m' + str().ljust(20) + "No model_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "model_areacell.shape = " + str(model_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model_areacell.axes = "\
              + str([ax.id for ax in model_areacell.getAxisList()]) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model_areacell.grid = " + str(model_areacell.getGrid()) + '\033[0m'
    if obs_areacell is None:
        print '\033[92m' + str().ljust(20) + "No obs_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "obs_areacell.shape = " + str(obs_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs_areacell.axes = "\
              + str([ax.id for ax in obs_areacell.getAxisList()]) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs_areacell.grid = " + str(obs_areacell.getGrid()) + '\033[0m'

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoLatRmse: the modeled time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoLatRmse: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average=['time'], compute_anom=False,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False, **kwargs)
    print '\033[92m' + str().ljust(15) + "after PreProcessTS" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.latitude = " + str(sst_model.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.latitude = " + str(sst_obs.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=box, **kwargs['regridding'])
        print '\033[92m' + str().ljust(15) + "after TwoVarRegrid" + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()])\
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'

    # Zonal average
    sst_model = AverageZonal(sst_model)
    sst_obs = AverageZonal(sst_obs)

    # Computes the root mean square difference
    sstRmse = RmsMeridional(sst_model, sst_obs, centered=centered_rmse)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return LatRmseMetric


def EnsoLonRmse(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoLonRmse() function computes the SST zonal (longitude) root mean square error (RMSE) in a 'box'
    (usually the Equatorial Pacific)

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
        name of box ('equatorial_pacific') for SST
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
    :return LonRmseMetric: dict
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
    Name = 'ENSO Zonal RMSE'
    Units = 'C'
    Method = 'Meridional root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    print '\033[92m' + str().ljust(10) + "EnsoLonRmse" + '\033[0m'
    print '\033[92m' + str().ljust(15) + "model var is " + str(sstnamemodel) + ", file is " \
          + str(sstfilemodel) + '\033[0m'
    print '\033[92m' + str().ljust(15) + "obs var is " + str(sstnameobs) + ", file is " + str(sstfileobs) + '\033[0m'
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadSelectRegionCheckUnits" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.timebounds = " + str(TimeBounds(sst_model)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.latitude = " + str(sst_model.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.longitude = " + str(sst_model.getLongitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.timebounds = " + str(TimeBounds(sst_obs)) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.latitude = " + str(sst_obs.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.longitude = " + str(sst_obs.getLongitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=box, **kwargs)
    print '\033[92m' + str().ljust(15) + "after ReadAreaSelectRegion" + '\033[0m'
    if model_areacell is None:
        print '\033[92m' + str().ljust(20) + "No model_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "model_areacell.shape = " + str(model_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model_areacell.axes = " \
              + str([ax.id for ax in model_areacell.getAxisList()]) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model_areacell.grid = " + str(model_areacell.getGrid()) + '\033[0m'
    if obs_areacell is None:
        print '\033[92m' + str().ljust(20) + "No obs_areacell" + '\033[0m'
    else:
        print '\033[92m' + str().ljust(20) + "obs_areacell.shape = " + str(obs_areacell.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs_areacell.axes = " \
              + str([ax.id for ax in obs_areacell.getAxisList()]) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs_areacell.grid = " + str(obs_areacell.getGrid()) + '\033[0m'

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoLonRmse: the modeled time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoLonRmse: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average=['time'], compute_anom=False,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False, **kwargs)
    print '\033[92m' + str().ljust(15) + "after PreProcessTS" + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.latitude = " + str(sst_model.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.latitude = " + str(sst_obs.getLatitude()) + '\033[0m'
    print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=box, **kwargs['regridding'])
        print '\033[92m' + str().ljust(15) + "after TwoVarRegrid" + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model.shape = " + str(sst_model.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "model.axes = " + str([ax.id for ax in sst_model.getAxisList()])\
              + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs.shape = " + str(sst_obs.shape) + '\033[0m'
        print '\033[92m' + str().ljust(20) + "obs.axes = " + str([ax.id for ax in sst_obs.getAxisList()]) + '\033[0m'

    # Meridional average
    sst_model = AverageMeridional(sst_model)
    sst_obs = AverageMeridional(sst_obs)

    # Computes the root mean square difference
    sstRmse = RmsZonal(sst_model, sst_obs, centered=centered_rmse)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return LonRmseMetric


def EnsoPrRmse(prfilemodel, prnamemodel, prfileobs, prnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoPrRmse() function computes the PR (precipitation) spatial root mean square error (RMSE) in a 'box' (usually
    the tropical Pacific)

    Inputs:
    ------
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PR variable (pr, precip) in 'prfilemodel'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param box: string
        name of box ('tropical_pacific') for PR
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
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return PrRmseMetric: dict
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
    Name = 'ENSO Pr RMSE'
    Units = 'mm/day'
    Method = 'Spatial root mean square error of ' + box + ' Pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=box,
                                          time_bounds=kwargs['time_bounds_model'], **kwargs)
    pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=box,
                                        time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(prfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(prfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(pr_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrRmse: the modeled time-period is too short: " + str(len(pr_model))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(pr_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrRmse: the observed time-period is too short: " + str(len(pr_obs))
                            + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = pr_model.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(pr_model)
    actualtimeboundsobs = TimeBounds(pr_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, average='time', compute_anom=False,
                                    **kwargs)
    pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='time', compute_anom=False, **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        pr_model, pr_obs, Method = TwoVarRegrid(pr_model, pr_obs, Method, region=box, **kwargs['regridding'])

    # Computes the root mean square difference
    pr_rmse = RmsHorizontal(pr_model, pr_obs, centered=centered_rmse)

    # Create output
    PrRmseMetric = {
        'name': Name, 'value': pr_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return PrRmseMetric


def EnsoPrLatRmse(prfilemodel, prnamemodel, prfileobs, prnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoPrLatRmse() function computes the PR (zonal wind stress) meridional (latitude) root mean square error (RMSE)
    in a 'box' (usually 'nino3_reduced')

    Inputs:
    ------
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PR variable (pr, precip) in 'prfilemodel'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param box: string
        name of box ('nino3_reduced') for PR
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
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return PrLatRmseMetric: dict
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
    Name = 'ENSO Pr Meridional RMSE'
    Units = 'mm/day'
    Method = 'Meridional root mean square error of ' + box + ' Pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=box,
                                          time_bounds=kwargs['time_bounds_model'], **kwargs)
    pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=box,
                                        time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(prfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(prfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(pr_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrLatRmse: the modeled time-period is too short: "
                            + str(len(pr_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(pr_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrLatRmse: the observed time-period is too short: "
                            + str(len(pr_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = pr_model.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(pr_model)
    actualtimeboundsobs = TimeBounds(pr_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, average=['time'], compute_anom=False,
                                    **kwargs)
    pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False, **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        pr_model, pr_obs, Method = TwoVarRegrid(pr_model, pr_obs, Method, region=box, **kwargs['regridding'])

    # Zonal average
    pr_model = AverageZonal(pr_model)
    pr_obs = AverageZonal(pr_obs)

    # Computes the root mean square difference
    pr_rmse = RmsMeridional(pr_model, pr_obs, centered=centered_rmse)

    # Create output
    PrLatRmseMetric = {
        'name': Name, 'value': pr_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return PrLatRmseMetric


def EnsoPrLonRmse(prfilemodel, prnamemodel, prfileobs, prnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoPrLonRmse() function computes the PR (zonal wind stress) zonal (longitude) root mean square error (RMSE) in
    a 'box' (usually the Equatorial Pacific)

    Inputs:
    ------
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PR variable (pr, precip) in 'prfilemodel'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param box: string
        name of box ('equatorial_pacific') for PR
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
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return PrLonRmseMetric: dict
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
    Name = 'ENSO Pr Zonal RMSE'
    Units = 'mm/day'
    Method = 'Meridional root mean square error of ' + box + ' Pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=box,
                                          time_bounds=kwargs['time_bounds_model'], **kwargs)
    pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=box,
                                        time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(prfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(prfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(pr_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrLonRmse: the modeled time-period is too short: "
                            + str(len(pr_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(pr_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrLonRmse: the observed time-period is too short: "
                            + str(len(pr_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = pr_model.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(pr_model)
    actualtimeboundsobs = TimeBounds(pr_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, average=['time'], compute_anom=False,
                                    **kwargs)
    pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False, **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        pr_model, pr_obs, Method = TwoVarRegrid(pr_model, pr_obs, Method, region=box, **kwargs['regridding'])

    # Meridional average
    pr_model = AverageMeridional(pr_model)
    pr_obs = AverageMeridional(pr_obs)

    # Computes the root mean square difference
    pr_rmse = RmsZonal(pr_model, pr_obs, centered=centered_rmse)

    # Create output
    PrLonRmseMetric = {
        'name': Name, 'value': pr_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return PrLonRmseMetric


def EnsoTauxRmse(tauxfilemodel, tauxnamemodel, tauxfileobs, tauxnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoTauxRmse() function computes the TAUX (zonal wind stress) spatial root mean square error (RMSE) in a 'box'
    (usually the tropical Pacific)

    Inputs:
    ------
    :param tauxfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemodel: string
        name of TAUX variable (taux, tauu) in 'tauxfilemodel'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param box: string
        name of box ('tropical_pacific') for TAUX
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
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return TauxRmseMetric: dict
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
    Name = 'ENSO Taux RMSE'
    Units = 'N/m2'
    Method = 'Spatial root mean square error of ' + box + ' Taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    taux_model = ReadSelectRegionCheckUnits(tauxfilemodel, tauxnamemodel, 'wind stress', box=box,
                                            time_bounds=kwargs['time_bounds_model'], **kwargs)
    taux_obs = ReadSelectRegionCheckUnits(tauxfileobs, tauxnameobs, 'wind stress', box=box,
                                          time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(tauxfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(tauxfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(taux_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxRmse: the modeled time-period is too short: "
                            + str(len(taux_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(taux_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxRmse: the observed time-period is too short: "
                            + str(len(taux_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = taux_model.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(taux_model)
    actualtimeboundsobs = TimeBounds(taux_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    taux_model, Method = PreProcessTS(taux_model, Method, areacell=model_areacell, average='time', compute_anom=False,
                                      **kwargs)
    taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average='time', compute_anom=False, **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        taux_model, taux_obs, Method = TwoVarRegrid(taux_model, taux_obs, Method, region=box, **kwargs['regridding'])

    # Computes the root mean square difference
    taux_rmse = RmsHorizontal(taux_model, taux_obs, centered=centered_rmse)

    # Create output
    TauxRmseMetric = {
        'name': Name, 'value': taux_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return TauxRmseMetric


def EnsoTauxLatRmse(tauxfilemodel, tauxnamemodel, tauxfileobs, tauxnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoTauxLatRmse() function computes the TAUX (zonal wind stress) meridional (latitude) root mean square error
    (RMSE) in a 'box' (usually 'nino3')

    Inputs:
    ------
    :param tauxfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemodel: string
        name of TAUX variable (taux, tauu) in 'tauxfilemodel'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param box: string
        name of box ('nino3') for TAUX
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
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return TauxLatRmseMetric: dict
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
    Name = 'ENSO Taux Meridional RMSE'
    Units = 'N/m2'
    Method = 'Meridional root mean square error of ' + box + ' Taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    taux_model = ReadSelectRegionCheckUnits(tauxfilemodel, tauxnamemodel, 'wind stress', box=box,
                                            time_bounds=kwargs['time_bounds_model'], **kwargs)
    taux_obs = ReadSelectRegionCheckUnits(tauxfileobs, tauxnameobs, 'wind stress', box=box,
                                          time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(tauxfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(tauxfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(taux_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxLatRmse: the modeled time-period is too short: "
                            + str(len(taux_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(taux_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxLatRmse: the observed time-period is too short: "
                            + str(len(taux_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = taux_model.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(taux_model)
    actualtimeboundsobs = TimeBounds(taux_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    taux_model, Method = PreProcessTS(taux_model, Method, areacell=model_areacell, average=['time'],
                                      compute_anom=False, **kwargs)
    taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False,
                                      **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        taux_model, taux_obs, Method = TwoVarRegrid(taux_model, taux_obs, Method, region=box, **kwargs['regridding'])

    # Zonal average
    taux_model = AverageZonal(taux_model)
    taux_obs = AverageZonal(taux_obs)

    # Computes the root mean square difference
    taux_rmse = RmsMeridional(taux_model, taux_obs, centered=centered_rmse)

    # Create output
    TauxLatRmseMetric = {
        'name': Name, 'value': taux_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return TauxLatRmseMetric


def EnsoTauxLonRmse(tauxfilemodel, tauxnamemodel, tauxfileobs, tauxnameobs, box, centered_rmse=0, **kwargs):
    """
    The EnsoTauxLonRmse() function computes the TAUX (zonal wind stress) zonal (longitude) root mean square error (RMSE)
    in a 'box' (usually the Equatorial Pacific)

    Inputs:
    ------
    :param tauxfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemodel: string
        name of TAUX variable (taux, tauu) in 'tauxfilemodel'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param box: string
        name of box ('equatorial_pacific') for TAUX
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
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return TauxLonRmseMetric: dict
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
    Name = 'ENSO Taux Zonal RMSE'
    Units = 'N/m2'
    Method = 'Meridional root mean square error of ' + box + ' Taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    taux_model = ReadSelectRegionCheckUnits(tauxfilemodel, tauxnamemodel, 'wind stress', box=box,
                                            time_bounds=kwargs['time_bounds_model'], **kwargs)
    taux_obs = ReadSelectRegionCheckUnits(tauxfileobs, tauxnameobs, 'wind stress', box=box,
                                          time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(tauxfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(tauxfileobs, box=box, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(taux_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxLonRmse: the modeled time-period is too short: "
                            + str(len(taux_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(taux_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoTauxLonRmse: the observed time-period is too short: "
                            + str(len(taux_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = taux_model.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(taux_model)
    actualtimeboundsobs = TimeBounds(taux_obs)

    # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    # here only the detrending (if applicable) and time averaging are performed
    taux_model, Method = PreProcessTS(taux_model, Method, areacell=model_areacell, average=['time'],
                                      compute_anom=False, **kwargs)
    taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average=['time'], compute_anom=False,
                                      **kwargs)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        taux_model, taux_obs, Method = TwoVarRegrid(taux_model, taux_obs, Method, region=box, **kwargs['regridding'])

    # Meridional average
    taux_model = AverageMeridional(taux_model)
    taux_obs = AverageMeridional(taux_obs)

    # Computes the root mean square difference
    taux_rmse = RmsZonal(taux_model, taux_obs, centered=centered_rmse)

    # Create output
    TauxLonRmseMetric = {
        'name': Name, 'value': taux_rmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return TauxLonRmseMetric


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
    # Read areacell
    sst_areacell = ReadAreaSelectRegion(sstfile, box=box, **kwargs)

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
    sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=False, **kwargs)

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


def NinaCompositeLon(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, event_definition, centered_rmse=0,
                     **kwargs):
    """
    The NinaCompositeLon() function computes a zonal composite of La Niña events during the peak of the event
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA > 'threshold' during 'season' are considered as La Niña events
    Then the zonal SSTA at the peak of the event is composited for each selected event

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
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
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
    :return NinaLonMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, events_model, events_observations,
        time_frequency, time_period_model, time_period_observations, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_model',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nina Composite Time Series'
    lat = ReferenceRegions('equatorial_pacific')['latitude']
    Method = 'Nina events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=region_ev, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=region_ev, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinaCompositeLon: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinaCompositeLon: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=True,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)

    # Lists event years
    event_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
    event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=False)

    # ------------------------------------------------
    # compute composite
    # ------------------------------------------------
    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=box, **kwargs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average=False, compute_anom=False,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False, **kwargs)

    # Seasonal mean
    sst_model = SeasonalMean(sst_model, season_ev, compute_anom=False)
    sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=False)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=box, **kwargs['regridding'])

    # Meridional average
    sst_model = AverageMeridional(sst_model)
    sst_obs = AverageMeridional(sst_obs)

    # composites
    composite_model = Composite(sst_model, event_years_model, kwargs['frequency'])
    composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'])

    # Computes the root mean square difference
    compositeRmse = RmsZonal(composite_model, composite_obs, centered=centered_rmse)

    # Create output
    NinaLonMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'events_model': event_years_model,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
    }
    return NinaLonMetric


def NinoCompositeLon(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, event_definition, centered_rmse=0,
                     **kwargs):
    """
    The NinoCompositeLon() function computes a zonal composite of El Niño events during the peak of the event
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA > 'threshold' during 'season' are considered as El Niño events
    Then the zonal SSTA at the peak of the event is composited for each selected event

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
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
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
    :return NinoLonMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, events_model, events_observations,
        time_frequency, time_period_model, time_period_observations, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_model',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino Composite Time Series'
    lat = ReferenceRegions('equatorial_pacific')['latitude']
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=region_ev, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=region_ev, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinoCompositeLon: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinoCompositeLon: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=True,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)

    # Lists event years
    event_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)

    # ------------------------------------------------
    # compute composite
    # ------------------------------------------------
    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=box, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=box, **kwargs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average=False, compute_anom=False,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False, **kwargs)

    # Seasonal mean
    sst_model = SeasonalMean(sst_model, season_ev, compute_anom=False)
    sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=False)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=box, **kwargs['regridding'])

    # Meridional average
    sst_model = AverageMeridional(sst_model)
    sst_obs = AverageMeridional(sst_obs)

    # composites
    composite_model = Composite(sst_model, event_years_model, kwargs['frequency'])
    composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'])

    # Computes the root mean square difference
    compositeRmse = RmsZonal(composite_model, composite_obs, centered=centered_rmse)

    # Create output
    NinoLonMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'events_model': event_years_model,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
    }
    return NinoLonMetric


def NinaCompositeTS(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, event_definition, nbr_years_window,
                    centered_rmse=0, **kwargs):
    """
    The NinaCompositeTS() function computes a time composite of La Niña events
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA < 'threshold' during 'season' are considered as La Niña events
        Then a 'nbr_years_window' long time series centered on selected events is composited for each selected event

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
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
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
    :return NinaTsMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, events_model, events_observations,
        time_frequency, time_period_model, time_period_observations, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_model',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nina Composite Time Series'
    Method = 'Nina events = ' + region_ev + ' sstA < ' + str(threshold) + ' during ' + season_ev + ', time series of '\
             + str(nbr_years_window) + ' years (centered on events)'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT rms (uncentered and biased) calculation'

    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=region_ev, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=region_ev, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinaCompositeTS: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinaCompositeTS: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average='horizontal',
                                     compute_anom=True, **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)

    # Lists event years
    event_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
    event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=False)

    # composites
    composite_model = Composite(sst_model, event_years_model, kwargs['frequency'], nbr_years_window=nbr_years_window)
    composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'], nbr_years_window=nbr_years_window)

    # Computes the root mean square difference
    compositeRmse = RmsTemporal(composite_model, composite_obs, centered=centered_rmse)

    # Create output
    NinaTsMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'events_model': event_years_model,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
    }
    return NinaTsMetric


def NinoCompositeTS(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, event_definition, nbr_years_window,
                    centered_rmse=0, **kwargs):
    """
    The NinoCompositeTS() function computes a time composite of El Niño events
    SSTA averaged in 'box' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA > 'threshold' during 'season' are considered as El Niño events
        Then a 'nbr_years_window' long time series centered on selected events is composited for each selected event

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
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
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
    :return NinoTsMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, events_model, events_observations,
        time_frequency, time_period_model, time_period_observations, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_model',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino Composite Time Series'
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', time series of '\
             + str(nbr_years_window) + ' years (centered on events)'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT rms (uncentered and biased) calculation'

    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_model'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    # Read areacell
    model_areacell = ReadAreaSelectRegion(sstfilemodel, box=region_ev, **kwargs)
    obs_areacell = ReadAreaSelectRegion(sstfileobs, box=region_ev, **kwargs)

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinoCompositeTS: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "NinoCompositeTS: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)

    # Lists event years
    event_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)

    # composites
    composite_model = Composite(sst_model, event_years_model, kwargs['frequency'], nbr_years_window=nbr_years_window)
    composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'], nbr_years_window=nbr_years_window)

    # Computes the root mean square difference
    compositeRmse = RmsTemporal(composite_model, composite_obs, centered=centered_rmse)

    # Create output
    NinoTsMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'events_model': event_years_model,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
    }
    return NinoTsMetric
# ---------------------------------------------------------------------------------------------------------------------#




