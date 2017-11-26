# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy import sqrt as NUMPYsqrt
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import defCollection
import EnsoErrorsWarnings
from EnsoUvcdatToolsLib import CheckTime, Multiply, SeasonalMean, TimeBounds,\
    TimeAnomaliesLinearRegressionAndNonlinearity, TimeAnomaliesStd, ReadSelectRegionCheckUnits,\
    SpatialRms


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO metrics
# These functions have file names and variable names as inputs and metric as output
#
def EnsoAlphaLhf(sstfile, lhffile, sstname, lhfname, sstbox, lhfbox, timebounds=None, frequency=None,
                 mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return alphaLhfMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'Latent feedback (alpha_lh)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + lhfbox + ' lhfA over ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)
    lhf = ReadSelectRegionCheckUnits(lhffile, lhfname, 'heat flux', box=lhfbox, time_bnds=timebounds,
                                     frequency=frequency)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lhf = CheckTime(sst, lhf, frequency=frequency, minimum_length=mintimesteps, metric_name='EnsoAlphaLhf') # -> plante

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the linear regression for all points, for SSTA >=0 and
    # for SSTA<=0
    alphaLhf, alphaLhfPos, alphaLhfNeg = TimeAnomaliesLinearRegressionAndNonlinearity(lhf, sst, return_stderr=True)

    # Create output
    alphaLhfMetric = {
        'name': Name, 'value': alphaLhf[0], 'value_error': alphaLhf[1], 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
        'nonlinearity': alphaLhfNeg[0] - alphaLhfPos[0], 'nonlinearity_error': alphaLhfNeg[1] + alphaLhfPos[1],
    }
    return alphaLhfMetric


def EnsoAlphaLwr(sstfile, lwrfile, sstname, lwrname, sstbox, lwrbox, timebounds=None, frequency=None,
                 mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return alphaLwrMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'Longwave feedback (alpha_lwr)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + lwrbox + ' lwrA over ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)
    if isinstance(lwrfile, list):
        rlds = ReadSelectRegionCheckUnits(lwrfile[0], lwrname[0], 'heat flux', box=lwrbox, time_bnds=timebounds,
                                          frequency=frequency)
        rlus = ReadSelectRegionCheckUnits(lwrfile[1], lwrname[1], 'heat flux', box=lwrbox, time_bnds=timebounds,
                                          frequency=frequency)
        lwr = rlds - rlus
    else:
        lwr = ReadSelectRegionCheckUnits(lwrfile, lwrname, 'heat flux', box=lwrbox, time_bnds=timebounds,
                                         frequency=frequency)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lwr = CheckTime(sst, lwr, frequency=frequency, minimum_length=mintimesteps, metric_name='EnsoAlphaLwr')

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the linear regression for all points, for SSTA >=0 and
    # for SSTA<=0
    alphaLwr, alphaLwrPos, alphaLwrNeg = TimeAnomaliesLinearRegressionAndNonlinearity(lwr, sst, return_stderr=True)

    # Create output
    alphaLwrMetric = {
        'name': Name, 'value': alphaLwr[0], 'value_error': alphaLwr[1], 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
        'nonlinearity': alphaLwrNeg[0] - alphaLwrPos[0], 'nonlinearity_error': alphaLwrNeg[1] + alphaLwrPos[1],
    }
    return alphaLwrMetric


def EnsoAlphaSwr(sstfile, swrfile, sstname, swrname, sstbox, swrbox, timebounds=None, frequency=None,
                 mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return alphaSwrMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'Longwave feedback (alpha_swr)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + swrbox + ' swrA over ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)
    if isinstance(swrfile, list):
        rsds = ReadSelectRegionCheckUnits(swrfile[0], swrname[0], 'heat flux', box=swrbox, time_bnds=timebounds,
                                          frequency=frequency)
        rsus = ReadSelectRegionCheckUnits(swrfile[1], swrname[1], 'heat flux', box=swrbox, time_bnds=timebounds,
                                          frequency=frequency)
        swr = rsds - rsus
    else:
        swr = ReadSelectRegionCheckUnits(swrfile, swrname, 'heat flux', box=swrbox, time_bnds=timebounds,
                                         frequency=frequency)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, swr = CheckTime(sst, swr, frequency=frequency, minimum_length=mintimesteps, metric_name='EnsoAlphaSwr')

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the linear regression for all points, for SSTA >=0 and
    # for SSTA<=0
    alphaSwr, alphaSwrPos, alphaSwrNeg = TimeAnomaliesLinearRegressionAndNonlinearity(swr, sst, return_stderr=True)

    # Create output
    alphaSwrMetric = {
        'name': Name, 'value': alphaSwr[0], 'value_error': alphaSwr[1], 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
        'nonlinearity': alphaSwrNeg[0] - alphaSwrPos[0], 'nonlinearity_error': alphaSwrNeg[1] + alphaSwrPos[1],
    }
    return alphaSwrMetric


def EnsoAlphaThf(sstfile, thffile, sstname, thfname, sstbox, thfbox, timebounds=None, frequency=None, mintimesteps=None):
    """
    The EnsoAlpha() function computes the regression of 'thfbox' thfA (total heat flux anomalies) over 'sstbox' sstA
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return alphaMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'Heat flux feedback (alpha)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + thfbox + ' thfA over ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)
    if isinstance(thffile, list):
        tmp = list()
        for ii in range(len(thffile)):
            val = ReadSelectRegionCheckUnits(thffile[ii], thfname[ii], 'heat flux', box=thfbox,
                                             time_bnds=timebounds, frequency=frequency)
            tmp.append(val)
        for ii in range(len(thffile)):
            try:
                thf
            except:
                thf = tmp[ii]
            else:
                thf = thf + tmp[ii]
    else:
        thf = ReadSelectRegionCheckUnits(thffile, thfname, 'heat flux', box=thfbox, time_bnds=timebounds,
                                         frequency=frequency)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, thf = CheckTime(sst, thf, frequency=frequency, minimum_length=mintimesteps, metric_name='EnsoAlpha')

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the linear regression for all points, for SSTA >=0 and
    # for SSTA<=0
    alpha, alphaPos, alphaNeg = TimeAnomaliesLinearRegressionAndNonlinearity(thf, sst, return_stderr=True)

    # Create output
    alphaMetric = {
        'name': Name, 'value': alpha[0], 'value_error': alpha[1], 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
        'nonlinearity': alphaNeg[0] - alphaPos[0], 'nonlinearity_error': alphaNeg[1] + alphaPos[1],
    }
    return alphaMetric


def EnsoAmpl(sstfile, sstname, sstbox, timebounds=None, frequency=None, mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return amplMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'ENSO amplitude'
    Units = 'C'
    Method = 'Standard deviation of ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)

    # checks if the time-period fulfills the minimum length criterion
    if mintimesteps is not None:
        if len(sst) < mintimesteps:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoAmpl', len(sst), mintimesteps, INSPECTstack())

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the standard deviation
    sstStd = TimeAnomaliesStd(sst)

    # Standard Error of the Standard Deviation (function of nyears)
    sstStdErr = sstStd / NUMPYsqrt(yearN)

    # Create output
    amplMetric = {
        'name': Name, 'value': sstStd, 'value_error': sstStdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
    }
    return amplMetric


def EnsoMu(sstfile, tauxfile, sstname, tauxname, sstbox, tauxbox, timebounds=None, frequency=None, mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return muMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, nonlinearity,
        nonlinearity_error

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'Bjerknes feedback (mu)'
    Units = '10-3 N/m2/C'
    Method = 'Regression of ' + tauxbox + ' tauxA over ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'

    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, time_bnds=timebounds,
                                     frequency=frequency)
    taux = ReadSelectRegionCheckUnits(tauxfile, tauxname, 'heat flux', box=tauxbox, time_bnds=timebounds,
                                      frequency=frequency)

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, taux = CheckTime(sst, taux, frequency=frequency, minimum_length=mintimesteps, metric_name='EnsoMu')

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Averages spatially, removes the annual cycle, computes the linear regression for all points, for SSTA >=0 and
    # for SSTA<=0
    mu, muPos, muNeg = TimeAnomaliesLinearRegressionAndNonlinearity(taux, sst, return_stderr=True)

    # Change units
    mu = [Multiply(mu[0], 1000.), Multiply(mu[1], 1000.)]
    muPos = [Multiply(muPos[0], 1000.), Multiply(muPos[1], 1000.)]
    muNeg = [Multiply(muNeg[0], 1000.), Multiply(muNeg[1], 1000.)]

    # Create output
    muMetric = {
        'name': Name, 'value': mu[0], 'value_error': mu[1], 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': muNeg[0] - muPos[0],
        'nonlinearity_error': muNeg[1] + muPos[1],
    }
    return muMetric


def EnsoRMSE(sstfilemodel, sstnamemodel, sstfileobs, sstnameobs, box, timeboundsmodel=None, timeboundsobs=None,
             centered_rmse=0, frequency=None, mintimesteps=None, regrid_model_on_obs=True):
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
    :param timeboundsmodel: tuple, optional
        tuple of the first and last dates to extract from the modeled SST file (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param timeboundsobs: tuple, optional
        tuple of the first and last dates to extract from the observed SST file (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param regrid_model_on_obs: boolean, optional
        default value = True, regrids the modeled SST on the observed SST grid
        True if you want to regrid, if you don't want it pass anything but true
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

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
    # Define metric attributes
    Name = 'ENSO RMSE'
    Units = 'C'
    Method = 'Spatial root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regriding and rms (uncentered and biased) calculation'

    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=box,
                                           time_bnds=timeboundsmodel, frequency=frequency)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=box,
                                         time_bnds=timeboundsobs, frequency=frequency)

    # checks if the time-period fulfills the minimum length criterion
    if mintimesteps is not None:
        if len(sst_model) < mintimesteps:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoRMSE: the modeled time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mintimesteps) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mintimesteps:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoRMSE: the observed time-period is too short: " + str(len(sst_model))
                            + " (minimum time-period: " + str(mintimesteps) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Averages temporally, regrids sst_model on sst_obs grid (if applicable) and computes the root mean square
    # difference
    sstRmse = SpatialRms(sst_model, sst_obs, centered=centered_rmse, regrid_tab_on_ref=regrid_model_on_obs)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': frequency,
        'time_period_model':actualtimeboundsmodel, 'time_period_observations':actualtimeboundsobs, 'ref': Ref,
    }
    return rmseMetric


def EnsoSeasonality(sstfile, sstname, box, timebounds=None, frequency=None, mintimesteps=None):
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
    :param timebounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., timebounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param mintimesteps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360

    Output:
    ------
    :return SeaMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # Define metric attributes
    Name = 'ENSO Seasonality'
    Units = ''
    Method = 'Ratio between NDJ and MAM standard deviation ' + box + ' sstA'
    Ref = 'Using CDAT std dev calculation'
    # Read file and select the right region
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=box, time_bnds=timebounds,
                                     frequency=frequency)

    # checks if the time-period fulfills the minimum length criterion
    if mintimesteps is not None:
        if len(sst) < mintimesteps:
            EnsoErrorsWarnings.TooShortTimePeriod('EnsoSeasonality', len(sst), mintimesteps, INSPECTstack())

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    # Seasonal ans Spatial average
    sst_NDJ = SeasonalMean(sst, 'NDJ', compute_anom=True)
    sst_MAM = SeasonalMean(sst, 'MAM', compute_anom=True)

    # Compute std dev and ratio
    sst_NDJ_std = TimeAnomaliesStd(sst_NDJ)
    sst_MAM_std = TimeAnomaliesStd(sst_MAM)
    ratioStd = float(sst_NDJ_std / sst_MAM_std)

    # Standard Error of the Standard Deviation (function of nyears)
    sst_NDJ_std_err = sst_NDJ_std / NUMPYsqrt(yearN - 1)
    sst_MAM_std_err = sst_MAM_std / NUMPYsqrt(yearN)

    # The error (dy) on ratio ('y = x/z'): dy = (z*dx + x*dz) / z2
    ratio_std_err = float((sst_MAM_std * sst_NDJ_std_err + sst_NDJ_std * sst_MAM_std_err) / NUMPYsquare(sst_MAM_std))

    # Create output
    seaMetric = {
        'name': Name, 'value': ratioStd, 'value_error': ratio_std_err, 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': frequency, 'time_period': actualtimebounds, 'ref': Ref,
    }
    return seaMetric
#----------------------------------------------------------------------------------------------------------------------#










#----------------------------------------------------------------------------------------------------------------------#
#
# This function has the metrics collection, model name, model file names and model variable names as inputs and metric
# as output
#


dict_oneVar_modelAndObs = {'EnsoRMSE': EnsoRMSE,
                           }

dict_oneVar = {'EnsoAmpl': EnsoAmpl, 'EnsoSeasonality': EnsoSeasonality,
               }

dict_twoVar = {'EnsoAlpha': EnsoAlpha, 'EnsoAlphaLhf': EnsoAlphaLhf, 'EnsoAlphaLwr': EnsoAlphaLwr,
               'EnsoAlphaSwr': EnsoAlphaSwr, 'EnsoMu': EnsoMu,
               }


def ComputeMetric(MetricCollection, metric, modelName, modelFile1, modelVarName1, obsName, obsFile1, obsVarName1,
                  regionVar1='', regionVar2='', modelFile2='', modelVarName2='', obsFile2='', obsVarName2=''):
    '''
    The ComputeMetric() function computes the given metric for the given model and observations

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    - ???		- ???

    Output:
    ------
    - ???

    Method:
    -------
    - ???

    '''
    MC_long_name = defCollection(MetricCollection)['long_name']
    CMIP_experiment = defCollection(MetricCollection)['CMIP_experiment']
    dict_regions = defCollection(MetricCollection)['dict_of_regions'][metric]
    frequency = defCollection(MetricCollection)['metrics_parameters']['frequency']
    mintimesteps = defCollection(MetricCollection)['metrics_parameters']['mintimesteps']
    observed_period = defCollection(MetricCollection)['metrics_parameters']['observed_period']
    modeled_period = defCollection(MetricCollection)['metrics_parameters']['modeled_period']
    var_names = defCollection(MetricCollection)['metrics_list'][metric]['variables']
    if not regionVar1:
        regionVar1 = dict_regions[var_names[0]]
    if not regionVar2:
        try: regionVar2 = dict_regions[var_names[2]]
        except: pass

    if metric in dict_oneVar_modelAndObs.keys():
        metric_mod_obs = dict_oneVar_modelAndObs[metric](modelFile1, modelVarName1, obsFile1, obsVarName1, regionVar1,
                                                         timeboundsmodel=modeled_period, timeboundsobs=observed_period,
                                                         centered_rmse=0, frequency=frequency,
                                                         mintimesteps=mintimesteps, regrid_model_on_obs=True)
        metric_val = {
            'MC_long_name': MC_long_name, 'M_name': metric_mod_obs['name'], 'metric': metric_mod_obs['value'],
            'metric_error': metric_mod_obs['value_error'],
            'comment': "The metric is the statistical value between the model and the observations",
            'CMIP_experiment': CMIP_experiment, 'model': modelName, 'nyears_model': metric_mod_obs['nyears_model'],
            'time_period_model': metric_mod_obs['time_period_model'], 'observations': obsName,
            'nyears_observations': metric_mod_obs['nyears_observations'],
            'time_period_observations': metric_mod_obs['time_period_observations'], 'units': metric_mod_obs['units'],
            'time_frequency': frequency, 'method': metric_mod_obs['method'], 'ref': metric_mod_obs['Ref'],
        }
    else:
        if metric in dict_oneVar.keys():
            metric_mod = dict_oneVar[metric](modelFile1, modelVarName1, regionVar1)
            metric_obs = dict_oneVar[metric](obsFile1, obsVarName1, regionVar1)
        elif metric in dict_twoVar.keys():
            metric_mod = dict_twoVar[metric](modelFile1, modelFile2, modelVarName1, modelVarName2, regionVar1,
                                             regionVar2)
            metric_obs = dict_twoVar[metric](obsFile1, obsFile2, obsVarName1, obsVarName2, regionVar1, regionVar2)
        v1, v2, err1, err2 = metric_mod['value'], metric_obs['value'], metric_mod['value_error'], metric_obs[
            'value_error']
        val1 = metric_mod['value'] / metric_obs['value']
        val1_err = (v1 * err2 + v2 * err1) / NUMPYsquare(v2)
        metric_val = {
            'MC_long_name': MC_long_name, 'M_name': metric_mod['name'], 'metric': val1, 'metric_error': val1_err,
            'comment': "The metric is the ratio value_model / value_observations", 'CMIP_experiment': CMIP_experiment,
            'model': modelName, 'nyears_model': metric['nyears_model'], 'time_period_model': metric_mod['time_period'],
            'value_model': v1, 'value_error_model': err1, 'observations': obsName,
            'nyears_observations':metric['nyears_obs'], 'time_period_observations': metric_obs['time_period'],
            'value_observations': v2, 'value_error_observations': err2, 'units': metric_mod['units'],
            'time_frequency': frequency, 'method': metric_mod['method'], 'ref': metric_mod['ref'],
        }
        try:
            metric_mod['nonlinearity']
        except:
            pass
        else:
            v1, v2 = metric_mod['nonlinearity'], metric_obs['nonlinearity']
            err1, err2 = metric_mod['nonlinearity_error'], metric_obs['nonlinearity_error']
            val2 = metric_mod['value'] / metric_obs['value']
            val2_err = (v1 * err2 + v2 * err1) / NUMPYsquare(v2)
            metric_val['nonlinearity_model'], metric_val['nonlinearity_error_model'] = v1, err1
            metric_val['nonlinearity_observations'], metric_val['nonlinearity_error_observations'] = v2, err2
            metric_val['metric_nonlinearity'], metric_val['metric_nonlinearity_error'] = val2, val2_err
    return metric_val
#----------------------------------------------------------------------------------------------------------------------#
