# -*- coding:UTF-8 -*-
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import sign as NUMPYsign
from numpy import sqrt as NUMPYsqrt
from numpy import square as NUMPYsquare

# ENSO_metrics package functions:
from EnsoCollectionsLib import ReferenceRegions
import EnsoErrorsWarnings
from EnsoToolsLib import percentage_val_eastward
from EnsoUvcdatToolsLib import ArrayOnes, ArrayToList, ApplyLandmask, ApplyLandmaskToArea, AverageMeridional,\
    AverageZonal, CheckTime, Composite, Composite_ev_by_ev, ComputeInterannualAnomalies, ComputePDF, Correlation,\
    DetectEvents, DurationAllEvent, FindXYMinMaxInTs, LinearRegressionAndNonlinearity, LinearRegressionTsAgainstMap,\
    MinMax, MyDerive, PreProcessTS, Read_data_mask_area, ReadAreaSelectRegion, ReadLandmaskSelectRegion,\
    ReadSelectRegionCheckUnits, Regrid, RmsAxis, RmsHorizontal, RmsMeridional, RmsZonal, SaveNetcdf,\
    SeasonalMean, SkewMonthly, SkewnessTemporal, Std, StdMonthly, TimeBounds, TwoVarRegrid
from KeyArgLib import DefaultArgValues


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO metrics
# These functions have file names and variable names as inputs and metric as output
#
def BiasSstRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False, netcdf_name='',
                **kwargs):
    """
    The BiasSstRmse() function computes the SST spatial root mean square error (RMSE) in a 'box' (usually the tropical
    Pacific)

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('tropical_pacific') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If you want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sst RMSE'
    Units = 'C'
    Method = 'Spatial root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasSstRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        sstRmse, sstRmseErr = None, None
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Computes the root mean square difference
        sstRmse = RmsHorizontal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        sstRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}

        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: sstRmse,
                     'metric_value_error_' + dataset2: sstRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sst_map__' + dataset1, var2=sst_obs,
                       var2_attributes=dict2, var2_name='sst_map__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstRmse), 'line2': 'metric value_error: ' + str(sstRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    rmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': sstRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return rmseMetric


def BiasSstLatRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                   sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                   centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False, netcdf_name='',
                   **kwargs):
    """
    The BiasSstLatRmse() function computes the SST meridional (latitude) root mean square error (RMSE) in a 'box'
    (usually 'nino3.3_LatExt')

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3.3_LatExt') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sst Meridional RMSE'
    Units = 'C'
    Method = 'Meridional root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasSstLatRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        sstRmse, sstRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Zonal average
        sst_mod = AverageZonal(sst_mod)
        sst_obs = AverageZonal(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageZonal', 15, **dict_debug)

        # Computes the root mean square difference
        sstRmse = RmsMeridional(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        sstRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                          'axis': list(sst_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafilemod, name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafileobs, name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time',
                                               compute_anom=False, **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '',
                                                          region='equatorial_pacific_LatExt2', **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: sstRmse,
                     'metric_value_error_' + dataset2: sstRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sst_lat__' + dataset1, var2=sst_obs,
                       var2_attributes=dict2, var2_name='sst_lat__' + dataset2, var3=map_mod, var3_attributes=dict3,
                       var3_name='sst_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='sst_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstRmse), 'line2': 'metric value_error: ' + str(sstRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': sstRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LatRmseMetric


def BiasSstLonRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod,
                   sstlandmasknamemod, sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs,
                   sstlandmasknameobs, box, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False,
                   netcdf=False, netcdf_name='', **kwargs):
    """
    The BiasSstLonRmse() function computes the SST zonal (longitude) root mean square error (RMSE) in a 'box'
    (usually the Equatorial Pacific)

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sst Zonal RMSE'
    Units = 'C'
    Method = 'Zonal root mean square error of ' + box + ' sst'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasSstLonRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        sstRmse, sstRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                       **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        sst_mod = AverageMeridional(sst_mod)
        sst_obs = AverageMeridional(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # Computes the root mean square difference
        sstRmse = RmsZonal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        sstRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                          'axis': list(sst_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafilemod, name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafileobs, name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time',
                                             compute_anom=False, **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '',
                                                        region='equatorial_pacific_LatExt2', **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: sstRmse,
                     'metric_value_error_' + dataset2: sstRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sst_lon__' + dataset1, var2=sst_obs,
                       var2_attributes=dict2, var2_name='sst_lon__' + dataset2, var3=map_mod, var3_attributes=dict3,
                       var3_name='sst_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='sst_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstRmse), 'line2': 'metric value_error: ' + str(sstRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': sstRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric


def BiasSstSkLonRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                     sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs,
                     box, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                     netcdf_name='', **kwargs):
    """
    The BiasSstSkLonRmse() function computes the SST zonal (longitude) skewness and then its root mean square error
    (RMSE) in a 'box' (usually the Equatorial Pacific)

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sstA Skewness Zonal RMSE'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Method = 'Zonal root mean square error of ' + box + ' sstA skewness'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasSstSkLonRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        skeRmse, skeRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average=False, compute_anom=True,
                                       **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=True,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        sst_mod = AverageMeridional(sst_mod)
        sst_obs = AverageMeridional(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # skewness
        sst_mod = SkewnessTemporal(sst_mod)
        sst_obs = SkewnessTemporal(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after SkewnessTemporal', 15, **dict_debug)

        # Computes the root mean square difference
        skeRmse = RmsZonal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        skeRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                          'axis': list(sst_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafilemod, name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafileobs, name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average=False, compute_anom=True,
                                             **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average=False, compute_anom=True,
                                             **kwargs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                              'line1': '(model) minmax' + str(MinMax(map_mod)),
                              'line2': '(obs) minmax' + str(MinMax(map_obs)),
                              'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS: netcdf', 15, **dict_debug)
            del mod_areacell, obs_areacell
            # skewness
            ske_map_mod = SkewnessTemporal(map_mod)
            ske_map_obs = SkewnessTemporal(map_obs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in ske_map_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in ske_map_obs.getAxisList()]),
                              'line1': '(model) minmax' + str(MinMax(ske_map_mod)),
                              'line2': '(obs) minmax' + str(MinMax(ske_map_obs)),
                              'shape1': '(model) ' + str(ske_map_mod.shape),
                              'shape2': '(obs) ' + str(ske_map_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SkewnessTemporal: netcdf', 15, **dict_debug)
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                ske_map_mod, ske_map_obs, unneeded = TwoVarRegrid(ske_map_mod, ske_map_obs, '',
                                                        region='equatorial_pacific_LatExt2', **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in ske_map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in ske_map_obs.getAxisList()]),
                                  'line1': '(model) minmax' + str(MinMax(ske_map_mod)),
                                  'line2': '(obs) minmax' + str(MinMax(ske_map_obs)),
                                  'shape1': '(model) ' + str(ske_map_mod.shape),
                                  'shape2': '(obs) ' + str(ske_map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: skeRmse,
                     'metric_value_error_' + dataset2: skeRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sstSke_lon__' + dataset1,
                       var2=sst_obs, var2_attributes=dict2, var2_name='sstSke_lon__' + dataset2, var3=ske_map_mod,
                       var3_attributes=dict3, var3_name='sstSke_map__' + dataset1, var4=ske_map_obs,
                       var4_attributes=dict4, var4_name='sstSke_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(skeRmse), 'line2': 'metric value_error: ' + str(skeRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': skeRmse, 'value_error': skeRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimebounds_mod, 'time_period_observations':actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric


def BiasPrRmse(prfilemod, prnamemod, prareafilemod, prareanamemod, prlandmaskfilemod, prlandmasknamemod,
               prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs, prlandmasknameobs, box,
               centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False, netcdf_name='',
               **kwargs):
    """
    The BiasPrRmse() function computes the PR spatial root mean square error (RMSE) in a 'box' (usually the tropical
    Pacific)

    Inputs:
    ------
    :param prfilemod: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemod: string
        name of PR variable (pr, precip) in 'prfilemod'
    :param prareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for PR
    :param prareanamemod: string
        name of areacell variable (areacella, areacello) in 'prareafilemod'
    :param prlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for PR
    :param prlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfilemod'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param prareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for PR
    :param prareanameobs: string
        name of areacell variable (areacella, areacello) in 'prareafileobs'
    :param prlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for PR
    :param prlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfileobs'
    :param box: string
        name of box ('tropical_pacific') for PR
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If you want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return rmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'pr RMSE'
    Units = 'mm/day'
    Method = 'Spatial root mean square error of ' + box + ' pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasPrRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    pr_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, box, file_area=prareafilemod,
                            name_area=prareanamemod, file_mask=prlandmaskfilemod, name_mask=prlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    pr_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, box, file_area=prareafileobs,
                            name_area=prareanameobs, file_mask=prlandmaskfileobs, name_mask=prlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = pr_mod.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(pr_mod)
    actualtimebounds_obs = TimeBounds(pr_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        prRmse, prRmseErr = None, None
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        pr_mod, Method = PreProcessTS(pr_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            pr_mod, pr_obs, Method = TwoVarRegrid(pr_mod, pr_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                              'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Computes the root mean square difference
        prRmse = RmsHorizontal(pr_mod, pr_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        prRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}

        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: prRmse,
                     'metric_value_error_' + dataset2: prRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=pr_mod, var1_attributes=dict1, var1_name='pr_map__' + dataset1, var2=pr_obs,
                       var2_attributes=dict2, var2_name='pr_map__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(prRmse), 'line2': 'metric value_error: ' + str(prRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    rmseMetric = {
        'name': Name, 'value': prRmse, 'value_error': prRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return rmseMetric


def BiasPrLatRmse(prfilemod, prnamemod, prareafilemod, prareanamemod, prlandmaskfilemod, prlandmasknamemod,
                  prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs, prlandmasknameobs, box,
                  centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False, netcdf_name='',
                  **kwargs):
    """
    The BiasPrLatRmse() function computes the PR meridional (latitude) root mean square error (RMSE) in a 'box'
    (usually 'nino3.3_LatExt')

    Inputs:
    ------
    :param prfilemod: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemod: string
        name of PR variable (pr, precip) in 'prfilemod'
    :param prareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for PR
    :param prareanamemod: string
        name of areacell variable (areacella, areacello) in 'prareafilemod'
    :param prlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for PR
    :param prlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfilemod'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param prareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for PR
    :param prareanameobs: string
        name of areacell variable (areacella, areacello) in 'prareafileobs'
    :param prlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for PR
    :param prlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfileobs'
    :param box: string
        name of box ('nino3.3_LatExt') for PR
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LatRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'pr Meridional RMSE'
    Units = 'mm/day'
    Method = 'Meridional root mean square error of ' + box + ' pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasPrLatRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    pr_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, box, file_area=prareafilemod,
                            name_area=prareanamemod, file_mask=prlandmaskfilemod, name_mask=prlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    pr_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, box, file_area=prareafileobs,
                            name_area=prareanameobs, file_mask=prlandmaskfileobs, name_mask=prlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = pr_mod.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(pr_mod)
    actualtimebounds_obs = TimeBounds(pr_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        prRmse, prRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        pr_mod, Method = PreProcessTS(pr_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            pr_mod, pr_obs, Method = TwoVarRegrid(pr_mod, pr_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                              'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Zonal average
        pr_mod = AverageZonal(pr_mod)
        pr_obs = AverageZonal(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageZonal', 15, **dict_debug)

        # Computes the root mean square difference
        prRmse = RmsMeridional(pr_mod, pr_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        prRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(pr_mod), 'observations': ArrayToList(pr_obs),
                          'axis': list(pr_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, 'equatorial_pacific_LatExt2',
                                    file_area=prareafilemod, name_area=prareanamemod, file_mask=prlandmaskfilemod,
                                    name_mask=prlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, 'equatorial_pacific_LatExt2',
                                    file_area=prareafileobs, name_area=prareanameobs, file_mask=prlandmaskfileobs,
                                    name_mask=prlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time',
                                               compute_anom=False, **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2',
                                                          **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            # change units
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: prRmse,
                     'metric_value_error_' + dataset2: prRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=pr_mod, var1_attributes=dict1, var1_name='pr_lat__' + dataset1, var2=pr_obs,
                       var2_attributes=dict2, var2_name='pr_lat__' + dataset2, var3=map_mod, var3_attributes=dict3,
                       var3_name='pr_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='pr_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(prRmse), 'line2': 'metric value_error: ' + str(prRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': prRmse, 'value_error': prRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LatRmseMetric


def BiasPrLonRmse(prfilemod, prnamemod, prareafilemod, prareanamemod, prlandmaskfilemod, prlandmasknamemod, prfileobs,
                  prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs, prlandmasknameobs, box, centered_rmse=0,
                  biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The BiasPrLonRmse() function computes the PR zonal (longitude) root mean square error (RMSE) in a 'box'
    (usually the Equatorial Pacific)

    Inputs:
    ------
    :param prfilemod: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemod: string
        name of PR variable (pr, precip) in 'prfilemod'
    :param prareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for PR
    :param prareanamemod: string
        name of areacell variable (areacella, areacello) in 'prareafilemod'
    :param prlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for PR
    :param prlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfilemod'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param prareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for PR
    :param prareanameobs: string
        name of areacell variable (areacella, areacello) in 'prareafileobs'
    :param prlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for PR
    :param prlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for PR
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LonRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'pr Zonal RMSE'
    Units = 'mm/day'
    Method = 'Zonal root mean square error of ' + box + ' pr'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasPrLonRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    pr_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, box, file_area=prareafilemod,
                            name_area=prareanamemod, file_mask=prlandmaskfilemod, name_mask=prlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    pr_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, box, file_area=prareafileobs,
                            name_area=prareanameobs, file_mask=prlandmaskfileobs, name_mask=prlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = pr_mod.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(pr_mod)
    actualtimebounds_obs = TimeBounds(pr_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        prRmse, prRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        pr_mod, Method = PreProcessTS(pr_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            pr_mod, pr_obs, Method = TwoVarRegrid(pr_mod, pr_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                              'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        pr_mod = AverageMeridional(pr_mod)
        pr_obs = AverageMeridional(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # Computes the root mean square difference
        prRmse = RmsZonal(pr_mod, pr_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        prRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(pr_mod), 'observations': ArrayToList(pr_obs),
                          'axis': list(pr_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, 'equatorial_pacific_LatExt2',
                                    file_area=prareafilemod, name_area=prareanamemod, file_mask=prlandmaskfilemod,
                                    name_mask=prlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, 'equatorial_pacific_LatExt2',
                                    file_area=prareafileobs, name_area=prareanameobs, file_mask=prlandmaskfileobs,
                                    name_mask=prlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time',
                                             compute_anom=False, **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2',
                                                          **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: prRmse,
                     'metric_value_error_' + dataset2: prRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=pr_mod, var1_attributes=dict1, var1_name='pr_lon__' + dataset1, var2=pr_obs,
                       var2_attributes=dict2, var2_name='pr_lon__' + dataset2, var3=map_mod, var3_attributes=dict3,
                       var3_name='pr_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='pr_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(prRmse), 'line2': 'metric value_error: ' + str(prRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': prRmse, 'value_error': prRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric


def BiasTauxRmse(tauxfilemod, tauxnamemod, tauxareafilemod, tauxareanamemod, tauxlandmaskfilemod, tauxlandmasknamemod,
                 tauxfileobs, tauxnameobs, tauxareafileobs, tauxareanameobs, tauxlandmaskfileobs, tauxlandmasknameobs,
                 box, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                 netcdf_name='', **kwargs):
    """
    The BiasTauxRmse() function computes the TAUX spatial root mean square error (RMSE) in a 'box' (usually the tropical
    Pacific)

    Inputs:
    ------
    :param tauxfilemod: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemod: string
        name of TAUX variable (taux, tauu) in 'tauxfilemod'
    :param tauxareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for TAUX
    :param tauxareanamemod: string
        name of areacell variable (areacella, areacello) in 'tauxareafilemod'
    :param tauxlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for TAUX
    :param tauxlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfilemod'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param tauxareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for TAUX
    :param tauxareanameobs: string
        name of areacell variable (areacella, areacello) in 'tauxareafileobs'
    :param tauxlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for TAUX
    :param tauxlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfileobs'
    :param box: string
        name of box ('tropical_pacific') for TAUX
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If you want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return rmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'taux RMSE'
    Units = '1e-3 N/m2'
    Method = 'Spatial root mean square error of ' + box + ' taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasTauxRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    taux_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(tauxfilemod, tauxnamemod, 'wind stress', metric, box, file_area=tauxareafilemod,
                            name_area=tauxareanamemod, file_mask=tauxlandmaskfilemod, name_mask=tauxlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    taux_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(tauxfileobs, tauxnameobs, 'wind stress', metric, box, file_area=tauxareafileobs,
                            name_area=tauxareanameobs, file_mask=tauxlandmaskfileobs, name_mask=tauxlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = taux_mod.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(taux_mod)
    actualtimebounds_obs = TimeBounds(taux_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        tauxRmse, tauxRmseErr = None, None
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        taux_mod, Method = PreProcessTS(taux_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                          'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            taux_mod, taux_obs, Method = TwoVarRegrid(taux_mod, taux_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                              'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # change units
        taux_mod = taux_mod * 1e3
        taux_obs = taux_obs * 1e3

        # Computes the root mean square difference
        tauxRmse = RmsHorizontal(taux_mod, taux_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        tauxRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': None, 'observations': None, 'axisLat': None, 'axisLon': None}

        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: tauxRmse,
                     'metric_value_error_' + dataset2: tauxRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=taux_mod, var1_attributes=dict1, var1_name='taux_map__' + dataset1,
                       var2=taux_obs, var2_attributes=dict2, var2_name='taux_map__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(tauxRmse), 'line2': 'metric value_error: ' + str(tauxRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    rmseMetric = {
        'name': Name, 'value': tauxRmse, 'value_error': tauxRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return rmseMetric


def BiasTauxLatRmse(tauxfilemod, tauxnamemod, tauxareafilemod, tauxareanamemod, tauxlandmaskfilemod,
                    tauxlandmasknamemod, tauxfileobs, tauxnameobs, tauxareafileobs, tauxareanameobs,
                    tauxlandmaskfileobs, tauxlandmasknameobs, box, centered_rmse=0, biased_rmse=1, dataset1='',
                    dataset2='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The BiasTauxLatRmse() function computes the TAUX meridional (latitude) root mean square error (RMSE) in a 'box'
    (usually 'nino3.3_LatExt')

    Inputs:
    ------
    :param tauxfilemod: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemod: string
        name of TAUX variable (taux, tauu) in 'tauxfilemod'
    :param tauxareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for TAUX
    :param tauxareanamemod: string
        name of areacell variable (areacella, areacello) in 'tauxareafilemod'
    :param tauxlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for TAUX
    :param tauxlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfilemod'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param tauxareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for TAUX
    :param tauxareanameobs: string
        name of areacell variable (areacella, areacello) in 'tauxareafileobs'
    :param tauxlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for TAUX
    :param tauxlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfileobs'
    :param box: string
        name of box ('nino3.3_LatExt') for TAUX
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LatRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'taux Meridional RMSE'
    Units = '1e-3 N/m2'
    Method = 'Meridional root mean square error of ' + box + ' taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasTauxLatRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    taux_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(tauxfilemod, tauxnamemod, 'wind stress', metric, box, file_area=tauxareafilemod,
                            name_area=tauxareanamemod, file_mask=tauxlandmaskfilemod, name_mask=tauxlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    taux_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(tauxfileobs, tauxnameobs, 'wind stress', metric, box, file_area=tauxareafileobs,
                            name_area=tauxareanameobs, file_mask=tauxlandmaskfileobs, name_mask=tauxlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = taux_mod.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(taux_mod)
    actualtimebounds_obs = TimeBounds(taux_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        tauxRmse, tauxRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        taux_mod, Method = PreProcessTS(taux_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                          'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            taux_mod, taux_obs, Method = TwoVarRegrid(taux_mod, taux_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                              'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Zonal average
        taux_mod = AverageZonal(taux_mod) * 1e3
        taux_obs = AverageZonal(taux_obs) * 1e3
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                          'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageZonal', 15, **dict_debug)

        # Computes the root mean square difference
        tauxRmse = RmsMeridional(taux_mod, taux_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        tauxRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(taux_mod), 'observations': ArrayToList(taux_obs),
                          'axis': list(taux_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(tauxfilemod, tauxnamemod, 'wind stress', metric, 'equatorial_pacific_LatExt2',
                                    file_area=tauxareafilemod, name_area=tauxareanamemod, file_mask=tauxlandmaskfilemod,
                                    name_mask=tauxlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(tauxfileobs, tauxnameobs, 'wind stress', metric, 'equatorial_pacific_LatExt2',
                                    file_area=tauxareafileobs, name_area=tauxareanameobs, file_mask=tauxlandmaskfileobs,
                                    name_mask=tauxlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time', compute_anom=False,
                                             **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2',
                                                          **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            # change units
            map_mod = map_mod * 1e3
            map_obs = map_obs * 1e3
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: tauxRmse,
                     'metric_value_error_' + dataset2: tauxRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=taux_mod, var1_attributes=dict1, var1_name='taux_lat__' + dataset1,
                       var2=taux_obs, var2_attributes=dict2, var2_name='taux_lat__' + dataset2, var3=map_mod,
                       var3_attributes=dict3, var3_name='taux_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='taux_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(tauxRmse), 'line2': 'metric value_error: ' + str(tauxRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': tauxRmse, 'value_error': tauxRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LatRmseMetric


def BiasTauxLonRmse(tauxfilemod, tauxnamemod, tauxareafilemod, tauxareanamemod, tauxlandmaskfilemod,
                    tauxlandmasknamemod, tauxfileobs, tauxnameobs, tauxareafileobs, tauxareanameobs,
                    tauxlandmaskfileobs, tauxlandmasknameobs, box, centered_rmse=0, biased_rmse=1, dataset1='',
                    dataset2='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The BiasTauxLonRmse() function computes the TAUX zonal (longitude) root mean square error (RMSE) in a 'box'
    (usually the Equatorial Pacific)

    Inputs:
    ------
    :param tauxfilemod: string
        path_to/filename of the file (NetCDF) of the modeled TAUX
    :param tauxnamemod: string
        name of TAUX variable (taux, tauu) in 'tauxfilemod'
    :param tauxareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for TAUX
    :param tauxareanamemod: string
        name of areacell variable (areacella, areacello) in 'tauxareafilemod'
    :param tauxlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for TAUX
    :param tauxlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfilemod'
    :param tauxfileobs: string
        path_to/filename of the file (NetCDF) of the observed TAUX
    :param tauxnameobs: string
        name of TAUX variable (taux, tauu) in 'tauxfileobs'
    :param tauxareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for TAUX
    :param tauxareanameobs: string
        name of areacell variable (areacella, areacello) in 'tauxareafileobs'
    :param tauxlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for TAUX
    :param tauxlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for TAUX
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed TAUX file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LonRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'taux Zonal RMSE'
    Units = '1e-3 N/m2'
    Method = 'Zonal root mean square error of ' + box + ' taux'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = "BiasTauxLonRmse"

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    taux_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(tauxfilemod, tauxnamemod, 'wind stress', metric, box, file_area=tauxareafilemod,
                            name_area=tauxareanamemod, file_mask=tauxlandmaskfilemod, name_mask=tauxlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    taux_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(tauxfileobs, tauxnameobs, 'wind stress', metric, box, file_area=tauxareafileobs,
                            name_area=tauxareanameobs, file_mask=tauxlandmaskfileobs, name_mask=tauxlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = taux_mod.shape[0] / 12
    yearN_obs = taux_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(taux_mod)
    actualtimebounds_obs = TimeBounds(taux_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        tauxRmse, tauxRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        taux_mod, Method = PreProcessTS(taux_mod, Method, areacell=mod_areacell, average='time', compute_anom=False,
                                         **kwargs)
        taux_obs, unneeded = PreProcessTS(taux_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                          'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            taux_mod, taux_obs, Method = TwoVarRegrid(taux_mod, taux_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                              'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        taux_mod = AverageMeridional(taux_mod) * 1e3
        taux_obs = AverageMeridional(taux_obs) * 1e3
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in taux_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in taux_obs.getAxisList()]),
                          'shape1': '(model) ' + str(taux_mod.shape), 'shape2': '(obs) ' + str(taux_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # Computes the root mean square difference
        tauxRmse = RmsZonal(taux_mod, taux_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        tauxRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(taux_mod), 'observations': ArrayToList(taux_obs),
                          'axis': list(taux_mod.getAxis(0)[:])}
        if netcdf is True:
            map_mod, mod_areacell, keyerror_mod = \
                Read_data_mask_area(tauxfilemod, tauxnamemod, 'wind stress', metric, 'equatorial_pacific_LatExt2',
                                    file_area=tauxareafilemod, name_area=tauxareanamemod, file_mask=tauxlandmaskfilemod,
                                    name_mask=tauxlandmasknamemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            map_obs, obs_areacell, keyerror_obs = \
                Read_data_mask_area(tauxfileobs, tauxnameobs, 'wind stress', metric, 'equatorial_pacific_LatExt2',
                                    file_area=tauxareafileobs, name_area=tauxareanameobs, file_mask=tauxlandmaskfileobs,
                                    name_mask=tauxlandmasknameobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
            map_mod, unneeded = PreProcessTS(map_mod, '', areacell=mod_areacell, average='time', compute_anom=False,
                                             **kwargs)
            map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average='time', compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            # Regridding
            if isinstance(kwargs['regridding'], dict):
                map_mod, map_obs, unneeded = TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2',
                                                          **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
            # change units
            map_mod = map_mod * 1e3
            map_obs = map_obs * 1e3
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict5 = {'metric_name': Name, 'metric_value_' + dataset2: tauxRmse,
                     'metric_value_error_' + dataset2: tauxRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=taux_mod, var1_attributes=dict1, var1_name='taux_lon__' + dataset1,
                       var2=taux_obs, var2_attributes=dict2, var2_name='taux_lon__' + dataset2, var3=map_mod,
                       var3_attributes=dict3, var3_name='taux_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                       var4_name='taux_map__' + dataset2, global_attributes=dict5)
            del dict1, dict2, dict3, dict4, dict5
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(tauxRmse), 'line2': 'metric value_error: ' + str(tauxRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': tauxRmse, 'value_error': tauxRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric


def EnsoAlphaLhf(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, lhffile, lhfname,
                 lhfareafile, lhfareaname, lhflandmaskfile, lhflandmaskname, lhfbox, dataset='', debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoAlphaLhf() function computes the regression of 'lhfbox' lhfA (latent heat flux anomalies) over 'sstbox' sstA
    (usually the regression of nino3 lhfA over nino3 sstA)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param lhffile: string
        path_to/filename of the file (NetCDF) of LHF
    :param lhfname: string
        name of LHF variable (lhf, hfls) in 'lhffile'
    :param lhfareafile: string
        path_to/filename of the file (NetCDF) of the areacell for LHF
    :param lhfareaname: string
        name of areacell variable (areacella, areacello) in 'lhfareafile'
    :param lhflandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for LHF
    :param lhflandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'lhflandmaskfile'
    :param lhfbox: string
        name of box (nino3') for LHF
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    metric = 'EnsoAlphaLhf'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(lhf) ' + lhffile, 'var1': '(sst) ' + sstname,
                      'var2': '(lhf) ' + lhfname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    lhf = ReadSelectRegionCheckUnits(lhffile, lhfname, 'heat flux', box=lhfbox, **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(lhf) ' + str([ax.id for ax in lhf.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(lhf) ' + str(lhf.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(lhf) ' + str(TimeBounds(lhf))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if lhfareafile:
        lhf_areacell = ReadAreaSelectRegion(lhfareafile, areaname=lhfareaname, box=lhfbox, **kwargs)
    else:
        lhf_areacell = ReadAreaSelectRegion(lhffile, areaname=lhfareaname, box=lhfbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if lhf_areacell is not None:
            dict_debug['axes2'] = '(lhf) ' + str([ax.id for ax in lhf_areacell.getAxisList()])
            dict_debug['shape2'] = '(lhf) ' + str(lhf_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if lhflandmaskfile:
        lhf_landmask = ReadLandmaskSelectRegion(lhflandmaskfile, landmaskname=lhflandmaskname, box=lhfbox, **kwargs)
    else:
        lhf_landmask = ReadLandmaskSelectRegion(lhffile, landmaskname=lhflandmaskname, box=lhfbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if lhf_landmask is not None:
            dict_debug['axes2'] = '(lhf) ' + str([ax.id for ax in lhf_landmask.getAxisList()])
            dict_debug['shape2'] = '(lhf) ' + str(lhf_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if lhf_landmask is not None:
        lhf = ApplyLandmask(lhf, lhf_landmask, maskland=True, maskocean=False)
        if lhf_areacell is None:
            lhf_areacell = ArrayOnes(lhf_landmask, id='areacell')
        lhf_areacell = ApplyLandmaskToArea(lhf_areacell, lhf_landmask, maskland=True, maskocean=False)
    del lhf_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lhf, keyerror = CheckTime(sst, lhf, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        alphaLhf, alphaLhfPos, alphaLhfNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        lhf, unneeded = PreProcessTS(lhf, '', areacell=lhf_areacell, average='horizontal', compute_anom=True,
                                     **kwargs)
        del sst_areacell, lhf_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(lhf) ' + str([ax.id for ax in lhf.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(lhf) ' + str(lhf.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(lhf) ' + str(TimeBounds(lhf))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        alphaLhf, alphaLhfPos, alphaLhfNeg = LinearRegressionAndNonlinearity(lhf, sst, return_stderr=True)

    # Create output
    alphaLhfMetric = {
        'name': Name, 'value': alphaLhf[0], 'value_error': alphaLhf[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaLhfNeg[0] - alphaLhfPos[0],
        'nonlinearity_error': alphaLhfNeg[1] + alphaLhfPos[1], 'keyerror': keyerror,
    }
    return alphaLhfMetric


def EnsoAlphaLwr(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, lwrfile, lwrname,
                 lwrareafile, lwrareaname, lwrlandmaskfile, lwrlandmaskname, lwrbox, dataset='', debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
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
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param lwrfile: string
        path_to/filename of the file (NetCDF) of LWR
    :param lwrname: string
        name of LWR variable (lwr, rlds - rlus) (may be a list of variables) in 'lwrfile'
    :param lwrareafile: string
        path_to/filename of the file (NetCDF) of the areacell for LWR
    :param lwrareaname: string
        name of areacell variable (areacella, areacello) in 'lwrareafile'
    :param lwrlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for LWR
    :param lwrlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'lwrlandmaskfile'
    :param lwrbox: string
        name of box (nino3') for LWR
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    metric = 'EnsoAlphaLwr'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(lwr) ' + lwrfile, 'var1': '(sst) ' + sstname,
                      'var2': '(lwr) ' + lwrname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(lwrfile, basestring):
        dict_var[lwrname] = ReadSelectRegionCheckUnits(lwrfile, lwrname, 'heat flux', box=lwrbox, **kwargs)
    else:
        for ii in range(len(lwrfile)):
            filename, varname = lwrfile[ii], lwrname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=lwrbox, **kwargs)
    lwr = MyDerive(kwargs['project_interpreter_var2'], 'lwr', dict_var)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(lwr) ' + str([ax.id for ax in lwr.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(lwr) ' + str(lwr.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(lwr) ' + str(TimeBounds(lwr))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if lwrareafile:
        if isinstance(lwrareafile, basestring):
            lwr_areacell = ReadAreaSelectRegion(lwrareafile, areaname=lwrareaname, box=lwrbox, **kwargs)
        else:
            for ii in range(len(lwrareafile)):
                lwr_areacell = ReadAreaSelectRegion(lwrareafile[ii], areaname=lwrareaname[ii], box=lwrbox, **kwargs)
                if lwr_areacell is not None:
                    break
    else:
        if isinstance(lwrfile, basestring):
            lwr_areacell = ReadAreaSelectRegion(lwrfile, areaname=lwrareaname, box=lwrbox, **kwargs)
        else:
            for ii in range(len(lwrfile)):
                lwr_areacell = ReadAreaSelectRegion(lwrfile[ii], areaname=lwrareaname[ii], box=lwrbox, **kwargs)
                if lwr_areacell is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if lwr_areacell is not None:
            dict_debug['axes2'] = '(lwr) ' + str([ax.id for ax in lwr_areacell.getAxisList()])
            dict_debug['shape2'] = '(lwr) ' + str(lwr_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if lwrlandmaskfile:
        if isinstance(lwrlandmaskfile, basestring):
            lwr_landmask = ReadAreaSelectRegion(lwrlandmaskfile, areaname=lwrlandmaskname, box=lwrbox, **kwargs)
        else:
            for ii in range(len(lwrlandmaskfile)):
                lwr_landmask = ReadAreaSelectRegion(lwrlandmaskfile[ii], areaname=lwrlandmaskname[ii], box=lwrbox,
                                                    **kwargs)
                if lwr_landmask is not None:
                    break
    else:
        if isinstance(lwrfile, basestring):
            lwr_landmask = ReadAreaSelectRegion(lwrfile, areaname=lwrlandmaskname, box=lwrbox, **kwargs)
        else:
            for ii in range(len(lwrfile)):
                lwr_landmask = ReadAreaSelectRegion(lwrfile[ii], areaname=lwrlandmaskname[ii], box=lwrbox, **kwargs)
                if lwr_landmask is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if lwr_landmask is not None:
            dict_debug['axes2'] = '(lwr) ' + str([ax.id for ax in lwr_landmask.getAxisList()])
            dict_debug['shape2'] = '(lwr) ' + str(lwr_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if lwr_landmask is not None:
        lwr = ApplyLandmask(lwr, lwr_landmask, maskland=True, maskocean=False)
        if lwr_areacell is None:
            lwr_areacell = ArrayOnes(lwr_landmask, id='areacell')
        lwr_areacell = ApplyLandmaskToArea(lwr_areacell, lwr_landmask, maskland=True, maskocean=False)
    del lwr_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, lwr, keyerror = CheckTime(sst, lwr, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        alphaLwr, alphaLwrPos, alphaLwrNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        lwr, unneeded = PreProcessTS(lwr, '', areacell=lwr_areacell, average='horizontal', compute_anom=True, **kwargs)
        del sst_areacell, lwr_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(lwr) ' + str([ax.id for ax in lwr.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(lwr) ' + str(lwr.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(lwr) ' + str(TimeBounds(lwr))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        alphaLwr, alphaLwrPos, alphaLwrNeg = LinearRegressionAndNonlinearity(lwr, sst, return_stderr=True)

    # Create output
    alphaLwrMetric = {
        'name': Name, 'value': alphaLwr[0], 'value_error': alphaLwr[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaLwrNeg[0] - alphaLwrPos[0],
        'nonlinearity_error': alphaLwrNeg[1] + alphaLwrPos[1], 'keyerror': keyerror,
    }
    return alphaLwrMetric


def EnsoAlphaShf(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, shffile, shfname,
                 shfareafile, shfareaname, shflandmaskfile, shflandmaskname, shfbox, dataset='', debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoAlphaShf() function computes the regression of 'shfbox' shfA (sensible heat flux anomalies) over 'sstbox'
    sstA (usually the regression of nino3 shfA over nino3 sstA)

    Author:	Yann Planton : yann.planton@locean-ipsl.upmc.fr
    Co-author:

    Created on Thu Oct  5 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param shffile: string
        path_to/filename of the file (NetCDF) of SHF
    :param shfname: string
        name of SHF variable (shf, hfss) in 'shffile'
    :param shfareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SHF
    :param shfareaname: string
        name of areacell variable (areacella, areacello) in 'shfareafile'
    :param shflandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SHF
    :param shflandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'shflandmaskfile'
    :param shfbox: string
        name of box (nino3') for SHF
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    Name = 'Sensible feedback (alpha_sh)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + shfbox + ' shfA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'
    metric = 'EnsoAlphaShf'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(shf) ' + shffile, 'var1': '(sst) ' + sstname,
                      'var2': '(shf) ' + shfname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    shf = ReadSelectRegionCheckUnits(shffile, shfname, 'heat flux', box=shfbox, **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(shf) ' + str([ax.id for ax in shf.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(shf) ' + str(shf.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(shf) ' + str(TimeBounds(shf))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if shfareafile:
        shf_areacell = ReadAreaSelectRegion(shfareafile, areaname=shfareaname, box=shfbox, **kwargs)
    else:
        shf_areacell = ReadAreaSelectRegion(shffile, areaname=shfareaname, box=shfbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if shf_areacell is not None:
            dict_debug['axes2'] = '(shf) ' + str([ax.id for ax in shf_areacell.getAxisList()])
            dict_debug['shape2'] = '(shf) ' + str(shf_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if shflandmaskfile:
        shf_landmask = ReadLandmaskSelectRegion(shflandmaskfile, landmaskname=shflandmaskname, box=shfbox, **kwargs)
    else:
        shf_landmask = ReadLandmaskSelectRegion(shffile, landmaskname=shflandmaskname, box=shfbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if shf_landmask is not None:
            dict_debug['axes2'] = '(shf) ' + str([ax.id for ax in shf_landmask.getAxisList()])
            dict_debug['shape2'] = '(shf) ' + str(shf_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if shf_landmask is not None:
        shf = ApplyLandmask(shf, shf_landmask, maskland=True, maskocean=False)
        if shf_areacell is None:
            shf_areacell = ArrayOnes(shf_landmask, id='areacell')
        shf_areacell = ApplyLandmaskToArea(shf_areacell, shf_landmask, maskland=True, maskocean=False)
    del shf_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, shf, keyerror = CheckTime(sst, shf, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        alphaShf, alphaShfPos, alphaShfNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        shf, unneeded = PreProcessTS(shf, '', areacell=shf_areacell, average='horizontal', compute_anom=True, **kwargs)
        del sst_areacell, shf_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(shf) ' + str([ax.id for ax in shf.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(shf) ' + str(shf.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(shf) ' + str(TimeBounds(shf))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        alphaShf, alphaShfPos, alphaShfNeg = LinearRegressionAndNonlinearity(shf, sst, return_stderr=True)

    # Create output
    alphaShfMetric = {
        'name': Name, 'value': alphaShf[0], 'value_error': alphaShf[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaShfNeg[0] - alphaShfPos[0],
        'nonlinearity_error': alphaShfNeg[1] + alphaShfPos[1], 'keyerror': keyerror,
    }
    return alphaShfMetric


def EnsoAlphaSwr(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, swrfile, swrname,
                 swrareafile, swrareaname, swrlandmaskfile, swrlandmaskname, swrbox, dataset='', debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
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
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param swrfile: string
        path_to/filename of the file (NetCDF) of SWR
    :param swrname: string
        name of SWR variable (swr, rsds - rsus) (may be a list of variables) in 'swrfile'
    :param swrareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SWR
    :param swrareaname: string
        name of areacell variable (areacella, areacello) in 'swrareafile'
    :param swrlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SWR
    :param swrlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'swrlandmaskfile'
    :param swrbox: string
        name of box (nino3') for SWR
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    Name = 'Shortwave feedback (alpha_swr)'
    Units = 'W/m2/C'
    Method = 'Regression of ' + swrbox + ' swrA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'
    metric = 'EnsoAlphaSwr'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(swr) ' + swrfile, 'var1': '(sst) ' + sstname,
                      'var2': '(swr) ' + swrname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(swrfile, basestring):
        dict_var[swrname] = ReadSelectRegionCheckUnits(swrfile, swrname, 'heat flux', box=swrbox, **kwargs)
    else:
        for ii in range(len(swrfile)):
            filename, varname = swrfile[ii], swrname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=swrbox, **kwargs)
    swr = MyDerive(kwargs['project_interpreter_var2'], 'swr', dict_var)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(swr) ' + str([ax.id for ax in swr.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(swr) ' + str(swr.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(swr) ' + str(TimeBounds(swr))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if swrareafile:
        if isinstance(swrareafile, basestring):
            swr_areacell = ReadAreaSelectRegion(swrareafile, areaname=swrareaname, box=swrbox, **kwargs)
        else:
            for ii in range(len(swrareafile)):
                swr_areacell = ReadAreaSelectRegion(swrareafile[ii], areaname=swrareaname[ii], box=swrbox, **kwargs)
                if swr_areacell is not None:
                    break
    else:
        if isinstance(swrfile, basestring):
            swr_areacell = ReadAreaSelectRegion(swrfile, areaname=swrareaname, box=swrbox, **kwargs)
        else:
            for ii in range(len(swrfile)):
                swr_areacell = ReadAreaSelectRegion(swrfile[ii], areaname=swrareaname[ii], box=swrbox, **kwargs)
                if swr_areacell is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if swr_areacell is not None:
            dict_debug['axes2'] = '(swr) ' + str([ax.id for ax in swr_areacell.getAxisList()])
            dict_debug['shape2'] = '(swr) ' + str(swr_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if swrlandmaskfile:
        if isinstance(swrlandmaskfile, basestring):
            swr_landmask = ReadAreaSelectRegion(swrlandmaskfile, areaname=swrlandmaskname, box=swrbox, **kwargs)
        else:
            for ii in range(len(swrlandmaskfile)):
                swr_landmask = ReadAreaSelectRegion(swrlandmaskfile[ii], areaname=swrlandmaskname[ii], box=swrbox,
                                                    **kwargs)
                if swr_landmask is not None:
                    break
    else:
        if isinstance(swrfile, basestring):
            swr_landmask = ReadAreaSelectRegion(swrfile, areaname=swrlandmaskname, box=swrbox, **kwargs)
        else:
            for ii in range(len(swrfile)):
                swr_landmask = ReadAreaSelectRegion(swrfile[ii], areaname=swrlandmaskname[ii], box=swrbox, **kwargs)
                if swr_landmask is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if swr_landmask is not None:
            dict_debug['axes2'] = '(swr) ' + str([ax.id for ax in swr_landmask.getAxisList()])
            dict_debug['shape2'] = '(swr) ' + str(swr_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if swr_landmask is not None:
        swr = ApplyLandmask(swr, swr_landmask, maskland=True, maskocean=False)
        if swr_areacell is None:
            swr_areacell = ArrayOnes(swr_landmask, id='areacell')
        swr_areacell = ApplyLandmaskToArea(swr_areacell, swr_landmask, maskland=True, maskocean=False)
    del swr_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, swr, keyerror = CheckTime(sst, swr, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        alphaSwr, alphaSwrPos, alphaSwrNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        swr, unneeded = PreProcessTS(swr, '', areacell=swr_areacell, average='horizontal', compute_anom=True, **kwargs)
        del sst_areacell, swr_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(swr) ' + str([ax.id for ax in swr.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(swr) ' + str(swr.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(swr) ' + str(TimeBounds(swr))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        alphaSwr, alphaSwrPos, alphaSwrNeg = LinearRegressionAndNonlinearity(swr, sst, return_stderr=True)

    # Create output
    alphaSwrMetric = {
        'name': Name, 'value': alphaSwr[0], 'value_error': alphaSwr[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaSwrNeg[0] - alphaSwrPos[0],
        'nonlinearity_error': alphaSwrNeg[1] + alphaSwrPos[1], 'keyerror': keyerror,
    }
    return alphaSwrMetric


def EnsoAlphaThf(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, thffile, thfname,
                 thfareafile, thfareaname, thflandmaskfile, thflandmaskname, thfbox, dataset='', debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
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
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param thffile: string
        path_to/filename of the file (NetCDF) of THF
    :param thfname: string
        name of THF variable (thf, netflux, thflx, thf + lwr + lhf + shf) (may be a list of variables) in 'thffile'
    :param thfareafile: string
        path_to/filename of the file (NetCDF) of the areacell for THF
    :param thfareaname: string
        name of areacell variable (areacella, areacello) in 'thfareafile'
    :param thflandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for THF
    :param thflandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'thflandmaskfile'
    :param thfbox: string
        name of box (nino3') for THF
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    metric = 'EnsoAlphaThf'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(thf) ' + thffile, 'var1': '(sst) ' + sstname,
                      'var2': '(thf) ' + thfname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    dict_var = dict()
    if isinstance(thffile, basestring):
        dict_var[thfname] = ReadSelectRegionCheckUnits(thffile, thfname, 'heat flux', box=thfbox, **kwargs)
    else:
        for ii in range(len(thffile)):
            filename, varname = thffile[ii], thfname[ii]
            dict_var[varname] = ReadSelectRegionCheckUnits(filename, varname, 'heat flux', box=thfbox, **kwargs)
    thf = MyDerive(kwargs['project_interpreter_var2'], 'thf', dict_var)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(thf) ' + str([ax.id for ax in thf.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(thf) ' + str(thf.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(thf) ' + str(TimeBounds(thf))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if thfareafile:
        if isinstance(thfareafile, basestring):
            thf_areacell = ReadAreaSelectRegion(thfareafile, areaname=thfareaname, box=thfbox, **kwargs)
        else:
            for ii in range(len(thfareafile)):
                thf_areacell = ReadAreaSelectRegion(thfareafile[ii], areaname=thfareaname[ii], box=thfbox, **kwargs)
                if thf_areacell is not None:
                    break
    else:
        if isinstance(thffile, basestring):
            thf_areacell = ReadAreaSelectRegion(thffile, areaname=thfareaname, box=thfbox, **kwargs)
        else:
            for ii in range(len(thffile)):
                thf_areacell = ReadAreaSelectRegion(thffile[ii], areaname=thfareaname[ii], box=thfbox, **kwargs)
                if thf_areacell is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if thf_areacell is not None:
            dict_debug['axes2'] = '(thf) ' + str([ax.id for ax in thf_areacell.getAxisList()])
            dict_debug['shape2'] = '(thf) ' + str(thf_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if thflandmaskfile:
        if isinstance(thflandmaskfile, basestring):
            thf_landmask = ReadAreaSelectRegion(thflandmaskfile, areaname=thflandmaskname, box=thfbox, **kwargs)
        else:
            for ii in range(len(thflandmaskfile)):
                thf_landmask = ReadAreaSelectRegion(thflandmaskfile[ii], areaname=thflandmaskname[ii], box=thfbox,
                                                    **kwargs)
                if thf_landmask is not None:
                    break
    else:
        if isinstance(thffile, basestring):
            thf_landmask = ReadAreaSelectRegion(thffile, areaname=thflandmaskname, box=thfbox, **kwargs)
        else:
            for ii in range(len(thffile)):
                thf_landmask = ReadAreaSelectRegion(thffile[ii], areaname=thflandmaskname[ii], box=thfbox, **kwargs)
                if thf_landmask is not None:
                    break
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if thf_landmask is not None:
            dict_debug['axes2'] = '(thf) ' + str([ax.id for ax in thf_landmask.getAxisList()])
            dict_debug['shape2'] = '(thf) ' + str(thf_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if thf_landmask is not None:
        thf = ApplyLandmask(thf, thf_landmask, maskland=True, maskocean=False)
        if thf_areacell is None:
            thf_areacell = ArrayOnes(thf_landmask, id='areacell')
        thf_areacell = ApplyLandmaskToArea(thf_areacell, thf_landmask, maskland=True, maskocean=False)
    del thf_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, thf, keyerror = CheckTime(sst, thf, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        alpha, alphaPos, alphaNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        thf, unneeded = PreProcessTS(thf, '', areacell=thf_areacell, average='horizontal', compute_anom=True, **kwargs)
        del sst_areacell, thf_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(thf) ' + str([ax.id for ax in thf.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(thf) ' + str(thf.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(thf) ' + str(TimeBounds(thf))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        alpha, alphaPos, alphaNeg = LinearRegressionAndNonlinearity(thf, sst, return_stderr=True)

    # Create output
    alphaMetric = {
        'name': Name, 'value': alpha[0], 'value_error': alpha[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': alphaNeg[0] - alphaPos[0],
        'nonlinearity_error': alphaNeg[1] + alphaPos[1], 'keyerror': keyerror,
    }
    return alphaMetric


def EnsoAmpl(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, dataset='',
             debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoAmpl() function computes the standard deviation of 'sstbox' sstA (usually the standard deviation of nino3
    sstA)

    Author:	Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Yann Planton : yann.planton@locean-ipsl.upmc.fr

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, keyerror, dive_down_diag

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
    metric = 'EnsoAmpl'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror =\
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, sstbox, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, maskland=True,
                            maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        sstStd, sstStdErr, dive_down_diag = None, None, {'value': None, 'axis': None}
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        del sst_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the standard deviation
        sstStd = float(Std(sst))

        # Standard Error of the Standard Deviation (function of nyears)
        sstStdErr = sstStd / NUMPYsqrt(yearN)

        # Dive down diagnostic
        sstStd_monthly = StdMonthly(sst)
        dive_down_diag = {'value': ArrayToList(sstStd_monthly), 'axis': list(sstStd_monthly.getAxis(0)[:])}
        if netcdf is True:
            # additional diagnostic
            # Read file and select the right region
            sst, sst_areacell, unneeded =\
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafile, name_area=sstareaname, file_mask=sstlandmaskfile,
                                    name_mask=sstlandmaskname, maskland=True, maskocean=False, debug=debug, **kwargs)
            # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, unneeded = PreProcessTS(sst, '', areacell=sst_areacell, compute_anom=True, **kwargs)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape),
                              'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 10, **dict_debug)
            # std
            sst = Std(sst)
            # Regridding
            if 'regridding' not in kwargs.keys():
                kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                        'newgrid_name': 'generic_1x1deg'}
            else:
                if not isinstance(kwargs['regridding'], dict):
                    kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                            'newgrid_name': 'generic_1x1deg'}
            sst = Regrid(sst, None, region='equatorial_pacific_LatExt2', **kwargs['regridding'])
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "monthly standard deviation of " + sstbox + " sstA",
                     'diagnostic_value': sstStd, 'diagnostic_value_error': sstStdErr}
            dict2 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "standard deviation of " + sstbox + " sstA (whole time series)"}
            dict3 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sstStd_monthly, var1_attributes=dict1, var1_name='sstStd_monthly__' + dataset,
                       var2=sst, var2_attributes=dict2, var2_name='sstStd_map__' + dataset, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstStd), 'line2': 'metric value_error: ' + str(sstStdErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    amplMetric = {
        'name': Name, 'value': sstStd, 'value_error': sstStdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref, 'keyerror': keyerror,
        'dive_down_diag': dive_down_diag,
    }
    return amplMetric


def EnsoDiversity(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, box, event_definition,
                  dataset='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoDiversity() function computes a zonal composite of El Nino and La Nina events during the peak of the event.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA > (<) 'threshold' during 'season' are considered as El Nino (La Nina) events
        2.) diversity of the zonal location of the maximum (minimum) SSTA
            2.1) zonal SSTA at the peak of the event is computed for each selected event
            2.2) find the zonal position of the maximum (minimum) SSTA for each selected event
            2.3) compute the percentage of EP events (maximum/minimum SSTA eastward of the given threshold)
            2.4) compute the ratio EP events during La Nina divided by EP events during El Nino

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param treshold_ep_ev: float, optional
        see EnsoToolsLib.percentage_val_eastward
        longitude, in degree east, of the westward boundary of eastern Pacific event
        default value is -140°E (i.e., 140°W)
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return NinaDivMetric: dict
        name, value, value_error, units, method, nyears, events, time_frequency, time_period, ref, keyerror,
        dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'treshold_ep_ev',
                    'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO Diversity (percentage of eastern Pacific El Nino / La Nina)'
    lat = ReferenceRegions(box)['latitude']
    lon = ReferenceRegions(box)['longitude']
    Method = 'Nino (Nina) events = ' + region_ev + ' sstA > (<) ' + str(threshold) + ' during ' + season_ev +\
             ', zonal SSTA ' + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) +\
             ']), westward boundary of EP events ' + str(kwargs['treshold_ep_ev']) + 'E'
    Units = '%'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'EnsoDiversity'

    # ------------------------------------------------
    # 1. detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, areacell, keyerror = \
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, region_ev, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname,
                            maskland=True, maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        ratioEP, StdErr, dive_down_diag, nino_years, nina_years = None, None, {'value': None, 'axis': None}, None, None
    else:
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, unneeded = PreProcessTS(sst, '', areacell=areacell, average='horizontal', compute_anom=False, **kwargs)
        del areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA > (<) 'threshold' during 'season' are considered as El Nino (La Nina) events
        # Lists event years
        nino_years = DetectEvents(sst, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        nina_years = DetectEvents(sst, season_ev, -threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nino1': 'nbr(' + str(len(nino_years)) + '): ' + str(nino_years),
                          'nina1': 'nbr(' + str(len(nina_years)) + '): ' + str(nina_years)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(nina_years) > 0 and len(nino_years) > 0:
            # ------------------------------------------------
            # 2. diversity of the zonal location of the minimum SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst, areacell, unneeded = \
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, box, file_area=sstareafile,
                                    name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskfile,
                                    maskland=True, maskocean=False, debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, Method = PreProcessTS(sst, Method, areacell=areacell, average=False, compute_anom=False, **kwargs)
            del areacell
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst = SeasonalMean(sst, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder', 'regridTool',
                              'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst = Regrid(sst, None, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                                  'shape1': '(sst) ' + str(sst.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst = AverageMeridional(sst)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sample_nino = Composite_ev_by_ev(sst, nino_years, kwargs['frequency'])
            sample_nina = Composite_ev_by_ev(sst, nina_years, kwargs['frequency'])

            # 2.2 find the zonal position of the maximum/minimum SSTA for each selected event
            lon_sstmax = FindXYMinMaxInTs(sample_nino, return_val='maxi', smooth=True, axis=0, window=5,
                                          method='triangle')
            lon_sstmin = FindXYMinMaxInTs(sample_nina, return_val='mini', smooth=True, axis=0, window=5,
                                          method='triangle')
            if debug is True:
                dict_debug = {'line1': 'longitude of the maximum SSTA (nino): ' + str(lon_sstmax),
                              'line2': 'longitude of the minimum SSTA (nina): ' + str(lon_sstmin)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after FindXYMinMaxInTs', 15, **dict_debug)

            # 2.3 compute the percentage of EP events (maximum/minimum SSTA eastward of the given threshold)
            ep_event_nino, keyerror_nino = percentage_val_eastward(lon_sstmax, metric, box,
                                                                   threshold=kwargs['treshold_ep_ev'])
            ep_event_nina, keyerror_nina = percentage_val_eastward(lon_sstmin, metric, box,
                                                                   threshold=kwargs['treshold_ep_ev'])
            if debug is True:
                dict_debug = {'nino1': 'percentage of EP event + ' + str(ep_event_nino),
                              'nina1': 'percentage of EP event + ' + str(ep_event_nina)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

            # 2.4 compute the ratio EP events during La Nina divided by EP events during El Nino
            if ep_event_nino == 0:
                keyerror_nino = "no EP El Nino during the period"
            else:
                ratioEP = float(ep_event_nina / ep_event_nino)
        else:
            if len(nina_years) == 0:
                keyerror_nina = "no La Nina event during the period"
            else:
                keyerror_nina = None
            if len(nino_years) == 0:
                keyerror_nino = "no El Nino event during the period"
            else:
                keyerror_nino = None

        if keyerror_nino is not None or keyerror_nina is not None:
            ratioEP, StdErr, dive_down_diag = None, None, {'value': None, 'axis': None}
            keyerror = ''
            if keyerror_nino is not None:
                keyerror = keyerror_nino
            if len(keyerror) > 0 and keyerror_nina is not None:
                keyerror += " ; "
            if keyerror_nina is not None:
                keyerror += keyerror_nina
        else:
            # Standard Error of the Standard Deviation (function of nyears)
            StdErr = None

            # Dive down diagnostic
            dive_down_diag = {'value': ArrayToList(lon_sstmax), 'axis': list(lon_sstmax.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': 'longitude (E)', 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nino_years': str(nino_years), 'diagnostic_value_' + dataset: ratioEP,
                         'diagnostic_value_error_' + dataset: StdErr}
                dict2 = {'units': 'longitude (E)', 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nina_years': str(nina_years), 'diagnostic_value_' + dataset: ratioEP,
                         'diagnostic_value_error_' + dataset: StdErr}
                dict3 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=lon_sstmax, var1_attributes=dict1,
                           var1_name='Nino_lon_pos_maxSSTA__' + dataset, var2=lon_sstmin, var2_attributes=dict2,
                           var2_name='Nina_lon_pos_minSSTA__' + dataset, global_attributes=dict3)
                del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(ratioEP), 'line2': 'metric value_error: ' + str(StdErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    EnsoDivMetric = {
        'name': Name, 'value': ratioEP, 'value_error': StdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'events': nino_years+nina_years, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds,
        'ref': Ref, 'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return EnsoDivMetric


def EnsoMu(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, tauxfile, tauxname,
           tauxareafile, tauxareaname, tauxlandmaskfile, tauxlandmaskname, tauxbox, dataset='', debug=False,
           netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoMu() function computes the regression of 'tauxbox' tauxA (surface downward zonal stress anomalies) over
    'sstbox' sstA (usually the regression of nino4 tauxA over nino3 sstA)

    Author:	Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Yann Planton : yann.planton@locean-ipsl.upmc.fr

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param tauxfile: string
        path_to/filename of the file (NetCDF) of TAUX
    :param tauxname: string
        name of TAUX variable (taux, tauu) in 'tauxfile'
    :param tauxareafile: string
        path_to/filename of the file (NetCDF) of the areacell for TAUX
    :param tauxareaname: string
        name of areacell variable (areacella, areacello) in 'tauxareafile'
    :param tauxlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for TAUX
    :param tauxlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'tauxlandmaskfile'
    :param tauxbox: string
        name of box (nino4') for TAUX
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    Units = '10e-3 N/m2/C'
    Method = 'Regression of ' + tauxbox + ' tauxA over ' + sstbox + ' sstA'
    Method_NL = 'The nonlinearity is the regression computed when sstA<0 minus the regression computed when sstA>0'
    Ref = 'Using CDAT regression calculation'
    metric = 'EnsoMu'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
        dict_debug = {'file1': '(sst) ' + sstfile, 'file2': '(taux) ' + tauxfile, 'var1': '(sst) ' + sstname,
                      'var2': '(taux) ' + tauxname}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files', 10, **dict_debug)
    sst = ReadSelectRegionCheckUnits(sstfile, sstname, 'temperature', box=sstbox, **kwargs)
    taux = ReadSelectRegionCheckUnits(tauxfile, tauxname, 'wind stress', box=tauxbox, **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                      'axes2': '(thf) ' + str([ax.id for ax in taux.getAxisList()]),
                      'shape1': '(sst) ' + str(sst.shape), 'shape2': '(thf) ' + str(taux.shape),
                      'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(thf) ' + str(TimeBounds(taux))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)

    # Read areacell
    if sstareafile:
        sst_areacell = ReadAreaSelectRegion(sstareafile, areaname=sstareaname, box=sstbox, **kwargs)
    else:
        sst_areacell = ReadAreaSelectRegion(sstfile, areaname=sstareaname, box=sstbox, **kwargs)
    if tauxareafile:
        taux_areacell = ReadAreaSelectRegion(tauxareafile, areaname=tauxareaname, box=tauxbox, **kwargs)
    else:
        taux_areacell = ReadAreaSelectRegion(tauxfile, areaname=tauxareaname, box=tauxbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_areacell is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_areacell.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_areacell.shape)
        if taux_areacell is not None:
            dict_debug['axes2'] = '(taux) ' + str([ax.id for ax in taux_areacell.getAxisList()])
            dict_debug['shape2'] = '(taux) ' + str(taux_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)

    # Read landmask
    if sstlandmaskfile:
        sst_landmask = ReadLandmaskSelectRegion(sstlandmaskfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    else:
        sst_landmask = ReadLandmaskSelectRegion(sstfile, landmaskname=sstlandmaskname, box=sstbox,
                                                **kwargs)
    if tauxlandmaskfile:
        taux_landmask = ReadLandmaskSelectRegion(tauxlandmaskfile, landmaskname=tauxlandmaskname, box=tauxbox, **kwargs)
    else:
        taux_landmask = ReadLandmaskSelectRegion(tauxfile, landmaskname=tauxlandmaskname, box=tauxbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if sst_landmask is not None:
            dict_debug['axes1'] = '(sst) ' + str([ax.id for ax in sst_landmask.getAxisList()])
            dict_debug['shape1'] = '(sst) ' + str(sst_landmask.shape)
        if taux_landmask is not None:
            dict_debug['axes2'] = '(taux) ' + str([ax.id for ax in taux_landmask.getAxisList()])
            dict_debug['shape2'] = '(taux) ' + str(taux_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)

    # Apply landmask
    if sst_landmask is not None:
        sst = ApplyLandmask(sst, sst_landmask, maskland=True, maskocean=False)
        if sst_areacell is None:
            sst_areacell = ArrayOnes(sst_landmask, id='areacell')
        sst_areacell = ApplyLandmaskToArea(sst_areacell, sst_landmask, maskland=True, maskocean=False)
        del sst_landmask
    if taux_landmask is not None:
        taux = ApplyLandmask(taux, taux_landmask, maskland=True, maskocean=False)
        if taux_areacell is None:
            taux_areacell = ArrayOnes(taux_landmask, id='areacell')
        taux_areacell = ApplyLandmaskToArea(taux_areacell, taux_landmask, maskland=True, maskocean=False)
    del taux_landmask

    # Checks if the same time period is used for both variables and if the minimum number of time steps is respected
    sst, taux, keyerror = CheckTime(sst, taux, metric_name=metric, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        mu, muPos, muNeg = None, None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smooths TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        taux, unneeded = PreProcessTS(taux, '', areacell=taux_areacell, average='horizontal', compute_anom=True,
                                      **kwargs)
        del sst_areacell, taux_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'axes2': '(taux) ' + str([ax.id for ax in taux.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'shape2': '(taux) ' + str(taux.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst)), 'time2': '(thf) ' + str(TimeBounds(taux))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the linear regression for all points, for SSTA >=0 and for SSTA<=0
        mu, muPos, muNeg = LinearRegressionAndNonlinearity(taux, sst, return_stderr=True)

        # Change units
        mu = [mu[0] * 1e3, mu[1] * 1e3]
        muPos = [muPos[0] * 1e3, muPos[1] * 1e3]
        muNeg = [muNeg[0] * 1e3, muNeg[1] * 1e3]

    # Create output
    muMetric = {
        'name': Name, 'value': mu[0], 'value_error': mu[1], 'units': Units, 'method': Method,
        'method_nonlinearity': Method_NL, 'nyears': yearN, 'time_frequency': kwargs['frequency'],
        'time_period': actualtimebounds, 'ref': Ref, 'nonlinearity': muNeg[0] - muPos[0],
        'nonlinearity_error': muNeg[1] + muPos[1], 'keyerror': keyerror,
    }
    return muMetric


def EnsoPrMap(sstfilemodel, sstnamemodel, sstareafilemodel, sstareanamemodel, sstlandmaskfilemodel,
              sstlandmasknamemodel, prfilemodel, prnamemodel, prareafilemodel, prareanamemodel,
              prlandmaskfilemodel, prlandmasknamemodel, sstfileobs, sstnameobs, sstareafileobs, sstareanameobs,
              sstlandmaskfileobs, sstlandmasknameobs, prfileobs, prnameobs, prareafileobs, prareanameobs,
              prlandmaskfileobs, prlandmasknameobs, sstbox, prbox, event_definition, centered_rmse=0,
              biased_rmse=1, debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoPrMap() function computes precipitation anomalies pattern associated with ENSO on the globe.
    First metric: rmse(observations vs model).
    Second metric: correlation(observations vs model).
    Third metric: std(model)/std(observations)
    These metrics can be used to compute a Taylor diagram.

    Inputs:
    ------
    :param sstfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemodel: string
        name of SST variable (tos, ts) in 'sstfilemodel'
    :param sstareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST areacell
    :param sstareanamemodel: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafilemodel'
    :param sstlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST landmask
    :param sstlandmasknamemodel: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfilemodel'
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PRvariable (pr) in 'prfilemodel'
    :param prareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR areacell
    :param prareanamemodel: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafilemodel'
    :param prlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR landmask
    :param prlandmasknamemodel: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfilemodel'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST areacell
    :param sstareanameobs: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafileobs'
    :param sstlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST landmask
    :param sstlandmasknameobs: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfileobs'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, precip) in 'prfileobs'
    :param prareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR areacell
    :param prareanameobs: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafileobs'
    :param prlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR landmask
    :param prlandmasknameobs: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfileobs'
    :param sstbox: string
        name of box (e.g. 'global') for SST
    :param prbox: string
        name of box (e.g. 'global') for PR
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return EnsoPrMapMetric: dict
        name, value (rms [obs;model]), value_error, units, method, value2 (corr [obs;model]),
        value_error2, units2, value3 (std_model / std_obs), value_error3, units3, nyears_model, nyears_observations,
        time_frequency, time_period_model, time_period_observations, ref, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO PRA pattern '
    Method = region_ev + 'SSTA during ' + season_ev + ' regressed against precipitation anomalies in ' + prbox
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'mm/day / C'
    Ref = 'Using CDAT regridding, correlation (centered and biased), std (centered and biased) and ' +\
          'rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # compute ENSO time series
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', 'EnsoPrMap', 10)
        dict_debug = {'file1': '(model) ' + sstfilemodel, 'file2': '(obs) ' + sstfileobs,
                      'var1': '(model) ' + sstnamemodel, 'var2': '(obs) ' + sstnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files ENSO', 10, **dict_debug)
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_mod'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if sstareafilemodel:
        model_areacell = ReadAreaSelectRegion(sstareafilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(sstfilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    if sstareafileobs:
        obs_areacell = ReadAreaSelectRegion(sstareafileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(sstfileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read landmask
    if sstlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(sstlandmaskfilemodel, landmaskname=sstlandmasknamemodel,
                                                  box=region_ev, **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(sstfilemodel, landmaskname=sstlandmasknamemodel, box=region_ev,
                                                  **kwargs)
    if sstlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(sstlandmaskfileobs, landmaskname=sstlandmasknameobs,
                                                  box=region_ev, **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(sstfileobs, landmaskname=sstlandmasknameobs,
                                                box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        sst_model = ApplyLandmask(sst_model, model_landmask, maskland=True, maskocean=False)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=True, maskocean=False)
        del model_landmask
    if obs_landmask is not None:
        sst_obs = ApplyLandmask(sst_obs, obs_landmask, maskland=True, maskocean=False)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=True, maskocean=False)
        del obs_landmask

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrMap: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrMap: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=False,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                     **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Seasonal mean and anomalies
    enso_model = SeasonalMean(sst_model, season_ev, compute_anom=True)
    enso_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
    del sst_model, sst_obs
    if season_ev == 'DJF':
        time_ax = enso_model.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
        enso_model.setAxis(0, time_ax)
        del time_ax
        time_ax = enso_obs.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
        enso_obs.setAxis(0, time_ax)
        del time_ax
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in enso_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in enso_obs.getAxisList()]),
                      'shape1': '(model) ' + str(enso_model.shape), 'shape2': '(obs) ' + str(enso_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(enso_model)), 'time2': '(obs) ' + str(TimeBounds(enso_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

    # ------------------------------------------------
    # regression map
    # ------------------------------------------------
    if debug is True:
        dict_debug = {
            'file1': '(model) ' + prfilemodel, 'file2': '(obs) ' + prfileobs, 'var1': '(model) ' + prnamemodel,
            'var2': '(obs) ' + prnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files regression map', 10, **dict_debug)
    # Read file and select the right region
    pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=prbox,
                                          time_bounds=kwargs['time_bounds_mod'], **kwargs)
    pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=prbox,
                                        time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                      'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if prareafilemodel:
        model_areacell = ReadAreaSelectRegion(prareafilemodel, areaname=prareanamemodel, box=prbox, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(prfilemodel, areaname=prareanamemodel, box=prbox, **kwargs)
    if prareafileobs:
        obs_areacell = ReadAreaSelectRegion(prareafileobs, areaname=prareanameobs, box=prbox, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(prfileobs, areaname=prareanameobs, box=prbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read if the given region is defined as a land region, an oceanic region, or both
    dict_reg = ReferenceRegions(prbox)
    if 'maskland' in dict_reg.keys():
        maskland = dict_reg['maskland']
    else:
        maskland = False
    if 'maskocean' in dict_reg.keys():
        maskocean = dict_reg['maskocean']
    else:
        maskocean = False
    # Read landmask
    if prlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(prlandmaskfilemodel, landmaskname=prlandmasknamemodel, box=prbox,
                                                  **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(prfilemodel, landmaskname=prlandmasknamemodel, box=prbox, **kwargs)
    if prlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(prlandmaskfileobs, landmaskname=prlandmasknameobs, box=prbox, **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(prfileobs, landmaskname=prlandmasknameobs, box=prbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        pr_model = ApplyLandmask(pr_model, model_landmask, maskland=maskland, maskocean=maskocean)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=maskland, maskocean=maskocean)
        del model_landmask
    if obs_landmask is not None:
        pr_obs = ApplyLandmask(pr_obs, obs_landmask, maskland=maskland, maskocean=maskocean)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=maskland, maskocean=maskocean)
        del obs_landmask

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, compute_anom=False, **kwargs)
    pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, compute_anom=False, **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                      'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Seasonal mean
    pr_model = SeasonalMean(pr_model, season_ev, compute_anom=True)
    pr_obs = SeasonalMean(pr_obs, season_ev, compute_anom=True)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                      'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        pr_model, pr_obs, Method = TwoVarRegrid(pr_model, pr_obs, Method, region=prbox, **kwargs['regridding'])
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

    # regression
    pr_model_slope, pr_model_stderr = LinearRegressionTsAgainstMap(pr_model, enso_model, return_stderr=True)
    pr_obs_slope, pr_obs_stderr = LinearRegressionTsAgainstMap(pr_obs, enso_obs, return_stderr=True)

    # Metric 1
    prRmse = float(RmsAxis(pr_model_slope, pr_obs_slope, axis='xy', centered=centered_rmse, biased=biased_rmse))
    # Metric 2
    prCorr = float(Correlation(pr_model_slope, pr_obs_slope, axis='xy', centered=1, biased=1))
    # Metric 3
    std_model = Std(pr_model_slope, weights=None, axis='xy', centered=1, biased=1)
    std_obs = Std(pr_obs_slope, weights=None, axis='xy', centered=1, biased=1)
    prStd = float(std_model) / float(std_obs)

    # Dive down diagnostic
    dive_down_diag = {'model': pr_model_slope, 'observations': pr_obs_slope,
                      'axisLat': list(pr_model_slope.getAxis(0)[:]), 'axisLon': list(pr_obs_slope.getAxis(1)[:])}

    # Create output
    EnsoPrMapMetric = {
        'name': Name, 'value': prRmse, 'value_error': None, 'units': Units, 'method': Method,
        'value2': prCorr, 'value_error2': None, 'units2': '', 'value3': prStd, 'value_error3': None, 'units3': Units,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
        'dive_down_diag': dive_down_diag,
    }
    return EnsoPrMapMetric


def EnsoPrJjaTel(sstfilemodel, sstnamemodel, sstareafilemodel, sstareanamemodel, sstlandmaskfilemodel,
                 sstlandmasknamemodel, prfilemodel, prnamemodel, prareafilemodel, prareanamemodel, prlandmaskfilemodel,
                 prlandmasknamemodel, sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs,
                 sstlandmasknameobs, prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs,
                 prlandmasknameobs, sstbox, prbox, event_definition, centered_rmse=0, biased_rmse=1, debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoPrJjaTel() function computes precipitations anomalies associated with El Nino and La Nina events in many AR5
        reference regions, then precipitations in JJA are composited for each selected event and the difference
        (El Nino PR - La Nina PR) is computed in each region.
    The first rmse(observations vs model) is the metric.
    The second metric is the number of regions where observations and models agree on the sign of the teleconnection

    Inputs:
    ------
    :param sstfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemodel: string
        name of SST variable (tos, ts) in 'sstfilemodel'
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PRvariable (pr) in 'prfilemodel'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, precip) in 'prfileobs'
    :param box: string
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param sstareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST areacell
    :param sstareanamemodel: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafilemodel'
    :param sstlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST landmask
    :param sstlandmasknamemodel: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfilemodel'
    :param prareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR areacell
    :param prareanamemodel: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafilemodel'
    :param prlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR landmask
    :param prlandmasknamemodel: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfilemodel'
    :param sstareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST areacell
    :param sstareanameobs: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafileobs'
    :param sstlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST landmask
    :param sstlandmasknameobs: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfileobs'
    :param prareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR areacell
    :param prareanameobs: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafileobs'
    :param prlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR landmask
    :param prlandmasknameobs: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfileobs'
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return EnsoPrTelMetric: dict
        name, value (rms [NinoPr-NinaPr]), value_error, units, method, value2 (sign agreement [NinoPr-NinaPr]),
        value_error2, units2, nyears_model, nyears_observations, nina_model, nino_model, nina_observations,
        nino_observations, time_frequency, time_period_model, time_period_observations, ref, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino composite minus Nina composite during JJA preceding the events in each region'
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ', Nina events = ' + region_ev + ' sstA < -'\
             + str(threshold) + ' during ' + season_ev + '; Precipitations associated with El Nino/La Nina' + \
             'events during the preceding JJA are composited and the difference (El Nino PR - La Nina PR) is' + \
             'computed in each region'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'mm/day'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', 'EnsoPrJjaTel', 10)
        dict_debug = {'file1': '(model) ' + sstfilemodel, 'file2': '(obs) ' + sstfileobs,
                      'var1': '(model) ' + sstnamemodel, 'var2': '(obs) ' + sstnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files ENSO', 10, **dict_debug)
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_mod'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if sstareafilemodel:
        model_areacell = ReadAreaSelectRegion(sstareafilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(sstfilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    if sstareafileobs:
        obs_areacell = ReadAreaSelectRegion(sstareafileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(sstfileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read landmask
    if sstlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(sstlandmaskfilemodel, landmaskname=sstlandmasknamemodel,
                                                  box=region_ev, **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(sstfilemodel, landmaskname=sstlandmasknamemodel, box=region_ev,
                                                  **kwargs)
    if sstlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(sstlandmaskfileobs, landmaskname=sstlandmasknameobs,
                                                  box=region_ev, **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(sstfileobs, landmaskname=sstlandmasknameobs,
                                                box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        sst_model = ApplyLandmask(sst_model, model_landmask, maskland=True, maskocean=False)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=True, maskocean=False)
        del model_landmask
    if obs_landmask is not None:
        sst_obs = ApplyLandmask(sst_obs, obs_landmask, maskland=True, maskocean=False)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=True, maskocean=False)
        del obs_landmask

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrJjaTel: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrJjaTel: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=False,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                     **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Lists event years
    nina_years_model = DetectEvents(sst_model, season_ev, -threshold, normalization=kwargs['normalization'], nino=False)
    nino_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    nina_years_obs = DetectEvents(sst_obs, season_ev, -threshold, normalization=kwargs['normalization'], nino=False)
    nino_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    del sst_model, sst_obs
    if debug is True:
        dict_debug = {'nina1': '(model) ' + str(nina_years_model), 'nina2': '(obs) ' + str(nina_years_obs),
                      'nino1': '(model) ' + str(nino_years_model), 'nino2': '(obs) ' + str(nino_years_obs)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

    # ------------------------------------------------
    # compute composite
    # ------------------------------------------------
    if debug is True:
        dict_debug = {
            'file1': '(model) ' + prfilemodel, 'file2': '(obs) ' + prfileobs, 'var1': '(model) ' + prnamemodel,
            'var2': '(obs) ' + prnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files Composite', 10, **dict_debug)
    # smoothing is not applied
    if 'smoothing' in kwargs.keys():
        smooth = deepcopy(kwargs['smoothing'])
        kwargs['smoothing'] = False
    if not isinstance(prbox, list):
        prbox = [prbox]
    prbox = sorted(prbox, key=str.lower)
    list_composite_model, list_composite_obs = list(), list()
    for reg in prbox:
        if debug is True:
            EnsoErrorsWarnings.DebugMode('\033[92m', 'region = '+str(reg), 10)
        # Read file and select the right region
        pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=reg,
                                              time_bounds=kwargs['time_bounds_mod'], **kwargs)
        pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=reg,
                                            time_bounds=kwargs['time_bounds_obs'], **kwargs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
        # Read areacell
        if prareafilemodel:
            model_areacell = ReadAreaSelectRegion(prareafilemodel, areaname=prareanamemodel, box=reg, **kwargs)
        else:
            model_areacell = ReadAreaSelectRegion(prfilemodel, areaname=prareanamemodel, box=reg, **kwargs)
        if prareafileobs:
            obs_areacell = ReadAreaSelectRegion(prareafileobs, areaname=prareanameobs, box=reg, **kwargs)
        else:
            obs_areacell = ReadAreaSelectRegion(prfileobs, areaname=prareanameobs, box=reg, **kwargs)
        if debug is True:
            dict_debug = {}
            if model_areacell is not None:
                dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
                dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
            if obs_areacell is not None:
                dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
                dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
        # Read if the given region is defined as a land region, an oceanic region, or both
        dict_reg = ReferenceRegions(reg)
        if 'maskland' in dict_reg.keys():
            maskland = dict_reg['maskland']
        else:
            maskland = False
        if 'maskocean' in dict_reg.keys():
            maskocean = dict_reg['maskocean']
        else:
            maskocean = False
        # Read landmask
        if prlandmaskfilemodel:
            model_landmask = ReadLandmaskSelectRegion(prlandmaskfilemodel, landmaskname=prlandmasknamemodel, box=reg,
                                                      **kwargs)
        else:
            model_landmask = ReadLandmaskSelectRegion(prfilemodel, landmaskname=prlandmasknamemodel, box=reg, **kwargs)
        if prlandmaskfileobs:
            obs_landmask = ReadLandmaskSelectRegion(prlandmaskfileobs, landmaskname=prlandmasknameobs, box=reg,
                                                    **kwargs)
        else:
            obs_landmask = ReadLandmaskSelectRegion(prfileobs, landmaskname=prlandmasknameobs, box=reg, **kwargs)
        if debug is True:
            dict_debug = {}
            if model_landmask is not None:
                dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
                dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
            if obs_landmask is not None:
                dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
                dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
        # Apply landmask
        if model_landmask is not None:
            pr_model = ApplyLandmask(pr_model, model_landmask, maskland=maskland, maskocean=maskocean)
            if model_areacell is None:
                model_areacell = ArrayOnes(model_landmask, id='areacell')
            model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=maskland, maskocean=maskocean)
            del model_landmask
        if obs_landmask is not None:
            pr_obs = ApplyLandmask(pr_obs, obs_landmask, maskland=maskland, maskocean=maskocean)
            if obs_areacell is None:
                obs_areacell = ArrayOnes(obs_landmask, id='areacell')
            obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=maskland, maskocean=maskocean)
            del obs_landmask

        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, average='horizontal',
                                        compute_anom=False, **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                        **kwargs)
        del model_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Seasonal mean
        pr_model = SeasonalMean(pr_model, 'JJA', compute_anom=False)
        pr_obs = SeasonalMean(pr_obs, 'JJA', compute_anom=False)

        # composites
        composite_nina_model = Composite(pr_model, nina_years_model, kwargs['frequency'])
        composite_nino_model = Composite(pr_model, nino_years_model, kwargs['frequency'])
        composite_nina_obs = Composite(pr_obs, nina_years_obs, kwargs['frequency'])
        composite_nino_obs = Composite(pr_obs, nino_years_obs, kwargs['frequency'])

        # list composites
        list_composite_model.append(float(composite_nino_model-composite_nina_model))
        list_composite_obs.append(float(composite_nino_obs-composite_nina_obs))
    if 'smoothing' in kwargs.keys():
        kwargs['smoothing'] = smooth
        del smooth

    # Computes the root mean square difference
    compositeRmse = RmsAxis(list_composite_model, list_composite_obs, centered=centered_rmse, biased=biased_rmse)

    # Computes the percentage of regions where observations and model agree on the sign of the teleconnection
    signAgreement = sum([1. for vmod,vobs in zip(list_composite_model,list_composite_obs)
                        if NUMPYsign(vmod)==NUMPYsign(vobs)])/len(list_composite_model)

    # Dive down diagnostic
    dive_down_diag = {'model': list_composite_model, 'observations': list_composite_obs, 'axis': prbox}

    # Create output
    EnsoPrTelMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'value2': signAgreement, 'value_error2': None, 'units2': '%', 'nyears_model': yearN_model,
        'nyears_observations': yearN_obs, 'nina_model': nina_years_model, 'nino_model': nino_years_model,
        'nina_observations': nina_years_obs, 'nino_observations': nino_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
        'dive_down_diag': dive_down_diag,
    }
    return EnsoPrTelMetric


def EnsoPrNdjTel(sstfilemodel, sstnamemodel, sstareafilemodel, sstareanamemodel, sstlandmaskfilemodel,
                 sstlandmasknamemodel, prfilemodel, prnamemodel, prareafilemodel, prareanamemodel, prlandmaskfilemodel,
                 prlandmasknamemodel, sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs,
                 sstlandmasknameobs, prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs,
                 prlandmasknameobs, sstbox, prbox, event_definition, centered_rmse=0, biased_rmse=1, debug=False,
                 netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoPrNdjTel() function computes precipitations anomalies associated with El Nino and La Nina events in many AR5
        reference regions, then precipitations in NDJ are composited for each selected event and the difference
        (El Nino PR - La Nina PR) is computed in each region.
    The first rmse(observations vs model) is the metric.
    The second metric is the number of regions where observations and models agree on the sign of the teleconnection

    Inputs:
    ------
    :param sstfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemodel: string
        name of SST variable (tos, ts) in 'sstfilemodel'
    :param prfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemodel: string
        name of PRvariable (pr) in 'prfilemodel'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, precip) in 'prfileobs'
    :param box: string
        name of box (e.g. 'nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param sstareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST areacell
    :param sstareanamemodel: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafilemodel'
    :param sstlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST landmask
    :param sstlandmasknamemodel: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfilemodel'
    :param prareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR areacell
    :param prareanamemodel: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafilemodel'
    :param prlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled PR landmask
    :param prlandmasknamemodel: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfilemodel'
    :param sstareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST areacell
    :param sstareanameobs: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafileobs'
    :param sstlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST landmask
    :param sstlandmasknameobs: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfileobs'
    :param prareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR areacell
    :param prareanameobs: string, optional
        name of areacell for the PR variable (areacella, areacello,...) in 'prareafileobs'
    :param prlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed PR landmask
    :param prlandmasknameobs: string, optional
        name of landmask for the PR variable (sftlf,...) in 'prlandmaskfileobs'
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return EnsoPrTelMetric: dict
        name, value (rms [NinoPr-NinaPr]), value_error, units, method, value2 (sign agreement [NinoPr-NinaPr]),
        value_error2, units2, nyears_model, nyears_observations, nina_model, nino_model, nina_observations,
        nino_observations, time_frequency, time_period_model, time_period_observations, ref, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino composite minus Nina composite during JJA preceding the events in each region'
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ', Nina events = ' + region_ev + ' sstA < -'\
             + str(threshold) + ' during ' + season_ev + '; Precipitations associated with El Nino/La Nina' + \
             'events during the preceding JJA are composited and the difference (El Nino PR - La Nina PR) is' + \
             'computed in each region'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'mm/day'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', 'EnsoPrNdjTel', 10)
        dict_debug = {'file1': '(model) ' + sstfilemodel, 'file2': '(obs) ' + sstfileobs,
                      'var1': '(model) ' + sstnamemodel, 'var2': '(obs) ' + sstnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files ENSO', 10, **dict_debug)
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_mod'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if sstareafilemodel:
        model_areacell = ReadAreaSelectRegion(sstareafilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(sstfilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    if sstareafileobs:
        obs_areacell = ReadAreaSelectRegion(sstareafileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(sstfileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read landmask
    if sstlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(sstlandmaskfilemodel, landmaskname=sstlandmasknamemodel,
                                                  box=region_ev, **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(sstfilemodel, landmaskname=sstlandmasknamemodel, box=region_ev,
                                                  **kwargs)
    if sstlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(sstlandmaskfileobs, landmaskname=sstlandmasknameobs,
                                                  box=region_ev, **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(sstfileobs, landmaskname=sstlandmasknameobs,
                                                box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        sst_model = ApplyLandmask(sst_model, model_landmask, maskland=True, maskocean=False)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=True, maskocean=False)
        del model_landmask
    if obs_landmask is not None:
        sst_obs = ApplyLandmask(sst_obs, obs_landmask, maskland=True, maskocean=False)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=True, maskocean=False)
        del obs_landmask

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrNdjTel: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoPrNdjTel: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=False,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                     **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Lists event years
    nina_years_model = DetectEvents(sst_model, season_ev, -threshold, normalization=kwargs['normalization'], nino=False)
    nino_years_model = DetectEvents(sst_model, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    nina_years_obs = DetectEvents(sst_obs, season_ev, -threshold, normalization=kwargs['normalization'], nino=False)
    nino_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
    del sst_model, sst_obs
    if debug is True:
        dict_debug = {'nina1': '(model) ' + str(nina_years_model), 'nina2': '(obs) ' + str(nina_years_obs),
                      'nino1': '(model) ' + str(nino_years_model), 'nino2': '(obs) ' + str(nino_years_obs)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

    # ------------------------------------------------
    # compute composite
    # ------------------------------------------------
    if debug is True:
        dict_debug = {
            'file1': '(model) ' + prfilemodel, 'file2': '(obs) ' + prfileobs, 'var1': '(model) ' + prnamemodel,
            'var2': '(obs) ' + prnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files Composite', 10, **dict_debug)
    # smoothing is not applied
    if 'smoothing' in kwargs.keys():
        smooth = deepcopy(kwargs['smoothing'])
        kwargs['smoothing'] = False
    if not isinstance(prbox, list):
        prbox = [prbox]
    prbox = sorted(prbox, key=str.lower)
    list_composite_model, list_composite_obs = list(), list()
    for reg in prbox:
        if debug is True:
            EnsoErrorsWarnings.DebugMode('\033[92m', 'region = '+str(reg), 10)
        # Read file and select the right region
        pr_model = ReadSelectRegionCheckUnits(prfilemodel, prnamemodel, 'precipitations', box=reg,
                                              time_bounds=kwargs['time_bounds_mod'], **kwargs)
        pr_obs = ReadSelectRegionCheckUnits(prfileobs, prnameobs, 'precipitations', box=reg,
                                            time_bounds=kwargs['time_bounds_obs'], **kwargs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
        # Read areacell
        if prareafilemodel:
            model_areacell = ReadAreaSelectRegion(prareafilemodel, areaname=prareanamemodel, box=reg, **kwargs)
        else:
            model_areacell = ReadAreaSelectRegion(prfilemodel, areaname=prareanamemodel, box=reg, **kwargs)
        if prareafileobs:
            obs_areacell = ReadAreaSelectRegion(prareafileobs, areaname=prareanameobs, box=reg, **kwargs)
        else:
            obs_areacell = ReadAreaSelectRegion(prfileobs, areaname=prareanameobs, box=reg, **kwargs)
        if debug is True:
            dict_debug = {}
            if model_areacell is not None:
                dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
                dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
            if obs_areacell is not None:
                dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
                dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
        # Read if the given region is defined as a land region, an oceanic region, or both
        dict_reg = ReferenceRegions(reg)
        if 'maskland' in dict_reg.keys():
            maskland = dict_reg['maskland']
        else:
            maskland = False
        if 'maskocean' in dict_reg.keys():
            maskocean = dict_reg['maskocean']
        else:
            maskocean = False
        # Read landmask
        if prlandmaskfilemodel:
            model_landmask = ReadLandmaskSelectRegion(prlandmaskfilemodel, landmaskname=prlandmasknamemodel, box=reg,
                                                      **kwargs)
        else:
            model_landmask = ReadLandmaskSelectRegion(prfilemodel, landmaskname=prlandmasknamemodel, box=reg, **kwargs)
        if prlandmaskfileobs:
            obs_landmask = ReadLandmaskSelectRegion(prlandmaskfileobs, landmaskname=prlandmasknameobs, box=reg,
                                                    **kwargs)
        else:
            obs_landmask = ReadLandmaskSelectRegion(prfileobs, landmaskname=prlandmasknameobs, box=reg, **kwargs)
        if debug is True:
            dict_debug = {}
            if model_landmask is not None:
                dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
                dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
            if obs_landmask is not None:
                dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
                dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
        # Apply landmask
        if model_landmask is not None:
            pr_model = ApplyLandmask(pr_model, model_landmask, maskland=maskland, maskocean=maskocean)
            if model_areacell is None:
                model_areacell = ArrayOnes(model_landmask, id='areacell')
            model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=maskland, maskocean=maskocean)
            del model_landmask
        if obs_landmask is not None:
            pr_obs = ApplyLandmask(pr_obs, obs_landmask, maskland=maskland, maskocean=maskocean)
            if obs_areacell is None:
                obs_areacell = ArrayOnes(obs_landmask, id='areacell')
            obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=maskland, maskocean=maskocean)
            del obs_landmask

        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        pr_model, Method = PreProcessTS(pr_model, Method, areacell=model_areacell, average='horizontal',
                                        compute_anom=False, **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                        **kwargs)
        del model_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_model.shape), 'shape2': '(obs) ' + str(pr_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(pr_model)), 'time2': '(obs) ' + str(TimeBounds(pr_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Seasonal mean
        pr_model = SeasonalMean(pr_model, 'NDJ', compute_anom=False)
        pr_obs = SeasonalMean(pr_obs, 'NDJ', compute_anom=False)

        # composites
        composite_nina_model = Composite(pr_model, nina_years_model, kwargs['frequency'])
        composite_nino_model = Composite(pr_model, nino_years_model, kwargs['frequency'])
        composite_nina_obs = Composite(pr_obs, nina_years_obs, kwargs['frequency'])
        composite_nino_obs = Composite(pr_obs, nino_years_obs, kwargs['frequency'])

        # list composites
        list_composite_model.append(float(composite_nino_model-composite_nina_model))
        list_composite_obs.append(float(composite_nino_obs-composite_nina_obs))
    if 'smoothing' in kwargs.keys():
        kwargs['smoothing'] = smooth
        del smooth

    # Computes the root mean square difference
    compositeRmse = RmsAxis(list_composite_model, list_composite_obs, centered=centered_rmse, biased=biased_rmse)

    # Computes the percentage of regions where observations and model agree on the sign of the teleconnection
    signAgreement = sum([1. for vmod,vobs in zip(list_composite_model,list_composite_obs)
                        if NUMPYsign(vmod)==NUMPYsign(vobs)])/len(list_composite_model)

    # Dive down diagnostic
    dive_down_diag = {'model': list_composite_model, 'observations': list_composite_obs, 'axis': prbox}

    # Create output
    EnsoPrTelMetric = {
        'name': Name, 'value': compositeRmse, 'value_error': None, 'units': Units, 'method': Method,
        'value2': signAgreement, 'value_error2': None, 'units2': '%', 'nyears_model': yearN_model,
        'nyears_observations': yearN_obs, 'nina_model': nina_years_model, 'nino_model': nino_years_model,
        'nina_observations': nina_years_obs, 'nino_observations': nino_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
        'dive_down_diag': dive_down_diag,
    }
    return EnsoPrTelMetric


def EnsoSeasonality(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, dataset='',
                    debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoSeasonality() function computes ratio between the November-December-January (NDJ) and March-April-May (MAM)
    average standard deviation of 'sstbox' sstA (usually nino3 sstA)

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, keyerror, dive_down_diag

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
    Name = 'ENSO seasonality'
    Units = 'C'
    Method = 'Ratio between NDJ and MAM standard deviation ' + sstbox + ' sstA'
    Ref = 'Using CDAT std dev calculation'
    metric = 'EnsoSeasonality'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror = \
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, sstbox, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, maskland=True,
                            maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        ratioStd, ratio_std_err, dive_down_diag = None, None, {'value': None, 'axis': None}
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_ts, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=False,
                                      **kwargs)
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst_ts.getAxisList()]),
                          'shape1': '(sst) ' + str(sst_ts.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst_ts))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Seasonal mean
        sst_NDJ = SeasonalMean(sst_ts, 'NDJ', compute_anom=True)
        sst_MAM = SeasonalMean(sst_ts, 'MAM', compute_anom=True)
        if debug is True:
            dict_debug = {'axes1': '(sst_NDJ) ' + str([ax.id for ax in sst_NDJ.getAxisList()]),
                          'shape1': '(sst_NDJ) ' + str(sst_NDJ.shape),
                          'axes2': '(sst_NDJ) ' + str([ax.id for ax in sst_MAM.getAxisList()]),
                          'shape2': '(sst_NDJ) ' + str(sst_MAM.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

        # Compute std dev and ratio
        sst_NDJ_std = Std(sst_NDJ)
        sst_MAM_std = Std(sst_MAM)
        ratioStd = float(sst_NDJ_std / sst_MAM_std)

        # Standard Error of the Standard Deviation (function of nyears)
        sst_NDJ_std_err = sst_NDJ_std / NUMPYsqrt(yearN - 1)
        sst_MAM_std_err = sst_MAM_std / NUMPYsqrt(yearN)
        # The error (dy) on ratio ('y = x/z'): dy = (z*dx + x*dz) / z2
        ratio_std_err = float(
            (sst_MAM_std * sst_NDJ_std_err + sst_NDJ_std * sst_MAM_std_err) / NUMPYsquare(sst_MAM_std))

        # Dive down diagnostic
        sst, unneeded = PreProcessTS(sst, '', areacell=sst_areacell, average='horizontal', compute_anom=True, **kwargs)
        del sst_areacell
        sstStd_monthly = StdMonthly(sst)
        dive_down_diag = {'value': ArrayToList(sstStd_monthly), 'axis': list(sstStd_monthly.getAxis(0)[:])}
        if netcdf is True:
            # additional diagnostic
            # Read file and select the right region
            sst, sst_areacell, unneeded = \
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafile, name_area=sstareaname, file_mask=sstlandmaskfile,
                                    name_mask=sstlandmaskname, maskland=True, maskocean=False, debug=debug, **kwargs)
            # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, unneeded = PreProcessTS(sst, '', areacell=sst_areacell, compute_anom=False, **kwargs)
            del sst_areacell
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape),
                              'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS: netcdf', 10, **dict_debug)
            # Seasonal mean
            sst_NDJ = SeasonalMean(sst, 'NDJ', compute_anom=True)
            sst_MAM = SeasonalMean(sst, 'MAM', compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(sst_NDJ) ' + str([ax.id for ax in sst_NDJ.getAxisList()]),
                              'shape1': '(sst_NDJ) ' + str(sst_NDJ.shape),
                              'axes2': '(sst_NDJ) ' + str([ax.id for ax in sst_MAM.getAxisList()]),
                              'shape2': '(sst_NDJ) ' + str(sst_MAM.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean: netcdf', 15, **dict_debug)
            # Compute std dev and ratio
            sst_NDJ_std = Std(sst_NDJ)
            sst_MAM_std = Std(sst_MAM)
            # Regridding
            if 'regridding' not in kwargs.keys():
                kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                        'newgrid_name': 'generic_1x1deg'}
            else:
                if not isinstance(kwargs['regridding'], dict):
                    kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                            'newgrid_name': 'generic_1x1deg'}
            sst_NDJ_std = Regrid(sst_NDJ_std, None, region='equatorial_pacific_LatExt2', **kwargs['regridding'])
            sst_MAM_std = Regrid(sst_MAM_std, None, region='equatorial_pacific_LatExt2', **kwargs['regridding'])
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "monthly standard deviation of " + sstbox + " sstA",
                     'diagnostic_value': ratioStd, 'diagnostic_value_error': ratio_std_err}
            dict2 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "standard deviation of " + sstbox + " sstA (during NDJ)"}
            dict3 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "standard deviation of " + sstbox + " sstA (during MAM)"}
            dict4 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sstStd_monthly, var1_attributes=dict1, var1_name='sstStd_monthly__' + dataset,
                       var2=sst_NDJ_std, var2_attributes=dict2, var2_name='sstStd_NDJ_map__' + dataset,
                       var3=sst_MAM_std, var3_attributes=dict3, var3_name='sstStd_MAM_map__' + dataset,
                       global_attributes=dict4)
            del dict1, dict2, dict3, dict4
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(ratioStd), 'line2': 'metric value_error: ' + str(ratio_std_err)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    seaMetric = {
        'name': Name, 'value': ratioStd, 'value_error': ratio_std_err, 'units': Units, 'method': Method,
        'nyears': yearN, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return seaMetric


def EnsoSstSkew(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, dataset='',
                debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The EnsoSstSkew() function computes the skewness of 'sstbox' sstA (usually the skewness of nino3 sstA)

    Author:	Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Yann Planton : yann.planton@locean-ipsl.upmc.fr

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param sstbox: string
        name of box (nino3') for SST
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return skewMetric: dict
        name, value, value_error, units, method, nyears, time_frequency, time_period, ref, keyerror, dive_down_diag

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
    Name = 'ENSO skewness'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Method = 'Standard deviation of ' + sstbox + ' sstA'
    Ref = 'Using CDAT regression calculation'
    metric = 'EnsoSstSkew'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror =\
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, sstbox, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, maskland=True,
                            maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        sstSke, sstSkeErr, dive_down_diag = None, None, {'value': None, 'axis': None}
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        del sst_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape),
                          'time1': '(sst) ' + str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the skewness
        sstSke = float(SkewnessTemporal(sst))

        # Standard Error of the skewness (function of nyears)
        sstSkeErr = sstSke / NUMPYsqrt(yearN)

        # Dive down diagnostic
        sstSke_monthly = SkewMonthly(sst)
        dive_down_diag = {'value': ArrayToList(sstSke_monthly), 'axis': list(sstSke_monthly.getAxis(0)[:])}
        if netcdf is True:
            # additional diagnostic
            # Read file and select the right region
            sst, sst_areacell, unneeded =\
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                    file_area=sstareafile, name_area=sstareaname, file_mask=sstlandmaskfile,
                                    name_mask=sstlandmaskname, maskland=True, maskocean=False, debug=debug, **kwargs)
            # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, unneeded = PreProcessTS(sst, '', areacell=sst_areacell, compute_anom=True, **kwargs)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape),
                              'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 10, **dict_debug)
            # std
            sst = SkewnessTemporal(sst)
            # Regridding
            if 'regridding' not in kwargs.keys():
                kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                        'newgrid_name': 'generic_1x1deg'}
            else:
                if not isinstance(kwargs['regridding'], dict):
                    kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                            'newgrid_name': 'generic_1x1deg'}
            sst = Regrid(sst, None, region='equatorial_pacific_LatExt2', **kwargs['regridding'])
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "monthly standard deviation of " + sstbox + " sstA",
                     'diagnostic_value': sstSke, 'diagnostic_value_error': sstSkeErr}
            dict2 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "standard deviation of " + sstbox + " sstA (whole time series)"}
            dict3 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sstSke_monthly, var1_attributes=dict1, var1_name='sstSke_monthly__' + dataset,
                       var2=sst, var2_attributes=dict2, var2_name='sstSke_map__' + dataset, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstSke), 'line2': 'metric value_error: ' + str(sstSkeErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    skewMetric = {
        'name': Name, 'value': sstSke, 'value_error': sstSkeErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref, 'keyerror': keyerror,
        'dive_down_diag': dive_down_diag,
    }
    return skewMetric


def EnsoSstMap(sstfilemodel, sstnamemodel, sstareafilemodel, sstareanamemodel, sstlandmaskfilemodel,
               sstlandmasknamemodel, sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs,
               sstlandmasknameobs, sstbox, event_definition, centered_rmse=0, biased_rmse=1, debug=False, netcdf=False,
               netcdf_name='', **kwargs):
    """
    The EnsoSstMap() function computes temperature anomalies pattern associated with ENSO on the globe.
    First metric: rmse(observations vs model).
    Second metric: correlation(observations vs model).
    Third metric: std(model)/std(observations)
    These metrics can be used to compute a Taylor diagram.

    Inputs:
    ------
    :param sstfilemodel: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemodel: string
        name of SST variable (tos, ts) in 'sstfilemodel'
    :param sstareafilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST areacell
    :param sstareanamemodel: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafilemodel'
    :param sstlandmaskfilemodel: string, optional
        path_to/filename of the file (NetCDF) of the modeled SST landmask
    :param sstlandmasknamemodel: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfilemodel'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST areacell
    :param sstareanameobs: string, optional
        name of areacell for the SST variable (areacella, areacello,...) in 'sstareafileobs'
    :param sstlandmaskfileobs: string, optional
        path_to/filename of the file (NetCDF) of the observed SST landmask
    :param sstlandmasknameobs: string, optional
        name of landmask for the SST variable (sftlf,...) in 'sstlandmaskfileobs'
    :param sstbox: string
        name of box (e.g. 'global') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return EnsoSstMapMetric: dict
        name, value (rms [obs;model]), value_error, units, method, value2 (corr [obs;model]),
        value_error2, units2, value3 (std_model / std_obs), value_error3, units3, nyears_model, nyears_observations,
        time_frequency, time_period_model, time_period_observations, ref, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'ENSO SSTA pattern '
    Method = region_ev + 'SSTA during ' + season_ev + ' regressed against temperature anomalies in ' + sstbox
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C / C'
    Ref = 'Using CDAT regridding, correlation (centered and biased), std (centered and biased) and ' +\
          'rms (uncentered and biased) calculation'

    # ------------------------------------------------
    # compute ENSO time series
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', 'EnsoSstMap', 10)
        dict_debug = {'file1': '(model) ' + sstfilemodel, 'file2': '(obs) ' + sstfileobs,
                      'var1': '(model) ' + sstnamemodel, 'var2': '(obs) ' + sstnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files ENSO', 10, **dict_debug)
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=region_ev,
                                           time_bounds=kwargs['time_bounds_mod'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=region_ev,
                                         time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if sstareafilemodel:
        model_areacell = ReadAreaSelectRegion(sstareafilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(sstfilemodel, areaname=sstareanamemodel, box=region_ev, **kwargs)
    if sstareafileobs:
        obs_areacell = ReadAreaSelectRegion(sstareafileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(sstfileobs, areaname=sstareanameobs, box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read landmask
    if sstlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(sstlandmaskfilemodel, landmaskname=sstlandmasknamemodel,
                                                  box=region_ev, **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(sstfilemodel, landmaskname=sstlandmasknamemodel, box=region_ev,
                                                  **kwargs)
    if sstlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(sstlandmaskfileobs, landmaskname=sstlandmasknameobs,
                                                  box=region_ev, **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(sstfileobs, landmaskname=sstlandmasknameobs,
                                                box=region_ev, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        sst_model = ApplyLandmask(sst_model, model_landmask, maskland=True, maskocean=False)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=True, maskocean=False)
        del model_landmask
    if obs_landmask is not None:
        sst_obs = ApplyLandmask(sst_obs, obs_landmask, maskland=True, maskocean=False)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=True, maskocean=False)
        del obs_landmask

    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        mini = kwargs['min_time_steps']
        if len(sst_model) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoSstMap: the modeled time-period is too short: "
                            + str(len(sst_model)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)
        if len(sst_obs) < mini:
            list_strings = ["ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": too short time-period",
                            str().ljust(5) + "EnsoSstMap: the observed time-period is too short: "
                            + str(len(sst_obs)) + " (minimum time-period: " + str(mini) + ")"]
            EnsoErrorsWarnings.MyError(list_strings)

    # Number of years
    yearN_model = sst_model.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimeboundsmodel = TimeBounds(sst_model)
    actualtimeboundsobs = TimeBounds(sst_obs)

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, unneeded = PreProcessTS(sst_model, '', areacell=model_areacell, average='horizontal', compute_anom=False,
                                       **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                     **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Seasonal mean and anomalies
    enso_model = SeasonalMean(sst_model, season_ev, compute_anom=True)
    enso_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
    del sst_model, sst_obs
    if season_ev == 'DJF':
        time_ax = enso_model.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
        enso_model.setAxis(0, time_ax)
        del time_ax
        time_ax = enso_obs.getTime()
        time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
        enso_obs.setAxis(0, time_ax)
        del time_ax
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in enso_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in enso_obs.getAxisList()]),
                      'shape1': '(model) ' + str(enso_model.shape), 'shape2': '(obs) ' + str(enso_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(enso_model)), 'time2': '(obs) ' + str(TimeBounds(enso_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

    # ------------------------------------------------
    # regression map
    # ------------------------------------------------
    if debug is True:
        dict_debug = {
            'file1': '(model) ' + sstfilemodel, 'file2': '(obs) ' + sstfileobs, 'var1': '(model) ' + sstnamemodel,
            'var2': '(obs) ' + sstnameobs}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'Files regression map', 10, **dict_debug)
    # Read file and select the right region
    sst_model = ReadSelectRegionCheckUnits(sstfilemodel, sstnamemodel, 'temperature', box=sstbox,
                                          time_bounds=kwargs['time_bounds_mod'], **kwargs)
    sst_obs = ReadSelectRegionCheckUnits(sstfileobs, sstnameobs, 'temperature', box=sstbox,
                                        time_bounds=kwargs['time_bounds_obs'], **kwargs)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadSelectRegionCheckUnits', 15, **dict_debug)
    # Read areacell
    if sstareafilemodel:
        model_areacell = ReadAreaSelectRegion(sstareafilemodel, areaname=sstareanamemodel, box=sstbox, **kwargs)
    else:
        model_areacell = ReadAreaSelectRegion(sstfilemodel, areaname=sstareanamemodel, box=sstbox, **kwargs)
    if sstareafileobs:
        obs_areacell = ReadAreaSelectRegion(sstareafileobs, areaname=sstareanameobs, box=sstbox, **kwargs)
    else:
        obs_areacell = ReadAreaSelectRegion(sstfileobs, areaname=sstareanameobs, box=sstbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_areacell is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_areacell.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_areacell.shape)
        if obs_areacell is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_areacell.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_areacell.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadAreaSelectRegion', 15, **dict_debug)
    # Read if the given region is defined as a land region, an oceanic region, or both
    dict_reg = ReferenceRegions(sstbox)
    if 'maskland' in dict_reg.keys():
        maskland = dict_reg['maskland']
    else:
        maskland = False
    if 'maskocean' in dict_reg.keys():
        maskocean = dict_reg['maskocean']
    else:
        maskocean = False
    # Read landmask
    if sstlandmaskfilemodel:
        model_landmask = ReadLandmaskSelectRegion(sstlandmaskfilemodel, landmaskname=sstlandmasknamemodel, box=sstbox,
                                                  **kwargs)
    else:
        model_landmask = ReadLandmaskSelectRegion(sstfilemodel, landmaskname=sstlandmasknamemodel, box=sstbox, **kwargs)
    if sstlandmaskfileobs:
        obs_landmask = ReadLandmaskSelectRegion(sstlandmaskfileobs, landmaskname=sstlandmasknameobs, box=sstbox,
                                                **kwargs)
    else:
        obs_landmask = ReadLandmaskSelectRegion(sstfileobs, landmaskname=sstlandmasknameobs, box=sstbox, **kwargs)
    if debug is True:
        dict_debug = {}
        if model_landmask is not None:
            dict_debug['axes1'] = '(model) ' + str([ax.id for ax in model_landmask.getAxisList()])
            dict_debug['shape1'] = '(model) ' + str(model_landmask.shape)
        if obs_landmask is not None:
            dict_debug['axes2'] = '(obs) ' + str([ax.id for ax in obs_landmask.getAxisList()])
            dict_debug['shape2'] = '(obs) ' + str(obs_landmask.shape)
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after ReadLandmaskSelectRegion', 15, **dict_debug)
    # Apply landmask
    if model_landmask is not None:
        sst_model = ApplyLandmask(sst_model, model_landmask, maskland=maskland, maskocean=maskocean)
        if model_areacell is None:
            model_areacell = ArrayOnes(model_landmask, id='areacell')
        model_areacell = ApplyLandmaskToArea(model_areacell, model_landmask, maskland=maskland, maskocean=maskocean)
        del model_landmask
    if obs_landmask is not None:
        sst_obs = ApplyLandmask(sst_obs, obs_landmask, maskland=maskland, maskocean=maskocean)
        if obs_areacell is None:
            obs_areacell = ArrayOnes(obs_landmask, id='areacell')
        obs_areacell = ApplyLandmaskToArea(obs_areacell, obs_landmask, maskland=maskland, maskocean=maskocean)
        del obs_landmask

    # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
    sst_model, Method = PreProcessTS(sst_model, Method, areacell=model_areacell, compute_anom=False, **kwargs)
    sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, compute_anom=False, **kwargs)
    del model_areacell, obs_areacell
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

    # Seasonal mean
    sst_model = SeasonalMean(sst_model, season_ev, compute_anom=True)
    sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
    if debug is True:
        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                      'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                      'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                      'time1': '(model) ' + str(TimeBounds(sst_model)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

    # Regridding
    if isinstance(kwargs['regridding'], dict):
        known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                      'regridTool', 'regridMethod'}
        extra_args = set(kwargs['regridding']) - known_args
        if extra_args:
            EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
        sst_model, sst_obs, Method = TwoVarRegrid(sst_model, sst_obs, Method, region=sstbox, **kwargs['regridding'])
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_model.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_model.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

    # regression
    sst_model_slope, sst_model_stderr = LinearRegressionTsAgainstMap(sst_model, enso_model, return_stderr=True)
    sst_obs_slope, sst_obs_stderr = LinearRegressionTsAgainstMap(sst_obs, enso_obs, return_stderr=True)

    # Metric 1
    sstRmse = float(RmsAxis(sst_model_slope, sst_obs_slope, axis='xy', centered=centered_rmse, biased=biased_rmse))
    # Metric 2
    sstCorr = float(Correlation(sst_model_slope, sst_obs_slope, axis='xy', centered=1, biased=1))
    # Metric 3
    std_model = Std(sst_model_slope, weights=None, axis='xy', centered=1, biased=1)
    std_obs = Std(sst_obs_slope, weights=None, axis='xy', centered=1, biased=1)
    sstStd = float(std_model) / float(std_obs)

    # Dive down diagnostic
    dive_down_diag = {'model': sst_model_slope, 'observations': sst_obs_slope,
                      'axisLat': list(sst_model_slope.getAxis(0)[:]), 'axisLon': list(sst_obs_slope.getAxis(1)[:])}

    # Create output
    EnsoSstMapMetric = {
        'name': Name, 'value': sstRmse, 'value_error': None, 'units': Units, 'method': Method,
        'value2': sstCorr, 'value_error2': None, 'units2': '', 'value3': sstStd, 'value_error3': None, 'units3': Units,
        'nyears_model': yearN_model, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimeboundsmodel, 'time_period_observations': actualtimeboundsobs, 'ref': Ref,
        'dive_down_diag': dive_down_diag,
    }
    return EnsoSstMapMetric


def NinaSstDiv(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, box, event_definition,
               dataset='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinaSstDiv() function computes a zonal composite of La Nina events during the peak of the event.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA < 'threshold' during 'season' are considered as La Nina events
        2.) diversity of the zonal location of the minimum SSTA
            2.1) zonal SSTA at the peak of the event is computed for each selected event
            2.2) find the zonal position of the minimum SSTA for each selected event
            2.3) compute the percentage of EP events (minimum SSTA eastward of the given threshold)

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param treshold_ep_ev: float, optional
        see EnsoToolsLib.percentage_val_eastward
        longitude, in degree east, of the westward boundary of eastern Pacific event
        default value is -140°E (i.e., 140°W)
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return NinaDivMetric: dict
        name, value, value_error, units, method, nyears, events, time_frequency, time_period, ref, keyerror,
        dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'treshold_ep_ev',
                    'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nina Diversity (percentage of eastern Pacific La Nina)'
    lat = ReferenceRegions(box)['latitude']
    lon = ReferenceRegions(box)['longitude']
    Method = 'Nina events = ' + region_ev + ' sstA < ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA ' +\
             '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']), westward boundary of EP events' +\
             str(kwargs['treshold_ep_ev']) + 'E'
    Units = '%'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinaSstDiv'

    # ------------------------------------------------
    # 1. detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, areacell, keyerror = \
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, region_ev, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname,
                            maskland=True, maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        ep_event, StdErr, dive_down_diag, event_years = None, None, {'value': None, 'axis': None}, None
    else:
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, unneeded = PreProcessTS(sst, '', areacell=areacell, average='horizontal', compute_anom=False, **kwargs)
        del areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as La Nina events
        # Lists event years
        event_years = DetectEvents(sst, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nina1': 'nbr(' + str(len(event_years)) + '): ' + str(event_years)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years) > 0:
            # ------------------------------------------------
            # 2. diversity of the zonal location of the minimum SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst, areacell, unneeded = \
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, box, file_area=sstareafile,
                                    name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskfile,
                                    maskland=True, maskocean=False, debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, Method = PreProcessTS(sst, Method, areacell=areacell, average=False, compute_anom=False, **kwargs)
            del areacell
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst = SeasonalMean(sst, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder', 'regridTool',
                              'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst = Regrid(sst, None, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                                  'shape1': '(sst) ' + str(sst.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst = AverageMeridional(sst)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sample = Composite_ev_by_ev(sst, event_years, kwargs['frequency'])

            # 2.2 find the zonal position of the minimum SSTA for each selected event
            lon_sstmax = FindXYMinMaxInTs(sample, return_val='mini', smooth=True, axis=0, window=5, method='triangle')
            if debug is True:
                dict_debug = {'line1': 'longitude of the minimum SSTA: ' + str(lon_sstmax)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after FindXYMinMaxInTs', 15, **dict_debug)

            # 2.3 compute the percentage of EP events (minimum SSTA eastward of the given threshold)
            ep_event, keyerror_metric = percentage_val_eastward(lon_sstmax, metric, box,
                                                                threshold=kwargs['treshold_ep_ev'])
            ep_event = float(ep_event)
        else:
            keyerror_metric = "no La Nina event during the period"

        if keyerror_metric is not None:
            StdErr, dive_down_diag = None, {'value': None, 'axis': None}
            keyerror = deepcopy(keyerror_metric)
        else:
            # Standard Error of the Standard Deviation (function of nyears)
            StdErr = None

            # Dive down diagnostic
            dive_down_diag = {'value': ArrayToList(lon_sstmax), 'axis': list(lon_sstmax.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': 'longitude (E)', 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nina_years': str(event_years), 'diagnostic_value_' + dataset: ep_event,
                         'diagnostic_value_error_' + dataset: StdErr}
                dict2 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=lon_sstmax, var1_attributes=dict1,
                           var1_name='Nina_lon_pos_minSSTA__' + dataset, global_attributes=dict2)
                del dict1, dict2

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(ep_event), 'line2': 'metric value_error: ' + str(StdErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinaDivMetric = {
        'name': Name, 'value': ep_event, 'value_error': StdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'events': event_years, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinaDivMetric


def NinaSstDivRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                   sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                   event_definition, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False,
                   netcdf=False, netcdf_name='', **kwargs):
    """
    The NinaSstDivRmse() function computes a zonal minimum of La Nina events during the peak of the event.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA < 'threshold' during 'season' are considered as La Nina events
        2.) diversity of the zonal location of the minimum SSTA
            2.1) zonal SSTA at the peak of the event is computed for each selected event
            2.2) find the zonal position of the minimum SSTA for each selected event and compute a pdf

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return NinaDivMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, events_model, events_observations,
        time_frequency, time_period_model, time_period_observations, ref, keyword, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'PDF of zonal min(SSTA) during Nina'
    lat = ReferenceRegions(box)['latitude']
    lon = ReferenceRegions(box)['longitude']
    Method = 'Nina events = ' + region_ev + ' sstA < ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    Units = 'density'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinaSstDivRmse'

    # ------------------------------------------------
    # 1. detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        pdfRmse, pdfRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                           **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as La Nina events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nina1': '(model) ' + str(event_years_mod), 'nina2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. diversity of the zonal location of the minimum SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst_mod, mod_areacell, unneeded = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                                    name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmaskfilemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            sst_obs, obs_areacell, unneeded = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                                    name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmaskfileobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average=False, compute_anom=False,
                                           **kwargs)
            sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                              'time1': '(model) ' + str(TimeBounds(sst_mod)),
                              'time2': '(obs) ' + str(TimeBounds(sst_obs))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst_mod = SeasonalMean(sst_mod, season_ev, compute_anom=True)
            sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                              'regridTool', 'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst_mod = AverageMeridional(sst_mod)
            sst_obs = AverageMeridional(sst_obs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sample_mod = Composite_ev_by_ev(sst_mod, event_years_mod, kwargs['frequency'])
            sample_obs = Composite_ev_by_ev(sst_obs, event_years_obs, kwargs['frequency'])

            # 2.2 find the zonal position of the minimum SSTA for each selected event and compute a pdf
            # longitude of the minimum SSTA for each selected event
            lon_min_mod = FindXYMinMaxInTs(sample_mod, return_val='mini', smooth=True, axis=0, window=5,
                                           method='triangle')
            lon_min_obs = FindXYMinMaxInTs(sample_obs, return_val='mini', smooth=True, axis=0, window=5,
                                           method='triangle')
            if debug is True:
                dict_debug = {'line1': '(model) longitude  of the maximum SSTA: ' + str(lon_min_mod),
                              'line2': '(obs) longitude  of the maximum SSTA: ' + str(lon_min_obs)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after FindXYMinMaxInTs', 15, **dict_debug)

            # compute PDFs
            if debug is True:
                dict_debug = {'line1': 'lon ' + str(lon) + '  ;  nbr_bins old = ' + str((lon[1] - lon[0]) / 10)
                                       + '  ;  nbr_bins new = ' + str(int((lon[1] - lon[0]) / 10))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'before ComputePDF', 15, **dict_debug)
            pdf_mod = ComputePDF(lon_min_mod, nbr_bins=int((lon[1] - lon[0]) / 10), interval=lon, axis_name='longitude')
            pdf_obs = ComputePDF(lon_min_obs, nbr_bins=int((lon[1] - lon[0]) / 10), interval=lon, axis_name='longitude')

            # Computes the root mean square difference
            pdfRmse = RmsZonal(pdf_mod, pdf_obs, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            pdfRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(pdf_mod), 'observations': ArrayToList(pdf_obs),
                              'axis': list(pdf_mod.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nina_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nina_years': str(event_years_obs)}
                dict3 = {'metric_name': Name, 'metric_value_' + dataset2: pdfRmse,
                         'metric_value_error_' + dataset2: pdfRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=pdf_mod, var1_attributes=dict1, var1_name='pdf__' + dataset1, var2=pdf_obs,
                           var2_attributes=dict2, var2_name='pdf__' + dataset2, global_attributes=dict3)
                del dict1, dict2, dict3
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated La Nina event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed La Nina event during the period"
            else:
                keyerror = "no simulated La Nina event during the period ; no observed La Nina event during the period"
            pdfRmse, pdfRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(pdfRmse), 'line2': 'metric value_error: ' + str(pdfRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)


    # Create output
    NinaDivMetric = {
        'name': Name, 'value': pdfRmse, 'value_error': pdfRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinaDivMetric


def NinaSstDur(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, box, event_definition,
               nbr_years_window, dataset='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinaSstDurRmse() function computes a duration of La Nina events.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA < 'threshold' during 'season' are considered as La Nina events
        2.) La Nina duration
            2.1) get a time series of 2 years before and 2 years after the La Nina peak (4 years time series)
            2.2) count the number of consecutive month bellow a threshold

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return NinaDurMetric: dict
        name, value, value_error, units, method, nyears, events, time_frequency, time_period, time_period, ref,
        keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nina Duration'
    Units = 'months'
    Method = 'Nina events = ' + region_ev + ' sstA < ' + str(threshold) + ' during ' + season_ev + \
             ', number of consecutive months when sstA < -0.5' + Units
    Ref = 'Using CDAT'
    metric = 'NinaSstDur'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror =\
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, region_ev, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, maskland=True,
                            maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        duration_mean, duration_err, dive_down_diag, event_years = None, None, {'value': None, 'axis': None}, None
    else:
        # ------------------------------------------------
        # 1. detect events
        # ------------------------------------------------
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        del sst_areacell
        if debug is True:
            dict_debug = {'axes1': str([ax.id for ax in sst.getAxisList()]), 'shape1': str(sst.shape),
                          'time1': str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as La Nina events
        # Lists event years
        event_years = DetectEvents(sst, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nina1': str(event_years)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years) > 0:
            # ------------------------------------------------
            # 2. La Nina duration
            # ------------------------------------------------
            # 2.1 get a time series of 2 years before and 2 years after the La Nina peak (4 years time series)
            # composites
            sample = Composite_ev_by_ev(sst, event_years, kwargs['frequency'], nbr_years_window=nbr_years_window)
            if debug is True:
                dict_debug = {'axes1': str([ax.id for ax in sample.getAxisList()]), 'shape1': str(sample.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after Composite', 15, **dict_debug)

            # 2.2 count the number of consecutive month bellow a threshold
            duration = DurationAllEvent(sample, -0.5, nino=False, debug=debug)

            duration_err = float(Std(duration) / NUMPYsqrt(len(duration)))
            duration_mean = float(duration.mean())

            # Dive down diagnostic
            dive_down_diag = {'value': ArrayToList(duration), 'axis': list(duration.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nina_years': str(event_years), 'description': "duration of Nina events",
                         'diagnostic_value': duration_mean, 'diagnostic_value_error': duration_err}
                dict2 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=duration, var1_attributes=dict1, var1_name='Nina_duration__' + dataset,
                           global_attributes=dict2)
                del dict1, dict2
        else:
            keyerror = "no La Nina event during the period"
            dive_down_diag = {'value': None, 'axis': None}
            duration_err = None
            duration_mean = None

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(duration_mean),
                      'line2': 'metric value_error: ' + str(duration_err)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinaDurMetric = {
        'name': Name, 'value': duration_mean, 'value_error': duration_err, 'units': Units, 'method': Method,
        'nyears': yearN, 'events': event_years, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds,
        'ref': Ref, 'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinaDurMetric


def NinaSstLonRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                   sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                   event_definition, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False,
                   netcdf=False, netcdf_name='', **kwargs):
    """
    The NinaSstLonRmse() function computes a zonal composite of La Nina events during the peak of the event
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA < 'threshold' during 'season' are considered as La Nina events
    Then the zonal SSTA at the peak of the event is composited for each selected event

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        time_frequency, time_period_model, time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nina Composite Time Series'
    lat = ReferenceRegions(box)['latitude']
    Method = 'Nina events = ' + region_ev + ' sstA < ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinaSstLonRmse'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        compRmse, compRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                           **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as La Nina events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nina1': '(model) ' + str(event_years_mod), 'nina2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. zonal composite of SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst_mod, mod_areacell, unneeded = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                                    name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmaskfilemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            sst_obs, obs_areacell, unneeded = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                                    name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmaskfileobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average=False, compute_anom=False,
                                           **kwargs)
            sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                              'time1': '(model) ' + str(TimeBounds(sst_mod)),
                              'time2': '(obs) ' + str(TimeBounds(sst_obs))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst_mod = SeasonalMean(sst_mod, season_ev, compute_anom=True)
            sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                              'regridTool', 'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst_mod = AverageMeridional(sst_mod)
            sst_obs = AverageMeridional(sst_obs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sst_mod = Composite(sst_mod, event_years_mod, kwargs['frequency'])
            sst_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after Composite', 15, **dict_debug)

            # Computes the root mean square difference
            compRmse = RmsZonal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            compRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                              'axis': list(sst_mod.getAxis(0)[:])}
            if netcdf is True:
                map_mod, mod_areacell, keyerror_mod = \
                    Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                        file_area=sstareafilemod, name_area=sstareanamemod,
                                        file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod, maskland=True,
                                        maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
                map_obs, obs_areacell, keyerror_obs = \
                    Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                        file_area=sstareafileobs, name_area=sstareanameobs,
                                        file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs, maskland=True,
                                        maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
                map_mod, Method = PreProcessTS(map_mod, Method, areacell=mod_areacell, average=False,
                                               compute_anom=False, **kwargs)
                map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                                 **kwargs)
                del mod_areacell, obs_areacell
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape),
                                  'time1': '(model) ' + str(TimeBounds(map_mod)),
                                  'time2': '(obs) ' + str(TimeBounds(map_obs))}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS: netcdf', 15, **dict_debug)
                # Seasonal mean
                map_mod = SeasonalMean(map_mod, season_ev, compute_anom=True)
                map_obs = SeasonalMean(map_obs, season_ev, compute_anom=True)
                # Regridding
                if isinstance(kwargs['regridding'], dict):
                    map_mod, map_obs, unneeded =\
                        TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2', **kwargs['regridding'])
                    if debug is True:
                        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                      'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                      'shape1': '(model) ' + str(map_mod.shape),
                                      'shape2': '(obs) ' + str(map_obs.shape)}
                        EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
                # samples
                map_mod = Composite(map_mod, event_years_mod, kwargs['frequency'])
                map_obs = Composite(map_obs, event_years_obs, kwargs['frequency'])
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nina_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nina_years': str(event_years_obs)}
                dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nina_years': str(event_years_mod)}
                dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nina_years': str(event_years_obs)}
                dict5 = {'metric_name': Name, 'metric_value_' + dataset2: compRmse,
                         'metric_value_error_' + dataset2: compRmseErr, 'metric_method': Method,
                         'metric_reference': Ref, 'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sstComp_lon__' + dataset1,
                           var2=sst_obs, var2_attributes=dict2, var2_name='sstComp_lon__' + dataset2, var3=map_mod,
                           var3_attributes=dict3, var3_name='sstComp_map__' + dataset1, var4=map_obs,
                           var4_attributes=dict4, var4_name='sstComp_map__' + dataset2, global_attributes=dict5)
                del dict1, dict2, dict3, dict4, dict5
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated La Nina event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed La Nina event during the period"
            else:
                keyerror = "no simulated La Nina event during the period ; no observed La Nina event during the period"
            compRmse, compRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(compRmse), 'line2': 'metric value_error: ' + str(compRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinaLonMetric = {
        'name': Name, 'value': compRmse, 'value_error': compRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinaLonMetric


def NinaSstTsRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                  sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs,
                  box, event_definition, nbr_years_window, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='',
                  debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinaSstTsRmse() function computes a time composite of La Nina events
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA < 'threshold' during 'season' are considered as La Nina events
        Then a 'nbr_years_window' long time series centered on selected events is composited for each selected event

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        time_frequency, time_period_model, time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
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
    metric = 'NinaSstTsRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        compRmse, compRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                           **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as La Nina events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=False)
        if debug is True:
            dict_debug = {'nina1': '(model) ' + str(event_years_mod), 'nina2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. temporal composite of SSTA
            # ------------------------------------------------
            # interannual anomalies
            sst_mod = ComputeInterannualAnomalies(sst_mod)
            sst_obs = ComputeInterannualAnomalies(sst_obs)

            # composites
            composite_mod = Composite(sst_mod, event_years_mod, kwargs['frequency'], nbr_years_window=nbr_years_window)
            composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'], nbr_years_window=nbr_years_window)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in composite_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in composite_obs.getAxisList()]),
                              'shape1': '(model) ' + str(composite_mod.shape),
                              'shape2': '(obs) ' + str(composite_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after Composite', 15, **dict_debug)

            # Computes the root mean square difference
            compRmse = RmsAxis(composite_mod, composite_obs, axis=0, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            compRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(composite_mod), 'observations': ArrayToList(composite_obs),
                              'axis': list(composite_mod.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nina_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nina_years': str(event_years_obs)}
                dict3 = {'metric_name': Name, 'metric_value_' + dataset2: compRmse,
                         'metric_value_error_' + dataset2: compRmseErr, 'metric_method': Method,
                         'metric_reference': Ref, 'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=composite_mod, var1_attributes=dict1, var1_name='sst_ts__' + dataset1,
                           var2=composite_obs, var2_attributes=dict2, var2_name='sst_ts__' + dataset2,
                           global_attributes=dict3)
                del dict1, dict2, dict3
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated La Nina event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed La Nina event during the period"
            else:
                keyerror = "no simulated La Nina event during the period ; no observed La Nina event during the period"
            compRmse, compRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(compRmse), 'line2': 'metric value_error: ' + str(compRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinaTsMetric = {
        'name': Name, 'value': compRmse, 'value_error': compRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinaTsMetric


def NinoSstDiv(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, box, event_definition,
               dataset='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinoSstDiv() function computes a zonal composite of El Nino events during the peak of the event.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA > 'threshold' during 'season' are considered as El Nino events
        2.) diversity of the zonal location of the maximum SSTA
            2.1) zonal SSTA at the peak of the event is computed for each selected event
            2.2) find the zonal position of the maximum SSTA for each selected event
            2.3) compute the percentage of EP events (maximum SSTA eastward of the given threshold)

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param treshold_ep_ev: float, optional
        see EnsoToolsLib.percentage_val_eastward
        longitude, in degree east, of the westward boundary of eastern Pacific event
        default value is -140°E (i.e., 140°W)
    :param time_bounds: tuple, optional
        tuple of the first and last dates to extract from the files (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return NinoDivMetric: dict
        name, value, value_error, units, method, nyears, events, time_frequency, time_period, ref, keyerror,
        dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'treshold_ep_ev',
                    'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino Diversity (percentage of eastern Pacific El Nino)'
    lat = ReferenceRegions(box)['latitude']
    lon = ReferenceRegions(box)['longitude']
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA ' +\
             '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']), westward boundary of EP events' +\
             str(kwargs['treshold_ep_ev']) + 'E'
    Units = '%'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinoSstDiv'

    # ------------------------------------------------
    # 1. detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, areacell, keyerror = \
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, region_ev, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname,
                            maskland=True, maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        ep_event, StdErr, dive_down_diag, event_years = None, None, {'value': None, 'axis': None}, None
    else:
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, unneeded = PreProcessTS(sst, '', areacell=areacell, average='horizontal', compute_anom=False, **kwargs)
        del areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as El Nino events
        # Lists event years
        event_years = DetectEvents(sst, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        if debug is True:
            dict_debug = {'nino1': 'nbr(' + str(len(event_years)) + '): ' + str(event_years)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years) > 0:
            # ------------------------------------------------
            # 2. diversity of the zonal location of the maximum SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst, areacell, unneeded = \
                Read_data_mask_area(sstfile, sstname, 'temperature', metric, box, file_area=sstareafile,
                                    name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskfile,
                                    maskland=True, maskocean=False, debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst, Method = PreProcessTS(sst, Method, areacell=areacell, average=False, compute_anom=False, **kwargs)
            del areacell
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape), 'time1': '(sst) ' + str(TimeBounds(sst))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst = SeasonalMean(sst, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder', 'regridTool',
                              'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst = Regrid(sst, None, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                                  'shape1': '(sst) ' + str(sst.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst = AverageMeridional(sst)
            if debug is True:
                dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                              'shape1': '(sst) ' + str(sst.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sample = Composite_ev_by_ev(sst, event_years, kwargs['frequency'])

            # 2.2 find the zonal position of the maximum SSTA for each selected event
            lon_sstmax = FindXYMinMaxInTs(sample, return_val='maxi', smooth=True, axis=0, window=5, method='triangle')
            if debug is True:
                dict_debug = {'line1': 'longitude of the maximum SSTA: ' + str(lon_sstmax)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after FindXYMinMaxInTs', 15, **dict_debug)

            # 2.3 compute the percentage of EP events (maximum SSTA eastward of the given threshold)
            ep_event, keyerror_metric = percentage_val_eastward(lon_sstmax, metric, box,
                                                                threshold=kwargs['treshold_ep_ev'])
            ep_event = float(ep_event)
        else:
            keyerror_metric = "no El Nino event during the period"

        if keyerror_metric is not None:
            StdErr, dive_down_diag = None, {'value': None, 'axis': None}
            keyerror = deepcopy(keyerror_metric)
        else:
            # Standard Error of the Standard Deviation (function of nyears)
            StdErr = None

            # Dive down diagnostic
            dive_down_diag = {'value': ArrayToList(lon_sstmax), 'axis': list(lon_sstmax.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': 'longitude (E)', 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nino_years': str(event_years), 'diagnostic_value_' + dataset: ep_event,
                         'diagnostic_value_error_' + dataset: StdErr}
                dict2 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=lon_sstmax, var1_attributes=dict1,
                           var1_name='Nino_lon_pos_maxSSTA__' + dataset, global_attributes=dict2)
                del dict1, dict2

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(ep_event), 'line2': 'metric value_error: ' + str(StdErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinoDivMetric = {
        'name': Name, 'value': ep_event, 'value_error': StdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'events': event_years, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinoDivMetric


def NinoSstDivRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                   sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                   event_definition, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False,
                   netcdf=False, netcdf_name='', **kwargs):
    """
    The NinoSstDivRmse() function computes a zonal composite of El Nino events during the peak of the event.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA > 'threshold' during 'season' are considered as El Nino events
        2.) diversity of the zonal location of the maximum SSTA
            2.1) zonal SSTA at the peak of the event is computed for each selected event
            2.2) find the zonal position of the maximum SSTA for each selected event and compute a pdf

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return NinoDivMetric: dict
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
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'PDF of zonal max(SSTA) during Nino'
    lat = ReferenceRegions(box)['latitude']
    lon = ReferenceRegions(box)['longitude']
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    Units = 'density'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinoSstDivRmse'

    # ------------------------------------------------
    # 1. detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        pdfRmse, pdfRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as El Nino events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        if debug is True:
            dict_debug = {'nino1': '(model) ' + str(event_years_mod), 'nino2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. diversity of the zonal location of the maximum SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst_mod, mod_areacell, unneeded = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                                    name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmaskfilemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            sst_obs, obs_areacell, unneeded = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                                    name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmaskfileobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average=False, compute_anom=False,
                                           **kwargs)
            sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                              'time1': '(model) ' + str(TimeBounds(sst_mod)),
                              'time2': '(obs) ' + str(TimeBounds(sst_obs))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst_mod = SeasonalMean(sst_mod, season_ev, compute_anom=True)
            sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                              'regridTool', 'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst_mod = AverageMeridional(sst_mod)
            sst_obs = AverageMeridional(sst_obs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sample_mod = Composite_ev_by_ev(sst_mod, event_years_mod, kwargs['frequency'])
            sample_obs = Composite_ev_by_ev(sst_obs, event_years_obs, kwargs['frequency'])

            # 2.2 find the zonal position of the maximum SSTA for each selected event and compute a pdf
            # longitude of the maximum SSTA for each selected event
            lon_min_mod =\
                FindXYMinMaxInTs(sample_mod, return_val='maxi', smooth=True, axis=0, window=5, method='triangle')
            lon_min_obs =\
                FindXYMinMaxInTs(sample_obs, return_val='maxi', smooth=True, axis=0, window=5, method='triangle')
            if debug is True:
                dict_debug = {'line1': '(model) longitude  of the maximum SSTA: ' + str(lon_min_mod),
                              'line2': '(obs) longitude  of the maximum SSTA: ' + str(lon_min_obs)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after FindXYMinMaxInTs', 15, **dict_debug)

            # compute PDFs
            if debug is True:
                dict_debug = {'line1': 'lon ' + str(lon) + '  ;  nbr_bins old = ' + str((lon[1] - lon[0]) / 10)
                                       + '  ;  nbr_bins new = ' + str(int((lon[1] - lon[0]) / 10))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'before ComputePDF', 15, **dict_debug)
            pdf_mod = ComputePDF(lon_min_mod, nbr_bins=int((lon[1] - lon[0]) / 10), interval=lon, axis_name='longitude')
            pdf_obs = ComputePDF(lon_min_obs, nbr_bins=int((lon[1] - lon[0]) / 10), interval=lon, axis_name='longitude')

            # Computes the root mean square difference
            pdfRmse = RmsZonal(pdf_mod, pdf_obs, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            pdfRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(pdf_mod), 'observations': ArrayToList(pdf_obs),
                              'axis': list(pdf_mod.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nino_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nino_years': str(event_years_obs)}
                dict3 = {'metric_name': Name, 'metric_value_' + dataset2: pdfRmse,
                         'metric_value_error_' + dataset2: pdfRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=pdf_mod, var1_attributes=dict1, var1_name='pdf__' + dataset1,
                           var2_attributes=dict2, var2=pdf_obs, var2_name='pdf__' + dataset2, global_attributes=dict3)
                del dict1, dict2, dict3
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated El Nino event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed El Nino event during the period"
            else:
                keyerror = "no simulated El Nino event during the period ; no observed El Nino event during the period"
            pdfRmse, pdfRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(pdfRmse), 'line2': 'metric value_error: ' + str(pdfRmse)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinoDivMetric = {
        'name': Name, 'value': pdfRmse, 'value_error': pdfRmse, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinoDivMetric


def NinoSstDur(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, box, event_definition,
               nbr_years_window, dataset='', debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinoSstDur() function computes a duration of El Nino events.
        1.) detect events
            1.1) SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
            1.2) SSTA > 'threshold' during 'season' are considered as El Nino events
        2.) El Nino duration
            2.1) get a time series of 2 years before and 2 years after the El Nino peak (4 years time series)
            2.2) count the number of consecutive month bellow a threshold

    Inputs:
    ------
    :param sstfile: string
        path_to/filename of the file (NetCDF) of the SST
    :param sstname: string
        name of SST variable (tos, ts) in 'sstfile'
    :param sstareafile: string
        path_to/filename of the file (NetCDF) of the areacell for SST
    :param sstareaname: string
        name of areacell variable (areacella, areacello) in 'sstareafile'
    :param sstlandmaskfile: string
        path_to/filename of the file (NetCDF) of the landmask for SST
    :param sstlandmaskname: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfile'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :return NinoDurMetric: dict
        name, value, value_error, units, method, nyears, events, time_frequency, time_period, time_period, ref

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino Duration'
    Units = 'months'
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + \
             ', number of consecutive months when sstA > 0.5' + Units
    Ref = 'Using CDAT'
    metric = 'NinoSstDur'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror = \
        Read_data_mask_area(sstfile, sstname, 'temperature', metric, region_ev, file_area=sstareafile,
                            name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, maskland=True,
                            maskocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = TimeBounds(sst)

    if keyerror is not None:
        duration_mean, duration_err, dive_down_diag, event_years = None, None, {'value': None, 'axis': None}, None
    else:
        # ------------------------------------------------
        # 1. detect events
        # ------------------------------------------------
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, Method = PreProcessTS(sst, Method, areacell=sst_areacell, average='horizontal', compute_anom=True,
                                   **kwargs)
        del sst_areacell
        if debug is True:
            dict_debug = {'axes1': str([ax.id for ax in sst.getAxisList()]), 'shape1': str(sst.shape),
                          'time1': str(TimeBounds(sst))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as El Nino events
        # Lists event years
        event_years = DetectEvents(sst, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        if debug is True:
            dict_debug = {'nino1': str(event_years)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years) > 0:
            # ------------------------------------------------
            # 2. El Nino duration
            # ------------------------------------------------
            # 2.1 get a time series of N years before and N years after the El Nino peak (2N years time series)
            # composites
            sample = Composite_ev_by_ev(sst, event_years, kwargs['frequency'], nbr_years_window=nbr_years_window)
            if debug is True:
                dict_debug = {'axes1': str([ax.id for ax in sample.getAxisList()]), 'shape1': str(sample.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after Composite', 15, **dict_debug)

            # 2.2 count the number of consecutive month bellow a threshold
            duration = DurationAllEvent(sample, 0.5, nino=True, debug=debug)

            duration_err = float(Std(duration) / NUMPYsqrt(len(duration)))
            duration_mean = float(duration.mean())

            # Dive down diagnostic
            dive_down_diag = {'value': ArrayToList(duration), 'axis': list(duration.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                         'nino_years': str(event_years), 'description': "duration of Nino events",
                         'diagnostic_value': duration_mean, 'diagnostic_value_error': duration_err}
                dict2 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                         'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=duration, var1_attributes=dict1, var1_name='Nino_duration__' + dataset,
                           global_attributes=dict2)
                del dict1, dict2
        else:
            keyerror = "no El Nino event during the period"
            dive_down_diag = {'value': None, 'axis': None}
            duration_err = None
            duration_mean = None

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(duration_mean),
                      'line2': 'metric value_error: ' + str(duration_err)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinoDurMetric = {
        'name': Name, 'value': duration_mean, 'value_error': duration_err, 'units': Units, 'method': Method,
        'nyears': yearN, 'events': event_years, 'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds,
        'ref': Ref, 'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinoDurMetric


def NinoSstLonRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                   sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs, box,
                   event_definition, centered_rmse=0, biased_rmse=1, dataset='', debug=False, netcdf=False,
                   netcdf_name='', **kwargs):
    """
    The NinoSstLonRmse() function computes a zonal composite of El Nino events during the peak of the event
    SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA > 'threshold' during 'season' are considered as El Nino events
    Then the zonal SSTA at the peak of the event is composited for each selected event

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        time_frequency, time_period_model, time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
                    'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'Nino Composite Time Series'
    lat = ReferenceRegions(box)['latitude']
    Method = 'Nino events = ' + region_ev + ' sstA > ' + str(threshold) + ' during ' + season_ev + ', zonal SSTA '\
             + '(meridional averaged [' + str(lat[0]) + ' ; ' + str(lat[1]) + ']'
    if kwargs['normalization']:
        Units = ''
    else:
        Units = 'C'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'NinoSstLonRmse'

    # ------------------------------------------------
    # detect events
    # ------------------------------------------------
    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        compRmse, compRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA < 'threshold' during 'season' are considered as El Nino events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        if debug is True:
            dict_debug = {'nino1': '(model) ' + str(event_years_mod), 'nino2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. zonal composite of SSTA
            # ------------------------------------------------
            # Read file and select the right region
            sst_mod, mod_areacell, unneeded = \
                Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                                    name_area=sstareanamemod, file_mask=sstlandmaskfilemod,
                                    name_mask=sstlandmaskfilemod, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
            sst_obs, obs_areacell, unneeded = \
                Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                                    name_area=sstareanameobs, file_mask=sstlandmaskfileobs,
                                    name_mask=sstlandmaskfileobs, maskland=True, maskocean=False,
                                    time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)

            # 2.1 zonal SSTA at the peak of the event is computed for each selected event
            # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst_mod, Method = PreProcessTS(sst_mod, Method, areacell=mod_areacell, average=False, compute_anom=False,
                                           **kwargs)
            sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                             **kwargs)
            del mod_areacell, obs_areacell
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                              'time1': '(model) ' + str(TimeBounds(sst_mod)),
                              'time2': '(obs) ' + str(TimeBounds(sst_obs))}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

            # Seasonal mean
            sst_mod = SeasonalMean(sst_mod, season_ev, compute_anom=True)
            sst_obs = SeasonalMean(sst_obs, season_ev, compute_anom=True)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after SeasonalMean', 15, **dict_debug)

            # Regridding
            if isinstance(kwargs['regridding'], dict):
                known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                              'regridTool', 'regridMethod'}
                extra_args = set(kwargs['regridding']) - known_args
                if extra_args:
                    EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
                sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

            # Meridional average
            sst_mod = AverageMeridional(sst_mod)
            sst_obs = AverageMeridional(sst_obs)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

            # samples
            sst_mod = Composite(sst_mod, event_years_mod, kwargs['frequency'])
            sst_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'])

            # Computes the root mean square difference
            compRmse = RmsZonal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            compRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                              'axis': list(sst_mod.getAxis(0)[:])}
            if netcdf is True:
                map_mod, mod_areacell, keyerror_mod = \
                    Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                        file_area=sstareafilemod, name_area=sstareanamemod,
                                        file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod, maskland=True,
                                        maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug, **kwargs)
                map_obs, obs_areacell, keyerror_obs = \
                    Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, 'equatorial_pacific_LatExt2',
                                        file_area=sstareafileobs, name_area=sstareanameobs,
                                        file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs, maskland=True,
                                        maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug, **kwargs)
                map_mod, Method = PreProcessTS(map_mod, Method, areacell=mod_areacell, average=False,
                                               compute_anom=False, **kwargs)
                map_obs, unneeded = PreProcessTS(map_obs, '', areacell=obs_areacell, average=False, compute_anom=False,
                                                 **kwargs)
                del mod_areacell, obs_areacell
                if debug is True:
                    dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                  'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                  'shape1': '(model) ' + str(map_mod.shape), 'shape2': '(obs) ' + str(map_obs.shape),
                                  'time1': '(model) ' + str(TimeBounds(map_mod)),
                                  'time2': '(obs) ' + str(TimeBounds(map_obs))}
                    EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS: netcdf', 15, **dict_debug)
                # Seasonal mean
                map_mod = SeasonalMean(map_mod, season_ev, compute_anom=True)
                map_obs = SeasonalMean(map_obs, season_ev, compute_anom=True)
                # Regridding
                if isinstance(kwargs['regridding'], dict):
                    map_mod, map_obs, unneeded =\
                        TwoVarRegrid(map_mod, map_obs, '', region='equatorial_pacific_LatExt2', **kwargs['regridding'])
                    if debug is True:
                        dict_debug = {'axes1': '(model) ' + str([ax.id for ax in map_mod.getAxisList()]),
                                      'axes2': '(obs) ' + str([ax.id for ax in map_obs.getAxisList()]),
                                      'shape1': '(model) ' + str(map_mod.shape),
                                      'shape2': '(obs) ' + str(map_obs.shape)}
                        EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid: netcdf', 15, **dict_debug)
                # samples
                map_mod = Composite(map_mod, event_years_mod, kwargs['frequency'])
                map_obs = Composite(map_obs, event_years_obs, kwargs['frequency'])
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nino_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nino_years': str(event_years_obs)}
                dict3 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nino_years': str(event_years_mod)}
                dict4 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nino_years': str(event_years_obs)}
                dict5 = {'metric_name': Name, 'metric_value_' + dataset2: compRmse,
                         'metric_value_error_' + dataset2: compRmseErr, 'metric_method': Method,
                         'metric_reference': Ref, 'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sstComp_lon__' + dataset1,
                           var2=sst_obs, var2_attributes=dict2, var2_name='sstComp_lon__' + dataset2, var3=map_mod,
                           var3_attributes=dict3, var3_name='sstComp_map__' + dataset1, var4=map_obs, var4_attributes=dict4,
                           var4_name='sstComp_map__' + dataset2, global_attributes=dict5)
                del dict1, dict2, dict3, dict4, dict5
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated El Nino event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed El Nino event during the period"
            else:
                keyerror = "no simulated El Nino event during the period ; no observed El Nino event during the period"
            compRmse, compRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(compRmse), 'line2': 'metric value_error: ' + str(compRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinoLonMetric = {
        'name': Name, 'value': compRmse, 'value_error': compRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinoLonMetric


def NinoSstTsRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                  sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs,
                  box, event_definition, nbr_years_window, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='',
                  debug=False, netcdf=False, netcdf_name='', **kwargs):
    """
    The NinoSstTsRmse() function computes a time composite of El Nino events
    SSTA averaged in 'box' are normalized / detrended / smoothed (running average) if applicable
        Then SSTA > 'threshold' during 'season' are considered as El Nino events
        Then a 'nbr_years_window' long time series centered on selected events is composited for each selected event

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3') for SST
    :param event_definition: dict
        dictionary providing the necessary information to detect ENSO events (region_ev, season_ev, threshold)
        e.g., event_definition = {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75}
    :param nbr_years_window: integer
        number of years used to compute the composite (e.g. 6)
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
        time_frequency, time_period_model, time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    """
    # setting variables
    region_ev = event_definition['region_ev']
    season_ev = event_definition['season_ev']
    threshold = event_definition['threshold']
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'smoothing', 'time_bounds_mod',
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
    metric = 'NinoSstTsRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, region_ev, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, region_ev, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        compRmse, compRmseErr, event_years_mod, event_years_obs = None, None, None, None
        dive_down_diag = {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # 1.1 SSTA averaged in 'region_ev' are normalized / detrended / smoothed (running average) if applicable
        # Preprocess sst (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst_mod, unneeded = PreProcessTS(sst_mod, '', areacell=mod_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', areacell=obs_areacell, average='horizontal', compute_anom=False,
                                         **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape),
                          'time1': '(model) ' + str(TimeBounds(sst_mod)), 'time2': '(obs) ' + str(TimeBounds(sst_obs))}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # 1.2 SSTA > 'threshold' during 'season' are considered as El Nino events
        # Lists event years
        event_years_mod = DetectEvents(sst_mod, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        event_years_obs = DetectEvents(sst_obs, season_ev, threshold, normalization=kwargs['normalization'], nino=True)
        if debug is True:
            dict_debug = {'nino1': '(model) ' + str(event_years_mod), 'nino2': '(obs) ' + str(event_years_obs)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after DetectEvents', 15, **dict_debug)

        if len(event_years_mod) > 0 and len(event_years_obs) > 0:
            # ------------------------------------------------
            # 2. temporal composite of SSTA
            # ------------------------------------------------
            # interannual anomalies
            sst_mod = ComputeInterannualAnomalies(sst_mod)
            sst_obs = ComputeInterannualAnomalies(sst_obs)

            # composites
            composite_mod = Composite(sst_mod, event_years_mod, kwargs['frequency'], nbr_years_window=nbr_years_window)
            composite_obs = Composite(sst_obs, event_years_obs, kwargs['frequency'], nbr_years_window=nbr_years_window)
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in composite_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in composite_obs.getAxisList()]),
                              'shape1': '(model) ' + str(composite_mod.shape),
                              'shape2': '(obs) ' + str(composite_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after Composite', 15, **dict_debug)

            # Computes the root mean square difference
            compRmse = RmsAxis(composite_mod, composite_obs, axis=0, centered=centered_rmse, biased=biased_rmse)

            # Error on the metric
            compRmseErr = None

            # Dive down diagnostic
            dive_down_diag = {'model': ArrayToList(composite_mod), 'observations': ArrayToList(composite_obs),
                              'axis': list(composite_mod.getAxis(0)[:])}
            if netcdf is True:
                if ".nc" in netcdf_name:
                    file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
                else:
                    file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
                dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod),
                         'nino_years': str(event_years_mod)}
                dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs),
                         'nino_years': str(event_years_obs)}
                dict3 = {'metric_name': Name, 'metric_value_' + dataset2: compRmse,
                         'metric_value_error_' + dataset2: compRmseErr, 'metric_method': Method,
                         'metric_reference': Ref, 'frequency': kwargs['frequency']}
                SaveNetcdf(file_name, var1=composite_mod, var1_attributes=dict1, var1_name='sst_ts__' + dataset1,
                           var2=composite_obs, var2_attributes=dict2, var2_name='sst_ts__' + dataset2,
                           global_attributes=dict3)
                del dict1, dict2, dict3
        else:
            if len(event_years_mod) == 0 and len(event_years_obs) > 0:
                keyerror = "no simulated El Nino event during the period"
            elif len(event_years_mod) > 0 and len(event_years_obs) == 0:
                keyerror = "no observed El Nino event during the period"
            else:
                keyerror = "no simulated El Nino event during the period ; no observed El Nino event during the period"
            compRmse, compRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}

    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(compRmse), 'line2': 'metric value_error: ' + str(compRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    NinoTsMetric = {
        'name': Name, 'value': compRmse, 'value_error': compRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'events_model': event_years_mod,
        'events_observations': event_years_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return NinoTsMetric


def SeasonalPrLatRmse(prfilemod, prnamemod, prareafilemod, prareanamemod, prlandmaskfilemod, prlandmasknamemod,
                      prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs, prlandmasknameobs, box,
                      centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                      netcdf_name='', **kwargs):
    """
    The SeasonalPrLatRmse() function computes the climatological (12 months) PR (precipitation) meridional (latitude)
    standard deviation root mean square error (RMSE) in a 'box' (usually the nino3.3_LatExt)

    Inputs:
    ------
    :param prfilemod: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemod: string
        name of PR variable (pr, precip) in 'prfilemod'
    :param prareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for PR
    :param prareanamemod: string
        name of areacell variable (areacella, areacello) in 'prareafilemod'
    :param prlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for PR
    :param prlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfilemod'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param prareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for PR
    :param prareanameobs: string
        name of areacell variable (areacella, areacello) in 'prareafileobs'
    :param prlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for PR
    :param prlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfileobs'
    :param box: string
        name of box ('nino3.3_LatExt') for PR
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LatRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'pr meridional seasonality RMSE'
    Units = 'mm/day'
    Method = 'Meridional root mean square error of ' + box + ' climatological pr STD'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'SeasonalPrLatRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    pr_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, box, file_area=prareafilemod,
                            name_area=prareanamemod, file_mask=prlandmaskfilemod, name_mask=prlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    pr_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, box, file_area=prareafileobs,
                            name_area=prareanameobs, file_mask=prlandmaskfileobs, name_mask=prlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = pr_mod.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(pr_mod)
    actualtimebounds_obs = TimeBounds(pr_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        prRmse, prRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        pr_mod, Method = PreProcessTS(pr_mod, Method, compute_sea_cycle=True, **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', compute_sea_cycle=True, **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # standard deviation computation
        pr_mod = Std(pr_mod)
        pr_obs = Std(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after Std', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            pr_mod, pr_obs, Method = TwoVarRegrid(pr_mod, pr_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                              'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Zonal average
        pr_mod = AverageZonal(pr_mod)
        pr_obs = AverageZonal(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageZonal', 15, **dict_debug)

        # Computes the root mean square difference
        prRmse = RmsMeridional(pr_mod, pr_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        prRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(pr_mod), 'observations': ArrayToList(pr_obs),
                          'axis': list(pr_mod.getAxis(0)[:])}
        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: prRmse,
                     'metric_value_error_' + dataset2: prRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=pr_mod, var1_attributes=dict1, var1_name='pr__' + dataset1, var2=pr_obs,
                       var2_attributes=dict2, var2_name='pr__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(prRmse), 'line2': 'metric value_error: ' + str(prRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': prRmse, 'value_error': prRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LatRmseMetric


def SeasonalPrLonRmse(prfilemod, prnamemod, prareafilemod, prareanamemod, prlandmaskfilemod, prlandmasknamemod,
                      prfileobs, prnameobs, prareafileobs, prareanameobs, prlandmaskfileobs, prlandmasknameobs, box,
                      centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                      netcdf_name='', **kwargs):
    """
    The SeasonalPrLonRmse() function computes the climatological (12 months) PR (precipitation) zonal (longitude)
    standard deviation root mean square error (RMSE) in a 'box' (usually the Equatorial Pacific)

    Inputs:
    ------
     :param prfilemod: string
        path_to/filename of the file (NetCDF) of the modeled PR
    :param prnamemod: string
        name of PR variable (pr, precip) in 'prfilemod'
    :param prareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for PR
    :param prareanamemod: string
        name of areacell variable (areacella, areacello) in 'prareafilemod'
    :param prlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for PR
    :param prlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfilemod'
    :param prfileobs: string
        path_to/filename of the file (NetCDF) of the observed PR
    :param prnameobs: string
        name of PR variable (pr, prec) in 'prfileobs'
    :param prareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for PR
    :param prareanameobs: string
        name of areacell variable (areacella, areacello) in 'prareafileobs'
    :param prlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for PR
    :param prlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'prlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for PR
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
        tuple of the first and last dates to extract from the modeled PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None
    :param time_bounds_obs: tuple, optional
        tuple of the first and last dates to extract from the observed PR file (strings)
        e.g., time_bounds=('1979-01-01T00:00:00', '2017-01-01T00:00:00')
        default value is None

    Output:
    ------
    :return LonRmseMetric: dict
        name, value, value_error, units, method, nyears_model, nyears_observations, time_frequency, time_period_model,
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'pr zonal seasonality RMSE'
    Units = 'mm/day'
    Method = 'Zonal root mean square error of ' + box + ' climatological pr STD'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'SeasonalPrLonRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    pr_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(prfilemod, prnamemod, 'precipitations', metric, box, file_area=prareafilemod,
                            name_area=prareanamemod, file_mask=prlandmaskfilemod, name_mask=prlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    pr_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(prfileobs, prnameobs, 'precipitations', metric, box, file_area=prareafileobs,
                            name_area=prareanameobs, file_mask=prlandmaskfileobs, name_mask=prlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = pr_mod.shape[0] / 12
    yearN_obs = pr_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(pr_mod)
    actualtimebounds_obs = TimeBounds(pr_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        prRmse, prRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        pr_mod, Method = PreProcessTS(pr_mod, Method, compute_sea_cycle=True, **kwargs)
        pr_obs, unneeded = PreProcessTS(pr_obs, '', compute_sea_cycle=True, **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # standard deviation computation
        pr_mod = Std(pr_mod)
        pr_obs = Std(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after Std', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            pr_mod, pr_obs, Method = TwoVarRegrid(pr_mod, pr_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                              'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        pr_mod = AverageMeridional(pr_mod)
        pr_obs = AverageMeridional(pr_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in pr_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in pr_obs.getAxisList()]),
                          'shape1': '(model) ' + str(pr_mod.shape), 'shape2': '(obs) ' + str(pr_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # Computes the root mean square difference
        prRmse = RmsZonal(pr_mod, pr_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        prRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(pr_mod), 'observations': ArrayToList(pr_obs),
                          'axis': list(pr_mod.getAxis(0)[:])}
        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: prRmse,
                     'metric_value_error_' + dataset2: prRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=pr_mod, var1_attributes=dict1, var1_name='pr__' + dataset1, var2=pr_obs,
                       var2_attributes=dict2, var2_name='pr__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(prRmse), 'line2': 'metric value_error: ' + str(prRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': prRmse, 'value_error': prRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric


def SeasonalSstLatRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                       sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs,
                       box, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                       netcdf_name='', **kwargs):
    """
    The SeasonalSstLatRmse() function computes the climatological (12 months) SST meridional (latitude) standard
    deviation root mean square error (RMSE) in a 'box' (usually the nino3.3_LatExt)

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('nino3.3_LatExt') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try: kwargs[arg]
        except: kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sst meridional seasonality RMSE'
    Units = 'C'
    Method = 'Meridional root mean square error of ' + box + ' climatological sst STD'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'SeasonalSstLatRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        sstRmse, sstRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, compute_sea_cycle=True, **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', compute_sea_cycle=True, **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # standard deviation computation
        sst_mod = Std(sst_mod)
        sst_obs = Std(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after Std', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Zonal average
        sst_mod = AverageZonal(sst_mod)
        sst_obs = AverageZonal(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageZonal', 15, **dict_debug)

        # Computes the root mean square difference
        sstRmse = RmsMeridional(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        sstRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                          'axis': list(sst_mod.getAxis(0)[:])}
        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: sstRmse,
                     'metric_value_error_' + dataset2: sstRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sst__' + dataset1, var2=sst_obs,
                       var2_attributes=dict2, var2_name='sst__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstRmse), 'line2': 'metric value_error: ' + str(sstRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LatRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': sstRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model':actualtimebounds_mod, 'time_period_observations':actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LatRmseMetric


def SeasonalSstLonRmse(sstfilemod, sstnamemod, sstareafilemod, sstareanamemod, sstlandmaskfilemod, sstlandmasknamemod,
                       sstfileobs, sstnameobs, sstareafileobs, sstareanameobs, sstlandmaskfileobs, sstlandmasknameobs,
                       box, centered_rmse=0, biased_rmse=1, dataset1='', dataset2='', debug=False, netcdf=False,
                       netcdf_name='', **kwargs):
    """
    The SeasonalSstLonRmse() function computes the climatological (12 months) SST zonal (longitude) standard
    deviation root mean square error (RMSE) in a 'box' (usually the Equatorial Pacific)

    Inputs:
    ------
    :param sstfilemod: string
        path_to/filename of the file (NetCDF) of the modeled SST
    :param sstnamemod: string
        name of SST variable (tos, ts) in 'sstfilemod'
    :param sstareafilemod: string
        path_to/filename of the file (NetCDF) of the model areacell for SST
    :param sstareanamemod: string
        name of areacell variable (areacella, areacello) in 'sstareafilemod'
    :param sstlandmaskfilemod: string
        path_to/filename of the file (NetCDF) of the model landmask for SST
    :param sstlandmasknamemod: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfilemod'
    :param sstfileobs: string
        path_to/filename of the file (NetCDF) of the observed SST
    :param sstnameobs: string
        name of SST variable (tos, ts) in 'sstfileobs'
    :param sstareafileobs: string
        path_to/filename of the file (NetCDF) of the observations areacell for SST
    :param sstareanameobs: string
        name of areacell variable (areacella, areacello) in 'sstareafileobs'
    :param sstlandmaskfileobs: string
        path_to/filename of the file (NetCDF) of the observations landmask for SST
    :param sstlandmasknameobs: string
        name of landmask variable (sftlf, lsmask, landmask) in 'sstlandmaskfileobs'
    :param box: string
        name of box ('equatorial_pacific') for SST
    :param centered_rmse: int, optional
        default value = 0 returns uncentered statistic (same as None). To remove the mean first (i.e centered statistic)
        set to 1. NOTE: Most other statistic functions return a centered statistic by default
    :param biased_rmse: int, optional
        default value = 1 returns biased statistic (number of elements along given axis)
        If want to compute an unbiased variance pass anything but 1 (number of elements along given axis minus 1)
    :param dataset: string, optional
        name of current dataset (e.g., 'model', 'obs', ...)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :param netcdf: boolean, optional
        default value = False dive_down are not saved in NetCDFs
        If you want to save the dive down diagnostics set it to True
    :param netcdf_name: string, optional
        default value = '' NetCDFs are saved where the program is ran without a root name
        the name of a metric will be append at the end of the root name
        e.g., netcdf_name='/path/to/directory/USER_DATE_METRICCOLLECTION_MODEL'
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
    :param time_bounds_mod: tuple, optional
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
        time_period_observations, ref, keyerror, dive_down_diag

    Method:
    -------
        uses tools from uvcdat library

    Notes:
    -----
        TODO: add error calculation to rmse (function of nyears)

    """
    # test given kwargs
    needed_kwarg = ['detrending', 'frequency', 'min_time_steps', 'normalization', 'regridding', 'smoothing',
                    'time_bounds_mod', 'time_bounds_obs']
    for arg in needed_kwarg:
        try:
            kwargs[arg]
        except:
            kwargs[arg] = DefaultArgValues(arg)

    # Define metric attributes
    Name = 'sst zonal seasonality RMSE'
    Units = 'C'
    Method = 'Zonal root mean square error of ' + box + ' climatological sst STD'
    Ref = 'Using CDAT regridding and rms (uncentered and biased) calculation'
    metric = 'SeasonalSstLonRmse'

    # Read file and select the right region
    if debug is True:
        EnsoErrorsWarnings.DebugMode('\033[92m', metric, 10)
    sst_mod, mod_areacell, keyerror_mod = \
        Read_data_mask_area(sstfilemod, sstnamemod, 'temperature', metric, box, file_area=sstareafilemod,
                            name_area=sstareanamemod, file_mask=sstlandmaskfilemod, name_mask=sstlandmasknamemod,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_mod'], debug=debug,
                            **kwargs)
    sst_obs, obs_areacell, keyerror_obs = \
        Read_data_mask_area(sstfileobs, sstnameobs, 'temperature', metric, box, file_area=sstareafileobs,
                            name_area=sstareanameobs, file_mask=sstlandmaskfileobs, name_mask=sstlandmasknameobs,
                            maskland=True, maskocean=False, time_bounds=kwargs['time_bounds_obs'], debug=debug,
                            **kwargs)

    # Number of years
    yearN_mod = sst_mod.shape[0] / 12
    yearN_obs = sst_obs.shape[0] / 12

    # Time period
    actualtimebounds_mod = TimeBounds(sst_mod)
    actualtimebounds_obs = TimeBounds(sst_obs)

    if keyerror_mod is not None or keyerror_obs is not None:
        sstRmse, sstRmseErr, dive_down_diag = None, None, {'model': None, 'observations': None, 'axis': None}
        keyerror = ''
        if keyerror_mod is not None:
            keyerror = keyerror_mod
        if len(keyerror) > 0 and keyerror_obs is not None:
            keyerror += " ; "
        if keyerror_obs is not None:
            keyerror += keyerror_obs
    else:
        keyerror = None
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        # here only the detrending (if applicable) and time averaging are performed
        sst_mod, Method = PreProcessTS(sst_mod, Method, compute_sea_cycle=True, **kwargs)
        sst_obs, unneeded = PreProcessTS(sst_obs, '', compute_sea_cycle=True, **kwargs)
        del mod_areacell, obs_areacell
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # standard deviation computation
        sst_mod = Std(sst_mod)
        sst_obs = Std(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after Std', 15, **dict_debug)

        # Regridding
        if isinstance(kwargs['regridding'], dict):
            known_args = {'model_orand_obs', 'newgrid', 'missing', 'order', 'mask', 'newgrid_name', 'regridder',
                          'regridTool', 'regridMethod'}
            extra_args = set(kwargs['regridding']) - known_args
            if extra_args:
                EnsoErrorsWarnings.UnknownKeyArg(extra_args, INSPECTstack())
            sst_mod, sst_obs, Method = TwoVarRegrid(sst_mod, sst_obs, Method, region=box, **kwargs['regridding'])
            if debug is True:
                dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                              'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                              'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
                EnsoErrorsWarnings.DebugMode('\033[92m', 'after TwoVarRegrid', 15, **dict_debug)

        # Meridional average
        sst_mod = AverageMeridional(sst_mod)
        sst_obs = AverageMeridional(sst_obs)
        if debug is True:
            dict_debug = {'axes1': '(model) ' + str([ax.id for ax in sst_mod.getAxisList()]),
                          'axes2': '(obs) ' + str([ax.id for ax in sst_obs.getAxisList()]),
                          'shape1': '(model) ' + str(sst_mod.shape), 'shape2': '(obs) ' + str(sst_obs.shape)}
            EnsoErrorsWarnings.DebugMode('\033[92m', 'after AverageMeridional', 15, **dict_debug)

        # Computes the root mean square difference
        sstRmse = RmsZonal(sst_mod, sst_obs, centered=centered_rmse, biased=biased_rmse)

        # Error on the metric
        sstRmseErr = None

        # Dive down diagnostic
        dive_down_diag = {'model': ArrayToList(sst_mod), 'observations': ArrayToList(sst_obs),
                          'axis': list(sst_mod.getAxis(0)[:])}
        if netcdf is True:
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metric + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metric + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN_mod, 'time_period': str(actualtimebounds_mod)}
            dict2 = {'units': Units, 'number_of_years_used': yearN_obs, 'time_period': str(actualtimebounds_obs)}
            dict3 = {'metric_name': Name, 'metric_value_' + dataset2: sstRmse,
                     'metric_value_error_' + dataset2: sstRmseErr, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            SaveNetcdf(file_name, var1=sst_mod, var1_attributes=dict1, var1_name='sst__' + dataset1, var2=sst_obs,
                       var2_attributes=dict2, var2_name='sst__' + dataset2, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstRmse), 'line2': 'metric value_error: ' + str(sstRmseErr)}
        EnsoErrorsWarnings.DebugMode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    LonRmseMetric = {
        'name': Name, 'value': sstRmse, 'value_error': sstRmseErr, 'units': Units, 'method': Method,
        'nyears_model': yearN_mod, 'nyears_observations': yearN_obs, 'time_frequency': kwargs['frequency'],
        'time_period_model': actualtimebounds_mod, 'time_period_observations': actualtimebounds_obs, 'ref': Ref,
        'keyerror': keyerror, 'dive_down_diag': dive_down_diag,
    }
    return LonRmseMetric
# ---------------------------------------------------------------------------------------------------------------------#




