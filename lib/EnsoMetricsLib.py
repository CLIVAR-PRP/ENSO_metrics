# -*- coding:UTF-8 -*-
from __future__ import print_function
from copy import deepcopy
from numpy import sqrt as NUMPYsqrt

# ENSO_metrics package functions:
import EnsoErrorsWarnings
from EnsoCdatToolsLib import average_meridional, pre_process_data, read_data, compute_regrid, save_netcdf, compute_std,\
    get_time_bounds
from KeyArgLib import DefaultArgValues


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to compute ENSO metrics
# These functions have file names and variable names as inputs and metric as output
#
def EnsoAmpl(sstfile, sstname, sstareafile, sstareaname, sstlandmaskfile, sstlandmaskname, sstbox, dataset='',
             debug=False, netcdf=False, netcdf_name='', metname='EnsoAmpl', **kwargs):
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
    :param metname: string, optional
        default value = 'EnsoAmpl', name of the metric
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
        EnsoErrorsWarnings.debug_mode('\033[92m', metric, 10)
    sst, sst_areacell, keyerror =\
        read_data(sstfile, sstname, 'temperature', metric, sstbox, file_area=sstareafile,
                  name_area=sstareaname, file_mask=sstlandmaskfile, name_mask=sstlandmaskname, mask_land=True,
                  mask_ocean=False, debug=debug, **kwargs)

    # Number of years
    yearN = sst.shape[0] / 12

    # Time period
    actualtimebounds = get_time_bounds(sst)

    if keyerror is not None:
        sstStd, sstStdErr = None, None
    else:
        # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
        sst, Method = pre_process_data(sst, Method, areacell=sst_areacell, average='horizontal', interannual_anomalies=True,
                                       region=sstbox, **kwargs)
        del sst_areacell
        if debug is True:
            dict_debug = {'axes1': '(sst) ' + str([ax.id for ax in sst.getAxisList()]),
                          'shape1': '(sst) ' + str(sst.shape),
                          'time1': '(sst) ' + str(get_time_bounds(sst))}
            EnsoErrorsWarnings.debug_mode('\033[92m', 'after PreProcessTS', 15, **dict_debug)

        # Computes the standard deviation
        sstStd = float(compute_std(sst))

        # Standard Error of the Standard Deviation (function of nyears)
        sstStdErr = sstStd / NUMPYsqrt(yearN)

        # Dive down diagnostic
        if netcdf is True:
            # additional diagnostic
            # Read file and select the right region
            sst1, sst_areacell1, unneeded =\
                read_data(sstfile, sstname, 'temperature', metric, 'equatorial_pacific_LatExt2',
                          file_area=sstareafile, name_area=sstareaname, file_mask=sstlandmaskfile,
                          name_mask=sstlandmaskname, mask_land=True, mask_ocean=False, debug=debug, **kwargs)
            sst2, sst_areacell2, unneeded = \
                read_data(sstfile, sstname, 'temperature', metric, 'equatorial_pacific',
                          file_area=sstareafile, name_area=sstareaname, file_mask=sstlandmaskfile,
                          name_mask=sstlandmaskname, mask_land=True, mask_ocean=False, debug=debug, **kwargs)
            # Preprocess variables (computes anomalies, normalizes, detrends TS, smoothes TS, averages horizontally)
            sst1, unneeded = pre_process_data(sst1, '', areacell=sst_areacell1, interannual_anomalies=True,
                                              region="equatorial_pacific_LatExt2", **kwargs)
            sst2, unneeded = pre_process_data(sst2, '', areacell=sst_areacell2, interannual_anomalies=True,
                                              region="equatorial_pacific", **kwargs)
            del sst_areacell1, sst_areacell2
            if debug is True:
                dict_debug = {'axes1': '(sst1) ' + str([ax.id for ax in sst1.getAxisList()]),
                              'axes2': '(sst2) ' + str([ax.id for ax in sst2.getAxisList()]),
                              'shape1': '(sst1) ' + str(sst1.shape), 'shape2': '(sst2) ' + str(sst2.shape),
                              'time1': '(sst1) ' + str(get_time_bounds(sst1)), 'time2': '(sst2) ' + str(get_time_bounds(sst2))}
                EnsoErrorsWarnings.debug_mode('\033[92m', 'after PreProcessTS', 10, **dict_debug)
            # std
            sst1 = compute_std(sst1)
            sst2 = compute_std(sst2)
            # Regridding
            if 'regridding' not in kwargs.keys():
                kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                        'newgrid_name': 'generic_1x1deg'}
            else:
                if not isinstance(kwargs['regridding'], dict):
                    kwargs['regridding'] = {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                            'newgrid_name': 'generic_1x1deg'}
            sst1 = compute_regrid(sst1, None, region='equatorial_pacific_LatExt2', **kwargs['regridding'])
            sst2 = compute_regrid(sst2, None, region='equatorial_pacific', **kwargs['regridding'])
            # Meridional average
            sst2 = average_meridional(sst2)
            # Dive down diagnostic
            if ".nc" in netcdf_name:
                file_name = deepcopy(netcdf_name).replace(".nc", "_" + metname + ".nc")
            else:
                file_name = deepcopy(netcdf_name) + "_" + metname + ".nc"
            dict1 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "zonal standard deviation of equatorial_pacific sstA",
                     'diagnostic_value': sstStd, 'diagnostic_value_error': sstStdErr}
            dict2 = {'units': Units, 'number_of_years_used': yearN, 'time_period': str(actualtimebounds),
                     'description': "standard deviation of equatorial_pacific sstA"}
            dict3 = {'metric_name': Name, 'metric_method': Method, 'metric_reference': Ref,
                     'frequency': kwargs['frequency']}
            save_netcdf(file_name, var1=sst2, var1_attributes=dict1, var1_name='sstStd_lon__' + dataset,
                        var2=sst1, var2_attributes=dict2, var2_name='sstStd_map__' + dataset, global_attributes=dict3)
            del dict1, dict2, dict3
    # metric value
    if debug is True:
        dict_debug = {'line1': 'metric value: ' + str(sstStd), 'line2': 'metric value_error: ' + str(sstStdErr)}
        EnsoErrorsWarnings.debug_mode('\033[92m', 'end of ' + metric, 10, **dict_debug)

    # Create output
    amplMetric = {
        'name': Name, 'value': sstStd, 'value_error': sstStdErr, 'units': Units, 'method': Method, 'nyears': yearN,
        'time_frequency': kwargs['frequency'], 'time_period': actualtimebounds, 'ref': Ref, 'keyerror': keyerror,
    }
    return amplMetric

