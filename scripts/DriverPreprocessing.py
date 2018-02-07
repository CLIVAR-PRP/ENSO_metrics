# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack

# ENSO_metrics package 'driver' functions:
from DriverPreprocessingUvcdatToolsLib import ComputeAnomalies, ComputeAverage, ComputeDetrend, ComputeNormalize,\
    ComputeRegrid, ComputeSmooth, ReadAndSelectRegion

# ENSO_metrics package functions:
import EnsoErrorsWarnings


dict_processes = {
    'anomalies': ComputeAnomalies, 'averaging': ComputeAverage, 'detrending': ComputeDetrend,
    'normalization': ComputeNormalize, 'regridding': ComputeRegrid, 'smoothing': ComputeSmooth,
}


def preprocess(file_name, var_name, preprocessing, region=None, period=None, frequency=None, units=None, **keyarg):
    info = 'Preprocessing:'
    # read 'var_name' in 'file_name' in the given 'region' and 'period'
    tab, info = ReadAndSelectRegion(file_name, var_name, info, box=region, time_bounds=period, frequency=frequency,
                                    units=units)
    if frequency == "daily":
        Nyears = len(tab) / 365
    elif frequency == "yearly":
        Nyears = len(tab) / 1
    else:
        Nyears = len(tab) / 12
    # processing data
    for num in preprocessing.keys()[1:]:
        process = preprocessing[num].keys()[0]
        print str(num)+' '+str(process)
        if preprocessing[num][process]:
            if process not in dict_processes.keys():
                EnsoErrorsWarnings.UnknownPreprocessing(process, dict_processes.keys(), INSPECTstack())
            tab, info = dict_processes[process](tab, info, region=region, period=period, frequency=frequency,
                                                **preprocessing[num][process])
    return tab, Nyears, info
