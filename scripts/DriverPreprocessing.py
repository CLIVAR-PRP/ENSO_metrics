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
    print '\033[93m' + str().ljust(5) + "DriverPreprocessing preprocess" + '\033[0m'
    print '\033[93m' + str().ljust(10) + "read " + str(var_name) + " from " + str(file_name) + '\033[0m'
    info = 'Preprocessing:'
    # read 'var_name' in 'file_name' in the given 'region' and 'period'
    tab, info = ReadAndSelectRegion(file_name, var_name, info, box=region, time_bounds=period, frequency=frequency,
                                    units=units)
    print '\033[93m' + str().ljust(10) + "after ReadAndSelectRegion" + '\033[0m'
    print '\033[93m' + str().ljust(15) + "tab.shape = " + str(tab.shape) + '\033[0m'
    print '\033[93m' + str().ljust(15) + "tab.axes = " + str([ax.id for ax in tab.getAxisList()]) + '\033[0m'
    print '\033[93m' + str().ljust(15) + "tab.units = " + str(tab.units) + '\033[0m'
    if frequency == "daily":
        Nyears = len(tab) / 365
    elif frequency == "yearly":
        Nyears = len(tab) / 1
    else:
        Nyears = len(tab) / 12
    # processing data
    for num in preprocessing.keys()[1:]:
        process = preprocessing[num].keys()[0]
        print '\033[93m' + str().ljust(10) + str(num)+' '+str(process) + '\033[0m'
        if preprocessing[num][process]:
            if process not in dict_processes.keys():
                EnsoErrorsWarnings.UnknownPreprocessing(process, dict_processes.keys(), INSPECTstack())
            tab, info = dict_processes[process](tab, info, region=region, period=period, frequency=frequency,
                                                **preprocessing[num][process])
            print '\033[93m' + str().ljust(15) + "tab.shape = " + str(tab.shape) + '\033[0m'
            print '\033[93m' + str().ljust(15) + "tab.axes = " + str([ax.id for ax in tab.getAxisList()]) + '\033[0m'
    return tab, Nyears, info
