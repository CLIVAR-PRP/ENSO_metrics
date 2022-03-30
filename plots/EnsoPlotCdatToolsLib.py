# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Provide functions based on CDAT (https://cdat.llnl.gov/) used for plots
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
from copy import deepcopy
from inspect import stack as INSPECTstack
from numpy import array as NUMPYarray
from numpy import bool_ as NUMPYbool_
from numpy import mean as NUMPYmean
from numpy.random import randint as NUMPYrandom__randint
# CDAT
from cdms2 import createAxis as CDMS2createAxis
from cdms2 import createVariable as CDMS2createVariable
from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
from genutil.statistics import percentiles as GENUTILstatistics__percentiles
from genutil.statistics import rms as GENUTILstatistics__rms
from MV2 import array as MV2array
from MV2 import average as MV2average
from MV2 import concatenate as MV2concatenate
from MV2 import masked_where as MV2masked_where
from MV2 import max as MV2max
from MV2 import min as MV2min
from MV2 import where as MV2where
from MV2 import zeros as MV2zeros
# ENSO_metrics package functions:
from EnsoMetrics import EnsoErrorsWarnings


# ---------------------------------------------------#


# ---------------------------------------------------#
# Parameters
# ---------------------------------------------------#
# Colors for prints
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def member_average(ldict_members, lvariable):
    """
    #################################################################################
    Description:
    Compute the average across ensemble members
    #################################################################################

    :param ldict_members: dictionary
        dictionary of ensemble members
        level 1: member names
        level 2: variable names
        level 3: 'array' or 'attributes'
        level 4: masked_array arr[region, x, y, z]
    :param lvariable: string
        name of the variable name to read (level 2 of the dictionary)

    :return: masked_array
        ensemble members average
    """
    lmem = sorted(list(ldict_members.keys()), key=lambda v: v.upper())
    list1 = MV2average(MV2array([ldict_members[mem][lvariable]["array"] for mem in lmem]), axis=0)
    ltab = ldict_members[lmem[0]][lvariable]["array"]
    list1 = CDMS2createVariable(list1, axes=ltab.getAxisList(), grid=ltab.getGrid(), mask=ltab.mask, id=ltab.id)
    return list1


def member_range(ldict_members, lvariable, lsample_size, lsig_level=90, lnum_samples=1000000, lmin_samples=100,
                 lmax_samples=10000, larray_mask=None):
    """
    #################################################################################
    Description:
    Generate 'lnum_samples' random composites, with the sample size corresponding to 'lsample_size' to find the
    'lsig_level'% interval. The minimum and maximum values of this interval is returned (i.e., the range)
    #################################################################################

    :param ldict_members: dictionary
        dictionary of ensemble members
        level 1: member names
        level 2: variable names
        level 3: 'array' or 'attributes'
        level 4: masked_array arr[events, region, x, y, z]
    :param lvariable: string
        name of the variable name to read (level 2 of the dictionary)
    :param lsample_size: integer
        number of events to select to create each composite
    :param lsig_level: integer, optional
        confidence level of the significance test (90 is 90% confidence level, the 5% and 95% percentiles are detected)
        default value is 90
    :param lnum_samples: integer, optional
        number of random selections (with replacement) of the given sample to use for the significance test
        default value is 100,000
    :param lmin_samples: integer, optional
        minimum number of composites created when a member is selected
        default value is 100
    :param lmax_samples: integer, optional
        maximum number of composites created when a member is selected
        default value is 10,000
    :param larray_mask: masked_array
        masked_array in which 0 means the values will be masked
        Default is None (data will not be masked)

    :return larray_o: masked_array
        masked_array of the shape arr[range, region, x, y, z], the first dimension is the minimum and maximum values of
        given confidence interval of the Monte Carlo resampling
    """
    if lmin_samples > lnum_samples / 100.:
        lmin_samples = int(round(lnum_samples / 100., 0))
    if lmax_samples > lnum_samples / 10.:
        lmax_samples = int(round(lnum_samples / 10., 0))
    s1 = (100 - lsig_level) / 2.
    s2 = lsig_level + s1
    # list the member names
    lmem = sorted(list(ldict_members.keys()), key=lambda v: v.upper())
    lnbr_mem = len(lmem)
    # Monte Carlo resampling
    larray_o = 1
    lnbr_samples_done = 0
    while lnbr_samples_done < lnum_samples:
        # first, randomly select a member
        smem = lmem[NUMPYrandom__randint(0, lnbr_mem)]
        # second, randomly select a number of composites that will be created
        lnbr_left_to_do = lnum_samples - lnbr_samples_done
        if lnbr_left_to_do < lmax_samples:
            if lnbr_left_to_do < min(lmin_samples * 10., lmax_samples / 10.) or lnbr_left_to_do <= lmin_samples:
                snum = deepcopy(lnbr_left_to_do)
            else:
                snum = NUMPYrandom__randint(lmin_samples, lnbr_left_to_do)
        else:
            snum = NUMPYrandom__randint(lmin_samples, lmax_samples)
        # random selection among events
        larray = ldict_members[smem][lvariable]["array"]
        idx = NUMPYrandom__randint(0, len(larray), (snum, lsample_size))
        samples = NUMPYarray(larray)[idx]
        # average randomly selected years (create lnum composites)
        stat = NUMPYmean(samples, axis=1)
        # concatenate array
        if lnbr_samples_done == 0:
            larray_o = deepcopy(stat)
        else:
            larray_o = MV2concatenate((larray_o, stat))
        # increment the number of composite created
        lnbr_samples_done += snum
        # delete
        del idx, larray, lnbr_left_to_do, samples, smem, snum, stat
    # compute lower and higher threshold of the given significance level
    if lsig_level == 100:
        low = MV2min(larray_o, axis=0)
        high = MV2max(larray_o, axis=0)
    else:
        low = GENUTILstatistics__percentiles(larray_o, percentiles=[s1], axis=0)[0]
        high = GENUTILstatistics__percentiles(larray_o, percentiles=[s2], axis=0)[0]
    # create array
    larray_o = MV2array([low, high])
    ax_o = CDMS2createAxis(list(range(2)), id="lowhigh")
    ax_o.short_name = str(s1) + "_and_" + str(s2) + "_precentiles"
    ax_o.long_name = str(s1) + " and " + str(s2) + " Precentiles of the Bootstrap"
    # create variable
    larray = ldict_members[lmem[0]][lvariable]["array"]
    axes = [ax_o] + larray.getAxisList()[1:]
    if type(larray.mask) == NUMPYbool_:
        larray_o = CDMS2createVariable(larray_o, axes=axes, grid=larray.getGrid(), id="range")
    else:
        lmask = MV2zeros(larray_o.shape)
        lmask[:] = larray[0].mask
        larray_o = CDMS2createVariable(larray_o, axes=axes, grid=larray.getGrid(), mask=lmask, id="range")
        del lmask
    # mask with given mask if
    if larray_mask is not None:
        lmask = MV2zeros(larray_o.shape)
        lmask[:] = larray_mask
        larray_o = my_mask(larray_o, lmask)
        del lmask
    return larray_o


def minimaxi(array, mini=1e20, maxi=-1e20):
    """
    #################################################################################
    Description:
    Find the minimum and the maximum
    #################################################################################

    :param array: masked_array
        masked_array for which the minimum and maximum values are wanted
    :param mini: float, optional
        preset minimum value
        default is 1e20
    :param maxi: float, optional
        preset maximum value
        default is -1e20

    :return: list
        minimum and maximum values
    """
    if isinstance(array, list) is True:
        for tmp in array:
            v1, v2 = float(MV2min(tmp)), float(MV2max(tmp))
            if v1 < mini:
                mini = deepcopy(v1)
            if v2 > maxi:
                maxi = deepcopy(v2)
    else:
        v1, v2 = float(MV2min(array)), float(MV2max(array))
        if v1 < mini:
            mini = deepcopy(v1)
        if v2 > maxi:
            maxi = deepcopy(v2)
    return mini, maxi


def my_mask(larray, larray_mask):
    """
    #################################################################################
    Description:
    Mask where given mask equal 1
    #################################################################################

    :param larray: masked_array
        masked_array to mask
    :param larray_mask: masked_array
        masked_array in which 0 means the values will be masked

    :return larray_out: masked_array
        masked_array with larray where larray_mask != 0 else masked
    """
    larray_out = MV2masked_where(larray_mask == 0, larray)
    return larray_out


def my_rmse(larray1, larray2=None):
    """
    #################################################################################
    Description:
    Mask where given mask equal 1
    #################################################################################

    :param larray1: masked_array
        masked_array to evaluate
    :param larray2: masked_array, optional
        masked_array to evaluate larray1 against
        default is None (larray1 evaluated against 0)

    :return rmse: float
        rmse between larray1 and 0 or larray2 if given
    """
    larray_i1 = MV2array(larray1)
    if larray2 is None:
        larray_i2 = MV2zeros(larray_i1.shape)
    else:
        larray_i2 = MV2array(larray2)
    rmse = float(GENUTILstatistics__rms(larray_i1, larray_i2, centered=0, biased=1))
    return rmse


def my_sig(larray1, larray2, larray3):
    tab = my_mask(NUMPYmean(larray2, axis=0), larray3).compressed()
    low, hig = larray1
    vlow = MV2zeros(low.shape)
    vlow = MV2where(tab < low, 1, vlow)
    vhigh = MV2zeros(hig.shape)
    vhigh = MV2where(tab > hig, 1, vhigh)
    significance = MV2zeros(low.shape)
    significance = MV2where(vlow + vhigh > 0, 1, significance).compressed()
    return significance


def significance_range(ldict_members, lvariable, lsample_size, lsig_level=90, lnum_samples=1000000, larray_mask=None):
    """
    #################################################################################
    Description:
    Generate 'lnum_samples' random ensemble mean, with the sample size corresponding to 'lsample_size' to find the
    'lsig_level'% interval. The minimum and maximum values of this interval is returned (i.e., the range)
    #################################################################################

    :param ldict_members:
        dictionary of ensemble members
        level 1: member names
        level 2: variable names
        level 3: 'array' or 'attributes'
        level 4: masked_array arr[region, x, y, z]
    :param lvariable: string
        name of the variable name to read (level 2 of the dictionary)
    :param lsample_size: integer
        number of events to select to create each composite
    :param lsig_level: integer, optional
        confidence level of the significance test (90 is 90% confidence level, the 5% and 95% percentiles are detected)
        default value is 90
    :param lnum_samples: integer, optional
        number of random selections (with replacement) of the given sample to use for the significance test
        default value is 100,000
    :param larray_mask: masked_array
        masked_array in which 0 means the values will be masked
        Default is None (data will not be masked)

    :return larray_o: masked_array
        masked_array of the shape arr[range, region, x, y, z], the first dimension is the minimum and maximum values of
        given confidence interval of the Monte Carlo resampling
    """
    # percentiles of the distribution to compute
    s1 = (100 - lsig_level) / 2.
    s2 = lsig_level + s1
    # create array from dict
    lmem = sorted(list(ldict_members.keys()), key=lambda v: v.upper())
    larray1 = ldict_members[lmem[0]][lvariable]["array"]
    larray2 = MV2array([ldict_members[mem][lvariable]["array"] for mem in lmem])
    # random selection among all years
    idx = NUMPYrandom__randint(0, len(larray2), (lnum_samples, lsample_size))
    samples = NUMPYarray(larray2)[idx]
    # average randomly selected years (create lnum_samples composites)
    stat2 = NUMPYmean(samples, axis=1)
    # compute lower and higher threshold of the given significance level
    if lsig_level == 100:
        low = MV2min(stat2, axis=0)
        high = MV2max(stat2, axis=0)
    else:
        low = GENUTILstatistics__percentiles(stat2, percentiles=[s1], axis=0)[0]
        high = GENUTILstatistics__percentiles(stat2, percentiles=[s2], axis=0)[0]
    # create array
    larray_o = MV2array([low, high])
    ax_o = CDMS2createAxis(list(range(2)), id="lowhigh")
    ax_o.short_name = str(s1) + "_and_" + str(s2) + "_precentiles"
    ax_o.long_name = str(s1) + " and " + str(s2) + " Precentiles of the Bootstrap"
    # create variable
    axes = [ax_o] + larray1.getAxisList()
    if type(larray1.mask) == NUMPYbool_:
        larray_o = CDMS2createVariable(larray_o, axes=axes, grid=larray1.getGrid(), id="range")
    else:
        lmask = MV2zeros(larray_o.shape)
        lmask[:] = larray1.mask
        larray_o = CDMS2createVariable(larray_o, axes=axes, grid=larray1.getGrid(), mask=lmask, id="range")
        del lmask
    # mask with given mask if
    if larray_mask is not None:
        lmask = MV2zeros(larray_o.shape)
        lmask[:] = larray_mask
        larray_o = my_mask(larray_o, lmask)
        del lmask
    return larray_o


def read_data(netcdf_file, netcdf_var, mod_name, ref_name, dict_diagnostic, dict_metric, member=None,
              netcdf_var_extra=None):
    """
    #################################################################################
    Description:
    Read and organize the given variables and associated metadata from the given netCDF file
    Select diagnostic and metric values
    #################################################################################

    :param netcdf_file: string
        path and file name of the netcdf file to read
        e.g., '/path/to/file/file_name.nc'
    :param netcdf_var: list
        list of the variable patterns to read
        e.g., ['var1__', 'var2__']
    :param mod_name: string
        name of the model
        e.g., 'ACCESS1-0'
    :param ref_name: string
        name of the reference observational dataset
        e.g., 'ERA-Interim'
    :param dict_diagnostic: dictionary
        dictionary containing the diagnostic values of the model and the observational datasets
        e.g., {'ACCESS1-0': 1., 'ERA-Interim': 1.}
    :param dict_metric: dictionary
        dictionary containing the metric values comparing the model to the observational datasets
        e.g., {'ERA-Interim': 1.}
    :param member: string, optional
        name of the model's member
        e.g., 'r1i1p1'
        default is None
    :param netcdf_var_extra: list
        list of the extra variable patterns to read if available
        e.g., ['extra_var1__', 'extra_var2__']

    :return:
    """
    # open netcdf file
    CDMS2setAutoBounds("on")
    ff = CDMS2open(netcdf_file)
    # global attributes
    att_glo = read_global_attributes(ff)
    # list the variables in the file
    list_variables = sorted(ff.listvariables(), key=lambda v: v.upper())
    # read netcdf
    dict_mod, dict_mod_extra, dict_ref, dict_ref_extra = dict(), dict(), dict(), dict()
    for var in netcdf_var:
        var_name1 = str(var) + str(mod_name)
        var_name2 = str(var_name1) + "_" + str(member) if isinstance(member, str) is True else None
        if var_name1 not in list_variables and (var_name2 is None or var_name2 not in list_variables):
            arr, att = None, None
            my_name = deepcopy(var_name1)
            if isinstance(member, str) is True:
                my_name += " or " + str(var_name2)
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": cannot find given variable",
                str().ljust(5) + "cannot find " + str(my_name) + " in " + str(netcdf_file),
                str().ljust(5) + "variables in file: " + str(list_variables)]
            EnsoErrorsWarnings.my_error(list_strings)
        elif var_name1 in list_variables:
            arr, att = read_variable_and_attributes(ff, var_name1)
        else:
            arr, att = read_variable_and_attributes(ff, var_name2)
        dict_mod[var] = {"array": arr, "attributes": att}
        del arr, att, var_name1, var_name2
        var_name1 = str(var) + str(ref_name)
        var_name2 = str(var_name1) + "_" + str(ref_name)
        if var_name1 not in list_variables and var_name2 not in list_variables:
            arr, att = None, None
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": cannot find given variable",
                str().ljust(5) + "cannot find " + str(var_name1) + " or " + str(var_name2) + " in " + str(netcdf_file),
                str().ljust(5) + "variables in file: " + str(list_variables)]
            EnsoErrorsWarnings.my_error(list_strings)
        elif var_name1 in list_variables:
            arr, att = read_variable_and_attributes(ff, var_name1)
        else:
            arr, att = read_variable_and_attributes(ff, var_name2)
        dict_ref[var] = {"array": arr, "attributes": att}
        del arr, att, var_name1, var_name2
    if isinstance(netcdf_var_extra, list) is True:
        for var in netcdf_var_extra:
            var_name1 = str(var) + str(mod_name)
            var_name2 = str(var_name1) + "_" + str(member) if isinstance(member, str) is True else None
            if var_name1 not in list_variables and (var_name2 is None or var_name2 not in list_variables):
                arr, att = None, None
            elif var_name1 in list_variables:
                arr, att = read_variable_and_attributes(ff, var_name1)
            else:
                arr, att = read_variable_and_attributes(ff, var_name2)
            if arr is not None:
                dict_mod_extra[var] = {"array": arr, "attributes": att}
            del arr, att, var_name1, var_name2
            var_name1 = str(var) + str(ref_name)
            var_name2 = str(var_name1) + "_" + str(ref_name)
            if var_name1 not in list_variables and var_name2 not in list_variables:
                arr, att = None, None
            elif var_name1 in list_variables:
                arr, att = read_variable_and_attributes(ff, var_name1)
            else:
                arr, att = read_variable_and_attributes(ff, var_name2)
            if arr is not None:
                dict_ref_extra[var] = {"array": arr, "attributes": att}
            del arr, att, var_name1, var_name2
    # read diagnostic value
    list_keys = sorted(list(dict_diagnostic.keys()), key=lambda v: v.upper())
    var_name1 = deepcopy(mod_name)
    var_name2 = str(var_name1) + "_" + str(member) if isinstance(member, str) is True else None
    if var_name1 not in list_keys and (var_name2 is None or var_name2 not in list_keys):
        dia_mod = None
        my_name = deepcopy(var_name1)
        if isinstance(member, str) is True:
            my_name += " or " + str(var_name2)
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": cannot find diagnostic for given model",
            str().ljust(5) + "cannot find " + str(my_name) + " in given dictionary",
            str().ljust(5) + "dictionary keys: " + str(list_keys)]
        EnsoErrorsWarnings.my_error(list_strings)
    elif var_name1 in list_keys:
        dia_mod = dict_diagnostic[var_name1]
    else:
        dia_mod = dict_diagnostic[var_name2]
    var_name3 = deepcopy(ref_name)
    var_name4 = str(var_name3) + "_" + str(ref_name)
    if var_name3 not in list_keys and var_name4 not in list_keys:
        dia_ref = None
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
            ": cannot find diagnostic for given reference",
            str().ljust(5) + "cannot find " + str(var_name3) + " or " + str(var_name4) + " in given dictionary",
            str().ljust(5) + "dictionary keys: " + str(list_keys)]
        EnsoErrorsWarnings.my_error(list_strings)
    elif var_name3 in list_keys:
        dia_ref = dict_diagnostic[var_name3]
    else:
        dia_ref = dict_diagnostic[var_name4]
    # read metric value
    list_keys = sorted(list(dict_metric.keys()), key=lambda v: v.upper())
    if var_name1 not in list_keys and var_name2 not in list_keys and var_name3 not in list_keys and \
            var_name4 not in list_keys:
        met_val = None
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": cannot find metric for given reference",
            str().ljust(5) + "cannot find " + str(var_name1) + " or " + str(var_name2) + " or " + str(var_name3) +
            " or " + str(var_name4) + " in given dictionary",
            str().ljust(5) + "dictionary keys: " + str(list_keys)]
        EnsoErrorsWarnings.my_error(list_strings)
    elif var_name1 in list_keys:
        met_val = dict_metric[var_name1]
    elif var_name2 in list_keys:
        met_val = dict_metric[var_name2]
    elif var_name3 in list_keys:
        met_val = dict_metric[var_name3]
    else:
        met_val = dict_metric[var_name4]
    return dict_mod, dict_ref, dia_mod, dia_ref, met_val, att_glo, dict_mod_extra, dict_ref_extra


def read_global_attributes(opened_file):
    attributes = dict()
    for att in opened_file.listglobal():
        tmp = opened_file.attributes[att]
        if isinstance(tmp, list) is True and len(tmp) == 1:
            tmp = tmp[0]
        attributes[att] = tmp
        del tmp
    return attributes


def read_variable_and_attributes(opened_file, variable_name):
    array = opened_file(variable_name)
    attributes = dict()
    for att in opened_file.listattribute(variable_name):
        tmp = opened_file.getattribute(variable_name, att)
        if isinstance(tmp, list) is True and len(tmp) == 1:
            tmp = tmp[0]
        attributes[att] = tmp
        del tmp
    return array, attributes
# ---------------------------------------------------------------------------------------------------------------------#
