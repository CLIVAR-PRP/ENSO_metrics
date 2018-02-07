# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack
from numpy.ma import array as NUMPYma__array
from numpy.ma import average as NUMPYma__average
from numpy.ma import count as NUMPYma__count
from numpy.ma import isarray as NUMPYma__isarray
from numpy.ma import masked_where as NUMPYma__masked_where
from numpy.ma import ones as NUMPYma__ones
from numpy.ma import sqrt as NUMPYma__sqrt
from numpy.ma import sum as NUMPYma__sum

# ENSO_metrics package functions:
import EnsoErrorsWarnings


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of functions without UVCDAT
#
def StringInDict(string_or_list, dictionary, inspect_stack):
    """
    #################################################################################
    Description:
    Tests if 'string_or_list' is in the given 'dictionary'
    #################################################################################

    :param string_or_list: string or list
        key or list of keys to look for in 'dictionary'
    :param dictionary: dict
        dictionary in which the given keys are looked for
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    # test input parameters
    if not isinstance(dictionary, dict):
        EnsoErrorsWarnings.ObjectTypeError('dictionary', 'dictionary', type(dictionary), INSPECTstack())

    if isinstance(string_or_list, basestring):
        # 'string_or_list' is a string
        if string_or_list not in dictionary.keys():
            # key 'string_or_list' is not in 'dictionary' -> raise error
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(inspect_stack) + ": item not included",
                str().ljust(5) + " key: " + str(string_or_list) + " is not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    elif isinstance(string_or_list, list):
        # 'string_or_list' is a list
        key_not_included = list()
        for key in string_or_list:
            if key not in dictionary.keys():
                # lists keys that are not in 'dictionary'
                key_not_included.append(key)
        if key_not_included:
            # if 'key_not_included' is not empty -> raise error
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.MessageFormating(inspect_stack) + ": item not included",
                str().ljust(5) + " key(s): " + str(key_not_included) + " are (is) not in the given dictionary",
                str().ljust(10) + "key(s) in the dictionary: " + str(sorted(dictionary.keys()))
            ]
            EnsoErrorsWarnings.MyError(list_strings)
    else:
        # 'string_or_list' is neither a string nor a list -> raise error
        EnsoErrorsWarnings.ObjectTypeError('string_or_list', '[string, list]', type(string_or_list), INSPECTstack())
    return


# ---------------------------------------------------------------------------------------------------------------------#
#
# statistics based on numpy
#
def __covariance(x, y, weights=None, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes the covariance
    Based on numpy.ma
    #################################################################################

    :param x: array_like
        1D array containing numbers whose covariance with 'y' is desired
        If a is not an array, a conversion is attempted
    :param y: array_like
        1D array containing numbers whose covariance with 'x' is desired
        If a is not an array, a conversion is attempted
    :param weights: array_like, optional
        the importance that each element has in the computation of the average. The weights array can either be 1-D (in
        which case its length must be the size of a along the given axis) or of the same shape as a. If weights=None,
        then all data in a are assumed to have a weight equal to one. If weights is complex, the imaginary parts are
        ignored.
        Default value is None
    :param centered: integer
        integer flag for whether the calculation should be centered (mean removed from array).
        A value of 1 indicates calculation should be centered, while any other value indicates that it should not.
        Default value is 1 (centered)
    :param biased: integer
        flag indicating whether calculation should be biased.
        A value of 1 indicates calculation should be biased, while any other value indicates that it should not.
        Default value is 1 (biased)

    :return: float
        covariance of x and y
    """
    # test if there is a contradiction between 'weights' and 'biased'
    if weights is not None and biased != 1:
        list_strings = [
            "ERROR " + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": contradiction between arguments",
            str().ljust(5) + "Error in covariance, you cannot have weights and unbiased together"
        ]
        EnsoErrorsWarnings.MyError(list_strings)
    # test if x and y are masked_array
    if not NUMPYma__isarray(x):
        x = NUMPYma__array(x)
    if not NUMPYma__isarray(y):
        y = NUMPYma__array(y)
    # remove mean from x and y if the calculation should be centered
    if centered == 1:
        xmean = NUMPYma__average(x, weights=weights, axis=0)
        ymean = NUMPYma__average(y, weights=weights, axis=0)
        x = x - xmean
        y = y - ymean
        del(xmean)
        del(ymean)
    # if 'weights' is None create an array of ones
    if weights is None:
        weights = NUMPYma__ones(x.shape, dtype=x.dtype.char)
    else:
        if not NUMPYma__isarray(weights):
            weights = NUMPYma__array(weights)
    # if x has a mask, apply the mask on weights
    if not (x.mask is None):
        weights = NUMPYma__masked_where(x.mask, weights)
    # if y has a mask, apply the mask on weights
    if not (y.mask is None):
            weights = NUMPYma__masked_where(y.mask, weights)
    # compute the covariance
    if biased == 1:
        cov = NUMPYma__sum(x * y * weights, axis=0) / NUMPYma__sum(weights, axis=0)
    else:
        cov = NUMPYma__sum(x * y, axis=0) / (NUMPYma__count(x * y, axis=0) - 1)
    return cov


def __rms(x, y, weights=None, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Computes the root mean square
    Based on numpy.ma
    #################################################################################

    :param x: array_like
        1D array containing numbers whose root mean square with 'y' is desired
        If a is not an array, a conversion is attempted
    :param y: array_like
        1D array containing numbers whose root mean square with 'x' is desired
        If a is not an array, a conversion is attempted
    :param weights: array_like, optional
        the importance that each element has in the computation of the average. The weights array can either be 1-D (in
        which case its length must be the size of a along the given axis) or of the same shape as a. If weights=None,
        then all data in a are assumed to have a weight equal to one. If weights is complex, the imaginary parts are
        ignored.
        Default value is None
    :param centered: integer
        integer flag for whether the calculation should be centered (mean removed from array).
        A value of 1 indicates calculation should be centered, while any other value indicates that it should not.
        Default value is 1 (centered)
    :param biased: integer
        flag indicating whether calculation should be biased.
        A value of 1 indicates calculation should be biased, while any other value indicates that it should not.
        Default value is 1 (biased)

    :return: float
        root mean square of x and y

    """
    return __std(x - y, centered=centered, biased=biased, weights=weights)


def __std(x, weights=None, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes the standard deviation
    Based on numpy.ma
    #################################################################################

    :param x: array_like
        1D array containing numbers whose standard deviation is desired
        If a is not an array, a conversion is attempted
    :param weights: array_like, optional
        the importance that each element has in the computation of the average. The weights array can either be 1-D (in
        which case its length must be the size of a along the given axis) or of the same shape as a. If weights=None,
        then all data in a are assumed to have a weight equal to one. If weights is complex, the imaginary parts are
        ignored.
        Default value is None
    :param centered: integer
        integer flag for whether the calculation should be centered (mean removed from array).
        A value of 1 indicates calculation should be centered, while any other value indicates that it should not.
        Default value is 1 (centered)
    :param biased: integer
        flag indicating whether calculation should be biased.
        A value of 1 indicates calculation should be biased, while any other value indicates that it should not.
        Default value is 1 (biased)

    :return: float
        standard deviation of x
    """
    return NUMPYma__sqrt(__variance(x, weights=weights, centered=centered, biased=biased))


def __variance(x, weights=None, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes the variance
    Based on numpy.ma
    #################################################################################

    :param x: array_like
        1D array containing numbers whose variance is desired
        If a is not an array, a conversion is attempted
    :param weights: array_like, optional
        the importance that each element has in the computation of the average. The weights array can either be 1-D (in
        which case its length must be the size of a along the given axis) or of the same shape as a. If weights=None,
        then all data in a are assumed to have a weight equal to one. If weights is complex, the imaginary parts are
        ignored.
        Default value is None
    :param centered: integer
        integer flag for whether the calculation should be centered (mean removed from array).
        A value of 1 indicates calculation should be centered, while any other value indicates that it should not.
        Default value is 1 (centered)
    :param biased: integer
        flag indicating whether calculation should be biased.
        A value of 1 indicates calculation should be biased, while any other value indicates that it should not.
        Default value is 1 (biased)

    :return: float
        variance of x
    """
    return __covariance(x, x, weights=weights, centered=centered, biased=biased)

