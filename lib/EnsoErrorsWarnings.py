# -*- coding:UTF-8 -*-
from sys import exit as sys_exit


#---------------------------------------------------#
# colors for printing
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of defined errors and warning functions used in EnsoUvcdatToolsLib.py
#
def PlusCommaSpace(string):
    """
    #################################################################################
    Description:
    Adds a comma and a space if the string is not empty or if the string is not composed of space only
    #################################################################################

    :param string: string
        string to which comma_space could be added
    :return string: string
        given string+', ' if applicable

    Examples
    ----------
    >>> string = string('')
    >>> print string
    ''
    # or
    >>> string = string('     ')
    >>> print string
    '     '
    # or
    >>> string = string('Where there’s a will')
    >>> print string
    'Where there’s a will, '
    """
    if string.isspace():
        return string
    elif not string:
        return string
    else:
        return string+', '


def MessageFormating(inspect_stack):
    """
    #################################################################################
    Description:
    Formats inspect.stack() as '   File filename, line n, in module'
    #################################################################################

    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return string: string
        formatted inspect.stack() in a string (using PlusCommaSpace)

    Examples
    ----------
    >>> string = message_formating(inspect_stack)
    >>> print string
    '   File filename, line n, in module'
    """
    string = '   '
    # adds file's name
    if inspect_stack[0][1] != '<stdin>':
        string = PlusCommaSpace(string)+'File '+str(inspect_stack[0][1])
    # adds line number
    string = PlusCommaSpace(string)+'line '+str(inspect_stack[0][2])
    # adds module's name
    if inspect_stack[0][3] != '<module>':
        string = PlusCommaSpace(string)+'in '+str(inspect_stack[0][3])
    return string


def MyWarning(list_strings):
    """
    #################################################################################
    Description:
    Prints the strings in 'list_strings' and continues
    #################################################################################

    :param list_strings: list
        list of strings to print
    :return:
    """
    for ii in range(2): print bcolors.WARNING + ""
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        print str().ljust(5) + str(string)
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for ii in range(2): print '' + bcolors.ENDC
    return


# ---------------------------------------------------------------------------------------------------------------------#
#
# ERRORS
#
def MyError(list_strings):
    """
    #################################################################################
    Description:
    Prints the strings in 'list_strings' and exits
    #################################################################################

    :param list_strings: list
        list of strings to print
    :return:
    """
    for ii in range(2): print bcolors.FAIL + ""
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        print str().ljust(5) + str(string)
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for ii in range(2): print '' + bcolors.ENDC
    sys_exit(1)
    return


def MismatchShapesError(tab1, tab2, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of array shape error
    Prints strings and exits
    #################################################################################

    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    try: name1 = tab1.name
    except: name1 = 'no_name'
    try: name2 = tab2.name
    except: name2 = 'no_name'
    list_strings = ["ERROR " + MessageFormating(inspect_stack) + ": array shape",
                    str().ljust(5) + "arrays shapes mismatch: " + str(name1) + " = " + str(tab1.shape) + "', and "
                    + str(name2) + " = " + str(tab2.shape)]
    MyError(list_strings)
    return


def ObjectTypeError(parameter_name, type_parameter, type_parameter_should_be, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of object type error
    Prints strings and exits
    #################################################################################

    :param parameter_name: string
        name of a parameter from which the error comes from
    :param type_parameter: string
        parameter's type
    :param type_parameter_should_be: string
        what the parameter's type should be
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR " + MessageFormating(inspect_stack) + ": object type",
                    str().ljust(5) + str(parameter_name) + ": should be '" + str(type_parameter_should_be) + "', not '"
                    + str(type_parameter) + "'"]
    MyError(list_strings)
    return


def TooShortTimePeriod(metric_name, length, minimum_length, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of a too short time-period
    Prints strings and exits
    #################################################################################

    :param metric_name: string
        name of the metric from which the error comes from
    :param length: integer
        length of the time axis of the variable
    :param minimum_length: integer
        minimum length of the time axis for the metric to make sens (defined in the metrics collection)
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR " + MessageFormating(inspect_stack) + ": too short time-period",
                    str().ljust(5) + str(metric_name) + ": the time-period is too short: " + str(length)
                    + " (minimum time-period: " + str(minimum_length) + ")"]
    MyError(list_strings)
    return


def UnlikelyUnits(var_name, name_in_file, units, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of unlikely units
    Prints strings and exits
    #################################################################################

    :param var_name: string
        generic name of the variable that has unlikely units
    :param name_in_file: string
        name of the variable in the file (usually the short_name) that has unlikely units
    :param units: string
        units of the variable
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR " + MessageFormating(inspect_stack) + ": units",
                    str().ljust(5) + "the file says that " + str(var_name) + " (" + str(name_in_file)
                    + ") is in " + str(units) + " but it seems unlikely"]
    MyError(list_strings)
    return


def UnknownFrequency(frequency, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of unknown frequency
    Prints strings
    #################################################################################

    :param frequency: string
        frequency of a dataset (should by daily, monthly or yearly)
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + MessageFormating(inspect_stack) + ": frequency",
                    str().ljust(5) + "unknown frequency: " + str(frequency)]
    MyError(list_strings)
    return


def UnknownKeyArg(arg, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of unknown argument
    Prints strings
    #################################################################################

    :param arg: string
        argument of a function
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + MessageFormating(inspect_stack) + ": argument",
                    str().ljust(5) + "unknown argument(s): " + str(arg)]
    MyError(list_strings)
    return


def UnknownUnits(var_name, name_in_file, units, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'MyError' in the case of unknown units
    Prints strings and exits
    #################################################################################

    :param var_name: string
        generic name of the variable that has unlikely units
    :param name_in_file: string
        name of the variable in the file (usually the short_name) that has unlikely units
    :param units: string
        units of the variable
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + MessageFormating(inspect_stack) + ": units",
                    str().ljust(5) + "unknown units: " + str(var_name) + " (" + str(name_in_file)
                    + ") is in " + str(units)]
    MyError(list_strings)
    return
# ---------------------------------------------------------------------------------------------------------------------#
