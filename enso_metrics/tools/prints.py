# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Classes and functions for errors and warnings
# ---------------------------------------------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------------------------------------------------#
# Classes
# ---------------------------------------------------------------------------------------------------------------------#
class BackgroundColors:
    blue = '\033[94m'
    green = '\033[92m'
    orange = '\033[93m'
    red = '\033[91m'
    normal = '\033[0m'
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
def plus_comma_space(string):
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
    string = string('')
    print string
    ''
    # or
    string = string('     ')
    print string
    '     '
    # or
    string = string('Where there’s a will')
    print string
    'Where there’s a will, '
    """
    if string.isspace():
        return string
    elif not string:
        return string
    else:
        return str(string) + ", "


def message_formating(inspect_stack):
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
    string = message_formating(inspect_stack)
    print string
    '   File filename, line n, in module'
    """
    string = "   "
    # adds file's name
    if inspect_stack[0][1] != "<stdin>":
        string = plus_comma_space(string) + "File " + str(inspect_stack[0][1])
    # adds line number
    string = plus_comma_space(string) + "line " + str(inspect_stack[0][2])
    # adds module's name
    if inspect_stack[0][3] != "<module>":
        string = plus_comma_space(string) + "in " + str(inspect_stack[0][3])
    return string


def my_warning(list_strings):
    """
    #################################################################################
    Description:
    Prints the strings in 'list_strings' and continues
    #################################################################################

    :param list_strings: list
        list of strings to print
    :return:
    """
    tmp = "\n\n" + str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        tmp += "\n" + str().ljust(5) + str(string)
    tmp += "\n" + str().ljust(5) + "%%%%%     -----     %%%%%\n\n"
    print(BackgroundColors.orange + tmp + BackgroundColors.normal)


def my_error(list_strings: list[str]):
    """
    #################################################################################
    Description:
    Prints the strings in 'list_strings' and exits
    #################################################################################

    :param list_strings: list
        list of strings to print
    :return:
    """
    tmp = "\n\n" + str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        tmp += "\n" + str().ljust(5) + str(string)
    tmp += "\n" + str().ljust(5) + "%%%%%     -----     %%%%%\n\n"
    ValueError(BackgroundColors.red + tmp + BackgroundColors.normal)


def mismatch_shapes_error(tab1, tab2, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_error' in the case of array shape error
    Prints strings and exits
    #################################################################################

    :param tab1: masked_array
    :param tab2: masked_array
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    try:
        name1 = tab1.name
    except:
        name1 = "no_name"
    try:
        name2 = tab2.name
    except:
        name2 = "no_name"
    list_strings = ["ERROR " + message_formating(inspect_stack) + ": array shape",
                    str().ljust(5) + "arrays shapes mismatch: " + str(name1) + " = " + str(tab1.shape) + "', and " +
                    str(name2) + " = " + str(tab2.shape)]
    my_warning(list_strings)


def object_type_error(parameter_name, type_parameter, type_parameter_should_be, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_error' in the case of object type error
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
    list_strings = ["ERROR " + message_formating(inspect_stack) + ": object type",
                    str().ljust(5) + str(parameter_name) + ": should be '" + str(type_parameter_should_be) +
                    "', not '" + str(type_parameter) + "'"]
    my_warning(list_strings)


def too_short_time_period(metric_name, length, minimum_length, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_warning' in the case of a too short time-period
    Prints strings and exits
    #################################################################################

    :param metric_name: string
        name of the metric from which the error comes from
    :param length: integer
        length of the time axis of the variable
    :param minimum_length: integer
        minimum length of the time axis for the metric to make sens (defined in the recipes collection)
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR " + message_formating(inspect_stack) + ": too short time-period",
                    str().ljust(5) + str(metric_name) + ": the time-period is too short: " + str(length) +
                    " (minimum time-period: " + str(minimum_length) + ")"]
    my_warning(list_strings)


def unlikely_units(var_name, name_in_file, units, min_max, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_warning' in the case of unlikely units
    Prints strings and exits
    #################################################################################

    :param var_name: string
        generic name of the variable that has unlikely units
    :param name_in_file: string
        name of the variable in the file (usually the short_name) that has unlikely units
    :param units: string
        units of the variable
    :param min_max: list
        minimum and maximum values of 'var_name'
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR " + message_formating(inspect_stack) + ": units",
                    str().ljust(5) + "the file says that " + str(var_name) + " (" + str(name_in_file) + ") is in " +
                    str(units) + " but it seems unlikely (" + str(min_max) + ")"]
    my_warning(list_strings)


def unknown_averaging(average, known_average, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_error' in the case of unknown frequency
    Prints strings
    #################################################################################

    :param average: string
        averaging method (axis) (should by horizontal, meridional, temporal or zonal)
    :param known_average: string
        list of defined averaging method (axis)
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + message_formating(inspect_stack) + ": averaging method",
                    str().ljust(5) + "unkwown averaging method (axis): " + str(average),
                    str().ljust(10) + "known averaging method: " + str(sorted(known_average))]
    my_error(list_strings)


def unknown_frequency(frequency, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_error' in the case of unknown frequency
    Prints strings
    #################################################################################

    :param frequency: string
        frequency of a dataset (should by daily, monthly or yearly)
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + message_formating(inspect_stack) + ": frequency",
                    str().ljust(5) + "unknown frequency: " + str(frequency)]
    my_error(list_strings)


def unknown_key_arg(arg, inspect_stack):
    """
    #################################################################################
    Description:
    Function 'my_error' in the case of unknown argument
    Prints strings
    #################################################################################

    :param arg: string
        argument of a function
    :param inspect_stack: array
        list of information about the program/module/line,... created using inspect.stack()
    :return:
    """
    list_strings = ["ERROR" + message_formating(inspect_stack) + ": argument",
                    str().ljust(5) + "unknown argument(s): " + str(arg)]
    my_warning(list_strings)


def unknown_units(var_name, name_in_file, units, inspect_stack):
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
    list_strings = [
        "ERROR" + message_formating(inspect_stack) + ": units",
        str().ljust(5) + "unknown units: " + str(var_name) + " (" + str(name_in_file) + ") is in " + str(units)]
    my_warning(list_strings)


def debug_mode(
        color, title, nbr_spaces, axes1="", axes2="", axes3="", axes4="", axes5="", axes6="", file1="", file2="",
        file3="", file4="", file5="", file6="", line1="", line2="", line3="", line4="", line5="", line6="", nina1="",
        nina2="", nina3="", nina4="", nina5="", nina6="", nino1="", nino2="", nino3="", nino4="", nino5="", nino6="",
        shape1="", shape2="", shape3="", shape4="", shape5="", shape6="", time1="", time2="", time3="", time4="",
        time5="", time6="", var1="", var2="", var3="", var4="", var5="", var6=""):
    """
    #################################################################################
    Description:
    Prints strings to ease debugging
    #################################################################################

    :param color: string
        color code (e.g. '\033[94m' is blue, '\033[92m' is green)
    :param title: string
        name of the section that is printed
    :param nbr_spaces: int
        number of leading spaces before printing the title
    :param axes1: string, optional
        axis list of variable 1
    :param axes2: string, optional
        axis list of variable 2
    :param axes3: string, optional
        axis list of variable 3
    :param axes4: string, optional
        axis list of variable 4
    :param axes5: string, optional
        axis list of variable 5
    :param axes6: string, optional
        axis list of variable 6
    :param file1: string, optional
        file name of variable 1
    :param file2: string, optional
        file name of variable 2
    :param file3: string, optional
        file name of variable 3
    :param file4: string, optional
        file name of variable 4
    :param file5: string, optional
        file name of variable 5
    :param file6: string, optional
        file name of variable 6
    :param line1: string, optional
        just a line to print 1
    :param line2: string, optional
        just a line to print 2
    :param line3: string, optional
        just a line to print 3
    :param line4: string, optional
        just a line to print 4
    :param line5: string, optional
        just a line to print 5
    :param line6: string, optional
        just a line to print 6
    :param nina1: string, optional
        list of nina years 1
    :param nina2: string, optional
        list of nina years 2
    :param nina3: string, optional
        list of nina years 3
    :param nina4: string, optional
        list of nina years 4
    :param nina5: string, optional
        list of nina years 5
    :param nina6: string, optional
        list of nina years 6
    :param nino1: string, optional
        list of nino years 1
    :param nino2: string, optional
        list of nino years 2
    :param nino3: string, optional
        list of nino years 3
    :param nino4: string, optional
        list of nino years 4
    :param nino5: string, optional
        list of nino years 5
    :param nino6: string, optional
        list of nino years 6
    :param shape1: string, optional
        shape of the array containing variable 1
    :param shape2: string, optional
        shape of the array containing variable 2
    :param shape3: string, optional
        shape of the array containing variable 3
    :param shape4: string, optional
        shape of the array containing variable 4
    :param shape5: string, optional
        shape of the array containing variable 5
    :param shape6: string, optional
        shape of the array containing variable 6
    :param time1: string, optional
        time bounds of variable 1
    :param time2: string, optional
        time bounds of variable 2
    :param time3: string, optional
        time bounds of variable 3
    :param time4: string, optional
        time bounds of variable 4
    :param time5: string, optional
        time bounds of variable 5
    :param time6: string, optional
        time bounds of variable 6
    :param var1: string, optional
        variable name 1
    :param var2: string, optional
        variable name 2
    :param var3: string, optional
        variable name 3
    :param var4: string, optional
        variable name 4
    :param var5: string, optional
        variable name 5
    :param var6: string, optional
        variable name 6

    :return:
    """
    # first variable
    print(color + str().ljust(nbr_spaces) + title + BackgroundColors.normal)
    if file1:
        print(color + str().ljust(nbr_spaces + 5) + "file name 1: " + file1 + BackgroundColors.normal)
    if var1:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 1: " + var1 + BackgroundColors.normal)
    if axes1:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 1: " + axes1 + BackgroundColors.normal)
    if time1:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 1: " + time1 + BackgroundColors.normal)
    if shape1:
        print(color + str().ljust(nbr_spaces + 5) + "shape 1: " + shape1 + BackgroundColors.normal)
    if nina1:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 1: " + nina1 + BackgroundColors.normal)
    if nino1:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 1: " + nino1 + BackgroundColors.normal)
    if line1:
        print(color + str().ljust(nbr_spaces + 5) + line1 + BackgroundColors.normal)
    # second variable
    if file2:
        print(color + str().ljust(nbr_spaces + 5) + "file name 2: " + file2 + BackgroundColors.normal)
    if var2:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 2: " + var2 + BackgroundColors.normal)
    if axes2:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 2: " + axes2 + BackgroundColors.normal)
    if time2:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 2: " + time2 + BackgroundColors.normal)
    if shape2:
        print(color + str().ljust(nbr_spaces + 5) + "shape 2: " + shape2 + BackgroundColors.normal)
    if nina2:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 2: " + nina2 + BackgroundColors.normal)
    if nino2:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 2: " + nino2 + BackgroundColors.normal)
    if line2:
        print(color + str().ljust(nbr_spaces + 5) + line2 + BackgroundColors.normal)
    # third variable
    if file3:
        print(color + str().ljust(nbr_spaces + 5) + "file name 3: " + file3 + BackgroundColors.normal)
    if var3:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 3: " + var3 + BackgroundColors.normal)
    if axes3:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 3: " + axes3 + BackgroundColors.normal)
    if time3:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 3: " + time3 + BackgroundColors.normal)
    if shape3:
        print(color + str().ljust(nbr_spaces + 5) + "shape 3: " + shape3 + BackgroundColors.normal)
    if nina3:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 3: " + nina3 + BackgroundColors.normal)
    if nino3:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 3: " + nino3 + BackgroundColors.normal)
    if line3:
        print(color + str().ljust(nbr_spaces + 5) + line3 + BackgroundColors.normal)
    # fourth variable
    if file4:
        print(color + str().ljust(nbr_spaces + 5) + "file name 4: " + file4 + BackgroundColors.normal)
    if var4:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 4: " + var4 + BackgroundColors.normal)
    if axes4:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 4: " + axes4 + BackgroundColors.normal)
    if time4:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 4: " + time4 + BackgroundColors.normal)
    if shape4:
        print(color + str().ljust(nbr_spaces + 5) + "shape 4: " + shape4 + BackgroundColors.normal)
    if nina4:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 4: " + nina4 + BackgroundColors.normal)
    if nino4:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 4: " + nino4 + BackgroundColors.normal)
    if line4:
        print(color + str().ljust(nbr_spaces + 5) + line4 + BackgroundColors.normal)
    # fifth variable
    if file5:
        print(color + str().ljust(nbr_spaces + 5) + "file name 5: " + file5 + BackgroundColors.normal)
    if var5:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 5: " + var5 + BackgroundColors.normal)
    if axes5:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 5: " + axes5 + BackgroundColors.normal)
    if time5:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 5: " + time5 + BackgroundColors.normal)
    if shape5:
        print(color + str().ljust(nbr_spaces + 5) + "shape 5: " + shape5 + BackgroundColors.normal)
    if nina5:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 5: " + nina5 + BackgroundColors.normal)
    if nino5:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 5: " + nino5 + BackgroundColors.normal)
    if line5:
        print(color + str().ljust(nbr_spaces + 5) + line5 + BackgroundColors.normal)
    # sixth variable
    if file6:
        print(color + str().ljust(nbr_spaces + 5) + "file name 6: " + file6 + BackgroundColors.normal)
    if var6:
        print(color + str().ljust(nbr_spaces + 5) + "variable name 6: " + var6 + BackgroundColors.normal)
    if axes6:
        print(color + str().ljust(nbr_spaces + 5) + "axes list 6: " + axes6 + BackgroundColors.normal)
    if time6:
        print(color + str().ljust(nbr_spaces + 5) + "time bounds 6: " + time6 + BackgroundColors.normal)
    if shape6:
        print(color + str().ljust(nbr_spaces + 5) + "shape 6: " + shape6 + BackgroundColors.normal)
    if nina6:
        print(color + str().ljust(nbr_spaces + 5) + "nina year 6: " + nina6 + BackgroundColors.normal)
    if nino6:
        print(color + str().ljust(nbr_spaces + 5) + "nino year 6: " + nino6 + BackgroundColors.normal)
    if line6:
        print(color + str().ljust(nbr_spaces + 5) + line6 + BackgroundColors.normal)
# ---------------------------------------------------------------------------------------------------------------------#
