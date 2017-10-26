# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Import packages
#---------------------------------------------------#
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
#---------------------------------------------------#
def PlusCommaSpace(string):
    """
    #################################################################################
	Description:  Adds a comma and a space if the string is not empty or if the
                  string is not composed of space only
    #################################################################################

    :param string (string):
        string to which comma_space could be added
    :return string (string):
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
    if string.isspace(): return string
    elif not string:     return string
    else:                return string+', '

#---------------------------------------------------#
def MessageFormating(inspect_stack):
    """
    #################################################################################
    Description:  Formats inspect.stack() as '   File filename, line n, in module'
    #################################################################################

    :param inspect_stack (array):
        list of informations about the program/module/line,... created using inspect.stack()
    :return string (string):
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
        string = plus_comma_space(string)+'File '+str(inspect_stack[0][1])
    # adds line number
    string = plus_comma_space(string)+'line '+str(inspect_stack[0][2])
    # adds module's name
    if inspect_stack[0][3] != '<module>':
        string = plus_comma_space(string)+'in '+str(inspect_stack[0][3])
    return string

#---------------------------------------------------#
def UnlikelyUnits(var_name, name_in_file, units, inspect_stack):
    list_strings = ["ERROR " + MessageFormating(inspect_stack) + ": units",
                    str().ljust(5) + "the file says that " + var_name + " (" + name_in_file + ") is in " + units
                    + " but it seems unlikely"]
    MyError(list_strings)
    return
#---------------------------------------------------#
def UnknownUnits(var_name, name_in_file, units, inspect_stack):
    list_strings = ["ERROR" + MessageFormating(inspect_stack) + ": units",
                    str().ljust(5) + "unknown units: " + var_name + " (" + name_in_file + ") is in " + units]
    MyError(list_strings)
    return
#---------------------------------------------------#
def MyWarning(list_strings):
    for ii in range(2): print bcolors.WARNING + ""
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        print str().ljust(5) + string
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for ii in range(2): print '' + bcolors.ENDC
    sys_exit(1)
    return
#---------------------------------------------------#
def MyError(list_strings):
    for ii in range(2): print bcolors.FAIL + ""
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for string in list_strings:
        print str().ljust(5) + string
    print str().ljust(5) + "%%%%%     -----     %%%%%"
    for ii in range(2): print '' + bcolors.ENDC
    sys_exit(1)
    return
