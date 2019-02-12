# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Parses arguments for main_driver.py
#---------------------------------------------------#


#---------------------------------------------------#
# import python packages
# usual python package
from copy import deepcopy
#---------------------------------------------------#


def test_type(arg_list, arg_type, arg_values, arg_default):
    if arg_type == 'int':
        type_in = int
    elif arg_type == 'float':
        type_in = float
    else:
        type_in = deepcopy(arg_type)
    for ii in range(len(arg_list)):
        tmp = arg_list[ii]
        if arg_type == 'int':
            try:
                tmp = int(tmp)
            except:
                pass
        elif arg_type == 'float':
            try:
                tmp = float(tmp)
            except:
                pass
        if isinstance(tmp, type_in):
            if tmp in arg_values:
                arg_out = deepcopy(tmp)
                break
    try:
        arg_out
    except:
        arg_out = deepcopy(arg_default)
    return arg_out


def my_arg(needed_arg, *kwargs):
    dict_out = dict()
    for key in needed_arg.keys():
        tmp = needed_arg[key]
        tmp_out = test_type(kwargs[0], tmp['type'], tmp['possible_values'], tmp['default'])
        if tmp_out == 'False':
            tmp_out = False
        elif tmp_out == 'True':
            tmp_out = True
        dict_out[key] = tmp_out
        print ''
    return dict_out

