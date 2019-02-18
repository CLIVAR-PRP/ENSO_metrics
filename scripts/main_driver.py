# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Main driver to compute the ENSO_metrics package
# This version of the package has been created to be able:
#      - to compute one metric at a time
#      - to compute the metric using only a specified number of simulated (model) years
#---------------------------------------------------#


#---------------------------------------------------#
# usual python package
from getpass import getuser as GETPASSgetuser
from os import makedirs as OSmakedirs
from os.path import isdir as OSpath__isdir
import sys
# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# ENSO_metrics package
# set new path where to find programs
sys.path.insert(0, "/home/" + user_name + "/ENSO_metrics/lib")
sys.path.insert(1, "/home/" + user_name + "/ENSO_metrics/scripts")
from EnsoCollectionsLib import defCollection
from Lib_CCollection_cmip_on_ciclad import main_compute
from Lib_CCollection_cmip_on_ciclad import nbryear_from_model
from my_arg import my_arg
# My (YYP) package
# set new path where to find programs
sys.path.insert(0, "/home/" + user_name + "/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import get_ensembles
from getfiles_sh_to_py import get_models
#---------------------------------------------------#


#---------------------------------------------------#
# general parameters (do not change this)
# where to save the files
path = "/data/" + user_name + "/ENSO_metrics"
# test if the directory exists, if not, creates it
if OSpath__isdir(path) is False:
    OSmakedirs(path)
# metric collection
mc_name = 'ENSO_perf'
# get metric collection definition
dict_mc = defCollection(mc_name)
# list of metrics of the metric collection
list_metrics = sorted(dict_mc['metrics_list'].keys())
#---------------------------------------------------#


#---------------------------------------------------#
# parameters
# At least in the beginning, you are not supposed to change the next 4 lines
# CMIP parameters
project = 'CMIP5'
experiment = 'pi'
frequency = 'mon'
realm = 'A'
# The next lines can be modified
# model
default_model = 'IPSL-CM5A-LR'
# metric (must be in list_metrics)
default_metric = 'EnsoAmpl'
# number of years used for the computation of the metric (must be an integer multiple of 5 between 5 and 500)
default_nyears = 100
# save output NetCDF files ('False' or 'True')
default_savenc = 'False'
# file name
file_name = user_name + "_" + mc_name + "_" + model + "_" + experiment
#---------------------------------------------------#


#---------------------------------------------------#
# arguments: metric, number or years used and save netcdf (do not change this)
needed_arg = {
    # define what arguments can be 'metric'
    'metric': {
        # 'metric' must be a string
        'type': basestring,
        # 'metric' must be define in EnsoCollectionsLib.defCollection(mc_name)
        'possible_values': list_metrics,
        # the default value of 'metric' is:
        'default': default_metric,
    },
    # define what arguments can be 'model'
    'model': {
        # 'model' must be a string
        'type': basestring,
        # 'model' must be define in getfiles_sh_to_py.get_models(project)
        'possible_values': sorted(get_models(project)),
        # the default value of 'model' is:
        'default': default_model,
    },
    # define what arguments can be 'nbr_years'
    'nbr_years': {
        # 'nbr_years' must be an integer
        'type': 'int',
        # 'nbr_years' must be a multiple of 5 between 5 and 500
        'possible_values': range(5, 505, 5),
        # the default value of 'nbr_years' is:
        # you can change this integer if you do not use main_driver.sh and if you do not want to give arguments when
        # running main_driver.py (e.g., python -i main_driver.py 10 True
        'default': default_nyears,
    },
    # define what arguments can be 'save_netcdf' (i.e., saving or not NetCDF files, used to create maps and curves)
    'save_netcdf': {
        # 'nbr_years' must be a string
        'type': basestring,
        # 'save_netcdf' must be either False or True
        'possible_values': ['False', 'True'],
        # the default value of 'save_netcdf' is:
        # False means NetCDF files are not saved, True means NetCDF files are saved
        'default': default_savenc,
    },
}
dict1 = my_arg(needed_arg, sys.argv[1:])
metric = dict1['metric']
nbr_years = dict1['nbr_years']
save_netcdf = dict1['save_netcdf']
#---------------------------------------------------#


#---------------------------------------------------#
# information
# list of CMIP models
list_models = sorted(get_models(project, experiment))
# list of ensembles
list_ensembles = sorted(get_ensembles(exp=experiment, fre=frequency, mod=model, pro=project, rea=realm))


# to get more information about the models or the metrics, you can uncomment the next lines

# to print the number of years of all ensembles for a given experiment/frequency/model/project/realm, uncomment next
# lines
# for ens in list_ensembles:
#     ens_nbr_years = nbryear_from_model(experiment, ens, frequency, model, project, realm, 'sst')
#     print ens + ': ' + str(ens_nbr_years).rjust(4) + ' years'

# to print all metrics for the given metric collection, uncomment next line
# print list_metrics

# to print all CMIP models available, uncomment next line
# print list_models

# to print all ensembles for the given experiment/frequency/model/project/realm, uncomment next line
# print list_ensembles

# to print all ensembles for all models for a given experiment/frequency/project/realm, uncomment next lines
# for mod in list_models:
#     print mod.ljust(14) + ': ' + str(get_ensembles(exp=experiment, fre=frequency, mod=mod, pro=project, rea=realm))

# to print the number of years of all ensembles for all models for a given experiment/frequency/project/realm, uncomment
# next lines
# for mod in list_models:
#     list_ens = sorted(get_ensembles(exp=experiment, fre=frequency, mod=mod, pro=project, rea=realm))
#     for ens in list_ens:
#         ens_nbr_years = nbryear_from_model(experiment, ens, frequency, mod, project, realm, 'sst')
#         print mod.ljust(14) + ' ' + ens + ': ' + str(ens_nbr_years).rjust(4) + ' years'
#---------------------------------------------------#


#---------------------------------------------------#
# to compute the metric for a given experiment/frequency/model/project/realm
main_compute(mc_name, metric, nbr_years, path, file_name, experiment, frequency, model, project, realm,
             save_netcdf=save_netcdf)
