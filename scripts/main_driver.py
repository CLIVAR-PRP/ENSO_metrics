# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Main driver to compute the ENSO_metrics package
# This version of the package has been created to be able:
#      - to compute one metric at a time
#      - to compute the metric using only a specified number of simulated (model) years
#---------------------------------------------------#


#---------------------------------------------------#
# usual python package
import datetime
import sys
# ENSO_metrics package
# set new path where to find programs
sys.path.insert(0, "/home/yplanton/ENSO_metrics/lib")
from EnsoCollectionsLib import defCollection
from Lib_CCollection_cmip_on_ciclad import main_compute
from Lib_CCollection_cmip_on_ciclad import nbryear_from_model
# My (YYP) package
# set new path where to find programs
sys.path.insert(0, "/home/yplanton/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import get_ensembles
from getfiles_sh_to_py import get_models
#---------------------------------------------------#


#---------------------------------------------------#
# At least in the beginning, you are not supposed to change the next lines
# metric collection
mc_name = 'ENSO_perf'
dict_mc = defCollection(mc_name)
# list of metrics
list_metrics = sorted(dict_mc['metrics_list'].keys())
# CMIP parameters
project = 'CMIP5'
experiment = 'pi'
frequency = 'mon'
realm = 'A'
list_models = sorted(get_models(project))
# today's date: used to save files
today = str(datetime.date.today())
today = today.replace("-", "")
#---------------------------------------------------#


#---------------------------------------------------#
# These lines can be modified
# model
model = 'IPSL-CM5A-LR'
# metric
metric = 'EnsoAmpl'
# list of ensembles
list_ensembles = sorted(get_ensembles(exp=experiment, fre=frequency, mod=model, pro=project, rea=realm))
# number of years used
nbr_years = 10 # you can test from 5 to 200 years
# where to save the files
path = '/data/yplanton/ENSO_metrics' # change that to your own environment
# file name
file_name = today + '_USER_NAME_' + mc_name + '_' + model + '_' + experiment

# to print all metrics for the given metric collection, uncomment next line
# print list_metrics

# to print all CMIP models available, uncomment next line
# print list_models

# to print all ensembles for all models for a given experiment/frequency/project/realm, uncomment next lines
# for mod in list_models:
#     print mod + ': ' + str(get_ensembles(exp=experiment, fre=frequency, mod=mod, pro=project, rea=realm))

# to print all ensembles for the given experiment/frequency/model/project/realm, uncomment next line
# print list_ensembles

# to print the number of years of all ensembles for a given experiment/frequency/model/project/realm
# for ens in list_ensembles:
#     ens_nbr_years = nbryear_from_model(experiment, ens, frequency, model, project, realm, 'sst')
#     print ens + ': ' + str(ens_nbr_years).rjust(4) + ' years'
#---------------------------------------------------#


#---------------------------------------------------#
# to compute the metric for a given experiment/frequency/model/project/realm
main_compute(mc_name, metric, nbr_years, path, file_name, experiment, frequency, model, project, realm)

