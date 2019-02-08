# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Main driver to plot outputs from the ENSO_metrics package
#---------------------------------------------------#


#---------------------------------------------------#
# usual python package
from getpass import getuser as GETPASSgetuser
from os.path import join as OSpath__join
import sys
# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# ENSO_metrics package
# set new path where to find programs
sys.path.insert(0, "/home/" + user_name + "/ENSO_metrics/lib")
sys.path.insert(1, "/home/" + user_name + "/ENSO_metrics/scripts")
from EnsoCollectionsLib import defCollection
from Lib_plot_on_ciclad import metric_boxplot
# My (YYP) package
# set new path where to find programs
sys.path.insert(0, "/home/" + user_name + "/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import get_ensembles
#---------------------------------------------------#


#---------------------------------------------------#
# general parameters (do not change this)
# where to save the files
path = "/data/" + user_name + "/ENSO_metrics"
# metric collection
mc_name = 'ENSO_perf'
# get metric collection definition
dict_mc = defCollection(mc_name)
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
model = 'IPSL-CM5A-LR'
# metric
metric = 'EnsoAmpl'
# list of ensembles
list_ensembles = sorted(get_ensembles(exp=experiment, fre=frequency, mod=model, pro=project, rea=realm))
#---------------------------------------------------#


#---------------------------------------------------#
# plot
# path for the plots
path_out = "/data/" + user_name + "/Plots"  # create your own
# loop on ensembles
# only if you have computed all available ensembles! If not, set 'list_ensembles' with the list of ensembles you have
# computed: e.g., list_ensembles = ['r1i1p1']
for ens in list_ensembles:
    # input pattern or file name
    # this is an example. You can create plots for all 'nbr_years' you have computed using the '*' before 'years' or
    # specify the number of years: e.g., '_0010years_'
    # You can also give a file name 'in hard': 'yplanton_ENSO_perf_IPSL-CM5A-LR_pi_r1i1p1_0010years_EnsoAmpl'
    input_name = OSpath__join(
        path, user_name + "_" + mc_name + "_" + model + "_" + experiment + "_" + ens + "_*years_" + metric)
    # output file name
    output_name = OSpath__join(path_out, mc_name + "_" + model + "_" + experiment + "_" + ens + "_" + metric)
    #---------------------------------------------------#
    # boxplot (uses json files)
    for key_val in ['metric', 'metric_err', 'diagnostic', 'diagnostic_err']:
        print key_val
        plot_title = model + " " + key_val  # this is an example you can give the name you want
        local_output_name = output_name + "_" + key_val  # this is an example you can give the name you want
        # this must be changed depending on the metric. Read the metadata in the json file or ask me for ideas
        if 'metric' in key_val:
            plot_yname = "STDmodel / STDobs"
        else:
            plot_yname = "STD"
        metric_boxplot(input_name, key_val, local_output_name, plot_title, plot_yname)


