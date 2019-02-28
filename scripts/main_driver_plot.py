# -*- coding:UTF-8 -*-
#---------------------------------------------------#
# Main driver to plot outputs from the ENSO_metrics package
#---------------------------------------------------#


#---------------------------------------------------#
# usual python package
from copy import deepcopy
from getpass import getuser as GETPASSgetuser
from os.path import join as OSpath__join
import sys
# user (get your user name for the paths and to save the files)
user_name = GETPASSgetuser()
# ENSO_metrics package
# set new path where to find programs
sys.path.insert(0, "/home/" + user_name + "/ENSO_metrics/lib")
sys.path.insert(1, "/home/" + user_name + "/ENSO_metrics/scripts")
from Lib_plot_on_ciclad import metric_boxplot, myplot
# My (YYP) package
# set new path where to find programs
sys.path.insert(0, "/home/yplanton/New_programs/lib_cmip_bash")
from getfiles_sh_to_py import get_ensembles
#---------------------------------------------------#


#---------------------------------------------------#
# general parameters (do not change this)
# where to save the files
path = "/data/" + user_name + "/ENSO_metrics"
# path = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Test"
# metric collection
mc_name = 'ENSO_perf'
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
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/data/Plots"
# loop on ensembles
# only if you have computed all available ensembles! If not, set 'list_ensembles' with the list of ensembles you have
# computed: e.g., list_ensembles = ['r1i1p1']
for ens in list_ensembles:
    # input pattern or file name
    # this is an example. You can create plots for all 'nbr_years' you have computed using the '*' before 'years' or
    # specify the number of years: e.g., '_0010years_'
    # You can also give a file name 'in hard': 'yplanton_ENSO_perf_IPSL-CM5A-LR_pi_r1i1p1_0010years_EnsoAmpl'
    tmp_name = mc_name + "_" + model + "_" + experiment + "_" + ens + "_*years_" + metric
    input_name = OSpath__join(path, user_name + "_" + tmp_name)
    #input_name = OSpath__join(path, "yplanton_" + tmp_name)
    print input_name
    # output file name
    output_name = OSpath__join(
        path_out, user_name + "_" + mc_name + "_" + model + "_" + experiment + "_" + ens + "_" + metric)
    # general plot title
    title = model + " " + metric
    #---------------------------------------------------#
    # boxplot (uses json files)
    print str().ljust(5) + "boxplot"
    if "Rmse" in metric:
        listkeys = ['metric']
    else:
        listkeys = ['metric', 'metric_err', 'diagnostic', 'diagnostic_err']
    for key_val in listkeys:
        print str().ljust(10) + key_val
        plot_title = title + " (" + key_val + ")"  # this is an example you can give the name you want
        local_output_name = output_name + "_" + key_val  # this is an example you can give the name you want
        # this must be changed depending on the metric. Read the metadata in the json file or ask me for ideas
        if "Rmse" in metric:
            plot_yname = "RMSE (model, obs)"
        elif metric == 'EnsoAmpl':
            if 'metric' in key_val:
                plot_yname = "STDmodel / STDobs"
            else:
                plot_yname = "sst standard deviation"
        elif metric == "EnsoDiversity":
            if 'metric' in key_val:
                plot_yname = "DIVmodel / DIVobs"
            else:
                plot_yname = "% EP La Nina / % EP El Nino"
        elif metric == "EnsoSeasonality":
            if 'metric' in key_val:
                plot_yname = "SEAmodel / SEAobs"
            else:
                plot_yname = "STD(NDJ) / STD(MAM)"
        elif metric == 'EnsoSstSkew':
            if 'metric' in key_val:
                plot_yname = "SKEWmodel / SKEWobs"
            else:
                plot_yname = "sst skewness"
        elif metric == "NinaSstDiv":
            if 'metric' in key_val:
                plot_yname = "% EP La Nina model / obs"
            else:
                plot_yname = "% EP La Nina"
        elif metric == "NinaSstDur":
            if 'metric' in key_val:
                plot_yname = "duration La Nina model / obs"
            else:
                plot_yname = "duration La Nina"
        elif metric == "NinoSstDiv":
            if 'metric' in key_val:
                plot_yname = "% EP El Nino model / obs"
            else:
                plot_yname = "% EP El Nino"
        elif metric == "NinoSstDur":
            if 'metric' in key_val:
                plot_yname = "duration El Nino model / obs"
            else:
                plot_yname = "duration EP El Nino"
        metric_boxplot(input_name, key_val, local_output_name, plot_title, plot_yname)
        del local_output_name, plot_title, plot_yname
    #---------------------------------------------------#
    # curveplot (uses 1d.nc files)
    print str().ljust(5) + "curveplot"
    myplot(input_name, output_name, metric, title, '1d')
    # ---------------------------------------------------#
    # mapplot (uses 1d.nc files)
    print str().ljust(5) + "mapplot"
    myplot(input_name, output_name, metric, title, '2d', model_name=model, z_name='spread 90%')


