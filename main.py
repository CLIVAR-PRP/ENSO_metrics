# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Program to test metric computation
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from inspect import getmembers, ismodule
from json import dumps as json__dumps

# local functions
from enso_metrics import recipes
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# call metric
# ---------------------------------------------------------------------------------------------------------------------#
available_recipe = dict((k[0], k[1]) for k in getmembers(recipes, ismodule) if "__" not in k)


if __name__ == '__main__':
    # -- metric arguments
    recipe = "bias_sst_lon"  # "stat_box"
    region1_main = "equatorial_pacific"  # "nino3"
    statistic = "average"
    variable1 = "ts"
    path_output = "/Users/yplanton-admin/Documents/Data/Test"
    kwargs = {}
    # -- model data
    project = "cmip6"
    dataset = "CanESM5-1"
    experiment = "historical"
    member = "r1i1p1f1"
    grid = "gn"
    path = "/Users/yplanton-admin/Documents/Data/%s/*/%s" % (project.upper(), dataset)
    area1_name = "areacella"
    area1_file = path + "/%s_fx_%s_*_*_%s.nc" % (area1_name, dataset, grid)
    area2_name = "areacello"
    area2_file = path + "/%s_Ofx_%s_*_*_%s.nc" % (area2_name, dataset, grid)
    mask_name = "sftlf"
    mask_file = path + "/%s_fx_%s_*_*_%s.nc" % (mask_name, dataset, grid)
    vari_name_nhf = ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"]
    vari_file_nhf = [path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (k, dataset, experiment, member, grid) for k in vari_name_nhf]
    vari_name_ssh = "zos"
    vari_file_ssh = path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (vari_name_ssh, dataset, experiment, member, grid)
    vari_name_ts = "ts"
    vari_file_ts = path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (vari_name_ts, dataset, experiment, member, grid)
    # model data input dictionary
    dict_data = {
        "areacella": {
            "file_name": area1_file,
            "variable": area1_name,
            "variable_computation": None,
            "variable_offset": None,
            "variable_scaling": None,
        },
        "areacello": {
            "file_name": area2_file,
            "variable": area2_name,
            "variable_computation": None,
            "variable_offset": None,
            "variable_scaling": None,
        },
        "landmask": {
            "file_name": mask_file,
            "variable": mask_name,
            "variable_computation": "1e-2 * sftlf",
            "variable_offset": None,
            "variable_scaling": 1e-2,
        },
        "nhf": {
            "area": "areacella",
            "mask": "landmask",
            "file_name": vari_file_nhf,
            "variable": vari_name_nhf,
            "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
            "variable_offset": None,
            "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
        },
        "ssh": {
            "area": "areacello",
            "mask": None,
            "file_name": vari_file_ssh,
            "variable": vari_name_ssh,
            "variable_computation": "1e2 * zos",
            "variable_offset": None,
            "variable_scaling": 1e2,
        },
        "ts": {
            "area": "areacella",
            "mask": "landmask",
            "file_name": vari_file_ts,
            "variable": vari_name_ts,
            "variable_computation": "ts - 273.15",
            "variable_offset": -273.15,
            "variable_scaling": None,
        },
    }
    # -- Observations data
    project_r = "observations"
    dataset_r = "COBE2"  # "ERSSTv5"  # "HadISST"  #
    experiment_r = "historical"
    member_r = "r1i1p1f1"
    path_r = "/Users/yplanton-admin/Documents/Data/%s/%s" % (str(project_r.upper()[0]) + str(project_r[1:]), dataset_r)
    vari_name_ts = "sst"
    vari_file_ts = path_r + "/*%s*.nc" % vari_name_ts
    # observations data input dictionary
    dict_reference = {
        "ts": {
            "area": None,
            "mask": None,
            "file_name": vari_file_ts,
            "variable": vari_name_ts,
            "variable_computation": "ts",
            "variable_offset": None,
            "variable_scaling": None,
        },
    }
    # -- Compute diagnostic
    print(variable1)
    # print(json__dumps(dict_model, indent=4))
    print(list(available_recipe.keys()))
    if recipe in list(available_recipe.keys()):
        print(str().ljust(5), "diagnostic", recipe)
        available_recipe[recipe].diagnostic(
            dict_data, dataset=dataset, experiment=experiment, project=project, member=member,
            supplementary=True, variable1=variable1, kwargs_saver={"path": path_output}, **kwargs)
        # available_recipe[recipe].diagnostic(
        #     dict_reference, dataset=dataset_r, experiment=experiment_r, project=project_r, member=member_r,
        #     supplementary=True, variable1=variable1, kwargs_saver={"path": path_output}, **kwargs)
        print("computed")
    # stop
    # -- Compute metric
    if recipe in list(available_recipe.keys()):
        print(str().ljust(5), "metric", recipe)
        # model to evaluate
        input_dataset = path_output + "/%s_%s_%s_%s_%s.nc" % (project, dataset, experiment, member, recipe)
        # reference(s) to use
        input_reference = {
            "COBE2": path_output + "/%s_%s_%s_%s_%s.nc" % (project_r, "COBE2", experiment_r, member_r, recipe),
            "ERSSTv5": path_output + "/%s_%s_%s_%s_%s.nc" % (project_r, "ERSSTv5", experiment_r, member_r, recipe),
            "HadISST": path_output + "/%s_%s_%s_%s_%s.nc" % (project_r, "HadISST", experiment_r, member_r, recipe),
        }
        # dictionary to save the metric values
        metric_dict = {}
        # keys to store the metric values
        keys = (recipe, project, dataset, experiment, member)
        # compute
        metric_dict = available_recipe[recipe].metric(
            input_dataset, input_reference, metric_dictionary=metric_dict, metric_dictionary_keys=keys)
        print(json__dumps(metric_dict, indent=4))
    print("done")
# ---------------------------------------------------------------------------------------------------------------------#
