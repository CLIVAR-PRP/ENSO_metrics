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
    recipe = "stat_box"
    region = "nino3"
    statistic = "average"
    variable = "ts"
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
    dict_model = {
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
    kwargs = {}
    print(variable)
    print(json__dumps(dict_model, indent=4))
    print(list(available_recipe.keys()))
    if recipe in list(available_recipe.keys()):
        available_recipe[recipe].diagnostic(
            dict_model, dataset=dataset, experiment=experiment, project=project, member=member, region=region,
            statistic=statistic, variable=variable, **kwargs)
        print("computed")
    print("done")
# ---------------------------------------------------------------------------------------------------------------------#
