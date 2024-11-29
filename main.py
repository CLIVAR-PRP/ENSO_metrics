# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Program to test metric computation
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from inspect import getmembers, isfunction

# local functions
from enso_metrics import recipes
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# call metric
# ---------------------------------------------------------------------------------------------------------------------#
available_recipe = dict((k[0], k[1]) for k in getmembers(recipes, isfunction))

if __name__ == '__main__':
    path = "/Users/ypla0001/Documents/Data/CMIP6/E3SM-Project/E3SM-2-0"
    vari_file = path + "/ts_Amon_E3SM-2-0_historical_r1i1p1f1_gr_185001-201412.nc"
    vari_name = "ts"
    area_file = path + "/areacella_fx_E3SM-2-0_historical_r1i1p1f1_gr.nc"
    area_name = "areacella"
    mask_file = path + "/sftlf_fx_E3SM-2-0_historical_r1i1p1f1_gr.nc"
    mask_name = "sftlf"
    dict_model = {
        "areacella": {
            "file_name": area_file,
            "variable": area_name,
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
        "ts":{
            "area": "areacella",
            "mask": "landmask",
            "file_name": vari_file,
            "variable": vari_name,
            "variable_computation": "ts - 273.15",
            "variable_offset": -273.15,
            "variable_scaling": None,
        },
    }
    project = "cmip6"
    dataset = "E3SM-2-0"
    experiment = "historical"
    member = "r1i1p1f1"
    region = "nino3"
    statistic = "average"
    variable = "ts"
    recipe = "stat_box"
    kwargs = {}
    print(vari_name)
    if recipe in list(available_recipe.keys()):
        available_recipe[recipe](
            dict_model, dataset=dataset, experiment=experiment, project=project, member=member, region=region,
            statistic=statistic, variable=variable, **kwargs)
        print("computed")
    print("done")
# ---------------------------------------------------------------------------------------------------------------------#
