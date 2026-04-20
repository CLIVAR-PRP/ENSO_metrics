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
    recipe = "bias_sst_lon"
    region1_main = None  # "equatorial_pacific"  # -> this is not needed, use this to change the default region
    statistic = "average"
    variable1 = None  # "ts"  # -> this is not needed, use this to change the default variable
    path_output = "/Users/yplanton-admin/Documents/Data/Test"
    kwargs_saver = {
        # -> define output filename and path
        # Note that if filename_netcdf is not provided (must be unique each time the 'diagnostic' function is called),
        # the package will create a filename based on the keys project, dataset, experiment, member provided when
        # calling the 'diagnostic' function
        # Note that the output path can be provided in 'filename_netcdf' and 'path' set to None
        "filename_netcdf": None,
        "path": path_output}
    kwargs = {}
    # -- model data
    project = "cmip6"
    dataset = "CanESM5-1"
    experiment = "historical"
    grid = "gn"
    path = "/Users/yplanton-admin/Documents/Data/%s/*/%s" % (project.upper(), dataset)
    area1_name = "areacella"
    area1_file = path + "/%s_fx_%s_*_*_%s.nc" % (area1_name, dataset, grid)
    area2_name = "areacello"
    area2_file = path + "/%s_Ofx_%s_*_*_%s.nc" % (area2_name, dataset, grid)
    mask_name = "sftlf"
    mask_file = path + "/%s_fx_%s_*_*_%s.nc" % (mask_name, dataset, grid)
    # -> nhf and ssh not need for the recipe unless you change variable1
    vari_name_nhf = ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"]
    vari_name_ssh = "zos"
    vari_name_ts = "ts"
    # just an example of dictionary with multiple members
    dict_data = {project: {dataset: {}}, "observations": {}}
    for mem in ["r1i1p1f1", "r2i1p1f1"]:
        # filenames
        vari_file_nhf = [path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (k, dataset, experiment, mem, grid)
                         for k in vari_name_nhf]
        vari_file_ssh = path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (vari_name_ssh, dataset, experiment, mem, grid)
        vari_file_ts = path + "/%s_?mon_%s_%s_%s_%s_*.nc" % (vari_name_ts, dataset, experiment, mem, grid)
        # model data input dictionary
        # this dictionary is structured as required by the package, anything set to None is no required
        dict_t = {
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
        # save in dictionary
        dict_data["model"][dataset][mem] = dict_t
    # -- Observations data
    project_r = "observations"
    experiment_r = "historical"
    member_r = "r1i1p1f1"
    for dat in ["ERSSTv5", "COBE2"]:
        # path and filename
        path_r = "/Users/yplanton-admin/Documents/Data/%s/%s" % (str(project_r.upper()[0]) + str(project_r[1:]), dat)
        vari_name_ts = "sst"
        vari_file_ts = path_r + "/*%s*.nc" % vari_name_ts
        # observations data input dictionary
        dict_t = {
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
        # save in dictionary
        dict_data[project_r][dat] = {member_r: dict_t}
    # dict_data is structured like:
    # dict_data = {
    #     "cmip6": {
    #         "model_name_1": {
    #             "member_name_1": "<dictionary of variables structured for the package>",
    #             ..., -> e.g., other members
    #         },
    #         ..., -> e.g., other models
    #     },
    #     "observations": {
    #         "reference_name_1": {
    #             "reference_name_1": "<dictionary of variables structured for the package>",
    #         },
    #         ..., -> e.g., other reference
    #     },
    # }
    # The only part very important here is "<dictionary of variables structured for the package>"
    # All the rest is defined by the user as one wants
    #
    # -- Compute diagnostic
    #
    # GOAL: compute diagnostics and save them in a netCDF (one file per input)
    # loop on dataset_type, dataset, member to compute everything in dict_data as defined by the user
    print("diagnostic", recipe)
    for dataset_type, d1 in dict_data.items():
        for dataset, d2 in d1.items():
            for member, d3 in d2.items():
                print(str(dataset_type).rjust(20), str(dataset).ljust(15), member)
                available_recipe[recipe].diagnostic(
                    d3, dataset=dataset, experiment=experiment, project=dataset_type, member=member,
                    region1_main=region1_main, supplementary=True, variable1=variable1,
                    kwargs_saver=kwargs_saver, **kwargs)
    # -- Compute metric
    # GOAL: compute the metric (distance model-observations)
    # Here the user must provide the netCDF outputs from the diagnostic function
    # loop on model & member, they will be evaluated against the observations
    print("metric", recipe)
    metric_dict = {}
    for dataset_type, d1 in dict_data.items():
        if dataset_type != project:
            continue
        for dataset, d2 in d1.items():
            for member, d3 in d2.items():
                # Here is the difficulty, if you haven't defined filename_netcdf to the diagnostic function (in
                # kwargs_saver), you need to find the files created by the package, here is an example of filename
                # generated by the package
                # model file
                input_dataset = path_output + "/%s_%s_%s_%s_%s.nc" % (dataset_type, dataset, experiment, member, recipe)
                # reference(s) file: this dictionary contains one or more references for the evaluation
                input_reference = {}
                for dataset_type_r, di1 in dict_data.items():
                    if dataset_type_r != project_r:
                        continue
                    for dataset_r, di2 in di1.items():
                        for member_r, di3 in di2.items():
                            input_reference[dataset_r] = path_output + "/%s_%s_%s_%s_%s.nc" % (
                                dataset_type_r, dataset_r, experiment, member_r, recipe)
                # keys to store the metric values
                keys = (recipe, dataset_type, dataset, experiment, member)
                # compute metric values
                metric_dict = available_recipe[recipe].metric(
                    input_dataset, input_reference, metric_dictionary=metric_dict, metric_dictionary_keys=keys)
    print(json__dumps(metric_dict, indent=4))
    #
    # -- Plot
    #
    # not done yet
# ---------------------------------------------------------------------------------------------------------------------#
