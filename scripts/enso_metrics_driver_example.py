# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Example of driver for the CLIVAR-PRP ENSO metrics package
#
# This examples uses the following datasets:
#     - CMIP6 MPI-ESM1-2-LR historical r1i1p1f1
# https://aims2.llnl.gov/search?project=CMIP6&resultType=originals+only&activeFacets=%7B%22activity_id%22%3A%22CMIP%22%2C%22source_id%22%3A%22MPI-ESM1-2-LR%22%2C%22experiment_id%22%3A%22historical%22%2C%22variant_label%22%3A%22r1i1p1f1%22%2C%22table_id%22%3A%5B%22Amon%22%2C%22Omon%22%2C%22fx%22%2C%22Ofx%22%5D%2C%22variable_id%22%3A%5B%22areacella%22%2C%22areacello%22%2C%22sftlf%22%2C%22hfls%22%2C%22hfss%22%2C%22pr%22%2C%22rlds%22%2C%22rlus%22%2C%22rsds%22%2C%22rsus%22%2C%22tauu%22%2C%22tauv%22%2C%22ts%22%2C%22zos%22%5D%7D
#     - CFSR
# https://aims2.llnl.gov/search?project=CREATE-IP&activeFacets=%7B%22experiment%22%3A%22CFSR%22%2C%22realm%22%3A%5B%22atmos%22%2C%22ocean%22%5D%2C%22time_frequency%22%3A%22mon%22%2C%22product%22%3A%22reanalysis%22%7D
#     - ERA5
# https://aims2.llnl.gov/search?project=CREATE-IP&activeFacets=%7B%22product%22%3A%22reanalysis%22%2C%22experiment%22%3A%22ERA5%22%7D
#     - GPCP-Monthly-3-2
# https://aims2.llnl.gov/search?project=obs4MIPs&activeFacets=%7B%22variable%22%3A%22pr%22%2C%22source_id%22%3A%22GPCP-Monthly-3-2%22%7D
#
# It assumes that a text representation of each dataset's variable has been create using the cdscan utility
# e.g.,
# cdscan -x <PROJECT>_<MODEL>_<EXPERIMENT>_<MEMBER>_<VARIABLE>.xml /path/to/data/*netcdf*filename*pattern*.nc
#
# The name of the text representation of each dataset's variable is
# <PROJECT>_<MODEL>_<EXPERIMENT>_<MEMBER>_<VARIABLE>.xml
# For CMIP6 MPI-ESM1-2-LR historical r1i1p1f1 it is therefore CMIP6_MPI-ESM1-2-LR_historical_r1i1p1f1_<VARIABLE>.xml
# For reference datasets, <PROJECT>, <EXPERIMENT> and <MEMBER> are respectively observation, historical and r1i1p1f1
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Import package
# ---------------------------------------------------------------------------------------------------------------------#
import json
# ENSO_metrics package
from EnsoMetrics.EnsoCollectionsLib import defCollection
from EnsoMetrics.EnsoComputeMetricsLib import ComputeCollection
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Parameters
# ---------------------------------------------------------------------------------------------------------------------#
# path where the input files (text representation of datasets created using the cdscan utility)
path_i = "/path/to/input"
# path where to save data
path_o = "/path/to/output"

# metric collection name
metric_collection = "ENSO_perf"
# to see all metric collections:
# print(sorted(list(defCollection().keys()), key=lambda v: v.upper()))

# project name
project = "CMIP6"
# model name
model = "MPI-ESM1-2-LR"
# experiment name
experiment = "historical"
# member name
member = "r1i1p1f1"

# name dataset (usualy model & member to properly store the dataset and member)
dataset = "_".join([model, member])

# output json file name (metrics values will be saved in a json file)
ouput_json_file = path_o + "/" + "_".join([metric_collection, project, model, experiment, member]) + ".json"

# output netCDF file name (dive-down diagnostics, i.e., additional data to understand the metric value, will be saved in
# a netCDF file)
ouput_netcdf_file = path_o + "/" + "_".join([metric_collection, project, model, experiment, member]) + ".nc"
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Definition of the data to analyze
# ---------------------------------------------------------------------------------------------------------------------#
# here is an example of a dictionary with model and observational datasets
# depending on the metric collection, the package will need some or all these variables
# if the variables are not given for the model or the observational datasets in this dictionary, the related metrics
# will be skipped
input_dataset_dictionary = {
    "model": {
        # MPI-ESM1-2-LR, please download:
        # https://aims2.llnl.gov/search?project=CMIP6&resultType=originals+only&activeFacets=%7B%22activity_id%22%3A%22CMIP%22%2C%22source_id%22%3A%22MPI-ESM1-2-LR%22%2C%22experiment_id%22%3A%22historical%22%2C%22variant_label%22%3A%22r1i1p1f1%22%2C%22table_id%22%3A%5B%22Amon%22%2C%22Omon%22%2C%22fx%22%2C%22Ofx%22%5D%2C%22variable_id%22%3A%5B%22areacella%22%2C%22areacello%22%2C%22sftlf%22%2C%22hfls%22%2C%22hfss%22%2C%22pr%22%2C%22rlds%22%2C%22rlus%22%2C%22rsds%22%2C%22rsus%22%2C%22tauu%22%2C%22tauv%22%2C%22ts%22%2C%22zos%22%5D%7D
        dataset: {
            "lhf": {
                # the CF standard_name is: surface_upward_latent_heat_flux
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_hfls.xml",
                "varname": "hfls",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "lwr": {
                # the CF standard_name is:
                # surface_downwelling_longwave_flux_in_air and surface_upwelling_longwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rlds.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rlus.xml"],
                "varname": ["rlds", "rlus"],
                # the CF standard_name is: cell_area
                # rlds and rlus could be on different grids so you need to give areacell/landmask for each variables 
                "path + filename_area": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml"],
                "areaname": ["areacella", "areacella"],
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml"],
                "landmaskname": ["sftlf", "sftlf"],
            },
            "pr": {
                # the CF standard_name is: precipitation_flux
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_pr.xml",
                "varname": "pr",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "shf": {
                # the CF standard_name is: surface_upward_sensible_heat_flux
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_hfss.xml",
                "varname": "hfss",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "ssh": {
                # the CF standard_name is: sea_surface_height_above_geoid
                #"path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_zos.xml",
                "path + filename": "/home/yplanton/New_XMLDIR/CMIP6/historical/CMIP6_MPI-M_MPI-ESM1-2-LR_historical_r1i1p1f1_Omon_zos_gn_latest.xml",
                "varname": "zos",
                # the CF standard_name is: cell_area
                # ocean grids are often not regular, this file becomes very important!
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacello.xml",
                "areaname": "areacello",
                # variables defined on ocean grid don't have a land-sea mask (or land area fraction)
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "sst": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "swr": {
                # the CF standard_name is:
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rsds.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rsus.xml"],
                "varname": ["rsds", "rsus"],
                # the CF standard_name is: cell_area
                # rsds and rsus could be on different grids so you need to give areacell/landmask for each variables
                "path + filename_area": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml"],
                "areaname": ["areacella", "areacella"],
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml"],
                "landmaskname": ["sftlf", "sftlf"],
            },
            "taux": {
                # the CF standard_name is: surface_downward_eastward_stress
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_tauu.xml",
                "varname": "tauu",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "tauy": {
                # the CF standard_name is: surface_downward_northward_stress
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_tauv.xml",
                "varname": "tauv",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
            "thf": {
                # the CF standard_name is:
                # surface_upward_latent_heat_flux, surface_upward_sensible_heat_flux,
                # surface_downwelling_longwave_flux_in_air, surface_upwelling_longwave_flux_in_air,
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_hfls.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_hfss.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rlds.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rlus.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rsds.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_rsus.xml"],
                "varname": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                # the CF standard_name is: cell_area
                # the variables could be on different grids so you need to give areacell/landmask for each variables
                "path + filename_area": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml"],
                "areaname": ["areacella", "areacella", "areacella", "areacella", "areacella", "areacella"],
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": [
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                    path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml"],
                "landmaskname": ["sftlf", "sftlf", "sftlf", "sftlf", "sftlf", "sftlf"],
            },
            "ts": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/" + "_".join([project, model, experiment, member]) + "_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": path_i + "/" + "_".join([project, model, experiment, member]) + "_areacella.xml",
                "areaname": "areacella",
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": path_i + "/" + "_".join([project, model, experiment, member]) + "_sftlf.xml",
                "landmaskname": "sftlf",
            },
        },
    },
    "observations": {
        # CFSR, please download:
        # https://aims2.llnl.gov/search?project=CREATE-IP&activeFacets=%7B%22experiment%22%3A%22CFSR%22%2C%22realm%22%3A%5B%22atmos%22%2C%22ocean%22%5D%2C%22time_frequency%22%3A%22mon%22%2C%22product%22%3A%22reanalysis%22%7D
        "CFSR": {
            "lhf": {
                # the CF standard_name is: surface_upward_latent_heat_flux
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_hfls.xml",
                "varname": "hfls",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "lwr": {
                # the CF standard_name is:
                # surface_downwelling_longwave_flux_in_air and surface_upwelling_longwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rlds.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rlus.xml"],
                "varname": ["rlds", "rlus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "pr": {
                # the CF standard_name is: precipitation_flux
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_pr.xml",
                "varname": "pr",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "shf": {
                # the CF standard_name is: surface_upward_sensible_heat_flux
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_hfss.xml",
                "varname": "hfss",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "ssh": {
                # the CF standard_name is: sea_surface_height_above_geoid
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_zos.xml",
                "varname": "zos",
                # the CF standard_name is: cell_area
                # ocean grids are often not regular, this file becomes very important!
                "path + filename_area": None,
                "areaname": None,
                # variables defined on ocean grid don't have a land-sea mask (or land area fraction)
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "sst": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "swr": {
                # the CF standard_name is:
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rsds.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rsus.xml"],
                "varname": ["rsds", "rsus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "taux": {
                # the CF standard_name is: surface_downward_eastward_stress
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_tauu.xml",
                "varname": "tauu",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "tauy": {
                # the CF standard_name is: surface_downward_northward_stress
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_tauv.xml",
                "varname": "tauv",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "thf": {
                # the CF standard_name is:
                # surface_upward_latent_heat_flux, surface_upward_sensible_heat_flux,
                # surface_downwelling_longwave_flux_in_air, surface_upwelling_longwave_flux_in_air,
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_CFSR_historical_r1i1p1f1_hfls.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_hfss.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rlds.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rlus.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rsds.xml",
                    path_i + "/observation_CFSR_historical_r1i1p1f1_rsus.xml"],
                "varname": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "ts": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/observation_CFSR_historical_r1i1p1f1_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
        },
        # ERA5, please download:
        # https://aims2.llnl.gov/search?project=CREATE-IP&activeFacets=%7B%22product%22%3A%22reanalysis%22%2C%22experiment%22%3A%22ERA5%22%7D
        "ERA5": {
            "lhf": {
                # the CF standard_name is: surface_upward_latent_heat_flux
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_hfls.xml",
                "varname": "hfls",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "lwr": {
                # the CF standard_name is:
                # surface_downwelling_longwave_flux_in_air and surface_upwelling_longwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rlds.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rlus.xml"],
                "varname": ["rlds", "rlus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "pr": {
                # the CF standard_name is: precipitation_flux
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_pr.xml",
                "varname": "pr",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "shf": {
                # the CF standard_name is: surface_upward_sensible_heat_flux
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_hfss.xml",
                "varname": "hfss",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "sst": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "swr": {
                # the CF standard_name is:
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rsds.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rsus.xml"],
                "varname": ["rsds", "rsus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "taux": {
                # the CF standard_name is: surface_downward_eastward_stress
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_tauu.xml",
                "varname": "tauu",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "tauy": {
                # the CF standard_name is: surface_downward_northward_stress
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_tauv.xml",
                "varname": "tauv",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "thf": {
                # the CF standard_name is:
                # surface_upward_latent_heat_flux, surface_upward_sensible_heat_flux,
                # surface_downwelling_longwave_flux_in_air, surface_upwelling_longwave_flux_in_air,
                # surface_downwelling_shortwave_flux_in_air and surface_upwelling_shortwave_flux_in_air
                # they need to be listed in this order!
                "path + filename": [
                    path_i + "/observation_ERA5_historical_r1i1p1f1_hfls.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_hfss.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rlds.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rlus.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rsds.xml",
                    path_i + "/observation_ERA5_historical_r1i1p1f1_rsus.xml"],
                "varname": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
            "ts": {
                # the CF standard_name is: surface_temperature
                "path + filename": path_i + "/observation_ERA5_historical_r1i1p1f1_ts.xml",
                "varname": "ts",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
        },
        # GPCP-Monthly-3-2, please download:
        # https://aims2.llnl.gov/search?project=obs4MIPs&activeFacets=%7B%22variable%22%3A%22pr%22%2C%22source_id%22%3A%22GPCP-Monthly-3-2%22%7D
        "GPCP-Monthly-3-2": {
            "pr": {
                # the CF standard_name is: precipitation_flux
                "path + filename": path_i + "/observation_GPCP-Monthly-3-2_historical_r1i1p1f1_pr.xml",
                "varname": "pr",
                # the CF standard_name is: cell_area
                "path + filename_area": None,
                "areaname": None,
                # the CF standard_name is: land_area_fraction
                "path + filename_landmask": None,
                "landmaskname": None,
            },
        },
    },
}

# variables have a nickname in the code:
#   lhf: Surface Latent Heat Flux (upward or downward), in CMIP this corresponds to hfls
#   lwr: Net Surface Downward Longwave Radiation, in CMIP this corresponds to rlds - rlus
#   pr: precipitation, in CMIP this corresponds to pr
#   shf: Surface Sensible Heat Flux (upward or downward), in CMIP this corresponds to hfss
#   ssh: Sea Surface Height Above Geoid, in CMIP this corresponds to zos
#   sst: Sea Surface Temperature, in CMIP this corresponds to ts
#   swr: Net Surface Downward Shortwave Radiation, in CMIP this corresponds to rsds - rsus
#   taux: Surface Downward Eastward Wind Stress, in CMIP this corresponds to tauu
#   tauy: Surface Downward Northward Wind Stress, in CMIP this corresponds to tauv
#   thf: Net Surface Downward Heat Flux, in CMIP this corresponds to hfls + hfss + rlds - rlus + rsds - rsus

# Note that the CMIP / observations computation of the fluxes is hard coded in the package
# If your data does not match the way it is defined in the package, you need to adapt the function
# ReferenceObservations in lib/EnsoCollectionsLib.py (for observation)
# CmipVariables in lib/EnsoCollectionsLib.py (for models)
# You can also add/modify observation datasets in ReferenceObservations
# ---------------------------------------------------------------------------------------------------------------------#




# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
# Compute metric collection
# dict_o, _ = ComputeCollection(
#     metric_collection, input_dataset_dictionary, dataset, netcdf=True, netcdf_name=ouput_netcdf_file)
# Note that the epoch used to compute the metric is defined in the metric collection, but you can override it!
# observed_fyear is the first year used for observational dataset
# observed_lyear is the last year used for observational dataset
# so the package will try to read the observations within the intervall [observed_fyear, observed_lyear], a shorter
# period will be read if a shorter period is available. The code will only keep full years (from Jan to Dec)
# the same applies for model data (modeled_fyear, modeled_lyear)
# e.g.:
dict_o, _ = ComputeCollection(
    metric_collection, input_dataset_dictionary, dataset, netcdf=True, netcdf_name=ouput_netcdf_file,
    observed_fyear=1980, observed_lyear=2014, modeled_fyear=1980, modeled_lyear=2014)
# save as json file
with open(ouput_json_file, "w") as outfile:
    json.dump({dataset: dict_o}, outfile, sort_keys=True)
# ---------------------------------------------------------------------------------------------------------------------#
