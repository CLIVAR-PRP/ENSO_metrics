# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Definition of climate variables
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def variables_cmip(variable=None):
    dict_o = {
        "file_name": "<variable>" + "_" + "<table_id>" + "_" + "<source_id>" + "_" + "<experiment_id>" + "_" +
                     "<variant_label>" + "_" + "<grid_label>" + "_*.nc",
        "website": "https://aims2.llnl.gov/search",
        # cell area
        "areacella": {
            "variable": "areacella",
            "variable_computation": "areacella",
            "variable_offset": None,
            "variable_scaling": None,
        },
        "areacello": {
            "variable": "areacello",
            "variable_computation": "areacello",
            "variable_offset": None,
            "variable_scaling": None,
        },
        # land-sea mask
        "landmask": {
            "variable": "sftlf",
            "variable_computation": "1e-2 * sftlf",
            "variable_offset": None,
            "variable_scaling": 1e-2,
        },
        # latent heat flux
        "lhf": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "hfls",
            "variable_computation": "- hfls",
            "variable_offset": None,
            "variable_scaling": -1,
        },
        # longwave radiation
        "lwr": {
            "area": "areacella",
            "mask": "landmask",
            "variable": ["rlds", "rlus"],
            "variable_computation": "rlds - rlus",
            "variable_offset": None,
            "variable_scaling": {"rlds": 1, "rlus": -1},
        },
        # net surface heat flux
        "nhf": {
            "area": "areacella",
            "mask": "landmask",
            "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
            "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
            "variable_offset": None,
            "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
        },
        "thf": {
            "area": "areacella",
            "mask": "landmask",
            "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
            "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
            "variable_offset": None,
            "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
        },
        # rainfall
        "pr": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "pr",
            "variable_computation": "86400 * pr",
            "variable_offset": None,
            "variable_scaling": 86400,
        },
        # sensible heat flux
        "shf": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "hfss",
            "variable_computation": "- hfss",
            "variable_offset": None,
            "variable_scaling": -1,
        },
        # sea level pressure
        "slp": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "psl",
            "variable_computation": "1e-2 * psl",
            "variable_offset": None,
            "variable_scaling": 1e-2,
        },
        # sea surface height
        "ssh": {
            "area": "areacello",
            "mask": None,
            "variable": "zos",
            "variable_computation": "1e2 * zos",
            "variable_offset": None,
            "variable_scaling": 1e2,
        },
        # sea surface temperature
        "sst": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "ts",
            "variable_computation": "ts - 273.15",
            "variable_offset": -273.15,
            "variable_scaling": None,
        },
        # shortwave radiation
        "swr": {
            "area": "areacella",
            "mask": "landmask",
            "variable": ["rsds", "rsus"],
            "variable_computation": "rsds - rsus",
            "variable_offset": None,
            "variable_scaling": {"rsds": 1, "rsus": -1},
        },
        # zonal surface wind stress
        "taux": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "tauu",
            "variable_computation": "1e3 * tauu",
            "variable_offset": None,
            "variable_scaling": 1e3,
        },
        # meridional surface wind stress
        "tauy": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "tauv",
            "cf_units": "Pa",
            "variable_computation": "1e3 * tauv",
            "variable_offset": None,
            "variable_scaling": 1e3,
        },
        # surface temperature
        "ts": {
            "area": "areacella",
            "mask": "landmask",
            "variable": "ts",
            "variable_computation": "ts - 273.15",
            "variable_offset": -273.15,
            "variable_scaling": None,
        },
    }
    if variable in list(dict_o.keys()):
        return dict_o[variable]
    else:
        return dict_o


def variables_observation(dataset=None, variable=None):
    dict_o = {
        "20CRv2c": {
            "file_name": "<variable>" + "_Amon_reanalysis_20CRv2c_185101-201212.nc",
            "website": "https://esgf-node.llnl.gov/search/create-ip/",
            # land-sea mask
            "landmask": {
                "variable": "land",
                "variable_computation": "land",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfls",
                "variable_computation": "- hfls",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rlds", "rlus"],
                "variable_computation": "rlds - rlus",
                "variable_offset": None,
                "variable_scaling": {"rlds": 1, "rlus": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "pr",
                "variable_computation": "86400 * pr",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfss",
                "variable_computation": "- hfss",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "psl",
                "variable_computation": "1e-2 * psl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rsds", "rsus"],
                "variable_computation": "rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"rsds": 1, "rsus": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "tauu",
                "variable_computation": "1e3 * tauu",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "tauv",
                "cf_units": "Pa",
                "variable_computation": "1e3 * tauv",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "20CRv3": {
            "file_name": "<variable>" + "*.mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html",
            # land-sea mask
            "landmask": {
                "variable": "land",
                "variable_computation": "land",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "lhtfl",
                "variable_computation": "- lhtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dlwrf", "ulwrf"],
                "variable_computation": "dlwrf - ulwrf",
                "variable_offset": None,
                "variable_scaling": {"dlwrf": 1, "ulwrf": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "prate",
                "variable_computation": "86400 * prate",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "shtfl",
                "variable_computation": "- shtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "prmsl",
                "variable_computation": "1e-2 * prmsl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dswrf", "uswrf"],
                "variable_computation": "dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"dswrf": 1, "uswrf": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "uflx",
                "variable_computation": "-1e3 * uflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "vflx",
                "cf_units": "Pa",
                "variable_computation": "-1e3 * vflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "AVISO": {
            "website": "https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/files?subdataset=" +
                       "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M-m_202112",
            "file_name": "dt_global_allsat_msla_h_y????_m??.nc",
            # sea surface height
            "ssh": {
                "area": None,
                "mask": None,
                "variable": "adt",
                "variable_computation": "1e2 * adt",
                "variable_offset": None,
                "variable_scaling": 1e2,
            },
        },
        "CFSR": {
            "file_name": "<variable>" + "_" + "<table_id>" + "_reanalysis_CFSR_197901-201912.nc",
            "website": "https://aims2.llnl.gov/search/create-ip/",
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "estimate",
                "variable": "hfls",
                "variable_computation": "- hfls",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "estimate",
                "variable": ["rlds", "rlus"],
                "variable_computation": "rlds - rlus",
                "variable_offset": None,
                "variable_scaling": {"rlds": 1, "rlus": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "estimate",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            "thf": {
                "area": None,
                "mask": "estimate",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "estimate",
                "variable": "pr",
                "variable_computation": "86400 * pr",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "estimate",
                "variable": "hfss",
                "variable_computation": "- hfss",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "estimate",
                "variable": "psl",
                "variable_computation": "1e-2 * psl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface height
            "ssh": {
                "area": None,
                "mask": None,
                "variable": "zos",
                "variable_computation": "1e2 * zos",
                "variable_offset": None,
                "variable_scaling": 1e2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "estimate",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "estimate",
                "variable": ["rsds", "rsus"],
                "variable_computation": "rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"rsds": 1, "rsus": -1},
            },
            # surface wind stress
            "taux": {
                "area": None,
                "mask": "estimate",
                "variable": "tauu",
                "variable_computation": "1e3 * tauu",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "estimate",
                "variable": "tauv",
                "cf_units": "Pa",
                "variable_computation": "1e3 * tauv",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "estimate",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "CMAP": {
            "file_name": "<variable>" + ".mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.cmap.html",
            # land-sea mask
            "landmask": {
                "variable": "lsmask",
                "variable_computation": "lsmask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "precip",
                "variable_computation": "precip",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "COBE1": {
            "file_name": "<variable>" + ".mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.cobe.html",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "COBE2": {
            "variable": "<var_name>" + ".mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.cobe2.html",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "ERA-Interim": {
            "file_name": "<variable>" + "_Amon_reanalysis_ERA-Interim_197901-201908.nc",
            "website": "https://aims2.llnl.gov/search/create-ip/",
            # land-sea mask
            "landmask": {
                "variable": "lsmask",
                "variable_computation": "lsmask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfls",
                "variable_computation": "- hfls",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rlds", "rlus"],
                "variable_computation": "rlds - rlus",
                "variable_offset": None,
                "variable_scaling": {"rlds": 1, "rlus": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "pr",
                "variable_computation": "86400 * pr",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfss",
                "variable_computation": "- hfss",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "psl",
                "variable_computation": "1e-2 * psl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rsds", "rsus"],
                "variable_computation": "rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"rsds": 1, "rsus": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "tauu",
                "variable_computation": "1e3 * tauu",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "tauv",
                "cf_units": "Pa",
                "variable_computation": "1e3 * tauv",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "ERA5": {
            "file_name": "<variable>" + "_Amon_reanalysis_ERA5_197901-202012.nc",
            "website": "https://aims2.llnl.gov/search/create-ip/",
            # land-sea mask
            "landmask": {
                "variable": "lsm",
                "variable_computation": "lsm",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfls",
                "variable_computation": "- hfls",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rlds", "rlus"],
                "variable_computation": "rlds - rlus",
                "variable_offset": None,
                "variable_scaling": {"rlds": 1, "rlus": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
                "variable_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"hfls": -1, "hfss": -1, "rlds": 1, "rlus": -1, "rsds": 1, "rsus": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "pr",
                "variable_computation": "86400 * pr",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "hfss",
                "variable_computation": "- hfss",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "psl",
                "variable_computation": "1e-2 * psl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["rsds", "rsus"],
                "variable_computation": "rsds - rsus",
                "variable_offset": None,
                "variable_scaling": {"rsds": 1, "rsus": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "tauu",
                "variable_computation": "1e3 * tauu",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "tauv",
                "cf_units": "Pa",
                "variable_computation": "1e3 * tauv",
                "variable_offset": None,
                "variable_scaling": 1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "ts",
                "variable_computation": "ts - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "ERSSTv3b": {
            "file_name": "<variable>" + ".mnmean.v3.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.noaa.ersst.v3.html",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "ERSSTv4": {
            "file_name": "<variable>" + ".mnmean.v4.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.noaa.ersst.v4.html",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "ERSSTv5": {
            "file_name": "<variable>" + ".mnmean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "GPCPv2.3": {
            "file_name": "<variable>" + ".mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.gpcp.html",
            # land-sea mask
            "landmask": {
                "variable": "lsmask",
                "variable_computation": "lsmask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "precip",
                "variable_computation": "precip",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "HadISST": {
            "file_name": "HadISST_" + "<variable>" + ".nc",
            "website": "https://www.metoffice.gov.uk/hadobs/hadisst/",
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
        },
        "NCEP1": {
            "file_name": "<variable>" + ".sfc.mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.html",
            # land-sea mask
            "landmask": {
                "variable": "lsmask",
                "variable_computation": "lsmask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "lhtfl",
                "variable_computation": "- lhtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dlwrf", "ulwrf"],
                "variable_computation": "dlwrf - ulwrf",
                "variable_offset": None,
                "variable_scaling": {"dlwrf": 1, "ulwrf": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "prate",
                "variable_computation": "86400 * prate",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "shtfl",
                "variable_computation": "- shtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "prmsl",
                "variable_computation": "1e-2 * prmsl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dswrf", "uswrf"],
                "variable_computation": "dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"dswrf": 1, "uswrf": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "uflx",
                "variable_computation": "-1e3 * uflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "vflx",
                "cf_units": "Pa",
                "variable_computation": "-1e3 * vflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "NCEP2": {
            "file_name": "<variable>" + ".sfc.mon.mean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html",
            # land-sea mask
            "landmask": {
                "variable": "lsmask",
                "variable_computation": "lsmask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": "landmask",
                "variable": "lhtfl",
                "variable_computation": "- lhtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dlwrf", "ulwrf"],
                "variable_computation": "dlwrf - ulwrf",
                "variable_offset": None,
                "variable_scaling": {"dlwrf": 1, "ulwrf": -1},
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            "thf": {
                "area": None,
                "mask": "landmask",
                "variable": ["lhtfl", "shtfl", "dlwrf", "ulwrf", "dswrf", "uswrf"],
                "variable_computation": "- lhtfl - shtfl + dlwrf - ulwrf + dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"lhtfl": -1, "shtfl": -1, "dlwrf": 1, "ulwrf": -1, "dswrf": 1, "uswrf": -1},
            },
            # rainfall
            "pr": {
                "area": None,
                "mask": "landmask",
                "variable": "prate",
                "variable_computation": "86400 * prate",
                "variable_offset": None,
                "variable_scaling": 86400,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": "landmask",
                "variable": "shtfl",
                "variable_computation": "- shtfl",
                "variable_offset": None,
                "variable_scaling": -1,
            },
            # sea level pressure
            "slp": {
                "area": None,
                "mask": "landmask",
                "variable": "prmsl",
                "variable_computation": "1e-2 * prmsl",
                "variable_offset": None,
                "variable_scaling": 1e-2,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": "landmask",
                "variable": ["dswrf", "uswrf"],
                "variable_computation": "dswrf - uswrf",
                "variable_offset": None,
                "variable_scaling": {"dswrf": 1, "uswrf": -1},
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": "landmask",
                "variable": "uflx",
                "variable_computation": "-1e3 * uflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": "landmask",
                "variable": "vflx",
                "cf_units": "Pa",
                "variable_computation": "-1e3 * vflx",
                "variable_offset": None,
                "variable_scaling": -1e3,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "OISSTv2": {
            "file_name": "<variable>" + ".mnmean.nc",
            "website": "https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html",
            # land-sea mask
            "landmask": {
                "variable": "mask",
                "variable_computation": "mask",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
            # surface temperature
            "ts": {
                "area": None,
                "mask": "landmask",
                "variable": "skt",
                "variable_computation": "skt - 273.15",
                "variable_offset": -273.15,
                "variable_scaling": None,
            },
        },
        "Tropflux": {
            "file_name": "<variable>" + "_tropflux_1m_*.nc",
            "website": "https://incois.gov.in/tropflux/tf_products.jsp",
            # latent heat flux
            "lhf": {
                "area": None,
                "mask": None,
                "variable": "lhf",
                "variable_computation": "lhf",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # longwave radiation
            "lwr": {
                "area": None,
                "mask": None,
                "variable": "lwr",
                "variable_computation": "lwr",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # net heat flux
            "nhf": {
                "area": None,
                "mask": None,
                "variable": "netflux",
                "variable_computation": "netflux",
                "variable_offset": None,
                "variable_scaling": None,
            },
            "thf": {
                "area": None,
                "mask": None,
                "variable": "netflux",
                "variable_computation": "netflux",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # sensible heat flux
            "shf": {
                "area": None,
                "mask": None,
                "variable": "shf",
                "variable_computation": "shf",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # sea surface temperature
            "sst": {
                "area": None,
                "mask": None,
                "variable": "sst",
                "variable_computation": "sst",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # shortwave radiation
            "swr": {
                "area": None,
                "mask": None,
                "variable": "swr",
                "variable_computation": "swr",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # zonal surface wind stress
            "taux": {
                "area": None,
                "mask": None,
                "variable": "taux",
                "variable_computation": "taux",
                "variable_offset": None,
                "variable_scaling": None,
            },
            # meridional surface wind stress
            "tauy": {
                "area": None,
                "mask": None,
                "variable": "tauy",
                "variable_computation": "tauy",
                "variable_offset": None,
                "variable_scaling": None,
            },
            
        },
    }
    # return
    if dataset in list(dict_o.keys()) and variable in list(dict_o[dataset].keys()):
        return dict_o[dataset][variable]
    elif dataset in list(dict_o.keys()) and variable not in list(dict_o[dataset].keys()):
        return dict_o[dataset]
    else:
        return dict_o


def variables_param(variable=None):
    # cf out_name / units for Amon
    # https://github.com/PCMDI/cmip6-cmor-tables/blob/main/Tables/CMIP6_Amon.json
    # cf out_name / units for Omon
    # https://github.com/PCMDI/cmip6-cmor-tables/blob/main/Tables/CMIP6_Omon.json
    dict_o = {
        # cell area
        "areacella": {
            "cf_computation": "areacella",
            "cf_out_name": "areacella",
            "cf_units": "m2",
            "long_name": "Grid-Cell Area for Atmospheric Variables",
            "short_name": "area",
            "units": "m**2",
            "variable_type": "area",
        },
        "areacello": {
            "cf_computation": "areacello",
            "cf_out_name": "areacello",
            "cf_units": "m2",
            "long_name": "Grid-Cell Area for Ocean Variables",
            "short_name": "area",
            "units": "m**2",
            "variable_type": "area",
        },
        # land-sea mask
        "landmask": {
            "cf_computation": "1e-2 * sftlf",
            "cf_out_name": "sftlf",
            "cf_units": "%",
            "long_name": "Grid-Cell Land Fraction for Atmospheric Variables",
            "short_name": "lsmask",
            "units": "1",
            "variable_type": "fraction",
        },
        # latent heat flux
        "lhf": {
            "cf_computation": "- hfls",
            "cf_out_name": "hfls",
            "cf_units": "W m-2",
            "long_name": "Surface Downward Latent Heat Flux (positive = heat into the ocean)",
            "short_name": "LHF",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        # longwave radiation
        "lwr": {
            "cf_computation": "rlds - rlus",
            "cf_out_name": ["rlds", "rlus"],
            "cf_units": "W m-2",
            "long_name": "Net Surface Downward Longwave Radiation (positive = heat into the ocean)",
            "short_name": "LWR",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        # net surface heat flux
        "nhf": {
            "cf_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
            "cf_out_name": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
            "cf_units": "W m-2",
            "long_name": "Net Surface Downward Heat Flux (positive = heat into the ocean)",
            "short_name": "NHF",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        "thf": {
            "cf_computation": "- hfls - hfss + rlds - rlus + rsds - rsus",
            "cf_out_name": ["hfls", "hfss", "rlds", "rlus", "rsds", "rsus"],
            "cf_units": "W m-2",
            "long_name": "Net Surface Downward Heat Flux (positive = heat into the ocean)",
            "short_name": "NHF",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        # rainfall
        "pr": {
            "cf_computation": "86400 * pr",
            "cf_out_name": "pr",
            "cf_units": "kg m-2 s-1",
            "long_name": "precipitation",
            "short_name": "PR",
            "units": "mm.day**-1",
            "variable_type": "precipitation",
        },
        # sea level pressure
        "slp": {
            "cf_computation": "1e-2 * psl",
            "cf_out_name": "psl",
            "cf_units": "Pa",
            "long_name": "Sea Level Pressure",
            "short_name": "SLP",
            "units": "hPa",
            "variable_type": "pressure",
        },
        # sensible heat flux
        "shf": {
            "cf_computation": "- hfss",
            "cf_out_name": "hfss",
            "cf_units": "W m-2",
            "long_name": "Surface Downward Sensible Heat Flux (positive = heat into the ocean)",
            "short_name": "SHF",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        # sea surface height
        "ssh": {
            "cf_computation": "1e2 * zos",
            "cf_out_name": "zos",
            "cf_units": "m",
            "long_name": "Sea Surface Height (Above Geoid or Anomaly Above Sea Level)",
            "short_name": "SSH",
            "units": "cm",
            "variable_type": "sea surface height",
        },
        # sea surface temperature
        "sst": {
            "cf_computation": "ts - 273.15",
            "cf_out_name": "ts",
            "cf_units": "K",
            "long_name": "Sea Surface Temperature",
            "short_name": "SST",
            "units": "degC",
            "variable_type": "temperature",
        },
        # shortwave radiation
        "swr": {
            "cf_computation": "rsds - rsus",
            "cf_out_name": ["rsds", "rsus"],
            "cf_units": "W m-2",
            "long_name": "Net Surface Downward Shortwave Radiation (positive = heat into the ocean)",
            "short_name": "SWR",
            "units": "W.m**-2",
            "variable_type": "heat flux",
        },
        # zonal surface wind stress
        "taux": {
            "cf_computation": "1e3 * tauu",
            "cf_out_name": "tauu",
            "cf_units": "Pa",
            "long_name": "Surface Downward Eastward Wind Stress",
            "short_name": "Taux",
            "units": "10**-3 Pa",
            "variable_type": "wind stress",
        },
        # meridional surface wind stress
        "tauy": {
            "cf_computation": "1e3 * tauv",
            "cf_out_name": "tauv",
            "cf_units": "Pa",
            "long_name": "Surface Downward Northward Wind Stress",
            "short_name": "Tauy",
            "units": "10**-3 Pa",
            "variable_type": "wind stress",
        },
        # surface temperature
        "ts": {
            "cf_computation": "ts - 273.15",
            "cf_out_name": "ts",
            "cf_units": "K",
            "long_name": "Surface Temperature (skin for open ocean)",
            "short_name": "TS",
            "units": "degC",
            "variable_type": "temperature",
        },
    }
    if variable in list(dict_o.keys()):
        return dict_o[variable]
    else:
        return deepcopy(dict_o)
# ---------------------------------------------------------------------------------------------------------------------#
