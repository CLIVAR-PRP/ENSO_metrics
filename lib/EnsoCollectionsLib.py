# -*- coding:UTF-8 -*-
#
# Define ENSO metrics collections as a function of science question/realm
#
# Draft version
#


# Define metrics collections
def defCollection(MC=True):
    # Name, list of metrics
    metrics_collection = {
        'MC1': {
            'long_name': 'Metrics Collection 1',
            'common_collection_parameters': {
                'frequency': 'monthly',
                'min_number_of_years': 30,
                'project_interpreter': 'CMIP',
            },
            'metrics_list': {
                'EnsoRmse': {
                    'variables': ['sst'],
                    'obs_name': {'sst': ['ERSSTv5', 'HadISST']},
                    'preprocessing': {
                        0: {
                            'selection_period_and_region': {
                                'period': {
                                    'model': ('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                                    'observations': ('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                                },
                                'regions': {
                                    'sst': 'tropical_pacific',
                                },
                                'units': {
                                    'sst': 'C',
                                },
                            },
                        },
                        1: {
                            'anomalies': False
                        },
                        2: {
                            'detrending': {'method': 'linear'},
                        },
                        3: {
                            'normalization': False,
                        },
                        4: {
                            'smoothing': False,
                        },
                        5: {
                            'averaging': {'average_dimension': ['time']},
                        },
                        6: {
                            'regridding': {
                                'newgrid': 'HadISST', 'newgrid_name': 'generic 1x1deg',
                                'mask_processing': 'model + observations', 'regridTool': 'esmf',
                                'regridMethod': 'linear',
                            },
                        },
                    },
                },
            },
            'description': 'Describe which science question this collection is about',
        },
    }
    if MC:
        return metrics_collection[MC]
    else:
        return metrics_collection


# List of reference observations for each variables
def ReferenceObservations(DATASET=True):
    dict_ref_obs = {
        'ERSSTv5': {
            'source': 'see https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/',
            'file_name': 'ersst.v5.' + '<YYYYMM>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'HadISST': {
            'source': 'see https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html',
            'file_name': 'HadISST_' + '<var_name>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'OISST': {
            'source': 'see https://www.earthsystemcog.org/search/obs4mips/?template=obs4mips&limit=200',
            'file_name': '<var_name>' + '_OISST_L4_AVHRR-only-v2_*-*.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
    }
    if DATASET:
        return dict_ref_obs[DATASET]
    else:
        return dict_ref_obs


def ReferenceRegions(AR=True):
    dict_reference_regions = {
        'tropical_pacific': {
            'long_name': 'Tropical Pacific (TP)', 'latitude': (-30., 30.), 'longitude': (120., 280.),
        },
        'equatorial_pacific': {
            'long_name': 'Equatorial Pacific (EP)', 'latitude': (-5., 5.), 'longitude': (120., 280.),
        },
        'eastern_equatorial_pacific': {
            'long_name': 'Eastern Equatorial Pacific (EEP)', 'latitude': (-5., 5.), 'longitude': (120., 205.),
        },
        'western_equatorial_pacific': {
            'long_name': 'Western Equatorial Pacific (WEP)', 'latitude': (-5., 5.), 'longitude': (205., 280.),
        },
        'nino1+2': {'long_name': 'Niño 1+2', 'latitude': (-10., 0.), 'longitude': (270., 280.)},
        'nino3': {'long_name': 'Niño 3', 'latitude': (-5., 5.), 'longitude': (210., 270.)},
        'nino3.4': {'long_name': 'Niño 3.4', 'latitude': (-5., 5.), 'longitude': (190., 240.)},
        'nino4': {'long_name': 'Niño 4', 'latitude': (-5., 5.), 'longitude': (160., 210.)},
        'nino3_reduced': {'long_name': 'Niño 3', 'latitude': (-5., 5.), 'longitude': (220., 250.)},
    }
    if AR:
        return dict_reference_regions[AR]
    else:
        return dict_reference_regions


def CmipVariables():
    dict_cmip_variables = {
        'reference': 'http://cfconventions.org/Data/cf-standard-names/46/build/cf-standard-name-table.html',
        'variable_name_in_file': {
            # line keys:
            # '<internal_metrics_variable_name>': {'var_name': '<var_name_in_file>', 'cf_name': <as per ref above>,
            # 'cf_unit': '<unit_in_file>'}
            #
            # latent heat flux (on ocean grid or ocean points only)
            'lhf': {'var_name': 'hfls', 'cf_name': 'surface_upward_latent_heat_flux', 'cf_units': 'W m-2'},
            # longwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # lwr = rlds - rlus
            # sometimes lwr is included in the datasets in a variable called 'rls'
            'lwr': {
                'var_name': ['rlds', 'rlus'],
                'cf_name': ['surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air',],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']},
            # Rainfall Flux (on ocean grid or ocean points only)
            'pr': {'var_name': 'pr', 'cf_name': 'rainfall_flux', 'cf_units': 'kg m-2 s-1'},
            # sensible heat flux (on ocean grid or ocean points only)
            'shf': {'var_name': 'hfss', 'cf_name': 'surface_upward_sensible_heat_flux', 'cf_units': 'W m-2'},
            # sea surface temperature (on ocean grid or ocean points only)
            'sst': {'var_name': 'tos', 'cf_name': 'sea_surface_temperature', 'cf_units': 'K'},
            # shortwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # swr = rsds - rsus
            # sometimes swr is included in the datasets in a variable called 'rss'
            'swr': {
                'var_name': ['rsds', 'rsus'],
                'cf_name': ['surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']
            },
            # zonal surface wind stress (on ocean grid or ocean points only)
            'taux': {'var_name': 'tauuo', 'cf_name': 'surface_downward_eastward_stress', 'cf_units': 'Pa'},
            # total heat flux computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # tfh = hfls + hfss + rlds - rlus + rsds - rsus
            # sometimes rls = rlds - rlus and rss = rsds - rsus
            # sometimes thf is included in the datasets in a variable called 'hfds', 'netflux', 'thflx',...
            'thf': {
                'var_name': ['hfls', 'hfss', 'rlds', 'rlus', 'rsds', 'rsus'],
                'cf_name': ['surface_upward_latent_heat_flux', 'surface_upward_sensible_heat_flux',
                            'surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air',
                            'surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'plus', 'plus', 'minus', 'plus', 'minus']
            },
        },
    }
    return dict_cmip_variables
