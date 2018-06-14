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
            'metrics_list': {
                'EnsoAmpl': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux']},
                    'metric_computation': 'difference', # i.e., (obs-model)/model
                    # the "science panel" will have to define the "metric", it could be the difference (obs-model), the
                    # ratio (model/model), the relative difference ([obs-model]/model),...

                    # the portrait plot will not be able to show the "real" metric value because of the normalization
                    # (to share a common colorbar) but we could dive down on a plot (like Bellenger et al. 2013, Fig. 1)
                    # showing the metric values and / or the diagnostic values for observations and model
                },
                'EnsoSeasonality': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux']},
                    'metric_computation': 'difference',
                },
                'BiasSstRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'tropical_pacific'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'EnsoAlphaLhf': {
                    'variables': ['sst','lhf'],
                    'regions': {'sst': 'nino3', 'lhf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'lwr': ['Tropflux']},
                    'metric_computation': 'difference',
                },
                'EnsoAlphaLwr': {
                    'variables': ['sst', 'lwr'],
                    'regions': {'sst': 'nino3', 'lwr': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'lwr': ['Tropflux']},
                    'metric_computation': 'difference',
                },
                'EnsoAlphaShf': {
                    'variables': ['sst', 'shf'],
                    'regions': {'sst': 'nino3', 'shf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'shf': ['Tropflux']},
                    'metric_computation': 'difference',
                },
                'EnsoAlphaSwr': {
                    'variables': ['sst', 'swr'],
                    'regions': {'sst': 'nino3', 'swr': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'swr': ['Tropflux']},
                    'metric_computation': 'difference',
                },
                'EnsoAlphaThf': {
                    'variables': ['sst', 'thf'],
                    'regions': {'sst': 'nino3', 'thf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'thf': ['Tropflux']},
                    'metric_computation': 'difference',
                },
                'EnsoMu': {
                    'variables': ['sst', 'taux'],
                    'regions': {'sst': 'nino3', 'taux': 'nino4'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'taux': ['Tropflux']},
                    'metric_computation': 'difference',
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 324,
                'project_interpreter': 'CMIP',
#                'observed_period': ('1979-01-01 00:00:00', '2016-12-31 23:59:60.0'),
#                'modeled_period': ('1979-01-01 00:00:00', '2005-12-31 23:59:60.0'),
            },
            'description': 'Describe which science question this collection is about',
        },
        'ENSO_perf': {
            'long_name': 'Metrics Collection for ENSO performance',
            'metrics_list': {
                'BiasSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'BiasPrLonRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'equatorial_pacific'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'BiasTauxLonRmse': {
                    'variables': ['taux'],
                    'regions': {'taux': 'equatorial_pacific'},
                    'obs_name': {'taux': ['ERA-Interim', 'Tropflux']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'BiasSstLatRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.3_LatExt'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'BiasPrLatRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'nino3.3_LatExt'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'BiasTauxLatRmse': {
                    'variables': ['taux'],
                    'regions': {'taux': 'equatorial_pacific_LatExt'},
                    'obs_name': {'taux': ['ERA-Interim', 'Tropflux']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'SeasonalSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'SeasonalSstLatRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.3_LatExt'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'SeasonalPrLatRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'nino3.3_LatExt'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'EnsoAmpl': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'metric_computation': 'difference',
                },
                'NinoSstTsRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['ERA-Interim']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                },
                'NinaSstTsRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                },
                'NinoSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': 0.75},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
                'NinaSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3', 'season_ev': 'DEC', 'threshold': -0.75},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear'},
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 204,
                'normalization': False,
                'project_interpreter': 'CMIP',
                'observed_period': ('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                'modeled_period': ('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
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
        'CFSR': {
            'source': 'see https://esgf.nccs.nasa.gov/search/create-ip/', # do not use 'source' keyword
            'file_name': '<var_name>' + '_Omon_reanalysis_CFSR_*.nc',
            'variable_name_in_file': {
                'ssh': {'var_name': 'zos'},
                'so': {'var_name': 'so'},
                'thetao': {'var_name': 'thetao'},
                'thf': {'var_name': 'hfds'}, # I'm not sure yet if it is the total heat flux
                'uo': {'var_name': 'uo'},
                'vo': {'var_name': 'vo'},
            },
        },
        'ERA-Interim': {
            'source': 'see https://esgf.nccs.nasa.gov/search/create-ip/', # do not use 'source' keyword
            'file_name': '<var_name>' + '_Amon_reanalysis_IFS-Cy31r2_*.nc',
            'variable_name_in_file': {
                'lhf': {'var_name': 'hfls'},
                # longwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
                # lwr = rlds - rlus
                # sometimes lwr is included in the datasets in a variable called 'rls'
                'lwr': {'var_name': ['rlds', 'rlus'], 'algebric_calculation': ['plus', 'minus']},
                'pr': {'var_name': 'pr'},
                'slp': {'var_name': 'psl'},
                'shf': {'var_name': 'hfss'},
                'sst': {'var_name': 'ts'},
                # shortwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
                # swr = rsds - rsus
                # sometimes swr is included in the datasets in a variable called 'rss'
                'swr': {'var_name': ['rsds', 'rsus'], 'algebric_calculation': ['plus', 'minus']},
                'taux': {'var_name': 'tauu'},
                'tauy': {'var_name': 'tauv'},
                # total heat flux computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
                # tfh = hfls + hfss + rlds - rlus + rsds - rsus
                # sometimes rls = rlds - rlus and rss = rsds - rsus
                # sometimes thf is included in the datasets in a variable called 'hfds', 'netflux', 'thflx',...
                'thf': {
                    'var_name': ['hfls', 'hfss', 'rlds', 'rlus', 'rsds', 'rsus'],
                    'algebric_calculation': ['plus', 'plus', 'plus', 'minus', 'plus', 'minus'],
                },
                'uas': {'var_name': 'uas'},
                'vas': {'var_name': 'vas'},
            },
        },
        'ERSSTv5': {
            'source': 'see https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/',
            'file_name': 'ersst.v5.' + '<YYYYMM>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'GODAS': {
            'source': 'see https://esgf.nccs.nasa.gov/search/create-ip/',
            'file_name': '<var_name>' + '_Omon_ORAreanalysis_GODAS_*.nc',
            'variable_name_in_file': {
                'so': {'var_name': 'so'},
                'thetao': {'var_name': 'thetao'},
                'uo': {'var_name': 'uo'},
                'vo': {'var_name': 'vo'},
                # many variables are missing from https://www.esrl.noaa.gov/psd/data/gridded/data.godas.html
                #   taux = uflx, tauy = vflx, thf = thflx
                #   and
                #   salt flux = sltfl, Sea Surface Height Relative to Geoid = sshg,
                #   Geometric Depth Below Sea Surface = dbss_obil, Geometric Depth Below Sea Surface = dbss_obml,
                #   Geometric vertical velocity (dz/dt) = dzdt
            },
        },
        'GPCPv2.3': {
            'source': 'see https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=GPCP+Version+2.3+'
                      + 'Combined+Precipitation+Dataset&group=0&submit=Search',
            'file_name': 'precip.mon.mean.nc',
            'variable_name_in_file': {
                'pr': {'var_name': 'precip'},
            },
        },
        'HadISST': {
            'source': 'see https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html',
            'file_name': 'HadISST_' + '<var_name>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'OAFlux': {
            'source': 'see ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/monthly/turbulence/',
            'file_name': '<???>' + '_oaflux_*.nc',
            'variable_name_in_file': {
                'lhf': {'var_name': 'lhtfl'},
                'shf': {'var_name': 'shtfl'},
                'sst': {'var_name': 'tmpsf'},
            },
        },
        'OISST': {
            'source': 'see https://www.earthsystemcog.org/search/obs4mips/?template=obs4mips&limit=200',
            'file_name': '<var_name>' + '_OISST_L4_AVHRR-only-v2_*-*.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'ORAS4': {
            'source': 'see https://esgf.nccs.nasa.gov/search/create-ip/',
            'file_name': '<var_name>' + '_Omon_ORAreanalysis_ORAS4_*.nc',
            'variable_name_in_file': {
                'so': {'var_name': 'so'},
                'thetao': {'var_name': 'thetao'},
                'uo': {'var_name': 'uo'},
                'vo': {'var_name': 'vo'},
            },
        },
        'Tropflux': {
            'source': 'see http://www.incois.gov.in/tropflux_datasets/data/monthly/',
            'file_name': '<var_name>' + '_tropflux_1m_*.nc',
            'variable_name_in_file': {
                'lhf': {'var_name': 'lhf'},
                'lwr': {'var_name': 'lwr'},
                'shf': {'var_name': 'shf'},
                'sst': {'var_name': 'sst'},
                'swr': {'var_name': 'swr'},
                'taux': {'var_name': 'taux'},
                'thf': {'var_name': 'netflux'},
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
        'equatorial_pacific_LatExt': {
            'long_name': 'Equatorial Pacific extended in latitude', 'latitude': (-15., 15.), 'longitude': (150., 270.),
        },
        'eastern_equatorial_pacific': {
            'long_name': 'Western Equatorial Pacific (WEP)', 'latitude': (-5., 5.), 'longitude': (205., 280.),
        },
        'western_equatorial_pacific': {
            'long_name': 'Eastern Equatorial Pacific (EEP)', 'latitude': (-5., 5.), 'longitude': (120., 205.),
        },
        'nino1+2': {'long_name': 'Niño 1+2', 'latitude': (-10., 0.), 'longitude': (270., 280.)},
        'nino3': {'long_name': 'Niño 3', 'latitude': (-5., 5.), 'longitude': (210., 270.)},
        'nino3.3_LatExt': {'long_name': 'small Niño 3 extended in latitude', 'latitude': (-15., 15.),
                           'longitude': (220., 250.)},
        'nino3.4': {'long_name': 'Niño 3.4', 'latitude': (-5., 5.), 'longitude': (190., 240.)},
        'nino4': {'long_name': 'Niño 4', 'latitude': (-5., 5.), 'longitude': (160., 210.)},
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
            # '<internal_metrics_variable_name>':{'var_name':'<var_name_in_file>','cf_name':<as per ref above>, 'cf_unit':'<unit_in_file>'}

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
            'sst': {'var_name': 'ts', 'cf_name': 'sea_surface_temperature', 'cf_units': 'K'},
            # shortwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # swr = rsds - rsus
            # sometimes swr is included in the datasets in a variable called 'rss'
            'swr': {
                'var_name': ['rsds', 'rsus'],
                'cf_name': ['surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']
            },
            # zonal surface wind stress (on ocean grid or ocean points only)
            'taux': {'var_name': 'tauu', 'cf_name': 'surface_downward_eastward_stress', 'cf_units': 'Pa'},
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

