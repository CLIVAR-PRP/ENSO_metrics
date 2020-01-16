# -*- coding:UTF-8 -*-
#
# Define ENSO metrics collections as a function of science question/realm
#
# Draft version
#


# Define metrics collections
def defCollection(mc=True):
    # Name, list of metrics
    metrics_collection = {
        'ENSO_perf': {
            'long_name': 'Metrics Collection for ENSO performance',
            'metrics_list': {
                'EnsoAmpl': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST', 'Tropflux']},
                    'metric_computation': 'abs_relative_difference',
                    'regridding': {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                   'newgrid_name': 'generic_1x1deg'},
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 204,
                'normalization': False,
                'project_interpreter': 'CMIP',
                'observed_period': ('1850-01-01 00:00:00', '2018-12-31 23:59:60.0'),
                'modeled_period': ('1850-01-01 00:00:00', '2015-12-31 23:59:60.0'),
            },
            'plot_order': ['EnsoAmpl'],
            'description': 'Describe which science question this collection is about',
        },
    }
    if mc is True:
        return metrics_collection
    else:
        return metrics_collection[mc]


# List of reference observations for each variables
def ReferenceObservations(dataset=True):
    dict_ref_obs = {
        'AVISO': {
            'website': 'see https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global.html',
            'file_name': 'dt_global_allsat_msla_h_y????_m??.nc',
            'variable_name_in_file': {'ssh': {'var_name': 'sla'}, },
        },
        'ERA-Interim': {
            'website': 'see https://esgf.nccs.nasa.gov/search/create-ip/',
            'file_name': '<var_name>' + '_Amon_reanalysis_IFS-Cy31r2_*.nc',
            'variable_name_in_file': {
                'landmask': {'var_name': 'lsmask'},
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
            'website': 'see https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/',
            'file_name': 'ersst.v5.' + '<YYYYMM>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'GPCPv2.3': {
            'website': 'see https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=GPCP+Version+2.3+' +
                       'Combined+Precipitation+Dataset&group=0&submit=Search',
            'file_name': 'precip.mon.mean.nc',
            'variable_name_in_file': {
                'landmask': {'var_name': 'lsmask'},
                'pr': {'var_name': 'precip'},
            },
        },
        'HadISST': {
            'website': 'see https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html',
            'file_name': 'HadISST_' + '<var_name>' + '.nc',
            'variable_name_in_file': {
                'sst': {'var_name': 'sst'},
            },
        },
        'SODA3.4.2': {
            'website': 'see https://www.atmos.umd.edu/~ocean/index_files/soda3.4.2_mn_download_b.htm',
            'file_name': 'soda3.4.2_mn_ocean_reg_????.nc',
            'variable_name_in_file': {
                'so': {'var_name': 'salt'},
                'ssh': {'var_name': 'ssh'},
                'taux': {'var_name': 'taux'},
                'thetao': {'var_name': 'temp'},
                'uo': {'var_name': 'u'},
                'vo': {'var_name': 'v'},
            },
        },
        'Tropflux': {
            'website': 'see https://incois.gov.in/tropflux/tf_products.jsp',
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
    if dataset is True:
        return dict_ref_obs
    else:
        return dict_ref_obs[dataset]


def ReferenceRegions(region=True):
    dict_reference_regions = {
        'global': {'long_name': 'Global 60S-60N', 'latitude': (-60., 60.), 'longitude': (0., 360.)},
        'equatorial_pacific': {
            'long_name': 'Equatorial Pacific (EP)', 'latitude': (-5., 5.), 'longitude': (150., 270.),
        },
        'equatorial_pacific_LatExt': {
            'long_name': 'Equatorial Pacific (EP)', 'latitude': (-15., 15.), 'longitude': (150., 270.),
        },
        'equatorial_pacific_LatExt2': {
            'long_name': 'Equatorial Pacific extended in latitude', 'latitude': (-15., 15.), 'longitude': (120., 285.),
        },
        'nino3': {'long_name': 'Niño 3', 'latitude': (-5., 5.), 'longitude': (210., 270.)},
        'nino3_LatExt': {
            'long_name': 'Niño 3 extended in latitude', 'latitude': (-15., 15.), 'longitude': (210., 270.)
        },
        'nino3.4': {'long_name': 'Niño 3.4', 'latitude': (-5., 5.), 'longitude': (190., 240.)},
        'nino4': {'long_name': 'Niño 4', 'latitude': (-5., 5.), 'longitude': (160., 210.)},
        # YYP regions
        'africaSE': {'long_name': 'South and East Africa', 'latitude': (-40, 15.), 'longitude': (0., 55.),
                     'maskland': False, 'maskocean': True},
        'americaN': {'long_name': 'North America', 'latitude': (10., 60.), 'longitude': (235., 300.),
                     'maskland': False, 'maskocean': True},
        'americaS': {'long_name': 'South America', 'latitude': (-60., 15.), 'longitude': (275., 330.),
                     'maskland': False, 'maskocean': True},
        'asiaS': {'long_name': 'South Asia', 'latitude': (-10., 30.), 'longitude': (65., 130.),
                  'maskland': False, 'maskocean': True},
        'oceania': {'long_name': 'oceania', 'latitude': (-50., 0.), 'longitude': (110., 180.), 'maskland': False,
                    'maskocean': True},
    }
    if region is True:
        return dict_reference_regions
    else:
        return dict_reference_regions[region]


def CmipVariables():
    dict_cmip_variables = {
        'reference': 'http://cfconventions.org/Data/cf-standard-names/46/build/cf-standard-name-table.html',
        'variable_name_in_file': {
            # line keys:
            # '<internal_metrics_variable_name>':{'var_name':'<var_name_in_file>','cf_name':<as per ref above>,
            #                                     'cf_unit':'<unit_in_file>'}
            # areacell
            'areacell': {'var_name': 'areacella', 'cf_name': 'cell_area', 'cf_units': 'm2'},
            # landmask
            'landmask': {'var_name': 'sftlf', 'cf_name': 'cell_area', 'cf_units': '1'},
            # latent heat flux (on ocean grid or ocean points only)
            'lhf': {'var_name': 'hfls', 'cf_name': 'surface_upward_latent_heat_flux', 'cf_units': 'W m-2'},
            # longwave radiation computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # lwr = rlds - rlus
            # sometimes lwr is included in the datasets in a variable called 'rls'
            'lwr': {
                'var_name': ['rlds', 'rlus'],
                'cf_name': ['surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']},
            # Rainfall Flux
            'pr': {'var_name': 'pr', 'cf_name': 'rainfall_flux', 'cf_units': 'kg m-2 s-1'},
            # Sea Level Pressure
            'slp': {'var_name': 'psl', 'cf_name': 'air_pressure_at_mean_sea_level', 'cf_units': 'Pa'},
            # sensible heat flux (on ocean grid or ocean points only)
            'shf': {'var_name': 'hfss', 'cf_name': 'surface_upward_sensible_heat_flux', 'cf_units': 'W m-2'},
            # sea surface height
            'ssh': {'var_name': 'zos', 'cf_name': 'sea_surface_height_above_geoid', 'cf_units': 'm'},
            # sea surface temperature
            'sst': {'var_name': 'ts', 'cf_name': 'sea_surface_temperature', 'cf_units': 'K'},
            # shortwave radiation computed from these variables IN THAT ORDER
            # swr = rsds - rsus
            # sometimes swr is included in the datasets in a variable called 'rss'
            'swr': {
                'var_name': ['rsds', 'rsus'],
                'cf_name': ['surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']
            },
            # zonal surface wind stress
            'taux': {'var_name': 'tauu', 'cf_name': 'surface_downward_eastward_stress', 'cf_units': 'Pa'},
            # total heat flux computed from these variables IN THAT ORDER
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
