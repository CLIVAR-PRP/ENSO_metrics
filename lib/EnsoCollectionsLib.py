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
                    'obs_name': {'sst': ['HadISST1.1','OISST']},
                },
                'EnsoSeasonality': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3'},
                    'obs_name': {'sst': ['HadISST1.1', 'OISST']},
                },
                'EnsoRMSE': {
                    'variables': ['sst'],
                    'regions': {'sst': 'tropical_pacific'},
                    'obs_name': {'sst': ['HadISST1.1', 'OISST']},
                },
                'EnsoAlphaSwr': {
                    'variables': ['sst','swr'],
                    'regions': {'sst': 'nino3', 'swr': 'nino3'},
                    'obs_name': {'sst': ['HadISST1.1', 'OISST'], 'swr': ['Tropflux','Tropflux']},
                },
                'EnsoMu': {
                    'variables': ['sst', 'taux'],
                    'regions': {'sst': 'nino3', 'taux': 'nino4'},
                    'obs_name': {'sst': ['HadISST1.1', 'OISST'], 'taux': ['Tropflux', 'ERA-Interim']},
                },
            },
            'common_collection_parameters': {
                'frequency': 'monthly',
                'minimum_number_of_time_steps': 324,
                'observed_period': ('1979-01-01 00:00:00', '2016-12-31 23:59:60.0'),
#                'modeled_period': ('1979-01-01 00:00:00', '2005-12-31 23:59:60.0'),
            },
            'description': 'Describe which science question this collection is about',
        },
    }
    if MC:
        return metrics_collection[MC]
    else:
        return metrics_collection


# List of metrics requirements (var name and reference obs)
# def metricReqs(VOR=True):
#     var_obs_requirements = {
#         'EnsoAlphaThf': {'nbvar': 2, 'var_names': ['sst', 'thf'], 'ref_obs': ['Tropflux', 'Tropflux']},
#         'EnsoAlphaLhf': {'nbvar': 2, 'var_names': ['sst', 'lhf'], 'ref_obs': ['Tropflux', 'Tropflux']},
#         'EnsoAlphaLwr': {'nbvar': 2, 'var_names': ['sst', 'lwr'], 'ref_obs': ['Tropflux', 'Tropflux']},
#         'EnsoAlphaSwr': {'nbvar': 2, 'var_names': ['sst', 'swr'], 'ref_obs': ['Tropflux', 'Tropflux']},
#         'EnsoAmpl': {'nbvar': 1, 'var_names': ['sst'], 'ref_obs': ['OISST']},
#         'EnsoMu': {'nbvar': 2, 'var_names': ['sst', 'taux'], 'ref_obs': ['Tropflux', 'Tropflux']},
#         'EnsoRMSE': {'nbvar': 1, 'var_names': ['sst'], 'ref_obs': ['OISST']},
#         'EnsoSeasonality': {'nbvar': 1, 'var_names': ['sst'], 'ref_obs': ['OISST']},
#     }
#     if VOR:
#         return var_obs_requirements[VOR]
#     else:
#         return var_obs_requirements


# List of reference observations for each variables
def ReferenceObservations(VAR=True, DATASET=True):
    dict_ref_obs = {
        'ERA-Interim': {
            'source': 'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/', # do not use 'source' keyword
            'file_name': '???',
            'variable_name_in_file': {
                'lhf':{'varname': 'hfls'},
                'lwr':{'var_name': 'rls'},
                'shf': {'varname': 'hfss'},
                'sst': {'varname': 'tos'},
                'swr': {'var_name': 'rss'},
                'taux': {'varname': 'taux'},
                # thf: 'algebric_calculation': [+1, +1, +1, +1] if nothing given the variables are summed
                'thf': {'var_name': ['hfls', 'hfss', 'rls', 'rss']},
            },
        },
        'HadISST1.1': {
            'source': 'see https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html',
            'file_name': 'HadISST_' + '<varname>' + '.nc',
            'variable_name_in_file': {
                'sst': {'varname': 'sst'},
            },
        },
        'OISST': {
            'source': 'see https://www.earthsystemcog.org/search/obs4mips/?template=obs4mips&limit=200',
            'file_name': '<varname>' + '_OISST_L4_AVHRR-only-v2_*-*.nc',
            'variable_name_in_file': {
                'sst': {'varname': 'sst'},
            },
        },
        'Tropflux': {
            'source': 'see http://www.incois.gov.in/tropflux_datasets/data/monthly/',
            'file_name': '<varname>' + '_tropflux_1m_*.nc',
            'variable_name_in_file': {
                'lhf': {'varname': 'lhf'},
                'lwr': {'var_name': 'lwr'},
                'shf': {'varname': 'shf'},
                'sst': {'varname': 'sst'},
                'swr': {'var_name': 'swr'},
                'taux': {'varname': 'taux'},
                'thf': {'var_name': 'netflux'},
            },
        },
    }
    if VAR:
        if DATASET:
            return dict_ref_obs[VAR][DATASET]
        else:
            return dict_ref_obs[VAR]
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
    }
    if AR:
        return dict_reference_regions[AR]
    else:
        return dict_reference_regions


def CmipVariables():
    dict_cmip_variables = {
        'reference':'http://cfconventions.org/Data/cf-standard-names/46/build/cf-standard-name-table.html',
        'variable_in_file':{
            # line keys:
            # '<internal_metrics_variable_name>':{'var_name':'<var_name_in_file>','cf_name':<as per ref above>, 'cf_unit':'<unit_in_file>'}

            # latent heat flux (on ocean grid or ocean points only)
            'lhf': {'varname': 'hfls', 'cf_name': 'surface_upward_latent_heat_flux', 'cf_units': 'W m-2'},
            # shortwave radiation (on ocean grid or ocean points only)
            # sometimes the rls = rlds - rlus
            'lwr': {
                'var_name': ['rlds', 'rlus'],
                'cf_name': ['surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air',],
                'cf_units': 'W m-2', 'algebric_calculation': [+1,-1]},
            # sensible heat flux (on ocean grid or ocean points only)
            'shf': {'varname': 'hfss', 'cf_name': 'surface_upward_sensible_heat_flux', 'cf_units': 'W m-2'},
            # sea surface temperature (on ocean grid or ocean points only)
            'sst': {'var_name': 'tos', 'cf_name': 'sea_surface_temperature', 'cf_units': 'K'},
            # shortwave radiation (on ocean grid or ocean points only)
            # sometimes the rss = rsds - rsus
            'swr': {
                'var_name': ['rsds', 'rsus'],
                'cf_name': ['surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': [+1, -1]
            },
            # zonal surface wind stress (on ocean grid or ocean points only)
            'taux': {'var_name': 'tauuo', 'cf_name': 'surface_downward_eastward_stress', 'cf_units': 'Pa'},
            # total heat flux computed from these variables IN THAT ORDER (on ocean grid or ocean points only)
            # tfh = hfls + hfss + rlds - rlus + rsds - rsus
            # sometimes rls = rlds - rlus and rss = rsds - rsus
            'tfh':{
                'var_name': ['hfls', 'hfss', 'rlds', 'rlus', 'rsds', 'rsus'],
                'cf_name': ['surface_upward_latent_heat_flux', 'surface_upward_sensible_heat_flux',
                            'surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air',
                            'surface_downwelling_shortwave_flux_in_air', 'surface_upwelling_shortwave_flux_in_air'],
                'cf_units': 'W m-2', 'algebric_calculation': [+1,+1,+1,-1,+1,-1]
            },
        },
    }
