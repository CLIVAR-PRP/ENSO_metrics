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
#                    'metric_computation': 'difference', # i.e., (obs-model)/model
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
#                    'metric_computation': 'difference',
                },
#                'BiasSstRmse': {
#                    'variables': ['sst'],
#                    'regions': {'sst': 'tropical_pacific'},
#                    'obs_name': {'sst': ['HadISST', 'Tropflux']},
#                    'regridding': {'model_orand_obs': 0, 'regridder': 'cdms', 'regridTool': 'esmf',
#                                   'regridMethod': 'linear'},
#                },
                'EnsoAlphaLhf': {
                    'variables': ['sst','lhf'],
                    'regions': {'sst': 'nino3', 'lhf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'lwr': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
                'EnsoAlphaLwr': {
                    'variables': ['sst', 'lwr'],
                    'regions': {'sst': 'nino3', 'lwr': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'lwr': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
                'EnsoAlphaShf': {
                    'variables': ['sst', 'shf'],
                    'regions': {'sst': 'nino3', 'shf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'shf': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
                'EnsoAlphaSwr': {
                    'variables': ['sst', 'swr'],
                    'regions': {'sst': 'nino3', 'swr': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'swr': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
                'EnsoAlphaThf': {
                    'variables': ['sst', 'thf'],
                    'regions': {'sst': 'nino3', 'thf': 'nino3'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'thf': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
                'EnsoMu': {
                    'variables': ['sst', 'taux'],
                    'regions': {'sst': 'nino3', 'taux': 'nino4'},
                    'obs_name': {'sst': ['HadISST', 'Tropflux'], 'taux': ['Tropflux']},
#                    'metric_computation': 'difference',
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 204,
                'project_interpreter': 'CMIP',
                'metric_computation': 'difference',
                'observed_period': ('1980-01-01 00:00:00', '2016-12-31 23:59:60.0'),#('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                'modeled_period': ('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
            },
            'description': 'Describe which science question this collection is about',
        },
        'ENSO_perf': {
            'long_name': 'Metrics Collection for ENSO performance',
            'metrics_list': {
                'BiasSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasPrLonRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'equatorial_pacific_LonRed'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasTauxLonRmse': {
                    'variables': ['taux'],
                    'regions': {'taux': 'equatorial_pacific_LonRed'},
                    'obs_name': {'taux': ['ERA-Interim', 'Tropflux']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasSstLatRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.3_LatExt'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasPrLatRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'nino3.3_LatExt'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasTauxLatRmse': {
                    'variables': ['taux'],
                    'regions': {'taux': 'equatorial_pacific_LatExt'},
                    'obs_name': {'taux': ['ERA-Interim', 'Tropflux']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'BiasSstSkLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'SeasonalSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'SeasonalSstLatRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.3_LatExt'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'SeasonalPrLatRmse': {
                    'variables': ['pr'],
                    'regions': {'pr': 'nino3.3_LatExt'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3']},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'EnsoAmpl': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'metric_computation': 'ratio',
                },
                # 'EnsoDiversity': {
                #     'variables': ['sst'],
                #     'regions': {'sst': 'equatorial_pacific_LonRed'},
                #     'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'treshold_ep_ev': -140,
                #     'regridding': {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                #                    'newgrid_name': 'generic_1x1deg'},
                #     'metric_computation': 'ratio',
                # },
                'EnsoSeasonality': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'metric_computation': 'ratio',
                },
                'EnsoSstSkew': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'metric_computation': 'ratio',
                },
                'NinaSstDiv': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    '': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75, 'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'treshold_ep_ev': -140,
                    'regridding': {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                   'newgrid_name': 'generic_1x1deg'},
                    'metric_computation': 'ratio',
                },
                'NinaSstDiv2': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    '': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75, 'normalization': True},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'treshold_ep_ev': -140,
                    'regridding': {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                   'newgrid_name': 'generic_1x1deg'},
                    'metric_computation': 'ratio',
                },
                # 'NinaSstDivRmse': {
                #     'variables': ['sst'],
                #     'regions': {'sst': 'equatorial_pacific_LonRed'},
                #     'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                #                    'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                # },
                'NinaSstDur': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'metric_computation': 'ratio',
                },
                'NinaSstDur2': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': True},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'metric_computation': 'ratio',
                },
                'NinaSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'NinaSstTsRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                },
                'NinoSstDiv': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'treshold_ep_ev': -140,
                    'regridding': {'regridder': 'cdms', 'regridTool': 'esmf', 'regridMethod': 'linear',
                                   'newgrid_name': 'generic_1x1deg'},
                    'metric_computation': 'ratio',
                },
                # 'NinoSstDivRmse': {
                #     'variables': ['sst'],
                #     'regions': {'sst': 'equatorial_pacific_LonRed'},
                #     'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                #                    'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                # },
                'NinoSstDur': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'nbr_years_window': 4,
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'metric_computation': 'ratio',
                },
                'NinoSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'NinoSstTsRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim']},
                    'nbr_years_window': 6,
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 204,
                'normalization': False,
                'project_interpreter': 'CMIP',
                'observed_period': ('1917-01-01 00:00:00', '2016-12-31 23:59:60.0'),#('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                'modeled_period': ('1906-01-01 00:00:00', '2005-12-31 23:59:60.0'),
            },
            'plot_groupings': {
                'plot_ratio': ['EnsoAmpl', 'EnsoDiversity', 'EnsoSeasonality', 'EnsoSstSkew', 'NinaSstDiv',
                               'NinaSstDur', 'NinoSstDiv', 'NinoSstDur'],
                'plot_rmse': ['BiasPrLatRmse', 'BiasPrLonRmse', 'BiasSstLatRmse', 'BiasSstLonRmse', 'BiasSstSkLonRmse',
                              'BiasTauxLatRmse', 'BiasTauxLonRmse', 'SeasonalPrLatRmse', 'SeasonalSstLatRmse',
                              'SeasonalSstLonRmse', 'NinaSstLonRmse', 'NinaSstTsRmse', 'NinoSstLonRmse',
                              'NinoSstTsRmse'],
                # 'plot_rmse': ['BiasPrLatRmse', 'BiasPrLonRmse', 'BiasSstLatRmse', 'BiasSstLonRmse', 'BiasSstSkLonRmse',
                #               'BiasTauxLatRmse', 'BiasTauxLonRmse', 'SeasonalPrLatRmse', 'SeasonalSstLatRmse',
                #               'SeasonalSstLonRmse', 'NinaSstDivRmse', 'NinaSstLonRmse', 'NinaSstTsRmse',
                #               'NinoSstDivRmse', 'NinoSstLonRmse', 'NinoSstTsRmse'],
            },
            'description': 'Describe which science question this collection is about',
        },
        'ENSO_tel': {
            'long_name': 'Metrics Collection for ENSO teleconnections',
            'metrics_list': {
                'EnsoAmpl': {
                    'variables': ['sst'],
                    'regions': {'sst': 'nino3.4'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'metric_computation': 'ratio',
                },
                # 'EnsoPrJjaTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'EnsoPrNdjTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'EnsoPrMap': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': 'global', 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                #                    'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'EnsoSlpMap': {
                #     'variables': ['sst', 'slp'],
                #     'regions': {'slp': 'global', 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                #                    'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'EnsoSstMap': {
                #     'variables': ['sst'],
                #     'regions': {'sst': 'global'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                #                    'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                #     'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'NinaPrNdjTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'NinaPrJjaTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                'NinaPrMap': {
                    'variables': ['sst', 'pr'],
                    'regions': {'pr': 'global', 'sst': 'nino3.4'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                },
                'NinaSlpMap': {
                    'variables': ['sst', 'slp'],
                    'regions': {'slp': 'global', 'sst': 'nino3.4'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                },
                'NinaSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'NinaSstMap': {
                    'variables': ['sst'],
                    'regions': {'sst': 'global'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': -0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                },
                # 'NinoPrNdjTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                # 'NinoPrJjaTel': {
                #     'variables': ['sst', 'pr'],
                #     'regions': {'pr': ['CAS', 'CEP', 'CNA', 'CNP', 'CSP', 'EAF', 'EAS', 'ENA', 'INO', 'MED', 'NAU',
                #                        'NEB', 'SAF', 'SAH', 'SAU', 'SEA', 'TIB', 'WAF', 'WAS', 'WNA'],
                #                 'sst': 'nino3.4'},
                #     'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                #                          'normalization': False},
                #     'smoothing': {'window': 5, 'method': 'triangle'},
                #     'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                # },
                'NinoPrMap': {
                    'variables': ['sst', 'pr'],
                    'regions': {'pr': 'global', 'sst': 'nino3.4'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                },
                'NinoSlpMap': {
                    'variables': ['sst', 'slp'],
                    'regions': {'slp': 'global', 'sst': 'nino3.4'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'pr': ['ERA-Interim', 'GPCPv2.3'], 'sst': ['ERA-Interim', 'HadISST']},
                },
                'NinoSstLonRmse': {
                    'variables': ['sst'],
                    'regions': {'sst': 'equatorial_pacific_LonRed'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                },
                'NinoSstMap': {
                    'variables': ['sst'],
                    'regions': {'sst': 'global'},
                    'event_definition': {'region_ev': 'nino3.4', 'season_ev': 'DEC', 'threshold': 0.75,
                                         'normalization': False},
                    'smoothing': {'window': 5, 'method': 'triangle'},
                    'regridding': {'model_orand_obs': 2, 'regridder': 'cdms', 'regridTool': 'esmf',
                                   'regridMethod': 'linear', 'newgrid_name': 'generic_1x1deg'},
                    'obs_name': {'sst': ['ERA-Interim', 'HadISST']},
                },
            },
            'common_collection_parameters': {
                'detrending': {'method': 'linear'},
                'frequency': 'monthly',
                'min_time_steps': 204,
                'normalization': False,
                'project_interpreter': 'CMIP',
                'observed_period': ('1917-01-01 00:00:00', '2016-12-31 23:59:60.0'),#('1905-01-01 00:00:00', '2005-12-31 23:59:60.0'),
                'modeled_period': ('1906-01-01 00:00:00', '2005-12-31 23:59:60.0'),
            },
            'plot_groupings': {
                # 'plot_corr': ['EnsoPrMapCorr', 'EnsoSlpMapCorr', 'EnsoSstMapCorr', 'NinaPrMapCorr', 'NinaSlpMapCorr',
                #               'NinaSstMapCorr', 'NinoPrMapCorr', 'NinoSlpMapCorr', 'NinoSstMapCorr'],
                # 'plot_ratio': ['EnsoPrMapStd', 'EnsoSlpMapStd', 'EnsoSstMapStd', 'NinaPrMapStd', 'NinaSlpMapStd',
                #                'NinaSstMapStd', 'NinoPrMapStd', 'NinoSlpMapStd', 'NinoSstMapStd'],
                # 'plot_rmse': ['EnsoPrJjaTelRmse', 'EnsoPrNdjTelRmse', 'EnsoPrMapRmse', 'EnsoSlpMapRmse',
                #               'EnsoSstMapRmse', 'NinaPrJjaTelRmse', 'NinaPrNdjTelRmse', 'NinaPrMapRmse',
                #               'NinaSlpMapRmse', 'NinaSstMapRmse', 'NinoPrJjaTelRmse', 'NinoPrNdjTelRmse',
                #               'NinoPrMapRmse', 'NinoSlpMapRmse', 'NinoSstMapRmse'],
                # 'plot_percentage': ['EnsoPrJjaTelRmse', 'EnsoPrNdjTelRmse', 'NinaPrJjaTelRmse', 'NinaPrNdjTelRmse',
                #                     'NinoPrJjaTelRmse', 'NinoPrNdjTelRmse']
                'plot_corr': ['NinaPrMapCorr', 'NinaSlpMapCorr', 'NinaSstMapCorr',
                              'NinoPrMapCorr', 'NinoSlpMapCorr', 'NinoSstMapCorr'],
                'plot_ratio': ['EnsoAmpl',
                               'NinaPrMapStd', 'NinaSlpMapStd', 'NinaSstMapStd',
                               'NinoPrMapStd', 'NinoSlpMapStd', 'NinoSstMapStd'],
                'plot_rmse': ['NinaPrMapRmse', 'NinaSlpMapRmse', 'NinaSstLonRmse', 'NinaSstMapRmse',
                              'NinoPrMapRmse', 'NinoSlpMapRmse', 'NinoSstLonRmse', 'NinoSstMapRmse'],
            },
            'description': 'Describe which science question this collection is about',
        },
    }
    if MC is True:
        return metrics_collection
    else:
        return metrics_collection[MC]


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
                'landmask': {'var_name': 'lsmask'},
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
    if DATASET is True:
        return dict_ref_obs
    else:
        return dict_ref_obs[DATASET]


def ReferenceRegions(AR=True):
    dict_reference_regions = {
        'global': {'long_name': 'Global 60S-60N', 'latitude': (-60., 60.), 'longitude': (0., 360.)},
        'global2': {'long_name': 'Global', 'latitude': (-90., 90.), 'longitude': (0., 360.)},
        'tropical_pacific': {
            'long_name': 'Tropical Pacific (TP)', 'latitude': (-30., 30.), 'longitude': (120., 280.),
        },
        'equatorial_pacific': {
            'long_name': 'Equatorial Pacific (EP)', 'latitude': (-5., 5.), 'longitude': (120., 280.),
        },
        'equatorial_pacific_LatExt': {
            'long_name': 'Equatorial Pacific (EP)', 'latitude': (-15., 15.), 'longitude': (150., 270.),
        },
        'equatorial_pacific_LatExt2': {
            'long_name': 'Equatorial Pacific extended in latitude', 'latitude': (-15., 15.), 'longitude': (120., 285.),
        },
        'equatorial_pacific_LonRed': {
            'long_name': 'Equatorial Pacific extended in latitude', 'latitude': (-5., 5.), 'longitude': (150., 270.),
        },
        'eastern_equatorial_pacific': {
            'long_name': 'Western Equatorial Pacific (WEP)', 'latitude': (-5., 5.), 'longitude': (205., 280.),
        },
        'western_equatorial_pacific': {
            'long_name': 'Eastern Equatorial Pacific (EEP)', 'latitude': (-5., 5.), 'longitude': (120., 205.),
        },
        'nino1+2': {'long_name': 'Niño 1+2', 'latitude': (-10., 0.), 'longitude': (270., 280.)},
        'nino3': {'long_name': 'Niño 3', 'latitude': (-5., 5.), 'longitude': (210., 270.)},
        'nino3.3_LatExt': {
            'long_name': 'small Niño 3 extended in latitude', 'latitude': (-15., 15.), 'longitude': (220., 250.)
        },
        'nino3.4': {'long_name': 'Niño 3.4', 'latitude': (-5., 5.), 'longitude': (190., 240.)},
        'nino4': {'long_name': 'Niño 4', 'latitude': (-5., 5.), 'longitude': (160., 210.)},
        # AR5 reference regions
        'ALA': {'long_name': 'Alaska/N.W. Canada', 'latitude': (60., 72.6), 'longitude': (192., 255.),
                'maskland': False, 'maskocean': True},
        # 'AMZ': {'long_name': 'Amazon', 'polygon shaped region, I do not know how to select it'},
        # 'CAM': {'long_name': 'Central America/Mexico', 'polygon shaped region, I do not know how to select it'},
        'CAS': {'long_name': 'Central Asia', 'latitude': (30., 50.), 'longitude': (60., 75.), 'maskland': False,
                'maskocean': True},
        # 'CEU': {'long_name': 'Central Europe', 'polygon shaped region, I do not know how to select it'},
        'CGI': {'long_name': 'Canada/Greenland/Iceland', 'latitude': (50., 85.), 'longitude': (255., 350.),
                'maskland': False, 'maskocean': True},
        'CNA': {'long_name': 'Central North America', 'latitude': (28.6, 50.), 'longitude': (255., 275.),
                'maskland': False, 'maskocean': True},
        'EAF': {'long_name': 'East Africa', 'latitude': (-11.4, 15.), 'longitude': (25., 52.), 'maskland': False,
                'maskocean': True},
        'EAS': {'long_name': 'East Asia', 'latitude': (20., 50.), 'longitude': (100., 145.), 'maskland': False,
                'maskocean': True},
        'ENA': {'long_name': 'East North America', 'latitude': (25., 50.), 'longitude': (275., 300.), 'maskland': False,
                'maskocean': True},
        'MED': {'long_name': 'South Europe/Mediterranean', 'latitude': (30., 45.), 'longitude': (350., 400.),
                'maskland': False, 'maskocean': True},
        'NAS': {'long_name': 'North Asia', 'latitude': (50., 70.), 'longitude': (40., 180.), 'maskland': False,
                'maskocean': True},
        'NAU': {'long_name': 'North Australia', 'latitude': (-30., -10.), 'longitude': (110., 155.), 'maskland': False,
                'maskocean': True},
        'NEB': {'long_name': 'North-East Brazil', 'latitude': (-20., 0.), 'longitude': (310., 326.), 'maskland': False,
                'maskocean': True},
        # 'NEU': {'long_name': 'North Europe', 'polygon shaped region, I do not know how to select it'},
        'SAF': {'long_name': 'Southern Africa', 'latitude': (-35., -11.4), 'longitude': (350., 412.), 'maskland': False,
                'maskocean': True},
        'SAH': {'long_name': 'Sahara', 'latitude': (15., 30.), 'longitude': (340., 400.), 'maskland': False,
                'maskocean': True},
        # 'SAS': {'long_name': 'South Asia', 'polygon shaped region, I do not know how to select it'},
        'SAU': {'long_name': 'South Australia/New Zealand', 'latitude': (-50., -30.), 'longitude': (110., 180.),
                'maskland': False, 'maskocean': True},
        'SEA': {'long_name': 'Southeast Asia', 'latitude': (-10., 20.), 'longitude': (95., 155.), 'maskland': False,
                'maskocean': False},
        # 'SSA': {'long_name': 'Southeastern South America', 'polygon shaped region, I do not know how to select it'},
        'TIB': {'long_name': 'Tibetan Plateau', 'latitude': (30., 50.), 'longitude': (75., 100.), 'maskland': False,
                'maskocean': True},
        'WAF': {'long_name': 'West Africa', 'latitude': (-11.4, 15.), 'longitude': (340., 385.), 'maskland': False,
                'maskocean': True},
        'WAS': {'long_name': 'West Asia', 'latitude': (15., 50.), 'longitude': (40., 60.), 'maskland': False,
                'maskocean': True},
        'WNA': {'long_name': 'West North America', 'latitude': (28.6, 60.), 'longitude': (230., 255.),
                'maskland': False, 'maskocean': True},
        # 'WSA': {'long_name': 'West Coast South America', 'polygon shaped region, I do not know how to select it'},
        # non-SREX reference regions
        'ANT': {'long_name': 'Antarctica', 'latitude': (-90., -50.), 'longitude': (0., 360.), 'maskland': False,
                'maskocean': False},
        'ARC': {'long_name': 'Arctic', 'latitude': (67.5, 90.), 'longitude': (0., 360.), 'maskland': False,
                'maskocean': False},
        # 'CAR': {
        #     'long_name': 'Caribbean', 'polygon shaped region, I do not know how to select it'
        # },
        'NTP': {'long_name': 'Northern Tropical Pacific', 'latitude': (5., 25.), 'longitude': (155., 210.)},
        'STP': {'long_name': 'Southern Topical Pacific', 'latitude': (-25., -5.), 'longitude': (155., 230.)},
        'ETP': {'long_name': 'Equatorial Tropical Pacific', 'latitude': (-5., 5.), 'longitude': (155., 210.)},
        'WIO': {'long_name': 'West Indian Ocean', 'latitude': (-25., 5.), 'longitude': (52., 75.), 'maskland': False,
                'maskocean': False},
        # Power and Delage's (2018) oceanic regions
        'CEP': {'long_name': 'Central Equatorial Pacific', 'latitude': (-5., 5.), 'longitude': (180., 220.),
                'maskland': True, 'maskocean': False},
        'CNP': {'long_name': 'Central Northern Tropical Pacific', 'latitude': (5., 15.), 'longitude': (180., 220.),
                'maskland': True, 'maskocean': False},
        'CSP': {'long_name': 'Central Southern Tropical Pacific', 'latitude': (-15., -5.), 'longitude': (180., 220.),
                'maskland': True, 'maskocean': False},
        'INO': {'long_name': 'Indian Ocean', 'latitude': (-25., 0.), 'longitude': (55., 95.), 'maskland': True,
                'maskocean': False},
    }
    if AR is True:
        return dict_reference_regions
    else:
        return dict_reference_regions[AR]


def CmipVariables():
    dict_cmip_variables = {
        'reference': 'http://cfconventions.org/Data/cf-standard-names/46/build/cf-standard-name-table.html',
        'variable_name_in_file': {
            # line keys:
            # '<internal_metrics_variable_name>':{'var_name':'<var_name_in_file>','cf_name':<as per ref above>, 'cf_unit':'<unit_in_file>'}
            # areacell
            'areacell': {'var_name': 'areacella', 'cf_name': 'cell_area', 'cf_units': 'm2'},
            # landmask
            'landmask': {'var_name': 'sftlf', 'cf_name': 'cell_area', 'cf_units': '1'},
            # latent heat flux
            'lhf': {'var_name': 'hfls', 'cf_name': 'land_area_fraction', 'cf_units': 'W m-2'},
            # longwave radiation computed from these variables IN THAT ORDER
            # lwr = rlds - rlus
            # sometimes lwr is included in the datasets in a variable called 'rls'
            'lwr': {
                'var_name': ['rlds', 'rlus'],
                'cf_name': ['surface_downwelling_longwave_flux_in_air', 'surface_upwelling_longwave_flux_in_air',],
                'cf_units': 'W m-2', 'algebric_calculation': ['plus', 'minus']},
            # Rainfall Flux
            'pr': {'var_name': 'pr', 'cf_name': 'rainfall_flux', 'cf_units': 'kg m-2 s-1'},
            # Sea Level Pressure
            'slp': {'var_name': 'psl', 'cf_name': 'air_pressure_at_mean_sea_level', 'cf_units': 'Pa'},
            # sensible heat flux (on ocean grid or ocean points only)
            'shf': {'var_name': 'hfss', 'cf_name': 'surface_upward_sensible_heat_flux', 'cf_units': 'W m-2'},
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

