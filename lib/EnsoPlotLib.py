# -*- coding:UTF-8 -*-
#
# Define ENSO metrics plots
#
from numpy import arange as NUMPYarange

dict_colorbar = {
    'amplitude': 'amp',
    'anomalies': 'balance',
    'PR': 'rain',
    'SST': 'thermal',
}

dict_label = {
    'amplitude': [round(ii, 1) for ii in NUMPYarange(0, 2.1, 0.5)],
    'PR': list(range(0, 13, 4)),
    'PRA': [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)],
    'SKEW': [round(ii, 1) for ii in NUMPYarange(-1.5, 1.6, 0.5)],
    'SST': list(range(22, 31, 2)),
    'SSTA': [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)],
    'TAUX': list(range(-9, 10, 3)),
}

dict_plot_parameters = {
    'BiasPrLatRmse': {
        'netcdf_variables': ['pr_lat__', 'pr_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasPrLatRmse: diagnostic',
            'varpattern': 'pr_lat__',
            'xname': 'latitude',
            'yname': 'PR',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['PR'],
            'label': dict_label['PR'],
            'title': 'Mean PR',
            'varpattern': 'pr_map__',
            'zname': 'PR',
        },
    },
    'BiasPrLonRmse': {
        'netcdf_variables': ['pr_lon__', 'pr_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasPrLonRmse: diagnostic',
            'varpattern': 'pr_lon__',
            'xname': 'latitude',
            'yname': 'PR',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['PR'],
            'label': dict_label['PR'],
            'title': 'Mean PR',
            'varpattern': 'pr_map__',
            'zname': 'PR',
        },
    },
    'BiasPrRmse': {
        'netcdf_variables': ['pr_map__'],
        'diagnostic': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['PR'],
            'label': dict_label['PR'],
            'title': 'BiasPrRmse: diagnostic',
            'varpattern': 'pr_map__',
            'zname': 'PR',
        },
    },
    'BiasSstLatRmse': {
        'netcdf_variables': ['sst_lat__', 'sst_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasSstLatRmse: diagnostic',
            'varpattern': 'sst_lat__',
            'xname': 'latitude',
            'yname': 'SST',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['SST'],
            'label': dict_label['SST'],
            'title': 'Mean SST',
            'varpattern': 'sst_map__',
            'zname': 'SST',
        },
    },
    'BiasSstLonRmse': {
        'netcdf_variables': ['sst_lon__', 'sst_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasSstLonRmse: diagnostic',
            'varpattern': 'sst_lon__',
            'xname': 'latitude',
            'yname': 'SST',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['SST'],
            'label': dict_label['SST'],
            'title': 'Mean SST',
            'varpattern': 'sst_map__',
            'zname': 'SST',
        },
    },
    'BiasSstRmse': {
        'netcdf_variables': ['sst_map__'],
        'diagnostic': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['SST'],
            'label': dict_label['SST'],
            'title': 'BiasSstRmse: diagnostic',
            'varpattern': 'sst_map__',
            'zname': 'SST',
        },
    },
    'BiasSstSkLonRmse': {
        'netcdf_variables': ['sstSke_lon__', 'sstSke_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasSstSkLonRmse: diagnostic',
            'varpattern': 'sstSke_lon__',
            'xname': 'latitude',
            'yname': 'SST skew',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['SKEW'],
            'title': 'SST skewness',
            'varpattern': 'sstSke_map__',
            'zname': 'SST skew',
        },
    },
    'BiasTauxLatRmse': {
        'netcdf_variables': ['taux_lat__', 'taux_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasTauxLatRmse: diagnostic',
            'varpattern': 'pr_lat__',
            'xname': 'latitude',
            'yname': 'TAUX',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['TAUX'],
            'title': 'Mean TAUX',
            'varpattern': 'pr_map__',
            'zname': 'TAUX',
        },
    },
    'BiasTauxLonRmse': {
        'netcdf_variables': ['taux_lon__', 'taux_map__'],
        'diagnostic': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'BiasTauxLonRmse: diagnostic',
            'varpattern': 'pr_lon__',
            'xname': 'latitude',
            'yname': 'TAUX',
        },
        'dive_down01': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['TAUX'],
            'title': 'Mean TAUX',
            'varpattern': 'pr_map__',
            'zname': 'TAUX',
        },
    },
    'BiasTauxRmse': {
        'netcdf_variables': ['taux_map__'],
        'diagnostic': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['TAUX'],
            'title': 'BiasTauxRmse: diagnostic',
            'varpattern': 'taux_map__',
            'zname': 'TAUX',
        },
    },
    'EnsoAmpl': {
        'netcdf_variables': ['sstStd_lon__', 'sstStd_map__'],
        'diagnostic': {
            'plot_type': 'dot',
            'nbr_panel': 1,
            'title': 'EnsoAmpl: diagnostic',
            'varpattern': 'diagnostic',
            'yname': 'SST std',
        },
        'dive_down01': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'SST standard deviation',
            'varpattern': 'sstStd_lon__',
            'xname': 'longitude',
            'yname': 'SST std',
        },
        'dive_down02': {
            'plot_type': 'map',
            'nbr_panel': 2,
            'colorbar': dict_colorbar['amplitude'],
            'label': dict_label['amplitude'],
            'title': 'SST standard deviation',
            'varpattern': 'sstStd_map__',
            'zname': 'SST std',
        },
    },
    'EnsodSstOce': {
        'netcdf_variables': ['dSST_ts__', 'dSSTthf_ts__', 'dSSToce_ts__', 'dSSTthf_lon__', 'dSSToce_lon__',
                             'dSST_hov__', 'dSSTthf_hov__', 'dSSToce_hov__'],
        'diagnostic': {
            'plot_type': 'dot',
            'nbr_panel': 1,
            'title': 'EnsodSstOce: diagnostic',
            'varpattern': 'diagnostic',
            'yname': 'normalized dSSToce',
        },
        'dive_down01': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'ENSO SST change',
            'varpattern': ['dSST_ts__', 'dSSTthf_ts__', 'dSSToce_ts__'],
            'colors': ['black', 'red', 'blue'],
            'legend': ['dSST', 'dSSTthf', 'dSSToce'],
            'xname': 'months',
            'yname': 'normalized dSST',
        },
        'dive_down02': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'ENSO SST change',
            'varpattern': ['dSSTthf_lon__', 'dSSToce_lon__'],
            'colors': ['red', 'blue'],
            'xname': 'longitude',
            'yname': 'normalized dSST',
        },
        'dive_down03': {
            'plot_type': 'hovmoeller',
            'nbr_panel': 6,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['SKEW'],
            'title': ['ENSO SST change', 'ENSO heat flux-driven SST change', 'ENSO ocean-driven SST change'],
            'varpattern': ['dSST_hov__', 'dSSTthf_hov__', 'dSSToce_hov__'],
            'xname': 'longitude',
            'yname': 'months',
            'zname': 'normalized dSST',
        },

    },
    'EnsoDuration': {
        'netcdf_variables': ['sst_over_sst_ts__', 'Nina_duration__', 'Nino_duration__'],
        'diagnostic': {
            'plot_type': 'dot',
            'nbr_panel': 1,
            'title': 'EnsoDuration: diagnostic',
            'varpattern': 'diagnostic',
            'yname': 'duration (reg>0)',
        },
        'dive_down01': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'ENSO life-cycle',
            'varpattern': 'sst_over_sst_ts__',
            'xname': 'months',
            'yname': 'reg(SSTA, SSTA)',
        },
        'dive_down02': {
            'plot_type': 'boxplot',
            'nbr_panel': 2,
            'title': ['La Nina duration', 'El Nino duration'],
            'varpattern': ['Nina_duration__', 'Nino_duration__'],
            'yname': ['duration (SSTA<-0.5)', 'duration (SSTA>0.5)'],
        },
    },
    'EnsoFbSshSst': {
        'netcdf_variables': ['ssh__', 'sst__', 'ssh_over_sst_lon__', 'sshPOS_over_sst_lon__', 'sshNEG_over_sst_lon__',
                             'ssh_over_sst_hov__', 'sshPOS_over_sst_hov__', 'sshNEG_over_sst_hov__'],
        'diagnostic': {
            'plot_type': 'scatterplot',
            'nbr_panel': 1,
            'title': 'EnsoFbSshSst: diagnostic',
            'varpattern': ['ssh__', 'sst__'],
            'xname': 'SSHA',
            'yname': 'SSTA',
        },
        'dive_down01': {
            'plot_type': 'scatterplot',
            'nbr_panel': 2,
            'title': 'EnsoFbSshSst: nonlinarity',
            'varpattern': ['ssh__', 'sst__'],
            'xname': 'SSHA',
            'yname': 'SSTA',
        },
        'dive_down02': {
            'plot_type': 'curve',
            'nbr_panel': 1,
            'title': 'Thermocline feedback',
            'varpattern': ['ssh_over_sst_lon__', 'sshPOS_over_sst_lon__', 'sshNEG_over_sst_lon__'],
            'colors': ['black', 'red', 'blue'],
            'legend': ['All', 'SSHA>0', 'SSHA<0'],
            'xname': 'longitude',
            'yname': 'reg(SSHA, SSTA)',
        },
        'dive_down03': {
            'plot_type': 'hovmoeller',
            'nbr_panel': 6,
            'colorbar': dict_colorbar['anomalies'],
            'label': dict_label['SKEW'],
            'title': ['reg(SSHA, SSTA)', 'reg(SSHA>0, SSTA)', 'reg(SSHA<0, SSTA)'],
            'varpattern': ['ssh_over_sst_hov__', 'sshPOS_over_sst_hov__', 'sshNEG_over_sst_hov__'],
            'xname': 'longitude',
            'yname': 'months',
            'zname': 'regression',
        },
    },

}
# EnsoDiversity


def plot_parameters(metric):
    return dict_plot_parameters[metric]


    # EnsoFbSshSst, EnsoFbSstLhf, EnsoFbSstLwr, EnsoFbSstShf, EnsoFbSstSwr, EnsoFbSstTaux, EnsoFbSstThf, EnsoFbTauxSsh,\
    # EnsoPrMap, EnsoPrJjaTel, EnsoPrNdjTel, EnsoPrTsRmse, EnsoSeasonality, EnsoSlpMap, EnsoSstLonRmse, EnsoSstMap,\
    # EnsoSstSkew, EnsoSstTsRmse, EnsoTauxTsRmse, NinaPrJjaTel, NinaPrNdjTel, NinaPrMap, NinaSlpMap, NinaSstDiv,\
    # NinaSstDivRmse, NinaSstDur, NinaSstLonRmse, NinaSstMap, NinaSstTsRmse, NinoPrJjaTel, NinoPrNdjTel, NinoPrMap,\
    # NinoSlpMap, NinoSstDiv, NinoSstDiversity, NinoSstDivRmse, NinoSstDur, NinoSstLonRmse, NinoSstMap, NinoSstTsRmse,\
    # SeasonalPrLatRmse, SeasonalPrLonRmse, SeasonalSstLatRmse, SeasonalSstLonRmse, SeasonalTauxLatRmse,\
    # SeasonalTauxLonRmse