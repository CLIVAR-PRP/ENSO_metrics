# -*- coding:UTF-8 -*-
#
# Define ENSO metrics plots
#
from numpy import arange as NUMPYarange
# ENSO_metrics functions
from .EnsoCollectionsLib import defCollection
from .KeyArgLib import default_arg_values


dict_colorbar = {
    "amplitude": "amp",
    "anomalies": "balance",
    "PR": "rain",
    "SST": "thermal",
}

dict_label = {
    "amplitude": [round(ii, 1) for ii in NUMPYarange(0, 2.1, 0.5)],
    "amplitude5": list(range(0, 6, 1)),
    "amplitude10": [round(ii, 1) for ii in NUMPYarange(0, 10.1, 2.5)],
    "amplitude15": list(range(0, 16, 5)),
    "amplitude60": list(range(0, 61, 20)),
    "PR": list(range(0, 13, 4)),
    "PRA": [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)],
    "REG03": [round(ii, 1) for ii in NUMPYarange(-0.3, 0.35, 0.1)],
    "REG05": [round(ii, 2) for ii in NUMPYarange(-0.5, 0.55, 0.25)],
    "REG12": [round(ii, 1) for ii in NUMPYarange(-1.2, 1.4, 0.6)],
    "REG2": list(range(-2, 3, 1)),
    "REG25": [round(ii, 1) for ii in NUMPYarange(-2.5, 2.6, 1.0)],
    "REG3": list(range(-3, 4, 1)),
    "REG4": list(range(-4, 5, 1)),
    "REG5": [round(ii, 1) for ii in NUMPYarange(-5, 6, 2.5)],
    "REG20": list(range(-20, 25, 10)),
    "REG30": list(range(-30, 35, 15)),
    "REG50": list(range(-50, 55, 25)),
    "REG60": list(range(-60, 65, 30)),
    "REG80": list(range(-80, 85, 40)),
    "SKEW": [round(ii, 1) for ii in NUMPYarange(-1.5, 1.6, 0.5)],
    "dSST": list(range(-2, 3, 1)),
    "SST": list(range(21, 31, 3)),
    "SSTA": [round(ii, 1) for ii in NUMPYarange(-1, 1.1, 0.5)],
    "TAUX": list(range(-100, 110, 50)),
}

plot_parameters = {
    "BiasPrLatRmse": {
        "netcdf_variables": ["pr_lat__", "pr_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean PR",  # "a) Mean meridional PR",  #
            "varpattern": "pr_lat__",
            "xname": "latitude",
            "yname": "PR",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["PR"],
            "label": dict_label["PR"],
            "maskland": True,
            "title": ["Mean PR", "Mean PR"],
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PR",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasPrLonRmse": {
        "netcdf_variables": ["pr_lon__", "pr_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean PR",
            "varpattern": "pr_lon__",
            "xname": "longitude",
            "yname": "PR",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["PR"],
            "label": dict_label["PR"],
            "maskland": True,
            "title": ["Mean PR", "Mean PR"],
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PR",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasPrRmse": {
        "netcdf_variables": ["pr_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["PR"],
            "label": dict_label["PR"],
            "maskland": True,
            "title": ["Mean PR", "Mean PR"],
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PR",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n\nMetric: RMSE$_{xy}$",
        },
    },
    "BiasSshLatRmse": {
        "netcdf_variables": ["ssh_lat__", "ssh_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean SSH",
            "varpattern": "ssh_lat__",
            "xname": "latitude",
            "yname": "SSH",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SSH", "Mean SSH"],
            "varpattern": "ssh_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSH",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasSshLonRmse": {
        "netcdf_variables": ["ssh_lon__", "ssh_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean SSH",
            "varpattern": "ssh_lon__",
            "xname": "longitude",
            "yname": "SSH",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SSH", "Mean SSH"],
            "varpattern": "ssh_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSH",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasSshRmse": {
        "netcdf_variables": ["ssh_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SSH", "Mean SSH"],
            "varpattern": "ssh_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSH",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n\nMetric: RMSE$_{xy}$",
        },
    },
    "BiasSstLatRmse": {
        "netcdf_variables": ["sst_lat__", "sst_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean SST",
            "varpattern": "sst_lat__",
            "xname": "latitude",
            "yname": "SST",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SST", "Mean SST"],
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SST",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasSstLonRmse": {
        "netcdf_variables": ["sst_lon__", "sst_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean SST",
            "varpattern": "sst_lon__",
            "xname": "longitude",
            "yname": "SST",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SST", "Mean SST"],
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SST",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasSstRmse": {
        "netcdf_variables": ["sst_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],
            "label": dict_label["SST"],
            "maskland": True,
            "title": ["Mean SST", "Mean SST"],
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SST",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n\nMetric: RMSE$_{xy}$",
        },
    },
    "BiasTauxLatRmse": {
        "netcdf_variables": ["taux_lat__", "taux_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean TAUX",
            "varpattern": "taux_lat__",
            "xname": "latitude",
            "yname": "TAUX",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["TAUX"],
            "maskland": True,
            "title": ["Mean TAUX", "Mean TAUX"],
            "varpattern": "taux_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TAUX",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasTauxLonRmse": {
        "netcdf_variables": ["taux_lon__", "taux_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Mean TAUX",
            "varpattern": "taux_lon__",
            "xname": "longitude",
            "yname": "TAUX",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n" +
                      "4) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["TAUX"],
            "maskland": True,
            "title": ["Mean TAUX", "Mean TAUX"],
            "varpattern": "taux_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TAUX",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°",
        },
    },
    "BiasTauxRmse": {
        "netcdf_variables": ["taux_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["TAUX"],
            "maskland": True,
            "title": ["Mean TAUX", "Mean TAUX"],
            "varpattern": "taux_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TAUX",
            "method": "1) Linearly detrended\n2) Temporal averaged\n3) Regridded to 1°x1°\n\nMetric: RMSE$_{xy}$",
        },
    },
    "EnsoAmpl": {
        "netcdf_variables": ["sstStd_lon__", "sstStd_map__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO amplitude",
            "varpattern": "diagnostic",
            "yname": "SSTA std",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) Standard deviation\n\nMetric: abs((STD$_{mod}$-STD$_{ref}$)/STD$_{ref}$)*100",

        },
        "dive_down01": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSTA standard deviation",
            "varpattern": "sstStd_lon__",
            "xname": "longitude",
            "yname": "SSTA std",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) 5S-5N meridional averaged",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["SSTA standard deviation", "SSTA standard deviation"],
            "varpattern": "sstStd_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA std",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
    },
    "EnsodSstOce": {
        "netcdf_variables": ["dSST_ts__", "dSSTthf_ts__", "dSSToce_ts__", "dSSTthf_lon__", "dSSToce_lon__",
                             "dSST_hov__", "dSSTthf_hov__", "dSSToce_hov__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO ocean-driven SST change",
            "varpattern": "diagnostic",
            "yname": "normalized dSSToce",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) dSST = REGION1 SSTA Dec. - Jul.\n3) REGION1 NHFA summed from Jul. to Dec.\n" +
                      "4) dSST = dSST/dSST and dSSTnhf = dSSTnhf/dSST\n" +
                      "5) dSSToce = dSST - dSSTnhf\n6) Mean dSSToce both El Nino and La Nina\n\n" +
                      "Metric: abs((dSSToce$_{mod}$-dSSToce$_{ref}$)/dSSToce$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO SST change",
            "varpattern": ["dSST_ts__", "dSSTthf_ts__", "dSSToce_ts__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["dSST", "dSSTnhf", "dSSToce"],
            "xname": "months",
            "yname": "normalized dSST",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) dSST$_i$ = REGION1 SSTA M$_i$ - M$_{i-1}$\n" +
                      "3) dSSTnhf$_i$ = REGION1 NHFA summed from M$_0$ to M$_i$\n" +
                      "4) dt = dSST$_{dec}$-dSST$_{jul}$,\n    dSST$_i$ = dSST$_i$/dt,\n    " +
                      "dSSTnhf$_i$ = dSSTnhf$_i$/dt\n" +
                      "5) dSSToce$_i$ = dSST$_i$ - dSSTnhf$_i$\n6) Mean dSSToce both El Nino and La Nina",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO SST change",
            "varpattern": ["dSSTthf_lon__", "dSSToce_lon__"],
            "colors": {"model": ["red", "blue"], "reference": ["red", "blue"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["dSSTnhf", "dSSToce"],
            "xname": "longitude",
            "yname": "normalized dSST",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) dSST = REGION1 SSTA Dec. - Jul.\n3) REGION1 NHFA summed from Jul. to Dec.\n" +
                      "4) dSST = dSST/dSST and dSSTnhf = dSSTnhf/dSST\n" +
                      "5) dSSToce = dSST - dSSTnhf\n6) Mean dSSToce both El Nino and La Nina",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["dSST"],
            "title": ["ENSO dSST", "ENSO heat flux dSST", "ENSO ocean dSST"],
            "varpattern": ["dSST_hov__", "dSSTthf_hov__", "dSSToce_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "normalized dSST",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) dSST$_i$ = REGION1 SSTA M$_i$ - M$_{i-1}$\n" +
                      "3) dSSTnhf$_i$ = REGION1 NHFA summed from M$_0$ to M$_i$\n" +
                      "4) dt = dSST$_{dec}$-dSST$_{jul}$,\n    dSST$_i$ = dSST$_i$/dt,\n    " +
                      "dSSTnhf$_i$ = dSSTnhf$_i$/dt\n" +
                      "5) dSSToce$_i$ = dSST$_i$ - dSSTnhf$_i$\n6) Mean dSSToce both El Nino and La Nina",
        },

    },
    "EnsoDuration": {
        "netcdf_variables": ["sst_against_sst_ts__", "Nina_duration__", "Nino_duration__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO duration",
            "varpattern": "diagnostic",
            "yname": "duration (reg>0.25)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) REGION1 averaged\n5) Dec. REGION1 SSTA regressed onto REGION1 SSTA\n" +
                      "6) Duration = nbr months > 0.25\n\nMetric: abs((DUR$_{mod}$-DUR$_{ref}$)/DUR$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO life-cycle",
            "varpattern": "sst_against_sst_ts__", #"sst_over_sst_ts__",
            "xname": "months",
            "yname": "reg(SSTA, SSTA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) REGION1 averaged\n5) Dec. REGION1 SSTA regressed onto REGION1 SSTA",
        },
        "dive_down02": {
            "plot_type": "boxplot",
            "nbr_panel": 2,
            "title": ["La Nina duration", "El Nino duration"],
            "varpattern": ["Nina_duration__", "Nino_duration__"],
            "yname": ["duration (SSTA<-0.5)", "duration (SSTA>0.5)"],
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 averaged\n6) Duration = nbr months > 0.5 ENSO STD",
        },
    },
    "EnsoFbSshSst": {
        "netcdf_variables": ["ssh__", "sst__", "ssh_over_sst_lon__", "sshPOS_over_sst_lon__", "sshNEG_over_sst_lon__",
                             "ssh_over_sst_hov__", "sshPOS_over_sst_hov__", "sshNEG_over_sst_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "SSH-to-SST coupling",
            "varpattern": ["ssh__", "sst__"],
            "xname": "SSHA",
            "yname": "SSTA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSHA regressed onto REGION1 SSTA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["ssh__", "sst__"],
            "xname": "SSHA",
            "yname": "SSTA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSHA>0 (SSHA<0) regressed onto REGION1 SSTA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Thermocline feedback",
            #"varpattern": ["ssh_over_sst_lon__", "sshPOS_over_sst_lon__", "sshNEG_over_sst_lon__"],
            "varpattern": ["reg_sst_over_ssh_lon__", "reg_sst_over_POSssh_lon__", "reg_sst_over_NEGssh_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSHA>0", "SSHA<0"],
            "xname": "longitude",
            "yname": "reg(SSHA, SSTA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSHA or SSHA>0 or SSHA<0 regressed onto SSTA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG03"],
            "title": ["reg(SSHA, SSTA)", "reg(SSHA>0, SSTA)", "reg(SSHA<0, SSTA)"],
            #"varpattern": ["ssh_over_sst_hov__", "sshPOS_over_sst_hov__", "sshNEG_over_sst_hov__"],
            "varpattern": ["reg_sst_over_ssh_hov__", "reg_sst_over_POSssh_hov__", "reg_sst_over_NEGssh_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSHA or SSHA>0 or SSHA<0 regressed onto SSTA",
        },
    },
    "EnsoFbSstLhf": {
        "netcdf_variables": ["sst__", "lhf__", "sst_over_lhf_lon__", "sstPOS_over_lhf_lon__", "sstNEG_over_lhf_lon__",
                             "sst_over_lhf_hov__", "sstPOS_over_lhf_hov__", "sstNEG_over_lhf_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Latent heat feedback",
            "varpattern": ["sst__", "lhf__"],
            "xname": "SSTA",
            "yname": "LHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION1 LHFA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "lhf__"],
            "xname": "SSTA",
            "yname": "LHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION1 LHFA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Latent heat feedback",
            #"varpattern": ["sst_over_lhf_lon__", "sstPOS_over_lhf_lon__", "sstNEG_over_lhf_lon__"],
            "varpattern": ["reg_lhf_over_sst_lon__", "reg_lhf_over_POSsst_lon__", "reg_lhf_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, LHFA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSTA or SSTA>0 or SSTA<0 regressed onto LHFA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG20"],
            "title": ["reg(SSTA, LHFA)", "reg(SSTA>0, LHFA)", "reg(SSTA<0, LHFA)"],
            #"varpattern": ["sst_over_lhf_hov__", "sstPOS_over_lhf_hov__", "sstNEG_over_lhf_hov__"],
            "varpattern": ["reg_lhf_over_sst_hov__", "reg_lhf_over_POSsst_hov__", "reg_lhf_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSTA or SSTA>0 or SSTA<0 regressed onto LHFA",
        },
    },
    "EnsoFbSstLwr": {
        "netcdf_variables": ["sst__", "lwr__", "sst_over_lwr_lon__", "sstPOS_over_lwr_lon__", "sstNEG_over_lwr_lon__",
                             "sst_over_lwr_hov__", "sstPOS_over_lwr_hov__", "sstNEG_over_lwr_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Longwave feedback",
            "varpattern": ["sst__", "lwr__"],
            "xname": "SSTA",
            "yname": "LWRA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION1 LWRA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "lwr__"],
            "xname": "SSTA",
            "yname": "LWRA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION1 LWRA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Longwave feedback",
            #"varpattern": ["sst_over_lwr_lon__", "sstPOS_over_lwr_lon__", "sstNEG_over_lwr_lon__"],
            "varpattern": ["reg_lwr_over_sst_lon__", "reg_lwr_over_POSsst_lon__", "reg_lwr_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, LWRA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSTA or SSTA>0 or SSTA<0 regressed onto LWRA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG20"],
            "title": ["reg(SSTA, LWRA)", "reg(SSTA>0, LWRA)", "reg(SSTA<0, LWRA)"],
            #"varpattern": ["sst_over_lwr_hov__", "sstPOS_over_lwr_hov__", "sstNEG_over_lwr_hov__"],
            "varpattern": ["reg_lwr_over_sst_hov__", "reg_lwr_over_POSsst_hov__", "reg_lwr_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSTA or SSTA>0 or SSTA<0 regressed onto LWRA",
        },
    },
    "EnsoFbSstShf": {
        "netcdf_variables": ["sst__", "shf__", "sst_over_shf_lon__", "sstPOS_over_shf_lon__", "sstNEG_over_shf_lon__",
                             "sst_over_shf_hov__", "sstPOS_over_shf_hov__", "sstNEG_over_shf_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Sensible heat feedback",
            "varpattern": ["sst__", "shf__"],
            "xname": "SSTA",
            "yname": "SHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION1 SHFA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "shf__"],
            "xname": "SSTA",
            "yname": "SHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION1 SHFA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Sensible heat feedback",
            #"varpattern": ["sst_over_shf_lon__", "sstPOS_over_shf_lon__", "sstNEG_over_shf_lon__"],
            "varpattern": ["reg_shf_over_sst_lon__", "reg_shf_over_POSsst_lon__", "reg_shf_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, SHFA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSTA or SSTA>0 or SSTA<0 regressed onto SHFA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG5"],
            "title": ["reg(SSTA, SHFA)", "reg(SSTA>0, SHFA)", "reg(SSTA<0, SHFA)"],
            #"varpattern": ["sst_over_shf_hov__", "sstPOS_over_shf_hov__", "sstNEG_over_shf_hov__"],
            "varpattern": ["reg_shf_over_sst_hov__", "reg_shf_over_POSsst_hov__", "reg_shf_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSTA or SSTA>0 or SSTA<0 regressed onto SHFA",
        },
    },
    "EnsoFbSstSwr": {
        "netcdf_variables": ["sst__", "swr__", "sst_over_swr_lon__", "sstPOS_over_swr_lon__", "sstNEG_over_swr_lon__",
                             "sst_over_swr_hov__", "sstPOS_over_swr_hov__", "sstNEG_over_swr_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Shortwave feedback",
            "varpattern": ["sst__", "swr__"],
            "xname": "SSTA",
            "yname": "SWRA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION1 SWRA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "swr__"],
            "xname": "SSTA",
            "yname": "SWRA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION1 SWRA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Shortwave feedback",
            #"varpattern": ["sst_over_swr_lon__", "sstPOS_over_swr_lon__", "sstNEG_over_swr_lon__"],
            "varpattern": ["reg_swr_over_sst_lon__", "reg_swr_over_POSsst_lon__", "reg_swr_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, SWRA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSTA or SSTA>0 or SSTA<0 regressed onto SWRA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG50"],
            "title": ["reg(SSTA, SWRA)", "reg(SSTA>0, SWRA)", "reg(SSTA<0, SWRA)"],
            #"varpattern": ["sst_over_swr_hov__", "sstPOS_over_swr_hov__", "sstNEG_over_swr_hov__"],
            "varpattern": ["reg_swr_over_sst_hov__", "reg_swr_over_POSsst_hov__", "reg_swr_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSTA or SSTA>0 or SSTA<0 regressed onto SWRA",
        },
    },
    "EnsoFbSstTaux": {
        "netcdf_variables": [
            "sst__", "taux__", "sst_over_taux_lon__", "sstPOS_over_taux_lon__", "sstNEG_over_taux_lon__",
            "sst_over_taux_hov__", "sstPOS_over_taux_hov__", "sstNEG_over_taux_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "SST-to-Taux coupling",
            "varpattern": ["sst__", "taux__"],
            "xname": "SSTA",
            "yname": "TAUXA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Horizontal averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION2 TAUXA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "taux__"],
            "xname": "SSTA",
            "yname": "TAUXA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Horizontal averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION2 TAUXA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Wind stress feedback",
            #"varpattern": ["sst_over_taux_lon__", "sstPOS_over_taux_lon__", "sstNEG_over_taux_lon__"],
            "varpattern": ["reg_taux_over_sst_lon__", "reg_taux_over_POSsst_lon__", "reg_taux_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, TAUXA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 SSTA averaged\n" +
                      "4) TAUXA 5S-5N meridional averaged\n5) TAUXA 30° zonal running ave.\n" +
                      "6) REGION1 SSTA or SSTA>0 or SSTA<0 regressed onto TAUXA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG20"],
            "title": ["reg(SSTA, TAUXA)", "reg(SSTA>0, TAUXA)", "reg(SSTA<0, TAUXA)"],
            #"varpattern": ["sst_over_taux_hov__", "sstPOS_over_taux_hov__", "sstNEG_over_taux_hov__"],
            "varpattern": ["reg_taux_over_sst_hov__", "reg_taux_over_POSsst_hov__", "reg_taux_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 SSTA averaged\n" +
                      "4) TAUXA 5S-5N meridional averaged\n5) TAUXA 30° zonal running ave.\n" +
                      "6) For each calendar month:\n    REGION1 SSTA or SSTA>0 or SSTA<0 regressed onto TAUXA",
        },
    },
    "EnsoFbSstThf": {
        "netcdf_variables": ["sst__", "thf__", "sst_over_thf_lon__", "sstPOS_over_thf_lon__", "sstNEG_over_thf_lon__",
                             "sst_over_thf_hov__", "sstPOS_over_thf_hov__", "sstNEG_over_thf_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Net heat flux feedback",
            "varpattern": ["sst__", "thf__"],
            "xname": "SSTA",
            "yname": "NHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA regressed onto REGION1 NHFA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["sst__", "thf__"],
            "xname": "SSTA",
            "yname": "NHFA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) REGION1 SSTA>0 (SSTA<0) regressed onto REGION1 NHFA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "Net heat flux feedback",
            #"varpattern": ["sst_over_thf_lon__", "sstPOS_over_thf_lon__", "sstNEG_over_thf_lon__"],
            "varpattern": ["reg_thf_over_sst_lon__", "reg_thf_over_POSsst_lon__", "reg_thf_over_NEGsst_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "SSTA>0", "SSTA<0"],
            "xname": "longitude",
            "yname": "reg(SSTA, NHFA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n5) SSTA or SSTA>0 or SSTA<0 regressed onto NHFA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG50"],
            "title": ["reg(SSTA, NHFA)", "reg(SSTA>0, NHFA)", "reg(SSTA<0, NHFA)"],
            #"varpattern": ["sst_over_thf_hov__", "sstPOS_over_thf_hov__", "sstNEG_over_thf_hov__"],
            "varpattern": ["reg_thf_over_sst_hov__", "reg_thf_over_POSsst_hov__", "reg_thf_over_NEGsst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5S-5N meridional averaged\n" +
                      "4) 30° zonal running ave.\n" +
                      "5) For each calendar month:\n    SSTA or SSTA>0 or SSTA<0 regressed onto NHFA",
        },
    },
    "EnsoFbTauxSsh": {
        "netcdf_variables": [
            "taux__", "ssh__", "taux_over_ssh_lon__", "tauxPOS_over_ssh_lon__", "tauxNEG_over_ssh_lon__",
            "taux_over_ssh_hov__", "tauxPOS_over_ssh_hov__", "tauxNEG_over_ssh_hov__"],
        "diagnostic": {
            "plot_type": "scatterplot",
            "nbr_panel": 1,
            "title": "Taux-to-SSH coupling",
            "varpattern": ["taux__", "ssh__"],
            "xname": "TAUXA",
            "yname": "SSHA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Horizontal averaged\n" +
                      "4) REGION1 TAUXA regressed onto REGION2 SSHA\n\n" +
                      "Metric: abs((Slope$_{mod}$-Slope$_{ref}$)/Slope$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "scatterplot",
            "nbr_panel": 2,
            "title": "nonlinarity",
            "varpattern": ["taux__", "ssh__"],
            "xname": "TAUXA",
            "yname": "SSHA",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Horizontal averaged\n" +
                      "4) REGION1 TAUXA>0 (TAUXA<0) regressed onto REGION2 SSHA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSH-Wind feedback",
            #"varpattern": ["taux_over_ssh_lon__", "tauxPOS_over_ssh_lon__", "tauxNEG_over_ssh_lon__"],
            "varpattern": ["reg_ssh_over_taux_lon__", "reg_ssh_over_POStaux_lon__", "reg_ssh_over_NEGtaux_lon__"],
            "colors": {"model": ["black", "red", "blue"], "reference": ["black", "red", "blue"]},
            "linestyles": {"model": ["-", "-", "-"], "reference": ["-.", "-.", "-."]},
            "legend": ["All", "TAUXA>0", "TAUXA<0"],
            "xname": "longitude",
            "yname": "reg(TAUXA, SSHA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 TAUXA averaged\n" +
                      "4) SSHA 5S-5N meridional averaged\n5) SSHA 30° zonal running ave.\n" +
                      "6) REGION1 TAUXA or TAUXA>0 or TAUXA<0 regressed onto SSHA",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 6,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG05"],
            "title": ["reg(TAUXA, SSHA)", "reg(TAUXA>0, SSHA)", "reg(TAUXA<0, SSHA)"],
            #"varpattern": ["taux_over_ssh_hov__", "tauxPOS_over_ssh_hov__", "tauxNEG_over_ssh_hov__"],
            "varpattern": ["reg_ssh_over_taux_hov__", "reg_ssh_over_POStaux_hov__", "reg_ssh_over_NEGtaux_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 TAUXA averaged\n" +
                      "4) SSHA 5S-5N meridional averaged\n5) SSHA 30° zonal running ave.\n" +
                      "6) For each calendar month:\n    REGION1 TAUXA or TAUXA>0 or TAUXA<0 regressed onto SSHA",
        },
    },
    "EnsoPrMap": {
        "netcdf_variables": ["reg_pr_over_sst_map__", "reg_pr_over_sst_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "reg_pr_over_sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            "varpattern": "reg_pr_over_sst_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Mask ocean",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            "varpattern": "reg_pr_over_sst_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            "varpattern": "reg_pr_over_sst_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            "varpattern": "reg_pr_over_sst_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            "varpattern": "reg_pr_over_sst_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. PRA\n" +
                      "7) Mask ocean",
        },
    },
    "EnsoPrMapDjf": {
        "netcdf_variables": ["reg_pr_over_sst_djf_map__", "reg_pr_over_sst_djf_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "reg_pr_over_sst_djf_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map__", "pr_nino_djf_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            "varpattern": "reg_pr_over_sst_djf_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            "varpattern": "reg_pr_over_sst_djf_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            "varpattern": "reg_pr_over_sst_djf_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            "varpattern": "reg_pr_over_sst_djf_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) DJF", "reg(ENSO SSTA, PRA) DJF"],
            "varpattern": "reg_pr_over_sst_djf_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF PRA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map_africaSE__", "pr_nino_djf_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map_americaN__", "pr_nino_djf_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map_americaS__", "pr_nino_djf_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map_asiaS__", "pr_nino_djf_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA DJF", "El Nino PRA DJF"],
            "varpattern": ["pr_nina_djf_map_oceania__", "pr_nino_djf_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina DJF PRA composited\n7) Mask ocean",
        },
    },
    "EnsoPrMapJja": {
        "netcdf_variables": ["reg_pr_over_sst_jja_map__", "reg_pr_over_sst_jja_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "reg_pr_over_sst_jja_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map__", "pr_nino_jja_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            "varpattern": "reg_pr_over_sst_jja_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            "varpattern": "reg_pr_over_sst_jja_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            "varpattern": "reg_pr_over_sst_jja_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            "varpattern": "reg_pr_over_sst_jja_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, PRA) JJA", "reg(ENSO SSTA, PRA) JJA"],
            "varpattern": "reg_pr_over_sst_jja_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) PRA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA PRA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map_africaSE__", "pr_nino_jja_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map_americaN__", "pr_nino_jja_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map_americaS__", "pr_nino_jja_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map_asiaS__", "pr_nino_jja_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina PRA JJA", "El Nino PRA JJA"],
            "varpattern": ["pr_nina_jja_map_oceania__", "pr_nino_jja_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged PR\n4) Temporal mean PR removed\n" +
                      "5) PRA regridded to 1°x1°\n6) El Nino and La Nina JJA PRA composited\n7) Mask ocean",
        },
    },
    "EnsoPrTsRmse": {
        "netcdf_variables": ["pr_over_sst_ts__", "pr_over_sst_hov__", "Nina_pr_ts__", "Nino_pr_ts__", "Nina_pr_hov__",
                             "Nino_pr_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO's PRA life-cycle",
            #"varpattern": "pr_over_sst_ts__",
            "varpattern": "sst_against_pr_ts__",
            "xname": "months",
            "yname": "reg(ENSO SSTA, PRA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) Horizontal averaged\n5) Dec. N3.4 SSTA regressed onto REGION1 PRA\n\nMetric: RMSE$_{t}$",
        },
        "dive_down01": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "title": ["reg(ENSO SSTA, PRA)", "reg(ENSO SSTA, PRA)"],
            #"varpattern": "pr_over_sst_hov__",
            "varpattern": "sst_against_pr_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) PRA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional PRA average\n7) Dec. N3.4 SSTA regressed onto PRA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO's PRA life-cycle",
            "varpattern": ["Nina_pr_ts__", "Nino_pr_ts__"],
            "colors": {"model": ["blue", "red"], "reference": ["blue", "red"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["La Nina", "El Nino"],
            "xname": "months",
            "yname": "ENSO PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 PRA averaged\n6) El Nino and La Nina PRA composited",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG5"],
            "title": ["La Nina PRA", "El Nino PRA"],
            "varpattern": ["Nina_pr_hov__", "Nino_pr_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "PRA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) PRA regridded to 1°x1°\n6) 5S-5N meridional PRA average\n" +
                      "7) El Nino and La Nina PRA composited",
        },
    },
    "EnsoSlpMap": {
        "netcdf_variables": ["reg_slp_over_sst_map__", "reg_slp_over_sst_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Mask ocean",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA)", "reg(ENSO SSTA, SLPA)"],
            "varpattern": "reg_slp_over_sst_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SLPA regridded to 1°x1°\n" +
                      "6) Dec. N3.4 SSTA regressed onto Dec. SLPA\n7) Mask ocean",
        },
    },
    "EnsoSlpMapDjf": {
        "netcdf_variables": ["reg_slp_over_sst_djf_map__", "reg_slp_over_sst_djf_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            # "varpattern": "sst_over_sst_map__",
            "varpattern": "reg_slp_over_sst_djf_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map__", "slp_nino_djf_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            "varpattern": "reg_slp_over_sst_djf_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            "varpattern": "reg_slp_over_sst_djf_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            "varpattern": "reg_slp_over_sst_djf_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            "varpattern": "reg_slp_over_sst_djf_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) DJF", "reg(ENSO SSTA, SLPA) DJF"],
            "varpattern": "reg_slp_over_sst_djf_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF SLPA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map_africaSE__", "slp_nino_djf_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map_americaN__", "slp_nino_djf_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map_americaS__", "slp_nino_djf_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map_asiaS__", "slp_nino_djf_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA DJF", "El Nino SLPA DJF"],
            "varpattern": ["slp_nina_djf_map_oceania__", "slp_nino_djf_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina DJF SLPA composited\n7) Mask ocean",
        },
    },
    "EnsoSlpMapJja": {
        "netcdf_variables": ["reg_slp_over_sst_jja_map__", "reg_slp_over_sst_jja_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            # "varpattern": "sst_over_sst_map__",
            "varpattern": "reg_slp_over_sst_jja_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map__", "slp_nino_jja_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            "varpattern": "reg_slp_over_sst_jja_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            "varpattern": "reg_slp_over_sst_jja_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            "varpattern": "reg_slp_over_sst_jja_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            "varpattern": "reg_slp_over_sst_jja_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, SLPA) JJA", "reg(ENSO SSTA, SLPA) JJA"],
            "varpattern": "reg_slp_over_sst_jja_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA SLPA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map_africaSE__", "slp_nino_jja_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map_americaN__", "slp_nino_jja_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map_americaS__", "slp_nino_jja_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map_asiaS__", "slp_nino_jja_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina SLPA JJA", "El Nino SLPA JJA"],
            "varpattern": ["slp_nina_jja_map_oceania__", "slp_nino_jja_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged SLP\n4) Temporal mean SLP removed\n" +
                      "5) SLPA regridded to 1°x1°\n6) El Nino and La Nina JJA SLPA composited\n7) Mask ocean",
        },
    },
    "EnsoSstLonRmse": {
        "netcdf_variables": ["sst_over_sst_lon__", "sst_over_sst_map__", "Nina_sst_lon__", "Nino_sst_lon__",
                             "Nina_sst_map__", "Nino_sst_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO pattern",
            #"varpattern": "sst_over_sst_lon__",
            "varpattern": "sst_against_sst_lon__",
            "xname": "longitude",
            "yname": "reg(ENSO SSTA, SSTA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA averaged\n" +
                      "7) Dec. N3.4 SSTA regressed onto Dec. SSTA\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": True,
            "title": ["reg(ENSO SSTA, SSTA)", "reg(ENSO SSTA, SSTA)"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "sst_against_sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SSTA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. SSTA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO's SSTA pattern",
            "varpattern": ["Nina_sst_lon__", "Nino_sst_lon__"],
            "colors": {"model": ["blue", "red"], "reference": ["blue", "red"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["La Nina", "El Nino"],
            "xname": "longitude",
            "yname": "ENSO SSTA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA averaged\n" +
                      "7) El Nino and La Nina Dec. SSTA composited",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG25"],
            "maskland": True,
            "title": ["La Nina SSTA", "El Nino SSTA"],
            "varpattern": ["Nina_sst_map__", "Nino_sst_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) El Nino and La Nina Dec. SSTA composited",
        },
    },
    "EnsoSstMap": {
        "netcdf_variables": ["reg_ts_over_sst_map__", "reg_ts_over_sst_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Mask ocean",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA)", "reg(ENSO SSTA, TSA)"],
            "varpattern": "reg_ts_over_sst_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TSA regridded to 1°x1°\n6) Dec. N3.4 SSTA regressed onto Dec. TSA\n" +
                      "7) Mask ocean",
        },
    },
    "EnsoSstMapDjf": {
        "netcdf_variables": ["reg_ts_over_sst_djf_map__", "reg_ts_over_sst_djf_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            # "varpattern": "sst_over_sst_map__",
            "varpattern": "reg_ts_over_sst_djf_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map__", "ts_nino_djf_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        # ["africaSE", "americaN", "americaS", "asiaS", "oceania"]
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            "varpattern": "reg_ts_over_sst_djf_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            "varpattern": "reg_ts_over_sst_djf_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            "varpattern": "reg_ts_over_sst_djf_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            "varpattern": "reg_ts_over_sst_djf_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) DJF", "reg(ENSO SSTA, TSA) DJF"],
            "varpattern": "reg_ts_over_sst_djf_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) DJF averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) DJF N3.4 SSTA regressed onto DJF TSA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map_africaSE__", "ts_nino_djf_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map_americaN__", "ts_nino_djf_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map_americaS__", "ts_nino_djf_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map_asiaS__", "ts_nino_djf_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA DJF", "El Nino TSA DJF"],
            "varpattern": ["ts_nina_djf_map_oceania__", "ts_nino_djf_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) DJF averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina DJF TSA composited\n7) Mask ocean",
        },
    },
    "EnsoSstMapJja": {
        "netcdf_variables": ["reg_ts_over_sst_jja_map__", "reg_ts_over_sst_jja_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            # "varpattern": "sst_over_sst_map__",
            "varpattern": "reg_ts_over_sst_jja_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n" +
                      "7) Equatorial Pacific masked\n\nMetric: RMSE$_{xy}$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map__", "ts_nino_jja_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n" +
                      "7) Equatorial Pacific masked",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            "varpattern": "reg_ts_over_sst_jja_map_africaSE__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n7) Mask ocean",
        },
        "dive_down03": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            "varpattern": "reg_ts_over_sst_jja_map_americaN__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n7) Mask ocean",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            "varpattern": "reg_ts_over_sst_jja_map_americaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n7) Mask ocean",
        },
        "dive_down05": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            "varpattern": "reg_ts_over_sst_jja_map_asiaS__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n7) Mask ocean",
        },
        "dive_down06": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["reg(ENSO SSTA, TSA) JJA", "reg(ENSO SSTA, TSA) JJA"],
            "varpattern": "reg_ts_over_sst_jja_map_oceania__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "regression",
            "method": "1) Linearly detrended\n2) N3.4 SST averaged\n3) JJA averaged\n4) Temporal mean removed\n" +
                      "5) TSA regridded to 1°x1°\n6) JJA N3.4 SSTA regressed onto JJA TSA\n7) Mask ocean",
        },
        "dive_down07": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map_africaSE__", "ts_nino_jja_map_africaSE__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n7) Mask ocean",
        },
        "dive_down08": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map_americaN__", "ts_nino_jja_map_americaN__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n7) Mask ocean",
        },
        "dive_down09": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map_americaS__", "ts_nino_jja_map_americaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n7) Mask ocean",
        },
        "dive_down10": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map_asiaS__", "ts_nino_jja_map_asiaS__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n7) Mask ocean",
        },
        "dive_down11": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["PRA"],
            "maskland": False,
            "maskocean": True,
            "title": ["La Nina TSA JJA", "El Nino TSA JJA"],
            "varpattern": ["ts_nina_jja_map_oceania__", "ts_nino_jja_map_oceania__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Linearly detrended\n3) JJA averaged TS\n4) Temporal mean TS removed\n" +
                      "5) TSA regridded to 1°x1°\n6) El Nino and La Nina JJA TSA composited\n7) Mask ocean",
        },
    },
    "EnsoSstTsRmse": {
        "netcdf_variables": ["sst_over_sst_ts__", "sst_over_sst_hov__", "Nina_sst_ts__", "Nino_sst_ts__",
                             "Nina_sst_hov__", "Nino_sst_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO life-cycle",
            #"varpattern": "sst_over_sst_ts__",
            "varpattern": "sst_against_sst_ts__",
            "xname": "months",
            "yname": "reg(ENSO SSTA, SSTA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) Horizontal averaged\n5) Dec. N3.4 SSTA regressed onto REGION1 SSTA\n\nMetric: RMSE$_{t}$",
        },
        "dive_down01": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG12"],
            "title": ["reg(ENSO SSTA, SSTA)", "reg(ENSO SSTA, SSTA)"],
            #"varpattern": "sst_over_sst_hov__",
            "varpattern": "sst_against_sst_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) SSTA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional SSTA average\n7) Dec. N3.4 SSTA regressed onto Dec. SSTA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO's SSTA life-cycle",
            "varpattern": ["Nina_sst_ts__", "Nino_sst_ts__"],
            "colors": {"model": ["blue", "red"], "reference": ["blue", "red"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["La Nina", "El Nino"],
            "xname": "months",
            "yname": "ENSO SSTA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 SSTA averaged\n6) El Nino and La Nina SSTA composited",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "title": ["La Nina SSTA", "El Nino SSTA"],
            "varpattern": ["Nina_sst_hov__", "Nino_sst_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "SSTA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA average\n" +
                      "7) El Nino and La Nina SSTA composited",
        },
    },
    "EnsoTauxTsRmse": {
        "netcdf_variables": ["taux_over_sst_ts__", "taux_over_sst_hov__", "Nina_taux_ts__", "Nino_taux_ts__",
                             "Nina_taux_hov__", "Nino_taux_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO life-cycle",
            #"varpattern": "taux_over_sst_ts__",
            "varpattern": "sst_against_taux_ts__",
            "xname": "months",
            "yname": "reg(ENSO SSTA, TAUXA)",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) Horizontal averaged\n5) Dec. N3.4 SSTA regressed onto REGION1 TAUXA\n\nMetric: RMSE$_{t}$",
        },
        "dive_down01": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG20"],
            "title": ["reg(ENSO SSTA, TAUXA)", "reg(ENSO SSTA, TAUXA)"],
            #"varpattern": "taux_over_sst_hov__",
            "varpattern": "sst_against_taux_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "regression",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) 5-month triangular running ave.\n" +
                      "4) N3.4 SSTA averaged\n5) TAUXA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional TAUXA average\n7) Dec. N3.4 SSTA regressed onto TAUXA",
        },
        "dive_down02": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "ENSO's TAUXA life-cycle",
            "varpattern": ["Nina_taux_ts__", "Nino_taux_ts__"],
            "colors": {"model": ["blue", "red"], "reference": ["blue", "red"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["La Nina", "El Nino"],
            "xname": "months",
            "yname": "ENSO TAUXA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 TAUXA averaged\n6) El Nino and La Nina TAUXA composited",
        },
        "dive_down03": {
            "plot_type": "hovmoeller",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG30"],
            "title": ["La Nina TAUXA", "El Nino TAUXA"],
            "varpattern": ["Nina_taux_hov__", "Nino_taux_hov__"],
            "xname": "longitude",
            "yname": "months",
            "zname": "TAUXA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) TAUXA regridded to 1°x1°\n6) 5S-5N meridional TAUXA average\n" +
                      "7) El Nino and La Nina TAUXA composited",
        },
    },
    "EnsoSeasonality": {
        "netcdf_variables": ["sstStd_monthly__", "sstStd_hov__", "sstStd_NDJ_lon__", "sstStd_MAM_lon__",
                             "sstStd_NDJ_map__", "sstStd_MAM_map__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO Seasonality",
            "varpattern": "diagnostic",
            "yname": "SSTA std (NDJ/MAM)",
            "method": "1) Linearly detrended\n2) REGION1 averaged\n3) NDJ and MAM averaged\n" +
                      "4) Temporal mean removed\n5) Standard deviation\n6) ratio = STD$_{NDJ}$/STD$_{MAM}$\n\n" +
                      "Metric: abs((Ratio$_{mod}$-Ratio$_{ref}$)/Ratio$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSTA standard deviation",
            "varpattern": "sstStd_monthly__",
            "xname": "months",
            "yname": "SSTA std",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n" +
                      "4) Standard deviation for each calendar month",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude"],
            "title": ["SSTA standard deviation", "SSTA standard deviation"],
            "varpattern": "sstStd_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "SSTA std",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n" +
                      "3) Standard deviation for each calendar month\n4) SSTA regridded to 1°x1°\n" +
                      "5) 5S-5N meridional averaged",
        },
        "dive_down03": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSTA standard deviation",
            "varpattern": ["sstStd_NDJ_lon__", "sstStd_MAM_lon__"],
            "colors": {"model": ["red", "blue"], "reference": ["red", "blue"]},
            "linestyles": {"model": ["-", "-"], "reference": ["-.", "-."]},
            "legend": ["NDJ", "MAM"],
            "xname": "longitude",
            "yname": "SSTA std",
            "method": "1) Linearly detrended\n2) REGION1 averaged\n3) NDJ and MAM averaged\n" +
                      "4) Temporal mean removed\n5) Standard deviation\n6) Regridded to 1°x1°\n" +
                      "7) 5S-5N meridional averaged",
        },
        "dive_down04": {
            "plot_type": "map",
            "nbr_panel": 4,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["NDJ", "MAM"], #["SSTA std NDJ", "SSTA std MAM"],
            "varpattern": ["sstStd_NDJ_map__", "sstStd_MAM_map__"],
            "xname": "longitude",
            "yname": "latitude",
            "zname": "monthly SSTA std",
            "method": "1) Linearly detrended\n2) REGION1 averaged\n3) NDJ and MAM averaged\n" +
                      "4) Temporal mean removed\n5) Standard deviation\n6) Regridded to 1°x1°",
        },
    },
    "EnsoSstDiversity": {
        "netcdf_variables": ["Enso_lon_pos_maxSSTA__", "Nina_lon_pos_minSSTA__", "Nino_lon_pos_maxSSTA__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO diversity",
            "varpattern": "diagnostic",
            "yname": "IQR of min/max SSTA",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA averaged\n7) 5° triangular running ave.\n" +
                      "8) Find zonal location of El Nino max(SSTA) and La Nina min(SSTA)\n" +
                      "7) IQR El Nino max(SSTA) with La Nina min(SSTA)\n\n" +
                      "Metric: abs((IQR$_{mod}$-IQR$_{ref}$)/IQR$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "boxplot",
            "nbr_panel": 3,
            "title": ["ENSO diversity", "La Nina diversity", "El Nino diversity"],
            "varpattern": ["Enso_lon_pos_maxSSTA__", "Nina_lon_pos_minSSTA__", "Nino_lon_pos_maxSSTA__"],
            "yname": ["longitude of min/max SSTA", "longitude of min SSTA", "longitude of max SSTA"],
            "custom_label": "longitude",
            "method": "1) Detect El Nino and La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA averaged\n7) 5° triangular running ave.\n" +
                      "8) Find zonal location of El Nino max(SSTA) and La Nina min(SSTA)\n" +
                      "7) IQR El Nino max(SSTA) with La Nina min(SSTA)\nand IQR El Nino and IQR La Nina\n\n" +
                      "Metric: abs((IQR$_{mod}$-IQR$_{ref}$)/IQR$_{ref}$)*100",
        },
    },
    "EnsoSstSkew": {
        "netcdf_variables": ["sstSke_lon__", "sstSke_map__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "ENSO skewness",
            "varpattern": "diagnostic",
            "yname": "SSTA skewness",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) REGION1 averaged\n4) Skewness\n\n" +
                      "Metric: abs((SKE$_{mod}$-SKE$_{ref}$)/SKE$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSTA skewness",
            "varpattern": "sstSke_lon__",
            "xname": "longitude",
            "yname": "SSTA skew",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Skewness\n4) Regridded to 1°x1°\n" +
                      "5) 5S-5N meridional averaged",
        },
        "dive_down02": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": True,
            "title": ["SSTA skew", "SSTA skew"],
            "varpattern": "sstSke_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA skew",
            "method": "1) Seasonal cycle removed\n2) Linearly detrended\n3) Skewness\n4) Regridded to 1°x1°",
        },
    },
    "NinaPrMap": {
        "netcdf_variables": ["prComp_map__", "prComp_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG2"],
            "maskland": False,
            "title": ["La Nina composite", "La Nina composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "prComp_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) PRA regridded to 1°x1°6) La Nina Dec. PRA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinaSlpMap": {
        "netcdf_variables": ["slp_map__", "slp_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "title": ["La Nina composite", "La Nina composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "slp_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SLPA regridded to 1°x1°6) La Nina Dec. SLPA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinaSstMap": {
        "netcdf_variables": ["ts_map__", "ts_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "title": ["La Nina composite", "La Nina composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "ts_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TSA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) TSA regridded to 1°x1°6) La Nina Dec. TSA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinaSstDur": {
        "netcdf_variables": ["Nina_duration__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "La Nina duration",
            "varpattern": "diagnostic",
            "yname": "duration (SSTA<-0.5)",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 averaged\n6) Duration = nbr months < -0.5 ENSO STD\n7) Mean La Nina duration\n\n" +
                      "Metric: abs((DUR$_{mod}$-DUR$_{ref}$)/DUR$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "boxplot",
            "nbr_panel": 1,
            "title": "La Nina duration",
            "varpattern": "Nina_duration__",
            "yname": "duration (SSTA<-0.5)",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 averaged\n6) Duration = nbr months < -0.5 ENSO STD",
        },
    },
    "NinaSstLonRmse": {
        "netcdf_variables": ["sst_lon__", "sst_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "La Nina pattern",
            "varpattern": "sst_lon__",
            "xname": "longitude",
            "yname": "SSTA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA averaged\n" +
                      "7) La Nina Dec. SSTA composited\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["dSST"],
            "maskland": True,
            "title": "La Nina SSTA",
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) La Nina Dec. SSTA composited",
        },
    },
    "NinaSstTsRmse": {
        "netcdf_variables": ["sst_ts__", "sst_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "La Nina life-cycle",
            "varpattern": "sst_ts__",
            "xname": "months",
            "yname": "SSTA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) REGION1 SSTA averaged\n6) La Nina SSTA composited\n\nMetric: RMSE$_{t}$",
        },
        "dive_down01": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "title": "La Nina SSTA",
            "varpattern": "sst_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "SSTA",
            "method": "1) Detect La Nina\n    (5-m. tri. ave. Dec. N3.4 SSTA < -0.75 STD)\n" +
                      "2) Seasonal cycle removed\n3) Linearly detrended\n4) 5-month triangular running ave.\n" +
                      "5) SSTA regridded to 1°x1°\n6) 5S-5N meridional SSTA average\n7) La Nina SSTA composited",
        },
    },
    "NinoPrMap": {
        "netcdf_variables": ["pr_map__", "pr_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "title": ["El Nino composite", "El Nino composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PRA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) PRA regridded to 1°x1°" +
                      "6) El Nino Dec. PRA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinoSlpMap": {
        "netcdf_variables": ["slp_map__", "slp_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG3"],
            "maskland": False,
            "title": ["El Nino composite", "El Nino composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "slp_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SLPA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SLPA regridded to 1°x1°" +
                      "6) El Nino Dec. SLPA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinoSstMap": {
        "netcdf_variables": ["ts_map__", "ts_map__"],
        "diagnostic": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "maskland": False,
            "title": ["El Nino composite", "El Nino composite"],
            #"varpattern": "sst_over_sst_map__",
            "varpattern": "ts_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°" +
                      "6) El Nino Dec. SSTA composited\n\nMetric: RMSE$_{xy}$",
        },
    },
    "NinoSstDiversity": {
        "netcdf_variables": ["Nina_lon_pos_minSSTA__", "Nino_lon_pos_maxSSTA__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "El Nino diversity",
            "varpattern": "diagnostic",
            "yname": "IQR of max SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional SSTA averaged\n7) 5° triangular running ave.\n" +
                      "8) Find zonal location of max(SSTA)\n7) IQR El Nino max(SSTA)\n\n" +
                      "Metric: abs((IQR$_{mod}$-IQR$_{ref}$)/IQR$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "boxplot",
            "nbr_panel": 2,
            "title": ["La Nina diversity", "El Nino diversity"],
            "varpattern": ["Nina_lon_pos_minSSTA__", "Nino_lon_pos_maxSSTA__"],
            "yname": ["longitude of min SSTA", "longitude of max SSTA"],
            "custom_label": "longitude",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional SSTA averaged\n7) 5° triangular running ave.\n" +
                      "8) Find zonal location of max(SSTA)",
        },
    },
    "NinoSstDur": {
        "netcdf_variables": ["Nino_duration__"],
        "diagnostic": {
            "plot_type": "dot",
            "nbr_panel": 1,
            "title": "El Nino duration",
            "varpattern": "diagnostic",
            "yname": "duration (SSTA>0.5)",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) REGION1 averaged\n" +
                      "6) Duration = nbr months > 0.5 ENSO STD\n7) Mean El Nino duration\n\n" +
                      "Metric: abs((DUR$_{mod}$-DUR$_{ref}$)/DUR$_{ref}$)*100",
        },
        "dive_down01": {
            "plot_type": "boxplot",
            "nbr_panel": 1,
            "title": "El Nino duration",
            "varpattern": "Nino_duration__",
            "yname": "duration (SSTA>0.5)",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) REGION1 averaged\n" +
                      "6) Duration = nbr months > 0.5 ENSO STD",
        },
    },
    "NinoSstLonRmse": {
        "netcdf_variables": ["sst_lon__", "sst_map__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "El Nino pattern",
            "varpattern": "sst_lon__",
            "xname": "longitude",
            "yname": "SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional SSTA averaged\n7) El Nino Dec. SSTA composited\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["dSST"],
            "maskland": True,
            "title": "El Nino SSTA",
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°\n" +
                      "6) El Nino Dec. SSTA composited",
        },
    },
    "NinoSstTsRmse": {
        "netcdf_variables": ["sst_ts__", "sst_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "El Nino life-cycle",
            "varpattern": "sst_ts__",
            "xname": "months",
            "yname": "SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) REGION1 SSTA averaged\n" +
                      "6) El Nino SSTA composited\n\nMetric: RMSE$_{t}$",
        },
        "dive_down01": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["SKEW"],
            "title": "El Nino SSTA",
            "varpattern": "sst_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "SSTA",
            "method": "1) Detect El Nino\n    (5-m. tri. ave. Dec. N3.4 SSTA > 0.75 STD)\n2) Seasonal cycle removed\n" +
                      "3) Linearly detrended\n4) 5-month triangular running ave.\n5) SSTA regridded to 1°x1°\n" +
                      "6) 5S-5N meridional SSTA average\n7) El Nino SSTA composited",
        },
    },
    "SeasonalPrLatRmse": {
        "netcdf_variables": ["pr_lat__", "pr_map__", "prMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "PR seasonal cycle std",
            "varpattern": "pr_lat__",
            "xname": "latitude",
            "yname": "PR std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude5"],
            "maskland": True,
            "title": ["PR seasonal cycle std"],
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PR std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["PR"],
            "label": dict_label["amplitude15"],
            "title": ["PR seasonal cycle", "PR seasonal cycle"],
            "varpattern": "prMac_hov__",
            "xname": "latitude",
            "yname": "months",
            "zname": "PR",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Zonal averaged (see box)",
        },
    },
    "SeasonalPrLonRmse": {
        "netcdf_variables": ["pr_lon__", "pr_map__", "prMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "PR seasonal cycle std",
            "varpattern": "pr_lon__",
            "xname": "longitude",
            "yname": "PR std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude5"],
            "maskland": True,
            "title": ["PR seasonal cycle std"],
            "varpattern": "pr_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "PR std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["PR"],
            "label": dict_label["amplitude10"],
            "title": ["PR seasonal cycle", "PR seasonal cycle"],
            "varpattern": "prMac_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "PR",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Meridional averaged (see box)",
        },
    },
    "SeasonalSshLatRmse": {
        "netcdf_variables": ["ssh_lat__", "ssh_map__", "sshMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSH seasonal cycle std",
            "varpattern": "ssh_lat__",
            "xname": "latitude",
            "yname": "SSH std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["SSH seasonal cycle std"],
            "varpattern": "ssh_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSH std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["SST"],
            "title": ["SSH seasonal cycle", "SSH seasonal cycle"],
            "varpattern": "sshMac_hov__",
            "xname": "latitude",
            "yname": "months",
            "zname": "SSH",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Zonal averaged (see box)",
        },
    },
    "SeasonalSshLonRmse": {
        "netcdf_variables": ["ssh_lon__", "ssh_map__", "sshMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SSH seasonal cycle std",
            "varpattern": "ssh_lon__",
            "xname": "longitude",
            "yname": "SSH std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["SSH seasonal cycle std"],
            "varpattern": "ssh_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SSH std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],  # YYP: I do not know yet the colobar / label needed
            "label": dict_label["SST"],
            "title": ["SSH seasonal cycle", "SSH seasonal cycle"],
            "varpattern": "sshMac_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "SSH",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Meridional averaged (see box)",
        },
    },
    "SeasonalSstLatRmse": {
        "netcdf_variables": ["sst_lat__", "sst_map__", "sstMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SST seasonal cycle std",
            "varpattern": "sst_lat__",
            "xname": "latitude",
            "yname": "SST std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["SST seasonal cycle std"],
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SST std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],
            "label": dict_label["SST"],
            "title": ["SST seasonal cycle", "SST seasonal cycle"],
            "varpattern": "sstMac_hov__",
            "xname": "latitude",
            "yname": "months",
            "zname": "SST",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Zonal averaged (see box)",
        },
    },
    "SeasonalSstLonRmse": {
        "netcdf_variables": ["sst_lon__", "sst_map__", "sstMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "SST seasonal cycle std",
            "varpattern": "sst_lon__",
            "xname": "longitude",
            "yname": "SST std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude"],
            "maskland": True,
            "title": ["SST seasonal cycle std"],
            "varpattern": "sst_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "SST std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["SST"],
            "label": dict_label["SST"],
            "title": ["SST seasonal cycle", "SST seasonal cycle"],
            "varpattern": "sstMac_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "SST",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Meridional averaged (see box)",
        },
    },
    "SeasonalTauxLatRmse": {
        "netcdf_variables": ["taux_lat__", "taux_map__", "tauxMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "TAUX seasonal cycle std",
            "varpattern": "taux_lat__",
            "xname": "latitude",
            "yname": "TAUX std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Zonal averaged (see box)\n\nMetric: RMSE$_y$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude60"],
            "maskland": True,
            "title": ["TAUX seasonal cycle std"],
            "varpattern": "taux_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TAUX std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["TAUX"],
            "title": ["TAUX seasonal cycle", "TAUX seasonal cycle"],
            "varpattern": "tauxMac_hov__",
            "xname": "latitude",
            "yname": "months",
            "zname": "TAUX",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Zonal averaged (see box)",

        },
    },
    "SeasonalTauxLonRmse": {
        "netcdf_variables": ["taux_lon__", "taux_map__", "tauxMac_hov__"],
        "diagnostic": {
            "plot_type": "curve",
            "nbr_panel": 1,
            "title": "TAUX seasonal cycle std",
            "varpattern": "taux_lon__",
            "xname": "longitude",
            "yname": "TAUX std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°\n5) Meridional averaged (see box)\n\nMetric: RMSE$_x$",
        },
        "dive_down01": {
            "plot_type": "map",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["amplitude"],
            "label": dict_label["amplitude60"],
            "maskland": True,
            "title": ["TAUX seasonal cycle std"],
            "varpattern": "taux_map__",
            "xname": "longitude",
            "yname": "latitude",
            "zname": "TAUX std",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Standard deviation\n" +
                      "4) Regridded to 1°x1°",
        },
        "dive_down02": {
            "plot_type": "hovmoeller",
            "nbr_panel": 2,
            "colorbar": dict_colorbar["anomalies"],
            "label": dict_label["REG80"],
            "title": ["TAUX seasonal cycle", "TAUX seasonal cycle"],
            "varpattern": "tauxMac_hov__",
            "xname": "longitude",
            "yname": "months",
            "zname": "TAUX",
            "method": "1) Linearly detrended\n2) Seasonal cycle computed\n3) Regridded to 1°x1°\n" +
                      "5) Meridional averaged (see box)",
        },
    },

}


reference_observations = {
    "ssh": "AVISO", "pr": "GPCPv2.3", "sst": "Tropflux", "lhf": "Tropflux", "lwr": "Tropflux", "shf": "Tropflux",
    "slp": "ERA-Interim", "swr": "Tropflux", "taux": "Tropflux", "thf": "Tropflux"
}


def plot_param(metric_collection, metric):
    dict_MC = defCollection(metric_collection)
    dict_MCm = dict_MC["metrics_list"][metric]
    # get plot parameters
    dict_out = plot_parameters[metric.replace("_1", "").replace("_2", "").replace("_3", "")]
    # get metric computation
    if metric_collection in ["ENSO_tel", "test_tel"] and "Map" in metric:
        computation = ["CORR", "RMSE"]
    elif "rmse" in metric.lower():
        computation = "RMSE"
    elif "corr" in metric.lower():
        computation = "CORR"
    else:
        try:
            computation = dict_MCm["metric_computation"]
        except:
            try:
                computation = dict_MC["common_collection_parameters"]["metric_computation"]
            except:
                computation = default_arg_values("metric_computation")
    dict_out["metric_computation"] = computation
    # get metric variables
    variables = dict_MCm["variables"]
    dict_out["metric_variables"] = variables
    # get metric reference
    references = dict((var, reference_observations[var]) for var in variables)
    if "sst" in variables and len(variables) > 1:
        references["sst"] = "Tropflux"
    refname = references[variables[0]]
    for var in variables[1:]:
        refname = refname + "_" + references[var]
    # if metric_collection == "ENSO_tel" and metric in ["EnsoPrMap", "EnsoSlpMap"]:
    #     refname = "ERA-Interim_ERA-Interim"
    # elif metric_collection == "ENSO_tel" and metric in ["EnsoSstMap", "NinaSstMap", "NinoSstMap"]:
    #     refname = "ERA-Interim"
    if metric_collection in ["ENSO_tel", "test_tel"] and metric in\
            ["EnsoSstMap", "NinaSstMap", "NinoSstMap", "EnsoSstMapDjf", "NinoSstMapDjf", "NinaSstMapDjf",
             "EnsoSstMapJja", "NinoSstMapJja", "NinaSstMapJja"]:
        refname = "ERA-Interim"
    dict_out["metric_reference"] = refname
    # get variable regions
    dict_out["metric_regions"] = dict_MCm["regions"]
    return dict_out
