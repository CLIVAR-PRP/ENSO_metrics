# -*- coding:UTF-8 -*-
from inspect import stack as INSPECTstack

# ENSO_metrics package functions:
from EnsoErrorsWarnings import UnknownKeyArg


# ---------------------------------------------------------------------------------------------------------------------#
#
# Library to ENSO metrics arguments (arg parser)
# These functions analyses given arguments and sets some arguments to their default value
#
def DefaultArgValues(arg):
    default = {
        'detrending': False, 'frequency': None, 'metric_computation': 'difference', 'min_time_steps': None,
        'normalization': False, 'project_interpreter': 'CMIP', 'regridding': False, 'smoothing': False,
        'treshold_ep_ev': -140, 'time_bounds': None, 'time_bounds_mod': None, 'time_bounds_obs': None,
    }
    try:
        default[arg]
    except:
        UnknownKeyArg(arg, INSPECTstack())
    return default[arg]
# ---------------------------------------------------------------------------------------------------------------------#