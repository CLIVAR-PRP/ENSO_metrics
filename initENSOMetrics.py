# Initialisation of ENSO metrics parameters
#
# name of variables in model data files 'internal_name','var_name in file'
# currently set up as standard CMIP short_name
def vars_names:
    var_names = {
    'sst':'tos',
    'taux':'tauuo',
    #...
    }
    return var_names

#
# list of reference observation (var name in file and file name)
# TODO: how to deal with several files ? use wild card * ? use obs4mips code/link ?
# TODO: how to deal with references which are not in obs4mips ?
def ref_obs(RO=True):
    reference_obs_file = {
        'HadiSST1.1':{'sst':{'var_name':'tos','file_name':'<file_name_sst>','source':'see <ENSO RF web site>'},'sie':{'var_name':'seaIce','file_name':'<file_name_sie>','source':'other'}},
        'ERA-interim':{'taux':{'var_name':'tauu','file_name':'<file_name_tauu>','source':'Obs4MIPS:<link>'}},
    }
    if RO:
        return reference_obs_file
    else:
        return reference_obs_file[RO]

#
# List of averaging regions
#
# TODO: how to deal with longitude bounds ?
# TODO See enso_bellenger/lib/PMP_rectangular_domains.py
# TODO add argument for periodicty or initial bounds ?
def averageRegion(AR=True):
    ave_region= {
        'nino3':{'long_name':'Niño 3','latitude':[-5.,5.],'longitude':[-150,-90]},
        'nino4':{'long_name':'Niño 4','latitude':[-5.,5.],'longitude':[160,210]},
    }
    if AR:
        return ave_region
    else:
        return ave_region[AR]
