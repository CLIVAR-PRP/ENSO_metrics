#
# Define ENSO metrics collections as a function of science question/realm
#
# Draft version
#
# Define metrics collections
def defCollection(MC=True):
# Name, list of metrics
    metrics_collection = {
        'MC1':{'long_name':'Metrics collection Q1','list_of_metrics':['EnsoAmpl','EnsoMu'],
               'description':'Describe which science question this collection is about'},
}
    if MC:
        return metrics_collection
    else:
        return metrics_collection[MC]

# List of metrics requirements (var name and reference obs)
def metricReqs(VOR=True):
    var_obs_requirements = {
        'EnsoAmpl':{'nbvar':1,'var_names':['sst'],'ref_obs':['HadiSST1.1']},
        'EnsoMu':{'nbvar':2,'var_names':['sst','taux'],'ref_obs':['HadiSST1.1','ERAint']},
    }
    if VOR:
        return var_obs_requirements
    else:
        return var_obs_requirements[VOR]

# list of reference observation (var name in file and file name)
# TODO: how to deal with several files ? use wild card * ? use obs4mips code/link ?
def ref_obs(RO=True):
    reference_obs_file = {
        'HadiSST1.1':{'sst':{'var_name':'tos','file_name':'<file_name_sst>'},'sie':{'var_name':'seaIce','file_name':'<file_name_sie>'}},
        'ERAint':{'taux':{'var_name':'tauu','file_name':'<file_name_tauu>'}},
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
