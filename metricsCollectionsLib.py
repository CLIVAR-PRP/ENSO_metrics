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
        'EnsoMu':{'nbvar':2,'var_names':['sst','taux'],'ref_obs':['HadiSST1.1','ERA-interim']},
    }
    if VOR:
        return var_obs_requirements
    else:
        return var_obs_requirements[VOR]

