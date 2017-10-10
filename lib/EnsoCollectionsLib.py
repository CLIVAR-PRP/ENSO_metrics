#
# Define ENSO metrics collections as a function of science question/realm
#
# Draft version
#





# Define metrics collections
def defCollection(MC=True):
# Name, list of metrics
    metrics_collection = {
        'MC1':{'long_name':'Metrics Collection 1',\
               'list_of_metrics':['EnsoAmpl','EnsoSeasonality','EnsoRMSE','EnsoMu','EnsoAlphaSwr'],\
		       'dict_of_regions':{'EnsoAmpl':{'sst':'nino3'},'EnsoSeasonality':{'sst':'nino3'},\
                                  'EnsoRMSE':{'sst':'tropical_pacific'},'EnsoMu':{'sst':'nino3','taux':'nino4'},\
                                  'EnsoAlphaSwr':{'sst':'nino3','swr':'nino3'}},\
               'description':'Describe which science question this collection is about'},\
    }
    if MC:
        return metrics_collection
    else:
        return metrics_collection[MC]


# List of metrics requirements (var name and reference obs)
def metricReqs(VOR=True):
    var_obs_requirements = {
		'EnsoAlphaSwr':   {'nbvar':2,'var_names':['sst','swr'], 'ref_obs':['Tropflux','Tropflux']},\
        'EnsoAmpl':       {'nbvar':1,'var_names':['sst'],       'ref_obs':['OISST']},\
        'EnsoMu':         {'nbvar':2,'var_names':['sst','taux'],'ref_obs':['Tropflux','Tropflux']},\
		'EnsoRMSE':       {'nbvar':1,'var_names':['sst'],       'ref_obs':['OISST']},\
		'EnsoSeasonality':{'nbvar':1,'var_names':['sst'],       'ref_obs':['OISST']},\
    }
    if VOR:
        return var_obs_requirements
    else:
        return var_obs_requirements[VOR]


# List of reference observations for each variables
def ReferenceObservations(VAR=True, DATASET=True):
    dict_reference_observations = {
		'lhf': {'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/lhf/',\
		                       'file_name':'lhf_tropflux_1m_*.nc',\
		                       'var_name_in_file':'lhf'},\
		        },\
		'lwr': {'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/lwr/',\
		                       'file_name':'lwr_tropflux_1m_*.nc',\
		                       'var_name_in_file':'lwr'},\
		        },\
		'shf': {'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/shf/',\
		                       'file_name':'shf_tropflux_1m_*.nc',\
		                       'var_name_in_file':'shf'},\
		        },\
		'sst': {'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'HadISST1.1': {'source':'see https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html',\
		                       'file_name':'HadISST_sst.nc',\
		                       'var_name_in_file':'sst'},\
		        'OISST':      {'source':'see https://www.earthsystemcog.org/search/obs4mips/?template=obs4mips&limit=200',\
		                       'file_name':'tos_OISST_L4_AVHRR-only-v2_*-*.nc',\
		                       'var_name_in_file':'tos'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/sst/',\
		                       'file_name':'sst_tropflux_1m_*.nc',\
		                       'var_name_in_file':'tos'},\
		        },\
		'swr': {'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/swr/',\
		                       'file_name':'swr_tropflux_1m_*.nc',\
		                       'var_name_in_file':'swr'},\
		        },\
		'taux':{'ERA-Interim':{'source':'see http://apps.ecmwf.int/datasets/data/interim-full-moda/levtype=sfc/',\
		                       'file_name':'???',\
		                       'var_name_in_file':'???'},\
		        'Tropflux':   {'source':'see http://www.incois.gov.in/tropflux_datasets/data/monthly/swr/',\
		                       'file_name':'swr_tropflux_1m_*.nc',\
		                       'var_name_in_file':'swr'},\
		        },\
    }
    if VAR:
        return dict_reference_observations
    else:
        if DATASET: return var_obs_requirements[VAR]
        else:       return var_obs_requirements[VAR][DATASET]


def ReferenceRegions(AR=True):
    dict_reference_regions = {
        'tropical_pacific':          {'long_name':'Tropical Pacific (TP)',           'latitude':(-30.,30.),'longitude':(120.,280.)},\
        'equatorial_pacific':        {'long_name':'Equatorial Pacific (EP)',         'latitude':( -5., 5.),'longitude':(120.,280.)},\
		'eastern_equatorial_pacific':{'long_name':'Eastern Equatorial Pacific (EEP)','latitude':( -5., 5.),'longitude':(120.,205.)},\
		'western_equatorial_pacific':{'long_name':'Western Equatorial Pacific (WEP)','latitude':( -5., 5.),'longitude':(205.,280.)},\
		'nino1+2':{'long_name':'Niño 1+2','latitude':(-10., 0.),'longitude':(270.,280.)},\
        'nino3':  {'long_name':'Niño 3',  'latitude':( -5., 5.),'longitude':(210.,270.)},\
		'nino3.4':{'long_name':'Niño 3.4','latitude':( -5., 5.),'longitude':(190.,240.)},\
        'nino4':  {'long_name':'Niño 4',  'latitude':( -5., 5.),'longitude':(160.,210.)},\
    }
    if AR:
        return dict_reference_regions
    else:
        return dict_reference_regions[AR]
