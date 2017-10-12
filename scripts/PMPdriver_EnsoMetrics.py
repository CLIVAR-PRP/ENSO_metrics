#!/usr/bin/env python
#=================================================
# Dependencies
#-------------------------------------------------
import cdms2
import collections
import copy
import json
import os
import pcmdi_metrics
import sys

from collections import defaultdict
from pcmdi_metrics.pcmdi.pmp_parser import PMPParser

from EnsoMetrics.EnsoCollectionsLib import defCollection 
from EnsoMetrics.EnsoMetricsLib import EnsoAlphaLhf, EnsoAlphaLwr, EnsoAlphaSwr, EnsoAlphaThf, EnsoAmpl, EnsoMu, EnsoRMSE, EnsoSeasonality
from EnsoMetrics.pmpParser import ReadOptions

#=================================================
# Collect user defined options
#-------------------------------------------------
param = ReadOptions()

#=================================================
# User input option
#-------------------------------------------------
mcdict = defCollection(param.metricsCollection)
metrics = mcdict['list_of_metrics']

#=================================================
# Prepare loop iteration
#-------------------------------------------------
# Setup where to output resulting ---
try:
    os.mkdir(param.outpathjson)
except BaseException:
    pass

# Insert observation at the beginning of the loop ---
models = copy.copy(param.modnames)
models.insert(0,'obs')

# Dictionary to save result ---
def tree(): return defaultdict(tree)
enso_stat_dic = tree() # Use tree dictionary to avoid declearing everytime

#=================================================
# Loop for Observation and Models 
#-------------------------------------------------
for mod in models:
    print ' ----- ', mod,' ---------------------'
  
    if mod == 'obs':
        sstName = param.sstNameObs
        tauxName = param.tauxNameObs

        sstFile = param.sstObsPath 
        tauxFile = param.tauxObsPath
    else:
        sstName = param.sstName
        tauxName = param.tauxName

        sstFile = (param.modpath.replace('MOD', mod)).replace('VAR',sstName) ## Will need land mask out at some point...!
        tauxFile = (param.modpath.replace('MOD', mod)).replace('VAR',tauxName)

    print sstFile
    print tauxFile
  
    try:
        for metric in metrics:

            print metric

            if metric == 'EnsoAmpl':
                nBox = mcdict['dict_of_regions'][metric]['sst']
                tmp_dict = EnsoAmpl(sstFile, sstName, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoSeasonality':
                nBox = mcdict['dict_of_regions'][metric]['sst']
                tmp_dict = EnsoSeasonality(sstFile, sstName, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoMu':
                sstBox = mcdict['dict_of_regions'][metric]['sst']
                tauxBox = mcdict['dict_of_regions'][metric]['taux']
                tmp_dict = EnsoMu(sstFile, tauxFile, sstName, tauxName, sstBox, tauxBox)
                tmp_dict['input_data'] = [sstFile, tauxFile]
                enso_stat_dic[mod][metric] = tmp_dict

            elif metric == 'EnsoRMSE' and mod != 'obs':
                nBox = mcdict['dict_of_regions'][metric]['sst']
                tmp_dict = EnsoRMSE(sstFile, sstName, param.sstObsPath, param.sstNameObs, nBox)
                tmp_dict['input_data'] = [sstFile]
                enso_stat_dic[mod][metric] = tmp_dict

    except:
        print 'failed for ', mod
  
#=================================================
#  OUTPUT METRICS TO JSON FILE
#-------------------------------------------------
OUT = pcmdi_metrics.io.base.Base(os.path.abspath(param.outpathjson), param.outnamejson)

disclaimer = open(
    os.path.join(
        sys.prefix,
        "share",
        "pmp",
        "disclaimer.txt")).read()

metrics_dictionary = collections.OrderedDict()
metrics_dictionary["DISCLAIMER"] = disclaimer
metrics_dictionary["REFERENCE"] = "The statistics in this file are based on Bellenger, H et al. Clim Dyn (2014) 42:1999-2018. doi:10.1007/s00382-013-1783-z"
metrics_dictionary["RESULTS"] = enso_stat_dic

OUT.write(
    metrics_dictionary,
    json_structure=["model", "metric", "item", "value or description"],
    indent=4,
    separators=(
        ',',
        ': '),
    sort_keys=True)

sys.exit('done')
