#!/usr/bin/env python

import logging
LOG_LEVEL = logging.INFO
logging.basicConfig(level=LOG_LEVEL)

import cdms2
import copy
import sys
import os
import json
import pcmdi_metrics
from pcmdi_metrics.pcmdi.pmp_parser import PMPParser
import collections
from collections import defaultdict

from EnsoMetrics.EnsoMetricsLib import EnsoAmpl, EnsoMu
from EnsoMetrics.pmpParser import ReadOptions

#########################################################
param = ReadOptions()

metrics = param.metrics
ninoBox = param.ninoBox

print metrics
print ninoBox
##########################################################

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
        sstFile = param.sstObsPath 
        tauxFile = param.tauxObsPath

        sstName = param.sstNameObs
        tauxName = param.tauxNameObs
    else:
        sstFile = (param.modpath.replace('MOD', mod)).replace('VAR',sstName) ## Will need land mask out at some point...!
        tauxFile = (param.modpath.replace('MOD', mod)).replace('VAR',tauxName)

        sstName = param.sstName
        tauxName = param.tauxName

    print sstFile
    print tauxFile
  
    #try:
    if 1:
        for metric in metrics:

            print metric

            if metric == 'EnsoAmpl':
                tmp_dict = EnsoAmpl(sstFile, sstName, ninoBox)
                tmp_dict['input_data'] = [sstFile]
            elif metric == 'EnsoMu':
                tmp_dict = EnsoMu(sstFile, tauxFile, sstName, tauxName, ninoBox, ninoBox)
                tmp_dict['input_data'] = [sstFile, tauxFile]
        
            enso_stat_dic[mod][metric] = tmp_dict
    else:    
    #except:
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
metrics_dictionary["RESULTS"] = enso_stat_dic  # collections.OrderedDict()

#OUT.var = var
OUT.write(
    metrics_dictionary,
    json_structure=["model", "metric", "item", "value or description"],
    indent=4,
    separators=(
        ',',
        ': '),
    sort_keys=True)

sys.exit('done')
