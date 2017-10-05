#!/Users/ericg/Projets/CMIP/Metrics/Metrics/bin/python
##!/usr/local/uvcdat/latest/bin/python

# run with python -i to stay in python

import cdms2 as cdm
import numpy as npy
import os, string, sys  
import cdutil as cdu
from genutil import statistics
from cdms2.selectors import Selector
import MV2 as mv

npy.set_printoptions(precision = 2)

var = 'tos'

pathin = '/work/guilyardi/metrics_cmip/cmip5clims/'

lst = os.listdir(pathin + var)

print ' Number of data files: ',len(lst)




# target grid
fileg='/work/guilyardi/metrics_cmip/obs/tos/UKMETOFFICE-HadISST-v1-1/ac/tos_pcmdi-metrics_Omon_UKMETOFFICE-HadISST-v1-1_198002-200501-clim.nc'
g   = cdm.open(fileg)
ref = g(var)
grd = ref.getGrid()
sel = Selector(latitude=(-5.,5.), longitude=(-150.,-90.))
refampl = mv.max(cdu.averager(ref(sel),axis='xy'),axis=0)-mv.min(cdu.averager(ref(sel),axis='xy'),axis=0)

#w=sys.stdin.readline() # stop the code here. [Ret] to keep going
print
print ' Annual mean SST RMSE in nino3 '
#for l in range(len(lst)):
for l in range(10):
  i=lst[l]
  mod = string.split(i,'.')[1]
  f = cdm.open(pathin + var + '/' + i)
  d = f(var)
  d = d.regrid(grd,regridTool='ESMF',regridMethod='linear')
  # build "annual" mean
  dy   = mv.average(d(sel), axis=0)(squeeze=1)
  refy = mv.average(ref(sel),axis=0)(squeeze=1)
  # 
  # compute RMSE between model and obs
  msstn3_ms = statistics.rms(dy,refy,axis='xy')


  #
  # AC amplitude in nino3 (max-min) try std dev as well (check whic one in Bellenger et al.)
  dampl = mv.max(cdu.averager(d(sel),axis='xy'),axis=0)-mv.min(cdu.averager(d(sel),axis='xy'),axis=0)
  msstn3_ac = dampl - refampl

  print '  ', mod, var,' ', msstn3_ms, msstn3_ac


#  aavg = cdu.averager(d(sel), axis = 'xy')
#  metric = aavg - ref(sel)

#  stdev=float(statistics.std(aavg))

#  print mod,'  ', aavg.shape, ' ',stdev 

  f.close()

#  aavg.id=var/Users/ericg/Projets/CMIP/Metrics/Metrics/
#  g=cdm.open(mod+"_out.nc","w+")
#  g.write(aavg)
g.close()



# d.info
# t=d.getTime() or t1=d.getAxis(0)
# print all: levs[:]
# import scipy
# import genutil (local dev)
# dir(genutil) doc for completion
# select: from genutil import statistics
