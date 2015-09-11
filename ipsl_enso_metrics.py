#!/bin/env python

# run with python -i to stay in python

# initialisation

import cdms2 as cdm
import numpy as npy
import os, string, sys  
import cdutil as cdu
from genutil import statistics
from cdms2.selectors import Selector
import MV2 as mv

npy.set_printoptions(precision = 2)

# Variable


# Directory and file names

pathin = '/prodigfs/esg/CMIP5/merge/CNRM-CERFACS/CNRM-CM5/piControl/mon/ocean/Omon/r1i1p1/latest/tos/'
file= 'tos_Omon_CNRM-CM5_piControl_r1i1p1_219001-219912.nc'
var = 'tos'

#pathin = '/home/ericglod/database/'
#file = 'HadISST1_1m_187001_199912_ORCA2_grid_T.nc'
var = 'sosstsst'

pathin = '/home/ericglod/IPSL_CM6/'
file = 'CPL6v5.17h_18500101_20191231_1M_tsol_oce.nc'
var = 'tsol_oce'

mod = pathin + file

print
print ' SST standard deviation in nino3 '

# Open file and read variable
# Select averaging region

f = cdm.open(mod)
data = f(var, latitude=(-5.,5.), longitude=(210., 270.))  # organised as data[time, latitude, longitude]
d_handle = f[var]
print data.shape

# Average in region

dataR = mv.average(mv.average(data, axis=2), axis=1)(squeeze=1)
print dataR.shape

nTime = dataR.shape[0] # number of time steps (e.g. months)
nYears = nTime / 12  # number of years

# Build mean annual cycle
 # init array
dataAC = npy.ma.ones([12], dtype='float32')*0.
for m in range(12):
  print m
  d = dataR[m::12] # select indices for month m every 12 months
  dataAC[m]   = mv.average(d) # average along time axis

print dataAC
  # Note: can be replace by [::12] directly in mv.average

# Remove mean annual cycle from data to obtain interannual anomalies

dataInter = dataR - npy.tile(dataAC, nYears)


# Compute standard deviation

dataStd = statistics.std(dataInter)


print '  ', file, var,' Standard deviation = ',dataStd


f.close()




