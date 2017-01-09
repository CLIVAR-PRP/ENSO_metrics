import cdms2 as cdm
import numpy as npy
import cdutil as cdu
from genutil import statistics
from cdms2.selectors import Selector
import MV2 as mv

#
# Libray to compute ENSO metrics
# These procedures have file names as inputs and metric as output
#

def EnsoAmpl (sstfile, varname, ninobox):
    '''
    The EnsoAmpl() function compute the standard deviation in a ninobox

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author:

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    - sstfile    - time/lat/lon array using Omon values
    - varname    - name of sst variable (tos, ts)
    - ninobox    - name of box ('nino3','nino3.4')

    Output:
    - amplMetric dict:
        - name, value, units, method, nyears, ref

    Notes:
    -----
    - could add error calculation

    '''
    cdm.setAutoBounds('on')

    # Define metric attributes
    Name   = 'ENSO amplitude'
    Units  = 'C'
    Method = 'Standard deviation of SST in '+ninobox
    Ref    = 'using CDAT calculation'

    # Open file and get time dimension
    fi = cdm.open(sstfile)
    ssth = fi[varname] # Create variable handle
    # Number of months and years
    timN = ssth.shape[0]
    yearN = timN / 12

    # define ninobox
    if ninobox =='nino3':
        #nbox = cdu.region.domain(latitude=(-5.,5.),longitude=(-150,-90))
        latBox = [-5,5]
        lonBox = [-150,-90]
    else:
        print '!!! ninobox not defined in EnsoAmpl', ninobox
    #print nbox
    # Read SST in box and average
    #sst = fi(varname, nbox)
    sst = fi(varname, latitude=(-5,5), longitude=(-150,-90))
    sstAveBox = cdu.averager(sst,axis='12',weights=cdu.area_weights(sst)).squeeze()

    # Compute anomaly wrt annual cycle and average
    #cdu.setTimeBoundsMonthly(sstAveBox)
    #sstAnom = cdu.ANNUALCYCLE.departures(sstAveBox)
    sstAnom = computeAnom(sstAveBox.data, yearN)

    # Compute standard deviation
    sstStd = statistics.std(sstAnom)

    # Create output
    amplMetric = {'name':Name, 'value':sstStd, 'units':Units, 'method':Method, 'ref':Ref, 'nyears':yearN}

    return amplMetric

def computeAnom(var1d, nYears):
    '''
    Compute internannual anomaly of 1D time serie as cdu.ANNUALCYCLE.departures complains about units
    :param var:
    :return:
    '''
    varAC = npy.ma.ones([12], dtype='float32')*0.
    for m in range(12):
        d = var1d[m::12] # select indices for month m every 12 months
        varAC[m]   = mv.average(d) # average along time axis
    print varAC-273.15
    varInter = var1d - npy.tile(varAC, nYears) # compute anomaly
    return varInter