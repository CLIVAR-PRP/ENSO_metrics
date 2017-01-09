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

def EnsoAmpl (sstfile, sstname, ninobox):
    '''
    The EnsoAmpl() function compute the SST standard deviation in a ninobox

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author:

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    - sstfile    - time/lat/lon SST array using Omon values
    - sstname    - name of sst variable (tos, ts)
    - ninobox    - name of box ('nino3','nino3.4')

    Output:
    - amplMetric dict:
        - name, value, units, method, nyears, ref, varname

    Notes:
    -----
    - TODO: add error calculation to stddev (function of nyears)

    '''
    cdm.setAutoBounds('on')

    # Define metric attributes
    Name   = 'ENSO amplitude'
    Units  = 'C'
    Method = 'Standard deviation of SST in '+ninobox
    Ref    = 'Using CDAT std dev calculation'

    # Open file and get time dimension
    fi = cdm.open(sstfile)
    ssth = fi[sstname] # Create variable handle
    # Number of months and years
    timN = ssth.shape[0]
    yearN = timN / 12

    # define ninobox
    if ninobox =='nino3':
        #nbox = cdu.region.domain(latitude=(-5.,5.),longitude=(-150,-90))
        latBox1 = -5  ; latBox2 = 5
        lonBox1 = 210; lonBox2 = 270
    else:
        print '!!! ninobox not defined in EnsoAmpl', ninobox
    #print nbox
    # Read SST in box and average
    #sst = fi(varname, nbox)
    sst = fi(sstname, latitude=(latBox1,latBox2), longitude=(lonBox1,lonBox2))
    sstAveBox = cdu.averager(sst,axis='12',weights=cdu.area_weights(sst)).squeeze()

    # Compute anomaly wrt annual cycle and average
    #cdu.setTimeBoundsMonthly(sstAveBox)
    #sstAnom = cdu.ANNUALCYCLE.departures(sstAveBox)
    sstAnom = computeAnom(sstAveBox.data, yearN)

    # Compute standard deviation
    sstStd = statistics.std(sstAnom)

    # Create output
    amplMetric = {'name':Name, 'value':sstStd, 'units':Units, 'method':Method, 'nyears':yearN, 'ref':Ref, 'varname':sstname}

    return amplMetric

def EnsoMu (sstfile, tauxfile, sstname, tauxname):
    '''
    The EnsoMu() function compute the regression of nino4 tauxA over nino3 sstA

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author:

    Created on Mon Jan  9 11:05:18 CET 2017

    Inputs:
    ------
    - sstfile    - time/lat/lon SST array using Omon values
    - tauxfile   - time/lat/lon Taux array using Omon values
    - sstname    - name of sst variable (tos, ts)
    - ninobox    - name of box ('nino3','nino3.4')

    Output:
    - amplMetric dict:
        - name, value, units, method, nyears, ref, varname

    Notes:
    -----
    - TODO: add error calculation to stddev (function of nyears)

    '''
    cdm.setAutoBounds('on')

    # Define metric attributes
    Name   = 'Bjerknes feedback (mu)'
    Units  = '10-3 N/m2/C'
    Method = 'Regression of nino4 tauxA over nino3 sstA'
    Ref    = 'Using CDAT regression calculation'

    # Open file and get time dimension
    fi = cdm.open(sstfile)
    ssth = fi[sstname] # Create variable handle
    # Number of months and years
    timN = ssth.shape[0]
    yearN = timN / 12

    # define nino boxes

    latn31 = -5  ; latn32 = 5
    lonn31 = 210; lonn32 = 270

    latn41 = -5  ; latn42 = 5
    lonn41 = 160; lonn42 = 210




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
    varInter = var1d - npy.tile(varAC, nYears) # compute anomaly
    return varInter