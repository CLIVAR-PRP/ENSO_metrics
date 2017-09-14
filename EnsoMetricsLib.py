import os
import cdms2 as cdm
import numpy as npy
import cdutil as cdu
from genutil import statistics
from cdms2.selectors import Selector
from monthly_variability_statistics import *
import MV2 as mv

#
# Libray to compute ENSO metrics
# These procedures have file names as inputs and metric as output
#
def Ensoampln3 (sstfile, sstname):
    ensoampl = Ensoampl(sstfile,sstname, nino3)
    return ensoampl

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
    ------
    - amplMetric dict:
        - name, value, units, method, nyears, ref, varname

    Method:
    -------
    - uses tools from PCMDI monthly_variability_statistics library

    Notes:
    -----
    - TODO: add error calculation to stddev (function of nyears)

    '''
    # Temp corrections for cdms2 to find the right axis
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
        nbox = cdu.region.domain(latitude=(-5.,5.),longitude=(-150,-90))
    else:
        print '!!! ninobox not defined in EnsoAmpl', ninobox
    # Read SST in box
    sst = fi(sstname, nbox)
    fi.close()

    # Average, in box, compute anomaly wrt annual cycle and std dev
    sstStd = interannual_variabilty_std_annual_cycle_removed(sst)

    # Create output
    amplMetric = {'name':Name, 'value':sstStd, 'units':Units, 'method':Method, 'nyears':yearN, 'ref':Ref}

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
    ------
    - muMetric dict:
        - name, value, units, method, nyears, ref, varname, intercept, nonlinearity

    Method:
    -------
    - uses tools from PCMDI monthly_variability_statistics library

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
    fsst = cdm.open(sstfile)
    ftaux = cdm.open(tauxfile)
    ssth = fsst[sstname] # Create variable handle
    # Number of months and years
    timN = ssth.shape[0]
    yearN = timN / 12

    # Read SST and Taux in boxes
    n3box = cdu.region.domain(latitude=(-5.,5.),longitude=(-150,-90))
    sst = fsst(sstname, n3box)

    n4box = cdu.region.domain(latitude=(-5.,5.),longitude=(160,210))
    taux = ftaux(tauxname, n4box)

    # Average and compute regression of interannual anomaly
    muSlope = get_slope_linear_regression_from_anomaly(taux,sst, 0) # (all points)
    muSlopePlus = get_slope_linear_regression_from_anomaly(taux,sst, 1) # (positive SSTA = El Nino)
    muSlopeNeg = get_slope_linear_regression_from_anomaly(taux,sst, -1) # (negative SSTA = La Nina)

    # Change units
    muSlope     = muSlope * 1000.
    muSlopePlus = muSlopePlus * 1000.
    muSlopeNeg  = muSlopeNeg * 1000.

    # Create output
    muMetric = {'name':Name, 'value':muSlope, 'units':Units, 'method':Method, 'nyears':yearN, 'ref':Ref, \
                'nonlinearity':muSlopeNeg-muSlopePlus}

    return muMetric


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