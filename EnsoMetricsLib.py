import cdms2 as cdm
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
    - sstfile    - time/lat/lon array
    - varname    - name of sst variable (tos, ts)
    - ninobox    - name of box ('nino3','nino3.4')

    Output:
    - amplMetric dict:
        - name, value, units, method, ref

    Notes:
    -----
    - could add error calculation

    '''
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
        nbox = cdu.region.domain(latitude=(-5.,5.),longitude=(210.,270.))
    else:
        print '!!! ninobox not defined in EnsoAmpl', ninobox

    # Read SST in box
    sst = fi(varname, nbox)
    cdu.setTimeBoundsMonthly(sst)

    # Compute mean annual cycle
    sstAnom = cdu.ANNUALCYCLE.departures(sst)

    # Compute standard deviation
    sstStd = statistics.std(sstAnom)

    # Create output
    amplMetric = {'name':Name, 'value':sstStd, 'units':Units, 'method':Method, ref:Ref}

    return amplMetric