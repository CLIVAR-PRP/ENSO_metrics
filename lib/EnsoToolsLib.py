from cdms2 import setAutoBounds as CDMS2setAutoBounds
from cdms2 import open as CDMS2open
import cdutil
from genutil.statistics import linearregression as GENUTILlinearregression
from genutil.statistics import rms as GENUTILrms
from genutil.statistics import std as GENUTILstd
from inspect import stack as INSPECTstack
from MV2 import minimum as MV2minimum
from MV2 import subtract as MV2subtract
from MV2 import multiply as MV2multiply
from numpy import array as NParray
from numpy import nonzero as NPnonzero


from EnsoCollectionsLib import defCollection, metricReqs, ReferenceObservations, ReferenceRegions
import EnsoErrorsWarnings


# ---------------------------------------------------------------------------------------------------------------------#
#
# Two functions, previously included in "monthly_variability_statistics.py", that can be called by the metrics functions
#


def interannual_variabilty_std_annual_cycle_removed(d):
    d_area_avg = cdutil.averager(d, axis='xy')
    d_area_avg_anom = cdutil.ANNUALCYCLE.departures(d_area_avg)
    d_area_avg_anom_sd = GENUTILstd(d_area_avg_anom)
    return float(d_area_avg_anom_sd)


def get_slope_linear_regression_from_anomaly(y, x, sign_x, return_intercept=False, return_stderr=False):
    y_area_avg = cdutil.averager(y, axis='xy')
    x_area_avg = cdutil.averager(x, axis='xy')
    x_area_avg_anom = NParray(cdutil.ANNUALCYCLE.departures(x_area_avg))
    y_area_avg_anom = NParray(cdutil.ANNUALCYCLE.departures(y_area_avg))
    if sign_x ==  1:
        idxpos = NPnonzero(x_area_avg_anom >= 0.)
        results = GENUTILlinearregression(y_area_avg_anom[idxpos], x=x_area_avg_anom[idxpos], error=1)
    elif sign_x == -1:
        idxneg = NPnonzero(x_area_avg_anom <= 0.)
        results = GENUTILlinearregression(y_area_avg_anom[idxneg], x=x_area_avg_anom[idxneg], error=1)
    else:
        results = GENUTILlinearregression(y_area_avg_anom, x=x_area_avg_anom, error=1)
    slope_intercept, std_errors = results
    slope = slope_intercept[0]
    intercept = slope_intercept[1]
    stderr = std_errors[0]
    if return_intercept:
        if return_stderr:
            return float(slope), float(intercept), float(stderr)
        else:
            return float(slope), float(intercept)
    else:
        if return_stderr:
            return float(slope), float(stderr)
        else:
            return float(slope)


# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Functions that can be called by the metrics functions
#


def CheckUnits(tab, var_name):
    """ Takes tab, checks its units, and change it if necessary
        Works for temperature, precipitation, wind stress, current/wind velocities
        Returns : tab with the new units

        tab        : variable, tab[tt,lat, lon] or tab[tt,j,i],...
        var_name   : name of the variable
    """
    units = tab.units
    try: name_in_file = tab.short_name
    except:
        try: name_in_file = tab.id
        except: name_in_file = ''
    if var_name in ['temperature']:
        if units == 'K':
            # check if the temperature units is really K
            if MV2minimum(tab) > 200:
                # unit change of the temperature: from K to degC
                tab = MV2subtract(tab, 273.15)
                units = "degC"
            else:
                EnsoErrorsWarnings.UnlikelyUnits(var_name, name_in_file, units, INSPECTstack())
        elif units in ['C', 'degree_Celsius', 'deg_Celsius', 'deg. C', 'degCelsius', 'degree_C', 'deg_C', 'degC',
                       'degrees C']:
            # check if the temperature units is really degC
            if MV2minimum(tab) < 50:
                units = "degC"
            else:
                EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['precipitation']:
        if units == 'kg m-2 s-1':
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = MV2multiply(tab, 86400)
            units = "mm/day"
        elif units == 'mm/day':
            pass
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['wind stress']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "N/m2"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['velocity']:
        if units in ['cm s-1', 'cm/s', 'cm s**-1']:
            # unit change of the velocity: from cm/s to m/s
            tab = MV2multiply(tab, 1e-2)
            units = "m/s"
        elif units in ['m s-1', 'm/s', 'm s**-1', 'm/sec']:
            units = "m/s"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['heat flux']:
        if units in ['W/m2', 'W m-2', 'W/m^2']:
            units = "W/m2"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    elif var_name in ['pressure']:
        if units in ['N/m^2', 'Pa', 'N m-2', 'N/m2']:
            units = "Pa"
        else:
            EnsoErrorsWarnings.UnknownUnits(var_name, name_in_file, units, INSPECTstack())
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.MyWarning(list_strings)
    tab.units = units
    tab.name = name_in_file
    return tab


def ComputeMultiply(tab, number):
    return MV2multiply(tab, number)


def ComputeStd(tab):
    return float(GENUTILstd(tab))


def ComputeRms(tab1, tab2, axis, weights, centered):
    return float(GENUTILrms(tab1, tab2, axis=axis, weights=weights, centered=centered))


def SpatialAverage(tab):
    return cdutil.averager(tab, axis='xy')


def TimeAverage(tab):
    return cdutil.averager(tab, axis='t')


# Reads file and selects the given region
def ReadAndSelectRegion(filename, varname, box):
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Open file and get time dimension
    fi = CDMS2open(filename)
    # define ninobox
    region_ref = ReferenceRegions(box)
    nbox = cdutil.region.domain(latitude=region_ref['latitude'], longitude=region_ref['longitude'])
    # Read SST in box
    var = fi(varname, nbox)
    fi.close()
    return var


# Dictionary of seasons
sea_dict = dict(JAN=cdutil.JAN, FEB=cdutil.FEB, MAR=cdutil.MAR, APR=cdutil.APR, MAY=cdutil.MAY, JUN=cdutil.JUN,
                JUL=cdutil.JUL, AUG=cdutil.AUG, SEP=cdutil.SEP, OCT=cdutil.OCT, NOV=cdutil.NOV, DEC=cdutil.DEC,
                JF=cdutil.times.Seasons('JF'), FM=cdutil.times.Seasons('FM'), MA=cdutil.times.Seasons('MA'),
                AM=cdutil.times.Seasons('AM'), MJ=cdutil.times.Seasons('MJ'), JJ=cdutil.times.Seasons('JJ'),
                JA=cdutil.times.Seasons('JA'), AS=cdutil.times.Seasons('AS'), SO=cdutil.times.Seasons('SO'),
                ON=cdutil.times.Seasons('ON'), ND=cdutil.times.Seasons('ND'), DJ=cdutil.times.Seasons('DJ'),
                JFM=cdutil.times.Seasons('JFM'), FMA=cdutil.times.Seasons('FMA'), MAM=cdutil.MAM,
                AMJ=cdutil.times.Seasons('AMJ'), MJJ=cdutil.times.Seasons('MJJ'), JJA=cdutil.JJA,
                JAS=cdutil.times.Seasons('JAS'), ASO=cdutil.times.Seasons('ASO'), SON=cdutil.SON,
                OND=cdutil.times.Seasons('OND'), NDJ=cdutil.times.Seasons('NDJ'), DJF=cdutil.DJF,
                JFMA=cdutil.times.Seasons('JFMA'),FMAM=cdutil.times.Seasons('FMAM'),MAMJ=cdutil.times.Seasons('MAMJ'),
                AMJJ=cdutil.times.Seasons('AMJJ'),MJJA=cdutil.times.Seasons('MJJA'),JJAS=cdutil.times.Seasons('JJAS'),
                JASO=cdutil.times.Seasons('JASO'),ASON=cdutil.times.Seasons('ASON'),SOND=cdutil.times.Seasons('SOND'),
                ONDJ=cdutil.times.Seasons('ONDJ'),NDJF=cdutil.times.Seasons('NDJF'),DJFM=cdutil.times.Seasons('DJFM'))


# Computes the given seasonal mean
def SeasonalMean(tab, season, compute_anom=False):
    # Temp corrections for cdms2 to find the right axis
    CDMS2setAutoBounds('on')
    # Checks if the season has been defined
    try:
        tab_sea = sea_dict[season]
    except:
        list_strings = ["ERROR" + EnsoErrorsWarnings.MessageFormating(INSPECTstack()) + ": season",
                        str().ljust(5) + "unknown season: " + str(season)]
        EnsoErrorsWarnings.MyError(list_strings)
    else:
        if season in ['DJ', 'NDJ', 'DJF', 'ONDJ', 'NDJF', 'NDJF']:
            # these 'seasons' are between two years
            # if I don't custom 'tab' cdutil will compute half season mean
            # (i.e., for NDJ the first element would be for J only and the last for ND only)
            time_ax_comp = tab.getTime().asComponentTime()
            ntime = len(time_ax_comp)
            ii, jj = 0, 0
            if season == 'DJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'NDJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 11: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'DJF':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 2: break
            elif season == 'ONDJ':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 10: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 1: break
            elif season == 'NDJF':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 11: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 2: break
            elif season == 'DJFM':
                for ii in range(ntime):
                    if time_ax_comp[ii].month == 12: break
                for jj in range(ntime):
                    if time_ax_comp[ntime - 1 - jj].month == 3: break
            tab = tab[ii:ntime - jj]
        if compute_anom:
            tab = sea_dict[season].departures(tab)  # extracts 'season' seasonal anomalies (from climatology)
        else:
            tab = sea_dict[season](tab)  # computes the 'season' climatology of a tab
    return tab


# ---------------------------------------------------------------------------------------------------------------------#