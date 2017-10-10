import cdutil, genutil, numpy, cdtime


def interannual_variabilty_std_annual_cycle_removed(d):
	d_area_avg = cdutil.averager(d,axis='xy')
	d_area_avg_anom = cdutil.ANNUALCYCLE.departures(d_area_avg)
	d_area_avg_anom_sd = genutil.statistics.std(d_area_avg_anom)
	return float(d_area_avg_anom_sd)


def get_slope_linear_regression_from_anomaly(y, x, sign_x, return_intercept=False, return_stderr=False):
	y_area_avg = cdutil.averager(y, axis='xy')
	x_area_avg = cdutil.averager(x, axis='xy')
	x_area_avg_anom = numpy.array( cdutil.ANNUALCYCLE.departures(x_area_avg) )
	y_area_avg_anom = numpy.array( cdutil.ANNUALCYCLE.departures(y_area_avg) )
	if sign_x == 0:
		results = genutil.statistics.linearregression(y_area_avg_anom, x=x_area_avg_anom, error=1)
	elif sign_x ==  1:
		idxpos  = numpy.nonzero(x_area_avg_anom >= 0.)
		results = genutil.statistics.linearregression(y_area_avg_anom[idxpos], x=x_area_avg_anom[idxpos], error=1)
	elif sign_x == -1:
		idxneg  = numpy.nonzero(x_area_avg_anom <= 0.)
		results = genutil.statistics.linearregression(y_area_avg_anom[idxneg], x=x_area_avg_anom[idxneg], error=1)
	slope, intercept, stderr = results
	if   return_intercept==False and return_stderr==False:
		return float(slope)
	elif return_intercept==True  and return_stderr==False:
		return float(slope), float(intercept)
	elif return_intercept==False and return_stderr==True:
		return float(slope), float(stderr)
	elif return_intercept==True  and return_stderr==True:
		return float(slope), float(intercept), float(stderr)





def interannual_variability_seasonal_std_mean_removed(d,season_string):
	d_area_avg = cdutil.averager(d,axis='xy')
	pre_defined_seasons = ['DJF', 'MAM', 'JJA', 'SON', 'YEAR']
	if season_string in pre_defined_seasons:
		d_area_avg_anom=getattr(cdutil,season_string).departures(d_area_avg)
	else:
		CustomSeason = cdutil.times.Seasons(season_string)
		d_area_avg_anom = CustomSeason.departures(d_area_avg)
	d_area_avg_anom_sd = genutil.statistics.std(d_area_avg_anom)
	return(float(d_area_avg_anom_sd))

def get_slope_linear_regression(y,x):
	results = genutil.statistics.linearregression(y,x=x)
	slope, intercept = results
	return(float(slope))

def get_area_avg_annual_cycle_removed(d):
	d_area_avg = cdutil.averager(d,axis='xy')
	d_area_avg_anom = cdutil.ANNUALCYCLE.departures(d_area_avg)
	return(d_area_avg_anom)

def get_axis_base_dataset(var, reg, path): # to be called from Atm Feedback driver
	f = cdms.open(path)
	if debug:
		reg_timeseries = f(var, regions_specs[reg]['domain'], time = slice(0,60)) # RUN CODE FAST ON 5 YEARS OF DATA
	else:
		reg_timeseries = f(var, regions_specs[reg]['domain'])
	f.close()
	# Get area averaged and annual cycle removed 1-D time series
	reg_timeseries_area_avg_anom = get_area_avg_annual_cycle_removed(reg_timeseries)
	return(reg_timeseries_area_avg_anom)


def MatchTimeDimension(a, b):
	cdutil.setTimeBoundsMonthly(a)
	cdutil.setTimeBoundsMonthly(b)

	stime1 = a.getTime().asComponentTime()[0]
	etime1 = a.getTime().asComponentTime()[-1]

	stime2 = b.getTime().asComponentTime()[0]
	etime2 = b.getTime().asComponentTime()[-1]

        stime = max(stime1, stime2)
        etime = min(etime1, etime2)

        stime_adjust = cdtime.comptime(stime.year, stime.mon, 1) 
        etime_adjust = cdtime.comptime(etime.year, etime.mon, 31) 

        a_sliced = a(time = (stime_adjust, etime_adjust))
        b_sliced = b(time = (stime_adjust, etime_adjust))

        return(a_sliced, b_sliced)
