import datetime
now = datetime.datetime.now()

#=================================================
# Observation
#-------------------------------------------------
reference_data_path = {
    'ERA-Interim': '/work/lee1043/DATA/reanalysis/ERAINT/mon/ERA-Interim_VAR_mo.xml',
    'HadISST': '/clim_obs/obs/ocn/mo/tos/UKMETOFFICE-HadISST-v1-1/130122_HadISST_sst.nc',
    'OISST': '/work/lee1043/DATA/OISST/xmls/OISST_tos_mo.xml',
    'Tropflux': '/work/lee1043/DATA/TropFlux/monthly/xmls/Tropflux_VAR_mo.xml',
    'OAFlux': '/work/lee1043/DATA/OAFlux/xmls/OAFlux_VAR_mo.xml',
    'GPCPv2.3': '/clim_obs/PMPObs/pmpobs1-5-1/atmos/mon/pr/GPCP-2-3/gn/v20180706/pr_mon_GPCP-2-3_BE_gn_197901-201803.nc', 
}

reference_data_lf_path = {
    'GPCPv2.3': '/work/lee1043/DATA/GPCP/gpcp_25_lsmask.nc'
}
#=================================================
# Models
#-------------------------------------------------
modpath = '/work/lee1043/ESGF/xmls/cmip5/historical/mon/%(variable)/cmip5.%(model).historical.%(realization).mon.%(variable).xml'
modpath_lf = '/work/lee1043/ESGF/xmls/cmip5/historical/fx/sftlf/cmip5.%(model).historical.r0i0p0.fx.sftlf.xml'

modnames = ['ACCESS1-0', 'ACCESS1-3', 'BCC-CSM1-1', 'BCC-CSM1-1-M', 'BNU-ESM',
            'CanCM4', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CESM1-FASTCHEM', 'CESM1-WACCM',
            'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CNRM-CM5-2', 'CSIRO-Mk3-6-0',
            'EC-EARTH', 'FGOALS-g2', 'FGOALS-s2', 'FIO-ESM',
            'GFDL-CM2p1', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',
            'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC',
            'HadCM3', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES',
            'INMCM4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',
            'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC4h', 'MIROC5',
            'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P', 'NorESM1-M', 'NorESM1-ME']

modnames = ['ACCESS1-0']
#modnames = ['IPSL-CM5A-LR']
#modnames = ['BCC-CSM1-1']

#=================================================
# Metrics Collection
#-------------------------------------------------
#metricsCollection = 'MC1'
metricsCollection = 'ENSO_perf'
#metricsCollection = 'ENSO_tel'

#=================================================
# Output
#-------------------------------------------------
nc_out = True
results_dir = '/work/lee1043/imsi/result_test/enso_metric' 
json_name = metricsCollection + '_' + now.strftime("%Y%m%d-%H%M")
netcdf_name = json_name + '_%(model)'

#=================================================
# Miscellaneous
#-------------------------------------------------
debug = True
