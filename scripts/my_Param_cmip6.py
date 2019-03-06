import datetime
import os

# =================================================
# Background Information
# -------------------------------------------------
mip = 'cmip6'
#exp = 'piControl'
exp = 'historical'

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
modpath = '/work/lee1043/ESGF/xmls/%(mip)/historical/mon/%(variable)/%(mip).%(model).historical.%(realization).mon.%(variable).xml'
modpath_lf = '/work/lee1043/ESGF/xmls/%(mip)/historical/fx/sftlf/%(mip).%(model).historical.r0i0p0.fx.sftlf.xml'

modnames = ['all']

#modnames = ['ACCESS1-0']
#modnames = ['BCC-CSM1-1']
#modnames = ['IPSL-CM5A-LR']

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
#case_id = "{:v%Y%m%d-%H%M}".format(datetime.datetime.now())
case_id = "{:v%Y%m%d}".format(datetime.datetime.now())
#results_dir = '/work/lee1043/imsi/result_test/enso_metric/' + case_id 
results_dir = os.path.join(
    '/p/user_pub/pmp/pmp_results/pmp_v1.1.2',
    '%(output_type)', 'enso_metric',
    mip, exp, case_id)
json_name = '_'.join([mip, exp, metricsCollection, case_id])
netcdf_name = json_name + '_%(model)'

#=================================================
# Miscellaneous
#-------------------------------------------------
debug = True
