import datetime
import os

# =================================================
# Background Information
# -------------------------------------------------
mip = 'cmip6'
#exp = 'piControl'
exp = 'historical'

#=================================================
# Miscellaneous
#-------------------------------------------------
debug = True
#debug = False
nc_out = True

#=================================================
# Observation
#-------------------------------------------------
reference_data_path = {
    'ERA-Interim': '/p/user_pub/PCMDIobs/PCMDIobs2.1/atmos/mon/VAR/ERA-INT/gn/v20190926/VAR_mon_ERA-INT_BE_gn_197901-201903.nc',
    'HadISST': '/work/lee1043/DATA/HadISSTv1.1/HadISSTv1.1.xml',
    'OISST': '/work/lee1043/DATA/OISST/xmls/OISST_tos_mo.xml',
    'Tropflux': '/work/lee1043/DATA/TropFlux/monthly/xmls/Tropflux_VAR_mo.xml',
    #'Tropflux': '/p/user_pub/PCMDIobs/PCMDIobs2.1/atmos/mon/VAR/TropFlux-1-0/gn/v20190912/VAR_mon_TropFlux-1-0_BE_gn_197901-201707.nc',
    #'OAFlux': '/work/lee1043/DATA/OAFlux/xmls/OAFlux_VAR_mo.xml',
    #'GPCPv2.3': '/p/user_pub/PCMDIobs/PCMDIobs2.1/atmos/mon/pr/GPCP-2-3/gn/v20190912/pr_mon_GPCP-2-3_BE_gn_197901-201803.nc', 
    'GPCPv2.3': '/p/user_pub/pmp/pmp_obs_preparation/orig/data/GPCP_v2.3_mon_jwl/precip.mon.mean.nc',
    #'AVISO': '/p/user_pub/PCMDIobs/PCMDIobs2.1/ocean/mon/zos/AVISO-1-0/gn/v20190912/zos_mon_AVISO-1-0_BE_gn_199210-201012.nc',
    'AVISO': '/work/lee1043/DATA/AVISO/sla_aviso_199301-201812.xml',
}

reference_data_lf_path = {
    'GPCPv2.3': '/work/lee1043/DATA/GPCP/gpcp_25_lsmask.nc'
}
#=================================================
# Models
#-------------------------------------------------
modpath = '/work/lee1043/ESGF/xmls/%(mip)/%(exp)/mon/%(variable)/%(mip).%(model).%(exp).%(realization).mon.%(variable).xml'
modpath_lf = '/work/lee1043/ESGF/xmls/%(mip)/historical/fx/sftlf/%(mip).%(model).historical.r0i0p0.fx.sftlf.xml'

modnames = ['all']

if debug:
    modnames = ['IPSL-CM5A-LR']

#=================================================
# Metrics Collection
#-------------------------------------------------
metricsCollection = 'ENSO_tel'

#=================================================
# Output
#-------------------------------------------------
case_id = "{:v%Y%m%d}".format(datetime.datetime.now())
pmprdir = '/p/user_pub/pmp/pmp_results/pmp_v1.1.2'

if debug:
    #case_id = "{:v%Y%m%d-%H%M}".format(datetime.datetime.now())
    pmprdir = '/work/lee1043/imsi/result_test'

results_dir = os.path.join(
    pmprdir,
    '%(output_type)', 'enso_metric',
    '%(mip)', '%(exp)', case_id, '%(metricsCollection)')
json_name = '_'.join(['%(mip)_%(exp)_%(metricsCollection)', case_id])
netcdf_name = json_name + '_%(model)' + '_%(realization)'
