#=================================================
# Observation
#-------------------------------------------------
sstObsPath = '/clim_obs/obs/ocn/mo/tos/UKMETOFFICE-HadISST-v1-1/130122_HadISST_sst.nc'
tauxObsPath = '/clim_obs/obs/atm/mo/tauu/ERAINT/tauu_ERAINT_198901-200911.nc'

sstNameObs = 'sst'
tauxNameObs = 'tauu'

#=================================================
# Models
#-------------------------------------------------
modpath = '/work/cmip5/historical/atm/mo/VAR/cmip5.MOD.historical.r1i1p1.mo.atm.Amon.VAR.ver-1.latestX.xml'

modnames = ['ACCESS1-0', 'ACCESS1-3', 
            'BNU-ESM', 
            'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 
            'CSIRO-Mk3-6-0', 'CanCM4', 
            'GISS-E2-H-CC', 'GISS-E2-H', 'GISS-E2-R-CC', 'GISS-E2-R', 
            'HadCM3', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES', 
            'IPSL-CM5A-LR', 
            'MIROC-ESM-CHEM', 'MIROC-ESM', 'MIROC4h', 'MIROC5', 
            'MPI-ESM-LR', 'MPI-ESM-MR', 
            'inmcm4'
           ]
modnames = ['IPSL-CM5A-LR']

# Variables
sstName = 'ts'
tauxName= 'tauu'

#=================================================
# Output
#-------------------------------------------------
outpathdata = '.' # e.g. '/user/directory/output/nc'
outpathjsons = '.' # e.g. '/user/directory/output/json'
outnamejson = 'test.json'

#=================================================
# Output
#-------------------------------------------------
# Metrics
metrics = ['EnsoAmpl', 'EnsoMu']

# Variable name and nino box
ninoBox = 'nino3'
