import sys,socket
from EnsoMetricsLib import EnsoAmpl
#
# Wrapper of EnsoMetricsLib for testing
#

# Define model and simulation
model = 'IPSL-CM5-LR'
simu = 'historical'


# define where we are working
hostname = socket.gethostname()
if 'locean-ipsl.upmc.fr' in hostname:
    baseDir = '/Volumes/hciclad/data/Density_binning/'
elif 'waippo.local' in hostname:
    baseDir = '/Volumes/hciclad/data/Density_binning/'
elif 'private.ipsl.fr' in hostname:
    baseDir = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5A-LR/historical/mon/ocean/Omon/r1i1p1/latest'
elif 'crunchy.llnl.gov' in hostname:
    baseDir = '/work/guilyardi/'
else:
    print hostname
    sys.exit('Unknown hostname')

# Variable name and nino box
varName = 'tos'
ninoBox = 'nino3'


sstFile = baseDir+'/tos/tos_Omon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc'

print sstFile

# Call metrics calculation
ensoAmpl = EnsoAmpl(sstFile, varName, ninoBox)

print ensoAmpl['name'] #+':'+ensoAmpl['value']+'('++ensoAmpl['units']+')'
print ensoAmpl['method']