
from EnsoMetrics import EnsoMetricsTable

Models =['IPSL-CM5A-LR','IPSL-CM5A-MR']
EnsoMetrics =[[0.82,4.1],
              [1.2,4.5]]

fig=EnsoMetricsTable(Models,EnsoMetrics, 'EnsoMetrics')
