#
# Use matplotlib to create a metrics table
#
import matplotlib.pyplot as plt


def EnsoMetricsTable(metrics, figName):

    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    colLabels=("Model", "ENSO Amplitude", "Mu")
    the_table = ax.table(cellText=metrics,
          colLabels=colLabels,loc='center')
    plt.savefig(figName+'.gif')

    #colLabels=("Model", "ENSO Amplitude", "Mu")
    #nrows, ncols = len(clust_data)+1, len(colLables)
    #hcell, wcell = 0.3, 1.
    #hpad, wpad = 0, 0
    #fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    #ax = fig.add_subplot(111)
    #ax.axis('off')
    ##do the table
    #the_table = ax.table(cellText=clust_data,
    #          colLabels=colLabels,
    #          loc='center')
    #plt.savefig("table.png")