#
# Use matplotlib to create a metrics table
#
import matplotlib.pyplot as plt


def EnsoMetricsTable(models,metrics, figName):

    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    # Labels
    colLabels = ("ENSO Amplitude", r'$\mu$')
    rowLabels = models
    # make the table
    the_table = ax.table(cellText=metrics,colWidths=[.2,.2],
          colLabels=colLabels,rowLabels=rowLabels,colLoc='center',loc='center')
    # height of rows
    table_props = the_table.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells: cell.set_height(0.08)
    # size of fonts
    the_table.set_fontsize(12)
    #the_table.scale(1.5, 1.5)
    plt.tight_layout(rect=[0.05,0.15,0.95,.95])

    plt.show()
    #plt.savefig(figName+'.jpeg')

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