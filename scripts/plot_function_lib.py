import math
import time


def opentemplate(window, *args) :
    if args[0] in window.listelements('template'):
        return window.gettemplate(args[0])
    else:
        return window.createtemplate(*args)


def organize(nplots) :
    """ Calculate a repartition of nplots between nrows and ncolumns
    """
    nrows = int(math.sqrt(nplots+2))
    ncolumns = int(math.ceil(float(nplots)/nrows))
    return nrows, ncolumns


def deftemplate(window, nrows=1, ncolumns=1, nplots=None, reftemplate="tpl_nu", left_margin=.0, right_margin=.0,
                bottom_margin=.05, top_margin=.05, horizontal_margin=.0, vertical_margin=.0, ratio='G',
                one_legend_only=1, axe_space=False):
    """
       Create a template with nrows x ncolumns plots
       if nplots is specified, calculates nrows and ncolumns dynamically
           (nrows and ncolumns are not taken into account)
       horizontal_margin=.05  # horizontal margin between plot
       vertical_margin=.05    # vertical margin between plot
       axe_space = T leave some space automatically for labels on x/y axes
    """
    if nplots is not None :
        nrows, ncolumns = organize(nplots)
    plot = opentemplate(window,reftemplate)
    tbn=""
    if axe_space is True:
        left_margin = 0.02
        horizontal_margin = 0.01
        vertical_margin = max(vertical_margin, 0.005)
    # define the plot area
    plot_area_x1 = left_margin
    plot_area_x2 = 1. - right_margin
    plot_area_y1 = bottom_margin
    plot_area_y2 = 1. - top_margin
    plot_area_width = plot_area_x2 - plot_area_x1
    plot_area_heigth = plot_area_y2 - plot_area_y1
    # computes the legend ratio compared to source template
    leg_ref_width = (plot.legend.x2 - plot.legend.x1)
    leg_ref_heigth = (plot.legend.y2 - plot.legend.y1)
    # maximize the size of each plot while keeping some space for one legend
    if one_legend_only is True:
        a_plot_area_width = (plot_area_width - (ncolumns - 1) * horizontal_margin) / ncolumns
        a_plot_area_heigth = (plot_area_heigth - (nrows - 1) * vertical_margin) / nrows
        a_legend_heigth = leg_ref_heigth * a_plot_area_heigth / 1.
        a_legend_width = leg_ref_width * a_plot_area_width / 1.
    else:
        a_legend_heigth = leg_ref_heigth * plot_area_heigth / 1.
        a_legend_width = leg_ref_width * plot_area_width / 1.
        if abs(plot.legend.y2 - plot.legend.y1) < abs(plot.legend.x2 - plot.legend.x1):
            # horizontal legend change y
            a_plot_area_heigth = (plot_area_heigth - a_legend_heigth - (nrows - 1) * vertical_margin) / nrows
            a_plot_area_width = (plot_area_width - (ncolumns - 1) * horizontal_margin) / ncolumns
        else :
            a_plot_area_width = (plot_area_width - a_legend_width - (ncolumns - 1) * horizontal_margin) / ncolumns
            a_plot_area_heigth = (plot_area_heigth - (nrows - 1) * vertical_margin) / nrows
    a_plot_data_x1 = plot.data.x1 * a_plot_area_width
    a_plot_data_y1 = plot.data.y1 * a_plot_area_heigth
    a_plot_data_x2 = plot.data.x2 * a_plot_area_width
    a_plot_data_y2 = plot.data.y2 * a_plot_area_heigth
    a_plot_leg_y1 = plot.legend.y1 * a_plot_area_heigth
    a_plot_leg_x1 = plot.legend.x1 * a_plot_area_width
    if one_legend_only is False:
        if abs(plot.legend.y2 - plot.legend.y1) < abs(plot.legend.x2 - plot.legend.x1):
            a_plot_data_y1 = (plot.data.y1 - leg_ref_heigth) * a_plot_area_heigth
        else:
            a_plot_data_x1 = (plot.data.x1 - leg_ref_width) *  a_plot_area_width
    natural_horizontal_margin = a_plot_area_width - (a_plot_data_x2 - a_plot_data_x1)
    natural_vertical_margin = a_plot_area_heigth - (a_plot_data_y2 - a_plot_data_y1)
    if one_legend_only is False:
        if abs(plot.legend.y2 - plot.legend.y1) < abs(plot.legend.x2 - plot.legend.x1):
            a_legend_width = a_legend_width + natural_horizontal_margin * (ncolumns - 1)
        else :
            a_legend_heigth = a_legend_heigth + natural_vertical_margin * (nrows - 1)
    t = opentemplate(window, 'tmp'+ '%5i' % (time.time() % 1 * 1e5), reftemplate)
    t.scale(min([0.95 / ncolumns, 0.95 / nrows]))
    if ratio is 'G':
        t.ratio_linear_projection(0, 360, -90, 90, box_and_ticks=1)
    # initialisation
    temp = {}
    ntpl = 0
    # loop through rows
    for ir in range(nrows):
        # loop through columns
        for ic in range(ncolumns):
            tnm = tbn + str(ir + 1) + "-" + str(ic + 1) + "_of_" + str(nrows) + "-" + str(ncolumns) +\
                  '%5i' % (time.time() % 1 * 1e5)
            t = opentemplate(window, tnm, t.name)
            x1_ref = left_margin + (a_plot_area_width + horizontal_margin) * ic
            y1_ref = 1. - top_margin - (a_plot_area_heigth + vertical_margin) * ir - a_plot_area_heigth
            t.moveto(x1_ref + a_plot_data_x1, y1_ref + a_plot_data_y1)
            # avoid having the plots in the left hand side if there is only on column (for example)
            if (t.data.x2 - t.data.x1) <= 0.7 * (a_plot_data_x2 - a_plot_data_x1):
                delta = (a_plot_data_x2 - a_plot_data_x1) - (t.data.x2 - t.data.x1)
                t.moveto(x1_ref + a_plot_data_x1+delta / 2., y1_ref + a_plot_data_y1)
            if (t.data.y2 - t.data.y1) <= 0.7 * (a_plot_data_y2 - a_plot_data_y1):
                delta = (a_plot_data_y2 - a_plot_data_y1) - (t.data.y2 - t.data.y1)
                t.moveto(x1_ref + a_plot_data_x1, y1_ref + a_plot_data_y1 + delta / 2.)
            if one_legend_only is False:
                if ir == 0 and ic == 0:
                    t.legend.x1 = left_margin + a_plot_leg_x1
                    t.legend.x2 = t.legend.x1 + a_legend_width
                    t.legend.y2 = t.legend.y1 + a_legend_heigth
                else:
                    t.legend.priority = 0
            else:
                t.legend.priority = 1
            temp[ir, ic] = t
            temp[ntpl] = t
            ntpl += 1
    return temp
