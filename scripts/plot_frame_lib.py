from vcs import init as VCSinit


def myboxplot(stats, col_stat=[], x_axis=[], x_dico={}, y_axis=[], y_dico={}, name_in_xlabel=False,
              draw_all_ylines=True, treshx=1, treshy=1, name='', titlec='', titlel='', titler='', xname='', yname='',
              writings=[], writings_pos_xy=[], writings_size=[], writings_halign=[], writings_color=[],
              writings_angle=[], linesb_x1x2=[], linesb_y1y2=[], linesb_width=[], linesb_type=[], linesb_color=[],
              linesa_x1x2=[], linesa_y1y2=[], linesa_width=[], linesa_type=[], linesa_color=[],  markers_xy=[],
              markers_size=[], markers_type=[], markers_color=[], colormap='', legend=True, path_plus_name_png='',
              draw_white_background=True, save_ps=False, shape='', bg=1):
    # define window and template
    x = VCSinit()
    if colormap != '' and isinstance(colormap, basestring):
        x.setcolormap(colormap)
    if shape == 'square':
        pl1 = PF.squared_template(x, nplots=1, reftemplate='tpl_axlegbig')
    else:
        pl1 = PF.right_pushed_template(x, nplots=1, reftemplate='tpl_axlegbig')
    tx2 = PF.opentextorientation(x, 'persox2', 'labx')
    tx2.height=30#25#
    tx2.angle=0
    tx2.halign = 'center'
    tx2.valign = 'top'
    if name_in_xlabel:
        tx2.height=25#30#
        tx2.angle=-25#-30#
        tx2.halign = 'right'
        tx2.valign = 'half'
        pl1[0].xlabel1.y = 0.16
    pl1[0].xlabel1.textorientation='persox2'
    ty2 = PF.opentextorientation(x, 'persoy2', 'laby')
    ty2.height=30
    pl1[0].ylabel1.textorientation='persoy2'
    # define x range and labels if user didn't
    try: 1./len(x_axis)
    except:
        x_axis, x_dico = PF.create_dico([0,len(stats)], fmin=0, fmax=len(stats))
    # define y range and labels if user didn't
    try: 1./len(y_axis)
    except:
        minmax = [ min([MV2minimum(stats[ii]) for ii in range(len(stats))]), max([MV2maximum(stats[ii]) for ii in range(len(stats))]) ]
        y_axis, y_dico = PF.create_dico(minmax)
    # param
    x1,x2,y1,y2 = x_axis[0],x_axis[1],y_axis[0],y_axis[1]
    x_dico2 = dict((elt,'') for elt in x_dico.keys())
    y_dico2 = dict((elt,'') for elt in y_dico.keys())
    param_plot = {'title':'','datawc_x1':x1,'datawc_x2':x2,'xticlabels1':x_dico,'xticlabels2':x_dico2,'datawc_y1':y1,'datawc_y2':y2,'yticlabels1':y_dico,'yticlabels2':y_dico2,'marker':'dot','markersize':1,'bg':bg}
    draw_plot = {'viewport':[pl1[0].box1.x1,pl1[0].box1.x2,pl1[0].box1.y1,pl1[0].box1.y2],'worldcoordinate':[x1,x2,y1,y2],'bg':bg}
    draw_plot2 = {'viewport':[pl1[0].box1.x1,pl1[0].box1.x2,pl1[0].box1.y1,pl1[0].box1.y2],'worldcoordinate':[0,100,0,100],'bg':bg}
    gra1 = PF.openxvsy(x, 'test1')
    gra1.datawc_x1, gra1.datawc_x2, gra1.datawc_y1, gra1.datawc_y2 = 0,100,0,100
    gra2 = PF.openxvsy(x, 'test2')
    gra2.datawc_x1, gra2.datawc_x2, gra2.datawc_y1, gra2.datawc_y2 = x1,x2,y1,y2
    # horizontal lines
    for ii in y_dico.keys():
        if y_dico[ii] == '': col = 'grey'
        else: col = 'black'
        if draw_all_ylines:
            x.drawline(name='l'+str('%2.2i'%(ii*100+2)), x=[x1,x2], y=[float(ii),float(ii)], width=1, color=x.match_color(col), **draw_plot)
        else:
            if ii>treshy: x.drawline(name='l'+str('%2.2i'%(ii*100+2)), x=[x1,x2-treshx], y=[float(ii),float(ii)], width=1, color=x.match_color(col), **draw_plot)
            else:         x.drawline(name='l'+str('%2.2i'%(ii*100+2)), x=[x1,x2], y=[float(ii),float(ii)], width=1, color=x.match_color(col), **draw_plot)
    # y zero line
    axis = CDMS2createAxis(MV2array(range(x_axis[0],x_axis[1]+1), dtype='int32'),id='events')
    zeros = MV2zeros(x_axis[1]-x_axis[0]+1)
    zeros.setAxis(0, axis)
    x.plot(zeros, pl1[0], linecolor='black', linewidth=4, markercolor='black', **param_plot)
    # lines given by the user
    for ii in range(max(len(lines_x1x2),len(lines_y1y2))):
        try:    x1x2 = lines_x1x2[ii]
        except: x1x2 = [x1,x2]
        try:    y1y2 = lines_y1y2[ii]
        except: y1y2 = [y1,y2]
        try:    widt = lines_width[ii]
        except: widt = 1
        try:    ltyp = lines_type[ii]
        except: ltyp = 'long-dash'
        try:    colo = lines_color[ii]
        except: colo = 'black'
        if isinstance(colo,basestring): colo = x.match_color(colo)
        x.drawline(name='l'+str('%2.2i'%(ii*100+3)), x=[x1x2[0],x1x2[1]], y=[y1y2[0],y1y2[1]], width=widt, ltype=ltyp, color=colo, **draw_plot)
        del colo, ltyp, widt, x1x2, y1y2
    # plot stats
    dia_size = 15 if len(stats)>30 else 30
    dx = 0.2 if len(stats)>30 else 0.1
    for ii in range(len(stats)):
        stats_boxplot = stats[ii]
        try: col = col_stat[ii]
        except: col = 'black'
        if isinstance(col,basestring): col = x.match_color(col)
        if isinstance(stats_boxplot, list):
            if len(stats_boxplot)==10:
                x.drawline(name='mini', x=[ii,    ii],       y=[stats_boxplot[0],stats_boxplot[3]], width=4, color=col, **draw_plot)
                x.drawline(name='c5.v', x=[ii-dx,ii+dx],     y=[stats_boxplot[1],stats_boxplot[1]], width=4, color=col, **draw_plot)
                x.drawline(name='d1.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[2],stats_boxplot[2]], width=4, color=col, **draw_plot)
                x.drawline(name='q1.h', x=[ii-2*dx,ii-2*dx], y=[stats_boxplot[3],stats_boxplot[5]], width=4, color=col, **draw_plot)
                x.drawline(name='q3.h', x=[ii+2*dx,ii+2*dx], y=[stats_boxplot[3],stats_boxplot[5]], width=4, color=col, **draw_plot)
                x.drawline(name='q1.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[3],stats_boxplot[3]], width=4, color=col, **draw_plot)
                x.drawline(name='medi', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[4],stats_boxplot[4]], width=4, color=col, **draw_plot)
                x.drawline(name='q3.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[5],stats_boxplot[5]], width=4, color=col, **draw_plot)
                x.drawline(name='d9.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[6],stats_boxplot[6]], width=4, color=col, **draw_plot)
                x.drawline(name='c95.v',x=[ii-dx,ii+dx],     y=[stats_boxplot[7],stats_boxplot[7]], width=4, color=col, **draw_plot)
                x.drawline(name='maxi', x=[ii,    ii],       y=[stats_boxplot[5],stats_boxplot[8]], width=4, color=col, **draw_plot)
                PF.drawmarkerlist(x, ['diamond_fill'], [col], [dia_size], [ii], [float(stats_boxplot[9])], names='%2.2i'%(31), **draw_plot)
                if ii == 0 and legend is True:
                    posx = [ii+2.5*dx for jj in range(len(stats_boxplot)-3)]
                    posy = [stats_boxplot[jj+1] for jj in range(len(stats_boxplot)-3)]
                    PF.ecritext(x, ['c2.5','D1','Q1','M','Q3','D9','c97.5'], pl1[0], gra2, ref='%2.2i'%(1), x=posx, y=posy, height=25, halign='left', bg=bg)
            elif len(stats_boxplot)==7:
                x.drawline(name='mini', x=[ii,    ii],       y=[stats_boxplot[0],stats_boxplot[2]], width=4, color=col, **draw_plot)
                x.drawline(name='d1.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[1],stats_boxplot[1]], width=4, color=col, **draw_plot)
                x.drawline(name='q1.h', x=[ii-2*dx,ii-2*dx], y=[stats_boxplot[2],stats_boxplot[4]], width=4, color=col, **draw_plot)
                x.drawline(name='q3.h', x=[ii+2*dx,ii+2*dx], y=[stats_boxplot[2],stats_boxplot[4]], width=4, color=col, **draw_plot)
                x.drawline(name='q1.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[2],stats_boxplot[2]], width=4, color=col, **draw_plot)
                x.drawline(name='medi', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[3],stats_boxplot[3]], width=4, color=col, **draw_plot)
                x.drawline(name='q3.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[4],stats_boxplot[4]], width=4, color=col, **draw_plot)
                x.drawline(name='d9.v', x=[ii-2*dx,ii+2*dx], y=[stats_boxplot[5],stats_boxplot[5]], width=4, color=col, **draw_plot)
                x.drawline(name='maxi', x=[ii,    ii],       y=[stats_boxplot[4],stats_boxplot[6]], width=4, color=col, **draw_plot)
                if ii == 0 and legend is True:
                    posx = [ii+2.5*dx for jj in range(len(stats_boxplot)-2)]
                    posy = [stats_boxplot[jj] for jj in range(len(stats_boxplot)-2)]
                    PF.ecritext(x, ['D1','Q1','M','Q3','D9'], pl1[0], gra2, ref='%2.2i'%(1), x=posx, y=posy, height=25, halign='left', bg=bg)
        del col, stats_boxplot
    # strings to write (given by the user)
    for ii in range(len(list_writings)):
        try:    height = list_writings_size[ii]
        except: height = 16
        try:    halign = list_writings_halign[ii]
        except: halign = 'left'
        try:    color = list_writings_color[ii]
        except: color = 'black'
        try:    angle = list_writings_angle[ii]
        except: angle = 0
        if isinstance(color,basestring): color = x.match_color(color)
        PF.ecritext(x, [list_writings[ii]], pl1[0], gra1, ref='%2.2i'%(ii*100+4), x=[list_writings_pos_xy[ii][0]], y=[list_writings_pos_xy[ii][1]], height=height, halign=halign, col=color, angle=angle, bg=bg)
        del angle, color, halign, height
    # lines given by the user
    for ii in range(max(len(linesa_x1x2),len(linesa_y1y2))):
        try:    x1x2 = linesa_x1x2[ii]
        except: x1x2 = [x1,x2]
        try:    y1y2 = linesa_y1y2[ii]
        except: y1y2 = [y1,y2]
        try:    widt = linesa_width[ii]
        except: widt = 1
        try:    ltyp = linesa_type[ii]
        except: ltyp = 'long-dash'
        try:    colo = linesa_color[ii]
        except: colo = 'black'
        if isinstance(colo,basestring): colo = x.match_color(colo)
        x.drawline(name='l' + str('%2.2i' % (ii * 100 + 3)), x=[x1x2[0], x1x2[1]], y=[y1y2[0], y1y2[1]], width=widt, ltype=ltyp, color=colo, **draw_plot)
        del colo, ltyp, widt, x1x2, y1y2
    # markers given by the user
    for ii in range(len(markers_xy)):
        try:    mcol = markers_color[ii]
        except: mcol = x.match_color('black')
        try:    mtyp = markers_type[ii]
        except: mtyp = 'dot'
        try:    msiz = markers_size[ii]
        except: msiz = 10
        if isinstance(mcol,basestring): mcol = x.match_color(mcol)
        x.drawmarker(name='dots'+str(ii), x=[markers_xy[ii][0]], y=[markers_xy[ii][1]], mtype=mtyp, color=mcol, size=msiz, **draw_plot)#, **draw_plot2)
        del mcol, mtyp, msiz
    # titles and axes names
    if shape == 'square':
        list_name, list_x, list_y, list_height, list_halign, list_angle = [name,title,title1,title2,ctitle,xname,yname], [0,40,25,75,50,50,-13], [103,102,103,103,103,-9,50], [35,25,35,35,30,30,30], ['left','left','center','center','center','center','center'], [0,0,0,0,0,0,-90]
        #list_name, list_x, list_y, list_height, list_halign, list_angle = [name,title,title1,title2,ctitle,xname,yname], [0,40,25,75,50,50,-13], [103,102,103,103,103,-9,25], [35,25,35,35,30,30,30], ['left','left','center','center','center','center','center'], [0,0,0,0,0,0,-90]
    else:
        list_name, list_x, list_y, list_height, list_halign, list_angle = [name,title,title1,title2,ctitle,xname,yname], [0,40,25,75,50,50,-10], [103,102,103,103,103,-9,50], [35,25,30,30,30,30,30], ['left','left','center','center','center','center','center'], [0,0,0,0,0,0,-90]
    for ii in range(len(list_name)): PF.ecritext(x, [list_name[ii]], pl1[0], gra1, ref='title'+'%2.2i'%(ii*100+5), x=[list_x[ii]], y=[list_y[ii]], height=list_height[ii], halign=list_halign[ii], angle=list_angle[ii], bg=bg)
    if save_ps:
        x.postscript(path_plus_name_png, mode='r', orientation='l', width=21, height=29.7, units='cm')
    else:
        x.png(path_plus_name_png, draw_white_background=draw_white_background)
    x.clear()
    return
