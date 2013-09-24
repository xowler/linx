import numpy 
import pylab
import matplotlib

def clines():
    ax = gca()
    while len(ax.lines)>0: a = ax.pop()


def lim(data,delta=.01):
    h,x = numpy.histogram(data,bins=1000)
    h = 1.0*numpy.cumsum(h)/sum(h)
    return x[numpy.where(h>delta)[0][0]], x[numpy.where(h>1-delta)[0][0]]

    

def h1d(data, hopts={}, popts={}, ret=False, f=None):
    h = {
        "bins":100,
        "normed":1,
    }   
    h.update(hopts)
    p = {}
    p.update(popts)

    y,x = numpy.histogram(data,**h)
    x =  (0.0 + x[:-1] + x[1:])/2
    if f==None:
        f = lambda x: x
    y = [f(i) for i in y]
    pylab.plot(x, y,**p)
    if ret:
        return x,y

import scipy.ndimage
def h2d(x,y, hopts={}, popts={},log=True, levels=10,clipped=False,gaussian=True,ret=False,lines=True):
    h={
        'bins':50,
    }
    if clipped:
        x0,x1 = lim(x)
        y0,y1 = lim(y)
        h.update({'range':[[y0,y1],[x0,x1]]})
    h.update(hopts)
    p = {}
    p.update(popts)
    h = numpy.histogram2d(y, x, **h)
    x = ( 0.0 + h[2][:-1] + h[2][1:] )/2
    y = ( 0.0 + h[1][:-1] + h[1][1:] )/2

    h = h[0]
    if gaussian:
        m = min(h[numpy.where(h>0)])
        h[numpy.where(h==0)] = m

    h = numpy.log10(h)
    if gaussian:
        h = scipy.ndimage.gaussian_filter(h,sigma=1,order=0)
    
    pylab.contourf(x,y,h, 10*levels)
    if lines: pylab.contour(x,y,h, levels, colors='k', linestyles='solid', linewidths=1,**p)
    if clipped:
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
    
    if ret: return (x,y,h)


def params(set=True):
    rcParams = {}
    #rcParams['text.usetex'] = True
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.serif'] = 'Avant Garde'
    rcParams['font.size'] = 25
    rcParams['figure.subplot.left'] = .15
    rcParams['figure.subplot.bottom'] = .15
    rcParams['legend.fancybox'] = 1
    rcParams['legend.loc'] = 0
    rcParams['legend.labelspacing'] = 0
    rcParams['legend.handletextpad'] = 0
    rcParams['figure.figsize'] = (10,7)
    rcParams['axes.grid'] = True
    #Problems when drawing contours
    rcParams['lines.linewidth'] = 5
    rcParams['axes.linewidth'] = 4
    rcParams['xtick.major.pad'] = 8
    rcParams['xtick.major.size'] = 8
    rcParams['ytick.major.pad'] = 8
    rcParams['ytick.major.size'] = 8
    rcParams['grid.linewidth'] = 3
    rcParams['grid.color'] = '#888888'

    if set:
        matplotlib.rcParams.update(rcParams)
    else:
        return rcParams


def multi_plot0(desc, left, bottom, width, height, sep, horizontal):
    res = []
    s = sum(desc)
    for d in desc:
        x0 = left
        y0 = bottom
        if horizontal:
            w = (width - sep*(len(desc)-1))*d/s
            h = height
            left = left + w + sep
        else:   
            w = width
            h = (height-sep*(len(desc)-1))*d/s
            bottom = bottom + h +sep
        res.append([x0,y0,w,h])
    return res
    

def multi_plot(desc, left, bottom, width, height, sep, horizontal=True):
    """
    ex:
            desc = {
                'w':(3,5),
                's':[
                    {'n':'histo'},
                    {
                        'w':(2,3,3),
                        's':[
                            {'n':'nrf'},
                            {'n':'dg'},
                            {'n':'ff'},
                        ]
                    }
                ]
            }
            
            figure(figsize=(14,8))
            plots = linx.multi_plot(desc, .1,.1,0.7,0.8,0.005)

    """
    res = {}
    if 'w' in desc.keys():
        boxes = multi_plot0(desc['w'], left, bottom, width, height, sep, horizontal)
        for i,box in enumerate(boxes):
            res.update(multi_plot(desc['s'][i],box[0],box[1],box[2],box[3],sep, not horizontal))
    else:
        if not desc['n'] == 'BLANK':
            return {desc['n']:pylab.axes([left,bottom,width,height])}
    
    return res

