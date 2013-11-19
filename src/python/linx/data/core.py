import numpy


def lims(data,delta=.01):
    h,x = numpy.histogram(data,bins=1000)
    h = 1.0*numpy.cumsum(h)/sum(h)
    return x[numpy.where(h>delta)[0][0]], x[numpy.where(h>1-delta)[0][0]]

def h1d(data, hopts={}):
    h = {
        "bins":100,
        "normed":1,
    }   
    h.update(hopts)
    y,x = numpy.histogram(data,**h)
    x =  (0.0 + x[:-1] + x[1:])/2
    return x,y

def smooth(data,window=-1):
    if window==-1: window = max(3,.05*len(data))
    if window%2==0: window = window + 1
    w = numpy.hanning(window)
    p = numpy.ones(numpy.floor((len(w)-1)/2))
    r = numpy.convolve(w/sum(w), hstack([p*data[0],data,data[-1]*p]), mode='valid')
    assert len(r)==len(data)
    return data

    
#def h2d(x,y, hopts={}, popts={},log=True, levels=10,clipped=False,gaussian=True,ret=False,lines=True):
#    h={
#        'bins':50,
#    }
#    if clipped:
#        x0,x1 = lim(x)
#        y0,y1 = lim(y)
#        h.update({'range':[[y0,y1],[x0,x1]]})
#    h.update(hopts)
#    p = {}
#    p.update(popts)
#    h = numpy.histogram2d(y, x, **h)
#    x = ( 0.0 + h[2][:-1] + h[2][1:] )/2
#    y = ( 0.0 + h[1][:-1] + h[1][1:] )/2
#
#    h = h[0]
#    if gaussian:
#        m = min(h[numpy.where(h>0)])
#        h[numpy.where(h==0)] = m
#
#    h = numpy.log10(h)
#    if gaussian:
#        h = scipy.ndimage.gaussian_filter(h,sigma=1,order=0)
#    
#    pylab.contourf(x,y,h, 10*levels)
#    if lines: pylab.contour(x,y,h, levels, colors='k', linestyles='solid', linewidths=1,**p)
#    if clipped:
#        pylab.xlim(x0,x1)
#        pylab.ylim(y0,y1)
#    
#    if ret: return (x,y,h)
#
