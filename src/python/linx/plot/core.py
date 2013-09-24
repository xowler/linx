import linx
import numpy
import pylab

def p1d(x,hopts={},popts={},f=None,ret=False):
    x,y = linx.h1d(x)
    if not f==None: y = 1/0
    pylab.plot(x,y)
    if ret: return x,y




#import numpy 

#import pylab
#import matplotlib
#
#def clines():
#    ax = gca()
#    while len(ax.lines)>0: a = ax.pop()
#
#
#
#def multi_plot0(desc, left, bottom, width, height, sep, horizontal):
#    res = []
#    s = sum(desc)
#    for d in desc:
#        x0 = left
#        y0 = bottom
#        if horizontal:
#            w = (width - sep*(len(desc)-1))*d/s
#            h = height
#            left = left + w + sep
#        else:   
#            w = width
#            h = (height-sep*(len(desc)-1))*d/s
#            bottom = bottom + h +sep
#        res.append([x0,y0,w,h])
#    return res
#    
#
#def multi_plot(desc, left, bottom, width, height, sep, horizontal=True):
#    """
#    ex:
#            desc = {
#                'w':(3,5),
#                's':[
#                    {'n':'histo'},
#                    {
#                        'w':(2,3,3),
#                        's':[
#                            {'n':'nrf'},
#                            {'n':'dg'},
#                            {'n':'ff'},
#                        ]
#                    }
#                ]
#            }
#            
#            figure(figsize=(14,8))
#            plots = linx.multi_plot(desc, .1,.1,0.7,0.8,0.005)
#
#    """
#    res = {}
#    if 'w' in desc.keys():
#        boxes = multi_plot0(desc['w'], left, bottom, width, height, sep, horizontal)
#        for i,box in enumerate(boxes):
#            res.update(multi_plot(desc['s'][i],box[0],box[1],box[2],box[3],sep, not horizontal))
#    else:
#        if not desc['n'] == 'BLANK':
#            return {desc['n']:pylab.axes([left,bottom,width,height])}
#    
#    return res
#
