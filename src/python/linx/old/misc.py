import numpy as np

def xpm(fname, flip=True):
    """
        Reads and parses a xpm file
        if flip then start from bottom left
        ex:
           xpm = linx.xpm('dm.xpm')
           imshow(xpm,interpolation='nearest',origin='lower')
    """
    data = [line for line in open(fname).read().split('\n')  if not line.strip()=='' and not line[0:2]=='/*' ]
    n,m,l,s = map(int, data[1][1:-2].split())
    #print n,m,l,s
    vals = {}
    for i in range(l):
        key,val = data[2+i].split('"')[1:4:2]
        key,val = key[:s], float(val)
        vals[key] = val
    #print vals

    res = np.zeros([n,m])
    for i,line in enumerate(data[l+2:]):
        line = line[1:-1]
        #print line
        for j in range(m):
            word = line[s*j:s*j+s]
            if not len(word) == s:
                continue
            if word[0] in ['"']:   
                continue

            if flip:
                res[n-1-i][j] = vals[word]
            else:
                res[i][j] = vals[word]
        #print res[n-1-i]

    return res

import sys
def kmeans(indata, k, iters=10,centers=None, derr=0.001, debug=True):
    """
        Performs kmeans on the input data, returns
            centers
            assignations
            error
            delta_error

    example:
        ai = vstack([
            random.normal([0,0],[.5,.5],(100,2)),
            random.normal([2,0],[.5,.5],(100,2)),
            random.normal([2.5,2.5],[.2,.5],(100,2))])
        centers, aa, err, derr = linx.kmeans(ai, 3)
        [plot(centers[i][0],centers[i][1],'o%s'%['r','g','b'][i],ms=40) for i in range(3)]
        [plot(ai[where(aa==i),0],ai[where(aa==i),1],'o%s'%['r','g','b'][i]) for i in range(3)]
    """
    data = indata
    if debug: print 'Initiating centers'
    if centers==None:
        centers = [ 0*data[0] ]
        ds = np.array([ np.linalg.norm(centers[0]-i) for i in data ])
        while len(centers)<k+1:
            if debug:   
                sys.stdout.write('.')
                sys.stdout.flush()
            centers.append(data[np.argmax(ds)])
            nds = np.array([ np.linalg.norm(centers[-1]-i) for i in data ])
            idx = np.where(nds<ds)
            ds[idx] = nds[idx]

        centers.pop(0)

    nerr = 9999999999999
    if debug: print '\nK-Means'
    for iter in range(iters):
        aa = np.array([ np.argmin([np.linalg.norm(xi-ci) for ci in centers]) for xi in data])
        centers = [ np.average(data[np.where(aa==i)],axis=0) for i in range(k)]
        err = nerr
        nerr = np.linalg.norm([ np.linalg.norm([np.linalg.norm(xi-centers[i]) for xi in data[np.where(aa==i)]]) for i in range(k)])
        assert nerr<2*err
        if debug: 
            sys.stdout.write('%s:%.3f,%.0f '%(iter,nerr,abs(nerr - err)/nerr)/derr)
            sys.stdout.flush()
        if abs(nerr - err)/nerr < derr:
            break

    return centers, aa, nerr, abs(nerr - err)/nerr

