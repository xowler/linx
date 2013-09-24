import numpy as np
from numpy.linalg import norm
import linx
import pylab
import sys

def load_data(temps,obs,verb=True):#{{{
    data = []
    for t in temps:
        data.append([])
        for o,cols in obs:
            f = '/home/jimenc/data/disk01/jimenc/bhp/analysis/data0/%s/data-%s.npz'%(o,t)
            if verb: print 'Loading %s'%f
            f = np.load(f)['arr_0']
            if f.shape == f.flatten().shape: f = f.reshape((f.shape[0],1))
            assert max(cols) < len(f[0])
            if len(f[0])>1: data[-1].append(f.T[cols])
            else: data[-1].append(f.T)

        data[-1] = np.vstack(data[-1]).T

    return data, [ "%s-%s"%(i,k) for i,j in obs for k in j ]#}}}

def trjs_old(temps,data):#{{{
    rex = np.load('/home/jimenc/data/disk01/jimenc/bhp/analysis/data0/rex/data.npz')['arr_0'].astype('int')[:,temps].T
    trjs = []
    blank = [np.nan]*len(data[0][0])
    for i in range(len(temps)):
        print 'trj-fying temp %s'%temps[i] 
        trj = []
        for c in range(len(rex[i])-1):
            trj.extend(data[i][4*c+1:4*c+5])
            if rex[i][c] != rex[i][c+1]: trj.append(blank)
    
        trj.append(blank)
        trj = np.vstack(trj)
        trjs.append(trj)

    return np.vstack(trjs)#}}}

def trjs(temps, data):#{{{
    rpt = np.load('/home/jimenc/data/disk01/jimenc/bhp/analysis/data0/rex/data.npz')['arr_0'].astype('int')
    tpr = rpt.argsort(axis=1)
    cmd = "&".join(["(tpr!=%s)"%i for i in temps])
    print 'Demuxing'
    exec('idx = np.argwhere(%s)'%cmd)
    for i,j in idx: tpr[i][j] = -1
    tpr = tpr.T

    lut = {}
    for i,j in enumerate(temps): lut[j] = i

    trjs = []
    blank = [np.nan]*len(data[0][0])

    sys.stdout.write('Trjfying replica ')
    for r in range(len(tpr)):
        sys.stdout.write('%s, '%r)
        sys.stdout.flush()
        for c in range(len(tpr[r])-1):
            if tpr[r][c]<0: continue
            trjs.extend(data[lut[tpr[r][c]]][4*c+1:4*c+5])
            if tpr[r][c] != tpr[r][c+1]: trjs.append(blank)

        trjs.append(blank)

    return np.array(trjs) #}}}

def FELs(data,lbls,sub=1,lims=None):#{{{
    ndata = np.vstack(data)[::sub].T
    if lims == None: lims = [ linx.lim(i,0.001) for i in ndata ]
    pairs = [[i,j] for i in range(len(lims)) for j in range(i+1,len(lims))]
    hs = [ np.histogram2d(ndata[j], ndata[i], bins=50, range=[lims[j],lims[i]] ) for i,j in pairs ] 
    for i,h in enumerate(hs):
        x = 0.5*( 2*h[2][:-1] + 0*h[2][1:] )
        y = 0.5*( 2*h[1][:-1] + 0*h[1][1:] )
        hs[i] = (x,y,h[0],lbls[pairs[i][0]],lbls[pairs[i][1]])

    return hs#}}}

def bin(trs):#{{{
    print 'Modifying trs IN PLACE'
    idx = np.where(np.isfinite(trs[:,0]))[0]
    lims = [linx.lim(trs[idx,i],0.005) for i in range(len(trs[0]))]
    for i,(m,M) in enumerate(lims):
        nidx = np.where(trs[idx,i]<m)[0]
        trs[idx[nidx],i] = m
        nidx = np.where(trs[idx,i]>M)[0]
        trs[idx[nidx],i] = M
        trs[idx,i] = 99.0*(trs[idx,i] - m)/(M-m)
        trs[idx,i] = (trs[idx,i]).astype('int') 

    return lims #}}}

def coords(l,lims):#{{{
    l = ((0.5+np.array(l))/100).T
    nl = []
    for i in range(len(l)):
        nl.append(l[i]*(lims[i][1]-lims[i][0]))
        nl[-1] = nl[-1]+lims[i][0]

    return np.array(nl).T#}}}

def Plot(h,ff=None):#{{{
    linx.params()
    if ff == None: 
        pylab.figure(figsize=(10,10))
        plots = linx.multi_plot({'w':[1],'s':[{'n':'p'}]}, .1,.1,0.9,0.9,0.005)
        ff = plots['p']
    ff.contourf(h[0],h[1],-np.log10(h[2]), 100)
    ff.contour(h[0],h[1],-np.log10(h[2]), 10, colors='k', lw=1,linestyles='solid', linewidths=1.5,alpha=.8)
    ff.set_xlim(h[0][[0,-1]])
    ff.set_ylim(h[1][[0,-1]])
    ff.set_xlabel(h[3])
    ff.set_ylabel(h[4])#}}}

def interpolate(l,N=2, rsum=False, rtype='int'):#{{{
    n = []
    for j in range(len(l[0])):
        n.append([])
        for i in range(len(l)-1):
            n[-1].extend(np.linspace(l[i][j],l[i+1][j],100))

    n = np.array(n).T
    ls = [ norm([ n[i][j]-n[i+1][j] for j in range(len(l[0]))])  for i in range(len(n)-1) ]
    nls = np.array(ls)/sum(ls)
    nls = np.cumsum(nls)
    nl = [l[0]]
    for i in np.linspace(0,1,N)[1:-1]:
        i = np.where(nls>i)[0][0]
        nl.append(n[i])

    nl.append(l[-1])
    if rsum: return np.array(nl).astype(rtype), sum(ls)
    return np.array(nl).astype(rtype)#}}}

def drift_point(trs,pt,d=2,D=8):#{{{
    s="&".join(["(np.abs(trs[:,%s]-%s)<%s)"%(i,pt[i],d) for i in range(len(pt))])
    exec('idx = np.where(%s)[0]'%s)

    if len(idx)<3:
        if (d+2)<D: return drift_point(trs,pt,d+2)
        else: return pt

    n = [ i for i in trs[idx+1] if np.isfinite(i[0]) ]
    n = np.average(n,axis=0)
    #dif = n-pt
    #norm = np.linalg.norm(n)
    #if norm>(len(pt)*D): 
    #    print n,
    #    dif = dif/norm*(len(pt)*D)
    #    n = pt+dif
    return n.astype('int')#}}}

def drift(trs, s):#{{{
    N = len(s)
    ns = np.array([drift_point(trs,i) for i in s])
    ns = interpolate(ns,N+1)
    ns = interpolate(ns,N+2)
    ns = interpolate(ns,N+1)
    ns = interpolate(ns,N)
    return ns#}}}

def relax(trs,ss,cb=None):#{{{
    ravg = ss[-1]
    for i in range(100): 
        ss.append( drift(trs,ss[-1]) )
        if not cb==None: cb[0](ss[-1],*cb[1])
        dif = np.linalg.norm(ravg-ss[-1])/len(ss[-1])
        if dif<.01*len(ss[-1][0]): break
        ravg = 0.55*ss[-1]+0.45*ravg
        if i%10==0: sys.stdout.write('%s'%i)
        sys.stdout.write('.')
        sys.stdout.flush()
        
    return ss[1:] #}}}

def do(trs, ss, t, m, cbp,lims):#{{{
    for n in range(len(ss[0]),100,2):
        ns, l = interpolate(ss[-1],n,1)
        ss.append(ns)
        if (1.0*l/n)<np.sqrt(2*len(ns[0])) : break
        print '%s images'%n
        relax(trs,ss,[cbp,[m,t]])
        np.savez('%s/strings.npz'%t,ss)
        np.savez('%s/strings_coords.npz'%t,[coords(i,lims) for i in ss ])#}}}

def quiver(trs,bins,lims):#{{{
    idx = np.where( (trs[:,0]>lims[0][0])& (trs[:,0]<lims[0][1])& (trs[:,1]>lims[1][0])& (trs[:,1]<lims[1][1]))[0]
    idx = idx[np.where(np.isfinite(trs[idx+1][:,0]))]

    ds = lims[:,1]-lims[:,0]
    nvf = [[[] for i in range(bins)] for j in range(bins)]
    x,y = ((trs[idx]-lims[:,0])/ds*(bins-1)).astype('int').T
    z = trs[idx+1]-trs[idx]
    for i,j,k in zip(x,y,z): nvf[i][j].append(k)

    for i in range(bins):
        for j in range(bins):
            if len(nvf[i][j])<10: nvf[i][j] = np.array([[0,0]])
            nvf[i][j] = np.average(nvf[i][j],axis=0)

    return np.linspace(lims[0][0],lims[0][1],bins),np.linspace(lims[1][0],lims[1][1],bins),np.array(nvf)
#}}}

def dg(s,trs,temp,d=1,rounds=0):#{{{
    n = []
    for pt in s:
        cmd = '&'.join([ "(np.fabs(trs[:,%s]-%s)<%s)"%(i,pt[i],d) for i in range(len(trs[0])) ])
        exec('n.append(len(np.where(%s)[0]))'%cmd)

    n = np.array(n)
    nn = n
    if rounds>0: nn = n[1:]+n[:-1]
    for i in range(rounds-1): nn = nn[1:]+nn[:-1]
    dg = -8.314*temp*np.log(nn)/1000
    return dg - min(dg), 1.0*nn/sum(nn)#}}}

def qs(s,trs,d=5):
    ds = (np.max(s,axis=0) - np.min(s,axis=0) )/100
    ts = np.vstack([s[1]-s[0],s[2:]-s[:-2],s[-1]-s[-2]])
    ts = [i/np.linalg.norm(i) for i in ts]
    ps = []
    Ds = []
    drifts = []
    for c,pt in enumerate(s):
        cmd = '&'.join(['( np.fabs(trs[:,%s]-%s)<%s )'%(i,pt[i],d*ds[i]) for i in range(len(trs[0]))])
        exec('idx = np.where(%s)[0]'%cmd)
        drift = trs[idx+1] - trs[idx]
        drift = drift[np.isfinite(drift[:,0])]
        drifts.append(np.average(drift,axis=0))
        Ds.append(np.corrcoef(drift.T))
        ps.append(len(drift))
        sys.stdout.write('%s '%c)
        sys.stdout.flush()

    qs = [ np.dot(np.dot(ts[i],Ds[i]),ts[i])/ps[i] for i in range(len(Ds))]
    qs = np.cumsum(qs)/np.sum(qs)
  
    return qs,ps,Ds,drifts










#for i in *png; do convert -trim $i $i; done
#mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc  -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1  "mf://?????.png" -mf type=png:fps=18 -o output.avi
#mencoder -mc 0 -vf harddup -ofps 25 -ovc lavc -of lavf -lavfopts format=flv -lavcopts vcodec=flv:vbitrate=236:mbd=2:mv0:trell:v4mv:cbp:last_pred=3 output.avi -o output.flv
def movie(T, M):#{{{
    import os
    for png in range(M): os.system('convert -trim %s/%05d.png %s/%05d.png'%(T,png,T,png))
    os.system('mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc  -lavcopts vcodec=mpeg4:vhq:trell:mbd=2:vmax_b_frames=1:v4mv:vb_strategy=0:vlelim=0:vcelim=0:cmp=6:subcmp=6:precmp=6:predia=3:dia=3:vme=4:vqscale=1  "mf://%s/?????.png" -mf type=png:fps=18 -o %s/output.avi'%(T,T))#}}}

def get_samples(s,obs,tol0=100,tol1=15,n=5):
    lims = np.array((np.min(s,axis=0), np.max(s,axis=0)))
    data,lbs = load_data(range(40),obs)

    sts = [[] for i in s]
    for pti,pt in enumerate(s):
        print 'Point %s of %s'%(pti,len(s))
        tol = tol0

        while tol>tol1 and len(sts[pti])==0:
            ds = (lims[1]-lims[0])/tol
            print '-- tol = %s'%tol
            for t in range(40)[::-1]:
                cmd = ['( np.fabs(data[%s][:,%s]-%s)<%s )'%(t,i,pt[i],ds[i]) for i in range(len(pt)) ]
                cmd = '&'.join(cmd)
                exec('idx = np.where( %s )[0]'%cmd)
                sts[pti].extend([(t,i) for i in idx])
                if len(sts[pti])>n: break

            tol = 0.85*tol

        if len(sts[pti])>n: 
            idx = (np.random.random(n)*len(sts[pti])).astype('int')
            sts[pti] = [ sts[pti][i] for i in idx]

    return sts




import commands
def get_pdbs(s,obs,n=5,dirname='render',tol0=80,tol1=15):
    sts = get_samples(s,obs,tol0,tol1,n)
    sts = np.vstack([[ [i,j,k] for (j,k) in sts[i] ] for i in range(len(sts))])
    ts = set(sts[:,1])

    gro = linx.gro('/home/jimenc/data/disk01/jimenc/bhp/data0/md.gro')
    [commands.getoutput('mkdir %s/%s'%(dirname,i)) for i,j in enumerate(s)]

    print 'Temperatures: ', ts
    for t in ts:
        print '';
        print 'Temperature %s'%t
        xtc = linx.xtc('/home/jimenc/data/disk01/jimenc/bhp/data0/md%s.xtc'%t)
        f = 0
        while xtc.next()==0:
            idx = np.where(sts[:,1]==t)[0]
            if f in sts[idx][:,2]:
                pt = sts[idx][np.where(sts[idx][:,2]==f)[0]][0][0]
                linx.writePDB(gro[0],xtc.x,'%s/%s/%s-%s.pdb'%(dirname,pt,t,f))
                
            f = f+1
            if f>max(sts[idx][:,2]): break
            if f%100000==0: sys.stdout.write(str(f))
            if f%10000==0:
                sys.stdout.write('.')
                sys.stdout.flush()
            
        
    


        


    

        
