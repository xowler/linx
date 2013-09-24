import numpy as np

#{{{ RMSD 
def rmsd_svd(a,b):
    assert a.shape == b.shape
    n = len(a)
    a = a-np.sum(a,axis = 0)/n
    b = b-np.sum(b,axis = 0)/n
    s = np.linalg.svd(np.dot(a.T,b),compute_uv=False)
    res =  np.sqrt( (np.sum(a*a+b*b) - 2*sum(s))/n )
    if res == np.nan:
        res = 0
    return res
        

def rmsd_quaternions(a,b):
    assert a.shape == b.shape
    n = len(a)
    a = a-np.sum(a,axis = 0)/n
    b = b-np.sum(b,axis = 0)/n
    R = np.dot(a.T,b)
    M = np.zeros((4,4))
    M[0][0] = R[0][0]+R[1][1]+R[2][2]
    M[1][1] = R[0][0]-R[1][1]-R[2][2]
    M[2][2] =-R[0][0]+R[1][1]-R[2][2]
    M[3][3] =-R[0][0]-R[1][1]+R[2][2]
    M[0][1] = M[1][0] = R[1][2]-R[2][1]
    M[0][2] = M[2][0] = R[2][0]-R[0][2]
    M[0][3] = M[3][0] = R[0][1]-R[1][0]
    M[1][2] = M[2][1] = R[0][1]+R[1][0]
    M[1][3] = M[3][1] = R[0][2]+R[2][0]
    M[2][3] = M[3][2] = R[1][2]+R[2][1]
    return np.sqrt( (np.sum(a*a+b*b) - 2*max(np.linalg.eig(M)[0]))/n )

rmsd = rmsd_svd
#}}}

#{{{ Dihedrals
from numpy import dot, sqrt, arccos, pi, cross, sin, cos
def dihedral(c1, c2, c3, c4):
    """ 
    Computes the dihedral angle (in degrees) between four atoms.

    c1, c2, c3, c4 are the positions of the atoms

    ex:
        gro = linx.read('../data0/md.gro')
        bb = linx.select(gro[0], 'bb')
        x = gro[1][bb]
        linx.dihedral(x[0],x[1],x[2],x[3])
    """
    v1 = c1-c2
    v2 = c3-c2
    v3 = c4-c3
    
    v1p =  v1 - v2*dot(v1,v2)/dot(v2,v2)
    v3p =  v3 - v2*dot(v3,v2)/dot(v2,v2)
    
    c = dot(v1p,v3p)/sqrt(dot(v1p,v1p)*dot(v3p,v3p))
    if c>1:
        c = 1
    elif c<-1:
        c = -1

    angle = arccos(c)
    
    if dot(cross(v1p,v3p),v2) < 0:
        angle = -angle
    
    return angle*180/3.1415926

 
# phi is C-N-CA-C and psi es N-CA-C-N
def dihedrals(backbone):
    """
    Computes the dihedrals phi and psi of the backbone. Assumes that the first atom 
    is a C and the order is C-N-CA-C

    ex:
        gro = linx.read('../data0/md.gro')
        bb = linx.select(gro[0], 'dihe')
        linx.dihedrals(gro[1][bb])

    """
    i = 0
    res = []
    while len(backbone)-i>4:
        phi = dihedral(*(backbone[i:i+4]))
        psi = dihedral(*(backbone[i+1:i+5]))
        res.append([phi,psi])
        i = i+3

    return np.array(res)
#}}}

#{{{ linear_size
def linear_size(coords, idx=False):
    """
    Returns the maximum linear size of the molecule (or a good approximate to it)
    """
    coords = coords - np.sum(coords, axis=0)/len(coords)
    n1 = np.argmax([ np.linalg.norm(x) for x in coords])
    coords = coords - coords[n1]
    n2 = np.argmax([ np.linalg.norm(x) for x in coords])
    d = np.linalg.norm(coords[n1]-coords[n2])

    if idx:
        return (d, n1, n2)
    else:
        return d
#}}}

#{{{ radius_of_gyration
def radius_of_gyration2(coords, mass=None):
    if mass==None:
        mass = np.array([1]*len(coords))
    assert len(coords) == len(mass)
    ncoords = coords - sum(coords)/len(coords)
    res = [ mass[i]*np.dot(ncoords[i],ncoords[i]) for i in range(len(coords)) ]
    return  sum(res)/sum(mass)

def radius_of_gyration(coords, mass=None):
    return np.sqrt(radius_of_gyration2(coords,mass))
#}}}

#{{{ inertia_tensor
def inertia_tensor(coords, mass=None):
    """ 
    Calculates the tensor of inertia of a distribution of points.
    If no mass, all the points are equally weighted 
    Ex:
        gro = linx.read('../data/folded.gro')
        ca = linx.select(gro[0], 'name', ['CA'])
        print linx.inertia_tensor(gro[1][ca])

    """
    print 'Not working!'
    assert 0
    if mass==None:
        mass = [1]*len(coords)

    coords = coords - sum(coords)/len(coords)

    xx =  np.dot(mass, [ (y*y+z*z) for (x,y,z) in coords])
    yy =  np.dot(mass, [ (x*x+z*z) for (x,y,z) in coords])
    zz =  np.dot(mass, [ (x*x+y*y) for (x,y,z) in coords])
    xy = -np.dot(mass, [ x*y for (x,y,z) in coords])
    xz = -np.dot(mass, [ x*z for (x,y,z) in coords])
    yz = -np.dot(mass, [ y*z for (x,y,z) in coords])
    return np.array([
    [xx,xy,xz],
    [xy,yy,yz],
    [xz,yz,zz],
    ])

    #
    #return numpy.dot(coords.T, coords)/len(coords)
#}}}

#{{{ principal_axes
def principal_axes(coords, mass=None, center=True):
    """
    Returns the eigenvalues and eigenvectors of the inertia tensor 
    """
    I = inertia_tensor(coords, mass, center)
    w,v = np.linalg.eig(I)
    return w,v
#}}}

#{{{ Contacts
def distances(coords1, coords2):
    pairs = [(i,j) for i,x in enumerate(coords1) for j,y in enumerate(coords2) ]
    d = lambda (x, y): np.linalg.norm(coords1[x]-coords2[y])    
    return map(d, pairs)

def residue_distances(atoms, coords, skip=2, what=np.min, sym=True):
    """
        Returns a contact map, if sym is false, returns a flat list
    """
    res = [[atoms[0]['resid'], []]]
    for atom,x in zip(atoms,coords):
        if not atom['resid']==res[-1][0]:
            res.append([atom['resid'],[]])
        res[-1][1].append(x)

    i = len(res)
    ds = np.zeros([i,i])
    for i,(r,x)  in enumerate(res):
        for j,(r,y) in enumerate(res):
            if j>i+skip:
                continue
            ds[i][j] = what(distances(x,y))
            ds[j][i] = ds[i][j]

    if sym:
        return ds
    lut = [(j,i) for j in range(len(res)) for i in range(len(res)) ] 
    lut = [(i,j) for i,j in lut if j>i+skip]
    return map(lambda (i,j): ds[i][j], lut), lut

def contacts(atoms, coords, skip=2, cutoff=.45):
    cs, idx = residue_distances(atoms, coords, skip=skip, sym=False)
    cs =  [i for c,i in zip(cs,idx) if c<cutoff]
    return cs

#}}}

def dists(res, coords, skip=2):
    ds = []
    for ri, ai in res.items():
        for rj, aj in res.items():
            if rj>ri+skip:
                ci = coords[ai]
                cj = coords[aj]
                d = np.min([np.linalg.norm(xi-xj) for xi in ci for xj in cj ])
                ds.append([[ri,rj],d])
    return zip(*ds)


    


# sasa Lee, B.; Richards, F. M.Interpretation of protein structures: Estimation of static accessibility J. Mol. Biol. 1971, 55,
# sasa Hasel, W.; Hendrickson, T. F.; Still, W. C.A rapid approximation to the solvent accessible surface areas of atoms Tetrahedron 
# dssp Kabsch, W.; Sander, C.Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features

