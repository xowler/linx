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

