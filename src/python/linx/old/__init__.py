from model import read, gro, pdb, select, writeGRO, writePDB, higgs, build_pdb, residues
from structure import rmsd, rmsd_svd, rmsd_quaternions
from structure import dihedral, dihedrals
from structure import radius_of_gyration2, radius_of_gyration, linear_size, principal_axes, inertia_tensor
from structure import residue_distances, contacts
from structure import dists

from pyplot import  h1d, h2d, params, multi_plot, lim
        
from misc import xpm, kmeans
from pyxtc import pyxtc as xtc, pytrr as trr

__doc__ = """
ex:
	gro = linx.read('../data0/md.gro')
	bb = linx.select(gro[0], 'name',['N','CA','C'])
	protein = linx.select(gro[0], 'protein', True)
	
	xtc = linx.xtc('../data0/md-0.xtc')
	dihedrals = []
	for i in range(10):
	    xtc.next()
	    dihedrals.append(linx.dihedrals(xtc.x[bb]))
	
"""
