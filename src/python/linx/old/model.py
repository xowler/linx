## ToDo: Check for multiple models
import numpy
import linx

#{{{ Aminoacid definitions
aminoacids3 = [ 
    'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 
    'MET', 'PRO', 'PHE', 'TRP', 'SER', 
    'THR', 'ASN', 'GLN', 'TYR', 'CYS', 
    'LYS', 'ARG', 'HIS', 'ASP', 'GLU',
    'NH2', 'ACE', 'LYP' ]

aminoacids1 = [ 
    'G', 'A', 'V', 'L', 'I', 
    'M', 'P', 'F', 'W', 'S', 
    'T', 'N', 'Q', 'Y', 'C', 
    'K', 'R', 'H', 'D', 'E',
    ]
#}}}

#{{{ Read PDB/GRO files 
def parsePDBline(line):
    atom = {
        'id':       int(line[4:11]),
        'name':     line[11:16].strip(),
        'resname':  line[16:20].strip(),
        'resid':    int(line[22:26]),
        'chain':    line[20:22].strip(),
    }
    x = [0.1*float(i) for i in line[30:].split()[:3]]
    return atom, x

def parseGROline(line):
    atom = {
        'resid':    int(line[:5].strip()),
        'resname':  line[5:10].strip(),
        'name':     line[10:15].strip(),
        'id':       int(line[15:21]),
    }
    x = map(float, line[21:].split())[:3]
    return atom, x

def gro(fname):
    """
    Reads a GRO file and returns atoms and coords
    """
    return read(fname)

def pdb(fname):
    """
    Reads a GRO file and returns atoms and coords
    """
    return read(fname)

def read(fname):
    """
    Reads a PDB or GRO file and returns two arrays

        atoms: dict with characteristics of the atoms
        coords: numpy array with the coordinates of the atoms

        THE VELOCITIES ARE DROPPED!

    ex:
        gro = linx.read('../data0/md.gro')
    """
    atoms = []
    coords = []
    lines = open(fname,'r').readlines()
    if fname[-3:] == "gro":
        parseline = parseGROline
        natoms = int(lines[1])
        lines = lines[2:][:-1]
    elif fname[-3:] == "pdb":
        parseline = parsePDBline
        nlines = []
        for line in lines:
            if line[:4]=='ATOM':
                nlines.append(line) 
            if 'TER' in line[:5]:
                break
        lines = nlines
    else:
        print "File format error: %s"%fname
        assert 0

    for line in lines:
        atom, x = parseline(line)
        atoms.append(atom)
        coords.append(x)
    return numpy.array(atoms), numpy.array(coords)
#}}}

#{{{ Write PDB
def writePDB(atoms, coords, fname):
    """
    Writes a PDB file, given the atom description and the coordinates
    """
    f = open(fname, 'w')
    f.write('TITLE NO-TITLE\n')
    f.write('MODEL        1\n')
    for i in range(len(coords)):
        if len(atoms[i]['name'])==4:
            r = atoms[i]['name'][0]
            n = atoms[i]['name'][1:]
        else:
            r = ' '
            n = atoms[i]['name'] 
        
        line = "ATOM  %5i %c%-3s %3s%2s%4i    %8.3f%8.3f%8.3f  1.00  0.00            \n"%(
            atoms[i]['id']%100000,
            r,
            n,
            atoms[i]['resname'],
            atoms[i].get('chain',' '),
            atoms[i]['resid'],
            10.0*coords[i][0], 10.0*coords[i][1], 10.0*coords[i][2])
        f.write(line)
    f.write('TER\nENDMDL')
#}}}

#{{{ Write GRO
def writeGRO(atoms, coords, fname, comment='No comments', size_string=None):
    """
    Writes a GRO file, given the atom description and the coordinates
    """
    if size_string == None:
        s = linx.linear_size(coords)
        size_string = "%10.5f%10.5f%10.5f"%(s,s,s)

    f = open(fname, 'w')
    f.write('%s\n'%comment)
    f.write('%5d\n'%len(atoms))
    for i in range(len(coords)):
        line = "%5d%-3s%7s%5s%8.3f%8.3f%8.3f\n"%(atoms[i]['resid'], atoms[i]['resname'], atoms[i]['name'], atoms[i]['id'], coords[i][0], coords[i][1], coords[i][2])
        f.write(line)

    #size = "".join(["%10.5f"%float(i) for i in size])
    #f.write(size)
    f.write(size_string)
#}}}

#{{{ Atom selection
def select(atoms, prop, sel=None, inv=False, heuristics=True):
    """
    Convenience selection for attributes among the atom list. Returns a list of integers.
    If sel==None, then sel = [ True ]
    
    ex:
        gro = linx.read('md.gro')
        bb = linx.select(gro[0], 'name',['N','CA','C'])
        protein = linx.select(gro[0],'protein')
        heavy = linx.select(gro[0],'element','H',inv=True)

    Special attributes:
        protein: Selects atoms from a known aminoacid
        bb: Selects N-CA-C atoms
        dihe: bb + makes sure it starts with N such that dihedrals can be computed

    if heuristics==True (default) it will try to check for the bb to be in order and for 
    the protein to not have gaps.
    """

    if sel==None:
        sel = [True]
    if not isinstance(sel, list):
        sel = [sel]

    if prop == 'protein':  
        idx = select(atoms, 'resname',aminoacids3)
        if heuristics:
            assert not idx[-1]==(len(idx)+1)
        return idx
    if prop == 'bb':
        return select(atoms, 'name', ['N','CA','C'])
    if prop == 'dihe':
        idx = select(atoms, 'bb')
        if atoms[idx[0]]['name'] == 'N':
            idx = idx[1:]
        if atoms[idx[0]]['name'] == 'CA':
            idx = idx[1:]
        if heuristics:
            assert False not in [ a['name']=='C'  for a in atoms[idx][0::3]]
            assert False not in [ a['name']=='N'  for a in atoms[idx][1::3]]
            assert False not in [ a['name']=='CA' for a in atoms[idx][2::3]]
        return idx
    if prop == 'element':
        f = lambda el: not el in ['C','H','N','O','S'] 
        elems = filter(f, [ atom['name'][0] for atom in atoms ])
        if not elems == []:
            print elems, 'Where not supported'
            assert 0
        if inv: 
            idx = [ i for i in range(len(atoms)) if not atoms[i]['name'][0] in sel]
        else:
            idx = [ i for i in range(len(atoms)) if atoms[i]['name'][0] in sel]

    if inv:
        idx = [ i for i in range(len(atoms)) if not atoms[i][prop] in sel] 
    else:
        idx = [ i for i in range(len(atoms)) if atoms[i][prop] in sel] 

    return idx
#}}}

def residues(atoms):
    """
        Returns a dictionary with a list of atoms per residue
    """
    residues = dict([ [j,[]] for j in  set([i['resid'] for i in atoms]) ])
    for i in range(len(atoms)):
        residues[atoms[i]['resid']].append(i)
    return residues


#{{{ Massify!
def higgs(atoms):
    """
        Gives mass to the atoms
    """
    masses = {
        'H': 1.0,
        'O': 15.9994,
        'N': 14.00674,
        'C': 12.0107,
        'S': 32.065,
        }
    mass = []
    for atom in atoms:
        mass.append(masses[atom['name'][0]])

    return numpy.array(mass)

#}}}

#{{{ Use pymol to build a pdb
def build_pdb(pdb, seq,ncap=None,ccap=None):
    """
    Builds a pdb given a sequence

    ex:
        import linx
        linx.build_pdb('/tmp/a.pdb','WGW','ACE','NH2')  
    """
    assert ncap in ['ACE','',None]
    assert ccap in ['NH2','',None]

    script = []
    if ncap=='ACE':
        script.append('cmd._alt("b")\n')

    for aa in seq:
        assert aa in aminoacids1
        script.append('cmd._alt("%s")\n'%aa.lower())

    if ccap=='NH2':
        script.append("editor.attach_amino_acid('pk1','nhh')\n")


    script.append('remove (hydro)\n')
    script.append('save %s\n'%pdb)
    script.append('quit\n')

    import tempfile
    f = tempfile.NamedTemporaryFile(suffix='.pymol', dir='/tmp', mode='w')
    f.file.writelines(script)
    f.file.flush()

    import os
    os.system('pymol -x -u %s'%f.name)
#}}}

    
