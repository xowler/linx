import numpy


def pdb(fname):
    def parseline(line):
        atom = {
            'id':       int(line[4:11]),
            'name':     line[11:16].strip(),
            'resname':  line[16:20].strip(),
            'resid':    int(line[22:26]),
            'chain':    line[20:22].strip(),
        }
        x = [0.1*float(i) for i in line[30:].split()[:3]]
        return atom, x
    atoms = []
    coords = []
    lines = open(fname,'r').readlines()
    nlines = []
    for line in lines:
        if line[:4]=='ATOM':
            nlines.append(line) 
        if 'TER' in line[:5]:
            break
    lines = nlines

    for line in lines:
        atom, x = parseline(line)
        atoms.append(atom)
        coords.append(x)
    return numpy.array(atoms), numpy.array(coords)

