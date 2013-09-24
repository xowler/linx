from core import p1d


import matplotlib
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
matplotlib.rcParams.update(rcParams)
