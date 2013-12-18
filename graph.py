__author__ = 'Дмитрий'

import matplotlib.pyplot as plt
import numpy
from numpy import *

NX = 50
NY = 50
max_files = 500

data_folder = 'd:/Work/kafedra-resheto-VS2012/data/'

for i in range(max_files):
    s = "%i" % i
    D = numpy.fromfile(data_folder+'Den/'+s+'.bin', dtype=float).reshape(NX, NY)
    plt.imshow(D, vmin=0.5, vmax=1.0, interpolation='nearest')
    plt.colorbar()
    plt.savefig(data_folder+'Den/Pic/'+s+'.png')
    plt.close()
    print("%i of %i" % (i, max_files))

