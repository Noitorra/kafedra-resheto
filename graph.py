__author__ = 'dimaxx'

import matplotlib.pyplot as plt
import numpy
from numpy import *

NX = 40
NY = 80
max_files = 100

data_folder = 'data/'

for i in range(max_files):
    s = "%i" % i
    D = numpy.fromfile(data_folder+'Den/'+s+'.bin', dtype=float).reshape(NX, NY)
    plt.imshow(D, vmin=0.0, vmax=1.0, interpolation='nearest')
    plt.colorbar()
    plt.savefig(data_folder+'Den/Pic/'+s+'.png')
    plt.close()
    print("%i of %i" % (i, max_files))

