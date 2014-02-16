__author__ = 'dimaxx'

import matplotlib.pyplot as plt
import numpy
from numpy import *

NX = 60
NY = 120
max_files = 50
gas_num = 2

for gas in range(gas_num):
  data_folder = 'data/Gas' + '%i' % gas + '/'
  for i in range(max_files):
    s = "%i" % i
    D = numpy.fromfile(data_folder+'Den/'+s+'.bin', dtype=float).reshape(NX, NY)
    plt.imshow(D, vmin=0.5, vmax=1.5, interpolation='nearest')
    plt.colorbar()
    plt.contour(D, colors='black')
    plt.savefig(data_folder+'Den/Pic/'+s+'.png', dpi=100)
    plt.close()
    print("%i of %i" % (i, max_files))
  for i in range(max_files):
    s = "%i" % i
    D = numpy.fromfile(data_folder+'Temp/'+s+'.bin', dtype=float).reshape(NX, NY)
    plt.imshow(D, interpolation='nearest')
    plt.colorbar()
    plt.contour(D, colors='black')
    plt.savefig(data_folder+'Temp/Pic/'+s+'.png', dpi=100)
    plt.close()
    print("%i of %i" % (i, max_files))

