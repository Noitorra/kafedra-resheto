"""
Pyplot animation example.

The method shown here is only for very simple, low-performance
use.  For more demanding applications, look at the animation
module and the examples that use it.
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import *

NX = 50
NY = 50
ITER = 500

for i in range(ITER):
    s = "Den%i" % (i)
    D = np.fromfile('Den/'+s+'.bin',dtype=float).reshape(NX,NY)
    plt.imshow(D, vmin=0.5, vmax=1.0, interpolation='nearest')
    plt.colorbar()
    plt.savefig('Den/Pic/'+s+'.png')
    plt.close()
    print("%i of %i" % (i, ITER))




