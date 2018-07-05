#! /usr/bin/python2.7

import ccgpack as ccg
import numpy as np
#import matplotlib.pylab as plt

from scipy import misc
m = misc.imread('../images/einstein.jpg')

cm = ccg.curvelet(m,4)

#plt.subplot(1,2,1)
#plt.imshow(m)
#plt.subplot(1,2,2)
#plt.imshow(cm)
#plt.show()

