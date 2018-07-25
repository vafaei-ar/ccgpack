from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import math
import pprint
import random
from time import gmtime, strftime


def ch_mkdir(directory):
    """
    ch_mkdir : This function creates a directory if it does not exist.

    Arguments:
        directory (string): Path to the directory.

    --------
    Returns:
        null.		
    """
    if not os.path.exists(directory):
          os.makedirs(directory)
          
def rank_table(mtx,xlabels,ylabels,plot_name=None,title=None):
    
    n_y = len(ylabels)
    n_x = len(xlabels)

    fig = plt.figure(figsize=(2*n_x,2*n_y))
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_aspect('auto')
    ax.imshow(mtx, cmap=plt.cm.jet,interpolation='nearest')

    width, height = mtx.shape

    rnk = np.prod(mtx.shape)-mtx.ravel().argsort().argsort().reshape(mtx.shape)

    for x in xrange(width):
        for y in xrange(height):
            ax.annotate('{:3.1f}\n rank: {:d}'.format(mtx[x][y],rnk[x][y]), xy=(y, x), 
                        horizontalalignment='center',
                        verticalalignment='center',fontsize=20);

    plt.xticks(np.arange(n_x),[str(i) for i in xlabels],fontsize=15,rotation=20)
    plt.yticks(np.arange(n_y),[str(i) for i in ylabels],fontsize=15)
    
    plt.xlabel('Disk size (Deg)',fontsize=20)
    plt.ylabel('N side',fontsize=20)
    if title is not None:
        plt.title(title,fontsize=25)
    if plot_name is not None :
        plt.savefig(plot_name+'.jpg',dpi=150,bbox_inches='tight')

