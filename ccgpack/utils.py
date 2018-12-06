from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import math
import pprint
import random
import urllib
import requests
from time import gmtime, strftime


def connected(url='http://www.google.com/'):
    timeout=5
    try:
        _ = requests.get(url, timeout=timeout)
        return True
    except requests.ConnectionError:
        print("Internet connection problem.")
        return False

#def internet_check():
#    try:
#        urllib.urlopen('http://216.58.192.142')
#        return True
#    except urllib2.URLError as err: 
#        return False

def report(count, blockSize, totalSize):
  	percent = int(count*blockSize*100/(totalSize))
  	sys.stdout.write("\r%d%%" % percent + ' complete')
  	sys.stdout.flush()

def download(getFile, saveFile=None):
    assert connected(),'Error! check your Internet connection.'
    if saveFile is None:
        saveFile = getFile.split('/')[-1]
    sys.stdout.write('\rFetching ' + saveFile + '...\n')
    try:
        urllib.urlretrieve(getFile, saveFile, reporthook=report)
    except:
        urllib.request.urlretrieve(getFile, saveFile, reporthook=report)
    sys.stdout.write("\rDownload complete, saved as %s" % (saveFile) + '\n\n')
    sys.stdout.flush()
    
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

def pop_percent(i,ntot):
    sys.stdout.write("\r{:3.1f}%".format(100.*(i+1)/ntot)+ ' complete')
    sys.stdout.flush()
