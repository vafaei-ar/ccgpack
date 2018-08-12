from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import urllib
import shutil
import tarfile
import requests
from setuptools import find_packages
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
    
extensions = Extension(name = 'routines',
#                 extra_compile_args = ['-O3'],
                 sources = ['ccgpack/f90_src/routines.f90'],
                 libraries=['gfortran'] 
                 )

def remove_dir(dirpath):
	if os.path.exists(dirpath) and os.path.isdir(dirpath):
		  shutil.rmtree(dirpath)

#def connected(host='http://google.com'):
#    try:
#        urllib.urlopen(host)
#        return True
#    except:
#        return False
        
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
     
def untar(fname):
    if (fname.endswith("tar.gz")):
        tar = tarfile.open(fname)
        tar.extractall()
        tar.close()
        print("Extracted in Current Directory")
    else:
        print('Not a tar.gz file: '+fname)

os.chdir("ccgpack/cpp_src")
if not os.path.exists('fftw-2.1.5'):
    download('http://www.fftw.org/fftw-2.1.5.tar.gz')
    untar('fftw-2.1.5.tar.gz')
    os.remove('fftw-2.1.5.tar.gz')
os.chdir("fftw-2.1.5")

os.system('./configure --enable-shared CFLAGS=-fPIC CPPFLAGS=-fPIC')
os.system('make')
os.chdir('../')

fftw_path = './fftw-2.1.5/fftw'
cmnd = 'g++ -g -Wall -W -w -shared -c -Wno-sign-compare -Wno-unused-label -MMD -fPIC -I'+fftw_path+' -DNDEBUG '
files = ['fdct_wrapping','ifdct_wrapping','fdct_wrapping_param','function'] 

objts = ''
for fil in files:
    print('Making '+fil+' ...')
    os.system(cmnd+fil+'.cpp')
    objts = objts+fil+'.o '

os.system('g++ -shared -Wl,-soname,curvelet.so -o curvelet.so '+objts+' -fPIC -L'+fftw_path+'/.libs -lfftw')

for fil in files:
    os.remove(fil+'.d')
    os.remove(fil+'.o')

print('Curvelet library is made.')
#print os.getcwd()
os.chdir('../../')

extensions.extra_f77_compile_args = []
extensions.extra_f90_compile_args = []

requires = [] #during runtime
tests_require=['pytest>=2.3'] #for testing

PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))

setup(
	name='ccgpack',
	version='0.1.0',
	description='Computational cosmology group package.',
	author='Alireza',
	url='https://github.com/vafaeiar/ccgpack',
	packages=find_packages(PACKAGE_PATH, "ccgpack"),
	package_dir={'ccgpack': 'ccgpack'},
	include_package_data=True,
    package_data={'': ['cpp_src/curvelet.so','cpp_src/fftw-2.1.5/fftw/.libs/*']},
	install_requires=requires,
	license='GPLv3',
	zip_safe=False,
	keywords='ccgpack',
	classifiers=[
		  'Development Status :: 2 - Pre-Alpha',
		  "Intended Audience :: Science/Research",
		  'Intended Audience :: Developers',
		  "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
		  'Natural Language :: English',
		  'Programming Language :: Python :: 2',
		  'Programming Language :: Python :: 2.6',
		  'Programming Language :: Python :: 2.7',
		  'Programming Language :: Python :: 3',
		  'Programming Language :: Python :: 3.3',
		  'Programming Language :: Fortran',
	],
	ext_modules=[extensions]
)

remove_dir('build')
remove_dir('ccgpack.egg-info')
remove_dir('dist')
remove_dir('ccgpack/cpp_src/fftw-2.1.5')
os.remove('ccgpack/cpp_src/curvelet.so')
