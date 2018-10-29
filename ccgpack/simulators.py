from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.interpolate import interp1d

class StochasticFieldSimulator(object):
    
    def __init__(self,cl,lmax=None):
        self.cl1d = interp1d(cl[:,0], cl[:,1])
        if lmax is None:
            self.lmax = cl[:,0].max()
        else:
            self.lmax = lmax
        self.clk = None
        
        assert cl[:,0].max()>=self.lmax,'Error!'
            
    def c_l(self,l_k):
        
        if l_k<self.lmax:
            try:
                return self.cl1d(l_k)
            except:
                return 0
        else:
            return 0

    def simulate(self,nside,size):
        delta_wavelength = 2*np.pi/(size*np.pi/180)
        m = np.zeros((nside,nside,4),dtype=np.complex)
        
        if self.clk is None or nside!=self.clk.shape[0]:
            delta_wavelength = 2*np.pi/(size*np.pi/180)
            X, Y = np.meshgrid(np.arange(0,nside)**2,np.arange(0,nside)**2)
            lk = np.sqrt(X+Y)*delta_wavelength
            self.clk = np.array([list(map(self.c_l, x)) for x in lk])
        
        for di in range(4):
            mr = np.random.normal(0,1,(nside,nside))
            mi = np.random.normal(0,1,(nside,nside))
#             for i in range(nside):
#                 for j in range(nside):
#                     l_k=np.sqrt((i*1.0)**2+(j*1.0)**2)*delta_wavelength    # this true when we use L as radian
#                     clk = self.c_l(l_k)
#                     m[i,j,di] = np.sqrt(clk/2)*(mr[i,j]+1j*mi[i,j])*(size/2*np.pi)**2

            m[:,:,di] = np.sqrt(self.clk/2)*(mr+1j*mi)*(size/2*np.pi)**2

            m[:,:,di] = np.rot90(np.real(np.fft.ifft2(m[:,:,di])),di)
        m = np.real(m)
        m = np.mean(m,axis=2)/2.

        return m*(nside*np.pi)
    
#class CMBSimulator(object):
#    
#    def __init__(self,cl,lmax=None):
#        self.cl1d = interp1d(cl[:,0], cl[:,1])
#        self.lmax = lmax
#        
#        assert cl[:,0].max()>lmax,'Error!'
#            
#    def c_l(self,l_k):
#        
#        if self.lmax is None or l_k<self.lmax:
#            try:
#                return self.cl1d(l_k)
#            except:
#                return 0
#        else:
#            return 0

#    def simulate(self,nside,size):
#        
#        delta_wavelength = 2*np.pi/(size*np.pi/180)
#        m = np.zeros((nside,nside,4),dtype=np.complex)

#        for di in range(4):
#            mr = np.random.normal(0,1,(nside,nside))
#            mi = np.random.normal(0,1,(nside,nside))
#            for i in range(nside):
#                for j in range(nside):
#                    l_k=np.sqrt((i*1.0)**2+(j*1.0)**2)*delta_wavelength    # this true when we use L as radian
#                    clk = self.c_l(l_k)
#                    m[i,j,di] = np.sqrt(clk/2)*(mr[i,j]+1j*mi[i,j])*(size/2*np.pi)**2

#            m[:,:,di] = np.rot90(np.real(np.fft.ifft2(m[:,:,di])),di)
#        m = np.real(m)
#        m = np.mean(m,axis=2)/2.

#        return m*(nside*np.pi)
