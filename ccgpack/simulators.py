

class StochasticFieldSimulator(object):
    
    def __init__(self,cl,lmax=None):
        self.cl1d = interp1d(cl[:,0], cl[:,1])
        self.lmax = lmax
        
        assert cl[:,0].max()>lmax,'Error!'
            
    def c_l(self,l_k):
        
        if self.lmax is None or l_k<self.lmax:
            try:
                return self.cl1d(l_k)
            except:
                return 0
        else:
            return 0

    def simulate(self,nside,size):
        
        delta_wavelength = 2*np.pi/(size*np.pi/180)
        m = np.zeros((nside,nside,4),dtype=np.complex)

        for di in range(4):
            mr = np.random.normal(0,1,(nside,nside))
            mi = np.random.normal(0,1,(nside,nside))
            for i in range(nside):
                for j in range(nside):
                    l_k=np.sqrt((i*1.0)**2+(j*1.0)**2)*delta_wavelength    # this true when we use L as radian
                    clk = self.c_l(l_k)
                    m[i,j,di] = np.sqrt(clk/2)*(mr[i,j]+1j*mi[i,j])*(size/2*np.pi)**2

            m[:,:,di] = np.rot90(np.real(np.fft.ifft2(m[:,:,di])),di)
        m = np.real(m)
        m = np.mean(m,axis=2)/2.

        return m*(nside*np.pi)
    
class CMBSimulator(object):
    
    def __init__(self,cl,lmax=None):
        self.cl1d = interp1d(cl[:,0], cl[:,1])
        self.lmax = lmax
        
        assert cl[:,0].max()>lmax,'Error!'
            
    def c_l(self,l_k):
        
        if self.lmax is None or l_k<self.lmax:
            try:
                return self.cl1d(l_k)
            except:
                return 0
        else:
            return 0

    def simulate(self,nside,size):
        
        delta_wavelength = 2*np.pi/(size*np.pi/180)
        m = np.zeros((nside,nside,4),dtype=np.complex)

        for di in range(4):
            mr = np.random.normal(0,1,(nside,nside))
            mi = np.random.normal(0,1,(nside,nside))
            for i in range(nside):
                for j in range(nside):
                    l_k=np.sqrt((i*1.0)**2+(j*1.0)**2)*delta_wavelength    # this true when we use L as radian
                    clk = self.c_l(l_k)
                    m[i,j,di] = np.sqrt(clk/2)*(mr[i,j]+1j*mi[i,j])*(size/2*np.pi)**2

            m[:,:,di] = np.rot90(np.real(np.fft.ifft2(m[:,:,di])),di)
        m = np.real(m)
        m = np.mean(m,axis=2)/2.

        return m*(nside*np.pi)
