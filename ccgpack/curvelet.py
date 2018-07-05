import ctypes
import numpy as np
import os

path = '/'.join(os.path.abspath(__file__).split('/')[:-1])
fftw_path = '/cpp_src/fftw-2.1.5/fftw/.libs/libfftw.so.2'
crv_path = '/cpp_src/curvelet.so'

print os.path.join(path,fftw_path)
print os.path.join(path,crv_path)

ctypes.cdll.LoadLibrary(os.path.join(path,fftw_path))
curlib = ctypes.cdll.LoadLibrary(os.path.join(path,crv_path))

def curvelet(m,r_scale,n_scales=7,n_wedges=10,ac=1):
    assert r_scale<=n_scales,'Incompatible scale request.'
    m = np.array(m, dtype=np.double)
    nx = m.shape[0]
    ny = m.shape[1]
    aptr = m.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))
    curlib.curvelet(aptr,nx,ny,n_scales,r_scale-1,n_wedges,ac)
    return m

