from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import cv2
import numpy as np
import routines as myr
from skimage import measure
from sky2face import sky_to_patch,patch_to_sky

def sky2face(m):
    print('sky2face does not exist any more, please use sky2patch.')
    return None

def sky2patch(m,npatch=1):
    numpa = 12*npatch**2
    lp = int(np.sqrt(m.shape[0]/12)/npatch)
    return sky_to_patch(m,npatch,numpa, lp)
    
def patch2sky(patch):
    npix = patch.size
    npatch = int(np.sqrt(patch.shape[0]/12))
    return patch_to_sky(patch,npix,npatch)

def fortranize(m):
    return np.array(m,order='F')
       
def correlarion_fucntion(m1,m2=None,n_p=None,mask_value=np.nan):
    m1 = fortranize(m1)
    if m2 is None:
        m2 = m1
    else:
        m2 = fortranize(m2)
    if n_p is None:
        n_p = 5*np.prod(m2.shape)
    (npixr,mean,var) = myr.make_standard(m1,0)
    (npixr,mean,var) = myr.make_standard(m2,0)
    (cor,vcor) = myr.cross_cf(m1,m2,n_p,mask_value)
    return (cor,vcor)
    
def ppcf(m,th,nfmax,rmax):
    m = fortranize(m)
    lg = m.shape[0]
    (npixr,mean,var) = myr.make_standard(m,0)
#    m = m-m.mean()
#    m = m/m.std()
    (nf1,fl1) = myr.peak_f(m,th,100)

    if nf1>nfmax:
        rlist = np.random.randint(0,nf1,size=nfmax)
        nf1 = nfmax
        fl1 = fl1[:,rlist]

    fl1 = fortranize(fl1[:,:nf1])
    (ksi,vksi) = myr.ffcf(1,lg,fl1,fl1,5*nf1,1,rmax)
    return ksi
    
def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    try:
          window_size = np.abs(np.int(window_size))
          order = np.abs(np.int(order))
    except ValueError:
          raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
          raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
          raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def power_spectrum(m,size):
    nside = m.shape[0]
    
    mk = np.fft.fft2(m)
    kmax = int(1.5*nside)
    power = np.zeros(kmax)
    nn = np.zeros(kmax)
    for i in range(nside):
        for j in range(nside):
            k = int(np.sqrt(i**2+j**2))
            power[k] += np.absolute(mk[i,j])**2
            nn[k] += 1
            
    filt = nn!=0
    power[filt] = power[filt]/nn[filt]
    ls = (np.arange(1,kmax)+1)*360./size
    return ls,power[1:]*(ls*(ls+1))/(2*np.pi)/((nside*np.pi)**2)
    
def der(data):
    lx = data.shape[0]
    ly = data.shape[1]

    datap = np.zeros((lx+1,ly+1))

    datap[:lx,:ly] = data
    datap[:lx,ly:] = data[:,:1]
    datap[lx:,:ly] = data[:1,:]
    datap[lx,ly] = data[0,0]

    dx = np.diff(datap,axis=0)
    dy = np.diff(datap,axis=1)
    
    strg = np.sqrt(dx[:,:ly]**2+dy[:lx,:]**2)
    mstrg = np.mean(strg)
    strg = np.where(strg>mstrg,strg,0)
#    strg = np.where(strg<2*mstrg,1/2,0)
    ornt = np.arctan2(dx[:,:ly],dy[:lx,:])
    return [strg,ornt]

def filters(d,edd_method='sch',R=0,smoothing='g'):

    if (R!=0):
        dt = np.fft.fft2(d)
        if smoothing=='g':
            for i in range(sz):
                for j in range(sz):
                    k2 = 1.*(i*i+j*j)/d.shape[0]
                    dt[i,j]=dt[i,j]*np.exp(-k2*R*R/2)

        elif smoothing=='tp':
            for i in range(sz):
                for j in range(sz):
                    k = np.sqrt(0.001+i*i+j*j)/sz
                    dt[i,j]=dt[i,j]* 3*(np.sin(k*R)-k*R*np.cos(k*R))/(k*R)**3

        d = np.fft.ifft2(dt)
        d = abs(d)

    if edd_method=='lap':
        d = cv2.Laplacian(d,cv2.CV_64F)

    elif edd_method=='sob':
        sobelx = cv2.Sobel(d,cv2.CV_64F,1,0,ksize=3)
        sobely = cv2.Sobel(d,cv2.CV_64F,0,1,ksize=3)
        d =np.sqrt(sobelx**2+sobely**2)

    elif edd_method=='sch':
        scharrx = cv2.Scharr(d,cv2.CV_64F,1,0)
        scharry = cv2.Scharr(d,cv2.CV_64F,0,1)
        d =np.sqrt(scharrx**2+scharry**2)
        
    elif edd_method=='der':
        (dx,dy,vx,vy) = myr.vdd(d)
        d = np.sqrt(dx**2+dy**2)
        
    else:
        print('The filter name is not recognized!')
        
    return d
        
        
def pdf(m,n_c=20):
    m = np.array(m)
    n, bin_edges = np.histogram(m.flatten('F'), n_c)
    dx = (bins[1]-bins[0])/2
    bins = np.array([bin_edges[i]+bin_edges[i+1] for i in range(n_c)])
    return bins,n 
    
def stat_describe(m,m_max=3):
    mean = np.mean(m)
    std = np.std(m)
    out = [mean,std]
    
    m = m - mean
    for n in range(3,m_max+1):
        m_n = np.mean(m**n)
        if m_n >= 0:
                m_n = m_n**(1./n)
        elif m_n < 0:
                m_n = -(abs(m_n)**(1./n))
        m_n = m_n/std
        out.append(m_n)

    return np.array(out)    
    
def hotspot(data,trsh):
    hotspots = measure.label(data>trsh)
    return hotspots.max()

def coldspot(data,trsh):
    coldspots = measure.label(data<trsh)
    return coldspots.max()
    
def genus(data,trshs):
    return [hotspot(data,trsh)-coldspot(data,trsh) for trsh in trshs]

def deform(image,df,inverse=False,verbose=True):
    
    cc = 1
    if inverse:
        cc = -1
    img_output = np.zeros(image.shape, dtype=image.dtype)
    f1, f2 = df.shape
    
    if image.ndim==2:
        rows, cols = image.shape
        if (f1==rows+1) & (f2==cols+1):
            pass
        else:
            if verbose:
                print ('Warning, deformation field dimensions are not compatible, it will be cropped!')
            df = df[:rows+1,:cols+1]
        mx = np.diff(df,1,axis=0)[:,:cols]
        my = np.diff(df,1,axis=1)[:rows,:]
        for i in range(rows):
            for j in range(cols):
                offset_x = int(mx[i,j])
                offset_y = int(my[i,j])
                img_output[i,j] = image[(i-cc*offset_x)%rows,(j-cc*offset_y)%cols]
                
    elif image.ndim==3:
        rows, cols, _ = image.shape
        if (f1==rows+1) & (f2==cols+1):
            pass
        else:
            if verbose:
                print ('Warning, deformation field dimensions are not compatible, it will be cropped!')
            df = df[:rows+1,:cols+1]
        mx = np.diff(df,1,axis=0)[:,:cols]
        my = np.diff(df,1,axis=1)[:rows,:]
        for i in range(rows):
            for j in range(cols):
                offset_x = int(mx[i,j])
                offset_y = int(my[i,j])
                img_output[i,j,:] = image[(i-cc*offset_x)%rows,(j-cc*offset_y)%cols,:]
    else:
        print('It is not an image!')
        return
                           
    return img_output,mx,my
