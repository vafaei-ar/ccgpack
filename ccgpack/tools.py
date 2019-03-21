from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import cv2
import numpy as np
import routines as myr
from skimage import measure
from sky2face import sky_to_patch,patch_to_sky

def cov_mat(data):
    if data.shape[0]<2 or data.max()==data.min():
        return np.eye(data.shape[1])
    c_m = np.mean(data,axis=0)
    n_obs = data.shape[0]
    n_data = data.shape[1]
    Cov = np.zeros((n_data,n_data))
    for i in range(n_obs):
          vec = (data[i,:]-c_m).reshape(n_data,1)
          Cov += np.matmul(vec,vec.T)
    return Cov/n_obs
    
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
    m1 = fortranize(m1.astype(np.float64))
    if m2 is None:
        m2 = m1
    else:
        m2 = fortranize(m2.astype(np.float64))
        
    assert m1.shape==m2.shape,'Both data sets have to be in same shape.'
    if n_p is None:
        n_p = 5*m2.size
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
    
#def ffcf(fl1,fl2,lg,1,rmax):
#    fl1 = fortranize(fl1[:,:nf1])
#    (ksi,vksi) = myr.ffcf(1,lg,fl1,fl1,5*nf1,1,rmax)
#    return ksi
    
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
        
        
def pdf(m,bins=20):
    m = np.array(m)
    hist, bin_edges = np.histogram(m, bins)
    bins = 0.5*(bin_edges[1:]+bin_edges[:1])
    return bins,hist
    
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
    
def genus(data,trshs,standard=False):
    if standard:
        data = data-data.mean()
        data = data/data.std()        
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
    
    
#contour function and exponent added by mehdi
def carachters(data,nu):    
    contour = measure.find_contours(data,nu)   
    num_c = len(contour)
    perimeter = []
    area = []
    radius = []
    for nc in range(num_c):
        n_points = contour[nc].shape[0]
        x = contour[nc][:,0]
        y = contour[nc][:,1]
        s_nc = 0
        x_c = np.mean(x)
        y_c = np.mean(y)
        r_nc=0
        for i in range(n_points):
            s_nc += np.sqrt((x[i-1]-x[i])**2+(y[i-1]-y[i])**2)
            r_nc += (((x[i]-x_c)**2)+((y_c-y[i])**2))
        A=np.abs(0.5*np.sum(y[:-1]*np.diff(x) - x[:-1]*np.diff(y)))
        perimeter.append(s_nc)
        area.append(A)
        radius.append(np.sqrt(r_nc/n_points))
    perimeter=np.asarray(perimeter)
    area=np.asarray(area)
    radius=np.asarray(radius)    
    return [contour,perimeter,area,radius]

def pdf_c(d,R):
    l_d=len(d)
    R=R-1
    min_d=np.min(d[:,0])
    d[:,0]-=min_d
    d_x=np.max(d[:,0])/R
    n=np.zeros(R+1)
    p=np.zeros(R+1)
    for i in range(l_d):
        k=int(d[i,0]/d_x)
        n[k]=n[k]+1
        p[k]=p[k]+d[i,1]        
    len_p=len(p[n>0])
    a=np.linspace(0.5,len(p)-1.5,len(p))
    a*=d_x
    a+=min_d
    out=np.zeros((len_p,3))
    out[:,0]=a[n>0]
    out[:,1]=n[n>0]
    out[:,2]=p[n>0]/n[n>0]
    return(out)

def D_f(data,th,R,r1,r2,plot):
    [contour,perimeter,area,radius] = carachters(data,th)
    zz=np.zeros((len(perimeter),2))
    k=0
    for i in range(len(radius)):
        if radius[i]>=r1 and radius[i]<=r2:
            zz[k,1]=perimeter[i]
            zz[k,0]=radius[i]
            k+=1
    zz.resize(k,2)
    ee=pdf_c(zz,R)
    new_c=[]
    for i in range(len(ee[:,1])):
        if ee[i,1]>0:
            new_c.append([ee[i,0],ee[i,2]])
    eee=np.array(new_c)
    x=np.log10(eee[:,0])
    y=np.log10(eee[:,1])
    fi2 = np.polyfit(x,y,1)
    if plot==1:
        print('D_f = ',fi2[0])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('R')
        plt.ylabel('<s>')    
        plt.scatter(10**x,10**y)
        plt.plot(10**x,10**(x*(fi2[0])+fi2[1]),'r--')
        plt.show()
    return fi2[0]

def etta(data,th,p_min,p_max,R,s_min,s_max,plot):
    [contour,perimeter1,area1,radius] = carachters(data,th)
    perimeter = perimeter1[perimeter1<p_max]
    per_filt = perimeter[perimeter>p_min]
    R=int((np.max(per_filt)-p_min)/((s_max-s_min)/R))
    P_s, bin_edges = np.histogram(per_filt,bins=R)
    P_s=P_s/P_s.sum()
    s=[]
    for i in range(R):    
        s.append(0.5*(bin_edges[i]+bin_edges[i+1]))
    c_out=[]
    for i in range(len(P_s)):
        if P_s[i]>0:   
            c_out.append([s[i],P_s[i]])
    c_out=np.array(c_out)
    s=c_out[:,0]
    P_s=c_out[:,1]
    for i in range(len(s)):
        if s[i]>s_min:
            s_min = i
            break
    for i in range(len(s)):
        if s[i]>s_max:
            s_max = i-1
            break      
    P_s=P_s[s_min:s_max]
    s=s[s_min:s_max]
    fit=np.polyfit(np.log(s),np.log(P_s),1)
    etta = (fit[0]*-1)+1
    if plot==1:
        print('etta = ',etta)
        plt.scatter(s,P_s)
        plt.xlabel('S')
        plt.ylabel('P(S)')
        plt.show()
        ys=np.exp((np.log(s)*fit[0])+fit[1])
        plt.loglog(s,P_s,'o')
        plt.loglog(s,ys,label='slope='+str(fit[0]))
        plt.xlabel('S')
        plt.ylabel('P(S) (log scale)')
        plt.show()
        plt.plot(s,P_s*(s**(etta-1)),'.')
        plt.xlabel('S')
        plt.ylabel('P(S)*s**(etta-1)')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(0.001,1000)
        plt.show()
    return etta

def kessi(data,th,area_min,area_max,R,A_min,A_max,plot):
    [contour,perimeter1,area1,radius] = carachters(data,th)
    arear = area1[area1<area_max]
    per_filt = arear[arear>area_min]
    R=int((np.max(per_filt)-area_min)/((A_max-A_min)/R))
    P_A, bin_edges = np.histogram(per_filt,bins=R)
    P_A=P_A/P_A.sum()
    A=[]
    for i in range(R):    
        A.append(0.5*(bin_edges[i]+bin_edges[i+1]))
    c_out=[]
    for i in range(len(P_A)):
        if P_A[i]>0:   
            c_out.append([A[i],P_A[i]])
    c_out=np.array(c_out)
    A=c_out[:,0]
    P_A=c_out[:,1]
    for i in range(len(A)):
        if A[i]>A_min:
            A_min = i
            break
    for i in range(len(A)):
        if A[i]>A_max:
            A_max = i-1
            break
    P_A=P_A[A_min:A_max]
    A=A[A_min:A_max]
    fit=np.polyfit(np.log(A),np.log(P_A),1)
    kessi = (fit[0]*-2)
    if plot==1:
        print('kessi = ',kessi)
        plt.scatter(A,P_A)
        plt.xlabel('A')
        plt.ylabel('P(A)')
        plt.show()
        yA=np.exp((np.log(A)*fit[0])+fit[1])
        plt.loglog(A,P_A,'o')
        plt.loglog(A,yA,label='slope='+str(fit[0]))
        plt.xlabel('A')
        plt.ylabel('P(A) (log scale)')
        plt.show()
        plt.plot(A,P_A*(A**(kessi/2)),'.')
        plt.xlabel('A')
        plt.ylabel('P(A)*s**(kessi/2)')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(0.001,1000)
        plt.show()
    return kessi 

#for row axis=0 - n is number of row - clean coloumn with zero in n_th row 
def sorter(data,axis,n,clean):
    if axis==0:
        data=data[:, data[n].argsort()]
        if clean==1:
            data=data[:,data[n,:]!=0]
    if axis==1:
        data=data[data[:,n].argsort(),:]
        if clean==1:
            data=data[data[:,n]!=0,:]
    return data
def Gr(data,th,points):
    import random
    data -= data.mean()
    data /= data.std()
    contour = measure.find_contours(data,th)
    l_c=len(contour)    
    l=[]
    for i in range(len(contour)):
        l.append([i,len(contour[i])])
    l=np.array(l)
    n_co=sorter(l,1,1,0)[:,0]
    nor=0   
    print('find ',l_c ,'conours and select Maximum',points,'points of any contour (;')
    per=0
    o=[]
    delta=1
    for k in (n_co):
        k=int(k)
        print("Progress {:2.1%}".format( (nor+1) / l_c), end="\r")
        R_max=0
        c=np.zeros(1)
        n=0
        len_contour=int(len(contour[k]))           
        n_c1=random.sample(range(len_contour),min(len_contour,points)) 
        nor+=1
        d=contour[k][[n_c1]]
        len_d=int(len(d))
        for i in range(len_d):
            for j in range(i+1,len_d):
                x=d[i,0]-d[j,0]
                y=d[i,1]-d[j,1]
                r=int(np.sqrt((x**2)+(y**2))/delta)
                R_max=max(r,R_max)
                c.resize(R_max+1)
                c[r]=c[r]+1
                n+=1    
        C=c/n
        G=[]
        for i in range(len(C)):
            if C[i]>0:
                G.append([i,C[i]])
        o.append(np.array(G))
    m=0
    for i in range(len(o)):
        if len(o[i])>0:           
            m=max(m,max(o[i][:,0]))    
    m=int(m)
    Cc=np.zeros((m+1,3))
    Cc[:,0]=np.arange(0,m+1,1)
    for i in range(len(o)):
        for j in range(len(o[i])):
            Cc[int(o[i][j,0]),1]+=o[i][j,1]
            Cc[int(o[i][j,0]),2]+=1
    Cc=Cc[Cc[:,2]!=0]       
    Cc[:,1]/=Cc[:,2]
    Cc=Cc[:,:2]
    Cc[:,1]/=Cc[:,1].sum()
    return Cc

def xl(data,th,points,r_min,r_max,plot):
    rr2=Gr(data,th,points)
#     rr2[:,1]=g[:,1]
    x=np.log(rr2[r_min:r_max,0])
    y=np.log(rr2[r_min:r_max,1])
    fit=np.polyfit(x,y,1)
    if plot==1:
        print('-2xl = ',fit[0])
        r=rr2[10:,0]
        gg=rr2[10:,1]
        #plt.ylim(-100,10)
        grr=gg*(r**(-1*fit[0]))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel('G(r)*r**(2xl)')
        plt.plot(r,grr,'r--')
    return (fit[0]/(-2))

#0_1 boxing
def zobox(data,th):
    [contour,perimeter,area,radius] = carachters(data,th)
    row=data.shape[0]
    col=data.shape[1]
    pos_con=np.zeros((row,col))
    for i in range(len(contour)):
        for j in range(len(contour[i])):
            pos_con[int(contour[i][j,0]),int(contour[i][j,1])]+=1
    return pos_con

def N_l(data,th):
    pos_con=zobox(data,th)
    len_d=len(data)
    N=np.zeros((int(np.log2(len_d))+1,2))
    for i in range(int(np.log2(len_d))+1):    
        nn=2**i
        len_d=len(data)
        sl=int(len_d/nn)
        N[i,0]=1/nn
        for j in range(nn):  
            for k in range(nn):
                if np.mean(pos_con[j*sl:(j+1)*sl,k*sl:(k+1)*sl])>0:
                    N[i,1]+=1
    return N
def d_d(data,th,plot):
    N=N_l(data,th)
    x=np.log(N[1:,0])
    y=np.log(N[1:,1])    
    fi2 = np.polyfit(x,y,1)
    if plot==1:
        print('d = ',fi2[0]*-1)
        plt.plot(x,y)
        plt.plot(x,x*(fi2[0])+fi2[1],'r--')
        plt.xlabel('log(l)')
        plt.ylabel('log N(l)')
        plt.show()
    return (fi2[0]*-1)

def p_l(data,th):
    pos_con=zobox(data,th)
    len_d=len(data)
    #N=np.zeros((int(np.log2(len_d))+1,2))
    N=[]
    for i in range(int(np.log2(len_d))+1):    
        nn=2**i
        len_d=len(data)
        sl=int(len_d/nn)        
        x=0        
        zp=0
        p=np.zeros((nn,nn))
        for j in range(nn):
            y=0            
            for k in range(nn):                
                s=pos_con[j*sl:(j+1)*sl,k*sl:(k+1)*sl].sum()
                zp=zp+s
                p[x,y]=s
                y=y+1
            x=x+1
        u=[]
        u.append([1/nn,p/zp])
        N.append(u)
    return N
# example output : p[3][0][1].sum() is 1

#D_q
def D_q(data,th,q1,q2,qq):
    k=0
    N=p_l(data,th)
    z_q=np.zeros((qq,len(N)-1,2))
    D_q=np.zeros((qq,len(N)-1,2)) 
    for q in np.linspace(q1,q2,qq):
        if q<=0:
            print(-1)
            for i in range(1,len(N)):
                yy=0
                for z in range(len(N[i][0][1])):
                    for e in range(len(N[i][0][1])):
                        if N[i][0][1][z,e]>0:
                            yy=yy+(N[i][0][1][z,e])**q 
                z_q[k,i-1,1]=yy       
                l=N[i][0][0]
                D_q[k,i-1,1]=(1/(q-1))*np.log(z_q[k,i-1,1])/np.log(l)
                z_q[k,i-1,0]=l
                D_q[k,i-1,0]=l

        if q==1:
            print(1)
            for i in range(1,len(N)):
                l=N[i][0][0]
                yy=0
                for z in range(len(N[i][0][1])):
                    for e in range(len(N[i][0][1])):
                        if N[i][0][1][z,e]>0:
                            yy=yy+(N[i][0][1][z,e])*np.log(N[i][0][1][z,e])  
                D_q[k,i-1,1]=yy/np.log(l)
                z_q[k,i-1,0]=l
                D_q[k,i-1,0]=l            
        if q>0 and q!=1:
            print(11)
            for i in range(1,len(N)):
                yy=0
                for z in range(len(N[i][0][1])):
                    for e in range(len(N[i][0][1])):
                        yy=yy+(N[i][0][1][z,e])**q 
                        #print(N[i][0][1][z,e])
                z_q[k,i-1,1]=yy       
                l=N[i][0][0]
                D_q[k,i-1,1]=(1/(q-1))*np.log(z_q[k,i-1,1])/np.log(l)
                z_q[k,i-1,0]=l
                D_q[k,i-1,0]=l
        k=k+1
    return D_q
