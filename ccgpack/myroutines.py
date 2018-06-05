import os
import cv2
import time
import ctypes
import numpy as np
import routines as myr

def fortranize(m):
	return np.array(m,order='F')

def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
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

def ccf(m,n_p):
	(npixr,mean,var) = myr.make_standard(m,0)
	m = fortranize(m)
	(cor,vcor) = myr.cross_cf(m,m,n_p,0)
	return cor


def ppcf(m,th,nfmax,rmax):
	lg = m.shape[0]
	(npixr,mean,var) = myr.make_standard(m,0)
	(nf1,fl1) = myr.peak_f(m,th,100)

	if nf1>nfmax:
		rlist = np.random.randint(0,nf1,size=nfmax)
		nf1 = nfmax
		fl1 = fl1[:,rlist]

	fl1 = fortranize(fl1[:,:nf1])
	(ksi,vksi) = myr.ffcf(1,lg,fl1,fl1,5*nf1,1,rmax)
	return ksi

def cccf(m,th,nfmax,rmax):
	lg = m.shape[0]
	(npixr,mean,var) = myr.make_standard(m,0)
	(nf1,nf2,fl1,fl2) = myr.cross_f(m,th,th+0.01)

	if nf1>nfmax:
		rlist = np.random.randint(0,nf1,size=nfmax)
		nf1 = nfmax
		fl1 = fl1[:,rlist]
	if nf2>nfmax:
		rlist = np.random.randint(0,nf2,size=nfmax)
		nf2 = nfmax
		fl2 = fl2[:,rlist]

	fl1 = fortranize(fl1[:,:nf1])
	fl2 = fortranize(fl2[:,:nf2])
	(ksi1,vksi) = myr.ffcf(1,lg,fl1,fl1,5*nf1,1,rmax)
	(ksi2,vksi) = myr.ffcf(1,lg,fl2,fl2,5*nf2,1,rmax)
	ksi = (ksi1+ksi2)/2
	return ksi

def pccf(m,th,nfmax,rmax):
	lg = m.shape[0]
	(npixr,mean,var) = myr.make_standard(m,0)
	(nf1,fl1) = myr.peak_f(m,th,100)
	(nf2,nf3,fl2,fl3) = myr.cross_f(m,th,th+0.01)

	if nf1>nfmax:
		rlist = np.random.randint(0,nf1,size=nfmax)
		nf1 = nfmax
		fl1 = fl1[:,rlist]
	if nf2>nfmax:
		rlist = np.random.randint(0,nf2,size=nfmax)
		nf2 = nfmax
		fl2 = fl2[:,rlist]
	if nf3>nfmax:
		rlist = np.random.randint(0,nf3,size=nfmax)
		nf3 = nfmax
		fl3 = fl3[:,rlist]

	fl1 = fortranize(fl1[:,:nf1])
	fl2 = fortranize(fl2[:,:nf2])
	fl3 = fortranize(fl3[:,:nf3])
	(ksi1,vksi) = myr.ffcf(1,lg,fl1,fl2,5*nf1,1,rmax)
	(ksi2,vksi) = myr.ffcf(1,lg,fl1,fl3,5*nf2,1,rmax)
	ksi = (ksi1+ksi2)/2
	return ksi

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
#	strg = np.where(strg<2*mstrg,1/2,0)
	ornt = np.arctan2(dx[:,:ly],dy[:lx,:])
	return [strg,ornt]

def canny(d,R,meth,edd):

	if (R!=0):
		dt = np.fft.fft2(d)
		if meth=='g':
			for i in range(sz):
				for j in range(sz):
					k2 = 1.*(i*i+j*j)/d.shape[0]
					dt[i,j]=dt[i,j]*np.exp(-k2*R*R/2)

		if meth=='tp':
			for i in range(sz):
				for j in range(sz):
					k = np.sqrt(0.001+i*i+j*j)/sz
					dt[i,j]=dt[i,j]* 3*(np.sin(k*R)-k*R*np.cos(k*R))/(k*R)**3

		d = np.fft.ifft2(dt)
		d = abs(d)

	if edd=='lap':
		d = cv2.Laplacian(d,cv2.CV_64F)

	if edd=='sob':
		sobelx = cv2.Sobel(d,cv2.CV_64F,1,0,ksize=3)
		sobely = cv2.Sobel(d,cv2.CV_64F,0,1,ksize=3)
		d =np.sqrt(sobelx**2+sobely**2)

	if edd=='sch':
		scharrx = cv2.Scharr(d,cv2.CV_64F,1,0)
		scharry = cv2.Scharr(d,cv2.CV_64F,0,1)
		d =np.sqrt(scharrx**2+scharry**2)

#	do = (d-np.min(d))/(np.max(d)-np.min(d))

#	dt = np.where((do < tl) | (do > tu), 0, 1)

#	if False:
#		for i in range(1,dt.shape[0]-1): 
#			for j in range(1,dt.shape[1]-1): 
#				if dt[i,j]==1:
#					if (dt[i+1,j]==0) and dt[i,j+1]==0 and dt[i,j-1]==0 and dt[i-1,j]==0:# and (dt[i+1,j+1]==0) and (dt[i-1,j+1]==0) and (dt[i+1,j-1]==0) and (dt[i-1,j-1]==0):
#						dt[i,j]==1

	return d

def badcurvelet(name,c1,c2):
	np.savetxt('curvelet_input',name)
	crv_path='/home/gf/work/pakages/curvelet'
	crn_path=os.getcwd()
	os.system("cat <<EOF | matlab -nodesktop -nosplash -nodisplay\n"+"mn="+str(c1)+";\n"+"mx="+str(c2)+";\n"+"CC = dlmread('curvelet_input');\n"+"cd "+crv_path+";\n"+"C = fdct_wrapping(CC);\n"+"CC=C;\n"+"for m=mn:mx\n"+"     C=CC;\n"+"%Making other components zero\n"+"    for s=1:length(C)\n"+"      for w=1:length(C{s})\n"+" C{s}{w}=C{s}{w}.*(s==m);\n"+"      end\n"+"    end\n"+"    y=ifdct_wrapping(C);\n"+"    out =['c' int2str(m) '_' 'curvelet_input'];\n"+"    y=real(y);\n"+"    cd "+crn_path+";\n"+"    dlmwrite(out,y,' ');\n"+"    cd "+crv_path+";\n"+"end % first\n"+"exit\n"+"EOF\n")

	print '\n'
	res = {i:None for i in range(c1,c2)}
	for i in range(c1,c2+1):
		res[i]=np.loadtxt('c'+str(i)+'_curvelet_input')
		os.remove('c'+str(i)+'_curvelet_input')
	return res

def imshow(ax,strg,tlt,rotation=False):
	im = ax.imshow(strg, cmap='spectral')
# Titling
	ax.set_title(tlt, y=1.04)
	if rotation:
		ax.set_title(ttl, rotation='vertical',x=-0.1,y=0.5)

	ax.axis('off')
#	ax1.set_xticks([])
#	ax1.set_yticks([])

def plot(ax,x,y,tlt,
clrs=['b','r','lime','c','indigo','gold','plum','k'],
xlab=[False,''],ylab=[False,''],
logx=False,logy=False,xlim=[False,0,1],
ylim=[False,0,1]):
# x, y and tlt must be inputed as list of plots
	for i in range(len(y)):
		ax.plot(x[i],y[i],clrs[i],label=tlt[i])

	if xlab[0]:
		ax.set_xlabel(xlab[1],fontsize=25)
	if ylab[0]:
		ax.set_ylabel(ylab[1],fontsize=25)
	if logx:
		ax.set_xscale("log", nonposx='clip')
	if logy:
		ax.set_yscale("log", nonposx='clip')
	if xlim[0]:
		ax.set_xlim(xlim[1],xlim[2])
	if ylim[0]:
		ax.set_ylim(ylim[1],ylim[2])
	
	ax.tick_params(axis='both', which='major', labelsize=18)

def histax(ax,m,n_c,tlt=[False,''],fc='b',alpha=0.75):
	n, bins, patches = ax.hist(m.flatten('F'), n_c, facecolor=fc,alpha=alpha)
#	bin_centers = 0.05 * np.diff(bins) + bins[:-1]
#	for count, x in zip(n, bin_centers):
		  # Label the raw counts
#		  ax1.annotate(str(count), xy=(x, 0), xycoords=('data', 'axes fraction'),
#		      xytext=(0, -18), textcoords='offset points', va='top', ha='center')

	ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	tks=[round(x,1) for x in np.linspace(np.min(m),np.max(m),5)]
	ax.set_xticks(tks)
	ax.tick_params(axis='both', which='major', labelsize=18)
	if tlt[0]:
		ax.set_title(tlt[1], fontsize=30,x=0.5,y=1.1)

def pdf_ax(ax,m,n_c,color='b',tlt=[False,'']):
	n, bins = np.histogram(m.flatten('F'), n_c)
	dx = (bins[1]-bins[0])/2
	bins = bins[:-1]+dx
	nt = 0.0001*np.size(m)

	bins = bins[n>nt]
	n = n[n>nt]

	ax.plot(bins,n,color=color)

	ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	tks=[round(x,1) for x in np.linspace(bins[0],bins[-1],5)]
	ax.set_xticks(tks)
	ax.tick_params(axis='both', which='major', labelsize=18)
	if tlt[0]:
		ax.set_title(tlt[1], fontsize=30,x=0.5,y=1.1)

def pdf(m,n_c):
	m = np.array(m)
	n, bins = np.histogram(m.flatten('F'), n_c)
	dx = (bins[1]-bins[0])/2
	bins = bins[:-1]+dx
#	nt = np.size(m)

#	bins = bins[n>nt]
#	n = n[n>nt]
	return [bins,n]

def curvelet(m,n_scales,r_scale,n_wedges,ac):
	dir_path = os.path.dirname(os.path.realpath(__file__))
	curlib = ctypes.cdll.LoadLibrary(dir_path+"/curvelet.so")
	m = np.array(m, dtype=np.double)
	nx = m.shape[0]
	ny = m.shape[1]

	aptr = m.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))

	curlib.curvelet(aptr,nx,ny,n_scales,r_scale-1,n_wedges,ac)
	return m

#def d_pdf(ax,m1,m2,n_c,tlt=[False,'']):
#	n1, bins1 = np.histogram(m1.flatten('F'), n_c)
#	dx1 = (bins1[1]-bins1[0])/2

#	n2, bins2 = np.histogram(m2.flatten('F'), n_c)
#	dx2 = (bins2[1]-bins2[0])/2

#	ax.plot(bins1[:-1]+dx1,n1)
#	ax.plot(bins2[:-1]+dx2,n2)

#	ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#	tks=[round(x,1) for x in np.linspace(np.min(m),np.max(m),5)]
#	ax.set_xticks(tks)
#	ax.tick_params(axis='both', which='major', labelsize=18)
#	if tlt[0]:
#		ax.set_title(tlt[1], fontsize=30,x=0.5,y=1.1)

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

def stat_describe(m,m_max):
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



