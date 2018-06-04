import numpy as np
import routines as myr

def fortranize(m):
	return np.array(m,order='F')
	
	
	
	
	
def cross_correlarion_fucntion(m1,m2,n_p,mask_value=np.nan):
#	(npixr,mean,var) = myr.make_standard(m1,0)
#	(npixr,mean,var) = myr.make_standard(m2,0)
	m1 = fortranize(m1)
	m2 = fortranize(m2)
	(cor,vcor) = myr.cross_cf(m1,m2,n_p,mask_value)
	return (cor,vcor)
	
def correlarion_fucntion(m,n_p,mask_value=np.nan):
#	(npixr,mean,var) = myr.make_standard(m,0)
	m = fortranize(m)
	(cor,vcor) = myr.cross_cf(m,m,n_p,mask_value)
	return (cor,vcor)
	
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
