from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from skimage.draw import (line, polygon,polygon_perimeter,
                          ellipse, ellipse_perimeter,
                          circle_perimeter, bezier_curve)
try:
    from skimage.draw import disk as circle
except:
    from skimage.draw import circle
def random_inside(img,num=1,margin=0):
    '''
    This function will return a random positon in an image.
    '''
    lx,ly = img.shape[:2]
    x = np.random.randint(margin,lx-margin,num)
    y = np.random.randint(margin,ly-margin,num)
    return np.c_[x,y]

def draw_line(imgx,begins=None,ends=None,num=10,value=1):
    if begins is None or ends is None:
        begins = random_inside(imgx,num=num)
        ends = random_inside(imgx,num=num)
    else:
        num = min(len(begins),len(ends))
    for i in range(num):
        rr, cc = line(begins[i][0],begins[i][1],ends[i][0],ends[i][1])
        imgx[rr, cc] += value
    return imgx

def draw_circle(img,centers=None,rs=None,num=10,rmin=10,rmax=20,fill=False,value=1):
    if fill:
        drawer = circle
    else:
        drawer = circle_perimeter
    if centers is None or rs is None:
        rs = np.random.randint(rmin,rmax,num)
        centers = random_inside(img,num=num,margin=rmax)
    else:
        num = min(len(centers),len(rs))
    for i in range(num):
        rr, cc = drawer(centers[i][0], centers[i][1], rs[i], shape=img.shape)
        img[rr, cc] += value
    return img

def draw_polygon(img,polys,fill=False,value=1):
    if fill:
        drawer = polygon
    else:
        drawer = polygon_perimeter
    for poly in polys:
        rr, cc = drawer(poly[:, 0], poly[:, 1],shape=img.shape)
        img[rr, cc] += value
    return img

def draw_square(img,centers,lxs,lys,fill=False,value=1):  
    polys = []
    num = len(centers)
    for i in range(num):
        exes = np.array([[-lxs[i],-lys[i]],[-lxs[i],lys[i]],[lxs[i],lys[i]],[lxs[i],-lys[i]]])
        polys.append([centers[i]+ex for ex in exes])
    polys = np.array(polys)        
    img = draw_polygon(img,polys=polys,fill=fill,value=value)
    return img

def draw_equitri(img,centers,rs,fill=False,value=1):  
    cc = np.sqrt(3.)/2.
    polys = []
    num = len(centers)
    for i in range(num):
        exes = np.array([[-rs[i],0],[0.5*rs[i],-cc*rs[i]],[0.5*rs[i],cc*rs[i]]])
        polys.append([centers[i]+ex for ex in exes])
    polys = np.array(polys)        
    img = draw_polygon(img,polys=polys,fill=fill,value=value)
    return img
    
    
# ellipse
# rr, cc = ellipse(300, 300, 100, 200, img.shape)
# img[rr, cc, 2] = 1
# ellipse
# rr, cc = ellipse_perimeter(120, 400, 60, 20, orientation=math.pi / 4.)
# img[rr, cc, :] = (1, 0, 1)
# # Bezier curve
# rr, cc = bezier_curve(70, 100, 10, 10, 150, 100, 1)
# img[rr, cc, :] = (1, 0, 0)
