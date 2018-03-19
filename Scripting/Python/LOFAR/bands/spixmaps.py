#!/usr/bin/python
#====================================================
#This script is used to create a spectral index map
#by fitting pixelwise a powerlaw function
#====================================================

#====================================================
#          Loading necessary packages
#====================================================

import astropy.io.fits as pyfits
import os
import numpy as np
from astropy.time import Time
import pylab as plt
import scipy
import scipy.optimize as opt
import matplotlib.colors as colors
import matplotlib
from matplotlib.mlab import bivariate_normal


#====================================================
#     Loading absolute and temporary variables 
#====================================================

path_out="spixmap.fits"
fromx=90
tox=160
fromy=90
toy=160

#=====================================================
#             Definitions for functions
#=====================================================

def extrapolate(matrix,y,x): #Function to correct dead pixels by extrapolating from N and NN neighbours
	matrix[y,x]=(((matrix[y+1,x]+matrix[y,x+1]+matrix[y-1,x]+matrix[y,x-1])/4)*2**0.5+((matrix[y+1,x+1]+matrix[y-1,x+1]+matrix[y+1,x-1]+matrix[y-1,x-1])/4))/(1+2**0.5)

	return matrix

def dist(x, y, xx, yy):
	x = x * 1.
	y = y * 1.
	xx = xx * 1.
	yy = yy * 1.

	return np.sqrt((yy-y)**2+(xx-x)**2) 

def median(a): # a: array of matrices
	if len(a) == 0:
		print "ERROR: No Data in Array!!!"
		return False
	y_len=np.shape(a[0])[0]
	x_len=np.shape(a[0])[1]

	for i in range(len(a)):
		test= (np.shape(a[0]) == np.shape(a[i]))
		if test == False:
			print "ERROR: Sizes are different!!!"
			return False

	c=np.zeros((y_len,x_len))
	for x in range(x_len):
		for y in range(y_len):
			vals=[]
			for i in range(len(a)):
				vals.append(a[i][y,x])
			c[y,x]= np.median(vals) * 1.
	
	return c

def linear1d((x),m,t):
    m = float(m)
    t = float(t)
    y = m*x+t
    return y

def powerlaw((x),a,k):
    a = float(a)
    k = float(k)
    y = a*x**(k)
    return y 


def twoD_Gaussian((x, y), amplitude,amplitude2, xo, yo,dx,dy, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    dx = float(dx)
    yo = float(yo)
    dy = float(dy) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo+dx)**2) + 2*bb*(x-xo+dx)*(y-yo+dy) 
                            + cc*((y-yo+dy)**2)))
    return g.ravel()

def twoD_Gaussian_norav((x, y), amplitude,amplitude2, xo, yo,dx,dy, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    dx = float(dx)
    yo = float(yo)
    dy = float(dy)  
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo+dx)**2) + 2*bb*(x-xo+dx)*(y-yo+dy) 
                            + cc*((y-yo+dy)**2)))
    return g

def twoD_Gaussian_norav1((x, y), amplitude,amplitude2, xo, yo,dx,dy, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    dx = float(dx)
    yo = float(yo)
    dy = float(dy) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g

def twoD_Gaussian_norav2((x, y), amplitude,amplitude2, xo, yo,dx,dy, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    dx = float(dx)
    yo = float(yo)
    dy = float(dy)  
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude2*np.exp( - (aa*((x-xo+dx)**2) + 2*bb*(x-xo+dx)*(y-yo+dy) 
                            + cc*((y-yo+dy)**2)))
    return g



#=====================================================
#                      Main
#=====================================================

#read in all your data:

fitslist_paths="fitslist.txt"
fitslist=list(open(fitslist_paths,"r"))
n_fitslist=len(fitslist)

datacubes=[]
headercubes=[]

for i in range(n_fitslist):  #Read-in the Fits-Files and sort images and headers in lists
    path=fitslist[i]
    path=path.rstrip()
    hdulist=pyfits.open(path) * 1
    dataset=hdulist[0].data * 1
    header=hdulist[0].header
    datacubes.append(dataset)
    headercubes.append(header)

datacubes_np=np.array(datacubes)


#build the output grid
x_size=headercubes[0]["NAXIS1"]
y_size=headercubes[0]["NAXIS2"]
xx=np.linspace(0,x_size-1,x_size)
yy=np.linspace(0,y_size-1,y_size)
spixmap=np.meshgrid(xx, yy)[0]
spixmap=spixmap * 0.0 + 20000.0





freqs=[]
for vals in range(datacubes_np.shape[0]): #build list of central frequencies pulled form the headers
    freqs.append(headercubes[vals]['CRVAL3'])

for y in range(datacubes_np.shape[1]): #go from pixel pos to pixel pos and fit powerlaw by freq.
    if (fromy<y<toy):
        for x in range(datacubes_np.shape[2]):
            if (fromx<x<tox):
                values=[]
                linvals=[]
                yvals=[]
                for vals in range(datacubes_np.shape[0]):
                    if datacubes_np[vals,y,x]>0:
                       values.append(datacubes_np[vals,y,x])
                       linvals.append(np.log10(datacubes_np[vals,y,x]))
                       yvals.append(freqs[vals])
                if len(yvals)>2:
                   popt, pcov = opt.curve_fit(powerlaw, yvals, values, maxfev=200000)
                   #popt, pcov = opt.curve_fit(powerlaw, freqs, values, maxfev=200000)
                   print y,x,popt[1]
                   print values
                   print linvals
                   spixmap[y,x] = popt[1]
                else:
                   spixmap[y,x] = 19000.0

if os.path.isfile(path_out): #Check if data already exists if yes, remove
   os.system("rm " + path_out)

    
pyfits.writeto(path_out, spixmap, headercubes[0])

#return 0






# hdulist=pyfits.open(path_in) * 1
# dataset=hdulist[0].data
# header=hdulist[0].header

# datasetr=dataset * 1
# datasetr=datasetr[0][0] * 1
# dataset=dataset.ravel()


# # We start to calculate median flats for before/after measurement

# x_size=header["NAXIS1"]
# y_size=header["NAXIS2"]
# freq=header["CRVAL3"]

# #-------------------------------------
# #build grid
# thet=(thet/360.)*scipy.pi
# thet_min=(thet_min/360.)*scipy.pi
# thet_max=(thet_max/360.)*scipy.pi


# x=np.linspace(0,x_size-1,x_size)
# y=np.linspace(0,y_size-1,y_size)
# x, y=np.meshgrid(x, y)

# #popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,p0=(amp,ampo,x_max,y_max,sigm1,sigm1o,sigm2,sigm2o,thet,tval))

# popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,maxfev=100000,p0=(amp,ampo,x_max,y_max,dx_g,dy_g,sigm1,sigm1o,sigm2,sigm2o,thet,tval),bounds=([amp_min,ampo_min, x_max_min, y_max_min,dx_min, dy_min, sigm1_min,sigm1o_min, sigm2_min,sigm2o_min, thet_min, tval_min],[amp_max,ampo_max, x_max_max, y_max_max, dx_max, dy_max, sigm1_max,sigm1o_max, sigm2_max,sigm2o_max, thet_max, tval_max]))

# data_fitted = twoD_Gaussian_norav((x, y), *popt)
# data_fitted1= twoD_Gaussian_norav1((x, y), *popt)
# data_fitted2= twoD_Gaussian_norav2((x, y), *popt)


# resid_img=datasetr-data_fitted1
# resid_img2=datasetr-data_fitted

# print popt
# print pcov
# plt.figure("Fitfunction"); plt.imshow(data_fitted,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
# plt.figure("Fitfunction1"); plt.imshow(data_fitted1,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
# plt.figure("Fitfunction2"); plt.imshow(data_fitted2,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
# plt.figure("Initial Data"); plt.imshow(datasetr,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
# plt.figure("Residualfull"); plt.imshow(resid_img2,cmap="magma",norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
# plt.figure("Residual"); plt.imshow(resid_img,cmap="magma",norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis();plt.show()

# if os.path.isfile(path_out):
# 	os.system("rm " + path_out)
# if os.path.isfile(path_out_comp):
# 	os.system("rm " + path_out_comp)
	
# pyfits.writeto(path_out, resid_img, header)
# pyfits.writeto(path_out_comp, data_fitted1, header)










