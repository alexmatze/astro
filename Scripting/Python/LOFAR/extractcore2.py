#!/usr/bin/python
#====================================================
#This script is used to successfully subtract gaussian
#components of a fits-file image
#====================================================

#====================================================
#          Loading necessary packages
#====================================================

import pyfits
import os
import numpy as np
from astropy.time import Time
import pylab as plt
import scipy
import scipy.optimize as opt

#====================================================
#     Loading absolute and temporary variables 
#====================================================

path_in="0836+710_stacked_hba.fits"
path_out="0836+710_resid.fits"
path_out_comp="0836+710_comp.fits"
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

def twoD_Gaussian((x, y), amplitude,amplitude2, xo, yo, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    yo = float(yo) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo)**2) + 2*bb*(x-xo)*(y-yo) 
                            + cc*((y-yo)**2)))
    return g.ravel()

def twoD_Gaussian_norav((x, y), amplitude,amplitude2, xo, yo, sigma_x,sigma_xo, sigma_y,sigma_yo, theta, offset):
    xo = float(xo)
    yo = float(yo) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo)**2) + 2*bb*(x-xo)*(y-yo) 
                            + cc*((y-yo)**2)))
    return g



#=====================================================
#                      Main
#=====================================================

#read in all your darks and flats:

hdulist=pyfits.open(path_in) * 1
dataset=hdulist[0].data
header=hdulist[0].header

datasetr=dataset * 1
dataset=dataset.ravel()


# We start to calculate median flats for before/after measurement

x_size=header["NAXIS1"]
y_size=header["NAXIS2"]

### Fitting Start Parameters
x_max=128 #x position of peak
y_max=129 #y position of peak
xo_max=128 #x position of second peak
yo_max=129 #y position of second peak


thet=45 #position angle of both 2d gauss in deg

amp=6 #Amplitude of max position
ampo=2 #Amplitude of max position
sigm1=6 #width of first 1d gauss
sigm2=4 #width of second 1d gauss
sigm1o=8 #width of first 1d gauss
sigm2o=6 #width of second 1d gauss

tval=0 #Value offset for fit (noise)
#-------------------------------------
#build grid
thet=(thet/360.)*scipy.pi


x=np.linspace(0,x_size-1,x_size)
y=np.linspace(0,y_size-1,y_size)
x, y=np.meshgrid(x, y)

popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,p0=(amp,ampo,x_max,y_max,sigm1,sigm1o,sigm2,sigm2o,thet,tval))
data_fitted = twoD_Gaussian_norav((x, y), *popt)

resid_img=datasetr-data_fitted

print popt
print pcov
plt.figure(); plt.imshow(data_fitted);plt.colorbar();plt.logspace();plt.gca().invert_yaxis();plt.show()
#plt.figure(); plt.imshow(datasetr);plt.colorbar();plt.gca().invert_yaxis();plt.show()
#plt.figure(); plt.imshow(resid_img);plt.colorbar();plt.gca().invert_yaxis();plt.show()
if os.path.isfile(path_out):
	os.system("rm " + path_out)
if os.path.isfile(path_out_comp):
	os.system("rm " + path_out_comp)
	
#pyfits.writeto(path_out, dataset_out, header)
#pyfits.writeto(path_out_comp, dataset_out_comp, header)










