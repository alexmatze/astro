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
import matplotlib.colors as colors
import matplotlib
from matplotlib.mlab import bivariate_normal

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

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def twoD_Gaussian_norav((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
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

x_1=0 #x start position of fitting box
y_1=0 #y start position of fitting box
x_2=0 #x end position of fitting box
y_2=0 #y end position of fitting box
thet=45 #position angle of 2d gauss in deg
sigm1=6 #width of first 1d gauss
sigm2=4 #width of second 1d gauss
tval=0 #Value offset for fit (noise)
amp=8 #Amplitude of max position
#-------------------------------------
#build grid
thet=(thet/360.)*scipy.pi


x=np.linspace(0,x_size-1,x_size)
y=np.linspace(0,y_size-1,y_size)
x, y=np.meshgrid(x, y)

popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,p0=(amp,x_max,y_max,sigm1,sigm2,thet,tval))
data_fitted = twoD_Gaussian_norav((x, y), *popt)

resid_img=datasetr-data_fitted

Z1 = bivariate_normal(x,y, 0.1, 0.2, 1.0, 1.0) + 0.1 * bivariate_normal(x, y, 1.0, 1.0, 0.0, 0.0)
print popt
print pcov

#fig, ax = plt.subplots(1, 1)
plt.figure("Fitfunction"); plt.imshow(data_fitted,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
plt.figure("Initial Data"); plt.imshow(datasetr,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
plt.figure("Residual"); plt.imshow(resid_img,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis();plt.show()
if os.path.isfile(path_out):
	os.system("rm " + path_out)
if os.path.isfile(path_out_comp):
	os.system("rm " + path_out_comp)
	
#pyfits.writeto(path_out, dataset_out, header)
#pyfits.writeto(path_out_comp, dataset_out_comp, header)










