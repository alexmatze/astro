#!/usr/bin/python
#====================================================
#This script is used to successfully subtract gaussian
#components of a fits-file image
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

path_in="0836+710_stacked_hba.fits"
path_in="./bands/126/band7.fits"
path_out="0836+710_resid1.fits"
path_out_comp="0836+710_comp1.fits"
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

#read in all your file:

hdulist=pyfits.open(path_in) * 1
dataset=hdulist[0].data
header=hdulist[0].header

dataset_norav=dataset * 1

datasetr=dataset * 1
dataset=dataset.ravel()

x_size=header["NAXIS1"]
y_size=header["NAXIS2"]

### Parameters for Modelcomponent
x_max=127 #x position of peak
y_max=128 #y position of peak

thet=-38.374932908286 #position angle of 2d gauss in deg
FWHM1=0.00013741709347007 #fwhm of first 1d gauss in deg
FWHM2=9.5506462157887E-05 #fwhm of second 1d gauss in deg
degperpix=2.6801750408464E-05 #degree per pixel
tval=0 #Value offset for fit (noise)
amp=1.88210 #Amplitude of max position
#-------------------------------------
#apply corrections to meet present definitions
thet=thet - 90. #correct position angle for needed definition
FWHM1=FWHM1/degperpix #converse FWHM1 from deg to pix
FWHM2=FWHM2/degperpix #converse FWHM2 from deg to pix
sigm1=FWHM1/2.355 #converse from FWHM1 to 1sigma
sigm2=FWHM2/2.355 #converse from FWHM2 to 1sigma
#build grid
thet=(thet/360.)*scipy.pi

#modelimg=dataset *0.
x=np.linspace(0,x_size-1,x_size)
y=np.linspace(0,y_size-1,y_size)
x, y=np.meshgrid(x, y)

modelimg=twoD_Gaussian_norav((x,y),amp,x_max,y_max,sigm1,sigm2,thet,tval)

plt.figure("Modelcomponent"); plt.imshow(modelimg,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis();plt.show()

resid=dataset_norav - modelimg

path_out="modeltest.fits"
if os.path.isfile(path_out):
	os.system("rm " + path_out)
#if os.path.isfile(path_out_comp):
#	os.system("rm " + path_out_comp)
	
pyfits.writeto(path_out, resid, header)
#pyfits.writeto(path_out_comp, dataset_out_comp, header)










