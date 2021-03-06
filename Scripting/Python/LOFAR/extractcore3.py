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
path_out="0836+710_resid3.fits"
path_out_comp="0836+710_comp3.fits"
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

def twoD_Gaussian((x, y), amplitude,amplitude2,amplitude3, xo, yo, sigma_x,sigma_xo,sigma_xoo, sigma_y,sigma_yo,sigma_yoo, theta, offset):
    xo = float(xo)
    yo = float(yo) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    aaa = (np.cos(theta)**2)/(2*sigma_xoo**2) + (np.sin(theta)**2)/(2*sigma_yoo**2)
    bbb = -(np.sin(2*theta))/(4*sigma_xoo**2) + (np.sin(2*theta))/(4*sigma_yoo**2)
    ccc = (np.sin(theta)**2)/(2*sigma_xoo**2) + (np.cos(theta)**2)/(2*sigma_yoo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo)**2) + 2*bb*(x-xo)*(y-yo) 
                            + cc*((y-yo)**2)))+ amplitude3*np.exp( - (aaa*((x-xo)**2) + 2*bbb*(x-xo)*(y-yo) 
                            + ccc*((y-yo)**2)))
    return g.ravel()

def twoD_Gaussian_norav((x, y), amplitude,amplitude2,amplitude3, xo, yo, sigma_x,sigma_xo,sigma_xoo, sigma_y,sigma_yo,sigma_yoo, theta, offset):
    xo = float(xo)
    yo = float(yo) 
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    aa = (np.cos(theta)**2)/(2*sigma_xo**2) + (np.sin(theta)**2)/(2*sigma_yo**2)
    bb = -(np.sin(2*theta))/(4*sigma_xo**2) + (np.sin(2*theta))/(4*sigma_yo**2)
    cc = (np.sin(theta)**2)/(2*sigma_xo**2) + (np.cos(theta)**2)/(2*sigma_yo**2)
    aaa = (np.cos(theta)**2)/(2*sigma_xoo**2) + (np.sin(theta)**2)/(2*sigma_yoo**2)
    bbb = -(np.sin(2*theta))/(4*sigma_xoo**2) + (np.sin(2*theta))/(4*sigma_yoo**2)
    ccc = (np.sin(theta)**2)/(2*sigma_xoo**2) + (np.cos(theta)**2)/(2*sigma_yoo**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + amplitude2*np.exp( - (aa*((x-xo)**2) + 2*bb*(x-xo)*(y-yo) 
                            + cc*((y-yo)**2)))+ amplitude3*np.exp( - (aaa*((x-xo)**2) + 2*bbb*(x-xo)*(y-yo) 
                            + ccc*((y-yo)**2)))
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
xoo_max=128 #x position of third peak
yoo_max=129 #y position of third peak

thet=60 #position angle of both 2d gauss in deg

amp=6 #Amplitude of max position
ampo=-1 #Amplitude of max position
ampoo=-4
sigm1=1.8 #width of first 1d gauss
sigm2=2.3 #width of second 1d gauss
sigm1o=3 #width of first 1d gauss
sigm2o=5 #width of second 1d gauss
sigm1oo=3 #width of first 1d gauss
sigm2oo=5 #width of second 1d gauss

tval=0.1 #Value offset for fit (noise)

#Fitting min-max param values:
#min:
x_max_min=126 #x position of peak
y_max_min=127 #y position of peak
xo_max_min=126 #x position of second peak
yo_max_min=127 #y position of second peak
xoo_max_min=126 #x position of second peak
yoo_max_min=127 #y position of second peak


thet_min=40 #position angle of both 2d gauss in deg

amp_min=5 #Amplitude of max position
ampo_min=-4 #Amplitude of max position
ampoo_min=-4 #Amplitude of max position
sigm1_min=0 #width of first 1d gauss
sigm2_min=0 #width of second 1d gauss
sigm1o_min=2 #width of first 1d gauss
sigm2o_min=5 #width of second 1d gauss
sigm1oo_min=2 #width of first 1d gauss
sigm2oo_min=5 #width of second 1d gauss

tval_min=0 #Value offset for fit (noise)


#max:
x_max_max=130 #x position of peak
y_max_max=131 #y position of peak
xo_max_max=130 #x position of second peak
yo_max_max=131 #y position of second peak
xoo_max_max=130 #x position of second peak
yoo_max_max=131 #y position of second peak


thet_max=70 #position angle of both 2d gauss in deg

amp_max=9 #Amplitude of max position
ampo_max=-1 #Amplitude of max position
ampoo_max=-1 #Amplitude of max position
sigm1_max=5 #width of first 1d gauss
sigm2_max=5 #width of second 1d gauss
sigm1o_max=4 #width of first 1d gauss
sigm2o_max=6 #width of second 1d gauss
sigm1oo_max=4 #width of first 1d gauss
sigm2oo_max=6 #width of second 1d gauss

tval_max=0.5 #Value offset for fit (noise)

#-------------------------------------
#build grid
thet=(thet/360.)*scipy.pi
thet_min=(thet_min/360.)*scipy.pi
thet_max=(thet_max/360.)*scipy.pi


x=np.linspace(0,x_size-1,x_size)
y=np.linspace(0,y_size-1,y_size)
x, y=np.meshgrid(x, y)

#popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,p0=(amp,ampo,x_max,y_max,sigm1,sigm1o,sigm2,sigm2o,thet,tval))

popt, pcov = opt.curve_fit(twoD_Gaussian, (x,y), dataset,maxfev=100000,p0=(amp,ampo,ampoo,x_max,y_max,sigm1,sigm1o,sigm1oo,sigm2,sigm2o,sigm2oo,thet,tval),bounds=([amp_min,ampo_min,ampoo_min, x_max_min, y_max_min, sigm1_min,sigm1o_min,sigm1oo_min, sigm2_min,sigm2o_min,sigm2oo_min, thet_min, tval_min],[amp_max,ampo_max,ampoo_max, x_max_max, y_max_max, sigm1_max,sigm1o_max,sigm1oo_max, sigm2_max,sigm2o_max, sigm2oo_max, thet_max, tval_max]))
data_fitted = twoD_Gaussian_norav((x, y), *popt)

resid_img=datasetr-data_fitted

print popt
print pcov
plt.figure("Fitfunction"); plt.imshow(data_fitted,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
plt.figure("Initial Data"); plt.imshow(datasetr,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
plt.figure("Residual"); plt.imshow(resid_img,norm=matplotlib.colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis();plt.show()
if os.path.isfile(path_out):
	os.system("rm " + path_out)
if os.path.isfile(path_out_comp):
	os.system("rm " + path_out_comp)
	
pyfits.writeto(path_out, dataset_out, header)
pyfits.writeto(path_out_comp, dataset_out_comp, header)










