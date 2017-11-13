#!/usr/bin/python
#====================================================
#This script is used to successfully create light-
#curves of one  host-star in order to detect exoplanets
#via transit method. Prepare a list ("darks_b.txt","darks_a.txt") of
#your dark images before/after, a list ("flats_b.txt","flats_a.txt") 
#of your flat images before/after and a list (data.txt) of
#your target images.
#Also: create a subfolder (./corrected/) for corrected
#target images.
#Run this script in the folder of target images!!
#====================================================

#====================================================
#          Loading necessary packages
#====================================================

import pyfits
import os
import numpy as np
from astropy.time import Time

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




#=====================================================
#                      Main
#=====================================================

#read in all your darks and flats:

hdulist=pyfits.open(path_in) * 1
dataset=hdulist[0].data
header=hdulist[0].header



# We start to calculate median flats for before/after measurement

x_size=header["NAXIS1"]
y_size=header["NAXIS2"]
	
if os.path.isfile(path_out):
	os.system("rm " + path_out)
if os.path.isfile(path_out_comp):
	os.system("rm " + path_out_comp)
	
#pyfits.writeto(path_out, dataset_out, header)
#pyfits.writeto(path_out_comp, dataset_out_comp, header)










